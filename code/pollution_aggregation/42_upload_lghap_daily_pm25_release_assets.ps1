param(
  [string]$Repository = "petergraffy/environment_transplant_survival",
  [string]$Tag = "lghap-pm25-zcta-daily-v1",
  [string]$Target = "main",
  [string]$AssetDir = "data/release/lghap_pm25_zcta_daily_parquet",
  [switch]$ReplaceExisting
)

$ErrorActionPreference = "Stop"

function Get-GitHubToken {
  if ($env:GH_TOKEN) { return $env:GH_TOKEN }
  if ($env:GITHUB_TOKEN) { return $env:GITHUB_TOKEN }
  $credInput = Join-Path $env:TEMP "github-credential-input.txt"
  Set-Content -LiteralPath $credInput -Value "protocol=https`nhost=github.com`n`n" -NoNewline -Encoding ascii
  $output = cmd /c "git credential fill < `"%TEMP%\github-credential-input.txt`""
  $passwordLine = $output | Where-Object { $_ -like "password=*" } | Select-Object -First 1
  if (-not $passwordLine) { throw "Could not retrieve a GitHub token." }
  return $passwordLine.Substring("password=".Length)
}

function Invoke-GitHubJson {
  param([string]$Method, [string]$Uri, [hashtable]$Headers, [object]$Body = $null)
  if ($null -eq $Body) { return Invoke-RestMethod -Method $Method -Uri $Uri -Headers $Headers }
  return Invoke-RestMethod -Method $Method -Uri $Uri -Headers $Headers -Body ($Body | ConvertTo-Json -Depth 20)
}

$token = Get-GitHubToken
$headers = @{
  Authorization = "Bearer $token"
  Accept = "application/vnd.github+json"
  "X-GitHub-Api-Version" = "2022-11-28"
  "User-Agent" = "environment-transplant-survival-lghap-pm25-upload"
}
$jsonHeaders = $headers.Clone()
$jsonHeaders["Content-Type"] = "application/json"

$apiRoot = "https://api.github.com/repos/$Repository"
$release = $null
try {
  $release = Invoke-GitHubJson -Method GET -Uri "$apiRoot/releases/tags/$Tag" -Headers $headers
  Write-Host "Using existing release $Tag"
} catch {
  if ($_.Exception.Response.StatusCode.value__ -ne 404) { throw }
  Write-Host "Creating release $Tag"
  $body = @{
    tag_name = $Tag
    target_commitish = $Target
    name = "LGHAP daily PM2.5 ZCTA Parquet assets v1"
    body = @"
Daily LGHAP PM2.5 exposures aggregated to 2020 ZCTA5 polygons for CONUS, 2005-2021.

Assets:
- `lghap_pm25_zcta_daily_YYYY.parquet`: yearly daily ZCTA PM2.5 values, one row per ZCTA-date.
- `lghap_pm25_zcta_daily_parquet_manifest.csv`: row counts, date coverage, completeness, sizes, and SHA-256 checksums.
- `lghap_pm25_zcta_daily_qc_all_years.csv`: monthly completeness and distribution QC.
- `README.md`: variables and citation notes.

Citation:
Bai, K., Li, K., Shao, L., Li, X., Liu, C., Li, Z., Ma, M., Han, D., Sun, Y., Zheng, Z., Li, R., Chang, N.-B., and Guo, J.: LGHAP v2: a global gap-free aerosol optical depth and PM2.5 concentration dataset since 2000 derived via big Earth data analytics, Earth Syst. Sci. Data, 16, 2425-2448, https://doi.org/10.5194/essd-16-2425-2024, 2024.

The repository contains the reproducible daily NetCDF-to-ZCTA aggregation and release packaging scripts.
"@
    draft = $false
    prerelease = $false
  }
  $release = Invoke-GitHubJson -Method POST -Uri "$apiRoot/releases" -Headers $jsonHeaders -Body $body
}

$uploadRoot = $release.upload_url.Split("{")[0]
$existingAssets = Invoke-GitHubJson -Method GET -Uri $release.assets_url -Headers $headers
$assetFiles = Get-ChildItem -LiteralPath $AssetDir -File |
  Where-Object { $_.Name -match "\.(parquet|csv|md)$" } |
  Sort-Object Name

foreach ($file in $assetFiles) {
  $existing = $existingAssets | Where-Object { $_.name -eq $file.Name } | Select-Object -First 1
  if ($existing) {
    if (-not $ReplaceExisting -and [int64]$existing.size -eq [int64]$file.Length) {
      Write-Host "Skipping existing asset $($file.Name)"
      continue
    }
    if (-not $ReplaceExisting) {
      Write-Host "Skipping existing asset with different size $($file.Name); pass -ReplaceExisting to overwrite"
      continue
    }
    Write-Host "Deleting existing asset $($file.Name)"
    Invoke-GitHubJson -Method DELETE -Uri "$apiRoot/releases/assets/$($existing.id)" -Headers $headers | Out-Null
  }

  $encodedName = [System.Uri]::EscapeDataString($file.Name)
  $uploadUri = "${uploadRoot}?name=$encodedName"
  Write-Host "Uploading $($file.Name) ($([math]::Round($file.Length / 1MB, 1)) MB)"
  $uploaded = $false
  for ($attempt = 1; $attempt -le 4 -and -not $uploaded; $attempt++) {
    try {
      Invoke-RestMethod `
        -Method POST `
        -Uri $uploadUri `
        -Headers $headers `
        -ContentType "application/octet-stream" `
        -InFile $file.FullName | Out-Null
      $uploaded = $true
    } catch {
      Write-Host "Upload attempt $attempt failed for $($file.Name): $($_.Exception.Message)"
      Start-Sleep -Seconds ([math]::Min(60, 10 * $attempt))
      $existingAssets = Invoke-GitHubJson -Method GET -Uri $release.assets_url -Headers $headers
      $existing = $existingAssets | Where-Object { $_.name -eq $file.Name } | Select-Object -First 1
      if ($existing -and [int64]$existing.size -eq [int64]$file.Length) {
        Write-Host "Detected completed asset after retry check $($file.Name)"
        $uploaded = $true
      } elseif ($existing -and $ReplaceExisting) {
        Write-Host "Deleting incomplete asset $($file.Name)"
        Invoke-GitHubJson -Method DELETE -Uri "$apiRoot/releases/assets/$($existing.id)" -Headers $headers | Out-Null
      }
      if ($attempt -eq 4 -and -not $uploaded) { throw }
    }
  }
}

Write-Host "Release URL: $($release.html_url)"
