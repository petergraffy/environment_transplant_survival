param(
  [string]$Repository = "petergraffy/environment_transplant_survival",
  [string]$Tag = "gridmet-zcta-daily-v1",
  [string]$Target = "main",
  [string]$AssetDir = "data/release/gridmet_zcta_daily_parquet",
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
  if (-not $passwordLine) {
    throw "Could not retrieve a GitHub token from GH_TOKEN, GITHUB_TOKEN, or Git Credential Manager."
  }
  return $passwordLine.Substring("password=".Length)
}

function Invoke-GitHubJson {
  param(
    [string]$Method,
    [string]$Uri,
    [hashtable]$Headers,
    [object]$Body = $null
  )
  if ($null -eq $Body) {
    return Invoke-RestMethod -Method $Method -Uri $Uri -Headers $Headers
  }
  return Invoke-RestMethod -Method $Method -Uri $Uri -Headers $Headers -Body ($Body | ConvertTo-Json -Depth 20)
}

$token = Get-GitHubToken
$headers = @{
  Authorization = "Bearer $token"
  Accept = "application/vnd.github+json"
  "X-GitHub-Api-Version" = "2022-11-28"
  "User-Agent" = "environment-transplant-survival-gridmet-upload"
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
    name = "gridMET ZCTA daily Parquet assets v1"
    body = @"
Daily gridMET weather exposures aggregated to 2020 ZCTA5 polygons for 2005-2025.

Assets:
- `gridmet_zcta_daily_YYYY.parquet`: value-only yearly daily ZCTA exposures.
- `gridmet_zcta_daily_parquet_manifest.csv`: file sizes, row counts, and SHA-256 checksums.
- `gridmet_zcta_daily_qc_all_years.csv`: completeness and distribution QC.
- `README.md`: variable notes.

The repository contains the reproducible download/aggregation/conversion scripts.
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
