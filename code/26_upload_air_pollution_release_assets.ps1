param(
  [string]$Repository = "petergraffy/environment_transplant_survival",
  [string]$Tag = "air-pollution-zcta-v1",
  [string]$Target = "main",
  [string]$AssetDir = "data/release/air_pollution_zcta_parquet",
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
  "User-Agent" = "environment-transplant-survival-air-pollution-upload"
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
    name = "Air pollution ZCTA Parquet assets v1"
    body = @"
Final filled air-pollution exposure surfaces aggregated to 2020 ZCTA5 polygons.

Assets:
- `air_pollution_zcta_pm25_monthly_2005_2023.parquet`: monthly PM2.5, `pm25_ug_m3`.
- `air_pollution_zcta_o3_monthly_2005_2023.parquet`: monthly ozone, `o3_ppb`.
- `air_pollution_zcta_no2_annual_2005_2024.parquet`: annual NO2, `no2`.
- `air_pollution_zcta_parquet_manifest.csv`: row counts, completeness, sizes, and SHA-256 checksums.
- `fill_summary.csv` and `fill_audit_summary.csv`: nearest-fill audit metadata.

PM2.5 and ozone currently cover 2005-2023 in the local source files; NO2 covers 2005-2024.
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
  Invoke-RestMethod `
    -Method POST `
    -Uri $uploadUri `
    -Headers $headers `
    -ContentType "application/octet-stream" `
    -InFile $file.FullName | Out-Null
}

Write-Host "Release URL: $($release.html_url)"
