param(
  [string]$Repository = "petergraffy/environment_transplant_survival",
  [string]$Tag = "lghap-pm25-zcta-daily-v1",
  [string]$ReadmePath = "data/release/lghap_pm25_zcta_daily_parquet/README.md"
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
  "User-Agent" = "environment-transplant-survival-lghap-pm25-metadata"
}
$jsonHeaders = $headers.Clone()
$jsonHeaders["Content-Type"] = "application/json"

$apiRoot = "https://api.github.com/repos/$Repository"
$release = Invoke-GitHubJson -Method GET -Uri "$apiRoot/releases/tags/$Tag" -Headers $headers

$body = @"
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

Invoke-GitHubJson -Method PATCH -Uri "$apiRoot/releases/$($release.id)" -Headers $jsonHeaders -Body @{ body = $body } | Out-Null
Write-Host "Updated release notes for $Tag"

$assets = Invoke-GitHubJson -Method GET -Uri $release.assets_url -Headers $headers
$existingReadme = $assets | Where-Object { $_.name -eq "README.md" } | Select-Object -First 1
if ($existingReadme) {
  Invoke-GitHubJson -Method DELETE -Uri "$apiRoot/releases/assets/$($existingReadme.id)" -Headers $headers | Out-Null
  Write-Host "Deleted existing README.md asset"
}

$uploadRoot = $release.upload_url.Split("{")[0]
$encodedName = [System.Uri]::EscapeDataString("README.md")
$uploadUri = "${uploadRoot}?name=$encodedName"
Invoke-RestMethod `
  -Method POST `
  -Uri $uploadUri `
  -Headers $headers `
  -ContentType "application/octet-stream" `
  -InFile $ReadmePath | Out-Null

Write-Host "Uploaded README.md asset"
Write-Host "Release URL: $($release.html_url)"
