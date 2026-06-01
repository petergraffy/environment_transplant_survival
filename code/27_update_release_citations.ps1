param(
  [string]$Owner = "petergraffy",
  [string]$Repo = "environment_transplant_survival"
)

$ErrorActionPreference = "Stop"

function Get-GitHubToken {
  if ($env:GH_TOKEN) { return $env:GH_TOKEN }
  if ($env:GITHUB_TOKEN) { return $env:GITHUB_TOKEN }

  $credInput = Join-Path $env:TEMP "github-credential-input.txt"
  Set-Content -LiteralPath $credInput -Value "protocol=https`nhost=github.com`n`n" -NoNewline -Encoding ascii
  $credential = cmd /c "git credential fill < `"%TEMP%\github-credential-input.txt`""
  $passwordLine = $credential | Where-Object { $_ -like "password=*" } | Select-Object -First 1
  if (-not $passwordLine) {
    throw "Could not retrieve a GitHub token from GH_TOKEN, GITHUB_TOKEN, or Git Credential Manager."
  }
  return $passwordLine.Substring("password=".Length)
}

function Invoke-GitHubJson {
  param(
    [string]$Method,
    [string]$Uri,
    [object]$Body = $null,
    [hashtable]$Headers
  )

  $params = @{
    Method = $Method
    Uri = $Uri
    Headers = $Headers
  }
  if ($null -ne $Body) {
    $params.Body = ($Body | ConvertTo-Json -Depth 12)
    $params.ContentType = "application/json"
  }
  return Invoke-RestMethod @params
}

function Update-ReadmeAsset {
  param(
    [object]$Release,
    [string]$ReadmePath,
    [hashtable]$Headers
  )

  $existing = $Release.assets | Where-Object { $_.name -eq "README.md" } | Select-Object -First 1
  if ($existing) {
    Invoke-RestMethod -Method Delete -Uri $existing.url -Headers $Headers | Out-Null
  }

  $uploadUrl = $Release.upload_url -replace "\{\?name,label\}", "?name=README.md"
  Invoke-RestMethod `
    -Method Post `
    -Uri $uploadUrl `
    -Headers $Headers `
    -ContentType "text/markdown" `
    -InFile $ReadmePath | Out-Null
}

$token = Get-GitHubToken
$headers = @{
  Authorization = "Bearer $token"
  Accept = "application/vnd.github+json"
  "X-GitHub-Api-Version" = "2022-11-28"
  "User-Agent" = "environment-transplant-survival-release-updater"
}

$gridmetBody = @'
Daily gridMET weather exposures aggregated to 2020 ZCTA5 polygons for 2005-2025.

Assets:
- One value-only yearly Parquet file per year.
- gridmet_zcta_daily_parquet_manifest.csv with row counts, sizes, and SHA-256 checksums.
- gridmet_zcta_daily_qc_all_years.csv with missingness and fill-distance audit summaries.
- README.md with schema and citation notes.

Source citation:
- Abatzoglou, J. T. (2013). Development of gridded surface meteorological data for ecological applications and modelling. International Journal of Climatology, 33(1), 121-131. https://doi.org/10.1002/joc.3413
- gridMET data portal: https://www.northwestknowledge.net/metdata/data/

The repository contains the reproducible download, aggregation, and conversion scripts.
'@

$airPollutionBody = @'
Final filled air-pollution exposure surfaces aggregated to 2020 ZCTA5 polygons.

Assets:
- air_pollution_zcta_pm25_monthly_2005_2023.parquet: monthly PM2.5, pm25_ug_m3.
- air_pollution_zcta_o3_monthly_2005_2023.parquet: monthly ozone, o3_ppb.
- air_pollution_zcta_no2_annual_2005_2024.parquet: annual NO2, no2.
- air_pollution_zcta_parquet_manifest.csv: row counts, coverage, completeness, sizes, and SHA-256 checksums.
- fill_summary.csv and fill_audit_summary.csv: nearest-fill audit metadata.
- README.md with schema and citation notes.

Source citations:
- PM2.5: Shen, S., Li, C., van Donkelaar, A., Jacobs, N., Wang, C., & Martin, R. V. (2024). Enhancing global estimation of fine particulate matter concentrations by including geophysical a priori information in deep learning. ACS ES&T Air. https://doi.org/10.1021/acsestair.3c00054
- PM2.5 data portal: WashU Atmospheric Composition Analysis Group SatPM2.5, https://sites.wustl.edu/acag/surface-pm2-5/
- Ozone: Liu, R., Chu, L., Deziel, N. C., & Chen, K. (2026). Four-decade (1980-2023) surface ozone concentrations across the contiguous United States: Fine-resolution estimates and health implications. Environmental Science & Technology. https://doi.org/10.1021/acs.est.5c16412
- Ozone data DOI: Yale Dataverse, https://doi.org/10.60600/YU/M1WT9R
- NO2: Mohegh, A., & Anenberg, S. (2021). Global surface NO2 concentrations 1990-2020. Figshare. https://doi.org/10.6084/m9.figshare.12968114
- NO2: Nawaz, M. O. (2025). Monthly and annual US TROPOMI surface NO2 estimates (~1 km x 1 km), version 1.01. Zenodo. https://doi.org/10.5281/zenodo.14646034

PM2.5 and ozone currently cover 2005-2023 in the local source files; NO2 covers 2005-2024.
'@

$updates = @(
  @{
    Tag = "gridmet-zcta-daily-v1"
    Body = $gridmetBody
    Readme = "data/release/gridmet_zcta_daily_parquet/README.md"
  },
  @{
    Tag = "air-pollution-zcta-v1"
    Body = $airPollutionBody
    Readme = "data/release/air_pollution_zcta_parquet/README.md"
  }
)

foreach ($update in $updates) {
  $tag = $update.Tag
  Write-Host "Updating release $tag"
  $releaseUri = "https://api.github.com/repos/$Owner/$Repo/releases/tags/$tag"
  $release = Invoke-GitHubJson -Method Get -Uri $releaseUri -Headers $headers
  $release = Invoke-GitHubJson -Method Patch -Uri $release.url -Headers $headers -Body @{ body = $update.Body }
  Update-ReadmeAsset -Release $release -ReadmePath $update.Readme -Headers $headers
  Write-Host "Updated release $tag"
}
