$ErrorActionPreference = 'Stop'
$specs = @(
    @{ Tag = 'o3-zcta-daily-v1'; Directory = 'o3_zcta_daily_parquet'; Prefix = 'o3_zcta_daily_' },
    @{ Tag = 'lghap-pm25-zcta-daily-v1'; Directory = 'lghap_pm25_zcta_daily_parquet'; Prefix = 'lghap_pm25_zcta_daily_' }
)
$records = foreach ($spec in $specs) {
    $release = Invoke-RestMethod -Uri "https://api.github.com/repos/petergraffy/environment_transplant_survival/releases/tags/$($spec.Tag)"
    foreach ($year in 2005..2024) {
        $name = "$($spec.Prefix)$year.parquet"
        $asset = @($release.assets | Where-Object { $_.name -eq $name })
        if ($asset.Count -ne 1) { throw "Release asset missing or ambiguous: $name" }
        $path = Join-Path "data/release/$($spec.Directory)" $name
        $hash = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
        if ("sha256:$hash" -ne $asset[0].digest) { throw "Published checksum mismatch: $path" }
        [pscustomobject]@{ release = $release.html_url; file = $name; year = $year; sha256 = $hash; verified = $true }
    }
}
$outDir = 'output/exposure_input_validation'
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
$records | Export-Csv -LiteralPath "$outDir/published_daily_release_checksums.csv" -NoTypeInformation
Write-Output "Verified all $($records.Count) annual Parquet files against GitHub release SHA-256 digests."
