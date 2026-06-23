param(
  [string]$Repository = "petergraffy/environment_transplant_survival",
  [string]$Tag = "gridmet-zcta-daily-v1"
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

$token = Get-GitHubToken
$headers = @{
  Authorization = "Bearer $token"
  Accept = "application/vnd.github+json"
  "X-GitHub-Api-Version" = "2022-11-28"
  "User-Agent" = "environment-transplant-survival-gridmet-verify"
}

$release = Invoke-RestMethod -Uri "https://api.github.com/repos/$Repository/releases/tags/$Tag" -Headers $headers
$assets = Invoke-RestMethod -Uri $release.assets_url -Headers $headers
$assets |
  Select-Object name, size, browser_download_url |
  Sort-Object name |
  Format-Table -AutoSize

$totalMb = [math]::Round((($assets | Measure-Object -Property size -Sum).Sum / 1MB), 1)
Write-Host "Assets=$($assets.Count) TotalMB=$totalMb"
Write-Host "Release URL: $($release.html_url)"
