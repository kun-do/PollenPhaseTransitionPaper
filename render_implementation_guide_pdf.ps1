param(
    [string]$HtmlPath = (Join-Path $PSScriptRoot 'IMPLEMENTATION_GUIDE_PRINT.html'),
    [string]$PdfPath = (Join-Path $PSScriptRoot 'IMPLEMENTATION_GUIDE.pdf')
)

$edgeCandidates = @(
    'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
    'C:\Program Files\Microsoft\Edge\Application\msedge.exe',
    'C:\Program Files\Google\Chrome\Application\chrome.exe',
    'C:\Program Files (x86)\Google\Chrome\Application\chrome.exe'
)

$browser = $edgeCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1

if (-not $browser) {
    throw 'No supported headless browser was found. Install Microsoft Edge or Google Chrome.'
}

$htmlFullPath = (Resolve-Path $HtmlPath).Path
$pdfFullPath = [System.IO.Path]::GetFullPath($PdfPath)
$fileUri = 'file:///' + ($htmlFullPath -replace '\\', '/')

& $browser `
    --headless `
    --disable-gpu `
    --allow-file-access-from-files `
    --run-all-compositor-stages-before-draw `
    --virtual-time-budget=15000 `
    --print-to-pdf-no-header `
    "--print-to-pdf=$pdfFullPath" `
    $fileUri

if (-not (Test-Path $pdfFullPath)) {
    throw "PDF generation failed: $pdfFullPath was not created."
}

Write-Output $pdfFullPath
