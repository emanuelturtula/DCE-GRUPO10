param (
    [Parameter(Mandatory=$true)]
    [Alias('f')]
    $Files
)

Foreach ($File in $Files) {
    $FFT = Get-Content $File -Raw -Encoding ASCII
    $FFT.replace("`t(", '|').replace("dB,",'|').replace("?)",'|') | Set-Content -Path "$File" 
}
