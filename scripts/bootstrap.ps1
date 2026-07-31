$ErrorActionPreference = "Stop"
$Python = if ($env:PYTHON_BIN) { $env:PYTHON_BIN } else { "py" }

if ($Python -eq "py") {
    & py -3.11 -m venv .venv
} else {
    & $Python -m venv .venv
}

& .\.venv\Scripts\python.exe -m pip install --upgrade pip setuptools wheel
& .\.venv\Scripts\python.exe -m pip install -e ".[dev]"
& .\.venv\Scripts\python.exe -m pytest -q
Write-Host "Environment ready. Select .venv\\Scripts\\python.exe in VS Code."
Write-Host "Next command: cml-repairnet doctor"
