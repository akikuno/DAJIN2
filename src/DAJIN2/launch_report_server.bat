@echo off
setlocal
cd /d "%~dp0"

where python >nul 2>nul
if errorlevel 1 (
    echo Python is not available in PATH.
    pause
    exit /b 1
)

set PORT=8000
start "" "http://127.0.0.1:%PORT%/report.html"
python -m http.server %PORT%

endlocal
