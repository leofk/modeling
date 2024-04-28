Write-Host "[TASK] Switching directories to build-dbg" -ForegroundColor Green
cd build-dbg
Write-Host "[SUCCESS]==="

Write-Host "[TASK] Building Debug" -ForegroundColor Green
cmake --build . --config debug
Write-Host "[SUCCESS]==="

Write-Host "[TASK] Starting Application" -ForegroundColor Green
.\bin\minimeshgui.exe
Write-Host "[SUCCESS]==="

cd ..