Write-Host "[TASK] Clearing contents of ./build-dbg/ and ./build-opt/" -ForegroundColor Green
rm -r ./build-dbg/*
rm -r ./build-opt/*
Write-Host "[SUCCESS]==="

Write-Host "[TASK] Switching directories to build-dbg" -ForegroundColor Green
cd build-dbg
Write-Host "[SUCCESS]==="

Write-Host "[TASK] Configuring CMAKE for DEBUG" -ForegroundColor Green
cmake .. -DCMAKE_BUILD_TYPE=debug -G"Visual Studio 15 2017 Win64" -DFREEGLUT_DIR=%cd%"\..\third-party\freeglut\bin-win64-msvc2017\debug" -DGLUI_DIR=%cd%"\..\third-party\glui\bin-win64-msvc2017\debug" -DEIGEN3_DIR=%cd%"\..\third-party\eigen"
Write-Host "[SUCCESS]==="

Write-Host "[TASK] building" -ForegroundColor Green
cmake --build . --config debug
Write-Host "[SUCCESS]==="

Write-Host "[TASK] copying target dlls" -ForegroundColor Green
cmake --build . --config debug --target COPY_DLLS
cd ..
Write-Host "[SUCCESS]"

Write-Host "[TASK] Configuring CMAKE for RELEASE" -ForegroundColor Green
cd build-opt
cmake .. -DCMAKE_BUILD_TYPE=relwithdebinfo  -G"Visual Studio 15 2017 Win64" -DFREEGLUT_DIR=%cd%"/../third-party/freeglut/bin-win64-msvc2017/release"  -DGLUI_DIR=%cd%"/../third-party/glui/bin-win64-msvc2017/release" -DEIGEN3_DIR=%cd%"/../third-party/eigen"
Write-Host "[SUCCESS]==="

Write-Host "[TASK] building" -ForegroundColor Green
cmake --build . --config relwithdebinfo
Write-Host "[SUCCESS]==="

Write-Host "[TASK] copying target dlls" -ForegroundColor Green
cmake --build . --config relwithdebinfo --target COPY_DLLS
Write-Host "[SUCCESS]"

Write-Host "[TASK] Starting Application [mode=release]" -ForegroundColor Green
cd ..
.\build-opt\bin\minimeshgui.exe
Write-Host "[SUCCESS]==="

Write-Host "[ACTION] Exiting" -ForegroundColor Green


