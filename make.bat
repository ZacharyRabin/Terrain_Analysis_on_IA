@echo off

rem This directory should contain the top-level CMakeLists.txt - it is assumed to reside in the same
rem directory as this script.
set CMAKE_DIR=%~dp0

rem Get the current directory - this is the default location for the build and install directory.
set CURRENT_DIR=%cd%

rem The build directory.
set BUILD_DIR=%CURRENT_DIR%\build

rem The install directory.
set INSTALL_DIR=%CURRENT_DIR%\install

if exist "%BUILD_DIR%" goto BUILD_DIR_CREATED
    mkdir "%BUILD_DIR%"
:BUILD_DIR_CREATED

if exist "%INSTALL_DIR%" goto INSTALL_DIR_CREATED
    mkdir "%INSTALL_DIR%"
:INSTALL_DIR_CREATED

set Eigen_DIR="C:\VESTEC\install\windows-externals-release"
rem set Eigen_DIR="C:\VESTEC\build\windows-externals-release\eigen"
set CMAKE_FLAGS=-G "Visual Studio 15"
set BUILD_TYPE=Release

cd "%BUILD_DIR%"
cmake %CMAKE_FLAGS% -DCMAKE_INSTALL_PREFIX="%INSTALL_DIR%"^
    -DCMAKE_BUILD_TYPE=%BUILD_TYPE%^
    -DEigen3_DIR="%Eigen_DIR%\share\eigen3\cmake"^
    -DEIGEN3_INCLUDE_DIR="%Eigen_DIR%\include\eigen3"^
    "%CMAKE_DIR%" || exit /b

rem cmake --build . --config %BUILD_TYPE% 
rem || exit /b
cmake --build . --config %BUILD_TYPE% --target ia_terrain_analysis --parallel 8  || exit /b
cd "%CURRENT_DIR%"
@echo on