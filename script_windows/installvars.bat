
:: modify the 3 following lines if necessary
call "c:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" amd64
set SOLUTION="Visual Studio 14 Win64"

:: example with VS 2017
:: call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat" -arch=amd64
:: set SOLUTION="Visual Studio 15 Win64"

:: add path for cmake and ninja if necessary
:: ninja is just one exec that can be copied with cmake
set PATH=%PATH%;C:\Program Files\CMake\bin\;

:: modify if necessary the 8 following variables (use / in path)
:: root folder where build directories will be created
set ROOT=D:/ASTex_Devel/
:: zlib source folder 
set ZLIB_SRC=%ROOT%/zlib-1.2.11
:: OpenEXR source folder 
set OPENEXR_SRC=%ROOT%/openexr-develop
:: Itk source folder 
set ITK_SRC=%ROOT%/InsightToolkit-4.11.0
:: Itk source folder 
set ASTEX_SRC=%ROOT%/ASTex

:: release install folder
set INSTALL_REL=%ROOT%/installed-Release
:: debug install folder
set INSTALL_DBG=%ROOT%/installed-Debug
:: itk version
set ITK_VER=4.11


set ITK_INCLUDES_R=%INSTALL_REL%/include/ITK-%ITK_VER%
set ITK_INCLUDES_D=%INSTALL_DBG%/include/ITK-%ITK_VER%
