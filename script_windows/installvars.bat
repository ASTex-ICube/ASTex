
:: modify the 3 following lines if necessary

:: example with VS 2015 (standard install)
::call "c:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" amd64
::set SOLUTION="Visual Studio 14 Win64"
::set IDE="C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\IDE\devenv.exe"

:: example with VS 2017 (standard install)
call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\Tools\VsDevCmd.bat" -arch=amd64
set SOLUTION="Visual Studio 15 Win64"
set IDE="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\IDE\devenv.exe"
set GIT="C:\Program Files\Git\bin\git.exe"


set JOMGEN="NMake Makefiles JOM"

:: add path for cmake and jom (qt) if necessary
set PATH=%PATH%;C:\Program Files\CMake\bin\;C:\Qt\Tools\QtCreator\bin\;

:: modify if necessary the 8 following variables (use / in path)
:: root folder where build directories will be created
set ROOT=D:/ASTex_Devel/

:: zlib source folder 
set ZLIB_SRC=%ROOT%/zlib
:: OpenEXR source folder 
set OPENEXR_SRC=%ROOT%/openexr
:: Itk source folder 
set ITK_SRC=%ROOT%/ITK
:: itk version
set ITK_VER=4.13

:: ASTex source folder 
set ASTEX_SRC=%ROOT%/ASTex

:: release install folder
set INSTALL_REL=%ROOT%/installed-Release
:: debug install folder
set INSTALL_DBG=%ROOT%/installed-Debug



set ITK_INCLUDES_R=%INSTALL_REL%/include/ITK-%ITK_VER%
set ITK_INCLUDES_D=%INSTALL_DBG%/include/ITK-%ITK_VER%
