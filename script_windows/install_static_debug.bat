call installvars.bat

set PATH=%PATH%;%INSTALL_DBG%\lib;%INSTALL_DBG%\bin

cd /d %ROOT%
if exist build-zlib-debug del /S /Q build-zlib-debug
mkdir build-zlib-debug
cd build-zlib-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG%  %ZLIB_SRC% || exit /b
cmake --build . --config Debug --target install || exit /b

cd /d %ROOT%
if exist build-openexr-static-debug del /S /Q build-openexr-static-debug
mkdir build-openexr-static-debug
cd build-openexr-static-debug
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DCMAKE_PREFIX_PATH=%INSTALL_DBG%^
 -DILMBASE_PACKAGE_PREFIX=%INSTALL_DBG% -DOPENEXR_NAMESPACE_VERSIONING=OFF -DOPENEXR_BUILD_STATIC=ON -DOPENEXR_BUILD_SHARED=OFF^
 -DOPENEXR_BUILD_ILMBASE=ON -DOPENEXR_BUILD_OPENEXR=ON -DOPENEXR_BUILD_PYTHON_LIBS=OFF -DOPENEXR_BUILD_UTILS=OFF^
 -DOPENEXR_BUILD_TESTS=OFF -DOPENEXR_BUILD_VIEWERS=OFF  %OPENEXR_SRC% || exit /b
cmake --build . --config Debug --target install || exit /b


cd /d %ROOT%
if exist build-itk-static-debug del /S /Q build-itk-static-debug
mkdir build-itk-static-debug
cd build-itk-static-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DITK_USE_SYSTEM_ZLIB=ON ^
 -DZLIB_INCLUDE_DIR=%INSTALL_DBG%/include -DZLIB_LIBRARY_DEBUG=%INSTALL_DBG%/lib/zlibstaticd.lib  %ITK_SRC% || exit /b
cmake --build . --config Debug --target install || exit /b


cd /d %ROOT%
if exist build-astex-static-debug del /S /Q build-astex-static-debug
mkdir build-astex-static-debug
cd build-astex-static-debug
cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Debug" -DCMAKE_PREFIX_PATH=%INSTALL_DBG% ^
 -DZLIB_INCLUDE_DIR=%INSTALL_DBG%/include -DZLIB_LIBRARY_DEBUG=%INSTALL_DBG%/lib/zlibstaticd.lib ^
 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_D%/itkpng -DPNG_LIBRARY_DEBUG=%INSTALL_DBG%/lib/itkpng-%ITK_VER%.lib ^
 -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% %ASTEX_SRC% || exit /b

copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_R:/=\%\itkpng\

echo You can remove build-zlib-debug build-ilmbase-static-debug build-openexr-static-debug build-itk-static-debug

cmake-gui . 

pause