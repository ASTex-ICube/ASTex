call installvars.bat

set PATH=%PATH%;%INSTALL_DBG%\lib;%INSTALL_DBG%\bin

cd /d %ROOT%
if exist build-zlib-debug del /S /Q build-zlib-debug
mkdir build-zlib-debug
cd build-zlib-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG%  %ZLIB_SRC% || exit /b
cmake --build . --config Debug --target install || exit /b

cd /d %ROOT%
if exist build-openexr-dyn-debug del /S /Q build-openexr-dyn-debug
mkdir build-openexr-dyn-debug
cd build-openexr-dyn-debug
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DCMAKE_PREFIX_PATH=%INSTALL_DBG%^
 -DILMBASE_PACKAGE_PREFIX=%INSTALL_DBG% -DOPENEXR_NAMESPACE_VERSIONING=OFF -DOPENEXR_BUILD_STATIC=OFF -DOPENEXR_BUILD_SHARED=ON^
 -DOPENEXR_BUILD_ILMBASE=ON -DOPENEXR_BUILD_OPENEXR=ON -DOPENEXR_BUILD_PYTHON_LIBS=OFF -DOPENEXR_BUILD_UTILS=OFF^
 -DOPENEXR_BUILD_TESTS=OFF -DOPENEXR_BUILD_VIEWERS=OFF  %OPENEXR_SRC%  || exit /b
cmake --build . --config Debug --target install  || exit /b

copy OpenEXR\IlmImf\*.dll %INSTALL_DBG:/=\%\bin 

cd /d %ROOT%
if exist build-itk-dyn-debug del /S /Q build-itk-dyn-debug
mkdir build-itk-dyn-debug
cd build-itk-dyn-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% ^
 -DCMAKE_PREFIX_PATH=%INSTALL_DBG% -DITK_USE_SYSTEM_ZLIB=ON -DBUILD_SHARED_LIBS=ON %ITK_SRC%  || exit /b
cmake --build . --config Debug --target install  || exit /b


cd /d %ROOT%
if exist build-astex-dyn-debug del /S /Q build-astex-dyn-debug
mkdir build-astex-dyn-debug
cd build-astex-dyn-debug
cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Debug" -DCMAKE_PREFIX_PATH=%INSTALL_DBG% ^
 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_R%/itkpng -DPNG_LIBRARY_DBGEASE=%INSTALL_DBG%/lib/itkpng-%ITK_VER%.lib ^
 -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DBUILD_SHARED_LIBS=ON %ASTEX_SRC%  || exit /b

copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_R:/=\%\itkpng\

echo You can remove build-zlib-debug build-ilmbase-dyn-debug build-openexr-dyn-debug build-itk-dyn-debug

cmake-gui . 

pause
