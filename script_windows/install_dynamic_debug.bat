call installvars.bat

set PATH=%PATH%;%INSTALL_DBG%\lib;%INSTALL_DBG%\bin

cd %ROOT%
if exist build-zlib-debug del /S /Q build-zlib-debug
mkdir build-zlib-debug
cd build-zlib-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG%  %ZLIB_SRC%
cmake --build . --config Debug --target install

cd %ROOT%
if exist build-ilmbase-dyn-debug del /S /Q build-ilmbase-dyn-debug
mkdir build-ilmbase-dyn-debug
cd build-ilmbase-dyn-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% ^
 -DBUILD_SHARED_LIBS=ON -DNAMESPACE_VERSIONING=OFF %OPENEXR_SRC%/IlmBase
cmake --build . --config Debug --target install

cd %ROOT%
if exist build-openexr-dyn-debug del /S /Q build-openexr-dyn-debug
mkdir build-openexr-dyn-debug
cd build-openexr-dyn-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DCMAKE_PREFIX_PATH=%INSTALL_DBG% ^
 -DBUILD_SHARED_LIBS=ON -DILMBASE_PACKAGE_PREFIX=%INSTALL_DBG% -DNAMESPACE_VERSIONING=OFF %OPENEXR_SRC%/OpenEXR
cmake --build . --config Debug --target install

copy %INSTALL_DBG:/=\%\lib\*.dll %INSTALL_DBG:/=\%\bin 


cd %ROOT%
if exist build-itk-dyn-debug del /S /Q build-itk-dyn-debug
mkdir build-itk-dyn-debug
cd build-itk-dyn-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% ^
 -DCMAKE_PREFIX_PATH=%INSTALL_DBG% -DITK_USE_SYSTEM_ZLIB=ON -DBUILD_SHARED_LIBS=ON %ITK_SRC%
cmake --build . --config Debug --target install

cd %ROOT%
if exist build-astex-dyn-debug del /S /Q build-astex-dyn-debug
mkdir build-astex-dyn-debug
cd build-astex-dyn-debug
cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Debug" -DCMAKE_PREFIX_PATH=%INSTALL_DBG% ^
 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_D%/itkpng -DPNG_LIBRARY_DEBUG=%INSTALL_DBG%/lib/itkpng-%ITK_VER%.lib ^
 -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DBUILD_SHARED_LIBS=ON %ASTEX_SRC%

copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_D:/=\%\itkpng\

cmake --build . --config Debug

echo You can remove build-zlib-debug build-ilmbase-dyn-debug build-openexr-dyn-debug build-itk-dyn-debug
pause

