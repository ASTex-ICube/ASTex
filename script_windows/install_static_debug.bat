call installvars.bat

cd %ROOT%
del /S /Q build-zlib-debug
mkdir build-zlib-debug
cd build-zlib-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG%  %ZLIB_SRC%
cmake --build . --config Debug --target install

cd %ROOT%
del /S /Q build-ilmbase-static-debug
mkdir build-ilmbase-static-debug
cd build-ilmbase-static-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% ^
 -DBUILD_SHARED_LIBS=OFF -DNAMESPACE_VERSIONING=OFF %OPENEXR_SRC%/IlmBase
cmake --build . --config Debug --target install

cd %ROOT%
del /S /Q build-openexr-static-debug
mkdir build-openexr-static-debug
cd build-openexr-static-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% ^
 -DZLIB_INCLUDE_DIR=%INSTALL_DBG%/include -DZLIB_LIBRARY_DEBUG=%INSTALL_DBG%/lib/zlibstaticd.lib ^
 -DBUILD_SHARED_LIBS=OFF -DILMBASE_PACKAGE_PREFIX=%INSTALL_DBG% -DNAMESPACE_VERSIONING=OFF %OPENEXR_SRC%/OpenEXR
cmake --build . --config Debug --target install


cd %ROOT%
del /S /Q build-itk-static-debug
mkdir build-itk-static-debug
cd build-itk-static-debug
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DITK_USE_SYSTEM_ZLIB=ON ^
 -DZLIB_INCLUDE_DIR=%INSTALL_DBG%/include -DZLIB_LIBRARY_DEBUG=%INSTALL_DBG%/lib/zlibstaticd.lib  %ITK_SRC%
cmake --build . --config Debug --target install

cd %ROOT%
del /S /Q build-astex-static-debug
mkdir build-astex-static-debug
cd build-astex-static-debug
cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Debug" -DCMAKE_PREFIX_PATH=%INSTALL_DBG% ^
 -DZLIB_INCLUDE_DIR=%INSTALL_DBG%/include -DZLIB_LIBRARY_DEBUG=%INSTALL_DBG%/lib/zlibstaticd.lib ^
 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_D%/itkpng -DPNG_LIBRARY_DEBUG=%INSTALL_DBG%/lib/itkpng-%ITK_VER%.lib ^
 -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% %ASTEX_SRC%

copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_D:/=\%\itkpng\

cmake --build . --config Debug

echo You can remove build-zlib-debug build-ilmbase-static-debug build-openexr-static-debug build-itk-static-debug
pause

