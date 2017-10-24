call installvars.bat

cd %ROOT%
del /S /Q build-zlib-release
mkdir build-zlib-release
cd build-zlib-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL%  %ZLIB_SRC%
cmake --build . --config Release --target install

cd %ROOT%
del /S /Q build-ilmbase-static-release
mkdir build-ilmbase-static-release
cd build-ilmbase-static-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% ^
 -DBUILD_SHARED_LIBS=OFF -DNAMESPACE_VERSIONING=OFF %OPENEXR_SRC%/IlmBase
cmake --build . --config Release --target install

cd %ROOT%
del /S /Q build-openexr-static-release
mkdir build-openexr-static-release
cd build-openexr-static-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% ^
 -DZLIB_INCLUDE_DIR=%INSTALL_REL%/include -DZLIB_LIBRARY_RELEASE=%INSTALL_REL%/lib/zlibstatic.lib ^
 -DBUILD_SHARED_LIBS=OFF -DILMBASE_PACKAGE_PREFIX=%INSTALL_REL% -DNAMESPACE_VERSIONING=OFF %OPENEXR_SRC%/OpenEXR
cmake --build . --config Release --target install

cd %ROOT%
del /S /Q build-itk-static-release
mkdir build-itk-static-release
cd build-itk-static-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% -DITK_USE_SYSTEM_ZLIB=ON ^
 -DZLIB_INCLUDE_DIR=%INSTALL_REL%/include -DZLIB_LIBRARY_RELEASE=%INSTALL_REL%/lib/zlibstatic.lib  %ITK_SRC%
cmake --build . --config Release --target install

cd %ROOT%
del /S /Q build-astex-static-release
mkdir build-astex-static-release
cd build-astex-static-release
cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Release" -DCMAKE_PREFIX_PATH=%INSTALL_REL% ^
 -DZLIB_INCLUDE_DIR=%INSTALL_REL%/include -DZLIB_LIBRARY_RELEASE=%INSTALL_REL%/lib/zlibstatic.lib ^
 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_R%/itkpng -DPNG_LIBRARY_RELEASE=%INSTALL_REL%/lib/itkpng-%ITK_VER%.lib ^
 -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% %ASTEX_SRC%
 
copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_R:/=\%\itkpng\

echo You can remove build-zlib-release build-ilmbase-static-release build-openexr-static-release build-itk-static-release

cmake-gui .

pause
