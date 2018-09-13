call installvars.bat

set PATH=%PATH%;%INSTALL_REL%\lib;%INSTALL_REL%\bin

cd /d %ROOT%
if exist build-zlib-release del /S /Q build-zlib-release
mkdir build-zlib-release
cd build-zlib-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL%  %ZLIB_SRC% || exit /b
cmake --build . --config Release --target install || exit /b

cd /d %ROOT%
if exist build-openexr-static-release del /S /Q build-openexr-static-release
mkdir build-openexr-static-release
cd build-openexr-static-release
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% -DCMAKE_PREFIX_PATH=%INSTALL_REL%^
 -DILMBASE_PACKAGE_PREFIX=%INSTALL_REL% -DOPENEXR_NAMESPACE_VERSIONING=OFF -DOPENEXR_BUILD_STATIC=ON -DOPENEXR_BUILD_SHARED=OFF^
 -DOPENEXR_BUILD_ILMBASE=ON -DOPENEXR_BUILD_OPENEXR=ON -DOPENEXR_BUILD_PYTHON_LIBS=OFF -DOPENEXR_BUILD_UTILS=OFF^
 -DOPENEXR_BUILD_TESTS=OFF -DOPENEXR_BUILD_VIEWERS=OFF  %OPENEXR_SRC% || exit /b
cmake --build . --config Release --target install || exit /b


cd /d %ROOT%
del /S /Q build-itk-static-release
mkdir build-itk-static-release
cd build-itk-static-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% -DITK_USE_SYSTEM_ZLIB=ON ^
 -DZLIB_INCLUDE_DIR=%INSTALL_REL%/include -DZLIB_LIBRARY_RELEASE=%INSTALL_REL%/lib/zlibstatic.lib  %ITK_SRC% || exit /b
cmake --build . --config Release --target install || exit /b

cd /d %ROOT%
del /S /Q build-astex-static-release
mkdir build-astex-static-release
cd build-astex-static-release
cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Release" -DCMAKE_PREFIX_PATH=%INSTALL_REL% ^
 -DZLIB_INCLUDE_DIR=%INSTALL_REL%/include -DZLIB_LIBRARY_RELEASE=%INSTALL_REL%/lib/zlibstatic.lib ^
 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_R%/itkpng -DPNG_LIBRARY_RELEASE=%INSTALL_REL%/lib/itkpng-%ITK_VER%.lib ^
 -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% %ASTEX_SRC% || exit /b
 
copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_R:/=\%\itkpng\

echo You can remove build-zlib-release build-ilmbase-static-release build-openexr-static-release build-itk-static-release

cmake-gui .

pause
