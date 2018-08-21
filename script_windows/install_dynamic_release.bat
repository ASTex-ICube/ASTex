call installvars.bat

set PATH=%PATH%;%INSTALL_REL%\lib;%INSTALL_REL%\bin

cd %ROOT%
if exist build-zlib-release del /S /Q build-zlib-release
mkdir build-zlib-release
cd build-zlib-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL%  %ZLIB_SRC%
cmake --build . --config Release --target install

cd %ROOT%
if exist build-openexr-dyn-release del /S /Q build-openexr-dyn-release
mkdir build-openexr-dyn-release
cd build-openexr-dyn-release
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% -DCMAKE_PREFIX_PATH=%INSTALL_REL%^
 -DILMBASE_PACKAGE_PREFIX=%INSTALL_REL% -DOPENEXR_NAMESPACE_VERSIONING=OFF -DOPENEXR_BUILD_STATIC=OFF -DOPENEXR_BUILD_SHARED=ON^
 -DOPENEXR_BUILD_ILMBASE=ON -DOPENEXR_BUILD_OPENEXR=ON -DOPENEXR_BUILD_PYTHON_LIBS=OFF -DOPENEXR_BUILD_UTILS=OFF^
 -DOPENEXR_BUILD_TESTS=OFF -DOPENEXR_BUILD_VIEWERS=OFF  %OPENEXR_SRC%
cmake --build . --config Release --target install

copy OpenEXR\IlmImf\Release\*.dll %INSTALL_REL:/=\%\bin 


cd %ROOT%
if exist build-itk-dyn-release del /S /Q build-itk-dyn-release
mkdir build-itk-dyn-release
cd build-itk-dyn-release
cmake -G %JOMGEN% -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% ^
 -DCMAKE_PREFIX_PATH=%INSTALL_REL% -DITK_USE_SYSTEM_ZLIB=ON -DBUILD_SHARED_LIBS=ON %ITK_SRC%
cmake --build . --config Release --target install


cd %ROOT%
if exist build-astex-dyn-release del /S /Q build-astex-dyn-release
mkdir build-astex-dyn-release
cd build-astex-dyn-release
cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Release" -DCMAKE_PREFIX_PATH=%INSTALL_REL% ^
 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_R%/itkpng -DPNG_LIBRARY_RELEASE=%INSTALL_REL%/lib/itkpng-%ITK_VER%.lib ^
 -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% -DBUILD_SHARED_LIBS=ON %ASTEX_SRC%

copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_R:/=\%\itkpng\

echo You can remove build-zlib-release build-ilmbase-dyn-release build-openexr-dyn-release build-itk-dyn-release

cmake-gui . 

pause
