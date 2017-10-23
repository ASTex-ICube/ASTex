@echo off
call ..\installvars.bat

REM RELEASE 

if exist INSTALL_REL:/=\% (
	copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_R:/=\%\itkpng\
	cd %ROOT%
	if exist build-astex-dyn-release del /S /Q build-astex-dyn-release
	mkdir build-astex-dyn-release
	cd build-astex-dyn-release
	cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Release" -DCMAKE_PREFIX_PATH=%INSTALL_REL% ^
	 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_R%/itkpng -DPNG_LIBRARY_RELEASE=%INSTALL_REL%/lib/itkpng-%ITK_VER%.lib ^
	 -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% -DBUILD_SHARED_LIBS=ON %ASTEX_SRC%
) else (
	echo ============================================
	echo no %INSTALL_REL:/=\%
	echo ============================================
)


REM DEBUG COMPIL

if exist %ITK_INSTALL_DBG:/=\% (
	copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_D:/=\%\itkpng\
	cd %ROOT%
	if exist build-astex-dyn-debug del /S /Q build-astex-dyn-debug
	mkdir build-astex-dyn-debug
	cd build-astex-dyn-debug
	cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Debug" -DCMAKE_PREFIX_PATH=%INSTALL_DBG% ^
	 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_D%/itkpng -DPNG_LIBRARY_DEBUG=%INSTALL_DBG%/lib/itkpng-%ITK_VER%.lib ^
	 -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% -DBUILD_SHARED_LIBS=ON %ASTEX_SRC%
) else (
	echo ============================================
	echo no %INSTALL_DBG:/=\%
	echo ============================================
)


pause
