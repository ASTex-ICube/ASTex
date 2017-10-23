@echo off
call ..\installvars.bat

REM RELEASE 

if exist INSTALL_REL:/=\% (
	copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_R:/=\%\itkpng\
	cd %ROOT%
	del /S /Q build-astex-static-release
	mkdir build-astex-static-release
	cd build-astex-static-release
	cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Release" -DCMAKE_PREFIX_PATH=%INSTALL_REL% ^
	 -DZLIB_INCLUDE_DIR=%INSTALL_REL%/include -DZLIB_LIBRARY_RELEASE=%INSTALL_REL%/lib/zlibstatic.lib ^
	 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_R%/itkpng -DPNG_LIBRARY_RELEASE=%INSTALL_REL%/lib/itkpng-%ITK_VER%.lib ^
	 -DCMAKE_INSTALL_PREFIX=%INSTALL_REL% %ASTEX_SRC%
) else (
	echo ============================================
	echo no %INSTALL_REL:/=\%
	echo ============================================
)


REM DEBUG COMPIL

if exist %ITK_INSTALL_DBG:/=\% (
	copy %ITK_SRC:/=\%\Modules\ThirdParty\PNG\src\itkpng\pnglibconf.h %ITK_INCLUDES_D:/=\%\itkpng\
	cd %ROOT%
	del /S /Q build-astex-static-debug
	mkdir build-astex-static-debug
	cd build-astex-static-debug
	cmake -G %SOLUTION% -DCMAKE_CONFIGURATION_TYPES="Debug" -DCMAKE_PREFIX_PATH=%INSTALL_DBG% ^
	 -DZLIB_INCLUDE_DIR=%INSTALL_DBG%/include -DZLIB_LIBRARY_DEBUG=%INSTALL_DBG%/lib/zlibstatic.lib ^
	 -DPNG_PNG_INCLUDE_DIR=%ITK_INCLUDES_D%/itkpng -DPNG_LIBRARY_DEBUG=%INSTALL_DBG%/lib/itkpng-%ITK_VER%.lib ^
	 -DCMAKE_INSTALL_PREFIX=%INSTALL_DBG% %ASTEX_SRC%
) else (
	echo ============================================
	echo no %INSTALL_DBG:/=\%
	echo ============================================
)


pause
