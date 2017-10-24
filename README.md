# ASTex

[ASTex](https://astex-icube.github.io/) is an open-source library for texture analysis and synthesis.
The term “texture” must be understood as in computer graphics.
The purpose of this library is to support collaborative coding and to allow for comparison with recently published algorithms. 

# Installation

## Dependencies:
To compile and use ASTex, you need some libraries:
- ITK 4.10 
- zlib
- openexr for saving/loading images in floating point format.

You need some classic development tools (minimal supported version)
- git
- cmake 3.0
- a _recent_ C++ compiler 
	- g++ 4.9
	- clang 3.3
	- Visual Studio C++ 2015

## Linux
Just install packages:
- libinsighttoolkit4-dev (4.10 min)
- libpng-dev
- libopenexr-dev 

## Mac OS/X
The most simple way to install dependencies is to use [homebrew](https://brew.sh/) package system.
Then you can install the dependencies:
- brew install insighttoolkit
- brew install openexr
- brew install libpng

## Windows

### Softwares:

- VisualStudio C++ (2015 min)
- CMake (3.0 min),
- jom (already installed if you have QtCreator) for multi-threaded compilation of dependencies 

### Sources to download:

- get sources of zlib
- get sources of OpenEXR (use github, archive on OpenEXR.com has bugs)
- get sources of InsightToolKit

### Compile dependencies

- use furnished compilation scripts:
  - make a local copy of folder __script_windows__ to avoid pushing your modifications on git server
  - edit your installvars.bat (check path for exec, src and installation)
  - **WARNING** due to a limitation in Visual-Studio, source and build (of itk) directory path should not be too long (50 char) !
  - double click on install\_static\_debug.bat for compil/install of zlib/openexr/itk and creation of ASTex solution in static debug mode
  - or/and double click on install\_static\_release.bat for compil/install of zlib/openexr/itk and creation of ASTex solution in static release mode
  - or/and double click on install\_dynamic\_debug.bat for compil/install of zlib/openexr/itk and creation of ASTex solution in dynamic debug mode
  - or/adn double click on install\_dynamic\_release.bat for compil/install of zlib/openexr/itk and creation of ASTex solution in dynamic release mode
  - you can reset the ASTex solutions with the reset\_static\_astex\_solu.bat and reset\_dynamic\_astex\_solu.bat scripts
  - after you can customize your build with cmake-gui (using the right build directory)
- __WARNINNG__: if you want to install dynamic and static do not forget to changethe INSTALL\_REL & INSTALL\_DBG vars in installvars.bat
- solutions are in folder build-astex-xxx-release / build-astex-xxx-debug
- all installed libs (including ASTex) are in the same folder:
  - to install ASTex just just generate the INSTALL target
  - installed-Release for release
  - installed-Debug for debug

## Data

Some tests, tutorials and algorithms use read example images and write some results.
In order to keep original Data directory of ASTex clean we use a copy which pass can be choosen 
at cmake configuration stage.

You have to copy yourself the contain of the Data directory in to the right place (see ASTEX\_TEMPO\_PATH clake variable).

## Compilation

### on Linux & Mac
Use CMake as usual:
* create a build directory as same level than ASTex (ASTex-build or ASTex-buildDebug for example)
* go inside build directory and do cmake ../ASTex (or use gui)
* or let (a recent) Qtcreator do the job !

### on Windows + VisualStudio

* build directory has been createde by install script
* use the script _reset\_astex\_solu_ in your copy of _script\_windows_ in case of trouble
* use CMake-gui to customize the build
* then launch Visual and load ASTex solution which has been generated in the build directory

### CMake Options
There are some original options/values that can be set at the cmake stage:

* ASTEX\_ALGO\_xxx choose to build the different implemented algorithms.
* ASTEX\_BUILD\_xxx choose to build bench/tuto/test
* ASTEX\_PERSO\_xxx for each directory added in ASTex that contain a CMakeLists.txt set this to ON to build. When you add a directory just relaunch cmake.
* ASTEX\_TEMPO\_PATH path of directory use to store images for test and tuto (copy ASTex/Data into it)
* ASTEX\_USE\_CPP14 set this to ON if VXL say that you are using a C++ standard version older than the one used ton compile the lib.

