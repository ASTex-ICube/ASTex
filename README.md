# ASTex

[ASTex](https://astex-icube.github.io/) is an open-source library for texture analysis and synthesis.
The term “texture” must be understood as in computer graphics.
The purpose of this library is to support collaborative coding and to allow for comparison with recently published algorithms. 

# Installation

## Dependencies:
To compile and use ASTex, you need some libraries:
- ITK 5. min 
- zlib
- openexr for saving/loading images in floating point format.

You need some classic development tools (minimal supported version)
- git
- cmake 3.0 min
- a _recent_ C++ compiler 
	- g++ 6 
	- clang 3.3
	- Visual Studio C++ 2017 /2019

## Linux
Just install packages:
- libpng-dev
- libopenexr-dev 
And get github version of InsightToolkit, compile & install

## Mac OS/X
The most simple way to install dependencies is to use [homebrew](https://brew.sh/) package system.
Then you can install the dependencies:
- brew install insighttoolkit
- brew install openexr
- brew install libpng

## Windows

### Softwares:

- VisualStudio C++ (2017 min)
- Git
- CMake (3.0 min)

### install VCPKG For ASTex
- open a power-shell
- git clone https://github.com/Microsoft/vcpkg.git 
- var env VCPKG_DEFAULT_TRIPLET=x64-windows
- cd  vcpkg
- .\vcpkg\bootstrap-vcpkg.bat
- .\vcpkg.exe install itk openexr

- Remarque: ATTENTION LA COMPILATION NECESSITE X.X Go D'ESPACE DISQUE ! Vous pouvez enlever les répertoires buildtrees, downloads et packages

### Compile ASTex
- launch CMake, chose src dir and binary dir
- Specify tool chain file:  `XXXX/vcpkg_ast/scripts/buildsystems/vcpkg.cmake`
- Configure again
- Generate
- Launch Solution 
- Enjoy

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


### CMake Options
There are some original options/values that can be set at the cmake stage:

* ASTEX\_ALGO\_xxx choose to build the different implemented algorithms.
* ASTEX\_BUILD\_xxx choose to build bench/tuto/test
* ASTEX\_PERSO\_xxx for each directory added in ASTex that contain a CMakeLists.txt set this to ON to build. When you add a directory just relaunch cmake.
* ASTEX\_TEMPO\_PATH path of directory use to store images for test and tuto (copy ASTex/Data into it)

# Contributing to ASTex

All contributions to ASTex are welcome and will be examined.

To contribute:
* you need a github account.
* fork ASTex repository on your account
* create a branch from the develop one
* develop your contribution 
* do pull-request on ASTex/develop branch
