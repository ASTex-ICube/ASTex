New easy install & compilation of ASTex (and dependencies) on Window

Requirements:
- Visual Studio C++ 2017 / 2015
- CMake (in the path)
- git (in the path)
- [ninja] optionnal faster compilation (used by default if found in path)
- [jom] optionnal faster compilation



Very easy way: double click on installer, answer to questions, and let omputer compile

A launchVS.bat is generated. It allows to launch VS with path included dll dir.

For CMake users, you can directly use CMake


Use this CMakeLists.txt to automatically
 - clone & compile zlib
 - clone & compile openEXR
 - clone & compile itk
 - generated ASTex build


Options:
- ASTEX_INSTALL_PREFIX : install dir for deps & ASTex
- ASTEX_SOURCE: the source dir of ASTex
- ASTEX_BUILD : the build dir of ASTex
- BUILD_SHARED  if checked build dynamic lib else static
- BUILD_RELEASE if checked build release mode else debug
- ADD_POST_TO_DIR if checked add _S for static & _rel or _dbg to created path

Usage
- launch CMake-gui
- where is the source code : this dir
- where to build the binary : you choose (WARNING clean dir between each build)
- Configure
- set vars
- Configure
- Generate
- Finally open the generated solution (in ASTEX_BUILD)



