Easy install & compilation of ASTex on Window

1 Copy this folder (and work on copy) to avoid unvolontary git commit
  It works whatever the place it is copied

2 Following scripts needs:
  - Zlib sources
  - OpenEXR sources (git version: IlmBase+OpenEXR)
  - ITK sources
  - Visual (2015/2017)
  - CMake
  - jomcompile tool (include in Qt)

3 You can download dependencies sources with dl_deps.bat

4 Check and modify vars in installvars.bat file:
  - Visual vcvarsall & solution type  (for cmake)
  - pathes of cmake/ninja (if not ib global path env)
  - root folder where build folder will be created
  - pathes of sources
  - pathes of installation
  - WARNING: due to a VS limitation pathes (source & build) should not be to long
    itk stop compilation at 50 chars !
 
5 Static script for static-lib installations:
  - release
  - debug

6 Dynamic script for dynamic-lib (dll) installations:
  - release
  - debug

If you want to install dynamic and static do not forget to change
the INSTALL_REL & INSTALL_DBG vars in installvars.bat


