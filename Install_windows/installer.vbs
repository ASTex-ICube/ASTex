
dim  objShell
Set objShell = wscript.createobject("wscript.shell")


err1 = objShell.Run("cmd.exe /C cmake --version", 0, True)
If err1 <> 0 Then
  Msgbox("install cmake !")
  wscript.Quit
End If

err2 = objShell.Run("cmd.exe /C git --version", 0, True)
If err2 <> 0 Then
  Msgbox("install git !")
  wscript.Quit
End If

err3 = objShell.Run("cmd.exe /C ninja --version", 0, True)
err4 = objShell.Run("cmd.exe /C jom --version", 0, True)
If (err3+err4 > 1) Then
    Msgbox "install ninja or jom", vbOkOnly, "ERROR"
   wscript.Quit
End If

build_deps_dir = replace(fnShellBrowseForFolderVB("Dependencies build dir, tempo (short path!)"),"\","/")
If (Len(build_deps_dir) > 25) Then
    Msgbox "path too long ", vbOkOnly, "ERROR"
   wscript.Quit
End If

inst_deps_dir = replace(fnShellBrowseForFolderVB("Dependencies Install dir"),"\","/")
solu_dir = replace(fnShellBrowseForFolderVB("AStex Solution dir"),"\","/")
script_dir = replace(objShell.CurrentDirectory,"\","/")

objShell.CurrentDirectory = build_deps_dir

dim bshared,brelease

Rep = MsgBox("Shared (Yes) Static (No) ?", vbYesNo + vbQuestion, "Shared / Static")
If Rep = vbYes Then
    bshared = "ON"
Else
    bshared = "OFF"
End If

Rep = MsgBox("Release (Yes) Debug (No) ?", vbYesNo + vbQuestion, "Release / Debug")
If Rep = vbYes Then
    brelease = "ON"
Else
    brelease = "OFF"
End If

dim oReg,vs2015,vs2017,gener
 
Set oReg=GetObject("winmgmts:{impersonationLevel=impersonate}!\\.\root\default:StdRegProv") 
oReg.CheckAccess &H80000000, "VisualStudio.DTE.13.0\CLSID", 1, vs2015
oReg.CheckAccess &H80000000, "VisualStudio.DTE.15.0\CLSID", 1, vs2017

If (Not(vs2017 Or vs2015)) Then
  msqbox("install Visual Studio !")
  wscript.Quit
End If


If (vs2017 And vs2015) Then
    Rep = MsgBox("Visual 2017(Yes) 2015(No) ?", vbYesNo + vbQuestion, "Visual Studio")
    If Rep = vbYes Then
        gener = """Visual Studio 15 2017 Win64"" "
    Else
        gener = """Visual Studio 13 2015 Win64"" "
    End If
else
    If (vs2017) Then 
        gener = """Visual Studio 15 2017 Win64"" "
    Else
        gener = """Visual Studio 13 2015 Win64"" "
    End If
End If

cmd = "cmake.exe -G " + gener + script_dir +_
" -DASTEX_INSTALL_PREFIX:PATH="+inst_deps_dir+_
" -DASTEX_SOURCE:PATH="+script_dir+"/.."+_
" -DASTEX_BUILD:PATH="+solu_dir+_
" -DBUILD_SHARED:BOOL="+bshared+" -DBUILD_RELEASE:BOOL="+brelease+" -DADD_POST_TO_DIR:BOOL=ON"

objShell.run "cmd.exe /C " + cmd , 1 , true
objShell.run "cmd.exe /K  cmake --build ." , 1 , true

dim objFSO
Set objFSO = CreateObject("Scripting.FileSystemObject")
Rep = MsgBox("Remove build dir ?", vbYesNo + vbQuestion, "Clear")
If Rep = vbYes Then
    objFSO.DeleteFolder(build_deps_dir)
End If


function fnShellBrowseForFolderVB(myTitle)
    dim objShell
    dim ssfWINDOWS
    dim objFolder
    set objShell = CreateObject("shell.application")
        set objFolder = objShell.BrowseForFolder(0, myTitle, 0, ssfDESKTOPDIRECTORY)
            If (not objFolder is nothing) Then
                fnShellBrowseForFolderVB = objFolder.Self.Path
            else
                fnShellBrowseForFolderVB = ""
            End If
        set objFolder = nothing
    set objShell = nothing
End function



