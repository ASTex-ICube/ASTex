
:: customize if necessary
set GIT="C:\Program Files\Git\bin\git.exe"

%GIT% clone https://github.com/madler/zlib.git -b v1.2.11 --depth 1
%GIT% clone https://github.com/openexr/openexr.git -b v2.3.0 --depth 1
%GIT% clone https://github.com/InsightSoftwareConsortium/ITK.git -b v4.13.1 --depth 1
