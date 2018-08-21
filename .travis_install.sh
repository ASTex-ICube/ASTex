#!/bin/bash
if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
 sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
 sudo apt-get update
 sudo apt-get install "g++-5"
 sudo apt-get install cmake3
 sudo apt-get install zlib1g-dev
 sudo apt-get install libpng12-dev
 sudo apt-get install libopenexr-dev 
 pushd /tmp/
 wget https://mycore.core-cloud.net/index.php/s/BfSTlzfQIDqvI8l/download -O itk.tgz
 cd /usr/
 sudo tar xf /tmp/itk.tgz
 popd
fi

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
 pushd /tmp/
 wget https://mycore.core-cloud.net/index.php/s/CSKCKkV1uPAWoJT/download -O itk_osx.tgz
 cd /usr/local
 sudo tar xf /tmp/itk_osx.tgz
 popd
fi
 
