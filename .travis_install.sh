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
 wget http://igg.unistra.fr/people/thery/itk.tgz
 cd /usr/
 sudo tar xf /tmp/itk.tgz
 popd
fi

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
 brew install insighttoolkit
 brew install openexr
# brew install libpng
fi
 
