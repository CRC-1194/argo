#! /usr/bin/env bash

BUILD_TYPE=$1

if [ -z $BUILD_TYPE ]
then
    BUILD_TYPE="Release"
fi

rm -rf build && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_EXPORT_COMPILE_COMMANDS=on ..
make -j4 && make install && ctest -V -L UNIT
