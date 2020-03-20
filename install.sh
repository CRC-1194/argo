#! /usr/bin/bash

rm -rf build && mkdir build && cd build 
cmake -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=on ..
make -j4 && make install && ctest -V -L UNIT
