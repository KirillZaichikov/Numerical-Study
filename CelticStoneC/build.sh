#!/bin/bash
cd build
cmake C:/Users/Kirill/Desktop/repos/Numerical-Study/СelticStoneC -G "MinGW Makefiles"
mingw32-make main
./main.exe