#!/bin/bash
cd build
cmake C:\\Users\\zki20\\Desktop\\Numerical-Study\\CelticStoneC -G "MinGW Makefiles"
mingw32-make main
./main.exe