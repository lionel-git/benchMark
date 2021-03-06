#!/bin/bash
#set WINE_PATH: export WINEPATH="/usr/x86_64-w64-mingw32/sys-root/mingw/bin/"
x86_64-w64-mingw32-g++ -O3 -static-libgcc -static-libstdc++ benchMark.cpp -o test_mingw64_static
x86_64-w64-mingw32-g++ -O3  benchMark.cpp -o test_mingw64_dynamic
