#!/bin/sh
g++ main.cpp -O2 -std=gnu++11 -pthread -lpthread -w -o mvs
#for visualization-enabled version use this:
#g++ main.cpp -O2 -std=gnu++11 -pthread -lpthread -w -o mvs -DUSE_VIS -lX11 
