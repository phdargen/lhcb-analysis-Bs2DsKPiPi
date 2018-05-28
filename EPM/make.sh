#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
	

#run "g++ -O3 -c mergeTrees.cpp `root-config --cflags` -I. -o mergeTrees.o"
run "g++ -o mergeTrees mergeTrees.cpp `root-config --cflags --glibs` "
