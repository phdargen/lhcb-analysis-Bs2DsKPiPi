#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
	
run "g++ -O3 -c pull.C `root-config --cflags` -I. -o pull.o"
run "g++ -o SystematicsMaker SystematicsMaker.C `root-config --cflags --glibs` -I. pull.o"
