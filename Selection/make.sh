#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
	

run "g++ -O3 -c DecayTree.cpp `root-config --cflags` -I. -o DecayTree.o"
#run "g++ `root-config --cflags --libs`  -lRooFit -lRooFitCore -lHtml -lMinuit -c DecayTree.cpp -o DecayTree.o"
#run "g++ `root-config --cflags --libs` -lRooFit -lRooFitCore -c SelectionMaker.cpp -o SelectionMaker.o"

#run "g++ -O3 -o DecayTree.o DecayTree.cpp -ggdb `root-config --cflags --glibs`" 

run "g++ -O3 -c MiniDecayTree.cpp `root-config --cflags` -I. -o MiniDecayTree.o"

#run "g++ -o MiniDecayTree.o MiniDecayTree.cpp `root-config --cflags --glibs` -I. DecayTree.o"


run "g++ -o MiniMaker MiniMaker.cpp `root-config --cflags --glibs` -I. DecayTree.o"
run "g++ -o SelectionMaker SelectionMaker.cpp `root-config --cflags --glibs` -I. DecayTree.o MiniDecayTree.o"

