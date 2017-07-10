#include "MiniDecayTree.h"
#include <iostream>
#include <stdlib.h>
#include <cstdlib>

int main(int argc, char** argv){
        
    Decay::Type decay;
    Year::Type year;
    Ds_finalState::Type finalState;
    DataType::Type dataType;
    
    if((string)argv[1] == "Signal") decay = Decay::signal ;
    else if((string)argv[1] == "Norm") decay = Decay::norm ;
    else {
		cout << "Missing or wrong config! I'll crash now, see commands.txt how to run me " << endl;
		throw "ERROR" ;
	}

    if((string)argv[2] == "Ds2KKpi") finalState = Ds_finalState::phipi; 
    else if((string)argv[2] == "Ds2pipipi") finalState = Ds_finalState::pipipi;
    else if((string)argv[2] == "Ds2Kpipi") finalState = Ds_finalState::Kpipi;
    else {
		cout << "Missing or wrong config! I'll crash now, see commands.txt how to run me " << endl;
		throw "ERROR" ;
	}

    if((string)argv[3] == "Data") dataType = DataType::data ;
    else if((string)argv[3] == "MC") dataType = DataType::mc ;
    else {
		cout << "Missing or wrong config! I'll crash now, see commands.txt how to run me " << endl;
		throw "ERROR" ;
	}

    if(atoi(argv[4]) == 11) year = Year::y11;
    else if(atoi(argv[4]) == 12) year = Year::y12;
    else if(atoi(argv[4]) == 15) year = Year::y15;
    else if(atoi(argv[4]) == 16) year = Year::y16;
    else {
		cout << "Missing or wrong config! I'll crash now, see commands.txt how to run me " << endl;
		throw "ERROR" ;
	}    

    MiniDecayTree d(decay, year, finalState, dataType);
    d.Loop();

    return 0;
}
