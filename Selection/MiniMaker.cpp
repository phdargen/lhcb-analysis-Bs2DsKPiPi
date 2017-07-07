#include "DecayTree.h"
#include <iostream>

int main(int argc, char** argv){
        
    Decay::Type decay;
    Year::Type year;
    Ds_finalState::Type finalState;
    DataType::Type dataType;
    
    if((string)argv[1] == "Signal") decay = Decay::signal ;
    else if((string)argv[1] == "Norm") decay = Decay::norm ;

    if((string)argv[2] == "Data") dataType = DataType::data ;
    if((string)argv[2] == "MC") dataType = DataType::mc ;
    
    if(atoi(argv[3]) == 0 || atoi(argv[3]) == 4) year = Year::y11;
    else if(atoi(argv[3]) == 1 || atoi(argv[3]) == 5) year = Year::y12;
    else if(atoi(argv[3]) == 2 || atoi(argv[3]) == 6) year = Year::y15;
    else if(atoi(argv[3]) == 3 || atoi(argv[3]) == 7) year = Year::y16;
    
    if(atoi(argv[3]) < 4) finalState = Ds_finalState::phipi;
    else if(atoi(argv[3]) > 4) finalState = Ds_finalState::pipipi;

    DecayTree d(decay, year, finalState, dataType);
    d.Loop();

    return 0;
}
