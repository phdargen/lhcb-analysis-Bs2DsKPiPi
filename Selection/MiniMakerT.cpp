#include "DecayTree.h"
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include "TThread.h"
//#include "TMethodCall.h"

DecayTree* d;

void *handle1(void *) {

    TThread::Lock();
    Decay::Type decay;
    Year::Type year;
    Ds_finalState::Type finalState;
    DataType::Type dataType;
    decay = Decay::signal ;
    dataType = DataType::data ;
    year = Year::y11;
    finalState = Ds_finalState::phipi;
    
    d = new DecayTree(decay, year, finalState, dataType);
    d->Init();  
    TThread::UnLock();

    d->Loop();

return (0);
}

void *handle2(void *) {

    return (0);

    TThread::Lock();

    Decay::Type decay;
    Year::Type year;
    Ds_finalState::Type finalState;
    DataType::Type dataType;
    
    
    decay = Decay::signal ;
    dataType = DataType::data ;
    year = Year::y12;
    finalState = Ds_finalState::phipi;

    
    //DecayTree d(decay, year, finalState, dataType);
    
    TThread::UnLock();

    //d.Loop();
    return (0);

}


int main(int argc, char** argv){
        
    Decay::Type decay;
    Year::Type year;
    Ds_finalState::Type finalState;
    DataType::Type dataType;
    /*    
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
     */
  
    //d->Loop();
    
    //return 0;
        

    printf("Starting Thread 0\n");
    TThread* thread1 = new TThread("t0", handle1, (void*) 0);
    thread1->Run();
    printf("Starting Thread 1\n");
    TThread*thread2 = new TThread("t1", handle2, (void*) 1);
    thread2->Run();
    //printf("Starting Joiner Thread \n");
    //threadj = new TThread("t4", joiner, (void*) 3);
    //threadj->Run();
    
    printf("Ps\n");
    TThread::Ps();


    return 0;
}
