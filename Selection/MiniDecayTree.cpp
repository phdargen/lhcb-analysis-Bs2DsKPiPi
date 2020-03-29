#define MiniDecayTree_cxx
#include "MiniDecayTree.h"
#include "DecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "TLorentzVector.h"
#include "TRandom3.h"

using namespace std; 


inline Bool_t MiniDecayTree::PhaseSpace_Cuts(){

    if(_decay == Decay::signal){
        if((BsDTF_K_plus + BsDTF_pi_plus + BsDTF_pi_minus).M() > 1950.) return false;
        if((BsDTF_K_plus + BsDTF_pi_minus).M()  > 1200.) return false;
        if((BsDTF_pi_plus + BsDTF_pi_minus).M() > 1200.) return false;

        if((BsDTF_K_plus + BsDTF_pi_plus + BsDTF_pi_minus).M() < massKaon + 2. * massPion) return false;
        if((BsDTF_K_plus + BsDTF_pi_minus).M()  < massKaon + massPion) return false;
        if((BsDTF_pi_plus + BsDTF_pi_minus).M() < 2. * massPion) return false;
    }
    else {
    	if((BsDTF_pi_plus1 + BsDTF_pi_plus2 + BsDTF_pi_minus).M() > 1950.) return false;
        if((BsDTF_pi_plus1 + BsDTF_pi_minus).M()  > 1200.) return false;
        if((BsDTF_pi_plus2 + BsDTF_pi_minus).M()  > 1200.) return false;

    	if((BsDTF_pi_plus1 + BsDTF_pi_plus2 + BsDTF_pi_minus).M() < 3. * massPion) return false;
        if((BsDTF_pi_plus1 + BsDTF_pi_minus).M()  < 2. * massPion) return false;
        if((BsDTF_pi_plus2 + BsDTF_pi_minus).M()  < 2. * massPion) return false;
    }
    return true;
}

inline Bool_t MiniDecayTree::PIDCalib_Cuts(){
	return true;
}

inline Bool_t MiniDecayTree::MC_Cuts(){

	if(Bs_BKGCAT != 20 && Bs_BKGCAT != 60)return false;
	if(abs(Ds_TRUEID)!= 431){
			vector<int> mother_ids;
			int count1,count2;
			if(_Ds_finalState == Ds_finalState::pipipi){
        			mother_ids.push_back(pi_plus_fromDs_MC_MOTHER_ID);
		        	mother_ids.push_back(pi_minus_fromDs_MC_MOTHER_ID);
        			mother_ids.push_back(pi_minus2_fromDs_MC_MOTHER_ID);
			}
			else if(_Ds_finalState == Ds_finalState::Kpipi){
        			mother_ids.push_back(K_minus_fromDs_MC_MOTHER_ID);
		        	mother_ids.push_back(pi_plus_fromDs_MC_MOTHER_ID);
        			mother_ids.push_back(pi_minus_fromDs_MC_MOTHER_ID);
			}
  			else {
        			mother_ids.push_back(K_plus_fromDs_MC_MOTHER_ID);
		        	mother_ids.push_back(K_minus_fromDs_MC_MOTHER_ID);
        			mother_ids.push_back(pi_minus_fromDs_MC_MOTHER_ID);
			}
			count1 = std::count(mother_ids.begin(),mother_ids.end(),431);
			count2 = std::count(mother_ids.begin(),mother_ids.end(),-431);
			if(count1 >= 2)Ds_TRUEID = 431;
			else if(count2 >= 2)Ds_TRUEID = -431;
			else Ds_TRUEID = 0;
	}	
	if(abs(Bs_TRUEID)!= 531){	
		vector<int> mother_ids;
		int count1,count2;
		if(_decay == Decay::signal){
			mother_ids.push_back(K_plus_MC_MOTHER_ID);
        		mother_ids.push_back(pi_plus_MC_MOTHER_ID);
        		mother_ids.push_back(pi_minus_MC_MOTHER_ID);
        		mother_ids.push_back(Ds_MC_MOTHER_ID);
		}
		else {
			mother_ids.push_back(pi_plus1_MC_MOTHER_ID);
        		mother_ids.push_back(pi_plus2_MC_MOTHER_ID);
        		mother_ids.push_back(pi_minus_MC_MOTHER_ID);
        		mother_ids.push_back(Ds_MC_MOTHER_ID);
		}
		count1 = std::count(mother_ids.begin(),mother_ids.end(),531);
		count2 = std::count(mother_ids.begin(),mother_ids.end(),-531);
		if(count1>count2 && count1 >= 2)Bs_TRUEID = 531;
		else if(count2 > count1 && count2 >= 2)Bs_TRUEID = -531;
		else Bs_TRUEID = 0;
	}
	/*
	if(abs(Bs_TRUEID)!= 531)return false;
    	if(_decay == Decay::signal){
		if(abs(K_plus_TRUEID)!= 321)return false;
		if(abs(pi_plus_TRUEID)!= 211)return false;
		if(abs(pi_minus_TRUEID)!= 211)return false;
	}
    	else{
		if(abs(pi_plus1_TRUEID)!= 211)return false;
		if(abs(pi_plus2_TRUEID)!= 211)return false;
		if(abs(pi_minus_TRUEID)!= 211)return false;
	}
        if(_Ds_finalState == Ds_finalState::pipipi){
        	if(abs(pi_plus_fromDs_TRUEID) != 211) return false; 
    		if(abs(pi_minus_fromDs_TRUEID) != 211) return false;
		if(abs(pi_minus2_fromDs_TRUEID) != 211) return false;
	}
  	else if(_Ds_finalState == Ds_finalState::Kpipi){
    	}
  	else {
        	if(abs(K_plus_fromDs_TRUEID) != 321) return false; 
    		if(abs(K_minus_fromDs_TRUEID) != 321) return false;
		if(abs(pi_minus_fromDs_TRUEID) != 211) return false;
	}	
	*/
	return true;
}

inline Bool_t MiniDecayTree::PID_Cuts(){

   if(_Ds_finalState == Ds_finalState::phipi || _Ds_finalState == Ds_finalState::KsK || _Ds_finalState == Ds_finalState::KKpi_NR){
	// remove events with no PID info
	if( K_plus_fromDs_hasRich == 0 ) return false;        
	if( K_minus_fromDs_hasRich == 0 ) return false;        
	if( pi_minus_fromDs_hasRich == 0 ) return false; 
	if( K_plus_fromDs_PIDK == 0 ) return false;        
	if( K_minus_fromDs_PIDK == 0 ) return false;        
	//if( pi_minus_fromDs_PIDK == 0 ) return false;         
    }

    if(_Ds_finalState == Ds_finalState::phipi){
        if(K_plus_fromDs_PIDK < -10) return false;
        else if(K_minus_fromDs_PIDK < -10) return false;
        else if(pi_minus_fromDs_PIDK > 20) return false;
    }

    else if(_Ds_finalState == Ds_finalState::KsK){
        if(K_plus_fromDs_PIDK < -10) return false;
        else if(K_minus_fromDs_PIDK < -5) return false;
        else if(pi_minus_fromDs_PIDK > 10) return false;        
    }

    else if(_Ds_finalState == Ds_finalState::KKpi_NR){
        if(K_plus_fromDs_PIDK < 5) return false;
        else if(K_minus_fromDs_PIDK < 5) return false;
        else if(pi_minus_fromDs_PIDK > 10) return false;
        //else if(pi_minus_fromDs_PIDp > 20) return false;
    }
    
    else if(_Ds_finalState == Ds_finalState::pipipi){
        if(pi_plus_fromDs_PIDK > 10) return false;
        else if(pi_minus_fromDs_PIDK > 10) return false;
        else if(pi_minus2_fromDs_PIDK > 10) return false;
        
        if(pi_plus_fromDs_PIDp > 20) return false;
        else if(pi_minus_fromDs_PIDp > 20) return false;
        else if(pi_minus2_fromDs_PIDp > 20) return false;
	// remove events with no PID info
	if( pi_plus_fromDs_hasRich == 0  ) return false;        
	if( pi_minus_fromDs_hasRich == 0  ) return false;        
	if( pi_minus2_fromDs_hasRich == 0  ) return false;  
	//if( pi_plus_fromDs_PIDK == 0  ) return false;        
	//if( pi_minus_fromDs_PIDK == 0  ) return false;        
	//if( pi_minus2_fromDs_PIDK == 0  ) return false;
    }
    
    else if(_Ds_finalState == Ds_finalState::Kpipi){
        if(pi_plus_fromDs_PIDK > 5) return false;
        else if(K_minus_fromDs_PIDK < 8) return false;
        else if(pi_minus_fromDs_PIDK > 5) return false;

        if(pi_plus_fromDs_PIDp > 20) return false;
        else if(pi_minus_fromDs_PIDp > 20) return false;
	// remove events with no PID info
	if( pi_plus_fromDs_hasRich == 0  ) return false;        
	if( K_minus_fromDs_hasRich == 0  ) return false;        
	if( pi_minus_fromDs_hasRich == 0  ) return false;     
	//if( pi_plus_fromDs_PIDK == 0  ) return false;        
	if( K_minus_fromDs_PIDK == 0  ) return false;        
	//if( pi_minus_fromDs_PIDK == 0  ) return false;                     
    }
    

    if(_decay == Decay::signal){
        if(K_plus_PIDK < 10) return false;
        else if(pi_plus_PIDK > 10) return false;
        else if(pi_minus_PIDK > 0) return false;
	// remove events with no PID info
	if( K_plus_hasRich == 0  ) return false;        
	if( pi_plus_hasRich == 0  ) return false;        
	if( pi_minus_hasRich == 0  ) return false;
	//if( K_plus_PIDK == 0  ) return false;        
	//if( pi_plus_PIDK == 0  ) return false;        
	//if( pi_minus_PIDK == 0  ) return false;        
    }

    else if(_decay == Decay::norm){
        if(pi_plus1_PIDK > 0) return false;
        else if(pi_plus2_PIDK > 0) return false;
        else if(pi_minus_PIDK > 10) return false;
	// remove events with no PID info
	if( pi_plus1_hasRich == 0  ) return false;        
	if( pi_plus2_hasRich == 0  ) return false;        
	if( pi_minus_hasRich == 0  ) return false;     
	//if( pi_plus1_PIDK == 0  ) return false;        
	//if( pi_plus2_PIDK == 0  ) return false;        
	//if( pi_minus_PIDK == 0  ) return false;     
    }

    return true;
}

inline Bool_t MiniDecayTree::Veto_Cuts(){

    if(_Ds_finalState == Ds_finalState::phipi || _Ds_finalState == Ds_finalState::KsK || _Ds_finalState == Ds_finalState::KKpi_NR){
        
        //D0 veto
        //if((K_plus_fromDs + K_minus_fromDs + pi_minus_fromDs).M() - (K_plus_fromDs + K_minus_fromDs).M() < 155.) return false;
        if((K_plus_fromDs + K_minus_fromDs).M() > 1840.) return false;
        
        if(_Ds_finalState == Ds_finalState::phipi){
            // Charmless veto
            if(!_ltu)if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < -1) return false;
            if(Ds_FDCHI2_ORIVX < 0) return false;
            //Lambda_c veto
            if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 40. && ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 5) ) return false;
	    //D- veto
            if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 40. && K_minus_fromDs_PIDK < 2 ) return false;
	}
        else if(_Ds_finalState == Ds_finalState::KsK){
            // Charmless veto
	    if(!_ltu)if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) return false;
            if(Ds_FDCHI2_ORIVX < 0) return false;
            //Lambda_c veto
            if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 40. && ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 5) ) return false;
	    //D- veto
            if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 40. && K_minus_fromDs_PIDK < 15 ) return false;
        }
        else if(_Ds_finalState == Ds_finalState::KKpi_NR){
	    // Charmless veto 
            if(!_ltu)if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) return false;
            if(Ds_FDCHI2_ORIVX < 4) return false;
	    //Lambda_c veto
            if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 40. && ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 5) ) return false;
	    //D- veto
            if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 40. && K_minus_fromDs_PIDK < 15 ) return false;
        }
    }

    else if(_Ds_finalState == Ds_finalState::pipipi){
        // D0 veto
        //if( ((pi_plus_fromDs + pi_minus_fromDs + pi_minus2_fromDs).M() - (pi_plus_fromDs + pi_minus_fromDs).M()) < 155. 
	//||  ((pi_plus_fromDs + pi_minus_fromDs + pi_minus2_fromDs).M() - (pi_plus_fromDs + pi_minus2_fromDs).M()) < 155.) return false;
        if( (pi_plus_fromDs + pi_minus_fromDs).M() > 1700. || (pi_plus_fromDs + pi_minus2_fromDs).M() > 1700.) return false;
        // Charmless veto
        if(!_ltu)if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) return false;
        if(Ds_FDCHI2_ORIVX < 6) return false;
        //Lambda_c veto
        //if( TMath::Abs((pi_plus_fromDs + piminus_fromDs_asProton_MissID + pi_minus2_fromDs).M() - massLambda_c) < 30. && pi_minus_fromDs_PIDp < 0 ) return false;
        //if( TMath::Abs((pi_plus_fromDs + piminus2_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 30. && pi_minus2_fromDs_PIDp < 0 ) return false;
    }
    
    else if(_Ds_finalState == Ds_finalState::Kpipi){
	// D0 veto
        //if((pi_plus_fromDs + K_minus_fromDs + pi_minus_fromDs).M() - (pi_plus_fromDs + K_minus_fromDs).M() < 155.) return false;
        if((pi_plus_fromDs + K_minus_fromDs).M() > 1750.) return false;
        // Charmless veto
        if(!_ltu)if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) return false;
        if(Ds_FDCHI2_ORIVX < 6) return false;
	//Lambda_c veto
        if( TMath::Abs((pi_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 40. && ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 5) ) return false;
	//D- veto
        if( TMath::Abs((pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 40. && K_minus_fromDs_PIDK < 15 ) return false;
    }
    
    if(_decay == Decay::signal){
        // Ds veto
        //if(TMath::Abs((BsDTF_K_plus + BsDTF_pi_plus + BsDTF_pi_minus).M() - massDs) < 20) return false;
        //if(TMath::Abs((K_plus + pi_plus + pi_minus_asK_MissID).M() - massDs) < 20 && pi_minus_PIDK > -5) return false;
	// Semi-lep. veto
	//if( K_plus_isMuon == 1 ) return false;
    }
    else if(_decay == Decay::norm){
        // Ds veto
        if(TMath::Abs((pi_plus1 + pi_plus2 + pi_minus).M() - massDs) < 20) return false;
        if(TMath::Abs((pi_plus1_asK_MissID + pi_plus2 + pi_minus).M() - massDs) < 20 && pi_plus1_PIDK > -5) return false;
        if(TMath::Abs((pi_plus2_asK_MissID + pi_plus1 + pi_minus).M() - massDs) < 20 && pi_plus2_PIDK > -5) return false;
	// Semi-lep. veto
	//if( pi_plus1_isMuon == 1 ) return false;
	//if( pi_plus2_isMuon == 1 ) return false;
    }

    // Wrong PV veto
    if(!_ltu)if(nPV > 1 && (Bs_MINIPCHI2NEXTBEST-Bs_IPCHI2_OWNPV) < 20) return false;

    return true;
}

inline Bool_t MiniDecayTree::Preselection_Cuts(){
    
	//if(Ds_PT < 1600) return false;
	if(_charmLess){
		if(fabs(Bs_PV_Dplus_M[0]-massDs) < 40 ) return false;
	}
	else {
        	if(_bkg)if(fabs(Ds_MM-massDs) > 25 ) return false;    
		if(!_bkg)if(fabs(Bs_PV_Dplus_M[0]-massDs) > 25 ) return false;
	}
        return true;
}

inline Bool_t MiniDecayTree::LTU_Cuts(){
    
	//if((Bs_ENDVERTEX_CHI2/Bs_ENDVERTEX_NDOF)> 5) return false;
        //if(Bs_IPCHI2_OWNPV>20) return false;
    	if((Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF)> 5) return false;
        //if(Ds_FDCHI2_ORIVX < 25) return false;
    	//if(Ds_FDCHI2_OWNPV < 25) return false;
        if(Ds_IPCHI2_OWNPV>9) return false;
    	//if(Ds_PT < 1800) return false;
        if(Ds_DIRA_OWNPV<0.99994) return false;

	if(_decay == Decay::signal){	
		//if(K_plus_PIDK < 15) return false;
		//if(pi_plus_PIDK > 0) return false;
		//else if(pi_minus_PIDK > 0) return false;
		//if(K_plus_isMuon == 1) return false;
		//if(pi_plus_isMuon == 1) return false;
		//if(pi_minus_isMuon == 1) return false;
	}
	
        return true;
}

inline Ds_finalState::Type MiniDecayTree::get_Ds_finalState(){

    if(_Ds_finalState == Ds_finalState::pipipi || _Ds_finalState == Ds_finalState::Kpipi) return _Ds_finalState;

    else {
    
        if( TMath::Abs(((DTF_K_plus_fromDs + DTF_K_minus_fromDs).M() - massPhi)) < 12) return Ds_finalState::phipi;
        
        else if (TMath::Abs(((DTF_pi_minus_fromDs + DTF_K_plus_fromDs).M() - massKstar)) < 75) return Ds_finalState::KsK;
        
        else return Ds_finalState::KKpi_NR;
    
    }
}

inline void MiniDecayTree::set_LorentzVectors(){
    
    if(_Ds_finalState == Ds_finalState::phipi || _Ds_finalState == Ds_finalState::KsK || _Ds_finalState == Ds_finalState::KKpi_NR){
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);        
        K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
        K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        Ds = pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs;
        
        DTF_pi_minus_fromDs=TLorentzVector(Bs_DTF_D_splus_piplus_PX[0],Bs_DTF_D_splus_piplus_PY[0],Bs_DTF_D_splus_piplus_PZ[0],Bs_DTF_D_splus_piplus_PE[0]  );  
        DTF_K_minus_fromDs=TLorentzVector(Bs_DTF_D_splus_Kplus_PX[0],Bs_DTF_D_splus_Kplus_PY[0],Bs_DTF_D_splus_Kplus_PZ[0],Bs_DTF_D_splus_Kplus_PE[0]);
        DTF_K_plus_fromDs=TLorentzVector(Bs_DTF_D_splus_Kplus_0_PX[0],Bs_DTF_D_splus_Kplus_0_PY[0],Bs_DTF_D_splus_Kplus_0_PZ[0],Bs_DTF_D_splus_Kplus_0_PE[0]);
	DTF_Ds = DTF_pi_minus_fromDs + DTF_K_plus_fromDs + DTF_K_minus_fromDs;

        BsDTF_pi_minus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_piplus_PX[0],Bs_BsDTF_D_splus_piplus_PY[0],Bs_BsDTF_D_splus_piplus_PZ[0],Bs_BsDTF_D_splus_piplus_PE[0]  );  
        BsDTF_K_minus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_Kplus_PX[0],Bs_BsDTF_D_splus_Kplus_PY[0],Bs_BsDTF_D_splus_Kplus_PZ[0],Bs_BsDTF_D_splus_Kplus_PE[0]);
        BsDTF_K_plus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_Kplus_0_PX[0],Bs_BsDTF_D_splus_Kplus_0_PY[0],Bs_BsDTF_D_splus_Kplus_0_PZ[0],Bs_BsDTF_D_splus_Kplus_0_PE[0]);
	BsDTF_Ds = BsDTF_pi_minus_fromDs + BsDTF_K_plus_fromDs + BsDTF_K_minus_fromDs;
        

        Kminus_fromDs_asProton_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);
        Kminus_fromDs_asPiminus_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);
    }
    
    else if(_Ds_finalState == Ds_finalState::pipipi){
        pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);        
        pi_minus2_fromDs.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massPion);
        Ds = pi_minus_fromDs + pi_plus_fromDs + pi_minus2_fromDs;

        DTF_pi_minus_fromDs=TLorentzVector(Bs_DTF_D_splus_piplus_PX[0],Bs_DTF_D_splus_piplus_PY[0],Bs_DTF_D_splus_piplus_PZ[0],Bs_DTF_D_splus_piplus_PE[0]  );  
        DTF_pi_plus_fromDs=TLorentzVector(Bs_DTF_D_splus_piplus_1_PX[0],Bs_DTF_D_splus_piplus_1_PY[0],Bs_DTF_D_splus_piplus_1_PZ[0],Bs_DTF_D_splus_piplus_1_PE[0]);
        DTF_pi_minus2_fromDs=TLorentzVector(Bs_DTF_D_splus_piplus_0_PX[0],Bs_DTF_D_splus_piplus_0_PY[0],Bs_DTF_D_splus_piplus_0_PZ[0],Bs_DTF_D_splus_piplus_0_PE[0]);
	DTF_Ds = DTF_pi_minus_fromDs + DTF_pi_plus_fromDs + DTF_pi_minus2_fromDs;

        BsDTF_pi_minus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_piplus_PX[0],Bs_BsDTF_D_splus_piplus_PY[0],Bs_BsDTF_D_splus_piplus_PZ[0],Bs_BsDTF_D_splus_piplus_PE[0]  );  
        BsDTF_pi_plus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_piplus_1_PX[0],Bs_BsDTF_D_splus_piplus_1_PY[0],Bs_BsDTF_D_splus_piplus_1_PZ[0],Bs_BsDTF_D_splus_piplus_1_PE[0]);
        BsDTF_pi_minus2_fromDs=TLorentzVector(Bs_BsDTF_D_splus_piplus_0_PX[0],Bs_BsDTF_D_splus_piplus_0_PY[0],Bs_BsDTF_D_splus_piplus_0_PZ[0],Bs_BsDTF_D_splus_piplus_0_PE[0]);
	BsDTF_Ds = BsDTF_pi_minus_fromDs + BsDTF_pi_plus_fromDs + BsDTF_pi_minus2_fromDs;
        
        piminus_fromDs_asProton_MissID = TLorentzVector(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massProton);
        piminus2_fromDs_asProton_MissID = TLorentzVector(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massProton);
    }
    
    else if(_Ds_finalState == Ds_finalState::Kpipi){
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);        
        pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
        K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        Ds = pi_minus_fromDs + pi_plus_fromDs + K_minus_fromDs;

        DTF_pi_minus_fromDs=TLorentzVector(Bs_DTF_D_splus_piplus_PX[0],Bs_DTF_D_splus_piplus_PY[0],Bs_DTF_D_splus_piplus_PZ[0],Bs_DTF_D_splus_piplus_PE[0]  );  
        DTF_pi_plus_fromDs=TLorentzVector(Bs_DTF_D_splus_piplus_0_PX[0],Bs_DTF_D_splus_piplus_0_PY[0],Bs_DTF_D_splus_piplus_0_PZ[0],Bs_DTF_D_splus_piplus_0_PE[0]);
        DTF_K_minus_fromDs=TLorentzVector(Bs_DTF_D_splus_Kplus_PX[0],Bs_DTF_D_splus_Kplus_PY[0],Bs_DTF_D_splus_Kplus_PZ[0],Bs_DTF_D_splus_Kplus_PE[0]);
	DTF_Ds = DTF_pi_minus_fromDs + DTF_pi_plus_fromDs + DTF_K_minus_fromDs;        

        BsDTF_pi_minus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_piplus_PX[0],Bs_BsDTF_D_splus_piplus_PY[0],Bs_BsDTF_D_splus_piplus_PZ[0],Bs_BsDTF_D_splus_piplus_PE[0]  );  
        BsDTF_pi_plus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_piplus_0_PX[0],Bs_BsDTF_D_splus_piplus_0_PY[0],Bs_BsDTF_D_splus_piplus_0_PZ[0],Bs_BsDTF_D_splus_piplus_0_PE[0]);
        BsDTF_K_minus_fromDs=TLorentzVector(Bs_BsDTF_D_splus_Kplus_PX[0],Bs_BsDTF_D_splus_Kplus_PY[0],Bs_BsDTF_D_splus_Kplus_PZ[0],Bs_BsDTF_D_splus_Kplus_PE[0]);
	BsDTF_Ds = BsDTF_pi_minus_fromDs + BsDTF_pi_plus_fromDs + BsDTF_K_minus_fromDs;        

        Kminus_fromDs_asProton_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);
        Kminus_fromDs_asPiminus_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);
    }
    
    if(_decay == Decay::signal){
        K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);        
    
        BsDTF_K_plus= TLorentzVector(Bs_BsDTF_K_1_1270_plus_Kplus_PX[0],Bs_BsDTF_K_1_1270_plus_Kplus_PY[0],Bs_BsDTF_K_1_1270_plus_Kplus_PZ[0],Bs_BsDTF_K_1_1270_plus_Kplus_PE[0]) ;
        BsDTF_pi_plus= TLorentzVector(Bs_BsDTF_K_1_1270_plus_piplus_0_PX[0], Bs_BsDTF_K_1_1270_plus_piplus_0_PY[0], Bs_BsDTF_K_1_1270_plus_piplus_0_PZ[0], Bs_BsDTF_K_1_1270_plus_piplus_0_PE[0]);
        BsDTF_pi_minus= TLorentzVector( Bs_BsDTF_K_1_1270_plus_piplus_PX[0], Bs_BsDTF_K_1_1270_plus_piplus_PY[0], Bs_BsDTF_K_1_1270_plus_piplus_PZ[0], Bs_BsDTF_K_1_1270_plus_piplus_PE[0]) ;  

        pi_minus_asK_MissID.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ, massKaon);
	K_plus_asMu_MissID.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massMuon);
        pi_plus_asMu_MissID.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massMuon);        
    }
    
    else if(_decay == Decay::norm){        
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);

        BsDTF_pi_plus1 = TLorentzVector(Bs_BsDTF_a_1_1260_plus_piplus_1_PX[0],Bs_BsDTF_a_1_1260_plus_piplus_1_PY[0],Bs_BsDTF_a_1_1260_plus_piplus_1_PZ[0], Bs_BsDTF_a_1_1260_plus_piplus_1_PE[0]) ;
        BsDTF_pi_plus2= TLorentzVector( Bs_BsDTF_a_1_1260_plus_piplus_PX[0],Bs_BsDTF_a_1_1260_plus_piplus_PY[0],Bs_BsDTF_a_1_1260_plus_piplus_PZ[0],Bs_BsDTF_a_1_1260_plus_piplus_PE[0]) ;
        BsDTF_pi_minus= TLorentzVector(Bs_BsDTF_a_1_1260_plus_piplus_0_PX[0],Bs_BsDTF_a_1_1260_plus_piplus_0_PY[0], Bs_BsDTF_a_1_1260_plus_piplus_0_PZ[0],Bs_BsDTF_a_1_1260_plus_piplus_0_PE[0]) ;

        pi_minus_asK_MissID.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ, massKaon);
        pi_plus1_asK_MissID.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ, massKaon);
        pi_plus2_asK_MissID.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ, massKaon);

        pi_plus1_asMu_MissID.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ, massMuon);
        pi_plus2_asMu_MissID.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ, massMuon);
    }
    
    return;
}


TTree* MiniDecayTree::GetInputTree(){
        
    cout << "Read from file " << _inFileName << endl;    
        
    TChain* chain= new TChain("DecayTree");
    chain->Add(_inFileName);
    //chain->Add("Signal/11D/*.root");
    if(chain==0){
        cout << "ERROR: No file found" << endl;
        throw "ERROR";
    }
    
    return (TTree*)chain;
}

void MiniDecayTree::Loop()
{

   time_t startTime = time(0);

   Init();
   if (fChain == 0) return;

   TFile* output = new TFile(_outFileName,"RECREATE");
   TTree* summary_tree = fChain->CloneTree(0);
    
   Long64_t nentries = fChain->GetEntries();
   cout << "Have " << nentries << " events" <<  endl;
   
    // Add new branches to output tree
    int Ds_finalState,Ds_finalState_mod,Ds_finalState_mod2, sample, TriggerCat,run,Bs_BKGCAT = 0;
    int year = (int)_year;
    if(year == 11 || year == 12) run = 1;
    else run = 2;
    double NTracks;

    summary_tree->Branch("Ds_finalState",&Ds_finalState,"Ds_finalState/I");    
    summary_tree->Branch("Ds_finalState_mod",&Ds_finalState_mod,"Ds_finalState_mod/I");    
    summary_tree->Branch("Ds_finalState_mod2",&Ds_finalState_mod2,"Ds_finalState_mod2/I");    
    summary_tree->Branch("sample", &sample, "sample/I");
    summary_tree->Branch("year", &year, "year/I");
    summary_tree->Branch("run", &run, "run/I");
    summary_tree->Branch("NTracks", &NTracks, "NTracks/D");
    summary_tree->Branch("TriggerCat", &TriggerCat, "TriggerCat/I");
    if(_data)summary_tree->Branch("Bs_BKGCAT", &Bs_BKGCAT, "Bs_BKGCAT/I");

    double weight = 1.;
    summary_tree->Branch("weight",&weight,"weight/D"); 
    
    // Variables for BDT 
    double DsDaughters_min_IPCHI2 = 0;
    double XsDaughters_min_IPCHI2 = 0;
    double DsDaughters_max_IPCHI2 = 0;
    double XsDaughters_max_IPCHI2 = 0;
    double DsDaughters_min_PT = 0;
    double XsDaughters_min_PT = 0;
    double DsDaughters_min_P = 0;
    double XsDaughters_min_P = 0;
    double Xs_max_DOCA = 0;
    double Ds_max_DOCA = 0;
    double max_TrackChi2 = 0;
    double min_TrackChi2 = 0;
    double Xs_max_ghostProb = 0;
    double Ds_max_ghostProb = 0;
    double max_ghostProb = 0;
    double Xs_max_ProbNNghost = 0;
    double Ds_max_ProbNNghost = 0;
    double max_ProbNNghost = 0;
    double track_min_PT,track_min_P,track_min_IPCHI2,Xs_ptasy,Xs_PT,Xs_ETA;    

    summary_tree->Branch("DsDaughters_min_IPCHI2",&DsDaughters_min_IPCHI2,"DsDaughters_min_IPCHI2/D");
    summary_tree->Branch("XsDaughters_min_IPCHI2",&XsDaughters_min_IPCHI2,"XsDaughters_min_IPCHI2/D");
    summary_tree->Branch("DsDaughters_max_IPCHI2",&DsDaughters_max_IPCHI2,"DsDaughters_max_IPCHI2/D");
    summary_tree->Branch("XsDaughters_max_IPCHI2",&XsDaughters_max_IPCHI2,"XsDaughters_max_IPCHI2/D");
    summary_tree->Branch("track_min_IPCHI2",&track_min_IPCHI2,"track_min_IPCHI2/D");
    summary_tree->Branch("DsDaughters_min_PT",&DsDaughters_min_PT,"DsDaughters_min_PT/D");
    summary_tree->Branch("XsDaughters_min_PT",&XsDaughters_min_PT,"XsDaughters_min_PT/D");
    summary_tree->Branch("track_min_PT",&track_min_PT,"track_min_PT/D");
    summary_tree->Branch("DsDaughters_min_P",&DsDaughters_min_P,"DsDaughters_min_P/D");
    summary_tree->Branch("XsDaughters_min_P",&XsDaughters_min_P,"XsDaughters_min_P/D");
    summary_tree->Branch("track_min_P",&track_min_P,"track_min_P/D");
    summary_tree->Branch("Xs_ptasy_1.00",&Xs_ptasy,"Xs_ptasy_1.00/D");
    summary_tree->Branch("Xs_PT",&Xs_PT,"Xs_PT/D");
    summary_tree->Branch("Xs_ETA",&Xs_ETA,"Xs_ETA/D");

    summary_tree->Branch("Xs_max_DOCA",&Xs_max_DOCA,"Xs_max_DOCA/D");
    summary_tree->Branch("Ds_max_DOCA",&Ds_max_DOCA,"Ds_max_DOCA/D");
    summary_tree->Branch("Xs_max_ghostProb",&Xs_max_ghostProb,"Xs_max_ghostProb/D");
    summary_tree->Branch("Ds_max_ghostProb",&Ds_max_ghostProb,"Ds_max_ghostProb/D");
    summary_tree->Branch("max_ghostProb",&max_ghostProb,"max_ghostProb/D");
    summary_tree->Branch("Xs_max_ProbNNghost",&Xs_max_ProbNNghost,"Xs_max_ProbNNghost/D");
    summary_tree->Branch("Ds_max_ProbNNghost",&Ds_max_ProbNNghost,"Ds_max_ProbNNghost/D");
    summary_tree->Branch("max_ProbNNghost",&max_ProbNNghost,"max_ProbNNghost/D");

    double angK,angPip, angPim, maxCos, maxGP;
    summary_tree->Branch("angK", &angK, "angK/D");
    summary_tree->Branch("angPip", &angPip, "angPip/D");
    summary_tree->Branch("angPim", &angPim, "angPim/D");
    summary_tree->Branch("maxCos", &maxCos, "maxCos/D");

    double Bs_RFD,Ds_RFD,Ds_FDsig,Ds_z;
    summary_tree->Branch("Bs_RFD", &Bs_RFD, "Bs_RFD/D");
    summary_tree->Branch("Ds_RFD", &Ds_RFD, "Ds_RFD/D");
    summary_tree->Branch("Ds_FDsig", &Ds_FDsig, "Ds_FDsig/D");
    summary_tree->Branch("Ds_z", &Ds_z, "Ds_z/D");

    // Tagging
    Int_t         OS_Muon_DEC;
    Double_t      OS_Muon_PROB;
    Int_t         OS_Electron_DEC;
    Double_t      OS_Electron_PROB;
    Int_t         OS_Kaon_DEC;
    Double_t      OS_Kaon_PROB;
    Int_t         OS_nnetKaon_DEC;
    Double_t      OS_nnetKaon_PROB;
    Int_t         OS_Charm_DEC;
    Double_t      OS_Charm_PROB;
    Int_t         OS_VtxCharge_DEC;
    Double_t      OS_VtxCharge_PROB;
    Int_t         SS_Kaon_DEC;
    Double_t      SS_Kaon_PROB;
    
    Int_t         OS_Muon_DEC_Run1;
    Double_t      OS_Muon_PROB_Run1;
    Int_t         OS_Electron_DEC_Run1;
    Double_t      OS_Electron_PROB_Run1;
    Int_t         OS_Kaon_DEC_Run1;
    Double_t      OS_Kaon_PROB_Run1;
    Int_t         SS_Kaon_DEC_Run1;
    Double_t      SS_Kaon_PROB_Run1;
 
    summary_tree->Branch("OS_Muon_DEC", &OS_Muon_DEC, "OS_Muon_DEC/I");
    summary_tree->Branch("OS_Electron_DEC", &OS_Electron_DEC, "OS_Electron_DEC/I");
    summary_tree->Branch("OS_Kaon_DEC", &OS_Kaon_DEC, "OS_Kaon_DEC/I");
    summary_tree->Branch("OS_Charm_DEC", &OS_Charm_DEC, "OS_Charm_DEC/I");
    summary_tree->Branch("OS_VtxCharge_DEC", &OS_VtxCharge_DEC, "OS_VtxCharge_DEC/I");
    summary_tree->Branch("SS_Kaon_DEC", &SS_Kaon_DEC, "SS_Kaon_DEC/I");

    summary_tree->Branch("OS_Muon_PROB", &OS_Muon_PROB, "OS_Muon_PROB/D");
    summary_tree->Branch("OS_Electron_PROB", &OS_Electron_PROB, "OS_Electron_PROB/D");
    summary_tree->Branch("OS_Kaon_PROB", &OS_Kaon_PROB, "OS_Kaon_PROB/D");
    summary_tree->Branch("OS_Charm_PROB", &OS_Charm_PROB, "OS_Charm_PROB/D");
    summary_tree->Branch("OS_VtxCharge_PROB", &OS_VtxCharge_PROB, "OS_VtxCharge_PROB/D");
    summary_tree->Branch("SS_Kaon_PROB", &SS_Kaon_PROB, "SS_Kaon_PROB/D");
    
    summary_tree->Branch("OS_nnetKaon_DEC", &OS_nnetKaon_DEC, "OS_nnetKaon_DEC/I");
    summary_tree->Branch("OS_nnetKaon_PROB", &OS_nnetKaon_PROB, "OS_nnetKaon_PROB/D");
    summary_tree->Branch("OS_Muon_DEC_Run1", &OS_Muon_DEC_Run1, "OS_Muon_DEC_Run1/I");
    summary_tree->Branch("OS_Electron_DEC_Run1", &OS_Electron_DEC_Run1, "OS_Electron_DEC_Run1/I");
    summary_tree->Branch("OS_Kaon_DEC_Run1", &OS_Kaon_DEC_Run1, "OS_Kaon_DEC_Run1/I");
    summary_tree->Branch("SS_Kaon_DEC_Run1", &SS_Kaon_DEC_Run1, "SS_Kaon_DEC_Run1/I");
    summary_tree->Branch("OS_Muon_PROB_Run1", &OS_Muon_PROB_Run1, "OS_Muon_PROB_Run1/D");
    summary_tree->Branch("OS_Electron_PROB_Run1", &OS_Electron_PROB_Run1, "OS_Electron_PROB_Run1/D");
    summary_tree->Branch("OS_Kaon_PROB_Run1", &OS_Kaon_PROB_Run1, "OS_Kaon_PROB_Run1/D");
    summary_tree->Branch("SS_Kaon_PROB_Run1", &SS_Kaon_PROB_Run1, "SS_Kaon_PROB_Run1/D");

    // Bkg studies
    double bkg_D_as_Ds_m, bkg_Lambdac_as_Ds_m;
    summary_tree->Branch("bkg_D_as_Ds_m", &bkg_D_as_Ds_m, "bkg_D_as_Ds_m/D");
    summary_tree->Branch("bkg_Lambdac_as_Ds_m", &bkg_Lambdac_as_Ds_m, "bkg_Lambdac_as_Ds_m/D");
    double bkg_D_as_Ds_Bs_m, bkg_Lambdac_as_Ds_Bs_m;
    summary_tree->Branch("bkg_D_as_Ds_Bs_m", &bkg_D_as_Ds_Bs_m, "bkg_D_as_Ds_Bs_m/D");
    summary_tree->Branch("bkg_Lambdac_as_Ds_Bs_m", &bkg_Lambdac_as_Ds_Bs_m, "bkg_Lambdac_as_Ds_Bs_m/D");
    
    double bkg_Dstar_as_Ds_dm1,bkg_Dstar_as_Ds_dm2,bkg_Dstar_as_Ds_dm3,bkg_Dstar_as_Ds_dm4;
    summary_tree->Branch("bkg_Dstar_as_Ds_dm1",&bkg_Dstar_as_Ds_dm1,"bkg_Dstar_as_Ds_dm1/D");
    summary_tree->Branch("bkg_Dstar_as_Ds_dm2",&bkg_Dstar_as_Ds_dm2,"bkg_Dstar_as_Ds_dm2/D");
    summary_tree->Branch("bkg_Dstar_as_Ds_dm3",&bkg_Dstar_as_Ds_dm3,"bkg_Dstar_as_Ds_dm3/D");
    summary_tree->Branch("bkg_Dstar_as_Ds_dm4",&bkg_Dstar_as_Ds_dm4,"bkg_Dstar_as_Ds_dm4/D");

    double bkg_Dstar_as_Xs_dm1,bkg_Dstar_as_Xs_dm2,bkg_Dstar_as_Xs_dm3;
    summary_tree->Branch("bkg_Dstar_as_Xs_dm1",&bkg_Dstar_as_Xs_dm1,"bkg_Dstar_as_Xs_dm1/D");
    summary_tree->Branch("bkg_Dstar_as_Xs_dm2",&bkg_Dstar_as_Xs_dm2,"bkg_Dstar_as_Xs_dm2/D");
    summary_tree->Branch("bkg_Dstar_as_Xs_dm3",&bkg_Dstar_as_Xs_dm3,"bkg_Dstar_as_Xs_dm3/D");
    
    double bkg_KKpi_as_Xs_m, bkg_3pi_as_Xs_m;
    summary_tree->Branch("bkg_KKpi_as_Xs_m", &bkg_KKpi_as_Xs_m, "bkg_KKpi_as_Xs_m/D");
    summary_tree->Branch("bkg_3pi_as_Xs_m", &bkg_3pi_as_Xs_m, "bkg_3pi_as_Xs_m/D");

    double bkg_Kppi_as_Xs_m,bkg_Kpip_as_Xs_m,bkg_ppipi_as_Xs_m;
    summary_tree->Branch("bkg_Kppi_as_Xs_m", &bkg_Kppi_as_Xs_m, "bkg_Kppi_as_Xs_m/D");
    summary_tree->Branch("bkg_Kpip_as_Xs_m", &bkg_Kpip_as_Xs_m, "bkg_Kpip_as_Xs_m/D");
    summary_tree->Branch("bkg_ppipi_as_Xs_m", &bkg_ppipi_as_Xs_m, "bkg_ppipi_as_Xs_m/D");

    double bkg_KKpi_as_Xs_Bs_m, bkg_3pi_as_Xs_Bs_m;
    summary_tree->Branch("bkg_KKpi_as_Xs_Bs_m", &bkg_KKpi_as_Xs_Bs_m, "bkg_KKpi_as_Xs_Bs_m/D");
    summary_tree->Branch("bkg_3pi_as_Xs_Bs_m", &bkg_3pi_as_Xs_Bs_m, "bkg_3pi_as_Xs_Bs_m/D");

    double bkg_Ks_as_KK_fromDs_m;
    summary_tree->Branch("bkg_Ks_as_KK_fromDs_m", &bkg_Ks_as_KK_fromDs_m, "bkg_Ks_as_KK_fromDs_m/D");

    double bkg_Ks_as_pipi1_m;
    summary_tree->Branch("bkg_Ks_as_pipi1_m", &bkg_Ks_as_pipi1_m, "bkg_Ks_as_pipi1_m/D");
    double bkg_Ks_as_pipi2_m;
    summary_tree->Branch("bkg_Ks_as_pipi2_m", &bkg_Ks_as_pipi2_m, "bkg_Ks_as_pipi2_m/D");
    double bkg_rho_as_Kpi_m;
    summary_tree->Branch("bkg_rho_as_Kpi_m", &bkg_rho_as_Kpi_m, "bkg_rho_as_Kpi_m/D");
    double bkg_phi_as_Kpi_m;
    summary_tree->Branch("bkg_phi_as_Kpi_m", &bkg_phi_as_Kpi_m, "bkg_phi_as_Kpi_m/D");

    double m_15,m_16,m_25,m_26,m_34;
    summary_tree->Branch("m_15", &m_15, "m_15/D");
    summary_tree->Branch("m_16", &m_16, "m_16/D");
    summary_tree->Branch("m_25", &m_25, "m_25/D");
    summary_tree->Branch("m_26", &m_26, "m_26/D");
    summary_tree->Branch("m_34", &m_34, "m_34/D");

    double m_135,m_156,m_345,m_125,m_245;
    summary_tree->Branch("m_135", &m_135, "m_135/D");
    summary_tree->Branch("m_156", &m_156, "m_156/D");
    summary_tree->Branch("m_345", &m_345, "m_345/D");
    summary_tree->Branch("m_125", &m_125, "m_125/D");
    summary_tree->Branch("m_245", &m_245, "m_245/D");

    double m_236,m_145,m_126,m_234,m_246,m_235,m_256,m_136,m_346;
    summary_tree->Branch("m_236", &m_236, "m_236/D");
    summary_tree->Branch("m_145", &m_145, "m_145/D");
    summary_tree->Branch("m_126", &m_126, "m_126/D");
    summary_tree->Branch("m_234", &m_234, "m_234/D");
    summary_tree->Branch("m_246", &m_246, "m_246/D");
    summary_tree->Branch("m_235", &m_235, "m_235/D");
    summary_tree->Branch("m_256", &m_256, "m_256/D");
    summary_tree->Branch("m_136", &m_136, "m_136/D");
    summary_tree->Branch("m_346", &m_346, "m_346/D");

    double m_1236,m_2346,m_1235,m_1256,m_2345,m_2456;
    summary_tree->Branch("m_1236", &m_1236, "m_1236/D");
    summary_tree->Branch("m_2346", &m_2346, "m_2346/D");
    summary_tree->Branch("m_1235", &m_1235, "m_1235/D");
    summary_tree->Branch("m_1256", &m_1256, "m_1256/D");
    summary_tree->Branch("m_2345", &m_2345, "m_2345/D");
    summary_tree->Branch("m_2456", &m_2456, "m_2456/D");

    double beta_K_plus;
    double beta_pi_plus;
    double beta_pi_minus;
    double beta_Ds;
    summary_tree->Branch("beta_K_plus", &beta_K_plus, "beta_K_plus/D");
    summary_tree->Branch("beta_pi_plus", &beta_pi_plus, "beta_pi_plus/D");
    summary_tree->Branch("beta_pi_minus", &beta_pi_minus, "beta_pi_minus/D");
    summary_tree->Branch("beta_Ds", &beta_Ds, "beta_Ds/D");
    
    double beta_K_plus_fromDs;
    double beta_pi_minus_fromDs;
    double beta_K_minus_fromDs;
    summary_tree->Branch("beta_K_plus_fromDs", &beta_K_plus_fromDs, "beta_K_plus_fromDs/D");
    summary_tree->Branch("beta_pi_minus_fromDs", &beta_pi_minus_fromDs, "beta_pi_minus_fromDs/D");
    summary_tree->Branch("beta_K_minus_fromDs", &beta_K_minus_fromDs, "beta_K_minus_fromDs/D");

    // Dalitz stuff
    double Ds_m12,Ds_m13;
    summary_tree->Branch("Ds_m12", &Ds_m12, "Ds_m12/D");
    summary_tree->Branch("Ds_m13", &Ds_m13, "Ds_m13/D");

    double m_DsK,m_Dspi,m_Dspipi,m_DsKpi,m_DsKpip,m_Kpipi,m_Kpi,m_pipi;
    summary_tree->Branch("m_DsK", &m_DsK, "m_DsK/D");
    summary_tree->Branch("m_Dspi", &m_Dspi, "m_Dspi/D");
    summary_tree->Branch("m_Dspipi", &m_Dspipi, "m_Dspipi/D");
    summary_tree->Branch("m_DsKpi", &m_DsKpi, "m_DsKpi/D");
    summary_tree->Branch("m_DsKpip", &m_DsKpip, "m_DsKpip/D");
    summary_tree->Branch("m_Kpipi", &m_Kpipi, "m_Kpipi/D");
    summary_tree->Branch("m_Kpi", &m_Kpi, "m_Kpi/D");
    summary_tree->Branch("m_pipi", &m_pipi, "m_pipi/D");

    // DTF stuff
    double PV_status, DTF_status, BsDTF_status;
    double PV_CHI2NDOF, DTF_CHI2NDOF, BsDTF_CHI2NDOF;
    summary_tree->Branch("PV_status", &PV_status, "PV_status/D");
    summary_tree->Branch("DTF_status", &DTF_status, "DTF_status/D");
    summary_tree->Branch("BsDTF_status", &BsDTF_status, "BsDTF_status/D");
    summary_tree->Branch("PV_CHI2NDOF", &PV_CHI2NDOF, "PV_CHI2NDOF/D");
    summary_tree->Branch("DTF_CHI2NDOF", &DTF_CHI2NDOF, "DTF_CHI2NDOF/D");
    summary_tree->Branch("BsDTF_CHI2NDOF", &BsDTF_CHI2NDOF, "BsDTF_CHI2NDOF/D");

    double Bs_PV_TAU,Bs_PV_TAUERR;
    double Bs_DTF_TAU,Bs_DTF_TAUERR;
    double Bs_BsDTF_TAU, Bs_BsDTF_TAUERR;
    summary_tree->Branch("Bs_PV_TAU", &Bs_PV_TAU, "Bs_PV_TAU/D");
    summary_tree->Branch("Bs_PV_TAUERR", &Bs_PV_TAUERR, "Bs_PV_TAUERR/D");
    summary_tree->Branch("Bs_DTF_TAU", &Bs_DTF_TAU, "Bs_DTF_TAU/D");
    summary_tree->Branch("Bs_DTF_TAUERR", &Bs_DTF_TAUERR, "Bs_DTF_TAUERR/D");
    summary_tree->Branch("Bs_BsDTF_TAU", &Bs_BsDTF_TAU, "Bs_BsDTF_TAU/D");
    summary_tree->Branch("Bs_BsDTF_TAUERR", &Bs_BsDTF_TAUERR, "Bs_BsDTF_TAUERR/D");
    
    double Bs_DTF_MM,Bs_DTF_MMERR;
    summary_tree->Branch("Bs_DTF_MM", &Bs_DTF_MM, "Bs_DTF_MM/D");
    summary_tree->Branch("Bs_DTF_MMERR", &Bs_DTF_MMERR, "Bs_DTF_MMERR/D");
    double Bs_PV_MM,Bs_PV_MMERR;
    summary_tree->Branch("Bs_PV_MM", &Bs_PV_MM, "Bs_PV_MM/D");
    summary_tree->Branch("Bs_PV_MMERR", &Bs_PV_MMERR, "Bs_PV_MMERR/D");    
    double Ds_PV_MM,Ds_PV_MMERR;
    summary_tree->Branch("Ds_PV_MM", &Ds_PV_MM, "Ds_PV_MM/D");
    summary_tree->Branch("Ds_PV_MMERR", &Ds_PV_MMERR, "Ds_PV_MMERR/D");    
    
    double BsDTF_Kplus_PX,BsDTF_Kplus_PY,BsDTF_Kplus_PZ,BsDTF_Kplus_PE;
    summary_tree->Branch("BsDTF_Kplus_PX", &BsDTF_Kplus_PX, "BsDTF_Kplus_PX/D");
    summary_tree->Branch("BsDTF_Kplus_PY", &BsDTF_Kplus_PY, "BsDTF_Kplus_PY/D");
    summary_tree->Branch("BsDTF_Kplus_PZ", &BsDTF_Kplus_PZ, "BsDTF_Kplus_PZ/D");
    summary_tree->Branch("BsDTF_Kplus_PE", &BsDTF_Kplus_PE, "BsDTF_Kplus_PE/D");
    
    double BsDTF_piplus_PX, BsDTF_piplus_PY,BsDTF_piplus_PZ,BsDTF_piplus_PE;
    summary_tree->Branch("BsDTF_piplus_PX", &BsDTF_piplus_PX, "BsDTF_piplus_PX/D");
    summary_tree->Branch("BsDTF_piplus_PY", &BsDTF_piplus_PY, "BsDTF_piplus_PY/D");
    summary_tree->Branch("BsDTF_piplus_PZ", &BsDTF_piplus_PZ, "BsDTF_piplus_PZ/D");
    summary_tree->Branch("BsDTF_piplus_PE", &BsDTF_piplus_PE, "BsDTF_piplus_PE/D");
    
    double BsDTF_piminus_PX,BsDTF_piminus_PY,BsDTF_piminus_PZ,BsDTF_piminus_PE;
    summary_tree->Branch("BsDTF_piminus_PX", &BsDTF_piminus_PX, "BsDTF_piminus_PX/D");
    summary_tree->Branch("BsDTF_piminus_PY", &BsDTF_piminus_PY, "BsDTF_piminus_PY/D");
    summary_tree->Branch("BsDTF_piminus_PZ", &BsDTF_piminus_PZ, "BsDTF_piminus_PZ/D");
    summary_tree->Branch("BsDTF_piminus_PE", &BsDTF_piminus_PE, "BsDTF_piminus_PE/D");
    
    double BsDTF_Ds_Kplus_PX,BsDTF_Ds_Kplus_PY,BsDTF_Ds_Kplus_PZ,BsDTF_Ds_Kplus_PE;
    summary_tree->Branch("BsDTF_Ds_Kplus_PX", &BsDTF_Ds_Kplus_PX, "BsDTF_Ds_Kplus_PX/D");
    summary_tree->Branch("BsDTF_Ds_Kplus_PY", &BsDTF_Ds_Kplus_PY, "BsDTF_Ds_Kplus_PY/D");
    summary_tree->Branch("BsDTF_Ds_Kplus_PZ", &BsDTF_Ds_Kplus_PZ, "BsDTF_Ds_Kplus_PZ/D");
    summary_tree->Branch("BsDTF_Ds_Kplus_PE", &BsDTF_Ds_Kplus_PE, "BsDTF_Ds_Kplus_PE/D");
    
    double BsDTF_Ds_Kminus_PX,BsDTF_Ds_Kminus_PY,BsDTF_Ds_Kminus_PZ,BsDTF_Ds_Kminus_PE;
    summary_tree->Branch("BsDTF_Ds_Kminus_PX", &BsDTF_Ds_Kminus_PX, "BsDTF_Ds_Kminus_PX/D");
    summary_tree->Branch("BsDTF_Ds_Kminus_PY", &BsDTF_Ds_Kminus_PY, "BsDTF_Ds_Kminus_PY/D");
    summary_tree->Branch("BsDTF_Ds_Kminus_PZ", &BsDTF_Ds_Kminus_PZ, "BsDTF_Ds_Kminus_PZ/D");
    summary_tree->Branch("BsDTF_Ds_Kminus_PE", &BsDTF_Ds_Kminus_PE, "BsDTF_Ds_Kminus_PE/D");    

    double BsDTF_Ds_piminus_PX,BsDTF_Ds_piminus_PY,BsDTF_Ds_piminus_PZ,BsDTF_Ds_piminus_PE;
    summary_tree->Branch("BsDTF_Ds_piminus_PX", &BsDTF_Ds_piminus_PX, "BsDTF_Ds_piminus_PX/D");
    summary_tree->Branch("BsDTF_Ds_piminus_PY", &BsDTF_Ds_piminus_PY, "BsDTF_Ds_piminus_PY/D");
    summary_tree->Branch("BsDTF_Ds_piminus_PZ", &BsDTF_Ds_piminus_PZ, "BsDTF_Ds_piminus_PZ/D");
    summary_tree->Branch("BsDTF_Ds_piminus_PE", &BsDTF_Ds_piminus_PE, "BsDTF_Ds_piminus_PE/D");   

    double BsDTF_Ds_PX,BsDTF_Ds_PY,BsDTF_Ds_PZ,BsDTF_Ds_PE;
    summary_tree->Branch("BsDTF_Ds_PX", &BsDTF_Ds_PX, "BsDTF_Ds_PX/D");
    summary_tree->Branch("BsDTF_Ds_PY", &BsDTF_Ds_PY, "BsDTF_Ds_PY/D");
    summary_tree->Branch("BsDTF_Ds_PZ", &BsDTF_Ds_PZ, "BsDTF_Ds_PZ/D");
    summary_tree->Branch("BsDTF_Ds_PE", &BsDTF_Ds_PE, "BsDTF_Ds_PE/D");   

    // MC PID
    double K_plus_PIDK_gen, K_plus_PIDK_corr, K_plus_PIDK_raw;
    double pi_plus_PIDK_gen, pi_plus_PIDK_corr, pi_plus_PIDK_raw;
    double pi_minus_PIDK_gen, pi_minus_PIDK_corr, pi_minus_PIDK_raw;    
    if(!_data){
	summary_tree->Branch("K_plus_PIDK_gen", &K_plus_PIDK_gen, "K_plus_PIDK_gen/D");
	summary_tree->Branch("K_plus_PIDK_corr", &K_plus_PIDK_corr, "K_plus_PIDK_corr/D");
	summary_tree->Branch("K_plus_PIDK_raw", &K_plus_PIDK_raw, "K_plus_PIDK_raw/D");
	summary_tree->Branch("pi_plus_PIDK_gen", &pi_plus_PIDK_gen, "pi_plus_PIDK_gen/D");
	summary_tree->Branch("pi_plus_PIDK_corr", &pi_plus_PIDK_corr, "pi_plus_PIDK_corr/D");
	summary_tree->Branch("pi_plus_PIDK_raw", &pi_plus_PIDK_raw, "pi_plus_PIDK_raw/D");
	summary_tree->Branch("pi_minus_PIDK_gen", &pi_minus_PIDK_gen, "pi_minus_PIDK_gen/D");
	summary_tree->Branch("pi_minus_PIDK_corr", &pi_minus_PIDK_corr, "pi_minus_PIDK_corr/D");
	summary_tree->Branch("pi_minus_PIDK_raw", &pi_minus_PIDK_raw, "pi_minus_PIDK_raw/D");
    }
    double pi_plus1_PIDK_gen, pi_plus1_PIDK_corr, pi_plus1_PIDK_raw;
    double pi_plus2_PIDK_gen, pi_plus2_PIDK_corr, pi_plus2_PIDK_raw;
    if(!_data){
	summary_tree->Branch("pi_plus1_PIDK_gen", &pi_plus1_PIDK_gen, "pi_plus1_PIDK_gen/D");
	summary_tree->Branch("pi_plus1_PIDK_corr", &pi_plus1_PIDK_corr, "pi_plus1_PIDK_corr/D");
	summary_tree->Branch("pi_plus1_PIDK_raw", &pi_plus1_PIDK_raw, "pi_plus1_PIDK_raw/D");
	summary_tree->Branch("pi_plus2_PIDK_gen", &pi_plus2_PIDK_gen, "pi_plus2_PIDK_gen/D");
	summary_tree->Branch("pi_plus2_PIDK_corr", &pi_plus2_PIDK_corr, "pi_plus2_PIDK_corr/D");
	summary_tree->Branch("pi_plus2_PIDK_raw", &pi_plus2_PIDK_raw, "pi_plus2_PIDK_raw/D");
    }
    double pi_minus_fromDs_PIDK_gen, pi_minus_fromDs_PIDK_corr, pi_minus_fromDs_PIDK_raw;
    double pi_minus2_fromDs_PIDK_gen, pi_minus2_fromDs_PIDK_corr, pi_minus2_fromDs_PIDK_raw;
    double K_minus_fromDs_PIDK_gen, K_minus_fromDs_PIDK_corr, K_minus_fromDs_PIDK_raw;
    double K_plus_fromDs_PIDK_gen, K_plus_fromDs_PIDK_corr, K_plus_fromDs_PIDK_raw;
    double pi_plus_fromDs_PIDK_gen, pi_plus_fromDs_PIDK_corr, pi_plus_fromDs_PIDK_raw;
    if(!_data){
	summary_tree->Branch("pi_minus_fromDs_PIDK_gen", &pi_minus_fromDs_PIDK_gen, "pi_minus_fromDs_PIDK_gen/D");
	summary_tree->Branch("pi_minus_fromDs_PIDK_corr", &pi_minus_fromDs_PIDK_corr, "pi_minus_fromDs_PIDK_corr/D");
	summary_tree->Branch("pi_minus_fromDs_PIDK_raw", &pi_minus_fromDs_PIDK_raw, "pi_minus_fromDs_PIDK_raw/D");
	summary_tree->Branch("pi_minus2_fromDs_PIDK_gen", &pi_minus2_fromDs_PIDK_gen, "pi_minus2_fromDs_PIDK_gen/D");
	summary_tree->Branch("pi_minus2_fromDs_PIDK_corr", &pi_minus2_fromDs_PIDK_corr, "pi_minus2_fromDs_PIDK_corr/D");
	summary_tree->Branch("pi_minus2_fromDs_PIDK_raw", &pi_minus2_fromDs_PIDK_raw, "pi_minus2_fromDs_PIDK_raw/D");
	summary_tree->Branch("K_minus_fromDs_PIDK_gen", &K_minus_fromDs_PIDK_gen, "K_minus_fromDs_PIDK_gen/D");
	summary_tree->Branch("K_minus_fromDs_PIDK_corr", &K_minus_fromDs_PIDK_corr, "K_minus_fromDs_PIDK_corr/D");
	summary_tree->Branch("K_minus_fromDs_PIDK_raw", &K_minus_fromDs_PIDK_raw, "K_minus_fromDs_PIDK_raw/D");
	summary_tree->Branch("K_plus_fromDs_PIDK_gen", &K_plus_fromDs_PIDK_gen, "K_plus_fromDs_PIDK_gen/D");
	summary_tree->Branch("K_plus_fromDs_PIDK_corr", &K_plus_fromDs_PIDK_corr, "K_plus_fromDs_PIDK_corr/D");
	summary_tree->Branch("K_plus_fromDs_PIDK_raw", &K_plus_fromDs_PIDK_raw, "K_plus_fromDs_PIDK_raw/D");
	summary_tree->Branch("pi_plus_fromDs_PIDK_gen", &pi_plus_fromDs_PIDK_gen, "pi_plus_fromDs_PIDK_gen/D");
	summary_tree->Branch("pi_plus_fromDs_PIDK_corr", &pi_plus_fromDs_PIDK_corr, "pi_plus_fromDs_PIDK_corr/D");
	summary_tree->Branch("pi_plus_fromDs_PIDK_raw", &pi_plus_fromDs_PIDK_raw, "pi_plus_fromDs_PIDK_raw/D");
    }
    /*
    double K_plus_PIDp_gen, K_plus_PIDp_corr, K_plus_PIDp_raw;
    double pi_plus_PIDp_gen, pi_plus_PIDp_corr, pi_plus_PIDp_raw;
    double pi_minus_PIDp_gen, pi_minus_PIDp_corr, pi_minus_PIDp_raw;    
    summary_tree->Branch("K_plus_PIDp_gen", &K_plus_PIDp_gen, "K_plus_PIDp_gen/D");
    summary_tree->Branch("K_plus_PIDp_corr", &K_plus_PIDp_corr, "K_plus_PIDp_corr/D");
    summary_tree->Branch("K_plus_PIDp_raw", &K_plus_PIDp_raw, "K_plus_PIDp_raw/D");
    summary_tree->Branch("pi_plus_PIDp_gen", &pi_plus_PIDp_gen, "pi_plus_PIDp_gen/D");
    summary_tree->Branch("pi_plus_PIDp_corr", &pi_plus_PIDp_corr, "pi_plus_PIDp_corr/D");
    summary_tree->Branch("pi_plus_PIDp_raw", &pi_plus_PIDp_raw, "pi_plus_PIDp_raw/D");
    summary_tree->Branch("pi_minus_PIDp_gen", &pi_minus_PIDp_gen, "pi_minus_PIDp_gen/D");
    summary_tree->Branch("pi_minus_PIDp_corr", &pi_minus_PIDp_corr, "pi_minus_PIDp_corr/D");
    summary_tree->Branch("pi_minus_PIDp_raw", &pi_minus_PIDp_raw, "pi_minus_PIDp_raw/D");

    double pi_plus1_PIDp_gen, pi_plus1_PIDp_corr, pi_plus1_PIDp_raw;
    double pi_plus2_PIDp_gen, pi_plus2_PIDp_corr, pi_plus2_PIDp_raw;
    summary_tree->Branch("pi_plus1_PIDp_gen", &pi_plus1_PIDp_gen, "pi_plus1_PIDp_gen/D");
    summary_tree->Branch("pi_plus1_PIDp_corr", &pi_plus1_PIDp_corr, "pi_plus1_PIDp_corr/D");
    summary_tree->Branch("pi_plus1_PIDp_raw", &pi_plus1_PIDp_raw, "pi_plus1_PIDp_raw/D");
    summary_tree->Branch("pi_plus2_PIDp_gen", &pi_plus2_PIDp_gen, "pi_plus2_PIDp_gen/D");
    summary_tree->Branch("pi_plus2_PIDp_corr", &pi_plus2_PIDp_corr, "pi_plus2_PIDp_corr/D");
    summary_tree->Branch("pi_plus2_PIDp_raw", &pi_plus2_PIDp_raw, "pi_plus2_PIDp_raw/D");

    double pi_minus_fromDs_PIDp_gen, pi_minus_fromDs_PIDp_corr, pi_minus_fromDs_PIDp_raw;
    double K_minus_fromDs_PIDp_gen, K_minus_fromDs_PIDp_corr, K_minus_fromDs_PIDp_raw;
    double K_plus_fromDs_PIDp_gen, K_plus_fromDs_PIDp_corr, K_plus_fromDs_PIDp_raw;
    summary_tree->Branch("pi_minus_fromDs_PIDp_gen", &pi_minus_fromDs_PIDp_gen, "pi_minus_fromDs_PIDp_gen/D");
    summary_tree->Branch("pi_minus_fromDs_PIDp_corr", &pi_minus_fromDs_PIDp_corr, "pi_minus_fromDs_PIDp_corr/D");
    summary_tree->Branch("pi_minus_fromDs_PIDp_raw", &pi_minus_fromDs_PIDp_raw, "pi_minus_fromDs_PIDp_raw/D");
    summary_tree->Branch("K_minus_fromDs_PIDp_gen", &K_minus_fromDs_PIDp_gen, "K_minus_fromDs_PIDp_gen/D");
    summary_tree->Branch("K_minus_fromDs_PIDp_corr", &K_minus_fromDs_PIDp_corr, "K_minus_fromDs_PIDp_corr/D");
    summary_tree->Branch("K_minus_fromDs_PIDp_raw", &K_minus_fromDs_PIDp_raw, "K_minus_fromDs_PIDp_raw/D");
    summary_tree->Branch("K_plus_fromDs_PIDp_gen", &K_plus_fromDs_PIDp_gen, "K_plus_fromDs_PIDp_gen/D");
    summary_tree->Branch("K_plus_fromDs_PIDp_corr", &K_plus_fromDs_PIDp_corr, "K_plus_fromDs_PIDp_corr/D");
    summary_tree->Branch("K_plus_fromDs_PIDp_raw", &K_plus_fromDs_PIDp_raw, "K_plus_fromDs_PIDp_raw/D");
    */

    TRandom3 r;
    for (Long64_t i=0; i<nentries;i++) {
        
        if(0ul == (i % 10000ul)) cout << "Read event " << i << "/" << nentries <<
        "  ( " << i/(double)nentries * 100. << " % )" << endl;

        fChain->GetEntry(i);           
        
        // Apply preselection cust
        if(!_ltu)if(!Preselection_Cuts()) continue;

	if(!_data){
		if(_decay == Decay::signal){
				K_plus_PIDK_raw = K_plus_PIDK;
				pi_plus_PIDK_raw = pi_plus_PIDK;
				pi_minus_PIDK_raw = pi_minus_PIDK;

				//K_plus_PIDp_raw = K_plus_PIDp;
				//pi_plus_PIDp_raw = pi_plus_PIDp;
				//pi_minus_PIDp_raw = pi_minus_PIDp;

				if(Polarity == 1){ 
					K_plus_PIDK_corr = K_plus_PIDK_corr_MagUp;	
					K_plus_PIDK_gen = K_plus_PIDK_gen_MagDown;
					pi_plus_PIDK_corr = pi_plus_PIDK_corr_MagUp;	
					pi_plus_PIDK_gen = pi_plus_PIDK_gen_MagDown;
					pi_minus_PIDK_corr = pi_minus_PIDK_corr_MagUp;	
					pi_minus_PIDK_gen = pi_minus_PIDK_gen_MagDown;

					//K_plus_PIDp_corr = K_plus_PIDp_corr_MagUp;	
					//K_plus_PIDp_gen = K_plus_PIDp_gen_MagUp;
					//pi_plus_PIDp_corr = pi_plus_PIDp_corr_MagUp;	
					//pi_plus_PIDp_gen = pi_plus_PIDp_gen_MagUp;
					//pi_minus_PIDp_corr = pi_minus_PIDp_corr_MagUp;	
					//pi_minus_PIDp_gen = pi_minus_PIDp_gen_MagUp;
				}	
				else { 
					K_plus_PIDK_corr = K_plus_PIDK_corr_MagDown;
					K_plus_PIDK_gen = K_plus_PIDK_gen_MagDown;
					pi_plus_PIDK_corr = pi_plus_PIDK_corr_MagDown;	
					pi_plus_PIDK_gen = pi_plus_PIDK_gen_MagDown;
					pi_minus_PIDK_corr = pi_minus_PIDK_corr_MagDown;	
					pi_minus_PIDK_gen = pi_minus_PIDK_gen_MagDown;

					//K_plus_PIDp_corr = K_plus_PIDp_corr_MagDown;
					//K_plus_PIDp_gen = K_plus_PIDp_gen_MagDown;
					//pi_plus_PIDp_corr = pi_plus_PIDp_corr_MagDown;	
					//pi_plus_PIDp_gen = pi_plus_PIDp_gen_MagDown;
					//pi_minus_PIDp_corr = pi_minus_PIDp_corr_MagDown;	
					//pi_minus_PIDp_gen = pi_minus_PIDp_gen_MagDown;
				}
				if(_usePIDvar == "Corr"){ 
					K_plus_PIDK = K_plus_PIDK_corr;
					pi_plus_PIDK = pi_plus_PIDK_corr;
					pi_minus_PIDK = pi_minus_PIDK_corr;

					//K_plus_PIDp = K_plus_PIDp_corr;
					//pi_plus_PIDp = pi_plus_PIDp_corr;
					//pi_minus_PIDp = pi_minus_PIDp_corr;
				}
				else if(_usePIDvar == "Gen") { 
					K_plus_PIDK = K_plus_PIDK_gen;
					pi_plus_PIDK = pi_plus_PIDK_gen;
					pi_minus_PIDK = pi_minus_PIDK_gen;
					//K_plus_PIDp = K_plus_PIDp_gen;
					//pi_plus_PIDp = pi_plus_PIDp_gen;
					//pi_minus_PIDp = pi_minus_PIDp_gen;
				}
		}
		else {
				pi_plus1_PIDK_raw = pi_plus1_PIDK;
				pi_plus2_PIDK_raw = pi_plus2_PIDK;
				pi_minus_PIDK_raw = pi_minus_PIDK;

				//pi_plus1_PIDp_raw = pi_plus1_PIDp;
				//pi_plus2_PIDp_raw = pi_plus2_PIDp;
				//pi_minus_PIDp_raw = pi_minus_PIDp;

				if(Polarity == 1){ 
					pi_plus1_PIDK_corr = pi_plus1_PIDK_corr_MagUp;	
					pi_plus1_PIDK_gen = pi_plus1_PIDK_gen_MagDown;
					pi_plus2_PIDK_corr = pi_plus2_PIDK_corr_MagUp;	
					pi_plus2_PIDK_gen = pi_plus2_PIDK_gen_MagDown;
					pi_minus_PIDK_corr = pi_minus_PIDK_corr_MagUp;	
					pi_minus_PIDK_gen = pi_minus_PIDK_gen_MagDown;

					//pi_plus1_PIDp_corr = pi_plus1_PIDp_corr_MagUp;	
					//pi_plus1_PIDp_gen = pi_plus1_PIDp_gen_MagUp;
					//pi_plus2_PIDp_corr = pi_plus2_PIDp_corr_MagUp;	
					//pi_plus2_PIDp_gen = pi_plus2_PIDp_gen_MagUp;
					//pi_minus_PIDp_corr = pi_minus_PIDp_corr_MagUp;	
					//pi_minus_PIDp_gen = pi_minus_PIDp_gen_MagUp;
				}	
				else { 
					pi_plus1_PIDK_corr = pi_plus1_PIDK_corr_MagDown;
					pi_plus1_PIDK_gen = pi_plus1_PIDK_gen_MagDown;
					pi_plus2_PIDK_corr = pi_plus2_PIDK_corr_MagDown;	
					pi_plus2_PIDK_gen = pi_plus2_PIDK_gen_MagDown;
					pi_minus_PIDK_corr = pi_minus_PIDK_corr_MagDown;	
					pi_minus_PIDK_gen = pi_minus_PIDK_gen_MagDown;

					//pi_plus1_PIDp_corr = pi_plus1_PIDp_corr_MagDown;
					//pi_plus1_PIDp_gen = pi_plus1_PIDp_gen_MagDown;
					//pi_plus2_PIDp_corr = pi_plus2_PIDp_corr_MagDown;	
					//pi_plus2_PIDp_gen = pi_plus2_PIDp_gen_MagDown;
					//pi_minus_PIDp_corr = pi_minus_PIDp_corr_MagDown;	
					//pi_minus_PIDp_gen = pi_minus_PIDp_gen_MagDown;
				}
				if(_usePIDvar == "Corr"){ 
					pi_plus1_PIDK = pi_plus1_PIDK_corr;
					pi_plus2_PIDK = pi_plus2_PIDK_corr;
					pi_minus_PIDK = pi_minus_PIDK_corr;
					//pi_plus1_PIDp = pi_plus1_PIDp_corr;
					//pi_plus2_PIDp = pi_plus2_PIDp_corr;
					//pi_minus_PIDp = pi_minus_PIDp_corr;
				}
				else if(_usePIDvar == "Gen") { 
					pi_plus1_PIDK = pi_plus1_PIDK_gen;
					pi_plus2_PIDK = pi_plus2_PIDK_gen;
					pi_minus_PIDK = pi_minus_PIDK_gen;
					//pi_plus1_PIDp = pi_plus1_PIDp_gen;
					//pi_plus2_PIDp = pi_plus2_PIDp_gen;
					//pi_minus_PIDp = pi_minus_PIDp_gen;
				}
		}

		if(_Ds_finalState == Ds_finalState::phipi){
				K_plus_fromDs_PIDK_raw = K_plus_fromDs_PIDK;
				K_minus_fromDs_PIDK_raw = K_minus_fromDs_PIDK;
				pi_minus_fromDs_PIDK_raw = pi_minus_fromDs_PIDK;

				//K_plus_fromDs_PIDp_raw = K_plus_fromDs_PIDp;
				//K_minus_fromDs_PIDp_raw = K_minus_fromDs_PIDp;
				//pi_minus_fromDs_PIDp_raw = pi_minus_fromDs_PIDp;

				if(Polarity == 1){ 
					K_plus_fromDs_PIDK_corr = K_plus_fromDs_PIDK_corr_MagUp;	
					K_plus_fromDs_PIDK_gen = K_plus_fromDs_PIDK_gen_MagDown;
					K_minus_fromDs_PIDK_corr = K_minus_fromDs_PIDK_corr_MagUp;	
					K_minus_fromDs_PIDK_gen = K_minus_fromDs_PIDK_gen_MagDown;
					pi_minus_fromDs_PIDK_corr = pi_minus_fromDs_PIDK_corr_MagUp;	
					pi_minus_fromDs_PIDK_gen = pi_minus_fromDs_PIDK_gen_MagDown;

					//K_plus_fromDs_PIDp_corr = K_plus_fromDs_PIDp_corr_MagUp;	
					//K_plus_fromDs_PIDp_gen = K_plus_fromDs_PIDp_gen_MagUp;
					//K_minus_fromDs_PIDp_corr = K_minus_fromDs_PIDp_corr_MagUp;	
					//K_minus_fromDs_PIDp_gen = K_minus_fromDs_PIDp_gen_MagUp;
					//pi_minus_fromDs_PIDp_corr = pi_minus_fromDs_PIDp_corr_MagUp;	
					//pi_minus_fromDs_PIDp_gen = pi_minus_fromDs_PIDp_gen_MagUp;
				}	
				else { 
					K_plus_fromDs_PIDK_corr = K_plus_fromDs_PIDK_corr_MagDown;
					K_plus_fromDs_PIDK_gen = K_plus_fromDs_PIDK_gen_MagDown;
					K_minus_fromDs_PIDK_corr = K_minus_fromDs_PIDK_corr_MagDown;	
					K_minus_fromDs_PIDK_gen = K_minus_fromDs_PIDK_gen_MagDown;
					pi_minus_fromDs_PIDK_corr = pi_minus_fromDs_PIDK_corr_MagDown;	
					pi_minus_fromDs_PIDK_gen = pi_minus_fromDs_PIDK_gen_MagDown;

					//K_plus_fromDs_PIDp_corr = K_plus_fromDs_PIDp_corr_MagDown;
					//K_plus_fromDs_PIDp_gen = K_plus_fromDs_PIDp_gen_MagDown;
					//K_minus_fromDs_PIDp_corr = K_minus_fromDs_PIDp_corr_MagDown;	
					//K_minus_fromDs_PIDp_gen = K_minus_fromDs_PIDp_gen_MagDown;
					//pi_minus_fromDs_PIDp_corr = pi_minus_fromDs_PIDp_corr_MagDown;	
					//pi_minus_fromDs_PIDp_gen = pi_minus_fromDs_PIDp_gen_MagDown;
				}
				if(_usePIDvar == "Corr"){ 
					K_plus_fromDs_PIDK = K_plus_fromDs_PIDK_corr;
					K_minus_fromDs_PIDK = K_minus_fromDs_PIDK_corr;
					pi_minus_fromDs_PIDK = pi_minus_fromDs_PIDK_corr;

					//K_plus_fromDs_PIDp = K_plus_fromDs_PIDp_corr;
					//K_minus_fromDs_PIDp = K_minus_fromDs_PIDp_corr;
					//pi_minus_fromDs_PIDp = pi_minus_fromDs_PIDp_corr;
				}
				else if(_usePIDvar == "Gen") { 
					K_plus_fromDs_PIDK = K_plus_fromDs_PIDK_gen;
					K_minus_fromDs_PIDK = K_minus_fromDs_PIDK_gen;
					pi_minus_fromDs_PIDK = pi_minus_fromDs_PIDK_gen;
					//K_plus_fromDs_PIDp = K_plus_fromDs_PIDp_gen;
					//K_minus_fromDs_PIDp = K_minus_fromDs_PIDp_gen;
					//pi_minus_fromDs_PIDp = pi_minus_fromDs_PIDp_gen;
				}
			
		}
	
		if(_Ds_finalState == Ds_finalState::Kpipi){
				pi_plus_fromDs_PIDK_raw = pi_plus_fromDs_PIDK;
				K_minus_fromDs_PIDK_raw = K_minus_fromDs_PIDK;
				pi_minus_fromDs_PIDK_raw = pi_minus_fromDs_PIDK;

				if(Polarity == 1){ 
					pi_plus_fromDs_PIDK_corr = pi_plus_fromDs_PIDK_corr_MagUp;	
					pi_plus_fromDs_PIDK_gen = pi_plus_fromDs_PIDK_gen_MagDown;
					K_minus_fromDs_PIDK_corr = K_minus_fromDs_PIDK_corr_MagUp;	
					K_minus_fromDs_PIDK_gen = K_minus_fromDs_PIDK_gen_MagDown;
					pi_minus_fromDs_PIDK_corr = pi_minus_fromDs_PIDK_corr_MagUp;	
					pi_minus_fromDs_PIDK_gen = pi_minus_fromDs_PIDK_gen_MagDown;
				}	
				else { 
					pi_plus_fromDs_PIDK_corr = pi_plus_fromDs_PIDK_corr_MagDown;
					pi_plus_fromDs_PIDK_gen = pi_plus_fromDs_PIDK_gen_MagDown;
					K_minus_fromDs_PIDK_corr = K_minus_fromDs_PIDK_corr_MagDown;	
					K_minus_fromDs_PIDK_gen = K_minus_fromDs_PIDK_gen_MagDown;
					pi_minus_fromDs_PIDK_corr = pi_minus_fromDs_PIDK_corr_MagDown;	
					pi_minus_fromDs_PIDK_gen = pi_minus_fromDs_PIDK_gen_MagDown;
				}
				if(_usePIDvar == "Corr"){ 
					pi_plus_fromDs_PIDK = pi_plus_fromDs_PIDK_corr;
					K_minus_fromDs_PIDK = K_minus_fromDs_PIDK_corr;
					pi_minus_fromDs_PIDK = pi_minus_fromDs_PIDK_corr;
				}
				else if(_usePIDvar == "Gen") { 
					pi_plus_fromDs_PIDK = pi_plus_fromDs_PIDK_gen;
					K_minus_fromDs_PIDK = K_minus_fromDs_PIDK_gen;
					pi_minus_fromDs_PIDK = pi_minus_fromDs_PIDK_gen;
				}
		}

		if(_Ds_finalState == Ds_finalState::pipipi){
				pi_plus_fromDs_PIDK_raw = pi_plus_fromDs_PIDK;
				pi_minus2_fromDs_PIDK_raw = pi_minus2_fromDs_PIDK;
				pi_minus_fromDs_PIDK_raw = pi_minus_fromDs_PIDK;

				if(Polarity == 1){ 
					pi_plus_fromDs_PIDK_corr = pi_plus_fromDs_PIDK_corr_MagUp;	
					pi_plus_fromDs_PIDK_gen = pi_plus_fromDs_PIDK_gen_MagDown;
					pi_minus2_fromDs_PIDK_corr = pi_minus2_fromDs_PIDK_corr_MagUp;	
					pi_minus2_fromDs_PIDK_gen = pi_minus2_fromDs_PIDK_gen_MagDown;
					pi_minus_fromDs_PIDK_corr = pi_minus_fromDs_PIDK_corr_MagUp;	
					pi_minus_fromDs_PIDK_gen = pi_minus_fromDs_PIDK_gen_MagDown;
				}	
				else { 
					pi_plus_fromDs_PIDK_corr = pi_plus_fromDs_PIDK_corr_MagDown;
					pi_plus_fromDs_PIDK_gen = pi_plus_fromDs_PIDK_gen_MagDown;
					pi_minus2_fromDs_PIDK_corr = pi_minus2_fromDs_PIDK_corr_MagDown;	
					pi_minus2_fromDs_PIDK_gen = pi_minus2_fromDs_PIDK_gen_MagDown;
					pi_minus_fromDs_PIDK_corr = pi_minus_fromDs_PIDK_corr_MagDown;	
					pi_minus_fromDs_PIDK_gen = pi_minus_fromDs_PIDK_gen_MagDown;
				}
				if(_usePIDvar == "Corr"){ 
					pi_plus_fromDs_PIDK = pi_plus_fromDs_PIDK_corr;
					pi_minus2_fromDs_PIDK = pi_minus2_fromDs_PIDK_corr;
					pi_minus_fromDs_PIDK = pi_minus_fromDs_PIDK_corr;
				}
				else if(_usePIDvar == "Gen") { 
					pi_plus_fromDs_PIDK = pi_plus_fromDs_PIDK_gen;
					pi_minus2_fromDs_PIDK = pi_minus2_fromDs_PIDK_gen;
					pi_minus_fromDs_PIDK = pi_minus_fromDs_PIDK_gen;
				}
		}


	}
	
        set_LorentzVectors();    
        _Ds_finalState = get_Ds_finalState();

	if(!PID_Cuts()) continue;    
	if(!Veto_Cuts()) continue;
 
	if(_ltu){
		if(!LTU_Cuts()) continue;
	}
	else {
		if(!_bkg)if(!PhaseSpace_Cuts()) continue;
 		if(!_data)if(!MC_Cuts()) continue;        
	}

        // Add new variables
        Ds_finalState = _Ds_finalState;
        Ds_finalState_mod = _Ds_finalState;
	if(Ds_finalState_mod == 4)Ds_finalState_mod = 2;
        Ds_finalState_mod2 = 1;

        Bs_RFD = sqrt(pow(Bs_ENDVERTEX_X-Bs_OWNPV_X,2)+pow(Bs_ENDVERTEX_Y-Bs_OWNPV_Y,2));
        Ds_RFD = sqrt(pow(Ds_ENDVERTEX_X-Ds_OWNPV_X,2)+pow(Ds_ENDVERTEX_Y-Ds_OWNPV_Y,2));
        Ds_FDsig = (Ds_ENDVERTEX_Z-Ds_ORIVX_Z)/sqrt(pow(Ds_ENDVERTEX_ZERR,2)+pow(Ds_ORIVX_ZERR,2));
        Ds_z = Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z;

    	vector< TLorentzVector > tv;
    	tv.push_back(TLorentzVector(0.,0.,0.,0.));
        
        if(_decay == Decay::signal){
            TVector3 v_Ds(Ds_PX,Ds_PY,0.);
            TVector3 v_K(K_plus_PX,K_plus_PY,0.);
            TVector3 v_pip(pi_plus_PX,pi_plus_PY,0.);
            TVector3 v_pim(pi_minus_PX,pi_minus_PY,0.);
            angK= v_Ds.Angle(v_K);
            angPip= v_Ds.Angle(v_pip);
            angPim= v_Ds.Angle(v_pim);
            maxCos = cos(max(angK,max(angPip,angPim)));
            
            XsDaughters_min_IPCHI2 = min(K_plus_IPCHI2_OWNPV,min(pi_minus_IPCHI2_OWNPV,pi_plus_IPCHI2_OWNPV));            
            XsDaughters_max_IPCHI2 = max(K_plus_IPCHI2_OWNPV,max(pi_minus_IPCHI2_OWNPV,pi_plus_IPCHI2_OWNPV));            
            XsDaughters_min_PT = min(K_plus_PT,min(pi_minus_PT,pi_plus_PT));       
            XsDaughters_min_P = min(K_plus_P,min(pi_minus_P,pi_plus_P));             
            Xs_max_DOCA = max(K_1_1270_plus_DOCA1,max(K_1_1270_plus_DOCA2,K_1_1270_plus_DOCA3));
            Xs_max_ghostProb = max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,pi_minus_TRACK_GhostProb)); 
            Xs_max_ProbNNghost = max(pi_plus_ProbNNghost,max(K_plus_ProbNNghost,pi_minus_ProbNNghost)); 
	    Xs_ptasy = K_1_1270_plus_ptasy_1_00;
	    if(_ltu){
	    	Xs_PT = (K_plus+pi_plus+pi_minus).Pt();
	    	Xs_ETA = (K_plus+pi_plus+pi_minus).PseudoRapidity();
	    }
	    else {
	    	Xs_PT = K_1_1270_plus_PT;
	    	Xs_ETA = K_1_1270_plus_ETA;
	    }
            m_DsK = (BsDTF_Ds+BsDTF_K_plus).M();
            m_Dspi = (BsDTF_Ds+BsDTF_pi_plus).M();
            m_Dspipi = (BsDTF_Ds+BsDTF_pi_plus+BsDTF_pi_minus).M();
            m_DsKpi = (BsDTF_Ds+BsDTF_K_plus+BsDTF_pi_minus).M();
            m_DsKpip = (BsDTF_Ds+BsDTF_K_plus+BsDTF_pi_plus).M();
            m_Kpipi = (BsDTF_K_plus+BsDTF_pi_plus+BsDTF_pi_minus).M();
            m_Kpi = (BsDTF_K_plus+BsDTF_pi_minus).M();
            m_pipi = (BsDTF_pi_plus+BsDTF_pi_minus).M();            

            TLorentzVector pi_plus_asK_MissID; 
            pi_plus_asK_MissID.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ, massKaon);
            TLorentzVector K_plus_asPi_MissID; 
            K_plus_asPi_MissID.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ, massPion);

            bkg_KKpi_as_Xs_m = (K_plus+pi_plus+pi_minus_asK_MissID).M();
            bkg_KKpi_as_Xs_Bs_m = (Ds+K_plus+pi_plus+pi_minus_asK_MissID).M();

            TLorentzVector pi_minus_asp_MissID; 
            pi_minus_asp_MissID.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ, massProton);
            bkg_Kpip_as_Xs_m = (K_plus+pi_plus+pi_minus_asp_MissID).M();

	    TLorentzVector pi_plus_asp_MissID; 
            pi_plus_asp_MissID.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ, massProton);
            bkg_Kppi_as_Xs_m = (K_plus+pi_minus+pi_plus_asp_MissID).M();
            
	    TLorentzVector K_plus_asp_MissID; 
            K_plus_asp_MissID.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ, massProton);
            bkg_ppipi_as_Xs_m = (pi_plus+pi_minus+K_plus_asp_MissID).M();

            bkg_Ks_as_pipi1_m = (pi_plus_asK_MissID + pi_minus).M(); 
	    bkg_Ks_as_pipi2_m = (pi_minus_asK_MissID + pi_plus).M(); 
	    bkg_rho_as_Kpi_m = (K_plus_asPi_MissID + pi_minus).M();           
	    bkg_phi_as_Kpi_m = (pi_minus_asK_MissID + K_plus).M();           

	    bkg_Dstar_as_Xs_dm1 = (K_plus + pi_minus_asK_MissID + pi_plus).M() - (K_plus + pi_minus_asK_MissID).M();
	    bkg_Dstar_as_Xs_dm2 = (K_plus + pi_minus + pi_plus).M() - (K_plus + pi_minus).M();
   	    bkg_Dstar_as_Xs_dm3 = (K_plus_asPi_MissID + pi_minus + pi_plus).M() - (pi_plus + pi_minus).M();
   
            beta_K_plus = (-K_plus_P + pi_plus_P + pi_minus_P + Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);
            beta_pi_plus = (K_plus_P - pi_plus_P + pi_minus_P + Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);
            beta_pi_minus = (K_plus_P + pi_plus_P - pi_minus_P + Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);
            beta_Ds = (K_plus_P + pi_plus_P + pi_minus_P - Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);
            
            BsDTF_Kplus_PX = Bs_BsDTF_K_1_1270_plus_Kplus_PX[0] ;
            BsDTF_Kplus_PY = Bs_BsDTF_K_1_1270_plus_Kplus_PY[0] ;
            BsDTF_Kplus_PZ = Bs_BsDTF_K_1_1270_plus_Kplus_PZ[0] ;
            BsDTF_Kplus_PE = Bs_BsDTF_K_1_1270_plus_Kplus_PE[0] ;
            
            BsDTF_piplus_PX = Bs_BsDTF_K_1_1270_plus_piplus_0_PX[0] ;
            BsDTF_piplus_PY = Bs_BsDTF_K_1_1270_plus_piplus_0_PY[0] ;
            BsDTF_piplus_PZ = Bs_BsDTF_K_1_1270_plus_piplus_0_PZ[0] ;
            BsDTF_piplus_PE = Bs_BsDTF_K_1_1270_plus_piplus_0_PE[0] ;
            
            BsDTF_piminus_PX = Bs_BsDTF_K_1_1270_plus_piplus_PX[0] ;
            BsDTF_piminus_PY = Bs_BsDTF_K_1_1270_plus_piplus_PY[0] ;
            BsDTF_piminus_PZ = Bs_BsDTF_K_1_1270_plus_piplus_PZ[0] ;
            BsDTF_piminus_PE = Bs_BsDTF_K_1_1270_plus_piplus_PE[0] ;    

	    tv.push_back(K_plus);
	    tv.push_back(pi_plus);
	    tv.push_back(pi_minus);
        }
        else {
            TVector3 v_Ds(Ds_PX,Ds_PY,0.);
            TVector3 v_K(pi_plus1_PX,pi_plus1_PY,0.);
            TVector3 v_pip(pi_plus2_PX,pi_plus2_PY,0.);
            TVector3 v_pim(pi_minus_PX,pi_minus_PY,0.);
            angK= v_Ds.Angle(v_K);
            angPip= v_Ds.Angle(v_pip);
            angPim= v_Ds.Angle(v_pim);
            maxCos = cos(max(angK,max(angPip,angPim)));
            
            XsDaughters_min_IPCHI2 = min(pi_plus1_IPCHI2_OWNPV,min(pi_minus_IPCHI2_OWNPV,pi_plus2_IPCHI2_OWNPV));            
            XsDaughters_max_IPCHI2 = max(pi_plus1_IPCHI2_OWNPV,max(pi_minus_IPCHI2_OWNPV,pi_plus2_IPCHI2_OWNPV));            
            XsDaughters_min_PT = min(pi_plus1_PT,min(pi_minus_PT,pi_plus2_PT));          
            XsDaughters_min_P = min(pi_plus1_P,min(pi_minus_P,pi_plus2_P));          
	    Xs_max_DOCA = max(a_1_1260_plus_DOCA1,max(a_1_1260_plus_DOCA2,a_1_1260_plus_DOCA3));
            Xs_max_ghostProb = max(pi_plus1_TRACK_GhostProb,max(pi_plus2_TRACK_GhostProb,pi_minus_TRACK_GhostProb));   
            Xs_max_ProbNNghost = max(pi_plus1_ProbNNghost,max(pi_plus2_ProbNNghost,pi_minus_ProbNNghost));          
	    Xs_ptasy = a_1_1260_plus_ptasy_1_00;
	    Xs_PT = a_1_1260_plus_PT;
	    Xs_ETA = a_1_1260_plus_ETA;

            m_DsK = (BsDTF_Ds+BsDTF_pi_plus1).M();
            m_Dspi = (BsDTF_Ds+BsDTF_pi_plus2).M();
            m_Dspipi = (BsDTF_Ds+BsDTF_pi_plus1+BsDTF_pi_minus).M();
            m_DsKpi = (BsDTF_Ds+BsDTF_pi_plus2+BsDTF_pi_minus).M();
            m_DsKpip = (BsDTF_Ds+BsDTF_pi_plus1+BsDTF_pi_plus2).M(); 
            m_Kpipi = (BsDTF_pi_minus+BsDTF_pi_plus1+BsDTF_pi_plus2).M(); 
     	    m_pipi = (BsDTF_pi_minus+BsDTF_pi_plus1).M();
	    m_Kpi = (BsDTF_pi_minus+BsDTF_pi_plus2).M();

            TLorentzVector pi_plus1_asK_MissID; 
            TLorentzVector pi_plus2_asK_MissID; 

            pi_plus1_asK_MissID.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ, massKaon);
            pi_plus2_asK_MissID.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ, massKaon);

            if(pi_plus1_PIDK > pi_plus2_PIDK ){
                bkg_3pi_as_Xs_m = (pi_minus+pi_plus2+pi_plus1_asK_MissID).M();
                bkg_3pi_as_Xs_Bs_m = (Ds+pi_minus+pi_plus2+pi_plus1_asK_MissID).M();
            }else { 
                bkg_3pi_as_Xs_m = (pi_minus+pi_plus1+pi_plus2_asK_MissID).M();
                bkg_3pi_as_Xs_Bs_m = (Ds+pi_minus+pi_plus1+pi_plus2_asK_MissID).M();
            }
            beta_K_plus = (-pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
            beta_pi_plus = (pi_plus1_P - pi_plus2_P + pi_minus_P + Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
            beta_pi_minus = (pi_plus1_P + pi_plus2_P - pi_minus_P + Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
            beta_Ds = (pi_plus1_P + pi_plus2_P + pi_minus_P - Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
            
            BsDTF_Kplus_PX = Bs_BsDTF_a_1_1260_plus_piplus_1_PX[0] ;
            BsDTF_Kplus_PY = Bs_BsDTF_a_1_1260_plus_piplus_1_PY[0] ;
            BsDTF_Kplus_PZ = Bs_BsDTF_a_1_1260_plus_piplus_1_PZ[0] ;
            BsDTF_Kplus_PE = Bs_BsDTF_a_1_1260_plus_piplus_1_PE[0] ;
            
            BsDTF_piplus_PX = Bs_BsDTF_a_1_1260_plus_piplus_PX[0] ;
            BsDTF_piplus_PY = Bs_BsDTF_a_1_1260_plus_piplus_PY[0] ;
            BsDTF_piplus_PZ = Bs_BsDTF_a_1_1260_plus_piplus_PZ[0] ;
            BsDTF_piplus_PE = Bs_BsDTF_a_1_1260_plus_piplus_PE[0] ;
            
            BsDTF_piminus_PX = Bs_BsDTF_a_1_1260_plus_piplus_0_PX[0] ;
            BsDTF_piminus_PY = Bs_BsDTF_a_1_1260_plus_piplus_0_PY[0] ;
            BsDTF_piminus_PZ = Bs_BsDTF_a_1_1260_plus_piplus_0_PZ[0] ;
            BsDTF_piminus_PE = Bs_BsDTF_a_1_1260_plus_piplus_0_PE[0] ;    

	    tv.push_back(pi_plus1);
	    tv.push_back(pi_plus2);
	    tv.push_back(pi_minus);
        }
        
        if(_Ds_finalState == Ds_finalState::pipipi){
            DsDaughters_min_IPCHI2 = min(pi_plus_fromDs_IPCHI2_OWNPV,min(pi_minus_fromDs_IPCHI2_OWNPV,pi_minus2_fromDs_IPCHI2_OWNPV));            
            DsDaughters_max_IPCHI2 = max(pi_plus_fromDs_IPCHI2_OWNPV,max(pi_minus_fromDs_IPCHI2_OWNPV,pi_minus2_fromDs_IPCHI2_OWNPV));            
            DsDaughters_min_PT = min(pi_plus_fromDs_PT,min(pi_minus_fromDs_PT,pi_minus2_fromDs_PT));            
            DsDaughters_min_P = min(pi_plus_fromDs_P,min(pi_minus_fromDs_P,pi_minus2_fromDs_P));            
            Ds_max_DOCA = max(Ds_DOCA1,max(Ds_DOCA2,Ds_DOCA3));
            Ds_max_ghostProb = max(pi_plus_fromDs_TRACK_GhostProb,max(pi_minus2_fromDs_TRACK_GhostProb,pi_minus_fromDs_TRACK_GhostProb));  
            Ds_max_ProbNNghost = max(pi_plus_fromDs_ProbNNghost,max(pi_minus2_fromDs_ProbNNghost,pi_minus_fromDs_ProbNNghost));  
            
            Ds_m12 = (DTF_pi_plus_fromDs + DTF_pi_minus_fromDs).M(); 
            Ds_m13 = (DTF_pi_plus_fromDs + DTF_pi_minus2_fromDs).M(); 
            
            beta_K_minus_fromDs = (pi_plus_fromDs_P + pi_minus_fromDs_P - pi_minus2_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P);
            beta_K_plus_fromDs = (-pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P);
            beta_pi_minus_fromDs = (pi_plus_fromDs_P - pi_minus_fromDs_P + pi_minus2_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P);
            
            BsDTF_Ds_Kplus_PX = Bs_BsDTF_D_splus_piplus_PX[0] ;
            BsDTF_Ds_Kplus_PY = Bs_BsDTF_D_splus_piplus_PY[0] ;
            BsDTF_Ds_Kplus_PZ = Bs_BsDTF_D_splus_piplus_PZ[0] ;
            BsDTF_Ds_Kplus_PE = Bs_BsDTF_D_splus_piplus_PE[0] ;
            
            BsDTF_Ds_Kminus_PX = Bs_BsDTF_D_splus_piplus_0_PX[0] ;
            BsDTF_Ds_Kminus_PY = Bs_BsDTF_D_splus_piplus_0_PY[0] ;
            BsDTF_Ds_Kminus_PZ = Bs_BsDTF_D_splus_piplus_0_PZ[0] ;
            BsDTF_Ds_Kminus_PE = Bs_BsDTF_D_splus_piplus_0_PE[0] ;
            
            BsDTF_Ds_piminus_PX = Bs_BsDTF_D_splus_piplus_1_PX[0] ;
            BsDTF_Ds_piminus_PY = Bs_BsDTF_D_splus_piplus_1_PY[0] ;
            BsDTF_Ds_piminus_PZ = Bs_BsDTF_D_splus_piplus_1_PZ[0] ;
            BsDTF_Ds_piminus_PE = Bs_BsDTF_D_splus_piplus_1_PE[0] ;

	    tv.push_back(pi_plus_fromDs);
	    tv.push_back(pi_minus_fromDs);
	    tv.push_back(pi_minus2_fromDs);
        }
        
        else if(_Ds_finalState == Ds_finalState::Kpipi){
            DsDaughters_min_IPCHI2 = min(pi_plus_fromDs_IPCHI2_OWNPV,min(pi_minus_fromDs_IPCHI2_OWNPV,K_minus_fromDs_IPCHI2_OWNPV));            
            DsDaughters_max_IPCHI2 = max(pi_plus_fromDs_IPCHI2_OWNPV,max(pi_minus_fromDs_IPCHI2_OWNPV,K_minus_fromDs_IPCHI2_OWNPV));            
            DsDaughters_min_PT = min(pi_plus_fromDs_PT,min(pi_minus_fromDs_PT,K_minus_fromDs_PT));            
            DsDaughters_min_P = min(pi_plus_fromDs_P,min(pi_minus_fromDs_P,K_minus_fromDs_P));            
	    Ds_max_DOCA = max(Ds_DOCA1,max(Ds_DOCA2,Ds_DOCA3));
            Ds_max_ghostProb = max(pi_plus_fromDs_TRACK_GhostProb,max(K_minus_fromDs_TRACK_GhostProb,pi_minus_fromDs_TRACK_GhostProb));  
            Ds_max_ProbNNghost = max(pi_plus_fromDs_ProbNNghost,max(K_minus_fromDs_ProbNNghost,pi_minus_fromDs_ProbNNghost));  
            
            Ds_m12 = (DTF_pi_plus_fromDs + DTF_K_minus_fromDs).M(); 
            Ds_m13 = (DTF_pi_plus_fromDs + DTF_pi_minus_fromDs).M(); 
            
            bkg_D_as_Ds_m= (pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus;

            bkg_Dstar_as_Ds_dm1= (pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - (pi_plus_fromDs + pi_minus_fromDs ).M();
            bkg_Dstar_as_Ds_dm3= (pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - (pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID ).M();
            bkg_Dstar_as_Ds_dm4= (pi_plus_fromDs + K_minus_fromDs + pi_minus_fromDs).M() - (pi_plus_fromDs + K_minus_fromDs ).M();

            bkg_Ks_as_KK_fromDs_m = (pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID).M();  

            bkg_Lambdac_as_Ds_m=  (pi_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c;  
            
            if(_decay == Decay::signal){
                bkg_D_as_Ds_Bs_m= (K_plus+pi_plus+pi_minus+pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M();
                bkg_Lambdac_as_Ds_Bs_m=  (K_plus+pi_plus+pi_minus+pi_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M();  
            }else {
                bkg_D_as_Ds_Bs_m= (pi_plus1+pi_plus2+pi_minus+pi_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M();
                bkg_Lambdac_as_Ds_Bs_m=  (pi_plus1+pi_plus2+pi_minus+pi_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M();  
            }
            beta_K_minus_fromDs = (pi_plus_fromDs_P + pi_minus_fromDs_P - K_minus_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            beta_K_plus_fromDs = (-pi_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            beta_pi_minus_fromDs = (pi_plus_fromDs_P - pi_minus_fromDs_P + K_minus_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            
            BsDTF_Ds_Kplus_PX = Bs_BsDTF_D_splus_piplus_0_PX[0] ;
            BsDTF_Ds_Kplus_PY = Bs_BsDTF_D_splus_piplus_0_PY[0] ;
            BsDTF_Ds_Kplus_PZ = Bs_BsDTF_D_splus_piplus_0_PZ[0] ;
            BsDTF_Ds_Kplus_PE = Bs_BsDTF_D_splus_piplus_0_PE[0] ;
            
            BsDTF_Ds_Kminus_PX = Bs_BsDTF_D_splus_Kplus_PX[0] ;
            BsDTF_Ds_Kminus_PY = Bs_BsDTF_D_splus_Kplus_PY[0] ;
            BsDTF_Ds_Kminus_PZ = Bs_BsDTF_D_splus_Kplus_PZ[0] ;
            BsDTF_Ds_Kminus_PE = Bs_BsDTF_D_splus_Kplus_PE[0] ;
            
            BsDTF_Ds_piminus_PX = Bs_BsDTF_D_splus_piplus_PX[0] ;
            BsDTF_Ds_piminus_PY = Bs_BsDTF_D_splus_piplus_PY[0] ;
            BsDTF_Ds_piminus_PZ = Bs_BsDTF_D_splus_piplus_PZ[0] ;
            BsDTF_Ds_piminus_PE = Bs_BsDTF_D_splus_piplus_PE[0] ;

	    tv.push_back(pi_plus_fromDs);
	    tv.push_back(K_minus_fromDs);
	    tv.push_back(pi_minus_fromDs);
        }
        
        else {
            DsDaughters_min_IPCHI2 = min(K_plus_fromDs_IPCHI2_OWNPV,min(pi_minus_fromDs_IPCHI2_OWNPV,K_minus_fromDs_IPCHI2_OWNPV));            
            DsDaughters_max_IPCHI2 = max(K_plus_fromDs_IPCHI2_OWNPV,max(pi_minus_fromDs_IPCHI2_OWNPV,K_minus_fromDs_IPCHI2_OWNPV));            
            DsDaughters_min_PT = min(K_plus_fromDs_PT,min(pi_minus_fromDs_PT,K_minus_fromDs_PT));            
            DsDaughters_min_P = min(K_plus_fromDs_P,min(pi_minus_fromDs_P,K_minus_fromDs_P));            
	    Ds_max_DOCA = max(Ds_DOCA1,max(Ds_DOCA2,Ds_DOCA3));
            Ds_max_ghostProb = max(K_plus_fromDs_TRACK_GhostProb,max(K_minus_fromDs_TRACK_GhostProb,pi_minus_fromDs_TRACK_GhostProb));  
            Ds_max_ProbNNghost = max(K_plus_fromDs_ProbNNghost,max(K_minus_fromDs_ProbNNghost,pi_minus_fromDs_ProbNNghost));  
            
            Ds_m12 = (DTF_K_plus_fromDs + DTF_K_minus_fromDs).M(); 
            Ds_m13 = (DTF_K_plus_fromDs + DTF_pi_minus_fromDs).M(); 
            
            bkg_D_as_Ds_m= (K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus;

	    TLorentzVector K_plus_fromDs_asPi_MissID; 
            K_plus_fromDs_asPi_MissID.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ, massPion);
            bkg_Dstar_as_Ds_dm1= (K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - (K_plus_fromDs + pi_minus_fromDs ).M();
            bkg_Dstar_as_Ds_dm2= (K_plus_fromDs_asPi_MissID + K_minus_fromDs + pi_minus_fromDs).M() - (K_minus_fromDs + K_plus_fromDs_asPi_MissID ).M();
            bkg_Dstar_as_Ds_dm3= (K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - (K_plus_fromDs + Kminus_fromDs_asPiminus_MissID ).M();
            bkg_Dstar_as_Ds_dm4= (K_plus_fromDs + K_minus_fromDs + pi_minus_fromDs).M() - (K_plus_fromDs + K_minus_fromDs ).M();

            bkg_Ks_as_KK_fromDs_m = (K_plus_fromDs + Kminus_fromDs_asPiminus_MissID).M();  

            bkg_Lambdac_as_Ds_m=  (K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c;  
            
            if(_decay == Decay::signal){
                bkg_D_as_Ds_Bs_m= (K_plus+pi_plus+pi_minus+K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M();
                bkg_Lambdac_as_Ds_Bs_m=  (K_plus+pi_plus+pi_minus+K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M();  
            }else {
                bkg_D_as_Ds_Bs_m= (pi_plus1+pi_plus2+pi_minus+K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M();
                bkg_Lambdac_as_Ds_Bs_m=  (pi_plus1+pi_plus2+pi_minus+K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M();  
            }
            beta_K_minus_fromDs = (K_plus_fromDs_P + pi_minus_fromDs_P - K_minus_fromDs_P) / (K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            beta_K_plus_fromDs = (-K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P) / (K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            beta_pi_minus_fromDs = (K_plus_fromDs_P - pi_minus_fromDs_P + K_minus_fromDs_P) / (K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            
            BsDTF_Ds_Kplus_PX = Bs_BsDTF_D_splus_Kplus_PX[0] ;
            BsDTF_Ds_Kplus_PY = Bs_BsDTF_D_splus_Kplus_PY[0] ;
            BsDTF_Ds_Kplus_PZ = Bs_BsDTF_D_splus_Kplus_PZ[0] ;
            BsDTF_Ds_Kplus_PE = Bs_BsDTF_D_splus_Kplus_PE[0] ;
            
            BsDTF_Ds_Kminus_PX = Bs_BsDTF_D_splus_Kplus_0_PX[0] ;
            BsDTF_Ds_Kminus_PY = Bs_BsDTF_D_splus_Kplus_0_PY[0] ;
            BsDTF_Ds_Kminus_PZ = Bs_BsDTF_D_splus_Kplus_0_PZ[0] ;
            BsDTF_Ds_Kminus_PE = Bs_BsDTF_D_splus_Kplus_0_PE[0] ;
            
            BsDTF_Ds_piminus_PX = Bs_BsDTF_D_splus_piplus_PX[0] ;
            BsDTF_Ds_piminus_PY = Bs_BsDTF_D_splus_piplus_PY[0] ;
            BsDTF_Ds_piminus_PZ = Bs_BsDTF_D_splus_piplus_PZ[0] ;
            BsDTF_Ds_piminus_PE = Bs_BsDTF_D_splus_piplus_PE[0] ;

	    tv.push_back(K_plus_fromDs);
	    tv.push_back(K_minus_fromDs);
	    tv.push_back(pi_minus_fromDs);
        }
        
        max_ghostProb = max(Xs_max_ghostProb,Ds_max_ghostProb);
        max_ProbNNghost = max(Xs_max_ProbNNghost,Ds_max_ProbNNghost);
        track_min_PT = min(XsDaughters_min_PT,DsDaughters_min_PT);
        track_min_P = min(XsDaughters_min_P,DsDaughters_min_P);  
        track_min_IPCHI2 = min(XsDaughters_min_IPCHI2,DsDaughters_min_IPCHI2);

	NTracks = (double)(nTracks);
	if(Bs_L0Global_TIS)TriggerCat = 2;
	else TriggerCat = 3;
        if(Bs_L0HadronDecision_TOS)TriggerCat = 0;
	else TriggerCat = 1;
	
        PV_status = Bs_PV_status[0];
        DTF_status = Bs_DTF_status[0];
        BsDTF_status = Bs_BsDTF_status[0];
        PV_CHI2NDOF = Bs_PV_chi2[0]/Bs_PV_nDOF[0];
        DTF_CHI2NDOF = Bs_DTF_chi2[0]/Bs_DTF_nDOF[0];
        BsDTF_CHI2NDOF = Bs_BsDTF_chi2[0]/Bs_BsDTF_nDOF[0];
 
        Bs_PV_MM = Bs_PV_M[0];
        Bs_DTF_MM = Bs_DTF_M[0];
        Ds_PV_MM = Bs_PV_Dplus_M[0];
        Bs_PV_MMERR = Bs_PV_MERR[0];
        Bs_DTF_MMERR = Bs_DTF_MERR[0];
        Ds_PV_MMERR = Bs_PV_Dplus_MERR[0];

        Bs_TAU = Bs_TAU*1000.;
        Bs_TAUERR = Bs_TAUERR*1000.;
        Ds_TAU = Ds_TAU*1000.;
        Ds_TAUERR = Ds_TAUERR*1000.;
        
        Bs_PV_TAU = Bs_PV_ctau[0] * 3.33564095; // 1e-3*1/(2.99792458*1e8)*1e12
        Bs_PV_TAUERR = Bs_PV_ctauErr[0] * 3.33564095;
        Bs_DTF_TAU = Bs_DTF_ctau[0] * 3.33564095;
        Bs_DTF_TAUERR = Bs_DTF_ctauErr[0] * 3.33564095;
        Bs_BsDTF_TAU = Bs_BsDTF_ctau[0] * 3.33564095;
        Bs_BsDTF_TAUERR = Bs_BsDTF_ctauErr[0] * 3.33564095;

        BsDTF_Ds_PX = BsDTF_Ds.Px() ;
        BsDTF_Ds_PY = BsDTF_Ds.Py() ;
        BsDTF_Ds_PZ = BsDTF_Ds.Pz() ;
        BsDTF_Ds_PE = BsDTF_Ds.E() ;

	if(run ==1){
		OS_Muon_DEC = (int) Bs_OS_Muon_DEC;
		OS_Muon_PROB = (double) Bs_OS_Muon_PROB;
		OS_Electron_DEC = (int) Bs_OS_Electron_DEC;
		OS_Electron_PROB= (double)Bs_OS_Electron_PROB;
		OS_Kaon_DEC = (int)Bs_OS_Kaon_DEC;
		OS_Kaon_PROB= (double)Bs_OS_Kaon_PROB;
		OS_VtxCharge_DEC = (int)Bs_VtxCharge_DEC;
		OS_VtxCharge_PROB = (double)Bs_VtxCharge_PROB;
		OS_nnetKaon_DEC= (int)Bs_OS_nnetKaon_DEC;
		OS_nnetKaon_PROB= (double)Bs_OS_nnetKaon_PROB;
		SS_Kaon_DEC= (int)Bs_SS_nnetKaon_DEC;
		SS_Kaon_PROB = (double)Bs_SS_nnetKaon_PROB;
		OS_Charm_DEC= (int)Bs_OS_Charm_DEC;
		OS_Charm_PROB = (double)Bs_OS_Charm_PROB;
	}
	else {
		OS_Muon_DEC = Bs_OS_Muon_TAGDEC;
		OS_Muon_PROB = Bs_OS_Muon_TAGETA;
		OS_Electron_DEC = Bs_OS_Electron_TAGDEC;
		OS_Electron_PROB= Bs_OS_Electron_TAGETA;
		OS_Kaon_DEC = Bs_OS_Kaon_TAGDEC;
		OS_Kaon_PROB= Bs_OS_Kaon_TAGETA;
		OS_VtxCharge_DEC = Bs_VtxCharge_TAGDEC;
		OS_VtxCharge_PROB = Bs_VtxCharge_TAGETA;
		OS_Kaon_DEC= Bs_OS_Kaon_TAGDEC;
		OS_Kaon_PROB=  Bs_OS_Kaon_TAGETA;
		SS_Kaon_DEC= Bs_SS_nnetKaon_TAGDEC;
		SS_Kaon_PROB = Bs_SS_nnetKaon_TAGETA;
		OS_Charm_DEC= Bs_OS_Charm_TAGDEC;
		OS_Charm_PROB = Bs_OS_Charm_TAGETA;
        
        OS_Muon_DEC_Run1 = Bs_OS_Muon_TAGDEC_Run1;
        OS_Muon_PROB_Run1 = Bs_OS_Muon_TAGETA_Run1;
        OS_Electron_DEC_Run1 = Bs_OS_Electron_TAGDEC_Run1;
        OS_Electron_PROB_Run1= Bs_OS_Electron_TAGETA_Run1;
        OS_Kaon_DEC_Run1= Bs_OS_Kaon_TAGDEC_Run1;
        OS_Kaon_PROB_Run1=  Bs_OS_Kaon_TAGETA_Run1;
        SS_Kaon_DEC_Run1= Bs_SS_nnetKaon_TAGDEC_Run1;
        SS_Kaon_PROB_Run1 = Bs_SS_nnetKaon_TAGETA_Run1;
	}

	//if(track_min_P < 2500) continue;
	if(!_ltu){
 		//if(Bs_BsDTF_TAU < 0.4 || Bs_BsDTF_TAU > 10.) continue;
		if(Bs_BsDTF_TAUERR < 0. || Bs_BsDTF_TAUERR > 0.1) continue;
	}
	else {
		//if(PV_CHI2NDOF > 10)continue;
		//if(max_ghostProb > 0.35) continue;
		if(track_min_PT < 500)continue;
// 		if(DsDaughters_min_IPCHI2 < 9)continue;
	}

	m_15 = (tv[1] + tv[5]).M();
        m_16 = (tv[1] + tv[6]).M();
	m_25 = (tv[2] + tv[5]).M();
	m_26 = (tv[2] + tv[6]).M();
	m_34 = (tv[3] + tv[4]).M();

	m_135 = (tv[1] + tv[3] + tv[5]).M();
	m_156 = (tv[1] + tv[5] + tv[6]).M();
	m_345 = (tv[3] + tv[4] + tv[5]).M();
	m_125 = (tv[1] + tv[2] + tv[5]).M();
	m_245 = (tv[2] + tv[4] + tv[5]).M();
	m_236 = (tv[2] + tv[3] + tv[6]).M();
	m_145 = (tv[1] + tv[4] + tv[5]).M();
	m_126 = (tv[1] + tv[2] + tv[6]).M();
	m_234 = (tv[2] + tv[3] + tv[4]).M();
	m_246 = (tv[2] + tv[4] + tv[6]).M();
	m_235 = (tv[2] + tv[3] + tv[5]).M();
	m_256 = (tv[2] + tv[5] + tv[6]).M();
	m_136 = (tv[1] + tv[3] + tv[6]).M();
	m_346 = (tv[3] + tv[4] + tv[6]).M();

	m_1236 = (tv[1] + tv[2] + tv[3]+ tv[6]).M();
	m_2346 = (tv[2] + tv[3] + tv[4]+ tv[6]).M();
	m_1235 = (tv[1] + tv[2] + tv[3]+ tv[5]).M();
	m_1256 = (tv[1] + tv[2] + tv[5]+ tv[6]).M();
	m_2345 = (tv[2] + tv[3] + tv[4]+ tv[5]).M();
	m_2456 = (tv[2] + tv[4] + tv[5]+ tv[6]).M();

        // Fill tree
        summary_tree->Fill();
        if(0ul == (i % 100000ul))summary_tree->AutoSave();
    }
   
    cout << "Selected " << summary_tree->GetEntries() << " events" <<  endl;
    cout << "Efficiency = " << summary_tree->GetEntries()/(double)nentries * 100. << " %" <<  endl;
    
    summary_tree->Write();
    output->Close();
    
    cout << "Created new file " << _outFileName << endl;
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;

}
