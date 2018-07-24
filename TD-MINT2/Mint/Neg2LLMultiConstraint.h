#ifndef NEG_TWO_LL_MultiConstraint_HH
#define NEG_TWO_LL_MultiConstraint_HH
// author: Philippe d'Argent (p.dargent@cern.ch)

#include "TMath.h"
#include "Mint/Minimisable.h"
#include "Mint/Neg2LLConstraint.h"
#include "Mint/NamedParameter.h"
#include "RooMultiVarGaussian.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooDataSet.h"
#include "TDecompChol.h"
#include <TROOT.h>
#include "TRandom3.h"
#include <vector>

using namespace std;
using namespace RooFit ;

namespace MINT{

class Neg2LLMultiConstraint: public Minimisable{
    protected:
    RooArgList* _x, *_mu;
    RooMultiVarGaussian* _gauss_cov;
    TMatrixTSym<double>* _cov;	
    MinuitParameterSet* _mps;
    TMatrixD* _UT;

    public:
        Neg2LLMultiConstraint(MinuitParameterSet* mps=0, string label = "")
        : Minimisable(mps), _mps(mps) {

			/// Parameters to constrain
			NamedParameter<string> constrain(("ConstrainMulti"+label).c_str(), (string)"" );
			string constrain_s = constrain;			
			vector<string> constrain_v = split(constrain_s, " ");

			/// Cov. matrix			
			NamedParameter<string> corr(("ConstrainMulti"+label+"_corr").c_str(),(string)"");
			string corr_s = corr;			
			vector<double> corr_v = splitValues(corr_s, " ");
			double* corr_a = &corr_v[0];
			TMatrixDSym cov(constrain_v.size(),corr_a);
			//cov.Print();

			/// Sanity checks
			if(constrain_v.size() != sqrt(corr_v.size())) { 
				cout << "ERROR in Neg2LLMultiConstraint::Inconsistent number of parameters and cov. matrix" << endl;
				throw "ERROR";			
			}

			_x = new RooArgList();
			_mu = new RooArgList();
			vector<double> sigma_v;

			for(int j= 0; j < constrain_v.size();j++){

				cout << "Adding gauss constraint for parameter: " << constrain_v[j] << endl;

				double mean = -99999;
				double sigma = -99999;

				for(unsigned int i=0; i < _mps->size(); i++){
				        if(0 == _mps->getParPtr(i)) continue;
					if(constrain_v[j] == _mps->getParPtr(i)->name()){	
						mean = _mps->getParPtr(i)->mean();
						sigma = ((FitParameter*)_mps->getParPtr(i))->stepInit();
					}
				}
				if(mean == -99999) {
					cout << "ERROR in Neg2LLMultiConstraint::Parameter not found" << endl;
					throw "ERROR";			
				}
				else {
 					RooRealVar* x = new RooRealVar((constrain_v[j]).c_str(), (constrain_v[j]).c_str(),mean);
 					_x->add(*x);
 					_mu->add(RooRealConstant::value(mean));
					sigma_v.push_back(sigma);
				}
			}
				
			for(int i=0; i < cov.GetNcols(); i++){
				for(int j=0; j < cov.GetNcols(); j++){    
					cov(i,j) = cov(i,j) * sigma_v[i] * sigma_v[j];
				}
			}
			_cov = new TMatrixDSym(cov);

			_x->Print();
			_mu->Print();
			_cov->Print();
			_gauss_cov = new RooMultiVarGaussian("gauss_cov","gauss_cov",*_x, *_mu, *_cov);

			TDecompChol tdc(*_cov);
			tdc.Decompose();
			TMatrixD U = tdc.GetU();
			_UT = new TMatrixD(TMatrixD::kTransposed,U);
			_UT->Print();

	};
                
        virtual void beginFit(){};
        virtual void parametersChanged(){
		for(int i=0; i < _x->getSize(); i++){
			((RooRealVar*)_x->at(i))->setVal(_mps->getParPtr(_x->at(i)->GetName())->mean());
		}
	};
        virtual void endFit(){};
        
        virtual double getVal(){
		return -2. * log(_gauss_cov->getVal(*_x));
        }
    
        virtual double getNewVal(){ 
            parametersChanged();
            return getVal();
        }
    
	int getNumberParams(){
		return _cov->GetNcols();
	}

	void smearInputValues(){
		
			cout << "Smearing input values " << endl;

			RooDataSet* data_cov = _gauss_cov->generate(*_x, 1);
			RooArgList* xvec_cov= (RooArgList*)data_cov->get(0);

			xvec_cov->Print();
		
			for(int i=0; i < xvec_cov->getSize(); i++){
				_mps->getParPtr(xvec_cov->at(i)->GetName())->setCurrentFitVal(((RooRealVar*)xvec_cov->at(i))->getVal());
				((FitParameter*)_mps->getParPtr(xvec_cov->at(i)->GetName()))->setInit(((RooRealVar*)xvec_cov->at(i))->getVal());
				cout << "Set parameter " << xvec_cov->at(i)->GetName() << " to " << _mps->getParPtr(xvec_cov->at(i)->GetName())->mean() << endl;   
			}
	}

	double smearInputValuesChol(int index = 0){
		
			cout << "Smearing input values for parameter " << index << endl;
			RooDataSet* data_cov = _gauss_cov->generate(*_x, 1);
			RooArgList* xvec_cov= (RooArgList*)data_cov->get(0);
 
			double val = 0;
			for(int i = 0 ; i < _UT->GetNcols(); i++){
				if(i != index)continue;
				val = _mps->getParPtr(xvec_cov->at(i)->GetName())->mean();
				for(int j = 0 ; j < _UT->GetNcols(); j++){
					// 			cout << UT(i,j) << " " ;
					// 			if(j == UT.GetNcols()-1) cout << endl;
					val += gRandom->Gaus(0.,1.) * (*_UT)(i,j);
				}	

				_mps->getParPtr(xvec_cov->at(i)->GetName())->setCurrentFitVal(val);
				((FitParameter*)_mps->getParPtr(xvec_cov->at(i)->GetName()))->setInit(val);
				cout << "Set parameter " << xvec_cov->at(i)->GetName() << " to " << val << endl;
			}
			return val;
	}

	RooDataSet* generateToys(int N = 100){
			cout << "Smearing input values " << endl;
			RooDataSet* data_cov = _gauss_cov->generate(*_x, N);
			return data_cov;
	}

	TMatrixDSym* getCovMatrix(){
		return _cov;
	}

	vector<string> split(string s, string delimiter){
		size_t pos = 0;
		string token;
		vector<string> list;
		while ((pos = s.find(delimiter)) != string::npos) {
			token = s.substr(0, pos);
			list.push_back(token);
			s.erase(0, pos + delimiter.length());
		}
		return list;
	}

	vector<double> splitValues(string s, string delimiter){
		size_t pos = 0;
		string token;
		vector<double> list;
		while ((pos = s.find(delimiter)) != string::npos) {
			token = s.substr(0, pos);
			list.push_back(atof(token.c_str()));
			s.erase(0, pos + delimiter.length());
		}
		return list;
	}


        virtual ~Neg2LLMultiConstraint(){}
        
};

}// namespace MINT
#endif
//
