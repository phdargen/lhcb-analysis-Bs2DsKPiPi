#define pull_cxx
#include "pull.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

TString pull::latexName(TString s){
    if(s == "D_bar") return "$\\bar{D}$";
    if(s == "S_bar") return "$\\bar{S}$";

    if(s == "p0_os_Run1") return "$p_{0}^{OS} \\, \\text{Run-I}$";
    if(s == "p1_os_Run1") return "$p_{1}^{OS} \\, \\text{Run-I}$";
    if(s == "delta_p0_os_Run1") return "$\\Delta p_{0}^{OS} \\, \\text{Run-I}$";
    if(s == "delta_p1_os_Run1") return "$\\Delta p_{1}^{OS} \\, \\text{Run-I}$";
    if(s == "tageff_os_Run1") return "$\\epsilon_{\\text{tag}}^{OS} \\, \\text{Run-I}$";
    if(s == "tageff_asym_os_Run1") return "$\\Delta \\epsilon_{\\text{tag}}^{OS} \\, \\text{Run-I}$";

    if(s == "p0_ss_Run1") return "$p_{0}^{SS} \\, \\text{Run-I}$";
    if(s == "p1_ss_Run1") return "$p_{1}^{SS} \\, \\text{Run-I}$";
    if(s == "delta_p0_ss_Run1") return "$\\Delta p_{0}^{SS} \\, \\text{Run-I}$";
    if(s == "delta_p1_ss_Run1") return "$\\Delta p_{1}^{SS} \\, \\text{Run-I}$";
    if(s == "tageff_ss_Run1") return "$\\epsilon_{\\text{tag}}^{SS} \\, \\text{Run-I}$";
    if(s == "tageff_asym_ss_Run1") return "$\\Delta \\epsilon_{\\text{tag}}^{SS} \\, \\text{Run-I}$";

    if(s == "p0_os_Run2") return "$p_{0}^{OS} \\, \\text{Run-II}$";
    if(s == "p1_os_Run2") return "$p_{1}^{OS} \\, \\text{Run-II}$";
    if(s == "delta_p0_os_Run2") return "$\\Delta p_{0}^{OS} \\, \\text{Run-II}$";
    if(s == "delta_p1_os_Run2") return "$\\Delta p_{1}^{OS} \\, \\text{Run-II}$";
    if(s == "tageff_os_Run2") return "$\\epsilon_{\\text{tag}}^{OS} \\, \\text{Run-II}$";
    if(s == "tageff_asym_os_Run2") return "$\\Delta \\epsilon_{\\text{tag}}^{OS} \\, \\text{Run-II}$";

    if(s == "p0_ss_Run2") return "$p_{0}^{SS} \\, \\text{Run-II}$";
    if(s == "p1_ss_Run2") return "$p_{1}^{SS} \\, \\text{Run-II}$";
    if(s == "delta_p0_ss_Run2") return "$\\Delta p_{0}^{SS} \\, \\text{Run-II}$";
    if(s == "delta_p1_ss_Run2") return "$\\Delta p_{1}^{SS} \\, \\text{Run-II}$";
    if(s == "tageff_ss_Run2") return "$\\epsilon_{\\text{tag}}^{SS} \\, \\text{Run-II}$";
    if(s == "tageff_asym_ss_Run2") return "$\\Delta \\epsilon_{\\text{tag}}^{SS} \\, \\text{Run-II}$";

    if(s == "production_asym_Run2") return "$A_{P} \\, \\text{Run-II}$";

    if(s == "dm") return "$\\Delta m_{s}$";

    return "$" + s + "$";
}


TString pull::latexNameMod(TString s){
    if(s == "p0_os_Run1") return "$p_{0}^{OS}$";
    if(s == "p1_os_Run1") return "$p_{1}^{OS}$";
    if(s == "delta_p0_os_Run1") return "$\\Delta p_{0}^{OS}$";
    if(s == "delta_p1_os_Run1") return "$\\Delta p_{1}^{OS}$";
    if(s == "tageff_os_Run1") return "$\\epsilon_{\\text{tag}}^{OS} \\, [\\%]$";
    if(s == "tageff_asym_os_Run1") return "$\\Delta \\epsilon_{\\text{tag}}^{OS} \\, [\\%]$";

    if(s == "p0_ss_Run1") return "$p_{0}^{SS}$";
    if(s == "p1_ss_Run1") return "$p_{1}^{SS}$";
    if(s == "delta_p0_ss_Run1") return "$\\Delta p_{0}^{SS}$";
    if(s == "delta_p1_ss_Run1") return "$\\Delta p_{1}^{SS}$";
    if(s == "tageff_ss_Run1") return "$\\epsilon_{\\text{tag}}^{SS} \\, [\\%]$";
    if(s == "tageff_asym_ss_Run1") return "$\\Delta \\epsilon_{\\text{tag}}^{SS} \\, [\\%]$";

    if(s == "p0_os_Run2") return "$p_{0}^{OS}$";
    if(s == "p1_os_Run2") return "$p_{1}^{OS}$";
    if(s == "delta_p0_os_Run2") return "$\\Delta p_{0}^{OS}$";
    if(s == "delta_p1_os_Run2") return "$\\Delta p_{1}^{OS}$";
    if(s == "tageff_os_Run2") return "$\\epsilon_{\\text{tag}}^{OS} \\, [\\%]$";
    if(s == "tageff_asym_os_Run2") return "$\\Delta \\epsilon_{\\text{tag}}^{OS} \\, [\\%]$";

    if(s == "p0_ss_Run2") return "$p_{0}^{SS}$";
    if(s == "p1_ss_Run2") return "$p_{1}^{SS}$";
    if(s == "delta_p0_ss_Run2") return "$\\Delta p_{0}^{SS}$";
    if(s == "delta_p1_ss_Run2") return "$\\Delta p_{1}^{SS}$";
    if(s == "tageff_ss_Run2") return "$\\epsilon_{\\text{tag}}^{SS} \\, [\\%]$";
    if(s == "tageff_asym_ss_Run2") return "$\\Delta \\epsilon_{\\text{tag}}^{SS}  \\, [\\%]$";
    if(s == "production_asym_Run2") return "$A_{P} \\, [\\%]$";

    if(s == "dm") return "$\\Delta m_{s} \\, [\\text{ps}^{-1}]$";

    return latexName(s);
}


TMatrixD pull::combineCov_maxVal(vector<TMatrixD*> vec){

	TMatrixD m(*vec[0]);

	for (int i = 0 ; i <  m.GetNcols(); i++) {
		for (int j = 0 ; j < m.GetNcols(); j++) {
			double max = 0.;
			for(int n = 0; n < vec.size(); n++){
				 if(abs((*vec[n])(i,j)) > abs(max)) max = (*vec[n])(i,j);
			}
			m(i,j) = max; 
		}
	}

	return m;
}

TMatrixD pull::sampleVariance(vector< vector<double> > vec_vals){

	TMatrixD m(vec_vals[0].size(),vec_vals[0].size());

	for(int i=0; i < vec_vals[0].size(); i++){
			double sum = 0.;
			double sum_square = 0.;
			for(int n = 0; n < vec_vals.size(); n++){
				sum += vec_vals[n][i];
				sum_square += pow(vec_vals[n][i],2);	
			}
			m(i,i) = (sum_square - pow(sum,2)/vec_vals.size())/(vec_vals.size()-1.);		
	}

	return m;
}


vector <double> pull::getVals(){
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }
    fChain->GetEntry(0);  
    vector<double> vals;
    for (int i = 0 ; i < _paraNames.size(); i++) 
	vals.push_back(*_means[i]);
    
    return vals;
}   

vector <double> pull::getErrs(){
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }
    fChain->GetEntry(0);  
    vector<double> vals;
    for (int i = 0 ; i < _paraNames.size(); i++) 
	vals.push_back(*_errs[i]);
    
    return vals;
}   


TMatrixD pull::getStatCov(TString label){
    
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }
    int N = fChain->GetEntries();
           
    TMatrixD cov(_paraNames.size(),_paraNames.size());
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
	for (int i = 0 ; i < _paraNames.size(); i++) {
		for (int j = 0 ; j < _paraNames.size(); j++) {
			cov[i][j] += (*_means[i] - *_inits[i]) * (*_means[j] - *_inits[j])/(N-1.);
		}
	}
    }

    return cov;
}

TMatrixD pull::getCov(TString label){
    
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }

    TMatrixD cov(_paraNames.size(),_paraNames.size());
    int N = fChain->GetEntries();

    if(N == 1){
        fChain->GetEntry(0);  
        for (int i = 0 ; i < _paraNames.size(); i++) 
            for (int j = 0 ; j < _paraNames.size(); j++) {
                if(i==j)cov[i][j] = pow((*_means[i] - *_inits[i]),2);
                else cov[i][j] = 0.;
            }
        return cov;
    }
    
    vector<TH1D*> h_pulls;
    for (int i = 0 ; i < _paraNames.size(); i++) 
        h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 40, -5.,5.));
    
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        for (int i = 0 ; i < _paraNames.size(); i++){
            h_pulls[i]->Fill(*_pulls[i]);
            for (int j = 0 ; j < _paraNames.size(); j++) {
                cov[i][j] += (*_means[i] - *_inits[i]) * (*_means[j] - *_inits[j])/(N-1.);
            }
        }
    }
    
    TCanvas* c = new TCanvas();
    TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
    gaussian->SetParameters(1.,0.,1.);
    gaussian->SetParLimits(1,-1., 1.);
    gaussian->SetParLimits(2, 0., 2.);
    gaussian->SetLineColor(kRed);
        
    ofstream SummaryFile;
    SummaryFile.open("pull_results/pull_table"+label+".tex",std::ofstream::trunc);
    SummaryFile << "\\begin{tabular}{l  c  c}" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    
    vector<double> fit_means,fit_sigmas;
    for(int i = 0 ; i < _paraNames.size(); i++) {
        h_pulls[i]->Fit(gaussian);
        fit_means.push_back(gaussian->GetParameter(1));
        fit_sigmas.push_back(gaussian->GetParameter(2));

        SummaryFile << std::fixed << std::setprecision(2) << latexName(_paraNames[i]) << " & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
        h_pulls[i]->Draw("");
        gaussian->Draw("SAME");
        c->Print("pull_results/pull_"+ _paraNames[i] + label + ".eps");
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";

    TMatrixD cov_prime(cov);
    for (int i = 0 ; i < _paraNames.size(); i++)
        for (int j = 0 ; j < _paraNames.size(); j++) 
            cov_prime[i][j] = cov[i][j]*abs(fit_means[i])*abs(fit_means[j]);
    
    return cov_prime;
}

TMatrixD pull::getDeltaCov(TString refFileName,TString label){
    
    int N = fChain->GetEntries();

    TChain* chain =  new TChain("MinuitParameterSetNtp");
    if(N>1)for(int i = 1; i <= N; i++){
	stringstream index;
	index << i;
	TString file = refFileName;
	chain->Add(file.ReplaceAll("*",index.str())); 
    }    
    else chain->Add(refFileName); 

    if(fChain == 0 || chain == 0){
        cout << "ERROR:: No file found" << endl;
        throw "ERROR";
    }

    if(N > chain->GetEntries()){
        cout << "ERROR:: Inconsistent number of entries" << endl;
        throw "ERROR";
    }
    
    vector<double*> means;
    vector<double*> inits;
    vector<double*> errs;
    vector<double*> pulls;
    
    for (int i = 0 ; i < _paraNames.size(); i++) {
        double * mean = new double[1];
        means.push_back(mean);
        chain->SetBranchAddress(_paraNames[i]+"_mean", mean);
        
        double * init = new double[1];
        inits.push_back(init);
        chain->SetBranchAddress(_paraNames[i]+"_init", init);
        
        double * err = new double[1];
        errs.push_back(err);
        chain->SetBranchAddress(_paraNames[i]+"_err", err);
        
        double * pull = new double[1];
        pulls.push_back(pull);
        chain->SetBranchAddress(_paraNames[i]+"_pull", pull);
    } 
    
    TMatrixD cov(_paraNames.size(),_paraNames.size());

    if(N == 1){
        fChain->GetEntry(0); 
        chain->GetEntry(0); 
        for (int i = 0 ; i < _paraNames.size(); i++) 
            for (int j = 0 ; j < _paraNames.size(); j++) {
                if(i==j)cov[i][j] = pow((*_means[i] - *means[i]),2);
                else cov[i][j] = 0.;
            }
        return cov;
    }

    double max = 0.;
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        chain->GetEntry(n);
        for (int i = 0 ; i < _paraNames.size(); i++){
	   double val = (*_means[i]-*means[i])/(*errs[i]);
           if(abs(val) > max) max = abs(val);
	}
    }
    max *= 1.1;
    
    vector<TH1D*> h_pulls;
    for (int i = 0 ; i < _paraNames.size(); i++) 
        h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 25, -max,max));
    
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        chain->GetEntry(n);
        for (int i = 0 ; i < _paraNames.size(); i++)
            h_pulls[i]->Fill((*_means[i]-*means[i])/(*errs[i]));
    }
    
    TCanvas* c = new TCanvas();
    TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
    gaussian->SetParameters(1.,0.,1.);
    gaussian->SetParLimits(1,-1., 1.);
    gaussian->SetParLimits(2, 0., 2.);
    gaussian->SetLineColor(kRed);
    
    ofstream SummaryFile;
    SummaryFile.open("pull_results/pull_table"+label+".tex",std::ofstream::trunc);
    SummaryFile << "\\begin{tabular}{l  c  c}" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    
    vector<double> fit_means,fit_sigmas;
    for (int i = 0 ; i < _paraNames.size(); i++) {
        h_pulls[i]->Fit(gaussian);
        fit_means.push_back(gaussian->GetParameter(1));
        fit_sigmas.push_back(gaussian->GetParameter(2));
        SummaryFile << std::fixed << std::setprecision(2) << latexName(_paraNames[i]) << " & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
        h_pulls[i]->Draw("");
        gaussian->Draw("SAME");
        c->Print("pull_results/pull_"+ _paraNames[i] + label + ".eps");
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";
    
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        chain->GetEntry(n);
        for (int i = 0 ; i < _paraNames.size(); i++)
            for (int j = 0 ; j < _paraNames.size(); j++) 
                cov[i][j] += (*_means[i] - *means[i]) * (*_means[j] - *means[j])/(N-1.);
    }
    
    TMatrixD cov_prime(cov);
    for (int i = 0 ; i < _paraNames.size(); i++)
        for (int j = 0 ; j < _paraNames.size(); j++) 
            cov_prime[i][j] = cov[i][j]*sqrt(pow(fit_means[i],2)+pow(fit_sigmas[i],2))*sqrt(pow(fit_means[j],2)+pow(fit_sigmas[j],2))/fit_sigmas[i]/fit_sigmas[j];
    
    return cov_prime;
}

TMatrixD pull::getDeltaCovChol(TString refFileName,TString label,int varPerParChol){
    
    int N = fChain->GetEntries();

    TChain* chain =  new TChain("MinuitParameterSetNtp");
    if(N>1)for(int i = 1; i <= N; i++){
	stringstream index;
	index << i;
	TString file = refFileName;
	chain->Add(file.ReplaceAll("*",index.str())); 
    }    
    else chain->Add(refFileName); 

    if(fChain == 0 || chain == 0){
        cout << "ERROR:: No file found" << endl;
        throw "ERROR";
    }

    if(N > chain->GetEntries()){
        cout << "ERROR:: Inconsistent number of entries" << endl;
        throw "ERROR";
    }
    
    vector<double*> means;
    vector<double*> inits;
    vector<double*> errs;
    vector<double*> pulls;
    
    for (int i = 0 ; i < _paraNames.size(); i++) {
        double * mean = new double[1];
        means.push_back(mean);
        chain->SetBranchAddress(_paraNames[i]+"_mean", mean);
        
        double * init = new double[1];
        inits.push_back(init);
        chain->SetBranchAddress(_paraNames[i]+"_init", init);
        
        double * err = new double[1];
        errs.push_back(err);
        chain->SetBranchAddress(_paraNames[i]+"_err", err);
        
        double * pull = new double[1];
        pulls.push_back(pull);
        chain->SetBranchAddress(_paraNames[i]+"_pull", pull);
    } 
    
    int N_chol = N/varPerParChol;
    
    TMatrixD cov_tot(_paraNames.size(),_paraNames.size());

    for(int ic = 0 ; ic < N_chol; ic ++ ){
    
        TMatrixD cov(_paraNames.size(),_paraNames.size());

        vector<TH1D*> h_pulls;
        for (int i = 0 ; i < _paraNames.size(); i++) 
            h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 40, -1.,1.));
        
        for (int n= ic*varPerParChol; n < (ic+1)*varPerParChol ;n++) {
            fChain->GetEntry(n);  
            chain->GetEntry(n);
            for (int i = 0 ; i < _paraNames.size(); i++)
                h_pulls[i]->Fill((*_means[i]-*means[i])/(*errs[i]));
        }
        
        TCanvas* c = new TCanvas();
        TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
        gaussian->SetParameters(1.,0.,1.);
        gaussian->SetParLimits(1,-1., 1.);
        gaussian->SetParLimits(2, 0., 2.);
        gaussian->SetLineColor(kRed);
        
        ofstream SummaryFile;
	stringstream number;
	number << ic;
        SummaryFile.open("pull_results/pull_table"+label+"_par_"+number.str()+".tex",std::ofstream::trunc);
        SummaryFile << "\\begin{tabular}{l  c  c}" << "\n";
        SummaryFile << "\\hline" << "\n";
        SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
        SummaryFile << "\\hline" << "\n";
        SummaryFile << "\\hline" << "\n";
        
        vector<double> fit_means,fit_sigmas;
        
        for (int i = 0 ; i < _paraNames.size(); i++) {
            h_pulls[i]->Fit(gaussian);
            fit_means.push_back(gaussian->GetParameter(1));
            fit_sigmas.push_back(gaussian->GetParameter(2));
            SummaryFile << std::fixed << std::setprecision(2) << latexName(_paraNames[i]) << " & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
            h_pulls[i]->Draw("");
            gaussian->Draw("SAME");
            c->Print("pull_results/pull_"+ _paraNames[i] + label + "_par_"+number.str()+ ".eps");
        }
        
        SummaryFile << "\\hline" << "\n";
        SummaryFile << "\\end{tabular}" << "\n";
        
        for (int n= ic*varPerParChol; n < (ic+1)*varPerParChol ;n++) {
            fChain->GetEntry(n);  
            chain->GetEntry(n);
            for (int i = 0 ; i < _paraNames.size(); i++)
                for (int j = 0 ; j < _paraNames.size(); j++) 
                    cov[i][j] += (*_means[i] - *means[i]) * (*_means[j] - *means[j])/(N-1.);
        }
        
        TMatrixD cov_prime(cov);
        for (int i = 0 ; i < _paraNames.size(); i++)
            for (int j = 0 ; j < _paraNames.size(); j++) 
                cov_prime[i][j] = cov[i][j]*sqrt(pow(fit_means[i],2)+pow(fit_sigmas[i],2))*sqrt(pow(fit_means[j],2)+pow(fit_sigmas[j],2))/fit_sigmas[i]/fit_sigmas[j];
        
        cov_tot += cov_prime;
    }    
        
    return cov_tot;
}
