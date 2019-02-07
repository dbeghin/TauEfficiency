#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TStyle.h"


using namespace std;


TH1F* MC_histo(TString prestring, TString dir_in, TString var, TFile* file_in, double xs, long Nevents, int rebin) {

  cout << file_in->GetName() << endl;

  double lumi = 41.529*pow(10,3); //luminosity in pb^-1

  double e_Nevents = pow(Nevents,0.5);
  double e_xs = 0.01*xs;

  //Weight
  double w = xs*lumi/Nevents;
  cout << "Events in data/events in MC " << xs*lumi/Nevents << endl;
  
  var = dir_in+prestring+var;
  cout << var << endl;
  TH1F* h = (TH1F*) file_in -> Get(var);
  TH1F* h1 = (TH1F*) h->Clone("tempname");
  h1->SetDirectory(0);
  //file_in->Close();

  h1 -> Scale(w);
  h1 -> Rebin(rebin);

  return h1;

}



int main(/*int argc, char** argv*/) {


  //string option = *(argv + 1);

  int rebin = 5;

  TString folder_in; 
  TString name_out;
  folder_in = "TauID2017/PzetaTest";
  name_out = "pzetatest";

  TFile* file_in_DY = new TFile(folder_in+"/Arranged_DY/DY.root", "R");
  TFile* file_in_WJets = new TFile(folder_in+"/Arranged_WJets/WJets.root", "R");
  TFile* file_in_QCD = new TFile(folder_in+"/Arranged_QCD/QCD.root", "R");


  TFile* file_in_TT_had = new TFile(folder_in+"/Arranged_TT/TT_had.root", "R");
  TFile* file_in_TT_semilep = new TFile(folder_in+"/Arranged_TT/TT_semilep.root", "R");
  TFile* file_in_TT_2l2nu = new TFile(folder_in+"/Arranged_TT/TT_2l2nu.root", "R");
  TFile* file_in_WW = new TFile(folder_in+"/Arranged_WW/WW.root", "R");
  TFile* file_in_WZ = new TFile(folder_in+"/Arranged_WZ/WZ.root", "R");
  TFile* file_in_ZZ = new TFile(folder_in+"/Arranged_ZZ/ZZ.root", "R");
  TFile* file_in_data = new TFile(folder_in+"/Arranged_data/data.root", "R");
  

  vector<TString> vars;
  //vars.push_back("tau_MVA");   
  vars.push_back("ev_DRmutau");
  //vars.push_back("ev_Mt_raw"); 
  vars.push_back("ev_Pzeta1"); 
  vars.push_back("ev_Pzeta2");

  

  //cross-sections
  double xs_DY = 6225.42;//5675.4; 
  double xs_WJets = 61526.7;


  double xs_TT_had = 831.76*0.457;//831.76;                 
  double xs_TT_semilep = 831.76*0.438;//831.76;                 
  double xs_TT_2l2nu = 831.76*0.105;//831.76;                 
  double xs_WW = 64.3;
  double xs_WZ = 23.4;
  double xs_ZZ = 10.16;

  //Nevents
  long N_DY = 141442384;//142161151;
  long N_WJets = 23219762;//25950745;

                  
  long N_TT_2l2nu = 8705576;
  long N_TT_semilep = 41221873;
  long N_TT_had = 42357944;
  long N_WW = 7791498;
  long N_WZ = 3928630;
  long N_ZZ = 1949768;


  TString var_in;
  TString dir_in;

  TFile* file_out = new TFile("Figures/"+name_out+".root", "RECREATE");
  file_out->cd();


  //options = is it the DY Sig?, variable name, which file to get the histo from, process cross-section
  for (unsigned int i = 0; i<vars.size(); ++i) {
    var_in = vars[i];
    dir_in = "NoTauID/";
    cout << var_in << endl;
	      
    TH1F* h_DYSig = MC_histo("Sig_", dir_in, var_in, file_in_DY, xs_DY, N_DY, rebin);
    h_DYSig->SetName("DYS_"+var_in);
    h_DYSig->SetTitle("DYS_"+var_in);
		
    h_DYSig->Write();

    TH1F* h_DYBkg = MC_histo("Bkg_", dir_in, var_in, file_in_DY, xs_DY, N_DY, rebin);
    h_DYBkg->SetName("DYB_"+var_in);
    
    h_DYBkg->Write();

    TH1F* h_WJets = MC_histo("Bkg_", dir_in, var_in, file_in_WJets, xs_WJets, N_WJets, rebin);
    h_WJets -> SetName("WJets_"+var_in);
    h_WJets->Write();
      	      
	      
    TH1F* h_TT_had = MC_histo("Bkg_", dir_in, var_in, file_in_TT_had, xs_TT_had, N_TT_had, rebin);
    TH1F* h_TT_semilep = MC_histo("Bkg_", dir_in, var_in, file_in_TT_semilep, xs_TT_semilep, N_TT_semilep, rebin);
    TH1F* h_TT_2l2nu = MC_histo("Bkg_", dir_in, var_in, file_in_TT_2l2nu, xs_TT_2l2nu, N_TT_2l2nu, rebin);
    TH1F* h_TTBkg = (TH1F*) h_TT_had->Clone("TTB_"+var_in);
    h_TTBkg->Add(h_TT_semilep);
    h_TTBkg->Add(h_TT_2l2nu);

    h_TTBkg->Write();
	      

    TH1F* h_WW = MC_histo("Bkg_", dir_in, var_in, file_in_WW, xs_WW, N_WW, rebin);
    TH1F* h_WZ = MC_histo("Bkg_", dir_in, var_in, file_in_WZ, xs_WZ, N_WZ, rebin);
    TH1F* h_ZZ = MC_histo("Bkg_", dir_in, var_in, file_in_ZZ, xs_ZZ, N_ZZ, rebin);
    TH1F* h_VV = (TH1F*) h_WW->Clone("VV_"+var_in);
    h_VV->Add(h_WZ);
    h_VV->Add(h_ZZ);

    h_VV->Write();
	      

    TH1F* h_data = (TH1F*) file_in_data -> Get(dir_in+"/Bkg_"+var_in);
    h_data->Rebin(rebin);
    h_data -> SetName("data_"+var_in);

    h_data->Write();
	      
  }//i
  file_out->Close();


  return 0;
}
