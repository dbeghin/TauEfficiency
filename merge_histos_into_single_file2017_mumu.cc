#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"


using namespace std;



TH1F* MC_histo(TString prestring, TString var, TFile* file_in, TFile* file_in_d, double xs, int rebin) {

  cout << file_in->GetName() << endl;

  TH1F* h_events_data = (TH1F*) file_in_d->Get("weighted_events");
  double full_data = 7.3968038e+08;
  double succ_data_ratio = h_events_data->Integral()/full_data;
  cout << "succesfull data ratio " << succ_data_ratio << endl;
  double lumi = 41.529 * pow(10,3) * succ_data_ratio; //luminosity in pb^-1

  TH1F* h_events = (TH1F*) file_in->Get("weighted_events");
  double Nevents = h_events->Integral();

  double e_Nevents = pow(Nevents,0.5);
  double e_xs = 0.01*xs;

  //Weight
  double w = 0;
  if (Nevents != 0) w = xs*lumi/Nevents;
  cout << "Events in data/events in MC " << w << endl;
  
  TH1F* h;
  h = (TH1F*) file_in -> Get(prestring+var);

  h -> Scale(w);
  h -> Rebin(rebin);
  
  return h;

}





int main(int argc, char** argv) {


  //int rebin = 2;
  //int rebin = 60;
  int rebin = 1;

  TString folder_in = "TauID2017/MuMu"; 
  TString name_out = "histos_mumu";


  TFile* file_out = new TFile("Figures/"+name_out+".root", "RECREATE");
  TFile* file_in_DY = new TFile(folder_in+"/Arranged_DY/DY.root", "R");
  TFile* file_in_WJets = new TFile(folder_in+"/Arranged_WJets/WJets.root", "R");
  TFile* file_in_TT = new TFile(folder_in+"/Arranged_TT/TT_2l2nu.root", "R");
  TFile* file_in_WW = new TFile(folder_in+"/Arranged_WW/WW.root", "R");
  TFile* file_in_WZ = new TFile(folder_in+"/Arranged_WZ/WZ.root", "R");
  TFile* file_in_ZZ = new TFile(folder_in+"/Arranged_ZZ/ZZ.root", "R");
  TFile* file_in_data = new TFile(folder_in+"/Arranged_data/data.root", "R");
  

  vector<TString> vars;                    vector<int> rrebin;
  vars.push_back("mu_pt");                 rrebin.push_back(2);
  vars.push_back("mu_eta");	 	   rrebin.push_back(2);
  vars.push_back("mu_phi");	 	   rrebin.push_back(2);
  vars.push_back("mu1_pt");	 	   rrebin.push_back(2);
  vars.push_back("mu1_eta");	 	   rrebin.push_back(2);
  vars.push_back("mu1_phi");	 	   rrebin.push_back(2);
  vars.push_back("mu2_pt");	 	   rrebin.push_back(2);
  vars.push_back("mu2_eta");	 	   rrebin.push_back(2);
  vars.push_back("mu2_phi");	 	   rrebin.push_back(2);
  vars.push_back("ev_DRmumu");	 	   rrebin.push_back(2);	 
  vars.push_back("ev_Mt_raw");	 	   rrebin.push_back(2);	 
  vars.push_back("ev_Mt");	 	   rrebin.push_back(2);
  vars.push_back("ev_Mvis");	 	   rrebin.push_back(2);//FIXME
  //vars.push_back("ev_Mvis_SS");	   rrebin.push_back(2);	 
  //vars.push_back("ev_METmumass");	   rrebin.push_back(2);
  vars.push_back("ev_MET");	 	   rrebin.push_back(2);
  vars.push_back("ev_METphi");	 	   rrebin.push_back(2);	 
  vars.push_back("ev_Nvertex");  	   rrebin.push_back(1);



  //cross-sections
  double xs_DY = 6225.42;
  //double xs_DY = 5675.4;
  double xs_WJets = 61526.7;
  double xs_TT = 87.31;//831.76;
  double xs_WW = 64.3;
  double xs_WZ = 23.4;
  double xs_ZZ = 10.16;



  TString var_in;

  file_out->cd();
  //options = is it the DY Sig?, variable name, which file to get the histo from, process cross-section
  for (unsigned int i = 0; i<vars.size(); ++i) {

    var_in = vars[i];
    int rebin = rrebin[i];
    
    TH1F* h_DYSig = MC_histo("Sig_", var_in, file_in_DY, file_in_data, xs_DY, rebin);
    h_DYSig -> SetName("DYB_"+var_in);
    h_DYSig->Write();
    TH1F* h_DYBkg = MC_histo("Bkg_", var_in, file_in_DY, file_in_data, xs_DY, rebin);
    h_DYBkg -> SetName("DYS_"+var_in);
    h_DYBkg->Write();
    
    //TH1F* h_WJets = MC_histo(true, var_in, file_in_WJets, xs_WJets, N_WJets, rebin);
    //h_WJets -> SetName("WJets_"+var_in);
    //h_WJets->Write();
      
    TH1F* h_TT = MC_histo("Sig_", var_in, file_in_TT, file_in_data, xs_TT, rebin);
    h_TT -> SetName("TTB_"+var_in);
    h_TT->Write();
    
    TH1F* h_WW = MC_histo("Sig_", var_in, file_in_WW, file_in_data, xs_WW, rebin);
    TH1F* h_WZ = MC_histo("Sig_", var_in, file_in_WZ, file_in_data, xs_WZ, rebin);
    TH1F* h_ZZ = MC_histo("Sig_", var_in, file_in_ZZ, file_in_data, xs_ZZ, rebin);
    cout << "before cloning" << endl;
    TH1F* h_VV = (TH1F*) h_WW->Clone("VV_"+var_in);
    cout << "after cloning" << endl;
    h_VV->Add(h_WZ);
    h_VV->Add(h_ZZ);
    //h_VV -> SetName("VV_"+var_in);
    h_VV->Write();
    
    cout << "?" << endl;

    TH1F* h_data = (TH1F*) file_in_data -> Get("Sig_"+var_in);//Data is, by definition, normalized
    h_data -> SetName("data_"+var_in);
    h_data->Rebin(rebin);
    h_data->Write();
  }
  file_out->Close();


  return 0;
}
