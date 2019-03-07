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


TH1F* MC_histo(TString var, TFile* file_in, double xs, long Nevents, int rebin) {

  cout << file_in->GetName() << endl;

  double lumi = 41529; //luminosity in pb^-1

  double e_Nevents = pow(Nevents,0.5);
  double e_xs = 0.01*xs;

  //Weight
  double w = xs*lumi/Nevents;
  cout << "Events in data/events in MC " << xs*lumi/Nevents << endl;
  
  TH1F* h;

  cout << var << endl;
  h = (TH1F*) file_in -> Get(var);

  h -> Sumw2();
  h -> Scale(w);
  //h -> GetXaxis()->SetRange(60,120);
  h -> Rebin(rebin);

  return h;

}



int main(int argc, char** argv) {


  int rebin = 2;//60;//2;

  TString folder_in = "TauID2017/FakeRate"; 
  TString name_out = "histos_fakerate";


  TFile* file_out = new TFile("Figures/"+name_out+".root", "RECREATE");
  TFile* file_in_DY = new TFile(folder_in+"/Arranged_DY/DY.root", "R");
  TFile* file_in_WJets = new TFile(folder_in+"/Arranged_WJets/WJets.root", "R");
  TFile* file_in_TT = new TFile(folder_in+"/Arranged_TT/TT_2l2nu.root", "R");
  TFile* file_in_WW = new TFile(folder_in+"/Arranged_WW/WW.root", "R");
  TFile* file_in_WZ = new TFile(folder_in+"/Arranged_WZ/WZ.root", "R");
  TFile* file_in_ZZ = new TFile(folder_in+"/Arranged_ZZ/ZZ.root", "R");
  TFile* file_in_data = new TFile(folder_in+"/Arranged_data/data.root", "R");
  
  //vector<TFile*> QCD_files;
  //TFile* file_in_QCD_20to30 = new TFile(folder_in+"/Arranged_QCD/QCD_20to30.root", "R");              QCD_files.push_back(file_in_QCD_20to30);
  //TFile* file_in_QCD_30to50 = new TFile(folder_in+"/Arranged_QCD/QCD_30to50.root", "R");              QCD_files.push_back(file_in_QCD_30to50);
  //TFile* file_in_QCD_50to80 = new TFile(folder_in+"/Arranged_QCD/QCD_50to80.root", "R");              QCD_files.push_back(file_in_QCD_50to80);
  //TFile* file_in_QCD_80to120 = new TFile(folder_in+"/Arranged_QCD/QCD_80to120.root", "R");            QCD_files.push_back(file_in_QCD_80to120);
  //TFile* file_in_QCD_120to170 = new TFile(folder_in+"/Arranged_QCD/QCD_120to170.root", "R");          QCD_files.push_back(file_in_QCD_120to170);
  //TFile* file_in_QCD_170to300 = new TFile(folder_in+"/Arranged_QCD/QCD_170to300.root", "R");          QCD_files.push_back(file_in_QCD_170to300);
  //TFile* file_in_QCD_300to470 = new TFile(folder_in+"/Arranged_QCD/QCD_300to470.root", "R");          QCD_files.push_back(file_in_QCD_300to470);
  //TFile* file_in_QCD_470to600 = new TFile(folder_in+"/Arranged_QCD/QCD_470to600.root", "R");          QCD_files.push_back(file_in_QCD_470to600);
  //TFile* file_in_QCD_600to800 = new TFile(folder_in+"/Arranged_QCD/QCD_600to800.root", "R");          QCD_files.push_back(file_in_QCD_600to800);
  //TFile* file_in_QCD_800to1000 = new TFile(folder_in+"/Arranged_QCD/QCD_800to1000.root", "R");        QCD_files.push_back(file_in_QCD_800to1000);
  //TFile* file_in_QCD_1000toInf = new TFile(folder_in+"/Arranged_QCD/QCD_1000toInf.root", "R");        QCD_files.push_back(file_in_QCD_1000toInf);
  

  vector<TString> vars;
  vars.push_back("ev_Mvis");
  vars.push_back("mu_pt");
  vars.push_back("mu_eta");
  vars.push_back("mu_phi");
  vars.push_back("ev_Nvertex");
  int tau_cutoff = vars.size();
  vars.push_back("taupt");

  vector<TString> HPS_WP;
  HPS_WP.push_back("cutbased_loose");
  HPS_WP.push_back("cutbased_medium");
  HPS_WP.push_back("cutbased_tight");

  HPS_WP.push_back("MVA_2017v2vvloose");
  HPS_WP.push_back("MVA_2017v2vloose");
  HPS_WP.push_back("MVA_2017v2loose");
  HPS_WP.push_back("MVA_2017v2medium");
  HPS_WP.push_back("MVA_2017v2tight");
  HPS_WP.push_back("MVA_2017v2vtight");
  HPS_WP.push_back("MVA_2017v2vvtight");

  HPS_WP.push_back("MVA_DBdR03vloose");
  HPS_WP.push_back("MVA_DBdR03loose");
  HPS_WP.push_back("MVA_DBdR03medium");
  HPS_WP.push_back("MVA_DBdR03tight");
  HPS_WP.push_back("MVA_DBdR03vtight");
  HPS_WP.push_back("MVA_DBdR03vvtight");

  HPS_WP.push_back("MVA_PWdR03vloose");
  HPS_WP.push_back("MVA_PWdR03loose");
  HPS_WP.push_back("MVA_PWdR03medium");
  HPS_WP.push_back("MVA_PWdR03tight");
  HPS_WP.push_back("MVA_PWdR03vtight");
  HPS_WP.push_back("MVA_PWdR03vvtight");

  vector<TString> passfail;
  passfail.push_back("pass");
  passfail.push_back("fail");


  vector<TString> dms;
  dms.push_back("DM0");
  dms.push_back("DM1");
  dms.push_back("DM10");

  vector<TString> eta;
  eta.push_back("barrel");
  eta.push_back("endcap");

  //cross-sections
  double xs_DY = 6225.42;//5675.4;
  double xs_WJets = 61526.7;
  double xs_TT = 87.31;//831.76;
  double xs_WW = 64.3;
  double xs_WZ = 23.4;
  double xs_ZZ = 10.16;

  vector<double> xs_QCD;
  double xs_QCD_20to30 = 562700000*0.007087;    xs_QCD.push_back(xs_QCD_20to30);
  double xs_QCD_30to50 = 139900000*0.01219;     xs_QCD.push_back(xs_QCD_30to50);
  double xs_QCD_50to80 = 19400000*0.02037;      xs_QCD.push_back(xs_QCD_50to80);
  double xs_QCD_80to120 = 2762000*0.0387;       xs_QCD.push_back(xs_QCD_80to120);
  double xs_QCD_120to170 = 479500*0.04958;      xs_QCD.push_back(xs_QCD_120to170);
  double xs_QCD_170to300 = 118100*0.07022;      xs_QCD.push_back(xs_QCD_170to300);
  double xs_QCD_300to470 = 7820.25*0.10196;     xs_QCD.push_back(xs_QCD_300to470);
  double xs_QCD_470to600 = 648.8*0.08722;       xs_QCD.push_back(xs_QCD_470to600);
  double xs_QCD_600to800 = 187.109*0.13412;     xs_QCD.push_back(xs_QCD_600to800);
  double xs_QCD_800to1000 = 32.3486*0.14552;    xs_QCD.push_back(xs_QCD_800to1000);
  double xs_QCD_1000toInf = 10.4305*0.15544;    xs_QCD.push_back(xs_QCD_1000toInf);

  //Nevents
  long N_DY = 142161151;//144230225;
  long N_WJets = 23133163;
  long N_TT = 8705576;//8615776;
  long N_WW = 3928630;//7528052;
  long N_WZ = 7791498;//3928630;
  long N_ZZ = 1949768;//1949768;

  vector<long> N_QCD;
  long N_QCD_20to30 = 24545944+3460458;      N_QCD.push_back(N_QCD_20to30);
  long N_QCD_30to50 = 24604628+4063290;      N_QCD.push_back(N_QCD_30to50);
  long N_QCD_50to80 = 23955198;              N_QCD.push_back(N_QCD_50to80);
  long N_QCD_80to120 = 23098427;             N_QCD.push_back(N_QCD_80to120);
  long N_QCD_120to170 = 20821535;            N_QCD.push_back(N_QCD_120to170);
  long N_QCD_170to300 = 24561531+21877145;   N_QCD.push_back(N_QCD_170to300);
  long N_QCD_300to470 = 17620456;            N_QCD.push_back(N_QCD_300to470);
  long N_QCD_470to600 = 18774218;            N_QCD.push_back(N_QCD_470to600);
  long N_QCD_600to800 = 16392140;            N_QCD.push_back(N_QCD_600to800);
  long N_QCD_800to1000 = 15694987;           N_QCD.push_back(N_QCD_800to1000);
  long N_QCD_1000toInf = 11464778;           N_QCD.push_back(N_QCD_1000toInf);

  TString var_in;

  file_out->cd();
  //options = is it the DY Sig?, variable name, which file to get the histo from, process cross-section
  for (unsigned int i = 0; i<vars.size(); ++i) {
    for (unsigned int j = 0; j<dms.size(); ++j) {
      if ( (i<tau_cutoff) && (j>0) ) break;
      for (unsigned int k = 0; k<eta.size(); ++k) {
	if ( (i<tau_cutoff) && (k>0) ) break;
	for (unsigned int l = 0; l<HPS_WP.size(); ++l) {
	  if ( (i<tau_cutoff) && (l>0) ) break;
	  for (unsigned int m = 0; m<passfail.size(); ++m) {

	    if (i < tau_cutoff) {
	      if (m>0) break;
	      var_in = vars[i];
	    }
	    else {
	      var_in = vars[i]+"_"+HPS_WP[l]+"_"+passfail[m]+"_"+dms[j]+"_"+eta[k];
	    }
	    cout << endl << endl <<var_in << endl;
    
            TH1F* h_DY = MC_histo(var_in, file_in_DY, xs_DY, N_DY, rebin);
            h_DY -> SetName("DY_"+var_in);
            h_DY->Write();
            
            //vector<TH1F*> h_QCD_vector;
            //for (unsigned int iBin = 0; iBin<QCD_files.size(); ++iBin) {
            //  h_QCD_vector.push_back( MC_histo(var_in, QCD_files[iBin], xs_QCD[iBin], N_QCD[iBin], rebin) );
            //}
            //TH1F* h_QCD = (TH1F*) h_QCD_vector[0]->Clone("QCD_"+var_in);
            //for (unsigned int iBin = 1; iBin<QCD_files.size(); ++iBin) {
            //  h_QCD->Add(h_QCD_vector[iBin]);
            //}
            //h_QCD->Write();
              
            /*TH1F* h_WJets = MC_histo(var_in, file_in_WJets, xs_WJets, N_WJets, rebin);
            h_WJets -> SetName("WJets_"+var_in);
            h_WJets->Write();
              
            TH1F* h_TT = MC_histo(var_in, file_in_TT, xs_TT, N_TT, rebin);
            h_TT -> SetName("TT_"+var_in);
            h_TT->Write();
            
            TH1F* h_WW = MC_histo(var_in, file_in_WW, xs_WW, N_WW, rebin);
            //TH1F* h_WZ = MC_histo(var_in, file_in_WZ, xs_WZ, N_WZ, rebin);
            //TH1F* h_ZZ = MC_histo(var_in, file_in_ZZ, xs_ZZ, N_ZZ, rebin);
            TH1F* h_VV = (TH1F*) h_WW->Clone("VV_"+var_in);
            //h_VV->Add(h_WZ);
            //h_VV->Add(h_ZZ);
            //h_VV -> SetName("VV_"+var_in);
            h_VV->Write();*/
            
            
            TH1F* h_data = (TH1F*) file_in_data -> Get(+var_in);//Data is, by definition, normalized
            h_data -> SetName("data_"+var_in);
            h_data->Rebin(rebin);
            h_data->Write();
	  }
	}
      }
    }
  }
  file_out->Close();


  return 0;
}
