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


TH1F* MC_histo(TString prestring, TString var, TFile* file_in, double xs, long Nevents, int rebin) {

  cout << file_in->GetName() << endl;

  double lumi = 41525.735; //luminosity in pb^-1

  double e_Nevents = pow(Nevents,0.5);
  double e_xs = 0.01*xs;

  //Weight
  double w = xs*lumi/Nevents;
  cout << "Events in data/events in MC " << xs*lumi/Nevents << endl;
  
  var = prestring+var;
  TH1F* h = (TH1F*) file_in -> Get(var);

  h -> Sumw2();
  h -> Scale(w);
  h -> Rebin(rebin);

  return h;

}



int main(int argc, char** argv) {


  string option = *(argv + 1);
  string option_name= option;
  //option = "QCD", "QCD_control", "final", "Wjets_pre", "Wjets_post", "quick_test"

  int rebin = 1;
  if (option_name == "Fakes") rebin=1;
  float Wjets_scale_factor_notauID = 1.22504;//1.16223;//1.24298;
  bool final = false;
  bool antimu = false;
  bool wjets_pre = false;
  if (option_name == "final" || option_name == "final_control") final = true;
  if (option_name == "AntiMu" || option_name == "AntiMu_control") antimu = true;
  if (option_name == "WJets_pre" || option_name == "WJets_post") wjets_pre = true;
  if (option_name == "final_control" || option_name == "quick_test" || option_name == "WJets_post" || option_name == "Wjets_post" || option_name == "QCD_control" || option_name == "AntiMu_control") rebin = 5;
  cout << rebin << endl;

  TString folder_in; 
  TString name_out;
  if (final || option_name == "quick_test") {
    folder_in = "TauID2017/Final";
    name_out = "allhistos";
    if (option_name == "final_control") name_out = "allhistos_control";
  }
  else if (option_name == "QCD" || option_name == "QCD_control") {
    folder_in = "TauID2017/QCD";
    if (option_name == "QCD_control") {
      name_out = "QCD_for_plotting";
    }
    else {
      name_out = "QCD_pre_estimation";
    }
  }
  else if (option_name == "Wjets_pre" || option_name == "WJets_pre") {
    folder_in = "TauID2017/WJets";
    name_out = "Wjets_pre_estimation";
  }
  else if (option_name == "Wjets_post" || option_name == "WJets_post") {
    folder_in = "TauID2017/WJets";
    name_out = "Wjets_post_estimation";
  }
  else if (antimu) {
    folder_in = "TauID2017/AntiMu";
    name_out = "antimu";
    if(option_name == "AntiMu_control") name_out = "antimu_control";
  }
  else if (option_name == "Fakes") {
    folder_in = "TauID2017/TauFakes";
    name_out = "TauFakes";
  }

  TFile* file_out = new TFile("Figures/"+name_out+".root", "RECREATE");
  TFile* file_in_DY = new TFile(folder_in+"/Arranged_DY/DY.root", "R");
  TFile* file_in_WJets = new TFile(folder_in+"/Arranged_WJets/WJets.root", "R");
  TFile* file_in_QCD = new TFile(folder_in+"/Arranged_QCD/QCD.root", "R");


  vector<TFile*> QCD_files;
  if (antimu) {
    TFile* file_in_QCD_15to30 = new TFile(folder_in+"/Arranged_QCD/QCD_15to30.root", "R");              QCD_files.push_back(file_in_QCD_15to30);
    TFile* file_in_QCD_30to50 = new TFile(folder_in+"/Arranged_QCD/QCD_30to50.root", "R");              QCD_files.push_back(file_in_QCD_30to50);
    TFile* file_in_QCD_50to80 = new TFile(folder_in+"/Arranged_QCD/QCD_50to80.root", "R");              QCD_files.push_back(file_in_QCD_50to80);
    TFile* file_in_QCD_80to120 = new TFile(folder_in+"/Arranged_QCD/QCD_80to120.root", "R");            QCD_files.push_back(file_in_QCD_80to120);
    TFile* file_in_QCD_120to170 = new TFile(folder_in+"/Arranged_QCD/QCD_120to170.root", "R");          QCD_files.push_back(file_in_QCD_120to170);
    TFile* file_in_QCD_170to300 = new TFile(folder_in+"/Arranged_QCD/QCD_170to300.root", "R");          QCD_files.push_back(file_in_QCD_170to300);
    TFile* file_in_QCD_300to470 = new TFile(folder_in+"/Arranged_QCD/QCD_300to470.root", "R");          QCD_files.push_back(file_in_QCD_300to470);
    TFile* file_in_QCD_470to600 = new TFile(folder_in+"/Arranged_QCD/QCD_470to600.root", "R");          QCD_files.push_back(file_in_QCD_470to600);
    TFile* file_in_QCD_600to800 = new TFile(folder_in+"/Arranged_QCD/QCD_600to800.root", "R");          QCD_files.push_back(file_in_QCD_600to800);
    TFile* file_in_QCD_800to1000 = new TFile(folder_in+"/Arranged_QCD/QCD_800to1000.root", "R");        QCD_files.push_back(file_in_QCD_800to1000);
    TFile* file_in_QCD_1000to1400 = new TFile(folder_in+"/Arranged_QCD/QCD_1000to1400.root", "R");      QCD_files.push_back(file_in_QCD_1000to1400);
    TFile* file_in_QCD_1400to1800 = new TFile(folder_in+"/Arranged_QCD/QCD_1400to1800.root", "R");      QCD_files.push_back(file_in_QCD_1400to1800);
    TFile* file_in_QCD_1800to2400 = new TFile(folder_in+"/Arranged_QCD/QCD_1800to2400.root", "R");      QCD_files.push_back(file_in_QCD_1800to2400);
    TFile* file_in_QCD_2400to3200 = new TFile(folder_in+"/Arranged_QCD/QCD_2400to3200.root", "R");      QCD_files.push_back(file_in_QCD_2400to3200);
  }
  else {
    TFile* file_in_QCD_20to30 = new TFile(folder_in+"/Arranged_QCD/QCD_20to30.root", "R");              QCD_files.push_back(file_in_QCD_20to30);
    TFile* file_in_QCD_30to50 = new TFile(folder_in+"/Arranged_QCD/QCD_30to50.root", "R");              QCD_files.push_back(file_in_QCD_30to50);
    TFile* file_in_QCD_50to80 = new TFile(folder_in+"/Arranged_QCD/QCD_50to80.root", "R");              QCD_files.push_back(file_in_QCD_50to80);
    TFile* file_in_QCD_80to120 = new TFile(folder_in+"/Arranged_QCD/QCD_80to120.root", "R");            QCD_files.push_back(file_in_QCD_80to120);
    TFile* file_in_QCD_120to170 = new TFile(folder_in+"/Arranged_QCD/QCD_120to170.root", "R");          QCD_files.push_back(file_in_QCD_120to170);
    TFile* file_in_QCD_170to300 = new TFile(folder_in+"/Arranged_QCD/QCD_170to300.root", "R");          QCD_files.push_back(file_in_QCD_170to300);
    TFile* file_in_QCD_300to470 = new TFile(folder_in+"/Arranged_QCD/QCD_300to470.root", "R");          QCD_files.push_back(file_in_QCD_300to470);
    TFile* file_in_QCD_470to600 = new TFile(folder_in+"/Arranged_QCD/QCD_470to600.root", "R");          QCD_files.push_back(file_in_QCD_470to600);
    TFile* file_in_QCD_600to800 = new TFile(folder_in+"/Arranged_QCD/QCD_600to800.root", "R");          QCD_files.push_back(file_in_QCD_600to800);
    TFile* file_in_QCD_800to1000 = new TFile(folder_in+"/Arranged_QCD/QCD_800to1000.root", "R");        QCD_files.push_back(file_in_QCD_800to1000);
    TFile* file_in_QCD_1000toInf = new TFile(folder_in+"/Arranged_QCD/QCD_1000toInf.root", "R");        QCD_files.push_back(file_in_QCD_1000toInf);
  }    

  TFile* file_in_TT_had = new TFile(folder_in+"/Arranged_TT/TT_had.root", "R");
  TFile* file_in_TT_semilep = new TFile(folder_in+"/Arranged_TT/TT_semilep.root", "R");
  TFile* file_in_TT_2l2nu = new TFile(folder_in+"/Arranged_TT/TT_2l2nu.root", "R");
  TFile* file_in_WW = new TFile(folder_in+"/Arranged_WW/WW.root", "R");
  TFile* file_in_WZ = new TFile(folder_in+"/Arranged_WZ/WZ.root", "R");
  TFile* file_in_ZZ = new TFile(folder_in+"/Arranged_ZZ/ZZ.root", "R");
  TFile* file_in_data = new TFile(folder_in+"/Arranged_data/data.root", "R");
  TFile* file_in_fakes = new TFile("TauID2017/Faketaus.root", "R");
  

  vector<TString> vars;
  //tau-iD independent
  //vars.push_back("tau_MVA");//0
  //vars.push_back("ev_DRmutau");//4
  //vars.push_back("ev_Mt_raw");//5
  //vars.push_back("ev_Pzeta");

  //tau-ID dependent
  int tau_dependent_cutoff = vars.size();
  vars.push_back("ev_Mt");
  vars.push_back("ev_Mvis");
  //vars.push_back("ev_Mtot");
  vars.push_back("tau_pt");
  vars.push_back("tau_eta");
  vars.push_back("tau_phi");
  vars.push_back("mu_pt");
  vars.push_back("mu_eta");
  vars.push_back("mu_phi");
  //vars.push_back("ev_METmumass");
  //vars.push_back("ev_MET");
  //vars.push_back("ev_METphi");
  vars.push_back("ev_Nvertex");
  if (option_name != "Wjets_pre" && option_name != "WJets_pre" && option_name != "Wjets_post" && option_name != "WJets_post" && option_name != "AntiMu") {
    vars.push_back("ev_Mvis_TESdown");
    vars.push_back("ev_Mvis_TESup");
  }
  if (antimu) {
    vars.push_back("ev_Mvis_SS");
    vars.push_back("tau_pt_SS");
    vars.push_back("mu_pt_SS");
  }

  vector<TString> tauIDs;
  //tauIDs.push_back("noTauID");
  //tauIDs.push_back("cutbased_loose");
  //tauIDs.push_back("cutbased_medium");
  //tauIDs.push_back("cutbased_tight");
  //tauIDs.push_back("MVA_veryloose");
  //tauIDs.push_back("MVA_loose");
  //tauIDs.push_back("MVA_medium");
  tauIDs.push_back("MVA_tight");
  //tauIDs.push_back("MVA_verytight");
  //tauIDs.push_back("MVA_veryverytight");
  //tauIDs.push_back("MVAnew_veryloose");
  //tauIDs.push_back("MVAnew_loose");
  //tauIDs.push_back("MVAnew_medium");
  //tauIDs.push_back("MVAnew_tight");
  //tauIDs.push_back("MVAnew_verytight");
  //tauIDs.push_back("MVAnew_veryverytight");

  vector<TString> categ;
  categ.push_back("pass");
  categ.push_back("fail");

  //cross-sections
  double xs_DY = 5675.4; 
  double xs_WJets = 61526.7;

  double xs_QCD_muenriched = 720648000*0.00042;//0.0003739 //Mu-enriched sample, this includes filter efficiency
  vector<double> xs_QCD;
  if (antimu) {
    double xs_QCD_15to30 = 1.822*pow(10,9);    xs_QCD.push_back(xs_QCD_15to30);
    double xs_QCD_30to50 = 1.387*pow(10,8);    xs_QCD.push_back(xs_QCD_30to50);
    double xs_QCD_50to80 = 1.913*pow(10,7);    xs_QCD.push_back(xs_QCD_50to80);
    double xs_QCD_80to120 = 2.736*pow(10,6);   xs_QCD.push_back(xs_QCD_80to120);
    double xs_QCD_120to170 = 4.663*pow(10,5);  xs_QCD.push_back(xs_QCD_120to170);
    double xs_QCD_170to300 = 1.172*pow(10,5);  xs_QCD.push_back(xs_QCD_170to300);
    double xs_QCD_300to470 = 7.76*pow(10,3);   xs_QCD.push_back(xs_QCD_300to470);
    double xs_QCD_470to600 = 640.5;            xs_QCD.push_back(xs_QCD_470to600);
    double xs_QCD_600to800 = 185.9;            xs_QCD.push_back(xs_QCD_600to800);
    double xs_QCD_800to1000 = 32;              xs_QCD.push_back(xs_QCD_800to1000);
    double xs_QCD_1000to1400 = 9.37;           xs_QCD.push_back(xs_QCD_1000to1400);
    double xs_QCD_1400to1800 = 0.838;          xs_QCD.push_back(xs_QCD_1400to1800);
    double xs_QCD_1800to2400 = 0.112;          xs_QCD.push_back(xs_QCD_1800to2400);
    double xs_QCD_2400to3200 = 6.74*pow(10,-3); xs_QCD.push_back(xs_QCD_2400to3200);
  }
  else {
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
  }  

  double xs_TT_had = 831.76*0.457;//831.76;                 
  double xs_TT_semilep = 831.76*0.438;//831.76;                 
  double xs_TT_2l2nu = 831.76*0.105;//831.76;                 
  double xs_WW = 64.3;
  double xs_WZ = 23.4;
  //double xs_WZTo3LNu = 1;
  //double xs_ZZTo2L2Nu = 1;
  double xs_ZZ = 10.16;

  //Nevents
  long N_DY = 96844363; //144230225;
  long N_WJets = 23219762;//25950745;

  long N_QCD_muenriched = 7373309;
  vector<long> N_QCD;
  if (antimu) {
    long N_QCD_15to30 = 37585689;              N_QCD.push_back(N_QCD_15to30);
    long N_QCD_30to50 = 9979945;               N_QCD.push_back(N_QCD_30to50);
    long N_QCD_50to80 = 9954259;               N_QCD.push_back(N_QCD_50to80);
    long N_QCD_80to120 = 7608728+6986638;      N_QCD.push_back(N_QCD_80to120);
    long N_QCD_120to170 = 5504047+6324260;     N_QCD.push_back(N_QCD_120to170);
    long N_QCD_170to300 = 6855630+6799704;     N_QCD.push_back(N_QCD_170to300);
    long N_QCD_300to470 = 4150323+14771394;    N_QCD.push_back(N_QCD_300to470);
    long N_QCD_470to600 = 3866704;             N_QCD.push_back(N_QCD_470to600);
    long N_QCD_600to800 = 3810427+9496102;     N_QCD.push_back(N_QCD_600to800);
    long N_QCD_800to1000 = 13715052;           N_QCD.push_back(N_QCD_800to1000);
    long N_QCD_1000to1400 = 2829635+6539554;   N_QCD.push_back(N_QCD_1000to1400);
    long N_QCD_1400to1800 = 312291+2310326;    N_QCD.push_back(N_QCD_1400to1800);
    long N_QCD_1800to2400 = 397083+1549881;    N_QCD.push_back(N_QCD_1800to2400);
    long N_QCD_2400to3200 = 398491+550068;     N_QCD.push_back(N_QCD_2400to3200);
  }
  else {
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
  }    
                  
  long N_TT_2l2nu = 8615776;
  long N_TT_semilep = 109715126;
  long N_TT_had = 127935682;
  long N_WW = 3928630;
  long N_WZ = 3928630;
  long N_ZZ = 1949768;


  TString var_in;

  file_out->cd();
  //options = is it the DY Sig?, variable name, which file to get the histo from, process cross-section
  for (unsigned int i = 0; i<vars.size(); ++i) {

    if (i < tau_dependent_cutoff) {
      var_in = vars[i];

      if (final || option_name == "quick_test" || option_name == "Fakes") {
	TH1F* h_DYSig = MC_histo("Sig_", var_in, file_in_DY, xs_DY, N_DY, rebin);
	h_DYSig->SetName("DYS_"+var_in);
	h_DYSig->Write();
      }
      

      TH1F* h_DYBkg = MC_histo("Bkg_", var_in, file_in_DY, xs_DY, N_DY, rebin);
      h_DYBkg->SetName("DYB_"+var_in);
      h_DYBkg->Write();


      TH1F* h_WJets = (TH1F*) file_in_WJets -> Get("Bkg_"+var_in); //already normed
      h_WJets->Rebin(rebin);
      if (option_name != "Wjets_pre" && option_name != "WJets_pre") h_WJets->Scale(Wjets_scale_factor_notauID);
      h_WJets -> SetName("WJets_"+var_in);
      h_WJets->Write();
      
      if (option_name == "Wjets_pre" || option_name == "Wjets_post" || option_name == "WJets_pre" || option_name == "WJets_post" || option_name == "quick_test" || option_name == "QCD_control" || antimu) {
	vector<TH1F*> h_QCD_vector;
	for (unsigned int iBin = 0; iBin<QCD_files.size(); ++iBin) {
	  h_QCD_vector.push_back( MC_histo("Bkg_", var_in, QCD_files[iBin], xs_QCD[iBin], N_QCD[iBin], rebin) );
	}
	TH1F* h_QCD = (TH1F*) h_QCD_vector[0]->Clone("QCD_"+var_in);
	for (unsigned int iBin = 1; iBin<QCD_files.size(); ++iBin) {
	  h_QCD->Add(h_QCD_vector[iBin]);
	}
	h_QCD->Write();
      }

      
      TH1F* h_TT_had = MC_histo("Bkg_", var_in, file_in_TT_had, xs_TT_had, N_TT_had, rebin);
      TH1F* h_TT_semilep = MC_histo("Bkg_", var_in, file_in_TT_semilep, xs_TT_semilep, N_TT_semilep, rebin);
      TH1F* h_TT_2l2nu = MC_histo("Bkg_", var_in, file_in_TT_2l2nu, xs_TT_2l2nu, N_TT_2l2nu, rebin);
      TH1F* h_TTBkg = (TH1F*) h_TT_had->Clone("TTB_"+var_in);
      h_TTBkg->Add(h_TT_semilep);
      h_TTBkg->Add(h_TT_2l2nu);
      h_TTBkg->Write();
      
      TH1F* h_WW = MC_histo("Bkg_", var_in, file_in_WW, xs_WW, N_WW, rebin);
      TH1F* h_WZ = MC_histo("Bkg_", var_in, file_in_WZ, xs_WZ, N_WZ, rebin);
      TH1F* h_ZZ = MC_histo("Bkg_", var_in, file_in_ZZ, xs_ZZ, N_ZZ, rebin);
      TH1F* h_VV = (TH1F*) h_WW->Clone("VV_"+var_in);
      h_VV->Add(h_WZ);
      h_VV->Add(h_ZZ);
      //h_VV -> SetName("VV_"+var_in);
      h_VV->Write();

      TH1F* h_data = (TH1F*) file_in_data -> Get("Bkg_"+var_in);//Data is, by definition, normalized
      h_data -> SetName("data_"+var_in);
      h_data->Rebin(rebin);
      h_data->Write();
    }//end tau-ID independent histos
    else {
      for (unsigned int j = 0; j<tauIDs.size(); ++j) {
	for (unsigned int k = 0; k<categ.size(); ++k) {
	  //if (j == 0 && k == 0) {
	  //  var_in = vars[i] + "_" + tauIDs[j];
	  //}
	  //else if (j == 0 && k == 1) {
	  //  continue;
	  //}
	  //else {
	  var_in = vars[i] + "_" + tauIDs[j] + "_" + categ[k];
	  //}
	  cout << var_in << endl;

	  if (final || option_name == "quick_test" || option_name == "Fakes") {
	    TH1F* h_DYSig = MC_histo("Sig_", var_in, file_in_DY, xs_DY, N_DY, rebin);
	    h_DYSig->SetName("DYS_"+var_in);
	    h_DYSig->Write();

	    //TH1F* h_TTSig = MC_histo("Sig_", var_in, file_in_TT, xs_TT, N_TT, rebin);
	    //h_TTSig->SetName("TTS_"+var_in);
	    //h_TTSig->Write();
	  }


	  TH1F* h_DYBkg = MC_histo("Bkg_", var_in, file_in_DY, xs_DY, N_DY, rebin);
	  /*TH1F* h_DY1Bkg = MC_histo("Bkg_", var_in, file_in_DY1, xs_DY1, N_DY1, rebin);
	  TH1F* h_DY2Bkg = MC_histo("Bkg_", var_in, file_in_DY2, xs_DY2, N_DY2, rebin);
	  TH1F* h_DY3Bkg = MC_histo("Bkg_", var_in, file_in_DY3, xs_DY3, N_DY3, rebin);
	  TH1F* h_DY4Bkg = MC_histo("Bkg_", var_in, file_in_DY4, xs_DY4, N_DY4, rebin);*/
	  
	  /*TH1F* h_DYBkg = (TH1F*) h_DY0Bkg->Clone("DYB_"+var_in);
	  h_DYBkg->Add(h_DY1Bkg);
	  h_DYBkg->Add(h_DY2Bkg);
	  h_DYBkg->Add(h_DY3Bkg);
	  h_DYBkg->Add(h_DY4Bkg);*/
	  h_DYBkg->SetName("DYB_"+var_in);
	  h_DYBkg->Write();


	  TH1F* h_WJets = (TH1F*) file_in_WJets -> Get("Bkg_"+var_in); //already normed
	  h_WJets->Rebin(rebin);
	  //if (option_name != "Wjets_pre" && option_name != "WJets_pre") h_WJets->Scale(Wjets_scale_factor[k][j]);
	  h_WJets -> SetName("WJets_"+var_in);
	  h_WJets->Write();
      
	  if (option_name == "Wjets_pre" || option_name == "Wjets_post" || option_name == "WJets_pre" || option_name == "WJets_post" || option_name == "quick_test" || option_name == "QCD_control" || antimu) {
	    vector<TH1F*> h_QCD_vector;
	    for (unsigned int iBin = 0; iBin<QCD_files.size(); ++iBin) {
	      h_QCD_vector.push_back( MC_histo("Bkg_", var_in, QCD_files[iBin], xs_QCD[iBin], N_QCD[iBin], rebin) );
	    }
	    TH1F* h_QCD = (TH1F*) h_QCD_vector[0]->Clone("QCD_"+var_in);
	    for (unsigned int iBin = 1; iBin<QCD_files.size(); ++iBin) {
	      h_QCD->Add(h_QCD_vector[iBin]);
	    }
	    h_QCD->Write();
	  }

	  TH1F* h_TT_had = MC_histo("Bkg_", var_in, file_in_TT_had, xs_TT_had, N_TT_had, rebin);
	  TH1F* h_TT_semilep = MC_histo("Bkg_", var_in, file_in_TT_semilep, xs_TT_semilep, N_TT_semilep, rebin);
	  TH1F* h_TT_2l2nu = MC_histo("Bkg_", var_in, file_in_TT_2l2nu, xs_TT_2l2nu, N_TT_2l2nu, rebin);
	  TH1F* h_TTBkg = (TH1F*) h_TT_had->Clone("TTB_"+var_in);
	  h_TTBkg->Add(h_TT_semilep);
	  h_TTBkg->Add(h_TT_2l2nu);
	  h_TTBkg->Write();

	  TH1F* h_WW = MC_histo("Bkg_", var_in, file_in_WW, xs_WW, N_WW, rebin);
	  TH1F* h_WZ = MC_histo("Bkg_", var_in, file_in_WZ, xs_WZ, N_WZ, rebin);
	  TH1F* h_ZZ = MC_histo("Bkg_", var_in, file_in_ZZ, xs_ZZ, N_ZZ, rebin);
	  TH1F* h_VV = (TH1F*) h_WW->Clone("VV_"+var_in);
	  h_VV->Add(h_WZ);
	  h_VV->Add(h_ZZ);
	  //h_VV -> SetName("VV_"+var_in);
	  h_VV->Write();

	  TH1F* h_data = (TH1F*) file_in_data -> Get("Bkg_"+var_in);
	  h_data->Rebin(rebin);
	  h_data -> SetName("data_"+var_in);
	  h_data->Write();

	  if (final && categ[k] == "pass") {
	    cout << file_in_fakes->GetName() << " " << "faketau_"+var_in <<  endl;
	    TH1F* h_fakes = (TH1F*) file_in_fakes -> Get("faketau_"+var_in);
	    cout << "b" << endl;
	    h_fakes->Rebin(rebin);
	    cout << "c" << endl;
	    h_fakes->Write();
	  }
	}
      }
    }
  }
  file_out->Close();


  return 0;
}
