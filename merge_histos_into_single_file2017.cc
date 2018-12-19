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

  double lumi = 40.76*pow(10,3); //luminosity in pb^-1

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

  h1 -> Sumw2();
  h1 -> Scale(w);
  h1 -> Rebin(rebin);

  return h1;

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
    name_out = "WJets";
  }
  else if (option_name == "Wjets_post" || option_name == "WJets_post") {
    folder_in = "TauID2017/WJets";
    name_out = "Wjets_post_estimation";
  }
  else if (option_name == "WjetsFake" || option_name == "WJetsFake") {
    folder_in = "TauID2017/WJetsFake";
    name_out = "WJetsFake";
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
  TFile* file_in_fakes;
  if (final) {
    file_in_fakes = new TFile("TauID2017/Faketaus.root", "R");
  }
  else if (wjets_pre) {
    file_in_fakes = new TFile("TauID2017/wjetsfake.root", "R");
  }
  

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
  vars.push_back("ev_Nvertex");
  vars.push_back("ev_Mvis_TESdown");
  vars.push_back("ev_Mvis_TESup");
  vars.push_back("tau_pt_TESdown");          
  vars.push_back("tau_pt_TESup");            
  vars.push_back("Z_pt");                    
  vars.push_back("tau_DM");                  
  vars.push_back("ev_Mvis_MESdown");         
  vars.push_back("ev_Mvis_MESup");           
  vars.push_back("tau_pt_MESdown");          
  vars.push_back("tau_pt_MESup");            
  vars.push_back("ev_Mvis_MinBiasdown");     
  vars.push_back("ev_Mvis_MinBiasup");       
  //vars.push_back("ev_Nvertex_MinBiasdown");  
  //vars.push_back("ev_Nvertex_Minbiasup");    
  vars.push_back("ev_Mvis_FRS_DM0_down");    
  vars.push_back("ev_Mvis_FRS_DM0_up");      
  vars.push_back("tau_pt_FRS_DM0_down");     
  vars.push_back("tau_pt_FRS_DM0_up");       
  vars.push_back("ev_Mvis_FRS_DM1_down");    
  vars.push_back("ev_Mvis_FRS_DM1_up");      
  vars.push_back("tau_pt_FRS_DM1_down");     
  vars.push_back("tau_pt_FRS_DM1_up");       
  vars.push_back("ev_Mvis_FRS_DM10_down");   
  vars.push_back("ev_Mvis_FRS_DM10_up");     
  vars.push_back("tau_pt_FRS_DM10_down");    
  vars.push_back("tau_pt_FRS_DM10_up");
  vars.push_back("ev_Mvis_antiisomu_down");  
  vars.push_back("ev_Mvis_antiisomu_up");    
  vars.push_back("tau_pt_antiisomu_down");   
  vars.push_back("tau_pt_antiisomu_up");     
  vars.push_back("ev_Mvis_antiisotau_down"); 
  vars.push_back("ev_Mvis_antiisotau_up");   
  vars.push_back("tau_pt_antiisotau_down");  
  vars.push_back("tau_pt_antiisotau_up");
  if (antimu) {
    vars.push_back("ev_Mvis_SS");
    vars.push_back("tau_pt_SS");
    vars.push_back("mu_pt_SS");
  }

  vector<TString> HPS_WP;
  //HPS_WP.push_back("cutbased_loose");
  //HPS_WP.push_back("cutbased_medium");
  //HPS_WP.push_back("cutbased_tight");
  //
  //HPS_WP.push_back("MVA_2017v1vvloose");
  //HPS_WP.push_back("MVA_2017v1vloose");
  //HPS_WP.push_back("MVA_2017v1loose");
  //HPS_WP.push_back("MVA_2017v1medium");
  //HPS_WP.push_back("MVA_2017v1tight");
  //HPS_WP.push_back("MVA_2017v1vtight");
  //HPS_WP.push_back("MVA_2017v1vvtight");
  //
  //HPS_WP.push_back("MVA_2017v2vvloose");
  //HPS_WP.push_back("MVA_2017v2vloose");
  //HPS_WP.push_back("MVA_2017v2loose");
  //HPS_WP.push_back("MVA_2017v2medium");
  HPS_WP.push_back("MVA_2017v2tight");
  //HPS_WP.push_back("MVA_2017v2vtight");
  //HPS_WP.push_back("MVA_2017v2vvtight");
  //
  //HPS_WP.push_back("MVA_DBdR03vloose"); 
  //HPS_WP.push_back("MVA_DBdR03loose");
  //HPS_WP.push_back("MVA_DBdR03medium");
  //HPS_WP.push_back("MVA_DBdR03tight");
  //HPS_WP.push_back("MVA_DBdR03vtight");
  //HPS_WP.push_back("MVA_DBdR03vvtight");
  //
  //HPS_WP.push_back("MVA_PWdR03vloose"); 
  //HPS_WP.push_back("MVA_PWdR03loose");
  //HPS_WP.push_back("MVA_PWdR03medium");
  //HPS_WP.push_back("MVA_PWdR03tight");
  //HPS_WP.push_back("MVA_PWdR03vtight");
  //HPS_WP.push_back("MVA_PWdR03vvtight");

  vector<TString> categ;
  categ.push_back("pass");
  categ.push_back("fail");

  vector<TString> dms;
  //dms.push_back("allDMs");
  dms.push_back("DM0");
  dms.push_back("DM1");
  dms.push_back("DM10");

  vector<TString> eta;
  //eta.push_back("fulletarange");
  eta.push_back("barrel");
  eta.push_back("endcap");

  vector<TString> ptrange;
  ptrange.push_back("pt_20_40");
  ptrange.push_back("pt_40_150");
  


  //cross-sections
  double xs_DY = 6225.42;//5675.4; 
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
  double xs_ZZ = 10.16;

  //Nevents
  long N_DY = 142161151;
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
                  
  long N_TT_2l2nu = 8705576;
  long N_TT_semilep = 41221873;
  long N_TT_had = 42357944;
  long N_WW = 7791498;
  long N_WZ = 3928630;
  long N_ZZ = 1949768;


  vector<TString> process;
  process.push_back("DYS");
  process.push_back("DYB");
  process.push_back("TTB");
  process.push_back("VV"); 
  process.push_back("data");
  if (final || wjets_pre) process.push_back("faketau");


  TString var_in;
  TString dir_in;

  TFile* file_out = new TFile("Figures/"+name_out+".root", "RECREATE");
  file_out->cd();


  //options = is it the DY Sig?, variable name, which file to get the histo from, process cross-section
  for (unsigned int l = 0; l<HPS_WP.size(); ++l) {
    TDirectory* work_dir =  file_out->mkdir(HPS_WP[l]);
    work_dir->cd();

    for (unsigned int i = 0; i<vars.size(); ++i) {
      map<TString, TH1F*> h_cons;
      map<TString, TH1F*> h_cons_bypt[ptrange.size()];
      map<TString, TH1F*> h_cons_byDMeta[dms.size()][eta.size()];

      for (unsigned int j = 0; j<dms.size(); ++j) {
	if ((i<tau_dependent_cutoff) && (j>0)) break; 
	for (unsigned int k = 0; k<eta.size(); ++k) {
	  if ((i<tau_dependent_cutoff) && (k>0)) break; 
	  for (unsigned int m = 0; m<categ.size(); ++m) {
	    if ((i<tau_dependent_cutoff) && (m>0)) break;
	    for (unsigned int p = 0; p<ptrange.size(); ++p) {
	      if (i<tau_dependent_cutoff) {
	        var_in = vars[i];
		dir_in = "NoTauID/";
	        if (p>0) break;
	      }
	      else {
	        var_in = vars[i] + "_" + dms[j] + "_" + eta[k] + "_" + ptrange[p] + "_" + HPS_WP[l] + "_" + categ[m];
		dir_in = HPS_WP[l] + "/";
	      }
	      cout << var_in << endl;
	      
	      if (final || option_name == "quick_test" || option_name == "Fakes") {
	        TH1F* h_DYSig = MC_histo("Sig_", dir_in, var_in, file_in_DY, xs_DY, N_DY, rebin);
		h_DYSig->SetName("DYS_"+var_in);
		h_DYSig->SetTitle("DYS_"+var_in);
		
	        h_DYSig->Write();

		if (m==0) {
		  if (j==0 && k==0 && p==0) h_cons["DYS"] = (TH1F*) h_DYSig->Clone("DYS_"+vars[i]+"_allDMs_alleta_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		  else h_cons["DYS"]->Add(h_DYSig);

		  if (j==0 && k==0) h_cons_bypt[p]["DYS"] = (TH1F*) h_DYSig->Clone("DYS_"+vars[i]+"_allDMs_alleta_"+ptrange[p]+"_"+ HPS_WP[l]+"_"+categ[m]);
		  else h_cons_bypt[p]["DYS"]->Add(h_DYSig);

		  if (p==0) h_cons_byDMeta[j][k]["DYS"] = (TH1F*) h_DYSig->Clone("DYS_"+vars[i]+"_"+dms[j]+"_"+eta[k]+"_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		  else h_cons_byDMeta[j][k]["DYS"]->Add(h_DYSig);
		}
	      }
	      
	      
	      TH1F* h_DYBkg = MC_histo("Bkg_", dir_in, var_in, file_in_DY, xs_DY, N_DY, rebin);
	      h_DYBkg->SetName("DYB_"+var_in);
	      
	      h_DYBkg->Write();

	      if (m==0) {
		if (j==0 && k==0 && p==0) h_cons["DYB"] = (TH1F*) h_DYBkg->Clone("DYB_"+vars[i]+"_allDMs_alleta_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons["DYB"]->Add(h_DYBkg);
		
		if (j==0 && k==0) h_cons_bypt[p]["DYB"] = (TH1F*) h_DYBkg->Clone("DYB_"+vars[i]+"_allDMs_alleta_"+ptrange[p]+"_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_bypt[p]["DYB"]->Add(h_DYBkg);
		
		if (p==0) h_cons_byDMeta[j][k]["DYB"] = (TH1F*) h_DYBkg->Clone("DYB_"+vars[i]+"_"+dms[j]+"_"+eta[k]+"_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_byDMeta[j][k]["DYB"]->Add(h_DYBkg);
	      }
	      

	      //TH1F* h_WJets = (TH1F*) file_in_WJets -> Get("Bkg_"+var_in); //already normed
	      //h_WJets->Rebin(rebin);
	      //if (option_name != "Wjets_pre" && option_name != "WJets_pre") h_WJets->Scale(Wjets_scale_factor[k][j]);
	      //h_WJets -> SetName("WJets_"+var_in);
	      //h_WJets->Write();
      	      
	      /*if (option_name == "Wjets_pre" || option_name == "Wjets_post" || option_name == "WJets_pre" || option_name == "WJets_post" || option_name == "quick_test" || option_name == "QCD_control" || antimu) {
	        vector<TH1F*> h_QCD_vector;
	        for (unsigned int iBin = 0; iBin<QCD_files.size(); ++iBin) {
	          h_QCD_vector.push_back( MC_histo("Bkg_", dir_in, var_in, QCD_files[iBin], xs_QCD[iBin], N_QCD[iBin], rebin) );
	        }
	        TH1F* h_QCD = (TH1F*) h_QCD_vector[0]->Clone("QCD_"+var_in);
	        for (unsigned int iBin = 1; iBin<QCD_files.size(); ++iBin) {
	          h_QCD->Add(h_QCD_vector[iBin]);
	        }
	        h_QCD->Write();
	        }*/
	      
	      TH1F* h_TT_had = MC_histo("Bkg_", dir_in, var_in, file_in_TT_had, xs_TT_had, N_TT_had, rebin);
	      TH1F* h_TT_semilep = MC_histo("Bkg_", dir_in, var_in, file_in_TT_semilep, xs_TT_semilep, N_TT_semilep, rebin);
	      TH1F* h_TT_2l2nu = MC_histo("Bkg_", dir_in, var_in, file_in_TT_2l2nu, xs_TT_2l2nu, N_TT_2l2nu, rebin);
	      TH1F* h_TTBkg = (TH1F*) h_TT_had->Clone("TTB_"+var_in);
	      h_TTBkg->Add(h_TT_semilep);
	      h_TTBkg->Add(h_TT_2l2nu);

	      h_TTBkg->Write();
	      
	      if (m==0) {
		if (j==0 && k==0 && p==0) h_cons["TTB"] = (TH1F*) h_TTBkg->Clone("TTB_"+vars[i]+"_allDMs_alleta_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons["TTB"]->Add(h_TTBkg);
		
		if (j==0 && k==0) h_cons_bypt[p]["TTB"] = (TH1F*) h_TTBkg->Clone("TTB_"+vars[i]+"_allDMs_alleta_"+ptrange[p]+"_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_bypt[p]["TTB"]->Add(h_TTBkg);
		
		if (p==0) h_cons_byDMeta[j][k]["TTB"] = (TH1F*) h_TTBkg->Clone("TTB_"+vars[i]+"_"+dms[j]+"_"+eta[k]+"_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_byDMeta[j][k]["TTB"]->Add(h_TTBkg);
	      }


	      TH1F* h_WW = MC_histo("Bkg_", dir_in, var_in, file_in_WW, xs_WW, N_WW, rebin);
	      TH1F* h_WZ = MC_histo("Bkg_", dir_in, var_in, file_in_WZ, xs_WZ, N_WZ, rebin);
	      TH1F* h_ZZ = MC_histo("Bkg_", dir_in, var_in, file_in_ZZ, xs_ZZ, N_ZZ, rebin);
	      TH1F* h_VV = (TH1F*) h_WW->Clone("VV_"+var_in);
	      h_VV->Add(h_WZ);
	      h_VV->Add(h_ZZ);
	      //h_VV -> SetName("VV_"+var_in);

	      h_VV->Write();
	      
	      if (m==0) {
		if (j==0 && k==0 && p==0) h_cons["VV"] = (TH1F*) h_VV->Clone("VV_"+vars[i]+"_allDMs_alleta_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons["VV"]->Add(h_VV);
		
		if (j==0 && k==0) h_cons_bypt[p]["VV"] = (TH1F*) h_VV->Clone("VV_"+vars[i]+"_allDMs_alleta_"+ptrange[p]+"_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_bypt[p]["VV"]->Add(h_VV);
		
		if (p==0) h_cons_byDMeta[j][k]["VV"] = (TH1F*) h_VV->Clone("VV_"+vars[i]+"_"+dms[j]+"_"+eta[k]+"_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_byDMeta[j][k]["VV"]->Add(h_VV);
	      }


	      TH1F* h_data = (TH1F*) file_in_data -> Get(HPS_WP[l]+"/Bkg_"+var_in);
	      h_data->Rebin(rebin);
	      h_data -> SetName("data_"+var_in);

	      h_data->Write();
	      
	      if (m==0) {
		if (j==0 && k==0 && p==0) h_cons["data"] = (TH1F*) h_data->Clone("data_"+vars[i]+"_allDMs_alleta_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons["data"]->Add(h_data);
		
		if (j==0 && k==0) h_cons_bypt[p]["data"] = (TH1F*) h_data->Clone("data_"+vars[i]+"_allDMs_alleta_"+ptrange[p]+"_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_bypt[p]["data"]->Add(h_data);
		
		if (p==0) h_cons_byDMeta[j][k]["data"] = (TH1F*) h_data->Clone("data_"+vars[i]+"_"+dms[j]+"_"+eta[k]+"_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		else h_cons_byDMeta[j][k]["data"]->Add(h_data);
	      }


	      if ((final || wjets_pre) && categ[m] == "pass") {
	        cout << file_in_fakes->GetName() << " " << "faketau_"+var_in <<  endl;
	        TH1F* h_fakes = (TH1F*) file_in_fakes -> Get(HPS_WP[l]+"/faketau_"+var_in);
	        h_fakes->Rebin(rebin);

	        h_fakes->Write();

		if (m==0) {
		  if (j==0 && k==0 && p==0) h_cons["faketau"] = (TH1F*) h_fakes->Clone("faketau_"+vars[i]+"_allDMs_alleta_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		  else h_cons["faketau"]->Add(h_fakes);
		  
		  if (j==0 && k==0) h_cons_bypt[p]["faketau"] = (TH1F*) h_fakes->Clone("faketau_"+vars[i]+"_allDMs_alleta_"+ptrange[p]+"_"+ HPS_WP[l]+"_"+categ[m]);
		  else h_cons_bypt[p]["faketau"]->Add(h_fakes);
		  
		  if (p==0) h_cons_byDMeta[j][k]["faketau"] = (TH1F*) h_fakes->Clone("faketau_"+vars[i]+"_"+dms[j]+"_"+eta[k]+"_allpt_"+ HPS_WP[l]+"_"+categ[m]);
		  else h_cons_byDMeta[j][k]["faketau"]->Add(h_fakes);
		}
	      }
	    }//p
	  }//m
	  for (unsigned int dd=0; dd<process.size(); ++dd) h_cons_byDMeta[j][k][process[dd]]->Write();
	}//k
      }//j
      for (unsigned int dd=0; dd<process.size(); ++dd) {
	h_cons[process[dd]]->Write();
	for (unsigned int p=0; p<ptrange.size(); ++p) h_cons_bypt[p][process[dd]]->Write();
      }
    }//i
    work_dir->Close();
  }
  file_out->Close();


  return 0;
}
