#define IIHEAnalysis_cxx
#include "IIHEAnalysis.h"
#include "meta.h"
#include "PU_reWeighting.cc"
#include "zptweight.cc"
//#include <TH1.h>
#include <TLorentzVector.h>
//#include <TCanvas.h>
#include "TString.h"
#include "TDirectory.h"
//#include "MC_pileup_weight2017.cc"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  string out = *(argv + 1);
  string out_name= out;
  string in = *(argv + 2);
  string inname= in;
  string mc_in = *(argv + 3);
  string mc_nickname= mc_in;
  string phase_in = *(argv + 4);
  string phase= phase_in;
  string type_in = *(argv + 5);
  string type= type_in;
  TFile *fIn = TFile::Open(inname.c_str());
  TTree* tree = (TTree*) fIn->Get("IIHEAnalysis");


  TTree* mmeta = (TTree*) fIn->Get("meta");
  meta* m = new meta(mmeta);
  Float_t nEvents = m->Loop(type);
  //delete m;

  IIHEAnalysis* a = new IIHEAnalysis(tree);
  a->Loop(phase, type, inname, out_name, mc_nickname, nEvents);
  fIn->Close();
  //delete a;
  return 0;
}



//Get weighted events
Float_t meta::Loop(string type_of_data) {
  if (fChain == 0) return -1;

  bool data;
  if (type_of_data == "Data" || type_of_data == "data" || type_of_data == "singlephoton" || type_of_data == "SinglePhoton" || type_of_data == "singlemu" || type_of_data == "SingleMu") {
    data = true;
  }
  else {
    data = false;
  }

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  Float_t nEvents = -1;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if (data) {
      nEvents = nEventsRaw;
    }
    else {
      nEvents = mc_nEventsWeighted;
    }
  }
  return nEvents;
}



//Main analysis code
void IIHEAnalysis::Loop(string phase, string type_of_data, string in_name, string out_name, string mc_nickname, Float_t nEvents) {
   if (fChain == 0) return;


   fstream event_file; //object of fstream class
    
   string out_event = out_name;
   int point_position = -1;
   if (out_event.find('.') != string::npos) {
     point_position = out_event.find('.');
   }
   else {
     cout << "output name error!!!" << endl;
   }
   out_event.erase(point_position);
   out_event += ".txt";
   cout << out_event << endl;

   //opening file in out(write) mode
   event_file.open(out_event,ios::out);
   event_file << "root file: " <<  in_name << endl;
   event_file << "run, LS, event" << endl;


   bool DY, data;
   if (type_of_data == "DYinc" || type_of_data == "DYhighM" || type_of_data == "DY") {
     DY = true;         
     data = false;
   }
   else if (type_of_data == "Data" || type_of_data == "data") {
     DY = false;
     data = true;
   }
   else {
     DY = false;
     data = false;
   }

   bool QCD=false;
   if (type_of_data == "QCD") QCD=true;

   bool WJets=false;
   if (type_of_data == "WJets" || type_of_data == "Wjets" || type_of_data == "WJetsinc" || type_of_data == "Wjetsinc") WJets=true;

   bool TT=false;
   if (type_of_data == "TT" || type_of_data == "TTinc") TT=true;

   bool isWJetsphase, isQCDphase, isfinalphase, isAntiMuphase, isFakesphase, isWJetsFakephase;
   cout << phase << endl;
   cout << DY << endl;
   if (phase == "QCD") {
     isQCDphase = true;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "WJets") {
     isQCDphase = false;
     isWJetsphase = true;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "WJetsFake") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = true;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "final") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = true;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "AntiMu") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = true;
     isFakesphase = false;
   }
   else if (phase == "Fakes") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = true;
   }
   else {
     cout << phase << " does not exist!!!!" << endl;
   }

   //string out_name = "out_"+type_of_data+".root";
   TFile* file_out = new TFile(out_name.c_str(),"RECREATE");

   const float mu_mass = 0.10565837;


   //list here the names and x-axis ranges of all gen-level histos we wish to create :
   vector<TString> histo_gen_names;                vector<int> nBins_gen;     vector<float> x_min_gen,    x_max_gen; 
   histo_gen_names.push_back("gen_tauh_vispt");    nBins_gen.push_back(100);  x_min_gen.push_back(0);     x_max_gen.push_back(100);
   histo_gen_names.push_back("gen_tauh_viseta");   nBins_gen.push_back(50);   x_min_gen.push_back(-2.5);  x_max_gen.push_back(2.5);
   histo_gen_names.push_back("gen_tauh_visphi");   nBins_gen.push_back(64);   x_min_gen.push_back(-3.2);  x_max_gen.push_back(3.2);
   histo_gen_names.push_back("gen_taumu_vispt");   nBins_gen.push_back(100);  x_min_gen.push_back(0);     x_max_gen.push_back(100);
   histo_gen_names.push_back("gen_taumu_viseta");  nBins_gen.push_back(50);   x_min_gen.push_back(-2.5);  x_max_gen.push_back(2.5);
   histo_gen_names.push_back("gen_taumu_visphi");  nBins_gen.push_back(64);   x_min_gen.push_back(-3.2);  x_max_gen.push_back(3.2);
   histo_gen_names.push_back("gen_ev_Mvis");       nBins_gen.push_back(150);  x_min_gen.push_back(0);     x_max_gen.push_back(150);
   histo_gen_names.push_back("gen_tau_decaymode");    nBins_gen.push_back(3);  x_min_gen.push_back(-0.5);     x_max_gen.push_back(2.5);

   vector<TH1F*> hgen;
   for (unsigned int i = 0; i<histo_gen_names.size(); ++i) hgen.push_back( new TH1F(histo_gen_names[i], histo_gen_names[i], nBins_gen[i], x_min_gen[i], x_max_gen[i]) ); 


   //list here the names and x-axis ranges of all reco-level histos we wish to create :
   vector<TString> histo_notauID_names;             vector<int> nBins;     vector<float> x_min,   x_max; 
   histo_notauID_names.push_back("tau_MVA");        nBins.push_back(200);  x_min.push_back(-1);    x_max.push_back(1);
   histo_notauID_names.push_back("ev_DRmutau");     nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(10);
   histo_notauID_names.push_back("ev_Mt_raw");      nBins.push_back(150);  x_min.push_back(0);    x_max.push_back(150);
   histo_notauID_names.push_back("ev_Pzeta");       nBins.push_back(300);  x_min.push_back(-150); x_max.push_back(150);


   vector<TString> bkg_or_sig;
   bkg_or_sig.push_back("Bkg");          int iBkg = bkg_or_sig.size()-1;
   bkg_or_sig.push_back("Jet");          int iJet = bkg_or_sig.size()-1;
   if (DY) bkg_or_sig.push_back("Sig");  int iSig = bkg_or_sig.size()-1;
   int taun = bkg_or_sig.size();


   vector<TH1F*> hnotauID[taun];
   for (unsigned int i = 0; i<histo_notauID_names.size(); ++i) {
     for (unsigned int k = 0; k<taun; ++k) {
       hnotauID[k].push_back( new TH1F(bkg_or_sig[k]+"_"+histo_notauID_names[i], bkg_or_sig[k]+"_"+histo_notauID_names[i], nBins[i], x_min[i], x_max[i]) ); 
       hnotauID[k][i]->Sumw2();
     }
   }

   vector<TString> htau_names;                     vector<int> nBins_tau;     vector<float> x_min_tau,   x_max_tau;
   htau_names.push_back("ev_Mt");                  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis");                nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mtot");                nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt");                 nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_eta");                nBins_tau.push_back(50);   x_min_tau.push_back(-2.5); x_max_tau.push_back(2.5);
   htau_names.push_back("tau_phi");                nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   htau_names.push_back("tau_DM");                 nBins_tau.push_back(11);   x_min_tau.push_back(0);    x_max_tau.push_back(11);
   htau_names.push_back("mu_pt");                  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("mu_eta");                 nBins_tau.push_back(50);   x_min_tau.push_back(-2.5); x_max_tau.push_back(2.5);
   htau_names.push_back("mu_phi");                 nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   htau_names.push_back("ev_METmumass");           nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_MET");                 nBins_tau.push_back(100);  x_min_tau.push_back(0);    x_max_tau.push_back(100);
   htau_names.push_back("ev_METphi");              nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   htau_names.push_back("ev_Nvertex");             nBins_tau.push_back(120);  x_min_tau.push_back(0);    x_max_tau.push_back(120);
   htau_names.push_back("Z_pt");                   nBins_tau.push_back(300);  x_min_tau.push_back(0);    x_max_tau.push_back(300);
						   
   htau_names.push_back("ev_Mvis_TES_DM0_down");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_TES_DM0_up");     nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_TES_DM0_down");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_TES_DM0_up");      nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
						   
   htau_names.push_back("ev_Mvis_TES_DM1_down");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_TES_DM1_up");     nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_TES_DM1_down");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_TES_DM1_up");      nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
						   
   htau_names.push_back("ev_Mvis_TES_DM10_down");  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_TES_DM10_up");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_TES_DM10_down");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_TES_DM10_up");     nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
						   
   htau_names.push_back("ev_Mvis_MESdown");        nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_MESup");          nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_MESdown");         nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_MESup");           nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);

   htau_names.push_back("ev_Mvis_MinBiasdown");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_MinBiasup");      nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Nvertex_MinBiasdown"); nBins_tau.push_back(120);  x_min_tau.push_back(0);    x_max_tau.push_back(120);
   htau_names.push_back("ev_Nvertex_MinBiasup");   nBins_tau.push_back(120);  x_min_tau.push_back(0);    x_max_tau.push_back(120);

   htau_names.push_back("ev_Mvis_FRS_DM0_down");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_FRS_DM0_up");     nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_FRS_DM0_down");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_FRS_DM0_up");      nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_FRS_DM1_down");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_FRS_DM1_up");     nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_FRS_DM1_down");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_FRS_DM1_up");      nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_FRS_DM10_down");  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_FRS_DM10_up");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_FRS_DM10_down");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_FRS_DM10_up");     nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);

   htau_names.push_back("ev_Mvis_antiisomu_down"); nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_antiisomu_up");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_antiisomu_down");  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_antiisomu_up");    nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);

   htau_names.push_back("ev_Mvis_antiisotau_down");nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_antiisotau_up");  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_antiisotau_down"); nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt_antiisotau_up");   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);

   htau_names.push_back("tau_pt_SS");              nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis_SS");             nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("mu_pt_SS");               nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);

   //Make the map relating the histo name to its position on the histo name list
   map<TString, int> htau_map;
   for (unsigned int iName=0; iName<htau_names.size(); ++iName) htau_map[htau_names[iName]] = iName;


   vector<TString> HPS_WP;                  vector<int> specialValues;		     
   //HPS_WP.push_back("cutbased_loose");	    specialValues.push_back(HPS_WP.size()-1);
   //HPS_WP.push_back("cutbased_medium");	                                                  
   //HPS_WP.push_back("cutbased_tight");	                                                  
   //		                      	                                                  
   HPS_WP.push_back("MVA_2017v2vvloose");   specialValues.push_back(HPS_WP.size()-1);
   //HPS_WP.push_back("MVA_2017v2vloose");                                                  
   //HPS_WP.push_back("MVA_2017v2loose");	                                                  
   //HPS_WP.push_back("MVA_2017v2medium");                                                  
   HPS_WP.push_back("MVA_2017v2tight");	                                                  
   //HPS_WP.push_back("MVA_2017v2vtight");                                                  
   //HPS_WP.push_back("MVA_2017v2vvtight");                                                 
   //		                      	                                                  
   //HPS_WP.push_back("MVA_DBdR03vloose");    specialValues.push_back(HPS_WP.size()-1);
   //HPS_WP.push_back("MVA_DBdR03loose");	                                                  
   //HPS_WP.push_back("MVA_DBdR03medium");                                                  
   //HPS_WP.push_back("MVA_DBdR03tight");	                                                  
   //HPS_WP.push_back("MVA_DBdR03vtight");                                                  
   //HPS_WP.push_back("MVA_DBdR03vvtight");                                                 
   //		                      	                                                  
   //HPS_WP.push_back("MVA_PWdR03vloose");    specialValues.push_back(HPS_WP.size()-1);
   //HPS_WP.push_back("MVA_PWdR03loose");
   //HPS_WP.push_back("MVA_PWdR03medium");
   //HPS_WP.push_back("MVA_PWdR03tight");
   //HPS_WP.push_back("MVA_PWdR03vtight");
   //HPS_WP.push_back("MVA_PWdR03vvtight");
   int nTauIDs = HPS_WP.size();

   vector<TString> passfail;   
   passfail.push_back("pass");
   passfail.push_back("fail");
   map<TString, int> pass_map;
   for (unsigned int iPass=0; iPass<passfail.size(); ++iPass) pass_map[passfail[iPass]] = iPass;

   vector<TString> dms;
   dms.push_back("DM0");
   dms.push_back("DM1");
   dms.push_back("DM10");
   map<TString, int> dms_map;
   for (unsigned int iDM=0; iDM<dms.size(); ++iDM) dms_map[dms[iDM]] = iDM;

   vector<TString> eta;
   eta.push_back("barrel");
   eta.push_back("endcap");
   map<TString, int> eta_map;
   for (unsigned int iEta=0; iEta<eta.size(); ++iEta) eta_map[eta[iEta]] = iEta;

   vector<TString> ptrange;
   ptrange.push_back("pt_20_40");
   ptrange.push_back("pt_40_150");
   map<TString, int> ptrange_map;
   for (unsigned int iPt=0; iPt<ptrange.size(); ++iPt) ptrange_map[ptrange[iPt]] = iPt;


   vector<TH1F*> htau[ptrange.size()][eta.size()][dms.size()][passfail.size()][taun][nTauIDs];
   for (unsigned int i = 0; i<htau_names.size(); ++i) {
     for (unsigned int j = 0; j<nTauIDs ; ++j) {
       for (unsigned int k = 0; k<taun; ++k) {
	 for (unsigned int l = 0; l<passfail.size(); ++l) {
	   for (unsigned int m = 0; m<dms.size(); ++m) {
	     for (unsigned int n = 0; n<eta.size(); ++n) {
	       for (unsigned int p = 0; p<ptrange.size(); ++p) {
		 TString nname_h = bkg_or_sig[k]+"_"+htau_names[i]+"_"+dms[m]+"_"+eta[n]+"_"+ptrange[p]+"_"+HPS_WP[j]+"_"+passfail[l];
		 htau[p][n][m][l][k][j].push_back( new TH1F(nname_h, nname_h, nBins_tau[i], x_min_tau[i], x_max_tau[i]) ); 
		 htau[p][n][m][l][k][j][i]->Sumw2();
	       }
	     }
	   }
	 }
       }
     }
   }


   TH1F* h_reweight = new TH1F("h_reweight", "h_reweight", 100, 0, 2.5);
   TH1F* h_weight_afterallsel = new TH1F("h_weight_afterallsel", "h_weight_afterallsel", 100, -2, 2);
   TH1F* h_debug = new TH1F("h_debug", "h_debug", 100, 0, 5);
   TH1F* h_genmass = new TH1F("h_genmass", "h_genmass", 120, 0, 120);
   TH1F* h_genpt = new TH1F("h_genpt", "h_genpt", 120, 0, 120);

   readZptFile(2017);

   Long64_t nEntries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int print_count = 0;
   long n1=0, n2=0;
   //start loop over all events
   for (Long64_t jEntry = 0; jEntry < nEntries; ++jEntry) {
      Long64_t iEntry = LoadTree(jEntry);
      if (iEntry < 0) break;
      if (jEntry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", jEntry, nEntries);
      //if (jEntry % 10 == 0) cout << endl << "Processed events:" << jEntry << " of " <<  nEntries;
      
      nb = fChain->GetEntry(jEntry);
      nbytes += nb;


      float pu_weight = 1, pu_weight_high = 1, pu_weight_low = 1;
      if (!data) {
        pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Data_METcorr_2017BtoF");
        pu_weight_high = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Data_METcorr_2017BtoF_high");
        pu_weight_low = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Data_METcorr_2017BtoF_low");
      }
      float first_weight = pu_weight;
      float reweight_njets = 1.0;
      float zpt_weight = 1.0;
      
      bool gen_mumu = false;
      vector<TLorentzVector> gentau_had_visp4, genmu_p4, anylepton_p4;
      gentau_had_visp4.clear(), genmu_p4.clear(), anylepton_p4.clear();
      if (!data) {
	
	//find Zpt weight
	if (DY) {
	  TLorentzVector p4, l1_p4, l2_p4, z_p4;
	  p4.SetPxPyPzE(0, 0, 0, 0);
	  l1_p4.SetPxPyPzE(0, 0, 0, 0);
	  l2_p4.SetPxPyPzE(0, 0, 0, 0);
	  z_p4.SetPxPyPzE(0, 0, 0, 0);
	  int l1_pdgid = -1, l2_pdgid = -1;
	  for (unsigned int iLHE = 0; iLHE < LHE_Pt->size(); ++iLHE) {
	    if (LHE_pdgid->at(iLHE) == 11 || LHE_pdgid->at(iLHE) == 13 || LHE_pdgid->at(iLHE) == 15) {
	      l1_pdgid = LHE_pdgid->at(iLHE);
	    }
	    else if (LHE_pdgid->at(iLHE) == -11 || LHE_pdgid->at(iLHE) == -13 || LHE_pdgid->at(iLHE) == -15) {
	      l2_pdgid = LHE_pdgid->at(iLHE);
	    }
	  }

	  bool found_1 = false, found_2 = false;
	  for (unsigned int iMC = 0; iMC < mc_pt->size(); ++iMC) {
	    p4.SetPtEtaPhiE(mc_pt->at(iMC),mc_eta->at(iMC),mc_phi->at(iMC),mc_energy->at(iMC));
	    bool isfromHardProcess = checkBit (mc_status_flags->at(iMC),8);
	    bool isGoodLep = (abs(mc_pdgId->at(iMC)) == 11 || abs(mc_pdgId->at(iMC)) == 13) && (mc_status->at(iMC) == 1 || mc_status->at(iMC) == 23);
	    bool isGoodTau = abs(mc_pdgId->at(iMC)) == 15 && mc_status->at(iMC) == 2;
	    bool isfromHardProcessFinalState = ( isGoodLep || isGoodTau ) && isfromHardProcess;

	    if (mc_pdgId->at(iMC) == l1_pdgid) {
	      if(!isfromHardProcessFinalState) continue;
	      l1_p4 = p4;
	      found_1 = true;
	    }
	    else if (mc_pdgId->at(iMC) == l2_pdgid) {
	      if(!isfromHardProcessFinalState) continue;
	      l2_p4 = p4;
	      found_2 = true;
	    }
	  }
	  if (found_1 && found_2) {
	    z_p4 = l1_p4 + l2_p4;
	    h_genmass->Fill(z_p4.M());
	    h_genpt->Fill(z_p4.Pt());
	    zpt_weight = getZpt(z_p4.M(), z_p4.Pt());
	  }
	}


	TLorentzVector p4, p4_2;
	//start loop over all simulated particules
	for (unsigned int iMC = 0; iMC < mc_tau_had_pt->size(); ++iMC) {
	  p4.SetPtEtaPhiE(mc_tau_had_pt->at(iMC), mc_tau_had_eta->at(iMC), mc_tau_had_phi->at(iMC), mc_tau_had_energy->at(iMC));
	  gentau_had_visp4.push_back( p4 );
	  anylepton_p4.push_back( p4 );
	}//1st loop over sim particules

	
	for (unsigned int iMC = 0; iMC < mc_pt->size(); ++iMC) {
	  if (abs(mc_pdgId->at(iMC)) == 11 || abs(mc_pdgId->at(iMC)) == 13) {
	    if (mc_pt->at(iMC) < 10) continue;
	    p4.SetPxPyPzE(mc_px->at(iMC), mc_py->at(iMC), mc_pz->at(iMC), mc_energy->at(iMC));
	    if (p4.Pt() < 10) continue;
	    if (p4.Pt() > 10000) continue;
	    anylepton_p4.push_back(p4);
	  }
	}

	
	for (unsigned int iMC1 = 0; iMC1 < mc_pt->size(); ++iMC1) {
	  if (abs(mc_pdgId->at(iMC1)) != 13) continue;
	  p4.SetPtEtaPhiE(mc_pt->at(iMC1), mc_eta->at(iMC1), mc_phi->at(iMC1), mc_energy->at(iMC1));
	  genmu_p4.push_back(p4);
	  for (unsigned int iMC2 = 0; iMC2 < iMC1; ++iMC2) {
	    if (abs(mc_pdgId->at(iMC2)) != 13) continue;
	    p4_2.SetPtEtaPhiE(mc_pt->at(iMC2), mc_eta->at(iMC2), mc_phi->at(iMC2), mc_energy->at(iMC2));
	    if (p4.DeltaR(p4_2) > 0.5) {
	      gen_mumu = true;
	      break;
	    }
	    if (gen_mumu) break;
	  }
	}

      }//end is this MC? condition


      if (WJets) {
	int njets = -2;
	for (unsigned int iLHE = 0; iLHE < LHE_Pt->size(); ++iLHE) {
          if (abs(LHE_pdgid->at(iLHE)) < 10 || abs(LHE_pdgid->at(iLHE)) == 21) ++njets;
	  //cout << LHE_pdgid->at(iLHE) << endl;
        }
	//cout << "njets  " << njets << endl;
	if (njets==0) {
	  reweight_njets = 111.024769762929;
	}
	else if (njets==1) {
	  reweight_njets = 14.68242648;
	}
	else if (njets==2) {
	  reweight_njets = 7.528580277;
	}
	else if (njets==3) {
	  reweight_njets = 2.436265125;
	}
	else if (njets==4) {
	  reweight_njets = 1.0599498;
	}
      }//end is this WJets? condition


      //Is one of the triggers fired?
      bool PassMuonTrigger = false;
      if (trig_HLT_IsoMu27_accept) PassMuonTrigger = true;
      if (!PassMuonTrigger) continue;

      if(!trig_Flag_goodVertices_accept) continue;
      if(!trig_Flag_globalSuperTightHalo2016Filter_accept) continue;
      if(!trig_Flag_HBHENoiseFilter_accept) continue;
      if(!trig_Flag_HBHENoiseIsoFilter_accept) continue;
      if(!trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept) continue;
      if(!trig_Flag_BadPFMuonFilter_accept) continue;
      //FIXME: add Ecal bad calibration when it's available


      //start muon counting loop
      int Nmu = 0;
      bool dimuon = false;
      for (unsigned int iMu = 0; iMu < mu_gt_pt->size(); ++iMu) {
	TLorentzVector mu1_p4, mu2_p4;
        /*if(mu_isPFMuon->at(iMu) && mu_gt_pt->at(iMu) > 10 && fabs(mu_gt_eta->at(iMu)) < 2.4 && fabs(mu_gt_dxy_firstPVtx->at(iMu)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu) < 0.3 && mu_isMediumMuon->at(iMu)) ++Nmu;
	  if (Nmu > 1) break;*/
	
	if (mu_gt_pt->at(iMu) > 15 && fabs(mu_gt_eta->at(iMu)) < 2.4 && mu_isGlobalMuon->at(iMu) && mu_isTrackerMuon->at(iMu) && mu_isPFMuon->at(iMu) &&fabs(mu_gt_dxy_firstPVtx->at(iMu)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu) < 0.3) {
	  mu1_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu), mu_gt_eta->at(iMu), mu_gt_phi->at(iMu), mu_mass);
	  for (unsigned int iMu2= 0; iMu2 < iMu; ++iMu2) {
	    if (mu_gt_pt->at(iMu2) > 15 && fabs(mu_gt_eta->at(iMu2)) < 2.4 && mu_isGlobalMuon->at(iMu2) && mu_isTrackerMuon->at(iMu2) && mu_isPFMuon->at(iMu2) &&fabs(mu_gt_dxy_firstPVtx->at(iMu2)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu2)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu2) < 0.3) {
	      mu2_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu2), mu_gt_eta->at(iMu2), mu_gt_phi->at(iMu2), mu_mass);
	      if (mu1_p4.DeltaR(mu2_p4) > 0.15) dimuon = true;
	    }
	  }
	}
	if (dimuon) break;
      }
      //if (Nmu > 1) continue; //2nd muon veto                             
      if (dimuon) continue;

      //electron veto
      bool electron = false;
      for (unsigned int iEle = 0; iEle < gsf_pt->size(); ++iEle) {
      if (gsf_VIDMVAMedium->at(iEle) && gsf_pt->at(iEle) > 10 && fabs(gsf_eta->at(iEle)) < 2.5 && fabs(gsf_dxy_firstPVtx->at(iEle)) < 0.045 && fabs(gsf_dz_firstPVtx->at(iEle)) < 0.2 && gsf_passConversionVeto->at(iEle) && gsf_nLostInnerHits->at(iEle) <= 1 && gsf_relIso->at(iEle) < 0.3) electron = true;
      if (electron) break;
      }
      if (electron) continue;//FIXME

      //bjet veto (medium WP for the bjet)                                                                                                                           
      /*bool bjet = false;
      for (unsigned int iJet = 0; iJet < jet_pt->size(); ++iJet) {
        if (jet_CSVv2->at(iJet) > 0.800 && jet_pt->at(iJet) > 20 && fabs(jet_eta->at(iJet)) < 2.4) bjet = true;
        if (bjet) break;
      }
      if (bjet) continue;*/



      //Sort muons, taus, by increasing isolation/decreasing pt
      float pt = -1000;
      int lowest = -1;
      vector<int> orderedMu, orderedTau;
      vector<int> rest, rest2;


      //sorting muons
      for (unsigned int ii = 0; ii < mu_gt_pt->size(); ++ii) {
	rest.push_back(ii);
      }
      while (rest.size()>0) {
	rest2.clear();
	lowest = -1;
	pt = -1000;
	for (unsigned int ii = 0; ii < rest.size(); ++ii) {
	  if (mu_gt_pt->at(rest[ii]) > pt) {
	    pt = mu_gt_pt->at(rest[ii]);
	    if (lowest > -1) rest2.push_back(lowest);
	    lowest = rest[ii];
	  }
	  else {
	    rest2.push_back(rest[ii]);
	  }
	}
	orderedMu.push_back(lowest);
	rest = rest2;
      }


      //sorting Taus
      rest.clear();
      for (unsigned int ii = 0; ii < tau_pt->size(); ++ii) {
	rest.push_back(ii);
      }
      while (rest.size()>0) {
	rest2.clear();
	lowest = -1;
	pt = -1000;
	for (unsigned int ii = 0; ii < rest.size(); ++ii) {
	  if (tau_pt->at(rest[ii]) > pt) {
	    pt = tau_pt->at(rest[ii]);
	    if (lowest > -1) rest2.push_back(lowest);
	    lowest = rest[ii];
	  }
	  else {
	    rest2.push_back(rest[ii]);
	  }
	}
	orderedTau.push_back(lowest);
	rest = rest2;
      }



      //start loop over reconstructed muons
      bool found_mutau_pair = false;
      for (unsigned int ii = 0; ii < orderedMu.size(); ++ii) {
	if (found_mutau_pair) break;
	int iMu = orderedMu[ii];
	if (mu_gt_pt->at(iMu) < 30.0) continue;
	if (fabs(mu_gt_eta->at(iMu)) > 2.4) continue;
	if (!mu_isPFMuon->at(iMu)) continue;
	if (!mu_isMediumMuon->at(iMu)) continue;
	if (fabs(mu_gt_dxy_firstPVtx->at(iMu)) >= 0.045) continue;
	if (fabs(mu_gt_dz_firstPVtx->at(iMu)) >= 0.2) continue;
	float reliso = mu_pfIsoDbCorrected04->at(iMu);
	if (isAntiMuphase) {
	  if (reliso < 0.15) continue;
	}
	else {
	  if (reliso > 0.15) continue;
	}


	TLorentzVector mu_p4;
	mu_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu), mu_gt_eta->at(iMu), mu_gt_phi->at(iMu), mu_mass);
	//start loop over reconstructed taus
	for (unsigned int jj = 0; jj < orderedTau.size(); ++jj) {
	  if (found_mutau_pair) break;

	  int iTau = orderedTau[jj];
	  //Move selections as far up as possible to speed things up...
	  //cuts which do not involve the tau pt or the MET
	  if (fabs(tau_eta->at(iTau)) > 2.3) continue;
	  if (tau_decayModeFinding->at(iTau) < 0.5) continue;
	  if (tau_againstMuonTight3->at(iTau) < 0.5) continue;
	  if (tau_againstElectronVLooseMVA6->at(iTau) < 0.5) continue;
	  if (fabs(tau_lead_dz->at(iTau)) > 0.2) continue;
	  if (fabs(tau_charge->at(iTau)) != 1) continue;

	
	  //misc. cuts
	  bool SS = false;
	  if (isfinalphase || isWJetsphase || isWJetsFakephase || isFakesphase) {
	    if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) > 0) continue; //SS veto (opposite for QCD estimation)
	    SS = false;
	  }
	  else if (isQCDphase) {
	    if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) < 0) continue; //OS veto
	    SS = true;
	  }
	  else {
	    if (isAntiMuphase) {
	      SS = false;
	      if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) > 0) SS = true;
	    }
	    else {
	      cout << "Error : non-existent phase!!!!" << endl;
	      break;
	    }
	  }


	  //Make the map relating the HPS WP name to its value (true/false)
	  map<TString, int> tauIDvalues_map;
	                                                                                                       
	  tauIDvalues_map["cutbased_loose"] = tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(iTau);        
	  tauIDvalues_map["cutbased_medium"] = tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(iTau);
	  tauIDvalues_map["cutbased_tight"] = tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(iTau);
			                    
	  tauIDvalues_map["MVA_2017v2vvloose"] = tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017->at(iTau);   
	  tauIDvalues_map["MVA_2017v2vloose"] = tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017->at(iTau);
	  tauIDvalues_map["MVA_2017v2loose"] = tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017->at(iTau);
	  tauIDvalues_map["MVA_2017v2medium"] = tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017->at(iTau);
	  tauIDvalues_map["MVA_2017v2tight"] = tau_byTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau);
	  tauIDvalues_map["MVA_2017v2vtight"] = tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau);
	  tauIDvalues_map["MVA_2017v2vvtight"] = tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau);
			                    
	  tauIDvalues_map["MVA_DBdR03vloose"] = tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->at(iTau);        
	  tauIDvalues_map["MVA_DBdR03loose"] = tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_DBdR03medium"] = tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_DBdR03tight"] = tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_DBdR03vtight"] = tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_DBdR03vvtight"] = tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT->at(iTau);
			                    
	  tauIDvalues_map["MVA_PWdR03vloose"] = tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT->at(iTau);        
	  tauIDvalues_map["MVA_PWdR03loose"] = tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_PWdR03medium"] = tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_PWdR03tight"] = tau_byTightIsolationMVArun2v1PWdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_PWdR03vtight"] = tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT->at(iTau);
	  tauIDvalues_map["MVA_PWdR03vvtight"] = tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT->at(iTau);



	  //check if if it's DY bkg or signal
	  TLorentzVector tau_p4, met_p4;
	  float met_px = MET_eefix_Px;
	  float met_py = MET_eefix_Py;
	  float met_pt = MET_eefix_Pt;
	  tau_p4.SetPtEtaPhiE(tau_pt->at(iTau), tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau));
	  met_p4.SetPxPyPzE(met_px, met_py, 0, met_pt);

	  bool DY_sig = false;
	  bool fakemu_match = false;
	  float fakemu_reweight = 1;
	  if (DY) {
	    for (unsigned int iGen = 0; iGen < gentau_had_visp4.size(); ++iGen) {
	      if (tau_p4.DeltaR(gentau_had_visp4[iGen]) < 0.5) {
		DY_sig = true;
		break;
	      }
	    }
	    if (!DY_sig) {
	      //reweighting for DY->mumu
	      if (gen_mumu) {
		for (unsigned int iGen = 0; iGen < genmu_p4.size(); ++iGen) {
		  if (tau_p4.DeltaR(genmu_p4[iGen]) < 0.5) {
		    fakemu_match = true;
		    break;
		  }
		}
		if (!fakemu_match) continue;
		
		if (fabs(tau_p4.Eta())<0.4) fakemu_reweight = 1.200;
		else if (fabs(tau_p4.Eta())<0.8) fakemu_reweight = 1.4;
		else if (fabs(tau_p4.Eta())<1.2) fakemu_reweight = 1.09;
		else if (fabs(tau_p4.Eta())<1.7) fakemu_reweight = 0.9;
		else if (fabs(tau_p4.Eta())<2.3) fakemu_reweight = 2.3;
	      }
	      //cout << tau_p4.Eta() << "  " << first_weight << "  " << final_weight << "  " << final_weight/first_weight << endl;
	    }
	  }
	  //final_weight = first_weight;
	  
	  
	  bool realtau = false;
	  if (!data) {
	    for (unsigned int iGen = 0; iGen < gentau_had_visp4.size(); ++iGen) {
	      if (tau_p4.DeltaR(gentau_had_visp4[iGen]) < 0.5) {
		realtau = true;
		break;
	      }
	    }
	  }
	  
	  bool is_not_jet = false;
	  if (anylepton_p4.size() != 0) {
	    for (unsigned int iGen = 0; iGen < anylepton_p4.size(); ++iGen) {
	      if (tau_p4.DeltaR(anylepton_p4[iGen]) < 0.5) {
		is_not_jet = true;
		break;
	      }
	    }
	  }

	  
	  int iBkgOrSig = -1;
	  if (DY_sig) {
	    iBkgOrSig = iSig;
	  }
	  else {
	    if (is_not_jet) {
	      iBkgOrSig = iBkg;
	    }
	    else {
	      iBkgOrSig = iJet;
	    }
	  }

	  for (int iTES = 0; iTES < 3; ++iTES) {
	    int iTES_down = 0; //0 : TES down
	    int iTES_up = 1; //1 : TES up
	    int iTES_nope = 2; //2 : no TES
	    //only vary the tau energy scale if it's a real tau
	    if (!realtau && iTES != iTES_nope) continue;
	    
	    map <TString, float> map_tau_TES;  map <TString, float> map_tau_TES_error;
	    map_tau_TES["DM0"] = 1.007;	       map_tau_TES_error["DM0"] = 0.008;	     
	    map_tau_TES["DM1"] = 0.998;	       map_tau_TES_error["DM1"] = 0.008;	     
	    map_tau_TES["DM10"] = 1.001;       map_tau_TES_error["DM10"] = 0.009;     

	    TLorentzVector tau_TES_p4, tau_MES_p4, vis_p4, total_p4;
	    //met_p4.SetPxPyPzE(MET_Px, MET_eefix_Py, 0, MET_eefix_Pt);


	    for (int iMES = 0; iMES < 3; ++iMES) {
	      int iMES_down = 0; //0 : MES down
	      int iMES_up = 1; //1 : MES up
	      int iMES_nope = 2; //2 : no MES
	      if (iTES !=iTES_nope && iMES !=iMES_nope) continue;
	      if (!fakemu_match && iMES != iMES_nope) continue;

	      

	      //determine decay mode and eta
	      TString str_dm = "", str_eta = "";
	      int m_dm = -1, n_eta = -1;
	      if (tau_decayMode->at(iTau) == 0) {
		str_dm = "DM0";
	      }
	      else if (tau_decayMode->at(iTau) == 1) {
		str_dm = "DM1";
	      }
	      else if (tau_decayMode->at(iTau) == 10) {
		str_dm = "DM10";
	      }
	      else {
		continue;
	      }
	      m_dm = dms_map[str_dm];

	      if (fabs(tau_eta->at(iTau)) < 1.5) {
		str_eta = "barrel";
	      }
	      else {
		str_eta = "endcap";
	      }
	      n_eta = eta_map[str_eta];


	      tau_p4.SetPtEtaPhiE(tau_pt->at(iTau), tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau));
	      if (realtau) {
		//correct tau energy scale central value
		tau_TES_p4.SetPtEtaPhiE(tau_pt->at(iTau)*map_tau_TES[str_dm], tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau)*map_tau_TES[str_dm]);
		met_p4 = met_p4 + tau_p4 - tau_TES_p4;
		tau_p4 = tau_TES_p4;

		//now vary tau energy scale down and up according to the uncertainty
		if (iTES == iTES_down) {
		  //TES down
		  tau_TES_p4.SetPtEtaPhiE(tau_p4.Pt()*(1-map_tau_TES_error[str_dm]), tau_p4.Eta(), tau_p4.Phi(), tau_p4.E()*(1-map_tau_TES_error[str_dm]));
		  met_p4 = met_p4 + tau_p4 - tau_TES_p4;
		  tau_p4 = tau_TES_p4;
	        }
		else if (iTES == iTES_up) {
		  //TES up
		  tau_TES_p4.SetPtEtaPhiE(tau_p4.Pt()*(1+map_tau_TES_error[str_dm]), tau_p4.Eta(), tau_p4.Phi(), tau_p4.E()*(1+map_tau_TES_error[str_dm]));
		  met_p4 = met_p4 + tau_p4 - tau_TES_p4;
		  tau_p4 = tau_TES_p4;
	        }
	      }


	      if (iMES == iMES_down) {
	        if (fakemu_match) {
		  //MES down by 3%
		  tau_MES_p4.SetPtEtaPhiE(tau_pt->at(iTau)*0.97, tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau)*0.97);
		  met_p4 = met_p4 + tau_p4 - tau_MES_p4;
		  tau_p4 = tau_MES_p4;
	        }
	      }
	      else if (iMES == iMES_up) {
	        if (fakemu_match) {
		  //MES up by 3%
		  tau_MES_p4.SetPtEtaPhiE(tau_pt->at(iTau)*1.03, tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau)*1.03);
		  met_p4 = met_p4 + tau_p4 - tau_MES_p4;
		  tau_p4 = tau_MES_p4;
	        }
	      }


	      vis_p4 = tau_p4 + mu_p4;
	      total_p4 = vis_p4;// + met_p4;
	      
	      if (tau_p4.Pt() < 20.0) continue;
	      
	      int p_pt = -1;
	      if (tau_p4.Pt() >= 20 && tau_p4.Pt() < 40) {
		p_pt = ptrange_map["pt_20_40"];
	      }
	      else if (tau_p4.Pt() >= 40) {
		p_pt = ptrange_map["pt_40_150"];
	      }

	      
	      //Pzeta calculation
	      float norm_zeta= norm_F( tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt(), tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt() );
	      //cout << norm_zeta << endl;
	      float x_zeta= (tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt())/norm_zeta;
	      float y_zeta= (tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt())/norm_zeta;
	      float p_zeta_mis=met_p4.Px()*x_zeta+met_p4.Py()*y_zeta;
	      float pzeta_vis=(tau_p4.Px()+mu_p4.Px())*x_zeta+(tau_p4.Py()+mu_p4.Py())*y_zeta;
	      //bool cut_zeta= p_zeta_mis-0.85*pzeta_vis>-25;
	      bool cut_zeta= p_zeta_mis+0.15*pzeta_vis>-25;
	      
	      

	      if (first_weight != first_weight) continue;
	      
	      float Mt;
	      bool Mt_accept = true;
	      if (2 * ( mu_gt_pt->at(iMu) * met_p4.Pt()  - mu_gt_px->at(iMu)*met_p4.Px() - mu_gt_py->at(iMu)*met_p4.Py() ) < 0) {
	        Mt = 0;
	        Mt_accept = false;
	      }
	      else {
	        Mt = sqrt(2 * ( mu_gt_pt->at(iMu) * met_p4.Pt()  - mu_gt_px->at(iMu)*met_p4.Px() - mu_gt_py->at(iMu)*met_p4.Py() ));
	      }
	      
	      
	      if (isfinalphase || isQCDphase || isAntiMuphase || isFakesphase) {
	        if (Mt > 50) Mt_accept = false;  //Mt cut, against Wjets (for Wjets CR, Mt>80)
	      }
	      else if (isWJetsphase || isWJetsFakephase) {
	        if (Mt < 80) Mt_accept = false;  //Mt cut, against Wjets (for Wjets CR, Mt>80)
	      }	      
	      else {
	        if (!isAntiMuphase) {
	      	cout << "Error : non-existent phase!!!!" << endl;
	      	break;
	        }
	      }	      
	      if (iTES==iTES_nope && Mt_accept) h_debug->Fill(2);
	      
	      float dR = tau_p4.DeltaR(mu_p4);
	      
	      float mc_weight = 1;
	      if (!data) mc_weight = mc_w_sign;
	      
	      float other_weights = 1;
	      if (!data) other_weights = GetTriggerMuonIDMuonIsoReweight(mu_p4.Pt(), mu_p4.Eta());
	      reweight_njets = 1;
	      float final_weight = pu_weight*other_weights*reweight_njets*fakemu_reweight*mc_weight*zpt_weight;
	      first_weight = final_weight;

	      if (Mt_accept && iTES == iTES_nope) hnotauID[iBkgOrSig][1]->Fill(dR, final_weight);
	      if (Mt_accept && iTES == iTES_nope) hnotauID[iBkgOrSig][3]->Fill(p_zeta_mis+0.15*pzeta_vis, final_weight);
	      //if (Mt_accept && iTES == 2) hnotauID[0][6]->Fill(x_zeta*x_zeta + y_zeta*y_zeta, final_weight);
	      
	      if (tau_p4.DeltaR(mu_p4) < 0.5) continue;
	      if (iTES==iTES_nope && Mt_accept) h_debug->Fill(3);
	      if (iTES == iTES_nope) hnotauID[iBkgOrSig][2]->Fill(Mt, final_weight);
	      if (!cut_zeta) continue;
	      
	      
	      
	      if(Mt_accept && iTES == iTES_nope && iMES == iMES_nope) {
	        found_mutau_pair = true;
	      }

	      
	      //mu histos
	      
	      if (Mt_accept && iTES == iTES_nope) {
		//tau MVA histo
		hnotauID[iBkgOrSig][0]->Fill(tau_byIsolationMVArun2v1DBoldDMwLTraw->at(iTau), final_weight);
	      }
	      
	      
	      int latestloosestID = -1, s_counter = 0;
	      float fakerate_weight = 1, fakerate_weight_low = 1, fakerate_weight_high = 1;
	      for (unsigned int iValue=0; iValue<HPS_WP.size(); ++iValue) {
	        int iPass = -1;

	        if (tauIDvalues_map[HPS_WP[iValue]] > 0.5) {
		  iPass = pass_map["pass"];
	        }
	        else {
		  iPass = pass_map["fail"];
	        }

	      
	        //fake rate method
	        if (isFakesphase || isWJetsFakephase) {
	      	  fakerate_weight = FakeRate(tau_p4.Pt(), HPS_WP[iValue], str_dm, str_eta);
	      	  fakerate_weight_high = FakeRateFlat(HPS_WP[iValue], str_dm);
		  fakerate_weight_low = 2*fakerate_weight - fakerate_weight_high;
	      	  
	      	  bool special = false;
		  if (iValue == specialValues[s_counter]) {
		    special = true;
		    latestloosestID = iValue;
		    ++s_counter;
		  }
		  
	      	  if (!special) {
	      	    if (tauIDvalues_map[HPS_WP[latestloosestID]] < 0.5) continue;
	      	  }
		}

		final_weight = first_weight*fakerate_weight;
		h_reweight->Fill(final_weight);
		if (final_weight != final_weight) continue;
	      
		if (iTES == iTES_nope) {
		  if (Mt_accept && (iMES == iMES_down || !fakemu_match)) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_MESdown"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept && (iMES == iMES_up || !fakemu_match))   htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_MESup"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept && (iMES == iMES_down || !fakemu_match)) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_MESdown"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept && (iMES == iMES_up || !fakemu_match))   htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_MESup"]]->Fill(tau_p4.Pt(), final_weight);
		}
		if (iMES != iMES_nope) continue;
	      
		//fill variations for the DM present in this event
	        if (Mt_accept && (iTES == iTES_down || !realtau)) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_TES_"+str_dm+"_down"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && (iTES == iTES_up || !realtau))   htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_TES_"+str_dm+"_up"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && (iTES == iTES_down || !realtau)) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_TES_"+str_dm+"_down"]]->Fill(tau_p4.Pt(), final_weight);
	        if (Mt_accept && (iTES == iTES_up || !realtau))   htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_TES_"+str_dm+"_up"]]->Fill(tau_p4.Pt(), final_weight);
		if (iTES != iTES_nope) continue;
		
		//iterate over other decay modes and fill the unmodified event
		for (map<TString,int>::iterator it=dms_map.begin(); it!=dms_map.end(); ++it) {
		  if (it->first == str_dm) continue;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_TES_"+it->first+"_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept)   htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_TES_"+it->first+"_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_TES_"+it->first+"_down"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept)   htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_TES_"+it->first+"_up"]]->Fill(tau_p4.Pt(), final_weight);
		}


		//fill text file with the run number etc.
		if (Mt_accept && iValue==0) event_file << ev_run << ", " << ev_luminosityBlock << ", " << ev_event << endl; 
		
		//normal histos
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mt"]]->Fill(Mt, final_weight);
	        if (!SS) if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis"]]->Fill(vis_p4.M(), final_weight);
	        if (SS) if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_SS"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mtot"]]->Fill(total_p4.M(), final_weight);
	        if (!SS) if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt"]]->Fill(tau_p4.Pt(), final_weight);
	        if (SS) if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_SS"]]->Fill(tau_p4.Pt(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_eta"]]->Fill(tau_p4.Eta(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_phi"]]->Fill(tau_p4.Phi(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_DM"]]->Fill(tau_decayMode->at(iTau), final_weight);
	        if (!SS) if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["mu_pt"]]->Fill(mu_gt_pt->at(iMu), final_weight);
	        if (SS) if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["mu_pt_SS"]]->Fill(mu_gt_pt->at(iMu), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["mu_eta"]]->Fill(mu_gt_eta->at(iMu), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["mu_phi"]]->Fill(mu_gt_phi->at(iMu), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["Z_pt"]]->Fill(total_p4.Pt(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_MET"]]->Fill(met_p4.Pt(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_METphi"]]->Fill(met_p4.Phi(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Nvertex"]]->Fill(pv_n, final_weight);

		//Min Bias variations, change pu weight
		final_weight = first_weight*fakerate_weight;
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_MinBiasdown"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Nvertex_MinBiasdown"]]->Fill(pv_n, final_weight);
		final_weight = first_weight*fakerate_weight;
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_MinBiasup"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Nvertex_MinBiasup"]]->Fill(pv_n, final_weight);


		//anti-isolation SFs, which could be different from one
		//tau is anti-iso in Fakes phase, tau could be either matched to a gen muon or a gen tau, or to some other object (mostly jets)
		double antiisomu_weight_low = 0.75;
		double antiisomu_weight_high = 1.25;
		if (fakemu_match) {
		  final_weight = first_weight*fakerate_weight*antiisomu_weight_low;
		}
		else {
		  final_weight = first_weight*fakerate_weight*1.0;
		}		  
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_antiisomu_down"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_antiisomu_down"]]->Fill(tau_p4.Pt(), final_weight);

		if (fakemu_match) {
		  final_weight = first_weight*fakerate_weight*antiisomu_weight_high;
		}
		else {
		  final_weight = first_weight*fakerate_weight*1.0;
		}		  
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_antiisomu_up"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_antiisomu_up"]]->Fill(tau_p4.Pt(), final_weight);


		double antiisotau_weight_low = 0.9;
		double antiisotau_weight_high = 1.1;
		if (realtau) {
		  final_weight = first_weight*fakerate_weight*antiisotau_weight_low;
		}
		else {
		  final_weight = first_weight*fakerate_weight*1.0;
		}		  
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_antiisotau_down"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_antiisotau_down"]]->Fill(tau_p4.Pt(), final_weight);

		if (realtau) {
		  final_weight = first_weight*fakerate_weight*antiisotau_weight_high;
		}
		else {
		  final_weight = first_weight*fakerate_weight*1.0;
		}		  
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_antiisotau_up"]]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_antiisotau_up"]]->Fill(tau_p4.Pt(), final_weight);




		if (!isFakesphase) continue;
		//fake rate variations
		if (tau_decayMode->at(iTau) == 0) {
		  //for the histos dedicated to DM0 variation, fill with a different weight
		  final_weight = first_weight*fakerate_weight_low;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM0_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM0_down"]]->Fill(tau_p4.Pt(), final_weight);

		  final_weight = first_weight*fakerate_weight_high;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM0_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM0_up"]]->Fill(tau_p4.Pt(), final_weight);

		  //other histos should be filled with the "normal" weight
		  final_weight = first_weight*fakerate_weight;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM1_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM1_down"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM1_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM1_up"]]->Fill(tau_p4.Pt(), final_weight);

		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_down"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_up"]]->Fill(tau_p4.Pt(), final_weight);
		}
		else if (tau_decayMode->at(iTau) == 1) {
		  //for the histos dedicated to DM1 variation, fill with a different weight
		  final_weight = first_weight*fakerate_weight_low;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM1_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM1_down"]]->Fill(tau_p4.Pt(), final_weight);

		  final_weight = first_weight*fakerate_weight_high;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM1_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM1_up"]]->Fill(tau_p4.Pt(), final_weight);

		  //other histos should be filled with the "normal" weight
		  final_weight = first_weight*fakerate_weight;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM0_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM0_down"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM0_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM0_up"]]->Fill(tau_p4.Pt(), final_weight);

		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_down"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_up"]]->Fill(tau_p4.Pt(), final_weight);
		}
		else if (tau_decayMode->at(iTau) == 10) {
		  //for the histos dedicated to DM10 variation, fill with a different weight
		  final_weight = first_weight*fakerate_weight_low;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_down"]]->Fill(tau_p4.Pt(), final_weight);

		  final_weight = first_weight*fakerate_weight_high;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_up"]]->Fill(tau_p4.Pt(), final_weight);

		  //other histos should be filled with the "normal" weight
		  final_weight = first_weight*fakerate_weight;
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM0_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM0_down"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM0_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM0_up"]]->Fill(tau_p4.Pt(), final_weight);

		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_down"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_down"]]->Fill(tau_p4.Pt(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["ev_Mvis_FRS_DM10_up"]]->Fill(vis_p4.M(), final_weight);
		  if (Mt_accept) htau[p_pt][n_eta][m_dm][iPass][iBkgOrSig][iValue][htau_map["tau_pt_FRS_DM10_up"]]->Fill(tau_p4.Pt(), final_weight);
		}

	      }//loop over tau IDs
	    }//loop over MES
	  }//loop over TES
	}//loop over taus
      }//loop over muons
   }//loop over events

   event_file << endl << endl << endl;
   event_file.close();

   TH1F* h_total_events =  new TH1F("weighted_events", "weighted_events", 1, 0, 1);
   h_total_events->Fill(0.5, nEvents);
   file_out->cd();
   h_genmass->Write();
   h_genpt->Write();
   h_total_events->Write();
   h_reweight->Write();
   h_weight_afterallsel->Write();
   h_debug->Write();

   TDirectory* NotauID_dir = file_out->mkdir("NoTauID");
   NotauID_dir->cd();
   for (unsigned int i = 0; i<histo_notauID_names.size(); ++i) for (unsigned int k = 0; k<taun; ++k) hnotauID[k][i]->Write();
   NotauID_dir->Close();

   TDirectory* Gen_dir = file_out->mkdir("GenLevel");
   Gen_dir->cd();
   if (DY) for (unsigned int i = 0; i<hgen.size(); ++i) hgen[i]->Write();
   Gen_dir->Close();

   vector<TDirectory*> tauID_dirs;
   for (unsigned int j = 0; j<nTauIDs; ++j) {
     tauID_dirs.push_back( file_out->mkdir( HPS_WP[j] ) );
     tauID_dirs[j]->cd();
     for (unsigned int i = 0; i<htau_names.size(); ++i) for (unsigned int k = 0; k<taun; ++k) for (unsigned int l = 0; l<passfail.size(); ++l) for (unsigned int m = 0; m<dms.size(); ++m) for (unsigned int n = 0; n<eta.size(); ++n) for (unsigned int p = 0; p<ptrange.size(); ++p) htau[p][n][m][l][k][j][i]->Write();
     tauID_dirs[j]->Close();
   }
   file_out->Close();
}
