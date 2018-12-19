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

using namespace std;
int main(int argc, char** argv) {
  string nature = *(argv + 1);

  TFile* file_out;
  TFile* file_in;

  bool final = false;

  if (nature == "final" || nature == "Final") {
    file_out = new TFile("TauID2017/Faketaus.root", "RECREATE");
    file_in = new TFile("Figures/TauFakes.root", "R");
    final = true;
  }
  else if (nature == "WJets" || nature == "Wjets" || nature == "wjets") {
    file_out = new TFile("TauID2017/wjetsfake.root", "RECREATE");
    file_in = new TFile("Figures/WJetsFake.root", "R");
  }
  

  vector<TString> names;
  names.push_back("data_");//0
  //names.push_back("WJets_");//1
  names.push_back("DYB_");
  if (final) names.push_back("DYS_");
  names.push_back("TTB_");
  names.push_back("VV_");

  vector<TString> vars;
  vars.push_back("tau_pt");
  //vars.push_back("tau_MVA");
  vars.push_back("ev_Mvis");
  vars.push_back("ev_Mt");
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

  vector<TString> HPS_WP;
  HPS_WP.push_back("cutbased_loose");
  HPS_WP.push_back("cutbased_medium");
  HPS_WP.push_back("cutbased_tight");

  HPS_WP.push_back("MVA_2017v1vvloose");
  HPS_WP.push_back("MVA_2017v1vloose");
  HPS_WP.push_back("MVA_2017v1loose");
  HPS_WP.push_back("MVA_2017v1medium");
  HPS_WP.push_back("MVA_2017v1tight");
  HPS_WP.push_back("MVA_2017v1vtight");
  HPS_WP.push_back("MVA_2017v1vvtight");

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



  //retrieve histograms from all control regions
  //only for CR4 (1) do we care to have all histos
  vector<TH1F*> h[names.size()][vars.size()][dms.size()][eta.size()][HPS_WP.size()];
  for (unsigned int i=0; i<names.size(); ++i) {
    for (unsigned int j=0; j<vars.size(); ++j) {
      for (unsigned int k=0; k<dms.size(); ++k) {
	for (unsigned int l=0; l<eta.size(); ++l) {
	  for (unsigned int m=0; m<HPS_WP.size(); ++m) {
	    for (unsigned int p=0; p<ptrange.size(); ++p) {
	      h[i][j][k][l][m].push_back( (TH1F*) file_in->Get(HPS_WP[m]+"/"+names[i]+vars[j]+"_"+dms[k]+"_"+eta[l]+"_"+ptrange[p]+"_"+HPS_WP[m]+"_fail") );
	      h[i][j][k][l][m][p]->SetName(names[i]+vars[j]+"_"+dms[k]+"_"+eta[l]+"_"+ptrange[p]+"_"+HPS_WP[m]);
	    }
	  }
	}
      }
    }
  }


  file_out->cd();
  for (unsigned int m=0; m<HPS_WP.size(); ++m) {
    TDirectory* work_dir = file_out->mkdir(HPS_WP[m]);
    work_dir->cd();
    for (unsigned int l=0; l<eta.size(); ++l) {
      for (unsigned int k=0; k<dms.size(); ++k) {
	for (unsigned int j=0; j<vars.size(); ++j) {
	  for (unsigned int p=0; p<ptrange.size(); ++p) {
	    TH1F* h_faketau = (TH1F*) h[0][j][k][l][m][p]->Clone("faketau_"+vars[j]+"_"+dms[k]+"_"+eta[l]+"_"+ptrange[p]+"_"+HPS_WP[m]+"_pass");
	    for (unsigned int i=1; i<names.size(); ++i) h_faketau->Add(h[i][j][k][l][m][p], -1);//subtract all real tau bg

	    for (unsigned int iBin = 0; iBin<h_faketau->GetNbinsX(); ++iBin) {
	      if (h_faketau->GetBinContent(iBin) < 0) h_faketau->SetBinContent(iBin,0);
	    }
	    h_faketau->Write();
	  }
	}
      }
    }
    work_dir->Close();
  }
  file_out->Close();


  return 0;
}
