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
int main(int argc, char** argv) {
  //string nature = *(argv + 1);

  TFile* file_out;
  TFile* file_in;

  file_out = new TFile("TauID2017/Faketaus.root", "RECREATE");
  file_in = new TFile("Figures/TauFakes.root", "R");
  

  vector<TString> names;
  names.push_back("data_");//0
  //names.push_back("WJets_");//1
  names.push_back("DYB_");
  names.push_back("DYS_");
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
  vars.push_back("ev_Nvertex_MinBiasdown");
  vars.push_back("ev_Nvertex_Minbiasup");
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

  HPS_WP.push_back("MVA_2016vloose");
  HPS_WP.push_back("MVA_2016loose");
  HPS_WP.push_back("MVA_2016medium");
  HPS_WP.push_back("MVA_2016tight");
  HPS_WP.push_back("MVA_2016vtight");
  HPS_WP.push_back("MVA_2016vvtight");

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

  //retrieve histograms from all control regions
  //only for CR4 (1) do we care to have all histos
  vector<TH1F*> h[names.size()][vars.size()];
  for (unsigned int j=0; j<names.size(); ++j) {
    for (unsigned int k=0; k<vars.size(); ++k) {
      for (unsigned int l=0; l<HPS_WP.size(); ++l) {
	h[j][k].push_back( (TH1F*) file_in->Get(names[j]+vars[k]+"_"+HPS_WP[l]+"_fail") );
	h[j][k][l]->SetName(names[j]+vars[k]);
      }
    }
  }


  file_out->cd();
  for (unsigned int l=0; l<HPS_WP.size(); ++l) {
    for (unsigned int k=0; k<vars.size(); ++k) {
      TH1F* h_faketau = (TH1F*) h[0][k][l]->Clone("faketau_"+vars[k]+"_"+HPS_WP[l]+"_pass");
      for (unsigned int j=1; j<names.size(); ++j) h_faketau->Add(h[j][k][l], -1);//subtract all real tau bg

      for (unsigned int iBin = 0; iBin<h_faketau->GetNbinsX(); ++iBin) {
	if (h_faketau->GetBinContent(iBin) < 0) h_faketau->SetBinContent(iBin,0);
      }
      h_faketau->Write();
    }
  }
  file_out->Close();


  return 0;
}
