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
  vars.push_back("ev_Mvis_TESdown");
  vars.push_back("ev_Mvis_TESup");
  vars.push_back("tau_eta");
  vars.push_back("tau_phi");
  vars.push_back("mu_pt");
  vars.push_back("mu_eta");
  vars.push_back("mu_phi");
  vars.push_back("ev_Nvertex");

  TString tauID = "_MVA_tight";

  //retrieve histograms from all control regions
  //only for CR4 (1) do we care to have all histos
  vector<TH1F*> h[names.size()];
  for (unsigned int j=0; j<names.size(); ++j) {
    for (unsigned int k=0; k<vars.size(); ++k) {
      h[j].push_back( (TH1F*) file_in->Get(names[j]+vars[k]+tauID+"_fail") );
      h[j][k]->SetName(names[j]+vars[k]);
    }
  }


  file_out->cd();
  for (unsigned int k=0; k<vars.size(); ++k) {
    TH1F* h_faketau = (TH1F*) h[0][k]->Clone("faketau_"+vars[k]+tauID+"_pass");
    for (unsigned int j=1; j<names.size(); ++j) h_faketau->Add(h[j][k], -1);//subtract all real tau bg

    for (unsigned int iBin = 0; iBin<h_faketau->GetNbinsX(); ++iBin) {
      if (h_faketau->GetBinContent(iBin) < 0) h_faketau->SetBinContent(iBin,0);
    }
    h_faketau->Write();
  }
  file_out->Close();


  return 0;
}
