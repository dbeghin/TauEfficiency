#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TDirectory.h>
#include "TString.h"
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  float firstbin = 50;
  float lastbin = 150;
  int Nbins = 10;
  int rebin = 10;//144-40=13*8

  TFile* file_in_mutau = new TFile("Figures/allhistos.root", "R");
  TFile* file_in_mumu = new TFile("Figures/histos_mumu.root", "R");


  vector<TString> tauID;
  tauID.push_back("cutbased_loose");
  tauID.push_back("cutbased_medium");
  tauID.push_back("cutbased_tight");

  tauID.push_back("MVA_2017v1vvloose");
  tauID.push_back("MVA_2017v1vloose");
  tauID.push_back("MVA_2017v1loose");
  tauID.push_back("MVA_2017v1medium");
  tauID.push_back("MVA_2017v1tight");
  tauID.push_back("MVA_2017v1vtight");
  tauID.push_back("MVA_2017v1vvtight");

  tauID.push_back("MVA_2017v2vvloose");
  tauID.push_back("MVA_2017v2vloose");
  tauID.push_back("MVA_2017v2loose");
  tauID.push_back("MVA_2017v2medium");
  tauID.push_back("MVA_2017v2tight");
  tauID.push_back("MVA_2017v2vtight");
  tauID.push_back("MVA_2017v2vvtight");

  tauID.push_back("MVA_DBdR03vloose");
  tauID.push_back("MVA_DBdR03loose");                                                  
  tauID.push_back("MVA_DBdR03medium");                                                  
  tauID.push_back("MVA_DBdR03tight");                                                  
  tauID.push_back("MVA_DBdR03vtight");                                                  
  tauID.push_back("MVA_DBdR03vvtight");                                                 

  tauID.push_back("MVA_PWdR03vloose");
  tauID.push_back("MVA_PWdR03loose");
  tauID.push_back("MVA_PWdR03medium");
  tauID.push_back("MVA_PWdR03tight");
  tauID.push_back("MVA_PWdR03vtight");
  tauID.push_back("MVA_PWdR03vvtight");


  vector<TString> dms;
  dms.push_back("allDMs");
  dms.push_back("DM0");
  dms.push_back("DM1");
  dms.push_back("DM10");

  vector<TString> eta;
  eta.push_back("alleta");
  eta.push_back("barrel");
  eta.push_back("endcap");

  vector<TString> ptrange;
  ptrange.push_back("allpt");
  ptrange.push_back("pt_20_40");
  ptrange.push_back("pt_40_150");


  vector<TString> in_names,                       out_names;
  in_names.push_back("faketau_ev_Mvis");	  out_names.push_back("faketau");
  in_names.push_back("DYS_ev_Mvis");              out_names.push_back("DYS");  
  in_names.push_back("DYB_ev_Mvis");  		  out_names.push_back("DYB");  
  in_names.push_back("VV_ev_Mvis");   		  out_names.push_back("VV");   
  in_names.push_back("TTB_ev_Mvis");   		  out_names.push_back("TTB");   
  in_names.push_back("data_ev_Mvis"); 		  out_names.push_back("data_obs");

  vector<TString> in_sys,      out_sys;
  in_sys.push_back("");        out_sys.push_back("");
  in_sys.push_back("TESup_");  out_sys.push_back("_tesUp");
  in_sys.push_back("TESdown_");out_sys.push_back("_tesDown");
  in_sys.push_back("MESup_");  out_sys.push_back("_mesUp");
  in_sys.push_back("MESdown_");out_sys.push_back("_mesDown");
  in_sys.push_back("MinBiasup_");  out_sys.push_back("_minbiasUp");
  in_sys.push_back("MinBiasdown_");out_sys.push_back("_minbiasDown");
  in_sys.push_back("FRS_DM0_up_");  out_sys.push_back("_frsdm0Up");
  in_sys.push_back("FRS_DM0_down_");out_sys.push_back("_frsdm0Down");
  in_sys.push_back("FRS_DM1_up_");  out_sys.push_back("_frsdm1Up");
  in_sys.push_back("FRS_DM1_down_");out_sys.push_back("_frsdm1Down");
  in_sys.push_back("FRS_DM10_up_");  out_sys.push_back("_frsdm10Up");
  in_sys.push_back("FRS_DM10_down_");out_sys.push_back("_frsdm10Down");
  in_sys.push_back("antiisomu_up_");  out_sys.push_back("_antiisomuUp");
  in_sys.push_back("antiisomu_down_");out_sys.push_back("_antiisomuDown");
  in_sys.push_back("antiisotau_up_");  out_sys.push_back("_antiisotauUp");
  in_sys.push_back("antiisotau_down_");out_sys.push_back("_antiisotauDown");


  for (unsigned int i=0; i<tauID.size(); ++i) {
    for (unsigned int l=0; l<dms.size(); ++l) {
      for (unsigned int m=0; m<eta.size(); ++m) {
	for (unsigned int p=0; p<ptrange.size(); ++p) {
	  TString fileID = dms[l]+"_"+eta[m]+"_"+ptrange[p]+"_"+tauID[i];
	  TH1F* h0 = (TH1F*) file_in_mutau->Get(tauID[i]+"/data_ev_Mvis_"+fileID+"_pass");
	  cout << tauID[i]+"/data_ev_Mvis_"+fileID+"_pass" << endl;
	  if (h0 == 0) continue;
	  TFile* file_out = new TFile("HistosForCombine/prefit_histos_"+fileID+".root", "RECREATE");
	  file_out->cd();


	  TDirectory* pass_dir = file_out->mkdir("pass");
	  pass_dir->cd();
	  for (unsigned int j=0; j<in_names.size(); ++j) {
	    for (unsigned int k=0; k<in_sys.size(); ++k) {
	      cout << in_names[j]+"_"+in_sys[k]+fileID+"_pass" << endl;
	      TH1F* h = (TH1F*) file_in_mutau->Get(tauID[i]+"/"+in_names[j]+"_"+in_sys[k]+fileID+"_pass");
	      h->Rebin(rebin);
	      TH1F* h2 = new TH1F(h->GetName(), h->GetTitle(), Nbins, firstbin, lastbin);
	      int ll=0;
	      for (unsigned int l=1; l<h->GetNbinsX()+1; ++l) {
		if (h->GetXaxis()->GetBinCenter(l) < firstbin) continue;
		ll += 1;
		float bin_content = h->GetBinContent(l);
		if (bin_content <= 0) {
		  h2->SetBinContent(ll, 0);
		}
		else {
		  h2->SetBinContent(ll, bin_content);
		}
		float bin_error = h->GetBinError(l);
		h2->SetBinError(ll, bin_error);
	      }
	      h2->SetName(out_names[j]+out_sys[k]);
	      if (out_names[j] == "data_obs" && k >0) continue;
	      h2->Write();
	    }
	  }
	  pass_dir->Close();



	  TDirectory* Zmumu_dir = file_out->mkdir("zmumu");
	  Zmumu_dir->cd();
	  for (unsigned int j=0; j<in_names.size(); ++j) {
	    TString inn = in_names[j];
	    if (j==0) inn = in_names[1];
	    TH1F* h = (TH1F*) file_in_mumu->Get(inn);
	    TH1F* h2 = new TH1F(h->GetName(), h->GetTitle(), 1, 60, 120);
	    int ll=0;
	    for (unsigned int l=1; l<h->GetNbinsX()+1; ++l) {
	      //if (h->GetXaxis()->GetBinCenter(l) < 35) continue;
	      ll += 1;
	      float bin_content = h->GetBinContent(l);
	      if (j==0) bin_content = 0;
	      if (bin_content <= 0) {
		h2->SetBinContent(ll, 0);
	      }
	      else {
		h2->SetBinContent(ll, bin_content);
	      }
	      float bin_error = h->GetBinError(l);
	      h2->SetBinError(ll, bin_error);
	    }
	    h2->SetName(out_names[j]);
	    h2->Write();
	  }
	  Zmumu_dir->Close();


	  file_out->Close();
	}
      }
    }
  }



  return 0;
}
