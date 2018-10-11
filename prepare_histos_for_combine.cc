#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TDirectory.h>
#include "TString.h"
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  float firstbin = 40;
  float lastbin = 150;
  int Nbins = 11;
  int rebin = 10;//144-40=13*8

  TFile* file_in_mutau = new TFile("Figures/allhistos.root", "R");
  TFile* file_in_mumu = new TFile("Figures/histos_mumu.root", "R");


  vector<TString> tauID;
  tauID.push_back("cutbased_loose");
  tauID.push_back("cutbased_medium");
  tauID.push_back("cutbased_tight");

  tauID.push_back("MVA_2016vloose");
  tauID.push_back("MVA_2016loose");
  tauID.push_back("MVA_2016medium");
  tauID.push_back("MVA_2016tight");
  tauID.push_back("MVA_2016vtight");
  tauID.push_back("MVA_2016vvtight");

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

  vector<TString> in_names,                       out_names;
  in_names.push_back("faketau_ev_Mvis");	  out_names.push_back("faketau");
  in_names.push_back("DYS_ev_Mvis");              out_names.push_back("DYS");  
  in_names.push_back("DYB_ev_Mvis");  		  out_names.push_back("DYB");  
  //in_names.push_back("QCD_ev_Mvis");  		  out_names.push_back("QCD");  
  in_names.push_back("VV_ev_Mvis");   		  out_names.push_back("VV");   
  //in_names.push_back("TTS_ev_Mvis");   		  out_names.push_back("TTS");   
  in_names.push_back("TTB_ev_Mvis");   		  out_names.push_back("TTB");   
  in_names.push_back("data_ev_Mvis"); 		  out_names.push_back("data_obs");

  vector<TString> in_TES,      out_TES;
  in_TES.push_back("");        out_TES.push_back("");
  in_TES.push_back("TESup_");  out_TES.push_back("_tesUp");
  in_TES.push_back("TESdown_");out_TES.push_back("_tesDown");


  for (unsigned int i=0; i<tauID.size(); ++i) {
    TFile* file_out = new TFile("prefit_histos_"+tauID[i]+".root", "RECREATE");
    file_out->cd();


    TDirectory* pass_dir = file_out->mkdir("pass");
    pass_dir->cd();
    for (unsigned int j=0; j<in_names.size(); ++j) {
      for (unsigned int k=0; k<in_TES.size(); ++k) {
	TH1F* h = (TH1F*) file_in_mutau->Get(in_names[j]+"_"+in_TES[k]+tauID[i]+"_pass");
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
	h2->SetName(out_names[j]+out_TES[k]);
	if (out_names[j] == "data_obs" && k >0) continue;
	h2->Write();
      }
    }
    pass_dir->Close();



    /*TDirectory* fail_dir = file_out->mkdir("fail");
    fail_dir->cd();
    for (unsigned int j=0; j<in_names.size(); ++j) {
      for (unsigned int k=0; k<in_TES.size(); ++k) {
	TH1F* h = (TH1F*) file_in_mutau->Get(in_names[j]+"_"+in_TES[k]+tauID[i]+"_fail");
	h->Rebin(rebin);
	TH1F* h2 = new TH1F(h->GetName(), h->GetTitle(), Nbins, firstbin, lastbin);
	int ll=0;
	float last_bin_error = 10;
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
	  if (bin_error == 0) {
	    bin_error = last_bin_error;
	  }
	  else {
	    last_bin_error = bin_error;
	  }
	  //if (bin_error = 10) cout << bin_error << endl;
	  h2->SetBinError(ll, bin_error);
	}
	h2->SetName(out_names[j]+out_TES[k]);
	if (out_names[j] == "data_obs" && k >0) continue;
	h2->Write();
      }
    }
    fail_dir->Close();*/


    TDirectory* Zmumu_dir = file_out->mkdir("zmumu");
    Zmumu_dir->cd();
    for (unsigned int j=0; j<in_names.size(); ++j) {
      TString inn = in_names[j];
      if (j==0) inn = "WJets_ev_Mvis";
      TH1F* h = (TH1F*) file_in_mumu->Get(inn);
      TH1F* h2 = new TH1F(h->GetName(), h->GetTitle(), 1, 60, 120);
      int ll=0;
      for (unsigned int l=1; l<h->GetNbinsX()+1; ++l) {
	//if (h->GetXaxis()->GetBinCenter(l) < 35) continue;
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
      h2->SetName(out_names[j]);
      h2->Write();
    }
    Zmumu_dir->Close();


    file_out->Close();
  }



  return 0;
}
