#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include <vector>
#include <iostream>
#include "TString.h"

using namespace std;

float norm_F(float x, float y){

  return sqrt(x*x+y*y);

}


double FakeRate(double taupt, TString HPS_WP, TString DM, TString eta) {
  if (taupt > 150) taupt = 149;

  double SF=0.5;

  TFile* fake_file = new TFile("Reweighting/fakerate.root","R");
  TString fake_string = "FakeRate_byTauPt_data_"+HPS_WP+"_"+DM+"_"+eta;
  TH1F* fake_histo = (TH1F*) fake_file->Get(fake_string);

  int iBin = -1;
  if (taupt >= 20) iBin = fake_histo->FindBin(taupt);
  if (iBin != -1) SF = fake_histo->GetBinContent(iBin);
  fake_file->Close();

  double reweight = 1;
  if (SF != 1) reweight = SF/(1-SF);

  return reweight;

}

double FakeRateFlat(TString HPS_WP, TString DM) {
  double SF=0.5;

  TFile* fake_file = new TFile("Reweighting/fakerate.root","R");
  TString fake_string = "ptratio_data_"+HPS_WP+"_"+DM+"_total";
  TH1F* fake_histo = (TH1F*) fake_file->Get(fake_string);

  int iBin = 1; //the first bin, with all low pt taus, is the interesting one                                                                                                                                                                 
  SF = fake_histo->GetBinContent(iBin);
  fake_file->Close();

  double reweight = 1;
  if (SF != 1) reweight = SF/(1-SF);

  return reweight;

}


pair<double,double> getSF (float mupt, float mueta) {
  TFile* id_file = new TFile("Reweighting/RunBCDEF_SF_ID.root","R");
  TH2F* id_histo = (TH2F*) id_file->Get("NUM_MediumID_DEN_genTracks_pt_abseta");
  int bin_in = id_histo->FindBin(mupt, fabs(mueta));
  double id_sf = id_histo->GetBinContent(bin_in);
  id_file->Close();

  TFile* iso_file = new TFile("Reweighting/RunBCDEF_SF_ISO.root","R");
  TH2F* iso_histo = (TH2F*) iso_file->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta");
  bin_in = iso_histo->FindBin(mupt, fabs(mueta));
  double iso_sf = iso_histo->GetBinContent(bin_in);
  iso_file->Close();

  if (id_sf == 0) id_sf = 1.0;
  if (iso_sf == 0) iso_sf = 1.0;

  pair<double, double> SF_pair;
  SF_pair.first =  id_sf;  SF_pair.second =iso_sf;
  return  SF_pair;
}


double GetTriggerMuonIDMuonIsoReweight(float mu_pt, float mu_eta) {
  //scale factor files that need to be open                                                                                                                                                                                                   
  TFile* tr_file = new TFile("Reweighting/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","R");
  TH2F* tr_histo = (TH2F*) tr_file->Get("IsoMu27_PtEtaBins/pt_abseta_ratio");
  int bin_in = tr_histo->FindBin(mu_pt, fabs(mu_eta));
  double tr_sf = tr_histo->GetBinContent(bin_in);
  tr_file->Close();
  if (tr_sf == 0) tr_sf = 1.0;

  float muID_sf = getSF(mu_pt, mu_eta).first, muIso_sf = getSF(mu_pt, mu_eta).second;

  double weight = tr_sf*muID_sf*muIso_sf;

  return weight;
}


double GetReweight_mumu(float mu1_pt, float mu1_eta, float mu2_pt, float mu2_eta) {
  //scale factor files that need to be open                                                                                                                                                                                                   
  TFile* tr_file = new TFile("Reweighting/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","R");
  TH2F* tr_data = (TH2F*) tr_file->Get("IsoMu27_PtEtaBins/efficienciesDATA/pt_abseta_DATA");
  TH2F* tr_MC = (TH2F*) tr_file->Get("IsoMu27_PtEtaBins/efficienciesMC/pt_abseta_MC");
  //Eff. of the trigger: COMPLEMENTARY of the prob. that NONE of the two muons will trigger                                                                                                                                                   

  int bin_in = tr_data->FindBin(mu1_pt, fabs(mu1_eta));
  double mu1_eff_data = tr_data->GetBinContent(bin_in);
  bin_in = tr_data->FindBin(mu2_pt, fabs(mu2_eta));
  double mu2_eff_data = tr_data->GetBinContent(bin_in);
  double tr_eff_data = 1 - (1 - mu1_eff_data)*(1 - mu2_eff_data);
  if (tr_eff_data == 0) tr_eff_data = 1.0;

  bin_in = tr_MC->FindBin(mu1_pt, fabs(mu1_eta));
  double mu1_eff_MC = tr_MC->GetBinContent(bin_in);
  bin_in = tr_data->FindBin(mu2_pt, fabs(mu2_eta));
  double mu2_eff_MC = tr_MC->GetBinContent(bin_in);
  double tr_eff_MC = 1 - (1 - mu1_eff_MC)*(1 - mu2_eff_MC);
  tr_file->Close();
  if (tr_eff_MC == 0) tr_eff_MC = 1.0;

  double tr_sf = tr_eff_data/tr_eff_MC;

  float mu1ID_sf = getSF(mu1_pt, mu1_eta).first, mu1Iso_sf = getSF(mu1_pt, mu1_eta).second;
  float mu2ID_sf = getSF(mu2_pt, mu2_eta).first, mu2Iso_sf = getSF(mu2_pt, mu2_eta).second;


  double weight=tr_sf*mu1ID_sf*mu1Iso_sf*mu2ID_sf*mu2Iso_sf;

  return weight;
}



