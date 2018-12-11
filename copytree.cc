#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "vector"
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <map>
#include <sys/stat.h>
#include <TVector3.h>
#include <TROOT.h>

using namespace std;

int main(int argc, char** argv) {
  string out = *(argv + 1);
  string in = *(argv + 2);
  TString out_name= out;
  TString in_name = in;

  TFile *file_in= new TFile(in_name, "R");
  TFile *newfile = new TFile(out_name,"RECREATE");

  TTree* oldtree = (TTree*) file_in->Get("IIHEAnalysis");


  //IIHEAnalysis *ev_event   = new IIHEAnalysis();
  //oldtree->SetBranchAddress("ev_event",&ev_event);

  oldtree->SetBranchStatus("*",0);
  oldtree->SetBranchStatus("ev_event",1);
  oldtree->SetBranchStatus("ev_run",1);
  oldtree->SetBranchStatus("ev_luminosityBlock",1);

  oldtree->SetBranchStatus("trig_HLT_IsoMu27_accept",1);
  oldtree->SetBranchStatus("pv_n",1);
  oldtree->SetBranchStatus("mc_mother_index",1);
  oldtree->SetBranchStatus("mc_pdgId",1);
  oldtree->SetBranchStatus("mc_px",1);
  oldtree->SetBranchStatus("mc_py",1);
  oldtree->SetBranchStatus("mc_pz",1);
  oldtree->SetBranchStatus("mc_energy",1);
  oldtree->SetBranchStatus("mc_eta",1);
  oldtree->SetBranchStatus("mc_phi",1);
  oldtree->SetBranchStatus("mc_pt",1);
  oldtree->SetBranchStatus("LHE_pdgid",1);
  oldtree->SetBranchStatus("LHE_Pt",1);
  oldtree->SetBranchStatus("mu_isPFMuon",1);
  oldtree->SetBranchStatus("mu_gt_pt",1);
  oldtree->SetBranchStatus("mu_gt_eta",1);
  oldtree->SetBranchStatus("mu_gt_phi",1);
  oldtree->SetBranchStatus("mu_gt_px",1);
  oldtree->SetBranchStatus("mu_gt_py",1);
  oldtree->SetBranchStatus("mu_gt_pz",1);
  oldtree->SetBranchStatus("mu_gt_dxy_firstPVtx",1);
  oldtree->SetBranchStatus("mu_gt_dz_firstPVtx",1);
  oldtree->SetBranchStatus("mu_pfIsoDbCorrected04",1);
  oldtree->SetBranchStatus("mu_isMediumMuon",1);
  oldtree->SetBranchStatus("mu_isGlobalMuon",1);
  oldtree->SetBranchStatus("mu_isTrackerMuon",1);
  oldtree->SetBranchStatus("mu_isPFMuon",1);
  oldtree->SetBranchStatus("mu_gt_charge",1);
  oldtree->SetBranchStatus("gsf_VIDMVAMedium",1);
  oldtree->SetBranchStatus("gsf_pt",1);
  oldtree->SetBranchStatus("gsf_eta",1);
  oldtree->SetBranchStatus("gsf_dxy_firstPVtx",1);
  oldtree->SetBranchStatus("gsf_dz_firstPVtx",1);
  oldtree->SetBranchStatus("gsf_passConversionVeto",1);
  oldtree->SetBranchStatus("gsf_nLostInnerHits",1);
  oldtree->SetBranchStatus("gsf_relIso",1);
  oldtree->SetBranchStatus("jet_CSVv2",1);
  oldtree->SetBranchStatus("jet_pt",1);
  oldtree->SetBranchStatus("jet_eta",1);
  oldtree->SetBranchStatus("tau_pt",1);
  oldtree->SetBranchStatus("tau_eta",1);
  oldtree->SetBranchStatus("tau_phi",1);
  oldtree->SetBranchStatus("tau_px",1);
  oldtree->SetBranchStatus("tau_py",1);
  oldtree->SetBranchStatus("tau_pz",1);
  oldtree->SetBranchStatus("tau_energy",1);
  oldtree->SetBranchStatus("tau_decayModeFinding",1);
  oldtree->SetBranchStatus("tau_againstMuonTight3",1);
  oldtree->SetBranchStatus("tau_againstElectronVLooseMVA6",1);
  oldtree->SetBranchStatus("tau_ptLeadChargedCand",1);
  oldtree->SetBranchStatus("tau_lead_dz",1);
  oldtree->SetBranchStatus("tau_charge",1);
  oldtree->SetBranchStatus("tau_decayMode",1);
  oldtree->SetBranchStatus("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits",1);
  oldtree->SetBranchStatus("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",1);
  oldtree->SetBranchStatus("tau_byTightCombinedIsolationDeltaBetaCorr3Hits",1);
  oldtree->SetBranchStatus("tau_byIsolationMVArun2v1DBoldDMwLTraw",1);
  oldtree->SetBranchStatus("tau_byIsolationMVArun2v1DBnewDMwLTraw",1);
  oldtree->SetBranchStatus("tau_byVLooseIsolationMVArun2v1DBoldDMwLT",1);
  oldtree->SetBranchStatus("tau_byLooseIsolationMVArun2v1DBoldDMwLT",1);
  oldtree->SetBranchStatus("tau_byMediumIsolationMVArun2v1DBoldDMwLT",1);
  oldtree->SetBranchStatus("tau_byTightIsolationMVArun2v1DBoldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVTightIsolationMVArun2v1DBoldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVVTightIsolationMVArun2v1DBoldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byTightIsolationMVArun2017v1DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byTightIsolationMVArun2017v2DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017",1);
  oldtree->SetBranchStatus("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byTightIsolationMVArun2v1PWdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT",1);
  oldtree->SetBranchStatus("tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT",1);

  oldtree->SetBranchStatus("MET_eefix_Px",1);
  oldtree->SetBranchStatus("MET_eefix_Py",1);
  oldtree->SetBranchStatus("MET_eefix_Pt",1);
  oldtree->SetBranchStatus("MET_nominal_Px",1);
  oldtree->SetBranchStatus("MET_nominal_Py",1);
  oldtree->SetBranchStatus("MET_nominal_Pt",1);

  oldtree->SetBranchStatus("mc_weight",1);
  oldtree->SetBranchStatus("mc_w_sign",1);
  oldtree->SetBranchStatus("mc_trueNumInteractions",1);
  oldtree->SetBranchStatus("gsf_energy",1);
  oldtree->SetBranchStatus("gsf_et",1);
  oldtree->SetBranchStatus("jet_px",1);
  oldtree->SetBranchStatus("jet_py",1);
  oldtree->SetBranchStatus("jet_pz",1);
  oldtree->SetBranchStatus("jet_energy",1);
  oldtree->SetBranchStatus("tau_decayModeFindingNewDMs",1);
  oldtree->SetBranchStatus("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits",1);
  oldtree->SetBranchStatus("tau_isPFTau",1);
  oldtree->SetBranchStatus("tau_byIsolationMVArun2017v1DBoldDMwLTraw2017",1);
  oldtree->SetBranchStatus("tau_byIsolationMVArun2017v2DBoldDMwLTraw2017",1);
  oldtree->SetBranchStatus("tau_lead_dxy",1);
  
  TTree *newtree = oldtree->CloneTree();
  file_in->Close();
  
  newfile->cd();
  newtree->AutoSave();
  newfile->Write();
  newfile->Close();

}
