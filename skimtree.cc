#define IIHEAnalysis_cxx
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
#include "IIHEAnalysis.h"

using namespace std;

int main(int argc, char** argv) {
  string out = *(argv + 1);
  string in = *(argv + 2);
  TString in_name = in;

  TFile *file_in= new TFile(in_name, "R");
  TTree* tree = (TTree*) file_in->Get("IIHEAnalysis");

  IIHEAnalysis* a = new IIHEAnalysis(tree);
  string phase = "eh", type = "eh", mc_nickname = "f";
  a->Loop(phase, type, out, mc_nickname);
  return 0;
}


void IIHEAnalysis::Loop(string phase, string type_of_data, string out_name, string whocares) {
  if (fChain == 0) return;

  TString proper_out_name = out_name;
  TFile *newfile = new TFile(proper_out_name,"RECREATE");
  TTree* skimmedTree = fChain->CloneTree(0);


  const float mu_mass = 0.10565837;

  Long64_t nEntries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int print_count = 0;
  //start loop over all events
  for (Long64_t jEntry = 0; jEntry < nEntries; ++jEntry) {
    Long64_t iEntry = LoadTree(jEntry);
    if (iEntry < 0) break;
    if (jEntry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", jEntry, nEntries);
      
    nb = fChain->GetEntry(jEntry);
    nbytes += nb;


    //Is one of the triggers fired?
    bool PassMuonTrigger = false;
    if (trig_HLT_IsoMu27_accept) PassMuonTrigger = true;
    if (!PassMuonTrigger) continue;

    bool keep_event = false;
    for (unsigned int iMu = 0; iMu < mu_gt_pt->size(); ++iMu) {
      if (keep_event) break;
      if (mu_gt_pt->at(iMu) < 27.0) continue;
      if (fabs(mu_gt_eta->at(iMu)) > 2.4) continue;
      if (!mu_isPFMuon->at(iMu)) continue;
      if (!mu_isMediumMuon->at(iMu)) continue;
      if (mu_pfIsoDbCorrected04->at(iMu) > 0.15) continue;
      if (fabs(mu_gt_dxy_firstPVtx->at(iMu)) > 0.045) continue;
      if (fabs(mu_gt_dz_firstPVtx->at(iMu)) > 0.2) continue;

      for (unsigned int iTau = 0; iTau < tau_pt->size(); ++iTau) {
	if (fabs(tau_eta->at(iTau)) > 2.3) continue;
	if (tau_againstMuonTight3->at(iTau) < 0.5) continue;
	if (tau_againstElectronVLooseMVA6->at(iTau) < 0.5) continue;
	keep_event = true;
	break;
      }
      for (unsigned int iMu2 = 0; iMu2 < mu_gt_pt->size(); ++iMu2) {
	if (iMu == iMu2) continue;
	if (mu_gt_pt->at(iMu2) < 20.0) continue;
	if (fabs(mu_gt_eta->at(iMu2)) > 2.4) continue;
	if (!mu_isPFMuon->at(iMu2)) continue;
	if (!mu_isMediumMuon->at(iMu2)) continue;
	if (fabs(mu_gt_dxy_firstPVtx->at(iMu2)) > 0.045) continue;
	if (fabs(mu_gt_dz_firstPVtx->at(iMu2)) > 0.2) continue;
	keep_event = true;
	break;
      }
    }

    if (keep_event) skimmedTree->Fill();
  }//loop over events

  newfile->cd();
  skimmedTree->AutoSave();
  newfile->Write();
  newfile->Close();
}
