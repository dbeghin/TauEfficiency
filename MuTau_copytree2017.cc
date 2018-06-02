#define IIHEAnalysis_cxx
#include "IIHEAnalysis.h"
//#include <TH1.h>
#include <TLorentzVector.h>
//#include <TCanvas.h>
#include "TString.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
  string out = *(argv + 1);
  string out_name= out;
  string in = *(argv + 2);
  string inname= in;
  string mc_in = *(argv + 3);
  string mc_nickname = mc_in;
  string phase_in = *(argv + 4);
  string phase= phase_in;
  string type_in = *(argv + 5);
  string type= type_in;
  TFile *fIn = TFile::Open(inname.c_str());
  TH1F* h1;//dummy TH1F to give as an argument to the loop function
  TH1F* h2;//dummy TH1F to give as an argument to the loop function
  TTree* tree = (TTree*) fIn->Get("IIHEAnalysis");


  IIHEAnalysis* a = new IIHEAnalysis(tree);
  a->Loop(phase, type, out_name, mc_nickname, h1, h2);
  return 0;
}


void IIHEAnalysis::Loop(string phase, string type_of_data, string out_name, string whocares, TH1F* h1, TH1F* h2) {
   if (fChain == 0) return;

   bool data;
   if (type_of_data == "Data" || type_of_data == "data") {
     data = true;
   }
   else {
     data = false;
   }

   //string out_name = "out_"+type_of_data+".root";
   TFile* file_out = new TFile(out_name.c_str(),"RECREATE");
   TH1F* hCounter = new TH1F("h1","h1",1,0,2); //Nevents
   TH1F* hCounter2 = new TH1F("h2","h2",1,0,2); //SumofWeights
   TTree* newTree = fChain->CloneTree(0);

   const float mu_mass = 0.10565837;


   //float event_weight = 1.0;
   //TBranch *newBranch = newTree->Branch("event_weight", &event_weight, "event_weight/F");


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


      //hCounter->Fill(1, event_weight);
      if (!trig_Flag_goodVertices_accept) continue;
      if (!trig_Flag_globalSuperTightHalo2016Filter_accept) continue;
      if (!trig_Flag_HBHENoiseFilter_accept) continue;
      if (!trig_Flag_HBHENoiseIsoFilter_accept) continue;
      if (!trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept) continue;
      if (!trig_Flag_BadPFMuonFilter_accept) continue;
      if (!trig_Flag_BadChargedCandidateFilter_accept) continue;
      if (data) if (!trig_Flag_eeBadScFilter_accept) continue; 
      //hCounter2->Fill(1, event_weight);


      
      //start loop over reconstructed muons
      bool tightmuon = false; //for the tag-and-probe and the High Mass LFV analysis, we need a tight muon
      TLorentzVector tightmu_p4;
      bool tauormu = false; //for the tag-and-probe, we need mumu and mutau events
      TLorentzVector tauormu_p4;
      for (unsigned int iMu = 0; iMu < mu_gt_pt->size(); ++iMu) {
	//is there a tight muon?
	if (mu_gt_pt->at(iMu) < 27.0) continue;
	if (fabs(mu_gt_eta->at(iMu)) > 2.4) continue;
	if (!mu_isPFMuon->at(iMu)) continue;
	if (!mu_isMediumMuon->at(iMu)) continue; //medium ID
	if (fabs(mu_gt_dxy_firstPVtx->at(iMu)) > 0.045) continue;
	if (fabs(mu_gt_dz_firstPVtx->at(iMu)) > 0.2) continue;
	
	tightmu_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu), mu_gt_eta->at(iMu), mu_gt_phi->at(iMu), mu_mass);
	tightmuon = true;

	//is there a 2nd, acceptable muon?
	//careful not to consider the same muon as in the first loop
	for (unsigned int iMu2 = 0; iMu2 < iMu; ++iMu2) {
	  if (mu_gt_pt->at(iMu2) < 20.0) continue;
	  if (fabs(mu_gt_eta->at(iMu2)) > 2.4) continue;
	  if (!mu_isPFMuon->at(iMu2)) continue;
	  if (!mu_isMediumMuon->at(iMu2)) continue; //medium ID
	  if (fabs(mu_gt_dxy_firstPVtx->at(iMu2)) > 0.045) continue;
	  if (fabs(mu_gt_dz_firstPVtx->at(iMu2)) > 0.2) continue;
	  tauormu = true;
	  if (tauormu) break;
	}//second muon loop
	if (tightmuon && tauormu) break;
	
	//start loop over reconstructed taus
	for (unsigned int iTau = 0; iTau < tau_pt->size(); ++iTau) {
	  if (tau_pt->at(iTau) < 17.0) continue;
	  if (fabs(tau_eta->at(iTau)) > 2.3) continue;
	  if (tau_againstMuonTight3->at(iTau) < 0.5) continue;
	  if (tau_againstElectronVLooseMVA6->at(iTau) < 0.5) continue;

	  tauormu = true;
	  if (tauormu) break;
	}//loop over taus
	if (tightmuon && tauormu) break;
	
      }//loop over muons
      if ( !(tightmuon && tauormu) ) continue;
      newTree -> Fill();

   }//loop over events

   file_out->cd();
   hCounter->Write();
   hCounter2->Write();
   newTree->AutoSave();
   file_out->Close();

}
