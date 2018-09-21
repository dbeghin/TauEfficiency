#define IIHEAnalysis_cxx
#include "IIHEAnalysis.h"
//#include <TH1.h>
#include <TLorentzVector.h>
//#include <TCanvas.h>
#include "TString.h"
#include <iostream>
#include "PU_reWeighting.cc"

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
  TH1F* hCounter = (TH1F*) fIn->Get("h1");
  TH1F* hCounter2 = (TH1F*) fIn->Get("h2");
  TTree* tree = (TTree*) fIn->Get("IIHEAnalysis");

  IIHEAnalysis* a = new IIHEAnalysis(tree);
  a->Loop(phase, type, out_name, mc_nickname, hCounter, hCounter2);
  return 0;
}

void IIHEAnalysis::Loop(string phase, string type_of_data, string out_name, string mc_nick, TH1F* hCounter, TH1F* hCounter2) {
   if (fChain == 0) return;

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

   //string out_name = "out_"+type_of_data+".root";
   TFile* file_out = new TFile(out_name.c_str(),"RECREATE");

   const float mu_mass = 0.10565837;


   //list here the names and x-axis ranges of all reco-level histos we wish to create :
   vector<TString> histo_names;             vector<int> nBins;     vector<float> x_min,   x_max; 
   histo_names.push_back("mu_pt");          nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(100);
   histo_names.push_back("mu_eta");         nBins.push_back(50);   x_min.push_back(-2.5); x_max.push_back(2.5);
   histo_names.push_back("mu_phi");         nBins.push_back(64);   x_min.push_back(-3.2); x_max.push_back(3.2);
   histo_names.push_back("mu1_pt");         nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(100);
   histo_names.push_back("mu1_eta");        nBins.push_back(50);   x_min.push_back(-2.5); x_max.push_back(2.5);
   histo_names.push_back("mu1_phi");        nBins.push_back(64);   x_min.push_back(-3.2); x_max.push_back(3.2);
   histo_names.push_back("mu2_pt");         nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(100);
   histo_names.push_back("mu2_eta");        nBins.push_back(50);   x_min.push_back(-2.5); x_max.push_back(2.5);
   histo_names.push_back("mu2_phi");        nBins.push_back(64);   x_min.push_back(-3.2); x_max.push_back(3.2);
   histo_names.push_back("ev_DRmumu");      nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(10);
   histo_names.push_back("ev_Mt_raw");      nBins.push_back(150);  x_min.push_back(0);    x_max.push_back(150);
   histo_names.push_back("ev_Mt");          nBins.push_back(150);  x_min.push_back(0);    x_max.push_back(150);
   histo_names.push_back("ev_Mvis");        nBins.push_back(150);  x_min.push_back(0);    x_max.push_back(150);
   histo_names.push_back("ev_METmumass");   nBins.push_back(150);  x_min.push_back(0);    x_max.push_back(150);
   histo_names.push_back("ev_MET");         nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(100);
   histo_names.push_back("ev_METphi");      nBins.push_back(64);   x_min.push_back(-3.2); x_max.push_back(3.2);
   histo_names.push_back("ev_Nvertex");     nBins.push_back(81);   x_min.push_back(-0.5); x_max.push_back(80.5);

   vector<TH1F*> h;
   for (unsigned int i = 0; i<histo_names.size(); ++i) {
     h.push_back( new TH1F(histo_names[i], histo_names[i], nBins[i], x_min[i], x_max[i]) ); 
   }

   vector<TString> htau_names;              vector<int> nBins_tau;     vector<float> x_min_tau,   x_max_tau; 
   htau_names.push_back("taupt_pass");      nBins_tau.push_back(1000); x_min_tau.push_back(0);    x_max_tau.push_back(1000);
   htau_names.push_back("taupt_fail");      nBins_tau.push_back(1000); x_min_tau.push_back(0);    x_max_tau.push_back(1000);
   htau_names.push_back("tau_MVA");         nBins_tau.push_back(200);  x_min_tau.push_back(-1);   x_max_tau.push_back(1);

   vector<TString> dms;
   dms.push_back("DM0");
   dms.push_back("DM1");
   dms.push_back("DM10");

   vector<TString> eta;
   eta.push_back("barrel");
   eta.push_back("endcap");

   vector<TH1F*> htau[histo_names.size()][dms.size()];
   for (unsigned int i = 0; i<htau_names.size(); ++i) {
     for (unsigned int j = 0; j<dms.size(); ++j) {
       for (unsigned int k = 0; k<eta.size(); ++k) {
	 htau[i][j].push_back( new TH1F(htau_names[i]+"_"+dms[j]+"_"+eta[k], htau_names[i]+"_"+dms[j]+"_"+eta[k], nBins_tau[i], x_min_tau[i], x_max_tau[i]) ); 
       }
     }
   }


   TH1F* h_reweight = new TH1F("h_r", "h_r", 100, -2, 2);

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

      float final_weight = 1;

      float pu_weight = 1;
      if (!data) {
        pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nick, "Data_2017BtoF");
      }
      


      //Is one of the triggers fired?
      bool PassMuonTrigger = false;
      if (trig_HLT_IsoMu27_accept) PassMuonTrigger = true;
      if (!PassMuonTrigger) continue;


      //start loop over reconstructed muons
      bool found_mumu_pair = false;
      for (unsigned int iMu1 = 0; iMu1 < mu_gt_pt->size(); ++iMu1) {
	//start 2nd loop over reconstructed mus
	for (unsigned int iMu2 = 0; iMu2 < iMu1; ++iMu2) {

	  if (mu_gt_pt->at(iMu1) < 30.0) continue;
	  if (fabs(mu_gt_eta->at(iMu1)) > 2.4) continue;
	  if (!mu_isPFMuon->at(iMu1)) continue;
	  if (!mu_isMediumMuon->at(iMu1)) continue; //medium ID
	  if (fabs(mu_gt_dxy_firstPVtx->at(iMu1)) > 0.045) continue;
	  if (fabs(mu_gt_dz_firstPVtx->at(iMu1)) > 0.2) continue;
	  float reliso = mu_pfIsoDbCorrected04->at(iMu1);
	  if (reliso > 0.15) continue;

	  if (mu_gt_pt->at(iMu2) < 30.0) continue;
	  if (fabs(mu_gt_eta->at(iMu2)) > 2.4) continue;
	  if (!mu_isPFMuon->at(iMu2)) continue; //medium ID
	  if (!mu_isMediumMuon->at(iMu2)) continue; //medium ID
	  if (fabs(mu_gt_dxy_firstPVtx->at(iMu2)) > 0.045) continue;
	  if (fabs(mu_gt_dz_firstPVtx->at(iMu2)) > 0.2) continue;
	  reliso = mu_pfIsoDbCorrected04->at(iMu2);
	  if (reliso > 0.15) continue;

	  TLorentzVector mu1_p4, mu2_p4, total_p4, met_p4, metmu_p4;
	  mu1_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu1), mu_gt_eta->at(iMu1), mu_gt_phi->at(iMu1), mu_mass);
	  mu2_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu2), mu_gt_eta->at(iMu2), mu_gt_phi->at(iMu2), mu_mass);
	  total_p4 = mu1_p4 + mu2_p4;

	  if (total_p4.M() < 60 || total_p4.M() > 120) continue;

	  if (!data) {
            final_weight = mc_w_sign*GetReweight_mumu(mu1_p4.Pt(), mu1_p4.Eta(), mu2_p4.Pt(), mu2_p4.Eta())*pu_weight;
	  }
	  
	  met_p4.SetPtEtaPhiM(MET_nominal_Pt, 0, MET_nominal_phi, 0);
	  metmu_p4 = met_p4 + mu1_p4;

	  if (final_weight != final_weight) continue;


	  float dR = mu1_p4.DeltaR(mu2_p4);
	  h[9]->Fill(dR, final_weight);
	  if (dR < 0.5) continue;

	  h_reweight->Fill(final_weight);


	  for (unsigned int iTau = 0; iTau < tau_pt->size(); ++iTau) {
	    if (tau_pt->at(iTau) < 20.0) continue;
	    if (fabs(tau_eta->at(iTau)) > 2.3) continue;
	    if (tau_decayModeFinding->at(iTau) < 0.5) continue;
	    if (tau_againstMuonTight3->at(iTau) < 0.5) continue;
	    if (tau_againstElectronVLooseMVA6->at(iTau) < 0.5) continue;
	    
	    TLorentzVector tau_p4;
	    tau_p4.SetPxPyPzE(tau_px->at(iTau), tau_py->at(iTau), tau_pz->at(iTau), tau_px->at(iTau));
	    if (tau_p4.DeltaR(mu1_p4) < 0.5) continue;
	    if (tau_p4.DeltaR(mu2_p4) < 0.5) continue;

	    //mu histos
	    h[0]->Fill(mu_gt_pt->at(iMu1), final_weight);
	    h[0]->Fill(mu_gt_pt->at(iMu2), final_weight);

	    h[1]->Fill(mu_gt_eta->at(iMu1), final_weight);
	    h[1]->Fill(mu_gt_eta->at(iMu2), final_weight);

	    h[2]->Fill(mu_gt_phi->at(iMu1), final_weight);
	    h[2]->Fill(mu_gt_phi->at(iMu2), final_weight);
	    

	    //Mu1 histos
	    h[3]->Fill(mu_gt_pt->at(iMu1), final_weight);
	    h[4]->Fill(mu_gt_eta->at(iMu1), final_weight);
	    h[5]->Fill(mu_gt_phi->at(iMu1), final_weight);


	    //Mu2 histos
	    h[6]->Fill(mu_gt_pt->at(iMu2), final_weight);
	    h[7]->Fill(mu_gt_eta->at(iMu2), final_weight);
	    h[8]->Fill(mu_gt_phi->at(iMu2), final_weight);

	    
	    //misc. hitos
	    h[12]->Fill(total_p4.M(), final_weight);
	    h[14]->Fill(MET_nominal_Pt, final_weight);
	    h[15]->Fill(MET_nominal_phi, final_weight);
	    h[16]->Fill(pv_n, final_weight);

	    int j_dm = -1, k_eta = -1;
	    if (tau_decayMode->at(iTau) == 0) {
	      j_dm = 0;
	    }
	    else if (tau_decayMode->at(iTau) == 1) {
	      j_dm = 1;
	    }
	    else if (tau_decayMode->at(iTau) == 10) {
	      j_dm = 2;
	    }
	    if (fabs(tau_eta->at(iTau)) < 1.5) {
	      k_eta = 0;
	    }
	    else {
	      k_eta = 1;
	    }

	    //Tau histos
	    if (tau_byTightIsolationMVArun2v1DBoldDMwLT->at(iTau) > 0.5) htau[0][j_dm][k_eta]->Fill(tau_pt->at(iTau), final_weight);
	    if ((tau_byTightIsolationMVArun2v1DBoldDMwLT->at(iTau) < 0.5) && (tau_byVLooseIsolationMVArun2v1DBoldDMwLT->at(iTau) > 0.5)) htau[1][j_dm][k_eta]->Fill(tau_pt->at(iTau), final_weight);
	    htau[2][j_dm][k_eta]->Fill(tau_byIsolationMVArun2v1DBoldDMwLTraw->at(iTau), final_weight);
	  }//loop over taus
	}//loop over mus
      }//loop over muons
   }//loop over events

   file_out->cd();
   //hCounter->Write();
   //hCounter2->Write();
   h_reweight->Write();
   for (unsigned int i = 0; i<histo_names.size(); ++i) h[i]->Write();
   for (unsigned int i=0; i<htau_names.size(); ++i) for (unsigned int j=0; j<dms.size(); ++j) for (unsigned int k=0; k<eta.size(); ++k) htau[i][j][k]->Write();
   file_out->Close();

}
