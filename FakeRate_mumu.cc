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
  TTree* tree = (TTree*) fIn->Get("IIHEAnalysis");

  IIHEAnalysis* a = new IIHEAnalysis(tree);
  a->Loop(phase, type, inname, out_name, mc_nickname);
  return 0;
}

void IIHEAnalysis::Loop(string phase, string type_of_data, string in_name, string out_name, string mc_nick) {
   if (fChain == 0) return;
   
   cout << type_of_data << endl;

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


   vector<TString> passfail;
   passfail.push_back("pass");
   passfail.push_back("fail");


   vector<TString> dms;
   dms.push_back("DM0");
   dms.push_back("DM1");
   dms.push_back("DM10");

   vector<TString> eta;
   eta.push_back("barrel");
   eta.push_back("endcap");

   int nBins_tau = 1000;  float x_min_tau = 0;    float x_max_tau = 1000;

   vector<TH1F*> htau[HPS_WP.size()][passfail.size()][dms.size()];
   for (unsigned int i1 = 0; i1<HPS_WP.size(); ++i1) {
     for (unsigned int i2 = 0; i2<passfail.size(); ++i2) {
       for (unsigned int j = 0; j<dms.size(); ++j) {
	 for (unsigned int k = 0; k<eta.size(); ++k) {
	   htau[i1][i2][j].push_back( new TH1F("taupt_"+HPS_WP[i1]+"_"+passfail[i2]+"_"+dms[j]+"_"+eta[k], "taupt_"+HPS_WP[i1]+"_"+passfail[i2]+"_"+dms[j]+"_"+eta[k], nBins_tau, x_min_tau, x_max_tau) ); 
	 }
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
      //if (jEntry % 10 == 0) cout << endl << jEntry;

      nb = fChain->GetEntry(jEntry);
      nbytes += nb;

      float final_weight = 1;

      float pu_weight = 1;
      if (!data) {
        pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nick, "Data_METcorr_2017BtoF");
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

	  final_weight = 1;
	  if (!data) {
            final_weight = mc_w_sign*GetReweight_mumu(mu1_p4.Pt(), mu1_p4.Eta(), mu2_p4.Pt(), mu2_p4.Eta())*pu_weight;
	  }

	  met_p4.SetPtEtaPhiM(MET_nominal_Pt, 0, MET_nominal_phi, 0);
	  metmu_p4 = met_p4 + mu1_p4;

	  if (final_weight != final_weight) continue;


	  float dR = mu1_p4.DeltaR(mu2_p4);
	  h[9]->Fill(dR, final_weight);
	  if (dR < 0.5) continue;



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
	    h_reweight->Fill(final_weight);
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
	    else if (tau_decayMode->at(iTau) == 1 || tau_decayMode->at(iTau) == 2) {
	      j_dm = 1;
	    }
	    else if (tau_decayMode->at(iTau) == 10 || tau_decayMode->at(iTau) == 11) {
	      j_dm = 2;
	    }
	    if (fabs(tau_eta->at(iTau)) < 1.5) {
	      k_eta = 0;
	    }
	    else {
	      k_eta = 1;
	    }

	    vector<float> tauIDvalues;                                                              vector<int> specialValues;
	    tauIDvalues.push_back(tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(iTau));        specialValues.push_back(tauIDvalues.size()-1);
	    tauIDvalues.push_back(tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(iTau));
	    tauIDvalues.push_back(tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(iTau));

	    tauIDvalues.push_back(tau_byVVLooseIsolationMVArun2017v1DBoldDMwLT2017->at(iTau));      specialValues.push_back(tauIDvalues.size()-1);
	    tauIDvalues.push_back(tau_byVLooseIsolationMVArun2017v1DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byMediumIsolationMVArun2017v1DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byTightIsolationMVArun2017v1DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byVTightIsolationMVArun2017v1DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byVVTightIsolationMVArun2017v1DBoldDMwLT2017->at(iTau));

	    tauIDvalues.push_back(tau_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017->at(iTau));      specialValues.push_back(tauIDvalues.size()-1);
	    tauIDvalues.push_back(tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau));
	    tauIDvalues.push_back(tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau));

	    tauIDvalues.push_back(tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->at(iTau));          specialValues.push_back(tauIDvalues.size()-1);
	    tauIDvalues.push_back(tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT->at(iTau));

	    tauIDvalues.push_back(tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT->at(iTau));          specialValues.push_back(tauIDvalues.size()-1);
	    tauIDvalues.push_back(tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byTightIsolationMVArun2v1PWdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT->at(iTau));
	    tauIDvalues.push_back(tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT->at(iTau));


	    //Tau histos
	    int latestloosestID = -1;
	    for (unsigned int iID=0; iID<tauIDvalues.size(); ++iID) {
	      bool special = false;
	      for (unsigned int s=0; s<specialValues.size(); ++s) {
		if (iID == specialValues[s]) {
		  special = true;
		  latestloosestID = iID;
		  break;
		}
	      }

	      //loosest IDs have to be treated differently
	      if (special) {
		if (tauIDvalues[iID] > 0.5) {
		  htau[iID][0][j_dm][k_eta]->Fill(tau_pt->at(iTau), final_weight);
		}
		else {
		  htau[iID][1][j_dm][k_eta]->Fill(tau_pt->at(iTau), final_weight);
		}
	      }
	      else {
		if (tauIDvalues[latestloosestID] < 0.5) continue;
		if (tauIDvalues[iID] > 0.5) {
		  htau[iID][0][j_dm][k_eta]->Fill(tau_pt->at(iTau), final_weight);
		}
		else {
		  htau[iID][1][j_dm][k_eta]->Fill(tau_pt->at(iTau), final_weight);
		}
	      }
	    }		
	  }//loop over taus
	}//loop over mus
      }//loop over muons
   }//loop over events

   file_out->cd();
   h_reweight->Write();
   for (unsigned int i = 0; i<histo_names.size(); ++i) h[i]->Write();
   for (unsigned int i1=0; i1<HPS_WP.size(); ++i1) for (unsigned int i2=0; i2<passfail.size(); ++i2) for (unsigned int j=0; j<dms.size(); ++j) for (unsigned int k=0; k<eta.size(); ++k) htau[i1][i2][j][k]->Write();
   file_out->Close();

}
