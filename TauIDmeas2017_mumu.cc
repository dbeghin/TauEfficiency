#define IIHEAnalysis_cxx
#include "IIHEAnalysis.h"
#include "PU_reWeighting.cc"
//#include <TH1.h>
#include <TLorentzVector.h>
//#include <TCanvas.h>
#include "TString.h"
#include <iostream>
#include "TRandom3.h"

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

void IIHEAnalysis::Loop(string phase, string type_of_data, string out_name, string mc_nickname, TH1F* hCounter, TH1F* hCounter2) {
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
   histo_names.push_back("ev_Mvis");        nBins.push_back(60);   x_min.push_back(60);    x_max.push_back(120);
   histo_names.push_back("ev_METmumass");   nBins.push_back(150);  x_min.push_back(0);    x_max.push_back(150);
   histo_names.push_back("ev_MET");         nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(100);
   histo_names.push_back("ev_METphi");      nBins.push_back(64);   x_min.push_back(-3.2); x_max.push_back(3.2);
   histo_names.push_back("ev_Nvertex");     nBins.push_back(81);   x_min.push_back(-0.5); x_max.push_back(80.5);
   histo_names.push_back("ev_Mvis_SS");     nBins.push_back(60);   x_min.push_back(60);   x_max.push_back(120);

   vector<TString> bkg_or_sig;
   bkg_or_sig.push_back("Bkg");
   bkg_or_sig.push_back("Sig");
   int two = bkg_or_sig.size();


   vector<TH1F*> h[two];
   for (unsigned int i = 0; i<histo_names.size(); ++i) {
     for (unsigned int k = 0; k<two; ++k) {
       h[k].push_back( new TH1F(bkg_or_sig[k]+"_"+histo_names[i], bkg_or_sig[k]+"_"+histo_names[i], nBins[i], x_min[i], x_max[i]) ); 
     }
   }


   TH1F* h_reweight = new TH1F("h_r", "h_r", 100, -2, 2);
   TH1F* h_events = new TH1F("h_events", "h_events", 1, 0, 1);

   Long64_t nEntries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int print_count = 0;
   long prev_jEntry = 0;
   //start loop over all events
   for (Long64_t jEntry = 0; jEntry < nEntries; ++jEntry) {
      Long64_t iEntry = LoadTree(jEntry);
      if (iEntry < 0) break;
      if (jEntry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", jEntry, nEntries);
      h_events->Fill(0);

      nb = fChain->GetEntry(jEntry);
      nbytes += nb;

      float pu_weight = 1;
      if (!data) {
	pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Data_2017BtoF");
	//pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Data_2017BtoF_high");
	//pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Data_2017BtoF_low");
	//pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Data_2017BtoF_80mb");
      }



      //Is one of the triggers fired?
      bool PassMuonTrigger = false;
      if (trig_HLT_IsoMu27_accept) PassMuonTrigger = true;
      if (!PassMuonTrigger) continue;


      //start muon counting loop
      int Nmu = 0;
      for (unsigned int iMu = 0; iMu < mu_gt_pt->size(); ++iMu) {
        if(mu_isPFMuon->at(iMu) && mu_gt_pt->at(iMu) > 10 && fabs(mu_gt_eta->at(iMu)) < 2.4 && fabs(mu_gt_dxy_firstPVtx->at(iMu)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu) < 0.3 && mu_isMediumMuon->at(iMu)) ++Nmu;
        if (Nmu > 2) break;
      }
      if (Nmu > 2) continue; //3rd muon veto

      //electron veto
      bool electron = false;
      for (unsigned int iEle = 0; iEle < gsf_pt->size(); ++iEle) {
        if (gsf_VIDLoose->at(iEle) && gsf_pt->at(iEle) > 10 && fabs(gsf_eta->at(iEle)) < 2.5 && fabs(gsf_dxy_firstPVtx->at(iEle)) < 0.045 && fabs(gsf_dz_firstPVtx->at(iEle)) < 0.2 && gsf_passConversionVeto->at(iEle) && gsf_nLostInnerHits->at(iEle) <= 1 && gsf_relIso->at(iEle) < 0.3) electron = true;
        if (electron) break;
      }
      if (electron) continue; //FIXME

      //bjet veto (medium WP for the bjet)                                                                                                                           
      /*bool bjet = false;
      for (unsigned int iJet = 0; iJet < jet_pt->size(); ++iJet) {
        if (jet_CSVv2->at(iJet) > 0.800 && jet_pt->at(iJet) > 20 && fabs(jet_eta->at(iJet)) < 2.4) bjet = true;
        if (bjet) break;
      }
      if (bjet) continue;*/



      //Sort muons, taus, by increasing isolation/decreasing pt
      float iso = 1.5, pt = 0.0;
      int lowest = -1;
      vector<int> orderedMu;
      vector<int> rest, rest2;

      //sorting muons
      for (unsigned int ii = 0; ii < mu_gt_pt->size(); ++ii) {
        rest.push_back(ii);
      }
      while (rest.size()>0) {
        rest2.clear();
        lowest = -1;
        iso = 10000000;
        pt = 0;
        for (unsigned int ii = 0; ii < rest.size(); ++ii) {
          if (mu_pfIsoDbCorrected04->at(rest[ii]) < iso) {
            iso = mu_pfIsoDbCorrected04->at(rest[ii]);
            pt = mu_gt_pt->at(rest[ii]);
            if (lowest > -1) rest2.push_back(lowest);
            lowest = rest[ii];
          }
          else if (mu_pfIsoDbCorrected04->at(rest[ii]) == iso && mu_gt_pt->at(rest[ii]) > pt) {
            pt = mu_gt_pt->at(rest[ii]);
            if (lowest > -1) rest2.push_back(lowest);
            lowest = rest[ii];
          }
          else {
            rest2.push_back(rest[ii]);
          }
        }
        orderedMu.push_back(lowest);
        rest = rest2;
      }
      //FIXME
      
      //start loop over reconstructed muons
      bool bPassedSel = false;
      for (unsigned int ii = 0; ii < orderedMu.size(); ++ii) {
      //for (unsigned int ii = 0; ii < mu_gt_pt->size(); ++ii) {
	if (bPassedSel) break;

	int iMuA = orderedMu[ii];
	//int iMu1 = ii;

	//start 2nd loop over reconstructed mus
	for (unsigned int jj = ii+1; jj < orderedMu.size(); ++jj) {
	  //for (unsigned int jj = ii+1; jj < mu_gt_pt->size(); ++jj) {
	  if (bPassedSel) break;

	  int iMuB = orderedMu[jj];
	  //int iMu2 = jj;
	  int iMu1 = -1, iMu2 = -1;
	  if (mu_gt_pt->at(iMuA) > mu_gt_pt->at(iMuB)) {
	      iMu1 = iMuA;
	      iMu2 = iMuB;
	  }
	  else {
	      iMu1 = iMuB;
	      iMu2 = iMuA;
	  }
	      
	  if (mu_gt_pt->at(iMu1) < 30.0) continue;
	  if (fabs(mu_gt_eta->at(iMu1)) > 2.4) continue;
	  if (!mu_isPFMuon->at(iMu1)) continue; //medium ID
	  if (!mu_isMediumMuon->at(iMu1)) continue; //medium ID
	  if (fabs(mu_gt_dxy_firstPVtx->at(iMu1)) > 0.045) continue;
	  if (fabs(mu_gt_dz_firstPVtx->at(iMu1)) > 0.2) continue;
	  float reliso = mu_pfIsoDbCorrected04->at(iMu1);
	  if (reliso > 0.15) continue;

	  if (mu_gt_pt->at(iMu2) < 30.0) continue;//20.0
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
	  
	  met_p4.SetPtEtaPhiM(MET_nominal_Pt, 0, MET_nominal_phi, 0);
	  metmu_p4 = met_p4 + mu1_p4;


	  float other_weights = 1;
	  float mc_wweight = 1;
	  //pu_weight = 1; //FIXME
	  if (!data) other_weights = GetReweight_mumu(mu1_p4.Pt(), mu1_p4.Eta(), mu2_p4.Pt(), mu2_p4.Eta()), mc_wweight = mc_w_sign;
	  float final_weight = pu_weight*other_weights*mc_wweight;
	  //final_weight = 1; //FIXME
	  if (final_weight != final_weight) continue;

	  float Mt = pow(2 * mu_gt_pt->at(iMu1) * MET_nominal_Pt * (1 - cos(mu_gt_phi->at(iMu1) - MET_nominal_phi) ), 0.5);
	  h[1][10]->Fill(Mt, final_weight);
	  h[1][13]->Fill(metmu_p4.M(), final_weight);

	  float dR = mu1_p4.DeltaR(mu2_p4);
	  h[1][9]->Fill(dR, final_weight);
	  if (dR < 0.5) continue;

	  //misc. cuts
	  if (mu_gt_charge->at(iMu1) * mu_gt_charge->at(iMu2) > 0) {
	    h[1][17]->Fill(total_p4.M(), final_weight);
	    continue; //SS veto
	  }

	  bPassedSel = true;
	  h_reweight->Fill(final_weight);
	  //h_reweight->Fill(pu_weight);
	  //h_reweight->Fill(mc_trueNumInteractions);

	  h[1][0]->Fill(mu_gt_pt->at(iMu1), final_weight);
	  h[1][0]->Fill(mu_gt_pt->at(iMu2), final_weight);

	  h[1][1]->Fill(mu_gt_eta->at(iMu1), final_weight);
	  h[1][1]->Fill(mu_gt_eta->at(iMu2), final_weight);

	  h[1][2]->Fill(mu_gt_phi->at(iMu1), final_weight);
	  h[1][2]->Fill(mu_gt_phi->at(iMu2), final_weight);
	 

	  //Mu1 histos
	  h[1][3]->Fill(mu_gt_pt->at(iMu1), final_weight);
	  h[1][4]->Fill(mu_gt_eta->at(iMu1), final_weight);
	  h[1][5]->Fill(mu_gt_phi->at(iMu1), final_weight);


	  //Mu2 histos
	  h[1][6]->Fill(mu_gt_pt->at(iMu2), final_weight);
	  h[1][7]->Fill(mu_gt_eta->at(iMu2), final_weight);
	  h[1][8]->Fill(mu_gt_phi->at(iMu2), final_weight);

	 
	  //misc. hitos
	  h[1][11]->Fill(Mt, final_weight);
	  h[1][12]->Fill(total_p4.M(), final_weight);
	  h[1][14]->Fill(MET_nominal_Pt, final_weight);
	  h[1][15]->Fill(MET_nominal_phi, final_weight);
	  h[1][16]->Fill(pv_n, final_weight);

	}//loop over mus
      }//loop over muons
   }//loop over events

   file_out->cd();
   //hCounter->Write();
   //hCounter2->Write();
   h_reweight->Write();
   h_events->Write();
   for (unsigned int i = 0; i<histo_names.size(); ++i) for (unsigned int k = 0; k<two; ++k) h[k][i]->Write();
   file_out->Close();

}
