#define IIHEAnalysis_cxx
#include "IIHEAnalysis.h"
#include "PU_reWeighting.cc"
//#include <TH1.h>
#include <TLorentzVector.h>
//#include <TCanvas.h>
#include "TString.h"
//#include "MC_pileup_weight2017.cc"
#include <iostream>
#include <vector>

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

   bool QCD=false;
   if (type_of_data == "QCD") QCD=true;

   bool WJets=false;
   if (type_of_data == "WJets" || type_of_data == "Wjets" || type_of_data == "WJetsinc" || type_of_data == "Wjetsinc") WJets=true;

   bool TT=false;
   if (type_of_data == "TT" || type_of_data == "TTinc") TT=true;

   bool isWJetsphase, isQCDphase, isfinalphase, isAntiMuphase;
   cout << phase << endl;
   cout << DY << endl;
   if (phase == "QCD") {
     isQCDphase = true;
     isWJetsphase = false;
     isfinalphase = false;
     isAntiMuphase = false;
   }
   else if (phase == "WJets") {
     isQCDphase = false;
     isWJetsphase = true;
     isfinalphase = false;
     isAntiMuphase = false;
   }
   else if (phase == "final") {
     isQCDphase = false;
     isWJetsphase = false;
     isfinalphase = true;
     isAntiMuphase = false;
   }
   else if (phase == "AntiMu") {
     isQCDphase = false;
     isWJetsphase = false;
     isfinalphase = false;
     isAntiMuphase = true;
   }
   else {
     cout << phase << " does not exist!!!!" << endl;
   }

   //string out_name = "out_"+type_of_data+".root";
   TFile* file_out = new TFile(out_name.c_str(),"RECREATE");

   const float mu_mass = 0.10565837;


   //list here the names and x-axis ranges of all gen-level histos we wish to create :
   vector<TString> histo_gen_names;                vector<int> nBins_gen;     vector<float> x_min_gen,    x_max_gen; 
   histo_gen_names.push_back("gen_tauh_vispt");    nBins_gen.push_back(100);  x_min_gen.push_back(0);     x_max_gen.push_back(100);
   histo_gen_names.push_back("gen_tauh_viseta");   nBins_gen.push_back(50);   x_min_gen.push_back(-2.5);  x_max_gen.push_back(2.5);
   histo_gen_names.push_back("gen_tauh_visphi");   nBins_gen.push_back(64);   x_min_gen.push_back(-3.2);  x_max_gen.push_back(3.2);
   histo_gen_names.push_back("gen_taumu_vispt");   nBins_gen.push_back(100);  x_min_gen.push_back(0);     x_max_gen.push_back(100);
   histo_gen_names.push_back("gen_taumu_viseta");  nBins_gen.push_back(50);   x_min_gen.push_back(-2.5);  x_max_gen.push_back(2.5);
   histo_gen_names.push_back("gen_taumu_visphi");  nBins_gen.push_back(64);   x_min_gen.push_back(-3.2);  x_max_gen.push_back(3.2);
   histo_gen_names.push_back("gen_ev_Mvis");       nBins_gen.push_back(150);  x_min_gen.push_back(0);     x_max_gen.push_back(150);
   histo_gen_names.push_back("gen_tau_decaymode");    nBins_gen.push_back(3);  x_min_gen.push_back(-0.5);     x_max_gen.push_back(2.5);

   vector<TH1F*> hgen;
   for (unsigned int i = 0; i<histo_gen_names.size(); ++i) hgen.push_back( new TH1F(histo_gen_names[i], histo_gen_names[i], nBins_gen[i], x_min_gen[i], x_max_gen[i]) ); 


   //list here the names and x-axis ranges of all reco-level histos we wish to create :
   vector<TString> histo_notauID_names;             vector<int> nBins;     vector<float> x_min,   x_max; 
   histo_notauID_names.push_back("tau_MVA");        nBins.push_back(200);  x_min.push_back(-1);    x_max.push_back(1);
   histo_notauID_names.push_back("ev_DRmutau");     nBins.push_back(100);  x_min.push_back(0);    x_max.push_back(10);
   histo_notauID_names.push_back("ev_Mt_raw");      nBins.push_back(150);  x_min.push_back(0);    x_max.push_back(150);
   histo_notauID_names.push_back("ev_Pzeta");       nBins.push_back(300);  x_min.push_back(-150); x_max.push_back(150);


   vector<TString> bkg_or_sig;
   bkg_or_sig.push_back("Sig");
   bkg_or_sig.push_back("Bkg");
   int two = bkg_or_sig.size();


   vector<TH1F*> hnotauID[two];
   for (unsigned int i = 0; i<histo_notauID_names.size(); ++i) {
     for (unsigned int k = 0; k<two; ++k) {
       hnotauID[k].push_back( new TH1F(bkg_or_sig[k]+"_"+histo_notauID_names[i], bkg_or_sig[k]+"_"+histo_notauID_names[i], nBins[i], x_min[i], x_max[i]) ); 
       hnotauID[k][i]->Sumw2();
     }
   }


   vector<TString> histo_taudependent_names;                      vector<int> nBins_tau;     vector<float> x_min_tau,   x_max_tau;
   histo_taudependent_names.push_back("ev_Mt");                   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("ev_Mvis");                 nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("ev_Mtot");                 nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("tau_pt");                  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("tau_eta");                 nBins_tau.push_back(50);   x_min_tau.push_back(-2.5); x_max_tau.push_back(2.5);
   histo_taudependent_names.push_back("tau_phi");                 nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   histo_taudependent_names.push_back("mu_pt");                   nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("mu_eta");                  nBins_tau.push_back(50);   x_min_tau.push_back(-2.5); x_max_tau.push_back(2.5);
   histo_taudependent_names.push_back("mu_phi");                  nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   histo_taudependent_names.push_back("ev_METmumass");            nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("ev_MET");                  nBins_tau.push_back(100);  x_min_tau.push_back(0);    x_max_tau.push_back(100);
   histo_taudependent_names.push_back("ev_METphi");               nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   histo_taudependent_names.push_back("ev_Nvertex");              nBins_tau.push_back(81);   x_min_tau.push_back(-0.5); x_max_tau.push_back(80.5);
   histo_taudependent_names.push_back("ev_Mvis_TESdown");         nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("ev_Mvis_TESup");           nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("tau_pt_TESdown");          nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("tau_pt_TESup");            nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   histo_taudependent_names.push_back("Z_pt");                    nBins_tau.push_back(300);  x_min_tau.push_back(0);    x_max_tau.push_back(300);
   histo_taudependent_names.push_back("n_jet");                   nBins_tau.push_back(20);   x_min_tau.push_back(0);    x_max_tau.push_back(20);
   if (isAntiMuphase) {
     histo_taudependent_names.push_back("tau_pt_SS");             nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
     histo_taudependent_names.push_back("ev_Mvis_SS");            nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
     histo_taudependent_names.push_back("mu_pt_SS");              nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   }

   vector<TString> tauIDs;
   tauIDs.push_back("noTauID");
   tauIDs.push_back("cutbased_loose_pass");
   tauIDs.push_back("cutbased_loose_fail");
   tauIDs.push_back("cutbased_medium_pass");
   tauIDs.push_back("cutbased_medium_fail");
   tauIDs.push_back("cutbased_tight_pass");
   tauIDs.push_back("cutbased_tight_fail");
   tauIDs.push_back("MVA_veryloose_pass");
   tauIDs.push_back("MVA_veryloose_fail");
   tauIDs.push_back("MVA_loose_pass");
   tauIDs.push_back("MVA_loose_fail");
   tauIDs.push_back("MVA_medium_pass");
   tauIDs.push_back("MVA_medium_fail");
   tauIDs.push_back("MVA_tight_pass");
   tauIDs.push_back("MVA_tight_fail");
   tauIDs.push_back("MVA_verytight_pass");
   tauIDs.push_back("MVA_verytight_fail");
   tauIDs.push_back("MVA_veryverytight_pass");
   tauIDs.push_back("MVA_veryverytight_fail");
   tauIDs.push_back("MVAnew_veryloose_pass");
   tauIDs.push_back("MVAnew_veryloose_fail");
   tauIDs.push_back("MVAnew_loose_pass");
   tauIDs.push_back("MVAnew_loose_fail");
   tauIDs.push_back("MVAnew_medium_pass");
   tauIDs.push_back("MVAnew_medium_fail");
   tauIDs.push_back("MVAnew_tight_pass");
   tauIDs.push_back("MVAnew_tight_fail");
   tauIDs.push_back("MVAnew_verytight_pass");
   tauIDs.push_back("MVAnew_verytight_fail");
   tauIDs.push_back("MVAnew_veryverytight_pass");
   tauIDs.push_back("MVAnew_veryverytight_fail");

   
   int nTauIDs = tauIDs.size();

   vector<TH1F*> htau[two][nTauIDs];
   for (unsigned int i = 0; i<histo_taudependent_names.size(); ++i) {
     for (unsigned int j = 0; j<nTauIDs ; ++j) {
       for (unsigned int k = 0; k<two; ++k) {
	 htau[k][j].push_back( new TH1F(bkg_or_sig[k]+"_"+histo_taudependent_names[i]+"_"+tauIDs[j], bkg_or_sig[k]+"_"+histo_taudependent_names[i]+"_"+tauIDs[j], nBins_tau[i], x_min_tau[i], x_max_tau[i]) ); 
	 htau[k][j][i]->Sumw2();
       }
     }
   }

   TH1F* h_reweight = new TH1F("h_reweight", "h_reweight", 100, 0, 2.5);
   TH1F* h_weight_afterallsel = new TH1F("h_weight_afterallsel", "h_weight_afterallsel", 100, -2, 2);
   TH1F* h_debug = new TH1F("h_debug", "h_debug", 100, 0, 5);

   Long64_t nEntries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int print_count = 0;
   long n1=0, n2=0;
   //start loop over all events
   for (Long64_t jEntry = 0; jEntry < nEntries; ++jEntry) {
      Long64_t iEntry = LoadTree(jEntry);
      if (iEntry < 0) break;
      if (jEntry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", jEntry, nEntries);

      nb = fChain->GetEntry(jEntry);
      nbytes += nb;

      float pu_weight = 1;
      if (!data) {
        pu_weight = PU_2017_Rereco::MC_pileup_weight(mc_trueNumInteractions, mc_nickname, "Rereco_all");
      }
      float first_weight = pu_weight;
      float reweight_njets = 1.0;
      
      bool gen_mutau = false, gen_mumu = false, unusualtau = false;
      vector<TLorentzVector> gentau_had_visp4, gentau_mu_visp4, genmu_p4;
      gentau_had_visp4.clear(), gentau_mu_visp4.clear(), genmu_p4.clear();
      if (!data) {
	
	vector<TLorentzVector> gentau_p4;
	vector<int> tau_ind, mu_ind;
	TLorentzVector p4;
	int moth_ind = -10;
	//start loop over all simulated particules
	for (unsigned int iMC = 0; iMC < mc_pt->size(); ++iMC) {
	  moth_ind = mc_mother_index->at(iMC).at(0);
	  if (moth_ind < 0) continue;
	  //Find out if the generated event is Z->mumu
	  if (abs(mc_pdgId->at(iMC)) == 13) {
	    if (abs(mc_pdgId->at(moth_ind)) != 15) {
	      mu_ind.push_back(iMC);
	      p4.SetPxPyPzE(mc_px->at(iMC), mc_py->at(iMC), mc_pz->at(iMC), mc_energy->at(iMC));
	      genmu_p4.push_back(p4);
	    }//end condition on mus' mothers
	  }//end condition on particle's id = mu


	  //Find taus whose mothers are either photons or Z bosons
	  if (abs(mc_pdgId->at(iMC)) == 15) {
	    tau_ind.push_back(iMC);
	    p4.SetPxPyPzE(mc_px->at(iMC), mc_py->at(iMC), mc_pz->at(iMC), mc_energy->at(iMC));
	    gentau_p4.push_back(p4);
	  }//end condition on particle's id = tau
	}//1st loop over sim particules

	
	int nTaus = tau_ind.size();
	int tau_dm[nTaus];
	TLorentzVector nutau_p4[nTaus], gentau_visp4[nTaus];
	for (unsigned int iTau = 0; iTau < nTaus; ++iTau) {
	  nutau_p4[iTau].SetPxPyPzE(0,0,0,0);
	  gentau_visp4[iTau].SetPxPyPzE(0,0,0,0);
	}
	for (unsigned int iTau = 0; iTau < nTaus; ++iTau) {
	  tau_dm[iTau] = 2;//we assume it's a hadronic tau by default
	  for (unsigned int iMC = 0; iMC < mc_pt->size(); ++iMC) {
	    moth_ind = mc_mother_index->at(iMC).at(0);
	    if (moth_ind < 0) continue;
	    if (moth_ind != tau_ind[iTau]) continue;
	    //now we now the mother of iMC's particle is a tau, we can classify that tau's decay mode
	    if (abs(mc_pdgId->at(iMC)) == 11) {
	      tau_dm[iTau] = 0;//electron tau
	      p4.SetPxPyPzE(mc_px->at(iMC), mc_py->at(iMC), mc_pz->at(iMC), mc_energy->at(iMC));
	      gentau_visp4[iTau] = p4;
	    }
	    else if (abs(mc_pdgId->at(iMC)) == 13) {
	      tau_dm[iTau] = 1;//muon tau
	      p4.SetPxPyPzE(mc_px->at(iMC), mc_py->at(iMC), mc_pz->at(iMC), mc_energy->at(iMC));
	      gentau_visp4[iTau] = p4;
	    }
	    else if (abs(mc_pdgId->at(iMC)) == 15) {
	      tau_dm[iTau] = -1;//tau becoming a tau!
	      p4.SetPxPyPzE(mc_px->at(iMC), mc_py->at(iMC), mc_pz->at(iMC), mc_energy->at(iMC));
	      gentau_visp4[iTau] = p4;
	    }
	    else if (abs(mc_pdgId->at(iMC)) < 11 || abs(mc_pdgId->at(iMC)) > 16) {
	      tau_dm[iTau] = 2;//hadronic tau
	      unusualtau = true;
	    }
	    else if (abs(mc_pdgId->at(iMC)) == 16) {
	      //cout << jEntry << " iTau " << iTau << "  Nu #" << iMC << "  Nu mother#" << moth_ind << "  Nu px " << mc_px->at(iMC) << endl; 
	      p4.SetPxPyPzE(mc_px->at(iMC), mc_py->at(iMC), mc_pz->at(iMC), mc_energy->at(iMC));
	      nutau_p4[iTau] = p4;
	    }
	  }
	}

	
	if (mu_ind.size() > 1) {
	  for (unsigned int iMu1 = 0; iMu1<mu_ind.size(); ++iMu1) {
	    for (unsigned int iMu2 = 0; iMu2<iMu1; ++iMu2) {
	      if (genmu_p4[iMu1].DeltaR(genmu_p4[iMu2]) > 0.5) {
		gen_mumu = true;
		break;
	      }
	    }
	    if (gen_mumu) break;
	  }
	}

	//getting vis tau_had p4 vector and tau_mu p4 vector
	for (unsigned int iTau1 = 0; iTau1 < nTaus; ++iTau1) {
	  if (print_count < 20) cout << jEntry << " iTau " << iTau1 << "  tau index " << tau_ind[iTau1] << "  decay mode " << tau_dm[iTau1] << endl;
	  if (tau_dm[iTau1] != 2) continue;//hadronic tau
	  gentau_visp4[iTau1] = gentau_p4[iTau1] - nutau_p4[iTau1];

	  for (unsigned int iTau2 = 0; iTau2 < nTaus; ++iTau2) {
	    if (tau_dm[iTau2] != 1) continue;//muon tau
	    gen_mutau = true;
	    gentau_had_visp4.push_back(gentau_visp4[iTau1]);
	    gentau_mu_visp4.push_back( gentau_visp4[iTau2]);
	  }
	}

	if (gentau_p4.size()>0 && print_count < 20) {
	  ++print_count;
	  cout << endl;
	  for (unsigned int iMC = 0; iMC < mc_pt->size(); ++iMC) cout << jEntry << " " << iMC << "  PDG ID " << mc_pdgId->at(iMC) << "  Mother PDG ID " << mc_mother_index->at(iMC).at(0) << "  px " << mc_px->at(iMC) << endl;
	  cout << gen_mutau << endl << endl;
	}


	//Now, depending on the tau decay modes, we fill the gen-level histograms
	for (unsigned int iMuTau = 0; iMuTau < gentau_had_visp4.size(); ++iMuTau) {
	  hgen[0]->Fill(gentau_had_visp4[iMuTau].Pt(), first_weight);
	  hgen[1]->Fill(gentau_had_visp4[iMuTau].Eta(), first_weight);
	  hgen[2]->Fill(gentau_had_visp4[iMuTau].Phi(), first_weight);
	  hgen[3]->Fill(gentau_mu_visp4[iMuTau].Pt(), first_weight);
	  hgen[4]->Fill(gentau_mu_visp4[iMuTau].Eta(), first_weight);
	  hgen[5]->Fill(gentau_mu_visp4[iMuTau].Phi(), first_weight);

	  p4 = gentau_had_visp4[iMuTau] + gentau_mu_visp4[iMuTau];
	  hgen[6]->Fill(p4.M(), first_weight);
	}
      }//end is this MC? condition

      if (WJets) {
	int njets = -2;
	for (unsigned int iLHE = 0; iLHE < LHE_Pt->size(); ++iLHE) {
          if (abs(LHE_pdgid->at(iLHE)) < 10 || abs(LHE_pdgid->at(iLHE)) == 21) ++njets;
	  //cout << LHE_pdgid->at(iLHE) << endl;
        }
	//cout << "njets  " << njets << endl;
	if (njets==0) {
	  reweight_njets = 111.024769762929;
	}
	else if (njets==1) {
	  reweight_njets = 14.68242648;
	}
	else if (njets==2) {
	  reweight_njets = 7.528580277;
	}
	else if (njets==3) {
	  reweight_njets = 2.436265125;
	}
	else if (njets==4) {
	  reweight_njets = 1.0599498;
	}
      }//end is this WJets? condition


      //Is one of the triggers fired?
      bool PassMuonTrigger = false;
      if (trig_HLT_IsoMu27_accept) PassMuonTrigger = true;
      if (!PassMuonTrigger) continue;


      //start muon counting loop
      int Nmu = 0;
      bool dimuon = false;
      for (unsigned int iMu = 0; iMu < mu_gt_pt->size(); ++iMu) {
	TLorentzVector mu1_p4, mu2_p4;
        if(mu_isPFMuon->at(iMu) && mu_gt_pt->at(iMu) > 10 && fabs(mu_gt_eta->at(iMu)) < 2.4 && fabs(mu_gt_dxy_firstPVtx->at(iMu)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu) < 0.3 && mu_isMediumMuon->at(iMu)) ++Nmu;
        if (Nmu > 1) break;
	
	if (mu_gt_pt->at(iMu) > 15 && fabs(mu_gt_eta->at(iMu)) < 2.4 && mu_isGlobalMuon->at(iMu) && mu_isTrackerMuon->at(iMu) && mu_isPFMuon->at(iMu) &&fabs(mu_gt_dxy_firstPVtx->at(iMu)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu) < 0.3) {
	  mu1_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu), mu_gt_eta->at(iMu), mu_gt_phi->at(iMu), mu_mass);
	  for (unsigned int iMu2= 0; iMu2 < mu_gt_pt->size(); ++iMu2) {
	    if (mu_gt_pt->at(iMu2) > 15 && fabs(mu_gt_eta->at(iMu2)) < 2.4 && mu_isGlobalMuon->at(iMu2) && mu_isTrackerMuon->at(iMu2) && mu_isPFMuon->at(iMu2) &&fabs(mu_gt_dxy_firstPVtx->at(iMu2)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu2)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu2) < 0.3) {
	      mu2_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu2), mu_gt_eta->at(iMu2), mu_gt_phi->at(iMu2), mu_mass);
	      if (mu1_p4.DeltaR(mu2_p4) > 0.15) dimuon = true;
	    }
	  }
	}
	if (dimuon) break;
      }
      if (Nmu > 1) continue; //2nd muon veto                                                                                                                                                                
      if (dimuon) continue;

      //electron veto
      bool electron = false;
      for (unsigned int iEle = 0; iEle < gsf_pt->size(); ++iEle) {
	if (gsf_VIDLoose->at(iEle) && gsf_pt->at(iEle) > 10 && fabs(gsf_eta->at(iEle)) < 2.5 && fabs(gsf_dxy_firstPVtx->at(iEle)) < 0.045 && fabs(gsf_dz_firstPVtx->at(iEle)) < 0.2 && gsf_passConversionVeto->at(iEle) && gsf_nLostInnerHits->at(iEle) <= 1 && gsf_relIso->at(iEle) < 0.3) electron = true;
        if (electron) break;
     }
      if (electron) continue;

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
      vector<int> orderedMu, orderedTau;
      vector<int> rest, rest2;
      bool interesting = false;


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
	    interesting = true;//just a flag to see if this actually happens
	  }
	  else {
	    rest2.push_back(rest[ii]);
	  }
	}
	orderedMu.push_back(lowest);
	rest = rest2;
      }


      //sorting Taus
      rest.clear();
      for (unsigned int ii = 0; ii < tau_pt->size(); ++ii) {
	rest.push_back(ii);
      }
      while (rest.size()>0) {
	rest2.clear();
	lowest = -1;
	iso = -1000;
	pt = 0;
	for (unsigned int ii = 0; ii < rest.size(); ++ii) {
	  if (tau_byIsolationMVArun2v1DBoldDMwLTraw->at(rest[ii]) > iso) {
	    iso = tau_byIsolationMVArun2v1DBoldDMwLTraw->at(rest[ii]);
	    pt = tau_pt->at(rest[ii]);
	    if (lowest > -1) rest2.push_back(lowest);
	    lowest = rest[ii];
	  }
	  else if (tau_byIsolationMVArun2v1DBoldDMwLTraw->at(rest[ii]) == iso && tau_pt->at(rest[ii]) > pt) {
	    pt = tau_pt->at(rest[ii]);
	    if (lowest > -1) rest2.push_back(lowest);
	    lowest = rest[ii];
	    interesting = true;
	  }
	  else {
	    rest2.push_back(rest[ii]);
	  }
	}
	orderedTau.push_back(lowest);
	rest = rest2;
      }



      /*if (print_count < 20 || interesting) {
	for (unsigned int ii = 0; ii < orderedMu.size(); ++ii) {
	  cout << "iso mu " << mu_pfIsoDbCorrected04->at(orderedMu[ii]) << "  pt mu " << mu_gt_pt->at(orderedMu[ii]) << endl;
	  for (unsigned int jj = 0; jj < orderedTau.size(); ++jj) cout << "                              iso tau " << tau_byIsolationMVArun2v1DBoldDMwLTraw->at(orderedTau[jj]) << "  pt tau " << tau_pt->at(orderedTau[jj]) << endl;
	}
	cout << endl << endl;
	}*/
	



      //start loop over reconstructed muons
      bool found_mutau_pair = false;
      for (unsigned int ii = 0; ii < orderedMu.size(); ++ii) {
	if (found_mutau_pair) break;
	int iMu = orderedMu[ii];
	if (mu_gt_pt->at(iMu) < 30.0) continue;
	if (fabs(mu_gt_eta->at(iMu)) > 2.4) continue;
	if (!mu_isPFMuon->at(iMu)) continue;
	if (!mu_isMediumMuon->at(iMu)) continue;
	if (fabs(mu_gt_dxy_firstPVtx->at(iMu)) >= 0.045) continue;
	if (fabs(mu_gt_dz_firstPVtx->at(iMu)) >= 0.2) continue;
	float reliso = mu_pfIsoDbCorrected04->at(iMu);
	if (isAntiMuphase) {
	  if (reliso < 0.15) continue;
	}
	else {
	  if (reliso > 0.15) continue;
	}

	TLorentzVector mu_p4;
	mu_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu), mu_gt_eta->at(iMu), mu_gt_phi->at(iMu), mu_mass);

	//start loop over reconstructed taus
	for (unsigned int jj = 0; jj < orderedTau.size(); ++jj) {
	  if (found_mutau_pair) break;
	  int iTau = orderedTau[jj];
	  //0 : TES down
	  //1 : TES up
	  //2 : no TES

	  vector<float> tauIDvalues;
	  tauIDvalues.push_back(tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(iTau));
	  tauIDvalues.push_back(tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(iTau));
	  tauIDvalues.push_back(tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(iTau));
	  tauIDvalues.push_back(tau_byVLooseIsolationMVArun2v1DBoldDMwLT->at(iTau));
	  tauIDvalues.push_back(tau_byLooseIsolationMVArun2v1DBoldDMwLT->at(iTau));
	  tauIDvalues.push_back(tau_byMediumIsolationMVArun2v1DBoldDMwLT->at(iTau));
	  tauIDvalues.push_back(tau_byTightIsolationMVArun2v1DBoldDMwLT->at(iTau));
	  tauIDvalues.push_back(tau_byVTightIsolationMVArun2v1DBoldDMwLT->at(iTau));
	  tauIDvalues.push_back(tau_byVVTightIsolationMVArun2v1DBoldDMwLT->at(iTau));
	  if (!QCD) {
	    tauIDvalues.push_back(tau_byVLooseIsolationMVArun2v1DBoldDMwLTNew->at(iTau));
	    tauIDvalues.push_back(tau_byLooseIsolationMVArun2v1DBoldDMwLTNew->at(iTau));
	    tauIDvalues.push_back(tau_byMediumIsolationMVArun2v1DBoldDMwLTNew->at(iTau));
	    tauIDvalues.push_back(tau_byTightIsolationMVArun2v1DBoldDMwLTNew->at(iTau));
	    tauIDvalues.push_back(tau_byVTightIsolationMVArun2v1DBoldDMwLTNew->at(iTau));
	    tauIDvalues.push_back(tau_byVVTightIsolationMVArun2v1DBoldDMwLTNew->at(iTau));
	  }

	  for (int iTES = 0; iTES < 3; ++iTES) {
	    TLorentzVector tau_p4, tau_TES_p4, vis_p4, met_p4, metmu_p4, total_p4;
	    //met_p4.SetPxPyPzE(MET_Px, MET_nominal_Py, 0, MET_nominal_Pt);
	    float met_px = MET_nominal_Px;
	    float met_py = MET_nominal_Py;
	    float met_pt = MET_nominal_Pt;
	    tau_p4.SetPtEtaPhiE(tau_pt->at(iTau), tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau));
	    met_p4.SetPxPyPzE(met_px, met_py, 0, met_pt);

	    //check if if it's DY bkg or signal
	    bool DY_sig = false;
	    float fakemu_reweight = 1;
	    if (DY) {
	      //we need to have a tau_mu tau_h pair at gen-level, the tau_h must have same sign and be DeltaR-compatible with the reco tau
	      if (gen_mutau) {
		for (unsigned int iGen = 0; iGen < gentau_had_visp4.size(); ++iGen) {
		  //cout << "dR  " << tau_p4.DeltaR(gentau_had_visp4[iGen]) << endl;
		  if (tau_p4.DeltaR(gentau_had_visp4[iGen]) < 0.5) {
		    DY_sig = true;
		    break;
		  }
		}
	      }
	      else {
		DY_sig = false;

	    	//reweighting for DY->mumu
		if (gen_mumu) {
		  bool fakemu_match = false;
		  for (unsigned int iGen = 0; iGen < genmu_p4.size(); ++iGen) {
		    if (tau_p4.DeltaR(genmu_p4[iGen]) < 0.5) {
		      fakemu_match = true;
		      break;
		    }
		  }
		  if (fabs(tau_p4.Eta())<0.4) fakemu_reweight = 1.200;
		  else if (fabs(tau_p4.Eta())<0.8) fakemu_reweight = 1.4;
		  else if (fabs(tau_p4.Eta())<1.2) fakemu_reweight = 1.09;
		  else if (fabs(tau_p4.Eta())<1.7) fakemu_reweight = 0.9;
		  else if (fabs(tau_p4.Eta())<2.3) fakemu_reweight = 2.3;
		}
		//cout << tau_p4.Eta() << "  " << first_weight << "  " << final_weight << "  " << final_weight/first_weight << endl;
	      }
	    }
	    //final_weight = first_weight;
	    
	    bool realtau = false;
	    bool TT_realtau = false;
	    if (!data) {
	      for (unsigned int iGen = 0; iGen < gentau_had_visp4.size(); ++iGen) {
		if (tau_p4.DeltaR(gentau_had_visp4[iGen]) > 0.5) continue;
		realtau = true;
	      }
	    }
	    if (TT && realtau) TT_realtau = true;

	    if (iTES == 0) {
	      if (realtau) {
		//TES down by 3%
		tau_TES_p4.SetPtEtaPhiE(tau_pt->at(iTau)*0.97, tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau)*0.97);
		met_p4 = met_p4 + tau_p4 - tau_TES_p4;
		tau_p4 = tau_TES_p4;
	      }
	    }
	    else if (iTES == 1) {
	      if (realtau) {
		//TES up by 3%
		tau_TES_p4.SetPtEtaPhiE(tau_pt->at(iTau)*1.03, tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau)*1.03);
		met_p4 = met_p4 + tau_p4 - tau_TES_p4;
		tau_p4 = tau_TES_p4;
	      }
	    }
	    vis_p4 = tau_p4 + mu_p4;
	    total_p4 = vis_p4;// + met_p4;
	    metmu_p4 = /*met_p4 +*/ mu_p4;
	    
	    if (tau_p4.Pt() < 20.0) continue;
	    if (fabs(tau_eta->at(iTau)) > 2.3) continue;
	    if (tau_decayModeFinding->at(iTau) < 0.5) continue;
	    if (tau_againstMuonTight3->at(iTau) < 0.5) continue;
	    if (tau_againstElectronVLooseMVA6->at(iTau) < 0.5) continue;
	    //if (tau_ptLeadChargedCand->at(iTau) < 5) continue;
	    if (!isAntiMuphase) {
	      if (fabs(tau_lead_dz->at(iTau)) > 0.2) continue;
	    }
	    if (fabs(tau_charge->at(iTau)) != 1) continue;
	    if (iTES==2) h_debug->Fill(0);



	    //Pzeta calculation
	    float norm_zeta= norm_F( tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt(), tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt() );
	    //cout << norm_zeta << endl;
	    float x_zeta= (tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt())/norm_zeta;
	    float y_zeta= (tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt())/norm_zeta;
	    float p_zeta_mis=met_p4.Px()*x_zeta+met_p4.Py()*y_zeta;
	    float pzeta_vis=(tau_p4.Px()+mu_p4.Px())*x_zeta+(tau_p4.Py()+mu_p4.Py())*y_zeta;
	    bool cut_zeta= p_zeta_mis-0.85*pzeta_vis>-25;

	    //misc. cuts
	    bool SS = false;
	    if (isfinalphase || isWJetsphase) {
	      if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) > 0) continue; //SS veto (opposite for QCD estimation)
	    }
	    else if (isQCDphase) {
	      if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) < 0) continue; //OS veto (opposite for QCD estimation)
	    }
	    else {
	      if (isAntiMuphase) {
		SS = false;
		if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) > 0) SS = true;
	      }
	      else {
		cout << "Error : non-existent phase!!!!" << endl;
		break;
	      }
	    }
	    if (iTES==2) h_debug->Fill(1);


	    if (first_weight != first_weight) continue;

	    float Mt;
	    bool Mt_accept = true;
	    if (2 * ( mu_gt_pt->at(iMu) * met_p4.Pt()  - mu_gt_px->at(iMu)*met_p4.Px() - mu_gt_py->at(iMu)*met_p4.Py() ) < 0) {
	      Mt = 0;
	      Mt_accept = false;
	    }
	    else {
	      Mt = sqrt(2 * ( mu_gt_pt->at(iMu) * met_p4.Pt()  - mu_gt_px->at(iMu)*met_p4.Px() - mu_gt_py->at(iMu)*met_p4.Py() ));
	    }


	    if (isfinalphase || isQCDphase || isAntiMuphase) {
	      if (Mt > 50) Mt_accept = false;  //Mt cut, against Wjets (for Wjets CR, Mt>80)
	    }
	    else if (isWJetsphase) {
	      if (Mt < 80) Mt_accept = false;  //Mt cut, against Wjets (for Wjets CR, Mt>80)
	    }	      
	    else {
	      if (!isAntiMuphase) {
		cout << "Error : non-existent phase!!!!" << endl;
		break;
	      }
	    }	      
	    if (iTES==2 && Mt_accept) h_debug->Fill(2);

	    float dR = tau_p4.DeltaR(mu_p4);

	    float other_weights = 1;
	    if (!data) other_weights = GetTriggerMuonIDMuonIsoReweight(mu_p4.Pt(), mu_p4.Eta());
	    float final_weight = pu_weight*other_weights*reweight_njets*fakemu_reweight;
	  
	    
	    if (DY && DY_sig) {
	      if (Mt_accept && iTES == 2) hnotauID[0][1]->Fill(dR, final_weight);
	      if (Mt_accept && iTES == 2) hnotauID[0][3]->Fill(p_zeta_mis-0.85*pzeta_vis, final_weight);
	      //if (Mt_accept && iTES == 2) hnotauID[0][6]->Fill(x_zeta*x_zeta + y_zeta*y_zeta, final_weight);
	    }
	    else {
	      if (Mt_accept && iTES == 2) hnotauID[1][1]->Fill(dR, final_weight);
	      if (Mt_accept && iTES == 2) hnotauID[1][3]->Fill(p_zeta_mis-0.85*pzeta_vis, final_weight);
	      }
	    
	    if (tau_p4.DeltaR(mu_p4) < 0.5) continue;
	    if (iTES==2 && Mt_accept) h_debug->Fill(3);
	    if (!cut_zeta) continue;



	    if(Mt_accept && iTES == 2 /*&& DY_sig*/) h_weight_afterallsel->Fill(final_weight), h_debug->Fill(4);	 
	    if(Mt_accept && iTES == 2) {
	      found_mutau_pair = true;
	    }
	    
	    bool fill_sig = false;
	    if ( isfinalphase && ((DY && DY_sig) || TT_realtau) ) fill_sig = true;

	    //fill "Sig" histograms
	    if (fill_sig) {
	      //mu histos
	      if (iTES == 2) hnotauID[0][2]->Fill(Mt, final_weight);
	      
	      //tau MVA histo
	      if (Mt_accept && iTES == 2) hnotauID[0][0]->Fill(tau_byIsolationMVArun2v1DBoldDMwLTraw->at(iTau), final_weight);
	    
	      //no tau ID
	      if (Mt_accept && iTES == 2) htau[0][0][0]->Fill(Mt, final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][1]->Fill(vis_p4.M(), final_weight);
	      if (Mt_accept && iTES == 0) htau[0][0][13]->Fill(vis_p4.M(), final_weight);
	      if (Mt_accept && iTES == 1) htau[0][0][14]->Fill(vis_p4.M(), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][2]->Fill(total_p4.M(), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][3]->Fill(tau_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 0) htau[0][0][15]->Fill(tau_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 1) htau[0][0][16]->Fill(tau_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][4]->Fill(tau_p4.Eta(), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][5]->Fill(tau_p4.Phi(), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][6]->Fill(mu_gt_pt->at(iMu), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][7]->Fill(mu_gt_eta->at(iMu), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][8]->Fill(mu_gt_phi->at(iMu), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][17]->Fill(vis_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][18]->Fill(jet_n, final_weight);
	      if (iTES == 2) htau[0][0][9]->Fill(metmu_p4.M(), final_weight);
	      if (iTES == 2) htau[0][0][10]->Fill(MET_nominal_Pt, final_weight);
	      if (iTES == 2) htau[0][0][11]->Fill(MET_nominal_phi, final_weight);
	      if (Mt_accept && iTES == 2) htau[0][0][12]->Fill(pv_n, final_weight);

	      for (unsigned int iValue=0; iValue<tauIDvalues.size(); ++iValue) {
		int iHisto = -1;

		if (tauIDvalues[iValue] > 0.5) {
		  iHisto = 2*iValue+1;
		}
		else {
		  iHisto = 2*iValue+2;
		}

	        if (Mt_accept && iTES == 2) htau[0][iHisto][0]->Fill(Mt, final_weight);
	        if (!SS) if (Mt_accept && iTES == 2) htau[0][iHisto][1]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && iTES == 0) htau[0][iHisto][13]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && iTES == 1) htau[0][iHisto][14]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && iTES == 2) htau[0][iHisto][2]->Fill(total_p4.M(), final_weight);
	        if (!SS) if (Mt_accept && iTES == 2) htau[0][iHisto][3]->Fill(tau_pt->at(iTau), final_weight);
		if (Mt_accept && iTES == 0) htau[0][iHisto][15]->Fill(tau_p4.Pt(), final_weight);
		if (Mt_accept && iTES == 1) htau[0][iHisto][16]->Fill(tau_p4.Pt(), final_weight);
	        if (Mt_accept && iTES == 2) htau[0][iHisto][4]->Fill(tau_p4.Eta(), final_weight);
	        if (Mt_accept && iTES == 2) htau[0][iHisto][5]->Fill(tau_p4.Phi(), final_weight);
		if (Mt_accept && iTES == 2) htau[0][iHisto][6]->Fill(mu_gt_pt->at(iMu), final_weight);
		if (Mt_accept && iTES == 2) htau[0][iHisto][7]->Fill(mu_gt_eta->at(iMu), final_weight);
		if (Mt_accept && iTES == 2) htau[0][iHisto][8]->Fill(mu_gt_phi->at(iMu), final_weight);
		if (Mt_accept && iTES == 2) htau[0][iHisto][17]->Fill(vis_p4.Pt(), final_weight);
		if (Mt_accept && iTES == 2) htau[0][iHisto][18]->Fill(jet_n, final_weight);
	        if (iTES == 2) htau[0][iHisto][9]->Fill(metmu_p4.M(), final_weight);
	        if (iTES == 2) htau[0][iHisto][10]->Fill(MET_nominal_Pt, final_weight);
	        if (iTES == 2) htau[0][iHisto][11]->Fill(MET_nominal_phi, final_weight);
	        if (Mt_accept && iTES == 2) htau[0][iHisto][12]->Fill(pv_n, final_weight);
	      }

	    }//end "Sig" histos
	    else {
	      //mu histos
	      if (iTES == 2) hnotauID[1][2]->Fill(Mt, final_weight);
	      
	      //tau MVA histo
	      if (Mt_accept && iTES == 2) hnotauID[1][0]->Fill(tau_byIsolationMVArun2v1DBoldDMwLTraw->at(iTau), final_weight);
	      
	      //various tau ID requirements
	      //no tau ID
	      if (Mt_accept && iTES == 2) htau[1][0][0]->Fill(Mt, final_weight);
	      if (!SS) if (Mt_accept && iTES == 2) htau[1][0][1]->Fill(vis_p4.M(), final_weight);
	      if (SS) if (Mt_accept && iTES == 2) htau[1][0][20]->Fill(vis_p4.M(), final_weight);
	      if (Mt_accept && iTES == 0) htau[1][0][13]->Fill(vis_p4.M(), final_weight);
	      if (Mt_accept && iTES == 1) htau[1][0][14]->Fill(vis_p4.M(), final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][2]->Fill(total_p4.M(), final_weight);
	      if (!SS) if (Mt_accept && iTES == 2) htau[1][0][3]->Fill(tau_p4.Pt(), final_weight);
	      if (SS) if (Mt_accept && iTES == 2) htau[1][0][19]->Fill(tau_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 0) htau[0][0][15]->Fill(tau_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 1) htau[0][0][16]->Fill(tau_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][4]->Fill(tau_p4.Eta(), final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][5]->Fill(tau_p4.Phi(), final_weight);
	      if (!SS) if (Mt_accept && iTES == 2) htau[1][0][6]->Fill(mu_gt_pt->at(iMu), final_weight);
	      if (SS) if (Mt_accept && iTES == 2) htau[1][0][21]->Fill(mu_gt_pt->at(iMu), final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][7]->Fill(mu_gt_eta->at(iMu), final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][8]->Fill(mu_gt_phi->at(iMu), final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][17]->Fill(vis_p4.Pt(), final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][18]->Fill(jet_n, final_weight);
	      if (iTES == 2) htau[1][0][9]->Fill(metmu_p4.M(), final_weight);
	      if (iTES == 2) htau[1][0][10]->Fill(MET_nominal_Pt, final_weight);
	      if (iTES == 2) htau[1][0][11]->Fill(MET_nominal_phi, final_weight);
	      if (Mt_accept && iTES == 2) htau[1][0][12]->Fill(pv_n, final_weight);
	      //if (Mt_accept && iTES == 2) htau[1][0][8]->Fill(tau_numberOfSignalChargedHadrCands->at(0), final_weight);
	      
	      for (unsigned int iValue=0; iValue<tauIDvalues.size(); ++iValue) {
		int iHisto = -1;
		if (tauIDvalues[iValue] > 0.5) {
		  iHisto = 2*iValue+1;
		}
		else {
		  iHisto = 2*iValue+2;
		}

		float w_weight = 1;
		if (WJets && !isWJetsphase) {
		  TFile* file_wjets = new TFile("Reweighting/Wjets_weights.root","R");
		  TH1F* Wjets_weights = (TH1F*) file_wjets->Get("ratio_"+tauIDs[iHisto]);
		  int w_bin = Wjets_weights->FindBin(vis_p4.M());
		  w_weight = Wjets_weights->GetBinContent(w_bin);
		  file_wjets->Close();
		}
		final_weight = pu_weight*other_weights*w_weight*reweight_njets;

	        if (Mt_accept && iTES == 2) htau[1][iHisto][0]->Fill(Mt, final_weight);
	        if (!SS) if (Mt_accept && iTES == 2) htau[1][iHisto][1]->Fill(vis_p4.M(), final_weight);
		if (SS) if (Mt_accept && iTES == 2) htau[1][iHisto][20]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && iTES == 0) htau[1][iHisto][13]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && iTES == 1) htau[1][iHisto][14]->Fill(vis_p4.M(), final_weight);
	        if (Mt_accept && iTES == 2) htau[1][iHisto][2]->Fill(total_p4.M(), final_weight);
	        if (!SS) if (Mt_accept && iTES == 2) htau[1][iHisto][3]->Fill(tau_pt->at(iTau), final_weight);
		if (SS) if (Mt_accept && iTES == 2) htau[1][iHisto][19]->Fill(tau_p4.Pt(), final_weight);
		if (Mt_accept && iTES == 0) htau[0][iHisto][15]->Fill(tau_p4.Pt(), final_weight);
		if (Mt_accept && iTES == 1) htau[0][iHisto][16]->Fill(tau_p4.Pt(), final_weight);
	        if (Mt_accept && iTES == 2) htau[1][iHisto][4]->Fill(tau_p4.Eta(), final_weight);
	        if (Mt_accept && iTES == 2) htau[1][iHisto][5]->Fill(tau_p4.Phi(), final_weight);
		if (!SS) if (Mt_accept && iTES == 2) htau[1][iHisto][6]->Fill(mu_gt_pt->at(iMu), final_weight);
		if (SS) if (Mt_accept && iTES == 2) htau[1][iHisto][21]->Fill(mu_gt_pt->at(iMu), final_weight);
		if (Mt_accept && iTES == 2) htau[1][iHisto][7]->Fill(mu_gt_eta->at(iMu), final_weight);
		if (Mt_accept && iTES == 2) htau[1][iHisto][8]->Fill(mu_gt_phi->at(iMu), final_weight);
		if (Mt_accept && iTES == 2) htau[0][iHisto][17]->Fill(vis_p4.Pt(), final_weight);
		if (Mt_accept && iTES == 2) htau[0][iHisto][18]->Fill(jet_n, final_weight);
	        if (iTES == 2) htau[1][iHisto][9]->Fill(metmu_p4.M(), final_weight);
	        if (iTES == 2) htau[1][iHisto][10]->Fill(MET_nominal_Pt, final_weight);
	        if (iTES == 2) htau[1][iHisto][11]->Fill(MET_nominal_phi, final_weight);
	        if (Mt_accept && iTES == 2) htau[1][iHisto][12]->Fill(pv_n, final_weight);
	      }

	    }//end background histos
	    
	  }//loop over TES
	}//loop over taus
      }//loop over muons
   }//loop over events

   file_out->cd();
   hCounter->Write();
   hCounter2->Write();
   h_reweight->Write();
   h_weight_afterallsel->Write();
   h_debug->Write();
   for (unsigned int i = 0; i<histo_notauID_names.size(); ++i) for (unsigned int k = 0; k<two; ++k) hnotauID[k][i]->Write();
   for (unsigned int i = 0; i<histo_taudependent_names.size(); ++i) for (unsigned int j = 0; j<nTauIDs; ++j) for (unsigned int k = 0; k<two; ++k) htau[k][j][i]->Write();
   if (DY) for (unsigned int i = 0; i<hgen.size(); ++i) hgen[i]->Write();
   file_out->Close();
}
