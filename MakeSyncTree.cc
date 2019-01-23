#define IIHEAnalysis_cxx
#include "IIHEAnalysis.h"
#include "PU_reWeighting.cc"
//#include <TH1.h>
#include <TLorentzVector.h>
//#include <TCanvas.h>
#include "TString.h"
//#include "MC_pileup_weight2017.cc"
#include <iostream>
#include <fstream>
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
  TTree* tree = (TTree*) fIn->Get("IIHEAnalysis");


  IIHEAnalysis* a = new IIHEAnalysis(tree);
  a->Loop(phase, type, inname, out_name, mc_nickname);
  return 0;
}

void IIHEAnalysis::Loop(string phase, string type_of_data, string in_name, string out_name, string mc_nickname) {
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

   bool isWJetsphase, isQCDphase, isfinalphase, isAntiMuphase, isFakesphase, isWJetsFakephase;
   cout << phase << endl;
   cout << DY << endl;
   if (phase == "QCD") {
     isQCDphase = true;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "WJets") {
     isQCDphase = false;
     isWJetsphase = true;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "WJetsFake") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = true;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "final") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = true;
     isAntiMuphase = false;
     isFakesphase = false;
   }
   else if (phase == "AntiMu") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = true;
     isFakesphase = false;
   }
   else if (phase == "Fakes") {
     isQCDphase = false;
     isWJetsphase = false;
     isWJetsFakephase = false;
     isfinalphase = false;
     isAntiMuphase = false;
     isFakesphase = true;
   }
   else {
     cout << phase << " does not exist!!!!" << endl;
   }

   TFile* file_out = new TFile(out_name.c_str(),"RECREATE");


   //newtree variables definitions
   TTree* newtree = new TTree("SyncTree", "SyncTree");
   // Declaration of leaf types
   Int_t           run;	      
   Int_t           lumi;      
   ULong64_t       evt;	      
   Float_t         pt_1;      
   Float_t         eta_1;     
   Float_t         phi_1;     
   Float_t         iso_mu;    
   Float_t         pt_2;      
   Float_t         phi_2;     
   Float_t         eta_2;     
   Float_t         charge_tau; 
   Float_t         dz_tau;     
   Float_t         TransverseMass;	       
   Float_t         PZeta;      
   Float_t         PZetaVis;      
   Float_t         dR;	       
   Float_t         met;	       
   Float_t         metphi;	       
   Float_t         weight;    

   // List of branches
   TBranch        *b_run        = newtree->Branch("run", &run);              
   TBranch        *b_lumi  	= newtree->Branch("lumi", &lumi);            
   TBranch        *b_evt   	= newtree->Branch("evt", &evt);              
   TBranch        *b_pt_1  	= newtree->Branch("pt_1", &pt_1);            
   TBranch        *b_eta_1 	= newtree->Branch("eta_1", &eta_1);          
   TBranch        *b_phi_1 	= newtree->Branch("phi_1", &phi_1);          
   TBranch        *b_iso_mu	= newtree->Branch("iso_mu", &iso_mu);        
   TBranch        *b_pt_2  	= newtree->Branch("pt_2", &pt_2);            
   TBranch        *b_phi_2 	= newtree->Branch("phi_2", &phi_2);          
   TBranch        *b_eta_2	= newtree->Branch("eta_2", &eta_2);          
   TBranch        *b_charge_tau = newtree->Branch("charge_tau", &charge_tau);
   TBranch        *b_dz_tau	= newtree->Branch("dz_tau", &dz_tau);        
   TBranch        *b_TransverseMass = newtree->Branch("TransverseMass", &TransverseMass);                
   TBranch        *b_PZeta	= newtree->Branch("PZeta", &PZeta);          
   TBranch        *b_PZetaVis	= newtree->Branch("PZetaVis", &PZetaVis);          
   TBranch        *b_dR	        = newtree->Branch("dR", &dR);                
   TBranch        *b_met	= newtree->Branch("met", &met);                
   TBranch        *b_metphi	= newtree->Branch("metphi", &metphi);                
   TBranch        *b_weight     = newtree->Branch("weight", &weight);        



   const float mu_mass = 0.10565837;



   Long64_t nEntries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int print_count = 0;
   long n1=0, n2=0;
   //start loop over all events
   for (Long64_t jEntry = 0; jEntry < nEntries; ++jEntry) {
      Long64_t iEntry = LoadTree(jEntry);
      if (iEntry < 0) break;
      if (jEntry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", jEntry, nEntries);
      //if (jEntry % 10 == 0) cout << endl << "Processed events:" << jEntry << " of " <<  nEntries;
      
      nb = fChain->GetEntry(jEntry);
      nbytes += nb;


      float pu_weight = 1, pu_weight_high = 1, pu_weight_low = 1;
      float first_weight = pu_weight;
      float reweight_njets = 1.0;
      


      //Is one of the triggers fired?
      bool PassMuonTrigger = false;
      if (trig_HLT_IsoMu27_accept) PassMuonTrigger = true;
      if (!PassMuonTrigger) continue;



      //start muon counting loop
      int Nmu = 0;
      bool dimuon = false;
      for (unsigned int iMu = 0; iMu < mu_gt_pt->size(); ++iMu) {
	TLorentzVector mu1_p4, mu2_p4;
        /*if(mu_isPFMuon->at(iMu) && mu_gt_pt->at(iMu) > 10 && fabs(mu_gt_eta->at(iMu)) < 2.4 && fabs(mu_gt_dxy_firstPVtx->at(iMu)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu) < 0.3 && mu_isMediumMuon->at(iMu)) ++Nmu;
	  if (Nmu > 1) break;*/
	
	if (mu_gt_pt->at(iMu) > 15 && fabs(mu_gt_eta->at(iMu)) < 2.4 && mu_isGlobalMuon->at(iMu) && mu_isTrackerMuon->at(iMu) && mu_isPFMuon->at(iMu) &&fabs(mu_gt_dxy_firstPVtx->at(iMu)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu) < 0.3) {
	  mu1_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu), mu_gt_eta->at(iMu), mu_gt_phi->at(iMu), mu_mass);
	  for (unsigned int iMu2= 0; iMu2 < iMu; ++iMu2) {
	    if (mu_gt_pt->at(iMu2) > 15 && fabs(mu_gt_eta->at(iMu2)) < 2.4 && mu_isGlobalMuon->at(iMu2) && mu_isTrackerMuon->at(iMu2) && mu_isPFMuon->at(iMu2) &&fabs(mu_gt_dxy_firstPVtx->at(iMu2)) < 0.045 && fabs(mu_gt_dz_firstPVtx->at(iMu2)) < 0.2 && mu_pfIsoDbCorrected04->at(iMu2) < 0.3) {
	      mu2_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu2), mu_gt_eta->at(iMu2), mu_gt_phi->at(iMu2), mu_mass);
	      if (mu1_p4.DeltaR(mu2_p4) > 0.15) dimuon = true;
	    }
	  }
	}
	if (dimuon) break;
      }
      //if (Nmu > 1) continue; //2nd muon veto                             
      if (dimuon) continue;

      //electron veto
      bool electron = false;
      for (unsigned int iEle = 0; iEle < gsf_pt->size(); ++iEle) {
	if (gsf_VIDMVAMedium->at(iEle) && gsf_pt->at(iEle) > 10 && fabs(gsf_eta->at(iEle)) < 2.5 && fabs(gsf_dxy_firstPVtx->at(iEle)) < 0.045 && fabs(gsf_dz_firstPVtx->at(iEle)) < 0.2 && gsf_passConversionVeto->at(iEle) && gsf_nLostInnerHits->at(iEle) <= 1 && gsf_relIso->at(iEle) < 0.3) electron = true;
	if (electron) break;
      }
      if (electron) continue;//FIXME

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
	  //Move selections as far up as possible to speed things up...
	  //cuts which do not involve the tau pt or the MET
	  if (fabs(tau_eta->at(iTau)) > 2.3) continue;
	  if (tau_decayModeFinding->at(iTau) < 0.5) continue;
	  if (tau_againstMuonTight3->at(iTau) < 0.5) continue;
	  if (tau_againstElectronVLooseMVA6->at(iTau) < 0.5) continue;
	  if (fabs(tau_lead_dz->at(iTau)) > 0.2) continue;
	  if (fabs(tau_charge->at(iTau)) != 1) continue;

	
	  //misc. cuts
	  bool SS = false;
	  if (isfinalphase || isWJetsphase || isWJetsFakephase || isFakesphase) {
	    if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) > 0) continue; //SS veto (opposite for QCD estimation)
	    SS = false;
	  }
	  else if (isQCDphase) {
	    if (tau_charge->at(iTau) * mu_gt_charge->at(iMu) < 0) continue; //OS veto
	    SS = true;
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



	  //check if if it's DY bkg or signal
	  TLorentzVector tau_p4, met_p4;
	  float met_px = MET_eefix_Px;
	  float met_py = MET_eefix_Py;
	  float met_pt = MET_eefix_Pt;
	  tau_p4.SetPtEtaPhiE(tau_pt->at(iTau), tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau));
	  met_p4.SetPxPyPzE(met_px, met_py, 0, met_pt);

	  
	  
	  map <TString, float> map_tau_TES;  map <TString, float> map_tau_TES_error;
	  map_tau_TES["DM0"] = 1.007;	       map_tau_TES_error["DM0"] = 0.008;	     
	  map_tau_TES["DM1"] = 0.998;	       map_tau_TES_error["DM1"] = 0.008;	     
	  map_tau_TES["DM10"] = 1.001;       map_tau_TES_error["DM10"] = 0.009;     
	  
	  TLorentzVector tau_TES_p4, tau_MES_p4, vis_p4, total_p4;
	  //met_p4.SetPxPyPzE(MET_Px, MET_eefix_Py, 0, MET_eefix_Pt);
	  

	  //determine decay mode and eta
	  TString str_dm = "", str_eta = "";
	  int m_dm = -1, n_eta = -1;
	  if (tau_decayMode->at(iTau) == 0) {
	    //m_dm = m_DM0;
	    str_dm = "DM0";
	  }
	  else if (tau_decayMode->at(iTau) == 1) {
	    //m_dm = m_DM1;
	    str_dm = "DM1";
	  }
	  else if (tau_decayMode->at(iTau) == 10) {
	    //m_dm = m_DM10;
	    str_dm = "DM10";
	  }
	  else {
	    continue;
	  }
	  if (fabs(tau_eta->at(iTau)) < 1.5) {
	    //n_eta = n_barrel;
	    str_eta = "barrel";
	  }
	  else {
	    //n_eta = n_endcap;
	    str_eta = "endcap";
	  }
	  n_eta = 0;
	  m_dm = 0;//FIXME



	  tau_p4.SetPtEtaPhiE(tau_pt->at(iTau), tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau));

	  vis_p4 = tau_p4 + mu_p4;
	  total_p4 = vis_p4;// + met_p4;
	  
	  if (tau_p4.Pt() < 20.0) continue;
	      
	      
	  //Pzeta calculation
	  float norm_zeta= norm_F( tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt(), tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt() );
	  //cout << norm_zeta << endl;
	  float x_zeta= (tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt())/norm_zeta;
	  float y_zeta= (tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt())/norm_zeta;
	  float p_zeta_mis=met_p4.Px()*x_zeta+met_p4.Py()*y_zeta;
	  float pzeta_vis=(tau_p4.Px()+mu_p4.Px())*x_zeta+(tau_p4.Py()+mu_p4.Py())*y_zeta;
	  bool cut_zeta= p_zeta_mis-0.85*pzeta_vis>-25;
	      
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
	      
	      
	  if (isfinalphase || isQCDphase || isAntiMuphase || isFakesphase) {
	    if (Mt > 50) Mt_accept = false;  //Mt cut, against Wjets (for Wjets CR, Mt>80)
	  }
	  else if (isWJetsphase || isWJetsFakephase) {
	    if (Mt < 80) Mt_accept = false;  //Mt cut, against Wjets (for Wjets CR, Mt>80)
	  }	      
	  else {
	    if (!isAntiMuphase) {
	      cout << "Error : non-existent phase!!!!" << endl;
	      break;
	    }
	  }	      
	  
	  float dR_mutau = tau_p4.DeltaR(mu_p4);
	  
	  float mc_weight = 1;
	  
	  float other_weights = 1;
	  reweight_njets = 1;
	  float final_weight = pu_weight*other_weights*reweight_njets*mc_weight;
	      
	      

	  if (tau_p4.DeltaR(mu_p4) < 0.5) continue;
	  if (!cut_zeta) continue;
	  
	  if(Mt_accept) {
	    if (tau_byTightIsolationMVArun2017v2DBoldDMwLT2017->at(iTau) < 0.5) continue;
	    found_mutau_pair = true;

	    run = ev_run;	      
	    lumi = ev_luminosityBlock;      
	    evt = ev_event;	      
	    pt_1 = mu_p4.Pt();      
	    eta_1 = mu_p4.Eta();     
	    phi_1 = mu_p4.Phi();     
	    iso_mu = reliso;    
	    pt_2 = tau_p4.Pt();      
	    eta_2 = tau_p4.Eta();     
	    phi_2 = tau_p4.Phi();     
	    charge_tau = tau_charge->at(iTau);
	    dz_tau = tau_lead_dz->at(iTau);    
	    TransverseMass = Mt;	      
	    PZeta = p_zeta_mis-0.85*pzeta_vis;     
	    PZetaVis = pzeta_vis;
	    dR = dR_mutau;	      
	    met = MET_eefix_Pt;
	    metphi = MET_eefix_phi;
	    weight = final_weight;    

	    newtree->Fill();
	  }

	}//loop over taus
      }//loop over muons
   }//loop over events

   file_out->cd();
   newtree->AutoSave();
   file_out->Close();
}
