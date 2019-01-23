#define IIHEAnalysis_cxx
#include "IIHEAnalysis.h"
#include "PU_reWeighting.cc"
#include <TLorentzVector.h>
#include "TString.h"
#include "TDirectory.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

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


   fstream event_file; //object of fstream class
    
   string out_event = out_name;
   int point_position = -1;
   if (out_event.find('.') != string::npos) {
     point_position = out_event.find('.');
   }
   else {
     cout << "output name error!!!" << endl;
   }
   out_event.erase(point_position);
   out_event += ".txt";
   cout << out_event << endl;

   //opening file in out(write) mode
   event_file.open(out_event,ios::out);
   event_file << "root file: " <<  in_name << endl;
   event_file << "run, LS, event" << endl;


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
   histo_notauID_names.push_back("ev_Pzeta");       nBins.push_back(300);  x_min.push_back(-150); x_max.push_back(150);



   vector<TH1F*> hnotauID;
   for (unsigned int i = 0; i<histo_notauID_names.size(); ++i) {
     hnotauID.push_back( new TH1F(histo_notauID_names[i], histo_notauID_names[i], nBins[i], x_min[i], x_max[i]) ); 
     hnotauID[i]->Sumw2();
   }

   vector<TString> htau_names;                     vector<int> nBins_tau;     vector<float> x_min_tau,   x_max_tau;
   htau_names.push_back("ev_Mt");                  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mvis");                nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_Mtot");                nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_pt");                 nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("tau_eta");                nBins_tau.push_back(50);   x_min_tau.push_back(-2.5); x_max_tau.push_back(2.5);
   htau_names.push_back("tau_phi");                nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   htau_names.push_back("tau_DM");                 nBins_tau.push_back(11);   x_min_tau.push_back(0);    x_max_tau.push_back(11);
   htau_names.push_back("mu_pt");                  nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("mu_eta");                 nBins_tau.push_back(50);   x_min_tau.push_back(-2.5); x_max_tau.push_back(2.5);
   htau_names.push_back("mu_phi");                 nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   htau_names.push_back("ev_METmumass");           nBins_tau.push_back(150);  x_min_tau.push_back(0);    x_max_tau.push_back(150);
   htau_names.push_back("ev_MET");                 nBins_tau.push_back(100);  x_min_tau.push_back(0);    x_max_tau.push_back(100);
   htau_names.push_back("ev_METphi");              nBins_tau.push_back(64);   x_min_tau.push_back(-3.2); x_max_tau.push_back(3.2);
   htau_names.push_back("ev_Nvertex");             nBins_tau.push_back(120);  x_min_tau.push_back(0);    x_max_tau.push_back(120);
   htau_names.push_back("Z_pt");                   nBins_tau.push_back(300);  x_min_tau.push_back(0);    x_max_tau.push_back(300);
   htau_names.push_back("mu_iso");                 nBins_tau.push_back(100);  x_min_tau.push_back(0);    x_max_tau.push_back(1);


   //Make the map relating the histo name to its position on the histo name list
   map<TString, int> htau_map;
   for (unsigned int iName=0; iName<htau_names.size(); ++iName) htau_map[htau_names[iName]] = iName;

   vector<TH1F*> htau;
   for (unsigned int i = 0; i<htau_names.size(); ++i) {
     TString nname_h = htau_names[i];
     htau.push_back( new TH1F(nname_h, nname_h, nBins_tau[i], x_min_tau[i], x_max_tau[i]) ); 
     htau[i]->Sumw2();
   }


   TH1F* h_reweight = new TH1F("h_reweight", "h_reweight", 100, 0, 2.5);
   TH1F* h_weight_afterallsel = new TH1F("h_weight_afterallsel", "h_weight_afterallsel", 100, -2, 2);
   TH1F* h_debug = new TH1F("h_debug", "h_debug", 100, 0, 5);


   vector<long> AndrewRun;
   vector<long> AndrewLS;
   vector<long> AndrewEvent;
   vector<double> AndrewMuPt;
   vector<double> AndrewMuEta;
   vector<double> AndrewMuPhi;
   vector<double> AndrewTauPt;
   vector<double> AndrewTauEta;
   vector<double> AndrewTauPhi;
   string line="";
   ifstream myfile ("different_totalinfo.csv");
   if (myfile.is_open()) {
     while ( getline (myfile,line) ) {
       vector<string> vec;
       boost::algorithm::split(vec, line, boost::is_any_of(","));
       AndrewRun.push_back( stol(vec[0]) );
       AndrewLS.push_back( stol(vec[1]) );
       AndrewEvent.push_back( stol(vec[2]) );
       AndrewMuPt.push_back( stod(vec[3]) );
       AndrewMuEta.push_back( stod(vec[4]) );
       AndrewMuPhi.push_back( stod(vec[5]) );
       AndrewTauPt.push_back( stod(vec[6]) );
       AndrewTauEta.push_back( stod(vec[7]) );
       AndrewTauPhi.push_back( stod(vec[8]) );
     }
     myfile.close();
   }
   else cout << "Unable to open file" << endl;


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

      
      //compare with Andrew
      ULong64_t runnbr = ev_run;
      ULong64_t LS = ev_luminosityBlock;
      ULong64_t eventnbr = ev_event;

      bool keep_event = false;
      long nAndrew = -1;
      for (unsigned long iAndrew=0; iAndrew < AndrewRun.size(); ++iAndrew) {
      	if (AndrewRun[iAndrew] == runnbr && AndrewLS[iAndrew] == LS && AndrewEvent[iAndrew] == eventnbr) {
      	  keep_event = true;
	  nAndrew=iAndrew;
      	  break;
      	}
      	if (AndrewRun[iAndrew] > runnbr) {
	  break;
	}
      }
      if (!keep_event) continue;
      cout << endl << "ev_run " << ev_run << endl << endl;




      //start loop over reconstructed muons
      bool found = false;
      for (unsigned int ii = 0; ii < mu_gt_pt->size(); ++ii) {
	if (found) break;
	int iMu = ii;
	TLorentzVector mu_p4;
	mu_p4.SetPtEtaPhiM(mu_gt_pt->at(iMu), mu_gt_eta->at(iMu), mu_gt_phi->at(iMu), mu_mass);

	bool matched_mu = false;
	if ((mu_p4.Pt() <= AndrewMuPt[nAndrew]+0.6) && (mu_p4.Pt() >= AndrewMuPt[nAndrew]-0.6)) {
	  if (mu_p4.Eta() <= AndrewMuEta[nAndrew]+0.01 && mu_p4.Eta() >= AndrewMuEta[nAndrew]-0.01) {
	    if (mu_p4.Phi() <= AndrewMuPhi[nAndrew]+0.01 && mu_p4.Phi() >= AndrewMuPhi[nAndrew]-0.01) {
	      matched_mu = true;
	    }
	  }
	}
	if (!matched_mu) continue;


	//start loop over reconstructed taus
	for (unsigned int jj = 0; jj < tau_pt->size(); ++jj) {
	  if (found) break;
	  int iTau = jj;


	  //check if if it's DY bkg or signal
	  TLorentzVector tau_p4, met_p4, vis_p4, total_p4;
	  float met_px = MET_eefix_Px;
	  float met_py = MET_eefix_Py;
	  float met_pt = MET_eefix_Pt;
	  met_p4.SetPxPyPzE(met_px, met_py, 0, met_pt);
	  tau_p4.SetPtEtaPhiE(tau_pt->at(iTau), tau_eta->at(iTau), tau_phi->at(iTau), tau_energy->at(iTau));

	  vis_p4 = tau_p4 + mu_p4;
	  total_p4 = vis_p4;// + met_p4;
	      
	  bool matched_tau = false;
	  if (tau_p4.Pt() <= AndrewTauPt[nAndrew]+0.6 && tau_p4.Pt() >= AndrewTauPt[nAndrew]-0.6)
	    if (tau_p4.Eta() <= AndrewTauEta[nAndrew]+0.01 && tau_p4.Eta() >= AndrewTauEta[nAndrew]-0.01)
	      if (tau_p4.Phi() <= AndrewTauPhi[nAndrew]+0.01 && tau_p4.Phi() >= AndrewTauPhi[nAndrew]-0.01)
		matched_tau = true;
	  if (!matched_tau) continue;
	  found = true;

	  //Pzeta calculation
	  float norm_zeta= norm_F( tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt(), tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt() );
	  //cout << norm_zeta << endl;
	  float x_zeta= (tau_p4.Px()/tau_p4.Pt()+mu_p4.Px()/mu_p4.Pt())/norm_zeta;
	  float y_zeta= (tau_p4.Py()/tau_p4.Pt()+mu_p4.Py()/mu_p4.Pt())/norm_zeta;
	  float p_zeta_mis=met_p4.Px()*x_zeta+met_p4.Py()*y_zeta;
	  float pzeta_vis=(tau_p4.Px()+mu_p4.Px())*x_zeta+(tau_p4.Py()+mu_p4.Py())*y_zeta;
	  bool cut_zeta= p_zeta_mis-0.85*pzeta_vis>-25;
	  if (!cut_zeta) continue;
	      
	  float first_weight=1;
	      
	  float Mt;
	  if (2 * ( mu_gt_pt->at(iMu) * met_p4.Pt()  - mu_gt_px->at(iMu)*met_p4.Px() - mu_gt_py->at(iMu)*met_p4.Py() ) < 0) {
	    Mt = 0;
	    continue;
	  }
	  else {
	    Mt = sqrt(2 * ( mu_gt_pt->at(iMu) * met_p4.Pt()  - mu_gt_px->at(iMu)*met_p4.Px() - mu_gt_py->at(iMu)*met_p4.Py() ));
	  }
	  if (Mt > 50) continue;
	      
	  float dR = tau_p4.DeltaR(mu_p4);
	      
	  float other_weights = 1;
	  if (!data) other_weights = GetTriggerMuonIDMuonIsoReweight(mu_p4.Pt(), mu_p4.Eta());
	  float final_weight = other_weights;
	      
	  final_weight = other_weights;
	  h_reweight->Fill(final_weight);

	      

	  hnotauID[1]->Fill(dR, final_weight);
	  hnotauID[2]->Fill(p_zeta_mis-0.85*pzeta_vis, final_weight);
	  //if (Mt_accept && iTES == 2) hnotauID[0][6]->Fill(x_zeta*x_zeta + y_zeta*y_zeta, final_weight);
	      
	  h_debug->Fill(3);
	      
	  //tau MVA histo
	  hnotauID[0]->Fill(tau_byIsolationMVArun2v1DBoldDMwLTraw->at(iTau), final_weight);
	  

	  //fill text file with the run number etc.
	  event_file << ev_run << ", " << ev_luminosityBlock << ", " << ev_event << endl; 
		

	  //normal histos
	  htau[htau_map["ev_Mt"]]->Fill(Mt, final_weight);
	  htau[htau_map["ev_Mvis"]]->Fill(vis_p4.M(), final_weight);
	  htau[htau_map["ev_Mtot"]]->Fill(total_p4.M(), final_weight);
	  htau[htau_map["tau_pt"]]->Fill(tau_p4.Pt(), final_weight);
	  htau[htau_map["tau_eta"]]->Fill(tau_p4.Eta(), final_weight);
	  htau[htau_map["tau_phi"]]->Fill(tau_p4.Phi(), final_weight);
	  htau[htau_map["tau_DM"]]->Fill(tau_decayMode->at(iTau), final_weight);
	  htau[htau_map["mu_pt"]]->Fill(mu_gt_pt->at(iMu), final_weight);
	  htau[htau_map["mu_eta"]]->Fill(mu_gt_eta->at(iMu), final_weight);
	  htau[htau_map["mu_phi"]]->Fill(mu_gt_phi->at(iMu), final_weight);
	  htau[htau_map["Z_pt"]]->Fill(total_p4.Pt(), final_weight);
	  htau[htau_map["ev_MET"]]->Fill(met_p4.Pt(), final_weight);
	  htau[htau_map["ev_METphi"]]->Fill(met_p4.Phi(), final_weight);
	  htau[htau_map["ev_Nvertex"]]->Fill(pv_n, final_weight);
	  htau[htau_map["mu_iso"]]->Fill(mu_pfIsoDbCorrected04->at(iMu), final_weight);

	}//loop over taus
      }//loop over muons
   }//loop over events

   event_file << endl << endl << endl;
   event_file.close();

   file_out->cd();
   h_reweight->Write();
   h_weight_afterallsel->Write();
   h_debug->Write();

   for (unsigned int i = 0; i<histo_notauID_names.size(); ++i) hnotauID[i]->Write();
   if (DY) for (unsigned int i = 0; i<hgen.size(); ++i) hgen[i]->Write();
   for (unsigned int i = 0; i<htau_names.size(); ++i) htau[i]->Write();
   file_out->Close();
}
