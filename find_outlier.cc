#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"


using namespace std;




int main(int argc, char** argv) {


  string file = *(argv + 1);
  TString file_name= file;

  TFile* file_in = new TFile(file_name, "R");
  TCanvas* c = new TCanvas("c","c",0,0,600,600);
  
  cout << file_name << endl;
  TString var = "Bkg_ev_Mvis_cutbased_medium_pass";
  TH1F* h = (TH1F*) file_in -> Get(var);
  c->cd();
  h->Draw();
  c->Modified();
  c->SaveAs("cccc.pdf");
  c->Close();
  if (h == 0) cout << "!!!!!!!" << endl;
  cout << h << endl << endl << endl;


  return 0;
}
