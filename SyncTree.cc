#define SyncTree_cxx
#include "SyncTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {
  string out = *(argv + 1);
  string out_name= out;
  string in = *(argv + 2);
  string inname= in;
  TFile *fIn = TFile::Open(inname.c_str());
  TTree* tree = (TTree*) fIn->Get("SyncTree");


  SyncTree* a = new SyncTree(tree);
  a->Loop(inname, out_name);
  return 0;
}

void SyncTree::Loop(string in_name, string out_name) {
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



   Long64_t nentries = fChain->GetEntriesFast();

   //start event loop
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if (jentry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", jentry, nentries);

      event_file << run << ", " << lumi << ", " << evt << endl;

   }
}
