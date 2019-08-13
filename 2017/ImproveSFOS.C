// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"
#include <sstream>
#include <iostream>
#include <fstream>

// CMS3
#include "Functions.h"
#include "CMS3_WWW100.cc"
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1) {

  vector<myevt> e;
  addeventtocheck(e,1, 2842, 1443084);
  addeventtocheck(e,1,17341, 8807395);
  addeventtocheck(e,1,15849, 8049420);
  addeventtocheck(e,1,21682,11012507);

  bool blindSR = false;
  bool btagreweighting = true;
  bool applylepSF      = false;
  bool applytrigSF     = false;
  bool applyPUrewgt    = true;
  bool getJECunc       = false;
  
  const char* json_file = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();
  histonames.push_back("NJets_SSee");                           hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("NJets_SSem");                           hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("NJets_SSmm");                           hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("NJetsLoose_SSee");                      hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("NJetsLoose_SSem");                      hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("NJetsLoose_SSmm");                      hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("MTmax_SSee");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_SSee");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTsum_SSee");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTmax_SSem");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_SSem");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTsum_SSem");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTmax_SSmm");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_SSmm");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTsum_SSmm");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("DPhill_SSee");                          hbins.push_back(13); hlow.push_back(    0); hup.push_back(3.25);
  histonames.push_back("DPhill_SSem");                          hbins.push_back(13); hlow.push_back(    0); hup.push_back(3.25);
  histonames.push_back("DPhill_SSmm");                          hbins.push_back(13); hlow.push_back(    0); hup.push_back(3.25);
  histonames.push_back("DPhillMET_SSee");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillMET_SSem");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillMET_SSmm");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilMETmax_SSee");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilMETmax_SSem");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilMETmax_SSmm");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilMETmin_SSee");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilMETmin_SSem");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilMETmin_SSmm");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("pTll_SSee");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTll_SSem");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTll_SSmm");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);

  histonames.push_back("MT3rd_1SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_0SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_1SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_2SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmax_0SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmax_1SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmax_2SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTsum_0SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTsum_1SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTsum_2SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTleastZ_2SFOS");                       hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTlll_0SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTlll_1SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTlll_2SFOS");                          hbins.push_back(25); hlow.push_back(    0); hup.push_back(250);

  histonames.push_back("MSFOS_1SFOS");                          hbins.push_back(15); hlow.push_back(   45); hup.push_back(120);
  histonames.push_back("MSFOSZlike_2SFOS");                     hbins.push_back(12); hlow.push_back(   60); hup.push_back(120);
  histonames.push_back("MSFOSnoZ_2SFOS");                       hbins.push_back(12); hlow.push_back(   60); hup.push_back(120);

  histonames.push_back("pTSFOS_1SFOS");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pT3rd_1SFOS");                          hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTW_1SFOS");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTSFOSmax_2SFOS");                      hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTSFOSmin_2SFOS");                      hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTSFOSZlike_2SFOS");                    hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTSFOSnoZ_2SFOS");                      hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTleastZ_2SFOS");                       hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTWleastZ_2SFOS");                      hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTlllMET_0SFOS");                       hbins.push_back(15); hlow.push_back(    0); hup.push_back(300);
  histonames.push_back("pTlllMET_1SFOS");                       hbins.push_back(15); hlow.push_back(    0); hup.push_back(300);
  histonames.push_back("pTlllMET_2SFOS");                       hbins.push_back(15); hlow.push_back(    0); hup.push_back(300);

  histonames.push_back("DPhillofZ_1SFOS");                      hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillofZlike_2SFOS");                  hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilloflessZ_2SFOS");                  hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillSFOSmin_2SFOS");                  hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillSFOSmax_2SFOS");                  hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillofZvsl_1SFOS");                   hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillofZlikevsl_2SFOS");               hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillSFOSvslmax_2SFOS");               hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillSFOSvslmin_2SFOS");               hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillofZvslMET_1SFOS");                hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillofZlikevslMET_2SFOS");            hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillSFOSvslMETmax_2SFOS");            hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhillSFOSvslMETmin_2SFOS");            hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);


  histonames.push_back("MET_0SFOS");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MET_1SFOS");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MET_2SFOS");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("DPhilllMET_0SFOS");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilllMET_1SFOS");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilllMET_2SFOS");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("pTlll_0SFOS");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("pTlll_1SFOS");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("pTlll_2SFOS");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("Mlll_0SFOS");                           hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("Mlll_1SFOS");                           hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("Mlll_2SFOS");                           hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);

  histonames.push_back("Mll_SSee");                             hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("Mll_SSem");                             hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("Mll_SSmm");                             hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MET_SSee");                             hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MET_SSem");                             hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MET_SSmm");                             hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("Mjj_SSee");                             hbins.push_back(32); hlow.push_back(    0); hup.push_back(160);
  histonames.push_back("Mjj_SSem");                             hbins.push_back(32); hlow.push_back(    0); hup.push_back(160);
  histonames.push_back("Mjj_SSmm");                             hbins.push_back(32); hlow.push_back(    0); hup.push_back(160);
  histonames.push_back("MjjL_SSee");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(600);
  histonames.push_back("MjjL_SSem");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(600);
  histonames.push_back("MjjL_SSmm");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(600);
  histonames.push_back("DetajjL_SSee");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(5);
  histonames.push_back("DetajjL_SSem");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(5);
  histonames.push_back("DetajjL_SSmm");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(5);


  histonames.push_back("Mjjll_SSee");                         hbins.push_back(15); hlow.push_back(    0); hup.push_back(1500);
  histonames.push_back("Mjjll_SSem");                         hbins.push_back(15); hlow.push_back(    0); hup.push_back(1500);
  histonames.push_back("Mjjll_SSmm");                         hbins.push_back(15); hlow.push_back(    0); hup.push_back(1500);
  histonames.push_back("MTjjllMET_SSee");                     hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTjjllMET_SSem");                     hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTjjllMET_SSmm");                     hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MljDR_SSee");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MljDR_SSem");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MljDR_SSmm");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("DPhiljmin_SSee");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhiljmin_SSem");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhiljmin_SSmm");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DRljmin_SSee");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DRljmin_SSem");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DRljmin_SSmm");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("RatMETvsll_SSee");                    hbins.push_back(16); hlow.push_back(    0); hup.push_back(4);
  histonames.push_back("RatMETvsll_SSem");                    hbins.push_back(16); hlow.push_back(    0); hup.push_back(4);
  histonames.push_back("RatMETvsll_SSmm");                    hbins.push_back(16); hlow.push_back(    0); hup.push_back(4);
  histonames.push_back("DPhijjvsllMET_SSee");                 hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhijjvsllMET_SSem");                 hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhijjvsllMET_SSmm");                 hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("pTlljjMET_SSee");                     hbins.push_back(10); hlow.push_back(    0); hup.push_back(300);
  histonames.push_back("pTlljjMET_SSem");                     hbins.push_back(10); hlow.push_back(    0); hup.push_back(300);
  histonames.push_back("pTlljjMET_SSmm");                     hbins.push_back(10); hlow.push_back(    0); hup.push_back(300);
  histonames.push_back("ST_SSee");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(1000);
  histonames.push_back("ST_SSem");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(1000);
  histonames.push_back("ST_SSmm");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(1000);

  histonames.push_back("DRllsum_0SFOS");                      hbins.push_back(16); hlow.push_back(    0); hup.push_back(16);
  histonames.push_back("DRllsum_1SFOS");                      hbins.push_back(16); hlow.push_back(    0); hup.push_back(16);
  histonames.push_back("DRllsum_2SFOS");                      hbins.push_back(16); hlow.push_back(    0); hup.push_back(16);
  histonames.push_back("pTllsum_0SFOS");                      hbins.push_back(20); hlow.push_back(    0); hup.push_back(1000);
  histonames.push_back("pTllsum_1SFOS");                      hbins.push_back(20); hlow.push_back(    0); hup.push_back(1000);
  histonames.push_back("pTllsum_2SFOS");                      hbins.push_back(20); hlow.push_back(    0); hup.push_back(1000);
  histonames.push_back("RatMETvslll_0SFOS");                  hbins.push_back(10); hlow.push_back(    0); hup.push_back(2.5);
  histonames.push_back("RatMETvslll_1SFOS");                  hbins.push_back(10); hlow.push_back(    0); hup.push_back(2.5);
  histonames.push_back("RatMETvslll_2SFOS");                  hbins.push_back(10); hlow.push_back(    0); hup.push_back(2.5);
  histonames.push_back("ST_0SFOS");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(1000);
  histonames.push_back("ST_1SFOS");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(1000);
  histonames.push_back("ST_2SFOS");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(1000);

  
  unsigned int nhistotemp = histonames.size();
  for(unsigned int i = 0; i<nhistotemp; ++i){
    histonames.push_back("ZPreselect_"+histonames[i]); hbins.push_back(hbins[i]); hlow.push_back(hlow[i]); hup.push_back(hup[i]);
  }
  //for(unsigned int i = 0; i<histonames.size(); ++i) hbins[i] *=2;
  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir);
  cout << "Loaded histograms" << endl;

  unsigned int nEventsRunning = 0;
  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    string fname = currentFile->GetTitle();

    TTree *tree = (TTree*)file->Get("t");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    cms3.Init(tree);
    cout << "File " <<  currentFile->GetTitle() << endl;
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
 
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      cms3.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      CMS3::progress( nEventsTotal, nEventsChain );

      
      double weight = evt_scale1fb()*35.9;
      //if(nEventsTotal==0) cout << weight << endl; 

      if(firstgoodvertex()!=0)   continue;
      if(nVert()<0)              continue;
      //if(nlep()<2)               continue;

      //weight = 1;
      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("www_2l_mia")!=string::npos)      weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")!=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(weight>100) cout << weight << " " << currentFile->GetTitle() << endl;
      if(isData()) weight = 1.;
      double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      float PUweight(1.), PUweightup(1.), PUweightdn(1.);
      if(applyPUrewgt&&!isData()){
	PUweight = getPUWeightAndError(PUweightdn,PUweightup);
	weight *= PUweight;
      }
      
      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      LorentzVector MET_up; MET_up.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*TMath::Cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*TMath::Sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
      LorentzVector MET_dn; MET_dn.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*TMath::Cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*TMath::Sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());

      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);

      int nj_up(0), nb_up(0), nj30_up(0);
      if(getJECunc) getalljetnumbers(nj_up,nj30_up,nb_up,1);
      float Mjj_up = -1;
      float MjjL_up = -1; float Detajj_up = -1;
      if(getJECunc) getMjjAndDeta(Mjj_up,MjjL_up,Detajj_up,1);

      int nj_dn(0), nb_dn(0), nj30_dn(0);
      if(getJECunc) getalljetnumbers(nj_dn,nj30_dn,nb_dn,-1);
      float Mjj_dn = -1;
      float MjjL_dn = -1; float Detajj_dn = -1;
      if(getJECunc) getMjjAndDeta(Mjj_dn,MjjL_dn,Detajj_dn,-1);

      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      //getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l);//oldID
      //getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l,1,25,20);//newID
      getleptonindices_v2(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l,25,20);
      float lepSF(1.), lepSFerr(0.);//i3l and iSS have same ID
      if(applylepSF&&!isData()){
	lepSF = getlepSFWeightandError(lepSFerr,i3l,ia3l);
	weight *= lepSF;
      }
      float trigSF(1.), trigSFerr(1.);
      if(applytrigSF&&!isData()){
	trigSF    = getTriggerWeightandError(trigSFerr, i3l,ia3l);
	weight *= trigSF;
      }
      float weight_lepSFup = weight;
      float weight_lepSFdn = weight;
      float weight_PUup    = weight;
      float weight_PUdn    = weight;
      float weight_bHFSFup = weight;
      float weight_bHFSFdn = weight;
      float weight_bLFSFup = weight;
      float weight_bLFSFdn = weight;
      if(!isData()&&btagreweighting&&weight_btagsf()!=0){
	weight_bHFSFup *= weight_btagsf_heavy_UP()/weight_btagsf();
	weight_bHFSFdn *= weight_btagsf_heavy_DN()/weight_btagsf();
	weight_bLFSFup *= weight_btagsf_light_UP()/weight_btagsf();
	weight_bLFSFdn *= weight_btagsf_light_DN()/weight_btagsf();
      }
      if(applyPUrewgt&&!isData()&&PUweight!=0){
	weight_PUup *= PUweightup/PUweight;
	weight_PUdn *= PUweightdn/PUweight;
      }
      if(applylepSF&&!isData()&&lepSF!=0){
	weight_lepSFup *= (lepSF+lepSFerr)/lepSF;
	weight_lepSFdn *= (lepSF-lepSFerr)/lepSF;
      }

      int nvetoSS = vSS.size();
      int nveto3l = v3l.size();
      int nSS = iSS.size();
      int n3l = i3l.size();
      int nvetoaSS = vaSS.size();
      int nvetoa3l = va3l.size();
      int naSS = iaSS.size();
      int na3l = ia3l.size();
      
      if((n3l+na3l)<2) continue;
      bool passofflineforTrigger = passofflineTriggers(i3l, ia3l);
      if(!passofflineforTrigger) continue;
      
      if(isData()){
	if(!passFilters()) continue;
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) { continue; }
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = passonlineTriggers(i3l, ia3l);//currently applied only to data
	if(!passonlineTrigger) continue;
      }

      string sample   = skimFilePrefix;
      string sn       = ((iSS.size()+iaSS.size())>=2) ? process(fname,true ,iSS,iaSS) : string("not2l");
      string sn2      = ((i3l.size()+ia3l.size())>=3) ? process(fname,false,i3l,ia3l) : string("not3l");
      bool isphotonSS = (sn =="photonfakes");
      bool isphoton3l = (sn2=="photonfakes");
      if(splitVH(fname)){ sample = "WHtoWWW"; }



      float MTmax = -1;
      if(iSS.size()==2) MTmax = calcMTmax(iSS,MET);
      else if(iSS.size()==1&&iaSS.size()>=1){
	vector<int> temp; temp.push_back(iSS[0]); temp.push_back(iaSS[0]);
	MTmax = calcMTmax(temp,MET);
      }
      float MTmax3l = calcMTmax(i3l,MET,true);
      float MTmax_up(-1), MTmax_dn(-1), MTmax3l_up(-1), MTmax3l_dn(-1);
      if(getJECunc) {
	MTmax_up   = calcMTmax(iSS,MET_up);
	//MTmax3l_up = calcMTmax(i3l,MET_up,true);
	MTmax_dn   = calcMTmax(iSS,MET_dn);
	//MTmax3l_dn = calcMTmax(i3l,MET_dn,true);
      }
      
      int SRSS[8]; bool selects3l[8];
      int SR3l[8];
      for(int i = 0; i<8; ++i) { SRSS[i] = -1; SR3l[i] = -1; selects3l[i] = false; }

      //SS
      //0: SR
      SRSS[0] = isSRSS(iSS,      vSS,false,MTmax,  nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      //1: SR preselect
      SRSS[1] = isSRSS(iSS,      vSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      //2: AR
      SRSS[2] = isARSS(iSS,iaSS,vaSS,false,MTmax,  nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      //3: AR preselect
      SRSS[3] = isARSS(iSS,iaSS,vaSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      //4: CR
      SRSS[4] = isCRSS(iSS,i3l,  v3l,false,MTmax3l,nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      //5: CR preselect
      SRSS[5] = isCRSS(iSS,i3l,  v3l,true ,MTmax3l,nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      //6: SR JEC up
      if(getJECunc) SRSS[6] = isSRSS(iSS, vSS,false,MTmax_up, nj30_up,nb_up,Mjj_up,MjjL_up,Detajj_up,MET_up,1);
      //7: SR JEC dn
      if(getJECunc) SRSS[7]= isSRSS(iSS, vSS,false,MTmax_dn, nj30_dn,nb_dn,Mjj_dn,MjjL_dn,Detajj_dn,MET_dn,-1);
      selects3l[4] = true; selects3l[5] = true;
      //3l
      //0: SR, 4: CR
      checkbothSRCR3l(SR3l[0],SR3l[4],i3l,false,nj,nb,MET,0,false,1);
      //1: SR preselect, 5: CR preselect
      checkbothSRCR3l(SR3l[1],SR3l[5],i3l,true ,nj,nb,MET,0,false,1);
      //2: AR
      SR3l[2] = isAR3l(i3l,ia3l,false,nj,nb,MET,0,false,1);
      //3: AR preselect
      SR3l[3] = isAR3l(i3l,ia3l,true ,nj,nb,MET,0,false,1);
      //6: SR JEC up
      if(getJECunc) SR3l[6] = isSR3l(i3l,false,nj_up,nb_up,MET_up,1);
      //7: SR JEC dn
      if(getJECunc) SR3l[7] = isSR3l(i3l,false,nj_dn,nb_dn,MET_dn,-1);
      

      for(int i = 0; i<8; ++i) {
	if(!selects3l[i]){
	  if(vetophotonprocess(fname,isphotonSS))    { SRSS[i] = -1; }
	}
	else if(vetophotonprocess(fname,isphoton3l)){ SRSS[i] = -1; }
	if(vetophotonprocess(fname,isphoton3l))     { SR3l[i] = -1; }
      }
      bool SFOS12(false), SFOS13(false), SFOS23(false); int nSFOS = 0;
      vector<LorentzVector> l;
      int mZ1(-1), mZ2(-1), mZ3(-1);
      int lZ1(-1), lZ2(-1), lZ3(-1);
      if(iSS.size()==2&&i3l.size()<3){ l.push_back(lep_p4()[iSS[0] ]); l.push_back(lep_p4()[iSS[1] ]); }
      if(i3l.size()==3){
	l.push_back(lep_p4()[i3l[0] ]); l.push_back(lep_p4()[i3l[1] ]); l.push_back(lep_p4()[i3l[2] ]);
	if(lep_pdgId()[i3l[0] ]==(-lep_pdgId()[i3l[1] ])) { SFOS12 = true; ++nSFOS; }
	if(lep_pdgId()[i3l[0] ]==(-lep_pdgId()[i3l[2] ])) { SFOS13 = true; ++nSFOS; }
	if(lep_pdgId()[i3l[1] ]==(-lep_pdgId()[i3l[2] ])) { SFOS23 = true; ++nSFOS; }
	if(nSFOS==0) { mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	if(nSFOS==1) {
	  if(SFOS12) { mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	  if(SFOS13) { mZ1 = 0; mZ2 = 2; mZ3 = 1; lZ1 = 0; lZ2 = 2; lZ3 = 1; }
	  if(SFOS23) { mZ1 = 1; mZ2 = 2; mZ3 = 0; lZ1 = 1; lZ2 = 2; lZ3 = 0; }
	}
	if(nSFOS==2) {
	  if(SFOS12&&SFOS13){
	    if(fabs((l[0]+l[1]).M()-MZ)<fabs((l[0]+l[2]).M()-MZ)){ mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 0; lZ2 = 2; lZ3 = 1; }
	    else                                                 { mZ1 = 0; mZ2 = 2; mZ3 = 1; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	  }
	  if(SFOS12&&SFOS23){
	    if(fabs((l[0]+l[1]).M()-MZ)<fabs((l[1]+l[2]).M()-MZ)){ mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 1; lZ2 = 2; lZ3 = 0; }
	    else                                                 { mZ1 = 1; mZ2 = 2; mZ3 = 0; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	  }
	  if(SFOS13&&SFOS23){
	    if(fabs((l[0]+l[2]).M()-MZ)<fabs((l[1]+l[2]).M()-MZ)){ mZ1 = 0; mZ2 = 2; mZ3 = 1; lZ1 = 1; lZ2 = 2; lZ3 = 0; }
	    else                                                 { mZ1 = 1; mZ2 = 2; mZ3 = 0; lZ1 = 0; lZ2 = 2; lZ3 = 1; }
	  }
	}
	if(nSFOS>=3) cout << "WTF nSFOS " << nSFOS << endl;
      }
      
      //0: SR, 1: SRpresel, 4: CR, 5: CRpresel
      vector<float> v;
      if(SRSS[0]>=0){
	string mysample = "";
	string mySS = "";
	if(SRSS[0]==0) { mySS = "SSee_"+sn; mysample = "SSee_"+sample; }
	if(SRSS[0]==1) { mySS = "SSem_"+sn; mysample = "SSem_"+sample; }
	if(SRSS[0]==2) { mySS = "SSmm_"+sn; mysample = "SSmm_"+sample; }
	fillhisto(histos, "NJets",       mysample, mySS, nj30,                weight);
	fillhisto(histos, "NJetsLoose",  mysample, mySS, nj,                  weight);
	fillhisto(histos, "MET",         mysample, mySS, met_pt(),            weight);
	fillhisto(histos, "Mjj",         mysample, mySS, Mjj,                 weight);
	fillhisto(histos, "MjjL",        mysample, mySS, MjjL,                weight);
	fillhisto(histos, "DetajjL",     mysample, mySS, Detajj,              weight);
	fillhisto(histos, "Mll",         mysample, mySS, (l[0]+l[1]).M(),     weight);
	v.clear();
	v.push_back(mT(l[0],MET));
	v.push_back(mT(l[1],MET));
	sort(v.begin(),v.end(),sortDecreasing);
	fillhisto(histos, "MTmax",       mysample, mySS, v[0],                weight);
	fillhisto(histos, "MTmin",       mysample, mySS, v[1],                weight);
	fillhisto(histos, "MTsum",       mysample, mySS, v[0]+v[1],           weight);
	fillhisto(histos, "DPhill",      mysample, mySS, dPhi(l[0],l[1]),     weight);
	fillhisto(histos, "DPhillMET",   mysample, mySS, dPhi(l[0]+l[1],MET), weight);
	v.clear();
	v.push_back(dPhi(l[0],MET));
	v.push_back(dPhi(l[1],MET));
	sort(v.begin(),v.end(),sortDecreasing);
	fillhisto(histos, "DPhilMETmax", mysample, mySS, v[0],                weight);
	fillhisto(histos, "DPhilMETmin", mysample, mySS, v[1],                weight);
	fillhisto(histos, "pTll",        mysample, mySS, (l[0]+l[1]).Pt(),    weight);
	//new 2017/11/22
	vector<LorentzVector> myjets; myjets.clear();
	for(unsigned int n = 0; n<jets_csv().size();++n){
	  if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5) myjets.push_back(jets_p4()[n]);
	}
	LorentzVector j1 = LorentzVector(0,0,0,0); LorentzVector j2 = j1;
	float minDR=999.;
	float ST = 0;
	for(unsigned int n1 = 0; n1<myjets.size();++n1){
	  ST += myjets[n1].Pt();
	  for(unsigned int n2 = n1+1; n2<myjets.size();++n2){
	    if(dR(myjets[n1], myjets[n2])<minDR){
	      minDR = dR(myjets[n1], myjets[n2]);
	      j1 = myjets[n1]; j2 = myjets[n2];
	    }
	  }
	}
	fillhisto(histos, "Mjjll",         mysample, mySS, (l[0]+l[1]+j1+j2).M(),           weight);
	fillhisto(histos, "MTjjllMET",     mysample, mySS, mT(l[0]+l[1]+j1+j2,MET),         weight);
	minDR=999.;
	float MljDR=-1;
	for(unsigned int n = 0; n<myjets.size();++n){
	  if(dR(myjets[n], l[0])<minDR){
	    minDR=dR(myjets[n],  l[0]);
	    MljDR = (myjets[n] + l[0]).M();
	  }
	  if(dR(myjets[n], l[1])<minDR){
	    minDR=dR(myjets[n],  l[1]);
	    MljDR = (myjets[n] + l[1]).M();
	  }
	}
	fillhisto(histos, "MljDR",         mysample, mySS, MljDR,                           weight);
	v.clear();
	for(unsigned int n = 0; n<myjets.size();++n){
	  v.push_back(dPhi(l[0],myjets[n]));
	  v.push_back(dPhi(l[1],myjets[n]));
	}
	sort(v.begin(),v.end());//increasing
	fillhisto(histos, "DPhiljmin",     mysample, mySS, v[0],                            weight);
	v.clear();
	for(unsigned int n = 0; n<myjets.size();++n){
	  v.push_back(dR(l[0],myjets[n]));
	  v.push_back(dR(l[1],myjets[n]));
	}
	sort(v.begin(),v.end());//increasing
	fillhisto(histos, "DRljmin",       mysample, mySS, v[0],                            weight);
	fillhisto(histos, "RatMETvsll",    mysample, mySS, MET.Pt()/(l[0]+l[1]).Pt(),       weight);
	fillhisto(histos, "DPhijjvsllMET", mysample, mySS, dPhi(j1+j2,l[0]+l[1]+MET),       weight);
	fillhisto(histos, "pTlljjMET",     mysample, mySS, (j1+j2+l[0]+l[1]+MET).Pt(),      weight);
	fillhisto(histos, "ST",            mysample, mySS, ST+l[0].Pt()+l[1].Pt()+met_pt(), weight);
      }
      
      if(SR3l[0]==1||SR3l[4]==1){
	fillhisto(histos, "MSFOS_1SFOS",      sample, sn2, (l[mZ1]+l[mZ2]).M(), weight);
      }
      if(SR3l[0]==2||SR3l[4]==2){
	fillhisto(histos, "MSFOSZlike_2SFOS", sample, sn2, (l[mZ1]+l[mZ2]).M(), weight);
	fillhisto(histos, "MSFOSnoZ_2SFOS",   sample, sn2, (l[lZ1]+l[lZ2]).M(), weight);
      }
      if(SR3l[0]>=0){
	string mysample = "";
	string my3l = "";
	if(SR3l[0]==0) { my3l = "0SFOS_"+sn2; mysample = "0SFOS_"+sample; }
	if(SR3l[0]==1) { my3l = "1SFOS_"+sn2; mysample = "1SFOS_"+sample; }
	if(SR3l[0]==2) { my3l = "2SFOS_"+sn2; mysample = "2SFOS_"+sample; }
	v.clear();
	v.push_back(mT(l[0],MET));
	v.push_back(mT(l[1],MET));
	v.push_back(mT(l[2],MET));
	sort(v.begin(),v.end(),sortDecreasing);
	fillhisto(  histos, "MTmax",               mysample, my3l, v[0],                           weight);
	fillhisto(  histos, "MTmin",               mysample, my3l, v[2],                           weight);
	fillhisto(  histos, "MTsum",               mysample, my3l, v[0]+v[1]+v[2],                 weight);
	fillhisto(  histos, "MTlll",               mysample, my3l, mT(l[0]+l[1]+l[2],MET),         weight);
	fillhisto(  histos, "pTlllMET",            mysample, my3l, (l[0]+l[1]+l[2]+MET).Pt(),      weight);
	fillhisto(  histos, "MET",                 mysample, my3l, met_pt(),                       weight);
	fillhisto(  histos, "Mlll",                mysample, my3l, (l[0]+l[1]+l[2]).M(),           weight);
	fillhisto(  histos, "pTlll",               mysample, my3l, (l[0]+l[1]+l[2]).Pt(),          weight);
	fillhisto(  histos, "DPhilllMET",          mysample, my3l, dPhi(l[0]+l[1]+l[2],MET),       weight);	
	if(SR3l[0]==1) {
	  fillhisto(histos, "MT3rd",               mysample, my3l, mT(l[mZ3],MET),                 weight);
	  fillhisto(histos, "pTSFOS",              mysample, my3l, (l[mZ1]+l[mZ2]).Pt(),           weight);
	  fillhisto(histos, "pT3rd",               mysample, my3l, (l[mZ3]).Pt(),                  weight);
	  fillhisto(histos, "pTW",                 mysample, my3l, (l[mZ3]+MET).Pt(),              weight);
	  fillhisto(histos, "DPhillofZ",           mysample, my3l, dPhi(l[mZ1],l[mZ2]),            weight);
	  fillhisto(histos, "DPhillofZvsl",        mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]),     weight);
	  fillhisto(histos, "DPhillofZvslMET",     mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]+MET), weight);
	}
	if(SR3l[0]==2) {
	  v.clear();
	  v.push_back((l[mZ1]+l[mZ2]).Pt());
	  v.push_back((l[lZ1]+l[lZ2]).Pt());
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "MTleastZ",            mysample, my3l, mT(l[mZ3],MET),                 weight);
	  fillhisto(histos, "pTSFOSmax",           mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "pTSFOSmin",           mysample, my3l, v[1],                           weight);
	  fillhisto(histos, "pTSFOSZlike",         mysample, my3l, (l[mZ1]+l[mZ2]).Pt(),           weight);
	  fillhisto(histos, "pTSFOSnoZ",           mysample, my3l, (l[lZ1]+l[lZ2]).Pt(),           weight);
	  fillhisto(histos, "pTleastZ",            mysample, my3l, (l[mZ3]).Pt(),                  weight);
	  fillhisto(histos, "pTWleastZ",           mysample, my3l, (l[mZ3]+MET).Pt(),              weight);
	  v.clear();
	  v.push_back(dPhi(l[mZ1],l[mZ2]));
	  v.push_back(dPhi(l[lZ1],l[lZ2]));
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "DPhillofZlike",       mysample, my3l, dPhi(l[mZ1],l[mZ2]),            weight);
	  fillhisto(histos, "DPhilloflessZ",       mysample, my3l, dPhi(l[lZ1],l[lZ2]),            weight);
	  fillhisto(histos, "DPhillSFOSmax",       mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "DPhillSFOSmin",       mysample, my3l, v[1],                           weight);
	  v.clear();
	  v.push_back(dPhi(l[mZ1]+l[mZ2],l[mZ3]));
	  v.push_back(dPhi(l[lZ1]+l[lZ2],l[lZ3]));
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "DPhillofZlikevsl",    mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]),     weight);
	  fillhisto(histos, "DPhillSFOSvslmax",    mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "DPhillSFOSvslmin",    mysample, my3l, v[1],                           weight);
	  v.clear();
	  v.push_back(dPhi(l[mZ1]+l[mZ2],l[mZ3]+MET));
	  v.push_back(dPhi(l[lZ1]+l[lZ2],l[lZ3]+MET));
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "DPhillofZlikevslMET", mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]+MET), weight);
	  fillhisto(histos, "DPhillSFOSvslMETmax", mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "DPhillSFOSvslMETmin", mysample, my3l, v[1],                           weight);
	}
	float ST = 0;
	for(unsigned int n = 0; n<jets_csv().size();++n){
	  if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5.0) ST += jets_p4()[n].Pt();
	}
	fillhisto(histos, "DRllsum",     mysample, my3l, dR(l[0],l[1])+dR(l[0],l[2])+dR(l[1],l[2]),          weight);
	fillhisto(histos, "pTllsum",     mysample, my3l, (l[0]+l[1]).Pt()+(l[0]+l[2]).Pt()+(l[1]+l[2]).Pt(), weight);
	fillhisto(histos, "RatMETvslll", mysample, my3l, MET.Pt()/(l[0]+l[1]+l[2]).Pt(),                     weight);
	fillhisto(histos, "ST",          mysample, my3l, ST+MET.Pt()+l[0].Pt()+l[1].Pt()+l[2].Pt(),          weight);

      }


      //PRESELECT
      if(SRSS[1]>=0){
	string mysample = "";
	string mySS = "";
	if(SRSS[1]==0) { mySS = "SSee_"+sn; mysample = "SSee_"+sample; }
	if(SRSS[1]==1) { mySS = "SSem_"+sn; mysample = "SSem_"+sample; }
	if(SRSS[1]==2) { mySS = "SSmm_"+sn; mysample = "SSmm_"+sample; }
	fillhisto(histos, "ZPreselect_NJets",       mysample, mySS, nj30,                weight);
	fillhisto(histos, "ZPreselect_NJetsLoose",  mysample, mySS, nj,                  weight);	
	fillhisto(histos, "ZPreselect_MET",         mysample, mySS, met_pt(),            weight);
	fillhisto(histos, "ZPreselect_Mjj",         mysample, mySS, Mjj,                 weight);
	fillhisto(histos, "ZPreselect_MjjL",        mysample, mySS, MjjL,                weight);
	fillhisto(histos, "ZPreselect_DetajjL",     mysample, mySS, Detajj,              weight);
	fillhisto(histos, "ZPreselect_Mll",         mysample, mySS, (l[0]+l[1]).M(),     weight);
	v.clear();
	v.push_back(mT(l[0],MET));
	v.push_back(mT(l[1],MET));
	sort(v.begin(),v.end(),sortDecreasing);
	fillhisto(histos, "ZPreselect_MTmax",       mysample, mySS, v[0],                weight);
	fillhisto(histos, "ZPreselect_MTmin",       mysample, mySS, v[1],                weight);
	fillhisto(histos, "ZPreselect_MTsum",       mysample, mySS, v[0]+v[1],           weight);
	fillhisto(histos, "ZPreselect_DPhill",      mysample, mySS, dPhi(l[0],l[1]),     weight);
	fillhisto(histos, "ZPreselect_DPhillMET",   mysample, mySS, dPhi(l[0]+l[1],MET), weight);
	v.clear();
	v.push_back(dPhi(l[0],MET));
	v.push_back(dPhi(l[1],MET));
	sort(v.begin(),v.end(),sortDecreasing);
	fillhisto(histos, "ZPreselect_DPhilMETmax", mysample, mySS, v[0],                weight);
	fillhisto(histos, "ZPreselect_DPhilMETmin", mysample, mySS, v[1],                weight);
	fillhisto(histos, "ZPreselect_pTll",        mysample, mySS, (l[0]+l[1]).Pt(),    weight);
	//new 2017/11/22
	vector<LorentzVector> myjets; myjets.clear();
	for(unsigned int n = 0; n<jets_csv().size();++n){
	  if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5) myjets.push_back(jets_p4()[n]);
	}
	LorentzVector j1 = LorentzVector(0,0,0,0); LorentzVector j2 = j1;
	float minDR=999.;
	float ST = 0;
	for(unsigned int n1 = 0; n1<myjets.size();++n1){
	  ST += myjets[n1].Pt();
	  for(unsigned int n2 = n1+1; n2<myjets.size();++n2){
	    if(dR(myjets[n1], myjets[n2])<minDR){
	      minDR = dR(myjets[n1], myjets[n2]);
	      j1 = myjets[n1]; j2 = myjets[n2];
	    }
	  }
	}
	fillhisto(histos, "ZPreselect_Mjjll",         mysample, mySS, (l[0]+l[1]+j1+j2).M(),           weight);
	fillhisto(histos, "ZPreselect_MTjjllMET",     mysample, mySS, mT(l[0]+l[1]+j1+j2,MET),         weight);
	minDR=999.;
	float MljDR=-1;
	for(unsigned int n = 0; n<myjets.size();++n){
	  if(dR(myjets[n], l[0])<minDR){
	    minDR=dR(myjets[n],  l[0]);
	    MljDR = (myjets[n] + l[0]).M();
	  }
	  if(dR(myjets[n], l[1])<minDR){
	    minDR=dR(myjets[n],  l[1]);
	    MljDR = (myjets[n] + l[1]).M();
	  }
	}
	fillhisto(histos, "ZPreselect_MljDR",         mysample, mySS, MljDR,                           weight);
	v.clear();
	for(unsigned int n = 0; n<myjets.size();++n){
	  v.push_back(dPhi(l[0],myjets[n]));
	  v.push_back(dPhi(l[1],myjets[n]));
	}
	sort(v.begin(),v.end());//increasing
	fillhisto(histos, "ZPreselect_DPhiljmin",     mysample, mySS, v[0],                            weight);
	v.clear();
	for(unsigned int n = 0; n<myjets.size();++n){
	  v.push_back(dR(l[0],myjets[n]));
	  v.push_back(dR(l[1],myjets[n]));
	}
	sort(v.begin(),v.end());//increasing
	fillhisto(histos, "ZPreselect_DRljmin",       mysample, mySS, v[0],                            weight);
	fillhisto(histos, "ZPreselect_RatMETvsll",    mysample, mySS, MET.Pt()/(l[0]+l[1]).Pt(),       weight);
	fillhisto(histos, "ZPreselect_DPhijjvsllMET", mysample, mySS, dPhi(j1+j2,l[0]+l[1]+MET),       weight);
	fillhisto(histos, "ZPreselect_pTlljjMET",     mysample, mySS, (j1+j2+l[0]+l[1]+MET).Pt(),      weight);
	fillhisto(histos, "ZPreselect_ST",            mysample, mySS, ST+l[0].Pt()+l[1].Pt()+met_pt(), weight);
      }
      if(SR3l[1]==1||SR3l[5]==1){
	fillhisto(histos, "ZPreselect_MSFOS_1SFOS",      sample, sn2, (l[mZ1]+l[mZ2]).M(), weight);
      }
      if(SR3l[1]==2||SR3l[5]==2){
	fillhisto(histos, "ZPreselect_MSFOSZlike_2SFOS", sample, sn2, (l[mZ1]+l[mZ2]).M(), weight);
	fillhisto(histos, "ZPreselect_MSFOSnoZ_2SFOS",   sample, sn2, (l[lZ1]+l[lZ2]).M(), weight);
      }
      if(SR3l[1]>=0){
	string mysample = "";
	string my3l = "";
	if(SR3l[1]==0) { my3l = "0SFOS_"+sn2; mysample = "0SFOS_"+sample; }
	if(SR3l[1]==1) { my3l = "1SFOS_"+sn2; mysample = "1SFOS_"+sample; }
	if(SR3l[1]==2) { my3l = "2SFOS_"+sn2; mysample = "2SFOS_"+sample; }
	v.clear();
	v.push_back(mT(l[0],MET));
	v.push_back(mT(l[1],MET));
	v.push_back(mT(l[2],MET));
	sort(v.begin(),v.end(),sortDecreasing);
	fillhisto(  histos, "ZPreselect_MTmax",               mysample, my3l, v[0],                           weight);
	fillhisto(  histos, "ZPreselect_MTmin",               mysample, my3l, v[2],                           weight);
	fillhisto(  histos, "ZPreselect_MTsum",               mysample, my3l, v[0]+v[1]+v[2],                 weight);
	fillhisto(  histos, "ZPreselect_MTlll",               mysample, my3l, mT(l[0]+l[1]+l[2],MET),         weight);
	fillhisto(  histos, "ZPreselect_pTlllMET",            mysample, my3l, (l[0]+l[1]+l[2]+MET).Pt(),      weight);
	fillhisto(  histos, "ZPreselect_MET",                 mysample, my3l, met_pt(),                       weight);
	fillhisto(  histos, "ZPreselect_Mlll",                mysample, my3l, (l[0]+l[1]+l[2]).M(),           weight);
	fillhisto(  histos, "ZPreselect_pTlll",               mysample, my3l, (l[0]+l[1]+l[2]).Pt(),          weight);
	fillhisto(  histos, "ZPreselect_DPhilllMET",          mysample, my3l, dPhi(l[0]+l[1]+l[2],MET),       weight);
	if(SR3l[1]==1) {
	  fillhisto(histos, "ZPreselect_MT3rd",               mysample, my3l, mT(l[mZ3],MET),                 weight);
	  fillhisto(histos, "ZPreselect_pTSFOS",              mysample, my3l, (l[mZ1]+l[mZ2]).Pt(),           weight);
	  fillhisto(histos, "ZPreselect_pT3rd",               mysample, my3l, (l[mZ3]).Pt(),                  weight);
	  fillhisto(histos, "ZPreselect_pTW",                 mysample, my3l, (l[mZ3]+MET).Pt(),              weight);
	  fillhisto(histos, "ZPreselect_DPhillofZ",           mysample, my3l, dPhi(l[mZ1],l[mZ2]),            weight);
	  fillhisto(histos, "ZPreselect_DPhillofZvsl",        mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]),     weight);
	  fillhisto(histos, "ZPreselect_DPhillofZvslMET",     mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]+MET), weight);
	}
	if(SR3l[1]==2) {
	  v.clear();
	  v.push_back((l[mZ1]+l[mZ2]).Pt());
	  v.push_back((l[lZ1]+l[lZ2]).Pt());
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "ZPreselect_MTleastZ",            mysample, my3l, mT(l[mZ3],MET),                 weight);
	  fillhisto(histos, "ZPreselect_pTSFOSmax",           mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "ZPreselect_pTSFOSmin",           mysample, my3l, v[1],                           weight);
	  fillhisto(histos, "ZPreselect_pTSFOSZlike",         mysample, my3l, (l[mZ1]+l[mZ2]).Pt(),           weight);
	  fillhisto(histos, "ZPreselect_pTSFOSnoZ",           mysample, my3l, (l[lZ1]+l[lZ2]).Pt(),           weight);
	  fillhisto(histos, "ZPreselect_pTleastZ",            mysample, my3l, (l[mZ3]).Pt(),                  weight);
	  fillhisto(histos, "ZPreselect_pTWleastZ",           mysample, my3l, (l[mZ3]+MET).Pt(),              weight);
	  v.clear();
	  v.push_back(dPhi(l[mZ1],l[mZ2]));
	  v.push_back(dPhi(l[lZ1],l[lZ2]));
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "ZPreselect_DPhillofZlike",       mysample, my3l, dPhi(l[mZ1],l[mZ2]),            weight);
	  fillhisto(histos, "ZPreselect_DPhilloflessZ",       mysample, my3l, dPhi(l[lZ1],l[lZ2]),            weight);
	  fillhisto(histos, "ZPreselect_DPhillSFOSmax",       mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "ZPreselect_DPhillSFOSmin",       mysample, my3l, v[1],                           weight);
	  v.clear();
	  v.push_back(dPhi(l[mZ1]+l[mZ2],l[mZ3]));
	  v.push_back(dPhi(l[lZ1]+l[lZ2],l[lZ3]));
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "ZPreselect_DPhillofZlikevsl",    mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]),     weight);
	  fillhisto(histos, "ZPreselect_DPhillSFOSvslmax",    mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "ZPreselect_DPhillSFOSvslmin",    mysample, my3l, v[1],                           weight);
	  v.clear();
	  v.push_back(dPhi(l[mZ1]+l[mZ2],l[mZ3]+MET));
	  v.push_back(dPhi(l[lZ1]+l[lZ2],l[lZ3]+MET));
	  sort(v.begin(),v.end(),sortDecreasing);
	  fillhisto(histos, "ZPreselect_DPhillofZlikevslMET", mysample, my3l, dPhi(l[mZ1]+l[mZ2],l[mZ3]+MET), weight);
	  fillhisto(histos, "ZPreselect_DPhillSFOSvslMETmax", mysample, my3l, v[0],                           weight);
	  fillhisto(histos, "ZPreselect_DPhillSFOSvslMETmin", mysample, my3l, v[1],                           weight);
	}
	float ST = 0;
	for(unsigned int n = 0; n<jets_csv().size();++n){
	  if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5.0) ST += jets_p4()[n].Pt();
	}
	fillhisto(histos, "ZPreselect_DRllsum",     mysample, my3l, dR(l[0],l[1])+dR(l[0],l[2])+dR(l[1],l[2]),          weight);
	fillhisto(histos, "ZPreselect_pTllsum",     mysample, my3l, (l[0]+l[1]).Pt()+(l[0]+l[2]).Pt()+(l[1]+l[2]).Pt(), weight);
	fillhisto(histos, "ZPreselect_RatMETvslll", mysample, my3l, MET.Pt()/(l[0]+l[1]+l[2]).Pt(),                     weight);
	fillhisto(histos, "ZPreselect_ST",          mysample, my3l, ST+MET.Pt()+l[0].Pt()+l[1].Pt()+l[2].Pt(),          weight);

      }

      
    }//event loop
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  SaveHistosToFile("rootfiles/ImproveSelection_newBabe.root",histos,true,true,(chainnumber==0));
  deleteHistograms(histos);

  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
