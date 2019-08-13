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
#define USE_CMS3_WWW100 

#include "Functions.h"
#ifdef USE_CMS3_WWW100
#include "CMS3_WWW106.cc"
#else
#include "CMS3_WWW0118.cc"
#endif
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1) {

  bool blindSR         = true;
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

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int>    hbins;      hbins.clear();
  vector<float>  hlow;       hlow.clear();
  vector<float>  hup;        hup.clear();

  histonames.push_back("SignalRegion");                        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegionPresel");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  //plot variables
  histonames.push_back("MET_SSee");                             hbins.push_back(15); hlow.push_back(    0); hup.push_back(150);
  histonames.push_back("MET_SSem");                             hbins.push_back(15); hlow.push_back(    0); hup.push_back(150);
  histonames.push_back("MET_SSmm");                             hbins.push_back(15); hlow.push_back(    0); hup.push_back(150);
  histonames.push_back("Mll_SSee");                             hbins.push_back(15); hlow.push_back(    0); hup.push_back(150);
  histonames.push_back("Mll_SSem");                             hbins.push_back(15); hlow.push_back(    0); hup.push_back(150);
  histonames.push_back("Mll_SSmm");                             hbins.push_back(15); hlow.push_back(    0); hup.push_back(150);
  histonames.push_back("Mjj_SSee");                             hbins.push_back(16); hlow.push_back(   40); hup.push_back(120);
  histonames.push_back("Mjj_SSem");                             hbins.push_back(16); hlow.push_back(   40); hup.push_back(120);
  histonames.push_back("Mjj_SSmm");                             hbins.push_back(16); hlow.push_back(   40); hup.push_back(120);
  histonames.push_back("Detajj_SSee");                          hbins.push_back(15); hlow.push_back(    0); hup.push_back(3.0);
  histonames.push_back("Detajj_SSem");                          hbins.push_back(15); hlow.push_back(    0); hup.push_back(3.0);
  histonames.push_back("Detajj_SSmm");                          hbins.push_back(15); hlow.push_back(    0); hup.push_back(3.0);
  
  histonames.push_back("NJets_SSee");                           hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("NJets_SSem");                           hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("NJets_SSmm");                           hbins.push_back( 4); hlow.push_back(    2); hup.push_back(6);
  histonames.push_back("MTmax_SSee");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_SSee");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTsum_SSee");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTmax_SSem");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_SSem");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTsum_SSem");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("MTmax_SSmm");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTmin_SSmm");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(250);
  histonames.push_back("MTsum_SSmm");                           hbins.push_back(10); hlow.push_back(    0); hup.push_back(500);
  histonames.push_back("DPhill_SSee");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhill_SSem");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhill_SSmm");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(3.2);
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
  histonames.push_back("DPhiljmin_SSee");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhiljmin_SSem");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhiljmin_SSmm");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DRljmin_SSee");                         hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DRljmin_SSem");                         hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DRljmin_SSmm");                         hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("RatMETvsll_SSee");                      hbins.push_back(16); hlow.push_back(    0); hup.push_back(4);
  histonames.push_back("RatMETvsll_SSem");                      hbins.push_back(16); hlow.push_back(    0); hup.push_back(4);
  histonames.push_back("RatMETvsll_SSmm");                      hbins.push_back(16); hlow.push_back(    0); hup.push_back(4);

  histonames.push_back("MET_0SFOS");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MET_1SFOS");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MET_2SFOS");                            hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("DPhilllMET_0SFOS");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilllMET_1SFOS");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhilllMET_2SFOS");                     hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("pTlll_0SFOS");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("pTlll_1SFOS");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("pTlll_2SFOS");                          hbins.push_back(12); hlow.push_back(    0); hup.push_back(120);
  histonames.push_back("MT3rd_1SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MT3rd_2SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTmin_0SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTmin_1SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTmin_2SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTmax_0SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTmax_1SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTmax_2SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTsum_0SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(480);
  histonames.push_back("MTsum_1SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(480);
  histonames.push_back("MTsum_2SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(480);
  histonames.push_back("MTlll_0SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTlll_1SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MTlll_2SFOS");                          hbins.push_back(24); hlow.push_back(    0); hup.push_back(240);
  histonames.push_back("MSFOSZ_1SFOS");                         hbins.push_back(18); hlow.push_back(   50); hup.push_back(140);
  histonames.push_back("MSFOSZ_2SFOS");                         hbins.push_back(18); hlow.push_back(   50); hup.push_back(140);
  histonames.push_back("pTSFOS_1SFOS");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTSFOS_2SFOS");                         hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pT3rd_1SFOS");                          hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pT3rd_2SFOS");                          hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTW_1SFOS");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTW_2SFOS");                            hbins.push_back(10); hlow.push_back(    0); hup.push_back(200);
  histonames.push_back("pTlllMET_0SFOS");                       hbins.push_back(14); hlow.push_back(    0); hup.push_back(280);
  histonames.push_back("pTlllMET_1SFOS");                       hbins.push_back(14); hlow.push_back(    0); hup.push_back(280);
  histonames.push_back("pTlllMET_2SFOS");                       hbins.push_back(14); hlow.push_back(    0); hup.push_back(280);
  histonames.push_back("DPhiZvsl_1SFOS");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhiZvsl_2SFOS");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhiZvsW_1SFOS");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("DPhiZvsW_2SFOS");                       hbins.push_back(16); hlow.push_back(    0); hup.push_back(3.2);
  histonames.push_back("RatMETvslll_0SFOS");                    hbins.push_back(10); hlow.push_back(    0); hup.push_back(2.5);
  histonames.push_back("RatMETvslll_1SFOS");                    hbins.push_back(10); hlow.push_back(    0); hup.push_back(2.5);
  histonames.push_back("RatMETvslll_2SFOS");                    hbins.push_back(10); hlow.push_back(    0); hup.push_back(2.5);          
  unsigned int nhistotemp = histonames.size();
  for(unsigned int i = 0; i < nhistotemp; ++i){ histonames.push_back("ZPresel_"+histonames[i]); hbins.push_back(hbins[i]); hlow.push_back(hlow[i]); hup.push_back(hup[i]); }

  //SR binning
  histonames.push_back("SR_addone_SSee");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMET_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMETMjj_SSee");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMETMll_SSee");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMETMjjMll_SSee");   hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMjj_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMll_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDeta_SSee");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDetaMjj_SSee");     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_SSem");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMET_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMETMjj_SSem");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMETMll_SSem");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMjj_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMll_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDeta_SSem");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDetaMjj_SSem");     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMTmax_SSem");       hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMTmaxMll_SSem");    hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMTmaxMjj_SSem");    hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_SSmm");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMET_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMETMjj_SSmm");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMETMll_SSmm");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMjj_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMll_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDeta_SSmm");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDetaMjj_SSmm");     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_0SFOS");                     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMET_0SFOS");             hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDPhilllMET_0SFOS");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_1SFOS");                     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMET_1SFOS");             hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_droppTlll_1SFOS");           hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDPhilllMET_1SFOS");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_droppTlllDPhilllMET_1SFOS"); hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMT3rd_1SFOS");           hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_2SFOS");                     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropMET_2SFOS");             hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_droppTlll_2SFOS");           hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_dropDPhilllMET_2SFOS");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SR_addone_droppTlllDPhilllMET_2SFOS"); hbins.push_back(75); hlow.push_back(0); hup.push_back(75);

  //SR binning: one jet
  histonames.push_back("SRonejet_addone_SSee");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMET_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMETMll_SSee");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMll_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_SSem");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMET_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMETMll_SSem");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMll_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMTmax_SSem");       hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMTmaxMll_SSem");    hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_SSmm");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMET_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMETMll_SSmm");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRonejet_addone_dropMll_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);

  //SR binning - soft jet
  histonames.push_back("SRsoftjet_addone_SSee");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMET_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMETMjj_SSee");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMETMll_SSee");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMETMjjMll_SSee");   hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMjj_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMll_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDeta_SSee");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDetaMjj_SSee");     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_SSem");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMET_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMETMjj_SSem");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMETMll_SSem");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMjj_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMll_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDeta_SSem");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDetaMjj_SSem");     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMTmax_SSem");       hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMTmaxMll_SSem");    hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMTmaxMjj_SSem");    hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_SSmm");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMET_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMETMjj_SSmm");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMETMll_SSmm");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMjj_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMll_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDeta_SSmm");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDetaMjj_SSmm");     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_0SFOS");                     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMET_0SFOS");             hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDPhilllMET_0SFOS");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_1SFOS");                     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMET_1SFOS");             hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_droppTlll_1SFOS");           hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDPhilllMET_1SFOS");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_droppTlllDPhilllMET_1SFOS"); hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMT3rd_1SFOS");           hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_2SFOS");                     hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropMET_2SFOS");             hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_droppTlll_2SFOS");           hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_dropDPhilllMET_2SFOS");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsoftjet_addone_droppTlllDPhilllMET_2SFOS"); hbins.push_back(75); hlow.push_back(0); hup.push_back(75);

  //SR binning
  histonames.push_back("SRsideMjj_addone_SSee");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMET_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMETMll_SSee");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMll_SSee");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropDeta_SSee");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_SSem");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMET_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMETMll_SSem");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMll_SSem");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropDeta_SSem");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMTmax_SSem");       hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMTmaxMll_SSem");    hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_SSmm");                 hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMET_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMETMll_SSmm");      hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropMll_SSmm");         hbins.push_back(75); hlow.push_back(0); hup.push_back(75);
  histonames.push_back("SRsideMjj_addone_dropDeta_SSmm");        hbins.push_back(75); hlow.push_back(0); hup.push_back(75);

  
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
      if(nLlep()<2)              continue;
      if(nTlepSS()<1&&nTlep()<2) continue;//preselection can be done already here
      if(nb()!=0)                continue;//preselection can be done already here
      //VETO DATA
      if(isData()) continue;

      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_") !=string::npos && gen_ht()>100.) continue;
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos && gen_ht()>100.) continue;
      //if(string(currentFile->GetTitle()).find("www_2l_mia")        !=string::npos) weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")   !=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(isData()) weight = 1.;
      //double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      if(!isData()&&applyPUrewgt)    weight *= purewgt();
      if(!isData()&&applylepSF)      weight *= lepsf();
      if(!isData()&&applytrigSF)     weight *= trigeff();
             
      if(isData()){
        if(!passFilters())                      continue;
        duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
        if( is_duplicate(id) )                  continue; 
        if( !goodrun(tas::run(), tas::lumi()) ) continue;
      } 
      if(!passTriggers()) continue;//pass trigger for data, and offline lepton kinematic cuts for data/simulation
      
      string sample   = skimFilePrefix;
      if(splitVH(fname)){ sample = "WHtoWWW"; }
      string sn = string(bkgtype().Data());
      if(vetophoton()) continue;
      
      int SRSS[20]; 
      int SR3l[20];
      for(int i = 0; i<20; ++i) { SRSS[i] = -1; SR3l[i] = -1;  }

      //SS
      SRSS[0] = isSRSS(false);//0: SR
      SRSS[1] = isSRSS(true );//1: SR presel
      SR3l[0] = isSR3l(false);//0: SR
      SR3l[1] = isSR3l(true );//1: SR presel

      int lepSR = -1;
      string SRlep = "";
      if(NtightSS()==2&&nVlep()==2&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0){
        if(passSSee() && fabs(MeeSS()-MZ)>10.)                                     { lepSR = 1; SRlep = "_SSee"; }
        if(passSSem())                                                             { lepSR = 2; SRlep = "_SSem"; }
        if(passSSmm())                                                             { lepSR = 3; SRlep = "_SSmm"; }
      }
      if(Ntight3l()==3&&nVlep()==3&&abs(lep_charge()[0]+lep_charge()[1]+lep_charge()[2])==1&&Mll3L ()>=20.&&fabs(M3l()-MZ)>=10.){
        if(nSFOS()==0&&fabs(Mee3L()-MZ)>=15.)                                        { lepSR = 4; SRlep = "_0SFOS"; }
        if(nSFOS()==1&&fabs(Mll3L()-MZ)>=20.)                                        { lepSR = 5; SRlep = "_1SFOS"; }
        if(nSFOS()==2&&fabs(Mll3L()-MZ)>=20.&&fabs(Mll3L1()-MZ)>=20.&&Mll3L1()>=20.) { lepSR = 6; SRlep = "_2SFOS"; }
      }
      if(lepSR<0) continue;


      if(!isData()||!blindSR){//SR is blinded
        //cout << __LINE__ << endl;
        fillSRhisto(histos, "SignalRegion",               sample, sn, SRSS[0], SR3l[0], weight);
        fillSRhisto(histos, "SignalRegionPresel",         sample, sn, SRSS[1], SR3l[1], weight);
        //cout << __LINE__ << endl;
      }

      bool passNJ    = false;
      bool passMjj   = false;
      bool passMjjL  = false;
      bool passDeta  = false;
      bool passMll   = false;
      bool passMET   = false;
      bool passMT    = false;
      bool passDPhi  = false;
      bool passPTlll = false;

      bool passMjjside = false;
      bool passOneJet  = false;
      bool passSoftJet = false;
      int njsoft = numJ(20.,3.0,-1,0) - nj30();
      if(njsoft<0) cout << "ERROR " << __LINE__ << " " << numJ(20.,3.0,-1,0) << " " << nj30() << endl;
      int j1(-1), j2(-1);
      if(lepSR>=1&&lepSR<=3){
        passNJ      = (nj30()>=2);
        passMjj     = (fabs(Mjj()-80.)<15.);
        passMjjL    = (fabs(MjjL())<=400.);
        passDeta    = (fabs(DetajjL())<=1.5);
        passMll     = (lepSR==2 ? MllSS()>30. : MllSS()>40.);
        passMET     = (lepSR==3 ? true : met_pt() > 60.);
        passMT      = (lepSR==2 ? MTmax()>90. : true);
        passDPhi    = true;
        passPTlll   = true;
        passMjjside = !passMjj;
        passOneJet  = (nj30()==1);
        passSoftJet = (nj30()==1&&njsoft>=1);
        if(passSoftJet){
          for(unsigned int n = 0; n<jets_csv().size();++n){
            if(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5) j1 = n;
            break;
          }
          for(unsigned int n = 0; n<jets_csv().size();++n){
            if(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5) continue;
            if(jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=3.0) {
              if(j2<0) j2 = n;
              else if(dR(jets_p4()[n],jets_p4()[j1])<dR(jets_p4()[j2],jets_p4()[j1])) j2 = n;
            }
          }
        }
      }
      else if(lepSR>=4){
        passNJ               = (nj()<=1);
        passMjj              = false;
        passMjjL             = false;
        passDeta             = false;
        passMll              = false;
        if(lepSR==4) passMET = (met_pt()>=30.);
        if(lepSR==5) passMET = (met_pt()>=40.);
        if(lepSR==6) passMET = (met_pt()>=55.);
        passMT               = (lepSR==5 ? MT3rd()>90. : true);
        passDPhi             = DPhi3lMET()>=2.5;
        passPTlll            = (lepSR==4 ? true : Pt3l()>=60.);
        passMjjside          = false;
        passOneJet           = false;
        passSoftJet          = !passNJ;
      }
      LorentzVector METlv; METlv.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      float MET         = met_pt();
      float Mll         = MllSS();
      float MassJJ      = Mjj();
      float Detajj      = DetajjL();
      float NJets       = nj30();
      float maxMT       = MTmax();
      float minMT       = MTmin();
      float MTsum       = MTmax()+MTmin();
      float DPhill      = -1;
      float DPhillMET   = -1;
      float DPhilMETmin = -1;
      float DPhilMETmax = -1;
      float DPhiljmin   = -1;
      float DRljmin     = -1;
      float RatMETvsll  = -1;
      float pTll        = -1;
      float DPhilllMET  = DPhi3lMET();
      float pTlll       = Pt3l();
      float thirdMT       = MT3rd();
      float MTlll       = -1;
      float MSFOS       = -1;
      float pTSFOS      = -1;
      float pT3rd       = -1;
      float pTW         = -1;
      float pTlllMET    = -1;
      float DPhiZvsl    = -1;
      float DPhiZvsW    = -1;
      float RatMETvslll = met_pt()/Pt3l();
      int id1(-1), id2(-1), id3(-1);
      for(unsigned int n = 0; n<lep_pass_VVV_cutbased_tight().size(); ++n){
        if(lep_p4()[n].Pt()<25.) continue;
        if((lep_p4()[n].Eta())>2.5) continue;
        if(lepSR>=1&&lepSR<=3&&lep_pass_VVV_cutbased_tight()[n]){
          if(id1<0)      id1 = n;
          else           id2 = n;
        } else if(lepSR>=4&&lep_pass_VVV_cutbased_3l_tight()[n]){
          if(id1<0)      id1 = n;
          else if(id2<0) id2 = n;
          else           id3 = n;
        }
      }
      if(lepSR>=1&&lepSR<=3){
        if(id2<0) cout << __LINE__ << endl;
        DPhill     = dPhi(lep_p4()[id1],lep_p4()[id2]);
        DPhillMET  = dPhi(lep_p4()[id1]+lep_p4()[id2],METlv);
        float dphi1 = dPhi(lep_p4()[id1],METlv);
        float dphi2 = dPhi(lep_p4()[id2],METlv);
        DPhilMETmin = (dphi1<dphi2 ? dphi1 : dphi2);
        DPhilMETmax = (dphi1<dphi2 ? dphi2 : dphi1);
        pTll = (lep_p4()[id1]+lep_p4()[id2]).Pt();
        RatMETvsll = (pTll>0 ? met_pt()/pTll : -1);
        if(passSoftJet&&MassJJ<=0&j1>=0&&j2>=0) { MassJJ = (jets_p4()[j1]+jets_p4()[j2]).M(); passMjj = fabs(MassJJ-80.)<15.; }
        for(unsigned int n = 0; n<jets_csv().size();++n){
          if(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5) {
            if(     DPhiljmin<0)                                DPhiljmin = dPhi(jets_p4()[n],lep_p4()[id1]);
            else if(DPhiljmin>dPhi(jets_p4()[n],lep_p4()[id1])) DPhiljmin = dPhi(jets_p4()[n],lep_p4()[id1]);
            if(     DPhiljmin>dPhi(jets_p4()[n],lep_p4()[id2])) DPhiljmin = dPhi(jets_p4()[n],lep_p4()[id2]);
            if(     DRljmin<0)                                DRljmin = dR(jets_p4()[n],lep_p4()[id1]);
            else if(DRljmin>dPhi(jets_p4()[n],lep_p4()[id1])) DRljmin = dR(jets_p4()[n],lep_p4()[id1]);
            if(     DRljmin>dPhi(jets_p4()[n],lep_p4()[id2])) DRljmin = dR(jets_p4()[n],lep_p4()[id2]);
          }
          if(passSoftJet&&jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=3.0) {
            if(Detajj<0&&j1>=0&&j1!=(int)n) { Detajj = dEta(jets_p4()[j1], jets_p4()[n]); passDeta = Detajj<=1.5; passMjjL = ((jets_p4()[j1]+jets_p4()[n]).M()<=400.); }
            if(     DPhiljmin<0)                                DPhiljmin = dPhi(jets_p4()[n],lep_p4()[id1]);
            else if(DPhiljmin>dPhi(jets_p4()[n],lep_p4()[id1])) DPhiljmin = dPhi(jets_p4()[n],lep_p4()[id1]);
            if(     DPhiljmin>dPhi(jets_p4()[n],lep_p4()[id2])) DPhiljmin = dPhi(jets_p4()[n],lep_p4()[id2]);
            if(     DRljmin<0)                                DRljmin = dR(jets_p4()[n],lep_p4()[id1]);
            else if(DRljmin>dPhi(jets_p4()[n],lep_p4()[id1])) DRljmin = dR(jets_p4()[n],lep_p4()[id1]);
            if(     DRljmin>dPhi(jets_p4()[n],lep_p4()[id2])) DRljmin = dR(jets_p4()[n],lep_p4()[id2]);
          }
        }
      }
      if(lepSR>=4){
        float MT1 = mT(lep_p4()[id1],METlv);
        float MT2 = mT(lep_p4()[id2],METlv);
        float MT3 = mT(lep_p4()[id3],METlv);
        LorentzVector thirdlv, SFOSlv;  
        if(nSFOS()==1){
          if(isSFOS01()) { thirdlv = lep_p4()[id3]; SFOSlv = lep_p4()[id1]+lep_p4()[id2]; }
          if(isSFOS02()) { thirdlv = lep_p4()[id2]; SFOSlv = lep_p4()[id1]+lep_p4()[id3]; }
          if(isSFOS12()) { thirdlv = lep_p4()[id1]; SFOSlv = lep_p4()[id2]+lep_p4()[id3]; }
        }
        if(nSFOS()==2){
          if(     !isSFOS01()&&fabs(M02()-MZ)<fabs(M12()-MZ)) { thirdlv = lep_p4()[id2]; SFOSlv = lep_p4()[id1]+lep_p4()[id3]; thirdMT = MT2; }
          else if(!isSFOS01()                               ) { thirdlv = lep_p4()[id1]; SFOSlv = lep_p4()[id2]+lep_p4()[id3]; thirdMT = MT1; }
          if(     !isSFOS02()&&fabs(M01()-MZ)<fabs(M12()-MZ)) { thirdlv = lep_p4()[id3]; SFOSlv = lep_p4()[id1]+lep_p4()[id2]; thirdMT = MT3; }
          else if(!isSFOS02()                               ) { thirdlv = lep_p4()[id1]; SFOSlv = lep_p4()[id2]+lep_p4()[id3]; thirdMT = MT1; }
          if(     !isSFOS12()&&fabs(M01()-MZ)<fabs(M02()-MZ)) { thirdlv = lep_p4()[id3]; SFOSlv = lep_p4()[id1]+lep_p4()[id2]; thirdMT = MT3; }
          else if(!isSFOS12()                               ) { thirdlv = lep_p4()[id2]; SFOSlv = lep_p4()[id1]+lep_p4()[id3]; thirdMT = MT2; } 
        }
        if(MT1>MT2 && MT1>MT3) { maxMT = MT1; if(MT2>MT3) minMT = MT3; else minMT = MT2; }
        if(MT2>MT1 && MT2>MT3) { maxMT = MT2; if(MT1>MT3) minMT = MT3; else minMT = MT1; }
        if(MT3>MT2 && MT3>MT2) { maxMT = MT3; if(MT2>MT1) minMT = MT1; else minMT = MT2; }
        MTsum    = MT1+MT2+MT3;
        MTlll    = mT(lep_p4()[id1]+lep_p4()[id2]+lep_p4()[id3],METlv);
        MSFOS = SFOSlv.M();
        pTSFOS   = SFOSlv.Pt();
        pT3rd    = thirdlv.Pt();
        pTW      = (thirdlv+METlv).Pt();
        pTlllMET = (lep_p4()[id1]+lep_p4()[id2]+lep_p4()[id3]+METlv).Pt();
        DPhiZvsl = dPhi(SFOSlv,thirdlv);
        DPhiZvsW = dPhi(SFOSlv,thirdlv+METlv);
      }



      

      vector<string> prestring;
      if(SRSS[0]>=0||SR3l[0]>=0) prestring.push_back("");
      if(SRSS[1]>=0||SR3l[1]>=0) prestring.push_back("ZPresel_");
      for(unsigned int n = 0; n<prestring.size(); ++n){
        if(SRSS[1]>=0){
          //cout << __LINE__ << endl;
          fillhisto(histos, prestring[n]+"MET"         +SRlep,     sample, sn, MET,         weight);
          fillhisto(histos, prestring[n]+"Mll"         +SRlep,     sample, sn, Mll,         weight);
          fillhisto(histos, prestring[n]+"Mjj"         +SRlep,     sample, sn, MassJJ,      weight);
          fillhisto(histos, prestring[n]+"Detajj"      +SRlep,     sample, sn, Detajj,      weight);
          fillhisto(histos, prestring[n]+"NJets"       +SRlep,     sample, sn, NJets,       weight);
          fillhisto(histos, prestring[n]+"MTmax"       +SRlep,     sample, sn, maxMT,       weight);
          fillhisto(histos, prestring[n]+"MTmin"       +SRlep,     sample, sn, minMT,       weight);
          fillhisto(histos, prestring[n]+"MTsum"       +SRlep,     sample, sn, MTsum,       weight);
          fillhisto(histos, prestring[n]+"DPhill"      +SRlep,     sample, sn, DPhill,      weight);
          fillhisto(histos, prestring[n]+"DPhillMET"   +SRlep,     sample, sn, DPhillMET,   weight);
          fillhisto(histos, prestring[n]+"DPhilMETmin" +SRlep,     sample, sn, DPhilMETmin, weight);
          fillhisto(histos, prestring[n]+"DPhilMETmax" +SRlep,     sample, sn, DPhilMETmax, weight);
          fillhisto(histos, prestring[n]+"pTll"        +SRlep,     sample, sn, pTll,        weight);
          fillhisto(histos, prestring[n]+"DPhiljmin"   +SRlep,     sample, sn, DPhiljmin,   weight);
          fillhisto(histos, prestring[n]+"DRljmin"     +SRlep,     sample, sn, DRljmin,     weight);
          fillhisto(histos, prestring[n]+"RatMETvsll"  +SRlep,     sample, sn, RatMETvsll,  weight);
          //cout << __LINE__ << endl;
        }
        if(SR3l[1]>=0){
          //cout << __LINE__ << " " << SRlep << endl;
          fillhisto(               histos, prestring[n]+"MET"         +SRlep,     sample, sn, MET,         weight);
          fillhisto(               histos, prestring[n]+"DPhilllMET"  +SRlep,     sample, sn, DPhilllMET,  weight);
          fillhisto(               histos, prestring[n]+"pTlll"       +SRlep,     sample, sn, pTlll,       weight);
          if(nSFOS()>=1) fillhisto(histos, prestring[n]+"MT3rd"       +SRlep,     sample, sn, thirdMT,     weight);
          fillhisto(               histos, prestring[n]+"MTlll"       +SRlep,     sample, sn, MTlll,       weight);
          fillhisto(               histos, prestring[n]+"MTmax"       +SRlep,     sample, sn, maxMT,       weight);
          fillhisto(               histos, prestring[n]+"MTmin"       +SRlep,     sample, sn, minMT,       weight);
          fillhisto(               histos, prestring[n]+"MTsum"       +SRlep,     sample, sn, MTsum,       weight);
          if(nSFOS()>=1) fillhisto(histos, prestring[n]+"MSFOSZ"      +SRlep,     sample, sn, MSFOS,       weight);
          if(nSFOS()>=1) fillhisto(histos, prestring[n]+"pTSFOS"      +SRlep,     sample, sn, pTSFOS,      weight);
          if(nSFOS()>=1) fillhisto(histos, prestring[n]+"pT3rd"       +SRlep,     sample, sn, pT3rd,       weight);
          if(nSFOS()>=1) fillhisto(histos, prestring[n]+"pTW"         +SRlep,     sample, sn, pTW,         weight);
          fillhisto(               histos, prestring[n]+"pTlllMET"    +SRlep,     sample, sn, pTlllMET,    weight);
          if(nSFOS()>=1) fillhisto(histos, prestring[n]+"DPhiZvsl"    +SRlep,     sample, sn, DPhiZvsl,    weight);
          if(nSFOS()>=1) fillhisto(histos, prestring[n]+"DPhiZvsW"    +SRlep,     sample, sn, DPhiZvsW,    weight);
          fillhisto(               histos, prestring[n]+"RatMETvslll" +SRlep,     sample, sn, RatMETvslll, weight);
          //cout << __LINE__ << endl;
        }
      }





      vector<bool> cuts;
      vector<string> selstring;
      if(lepSR>=1&&lepSR<=3){
        cuts.push_back(true);                                        //0
        cuts.push_back(MET>30.);                                     //1
        cuts.push_back(MET>40.);                                     //2
        cuts.push_back(MET>60.);                                     //3
        cuts.push_back(MET>80.);                                     //4
        cuts.push_back(MTsum>90.);                                   //5
        cuts.push_back(MTsum>150.);                                  //6
        cuts.push_back(maxMT>90.);                                   //7
        cuts.push_back(maxMT>120.);                                  //8
        cuts.push_back(maxMT>150.);                                  //9
        cuts.push_back(minMT>50.);                                   //10
        cuts.push_back(minMT>90.);                                   //11
        cuts.push_back(minMT>120.);                                  //12
        cuts.push_back(DPhill>0.5);                                  //13
        cuts.push_back(Mll>20.);                                     //14
        cuts.push_back(Mll>40.);                                     //15
        cuts.push_back(Mll>60.);                                     //16
        cuts.push_back(fabs(MassJJ-80.)<10.);                        //17
        cuts.push_back(fabs(MassJJ-80.)<20.);                        //18
        cuts.push_back(MassJJ<65.);                                  //19
        cuts.push_back(MassJJ<100.);                                 //20
        cuts.push_back(MassJJ<120.);                                 //21
        cuts.push_back(RatMETvsll>0.5);                              //22
        cuts.push_back(DPhiljmin<2.0);                               //23
        cuts.push_back(DPhiljmin<1.5);                               //24
        cuts.push_back(DPhiljmin>0.5);                               //25
        cuts.push_back(DPhiljmin<1.5&&DPhiljmin>0.5);                //26
        cuts.push_back(DRljmin<1.5);                                 //27
        cuts.push_back(DRljmin>0.5);                                 //28
        cuts.push_back(DRljmin<1.5&&DRljmin>0.5);                    //29
        cuts.push_back(Detajj<2.5);                                  //30
        cuts.push_back(nj30()>=3);                                   //31
        cuts.push_back(MET>80. && DPhiljmin<1.5);                    //32
        cuts.push_back(MET>80. && DPhiljmin>0.5);                    //33
        cuts.push_back(MET>80. && DPhiljmin<1.5&&DPhiljmin>0.5);     //34
        cuts.push_back(maxMT>90. && DPhiljmin<1.5);                  //35
        cuts.push_back(maxMT>90. && DPhiljmin>0.5);                  //36
        cuts.push_back(maxMT>90. && DPhiljmin<1.5&&DPhiljmin>0.5);   //37
        cuts.push_back(minMT>50. && DPhiljmin<1.5);                  //38
        cuts.push_back(minMT>50. && DPhiljmin>0.5);                  //39
        cuts.push_back(minMT>50. && DPhiljmin<1.5&&DPhiljmin>0.5);   //40
        cuts.push_back(MET>80. && MTsum>90.);                        //41
        cuts.push_back(MET>80. && maxMT>90.);                        //42
        cuts.push_back(MET>80. && minMT>50.);                        //43
        cuts.push_back(minMT>50. && maxMT>90.);                      //44
        cuts.push_back(minMT>50. && maxMT>120.);                     //45
        cuts.push_back(minMT>90. && maxMT>90.);                      //46
        cuts.push_back(minMT>90. && maxMT>120.);                     //47
        cuts.push_back(DPhill>0.5 && minMT>50.);                     //48
        cuts.push_back(DPhill>0.5 && maxMT>90.);                     //49
        if(passNJ && passMjj && passMjjL && passDeta && passMll  && passMET && passMT)   selstring.push_back("SR_addone");
        if(passNJ && passMjj && passMjjL && passDeta && passMll  &&            passMT)   selstring.push_back("SR_addone_dropMET");
        if(passNJ &&            passMjjL && passDeta && passMll  &&            passMT)   selstring.push_back("SR_addone_dropMETMjj");
        if(passNJ && passMjj && passMjjL && passDeta &&                        passMT)   selstring.push_back("SR_addone_dropMETMll");
        if(passNJ &&            passMjjL && passDeta && lepSR==1 &&            passMT)   selstring.push_back("SR_addone_dropMETMjjMll");
        if(passNJ &&            passMjjL && passDeta && passMll  && passMET && passMT)   selstring.push_back("SR_addone_dropMjj");
        if(passNJ && passMjj && passMjjL && passDeta &&             passMET && passMT)   selstring.push_back("SR_addone_dropMll");
        if(passNJ && passMjj && passMjjL &&             passMll  && passMET && passMT)   selstring.push_back("SR_addone_dropDeta");
        if(passNJ &&            passMjjL &&             passMll  && passMET && passMT)   selstring.push_back("SR_addone_dropDetaMjj");
        if(passNJ && passMjj && passMjjL && passDeta && passMll  && passMET && lepSR==2) selstring.push_back("SR_addone_dropMTmax"); 
        if(passNJ &&            passMjjL && passDeta && passMll  && passMET && lepSR==2) selstring.push_back("SR_addone_dropMTmaxMjj"); 
        if(passNJ && passMjj && passMjjL && passDeta &&             passMET && lepSR==2) selstring.push_back("SR_addone_dropMTmaxMll");
        if(passOneJet && passMll && passMET && passMT)   selstring.push_back("SRonejet_addone");
        if(passOneJet && passMll &&            passMT)   selstring.push_back("SRonejet_addone_dropMET");
        if(passOneJet &&                       passMT)   selstring.push_back("SRonejet_addone_dropMETMll");
        if(passOneJet &&            passMET && passMT)   selstring.push_back("SRonejet_addone_dropMll");
        if(passOneJet && passMll && passMET && lepSR==2) selstring.push_back("SRonejet_addone_dropMTmax");
        if(passOneJet &&            passMET && lepSR==2) selstring.push_back("SRonejet_addone_dropMTmaxMll");
        if(passSoftJet && passMjj && passMjjL && passDeta && passMll  && passMET && passMT)   selstring.push_back("SRsoftjet_addone");
        if(passSoftJet && passMjj && passMjjL && passDeta && passMll  &&            passMT)   selstring.push_back("SRsoftjet_addone_dropMET");
        if(passSoftJet &&            passMjjL && passDeta && passMll  &&            passMT)   selstring.push_back("SRsoftjet_addone_dropMETMjj");
        if(passSoftJet && passMjj && passMjjL && passDeta &&                        passMT)   selstring.push_back("SRsoftjet_addone_dropMETMll");
        if(passSoftJet &&            passMjjL && passDeta && lepSR==1 &&            passMT)   selstring.push_back("SRsoftjet_addone_dropMETMjjMll");
        if(passSoftJet &&            passMjjL && passDeta && passMll  && passMET && passMT)   selstring.push_back("SRsoftjet_addone_dropMjj");
        if(passSoftJet && passMjj && passMjjL && passDeta &&             passMET && passMT)   selstring.push_back("SRsoftjet_addone_dropMll");
        if(passSoftJet && passMjj && passMjjL &&             passMll  && passMET && passMT)   selstring.push_back("SRsoftjet_addone_dropDeta");
        if(passSoftJet &&            passMjjL &&             passMll  && passMET && passMT)   selstring.push_back("SRsoftjet_addone_dropDetaMjj");
        if(passSoftJet && passMjj && passMjjL && passDeta && passMll  && passMET && lepSR==2) selstring.push_back("SRsoftjet_addone_dropMTmax");
        if(passSoftJet &&            passMjjL && passDeta && passMll  && passMET && lepSR==2) selstring.push_back("SRsoftjet_addone_dropMTmaxMjj");
        if(passSoftJet && passMjj && passMjjL && passDeta &&             passMET && lepSR==2) selstring.push_back("SRsoftjet_addone_dropMTmaxMll");
        if(passNJ && passMjjside && passMjjL && passDeta && passMll && passMET && passMT)   selstring.push_back("SRsideMjj_addone");
        if(passNJ && passMjjside && passMjjL && passDeta && passMll &&            passMT)   selstring.push_back("SRsideMjj_addone_dropMET");
        if(passNJ && passMjjside && passMjjL && passDeta &&                       passMT)   selstring.push_back("SRsideMjj_addone_dropMETMll");
        if(passNJ && passMjjside && passMjjL && passDeta &&            passMET && passMT)   selstring.push_back("SRsideMjj_addone_dropMll");
        if(passNJ && passMjjside && passMjjL &&             passMll && passMET && passMT)   selstring.push_back("SRsideMjj_addone_dropDeta");
        if(passNJ && passMjjside && passMjjL && passDeta && passMll && passMET && lepSR==2) selstring.push_back("SRsideMjj_addone_dropMTmax");
        if(passNJ && passMjjside && passMjjL && passDeta &&            passMET && lepSR==2) selstring.push_back("SRsideMjj_addone_dropMTmaxMll");
      }//SS
      if(lepSR>=4&&lepSR<=6){
        cuts.push_back(true);                                        //0
        cuts.push_back(MET>30.);                                     //1
        cuts.push_back(MET>40.);                                     //2
        cuts.push_back(MET>60.);                                     //3
        cuts.push_back(MET>80.);                                     //4
        cuts.push_back(pTW>40.);                                     //5
        cuts.push_back(pTW>60.);                                     //6
        cuts.push_back(pT3rd>40.);                                   //7
        cuts.push_back(pT3rd>60.);                                   //8
        cuts.push_back(pTSFOS>40.);                                  //9
        cuts.push_back(pTSFOS>60.);                                  //10
        cuts.push_back(MTsum>100.);                                  //11
        cuts.push_back(MTsum>120.);                                  //12
        cuts.push_back(maxMT>90.);                                   //13
        cuts.push_back(minMT>50.);                                   //14
        cuts.push_back(thirdMT>50.);                                 //15
        cuts.push_back(thirdMT>90.);                                 //16
        cuts.push_back(thirdMT>120.);                                //17
        cuts.push_back(MTlll>90.);                                   //18
        cuts.push_back(DPhiZvsW<2.0);                                //19
        cuts.push_back(DPhiZvsW>2.0);                                //20
        cuts.push_back(pTlll>30.);                                   //21
        cuts.push_back(pTlll>45.);                                   //22
        cuts.push_back(pTlll>60.);                                   //23
        cuts.push_back(pTlll>90.);                                   //24
        cuts.push_back(DPhilllMET>2.1);                              //25
        cuts.push_back(DPhilllMET>2.7);                              //26
        cuts.push_back(RatMETvslll>0.5);                             //27
        cuts.push_back(RatMETvslll>1.0);                             //28
        cuts.push_back(MET>30. && maxMT>90.);                        //29
        cuts.push_back(MET>40. && maxMT>90.);                        //30
        cuts.push_back(MET>55. && maxMT>90.);                        //31
        cuts.push_back(MET>80. && maxMT>90.);                        //32
        cuts.push_back(MET>30. && minMT>50.);                        //33
        cuts.push_back(MET>40. && minMT>50.);                        //34
        cuts.push_back(MET>55. && minMT>50.);                        //35
        cuts.push_back(MET>80. && minMT>50.);                        //36
        cuts.push_back(pTW>40. && maxMT>90.);                        //37
        cuts.push_back(pTW>40. && minMT>50.);                        //38
        cuts.push_back(pTW>40. && thirdMT>90.);                      //39
        if(passNJ && passMET && passMT   && passDPhi && passPTlll)   selstring.push_back("SR_addone");
        if(passNJ &&            passMT   && passDPhi && passPTlll)   selstring.push_back("SR_addone_dropMET");
        if(passNJ && passMET && passMT   &&             passPTlll)   selstring.push_back("SR_addone_dropDPhilllMET");
        if(passNJ && passMET && passMT   && passDPhi && lepSR>=5 )   selstring.push_back("SR_addone_droppTlll");
        if(passNJ && passMET && passMT   &&             lepSR>=5 )   selstring.push_back("SR_addone_droppTlllDPhilllMET");
        if(passNJ && passMET && lepSR==5 && passDPhi && passPTlll)   selstring.push_back("SR_addone_dropMT3rd");
        if(passSoftJet && passMET && passMT   && passDPhi && passPTlll)   selstring.push_back("SRsoftjet_addone");
        if(passSoftJet &&            passMT   && passDPhi && passPTlll)   selstring.push_back("SRsoftjet_addone_dropMET");
        if(passSoftJet && passMET && passMT   &&             passPTlll)   selstring.push_back("SRsoftjet_addone_dropDPhilllMET");
        if(passSoftJet && passMET && passMT   && passDPhi && lepSR>=5 )   selstring.push_back("SRsoftjet_addone_droppTlll");
        if(passSoftJet && passMET && passMT   &&             lepSR>=5 )   selstring.push_back("SRsoftjet_addone_droppTlllDPhilllMET");
        if(passSoftJet && passMET && lepSR==5 && passDPhi && passPTlll)   selstring.push_back("SRsoftjet_addone_dropMT3rd");
      }//3l
      for(unsigned int ss = 0; ss<selstring.size(); ++ss){
        for(unsigned int c = 0; c<cuts.size(); ++c){
          if(cuts[c]) {
            //cout << __LINE__ << " " << selstring[ss]+SRlep+"_"+sample << " " << selstring[ss]+SRlep+"_"+sn << " " << c+0.5 << " " << weight << endl;
            fillhisto(histos, selstring[ss]+SRlep,     sample, sn, c+0.5, weight);//filling all histograms
          }
        }
      }
      //cout << __LINE__ << endl;

    }//event loop
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  SaveHistosToFile("rootfiles/InvestigateSR.root",histos,true,true,(chainnumber==0));
  deleteHistograms(histos);
  
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.02f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.02f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
