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
  int counterSR(0), counter(0);
  double timerSR(0), timer(0);
  vector<myevt> e;
  addeventtocheck(e,1, 2842, 1443084);

  TFile *fpu = TFile::Open("rootfiles/puWeights.root");
  map<string, TH1D*> hPU;
  hPU["xs80_observed"  ] = (TH1D*)fpu->Get("puWeight_80_observed");
  hPU["xs69p2_observed"] = (TH1D*)fpu->Get("puWeight_69p2_observed");
  hPU["xs63"] = (TH1D*)fpu->Get("puWeight_63");
  hPU["xs64"] = (TH1D*)fpu->Get("puWeight_64");
  hPU["xs65"] = (TH1D*)fpu->Get("puWeight_65");
  hPU["xs66"] = (TH1D*)fpu->Get("puWeight_66");
  hPU["xs67"] = (TH1D*)fpu->Get("puWeight_67");
  hPU["xs68"] = (TH1D*)fpu->Get("puWeight_68");
  hPU["xs69p2"] = (TH1D*)fpu->Get("puWeight_69p2");
  hPU["xs70"] = (TH1D*)fpu->Get("puWeight_70");
  hPU["xs71"] = (TH1D*)fpu->Get("puWeight_71");
  hPU["xs72"] = (TH1D*)fpu->Get("puWeight_72");
  hPU["xs73"] = (TH1D*)fpu->Get("puWeight_73");
  hPU["xs74"] = (TH1D*)fpu->Get("puWeight_74");
  hPU["xs75"] = (TH1D*)fpu->Get("puWeight_75");
  hPU["xs80"] = (TH1D*)fpu->Get("puWeight_80");

  bool blindSR         = true;
  bool btagreweighting = false;
  bool applylepSF      = false;
  bool applytrigSF     = false;
  bool applyPUrewgt    = false;
  bool getJECunc       = false;
  
  const char* json_file = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  bool storeeventnumbers = true;
  std::ostringstream*  SREE    ;
  std::ostringstream*  SREM    ;
  std::ostringstream*  SRMM    ;
  std::ostringstream*  SR0SFOS ;
  std::ostringstream*  SR1SFOS ;
  std::ostringstream*  SR2SFOS ;
  std::ostringstream*  AREE    ;
  std::ostringstream*  AREM    ;
  std::ostringstream*  ARMM    ;
  std::ostringstream*  AR0SFOS ;
  std::ostringstream*  AR1SFOS ;
  std::ostringstream*  AR2SFOS ;
  std::ostringstream*  CREE    ;
  std::ostringstream*  CREM    ;
  std::ostringstream*  CRMM    ;
  std::ostringstream*  CR0SFOS ;
  std::ostringstream*  CR1SFOS ;
  std::ostringstream*  CR2SFOS ;
  if(storeeventnumbers){
    SREE    = new std::ostringstream();
    SREM    = new std::ostringstream();
    SRMM    = new std::ostringstream();
    SR0SFOS = new std::ostringstream();
    SR1SFOS = new std::ostringstream();
    SR2SFOS = new std::ostringstream();
    AREE    = new std::ostringstream();
    AREM    = new std::ostringstream();
    ARMM    = new std::ostringstream();
    AR0SFOS = new std::ostringstream();
    AR1SFOS = new std::ostringstream();
    AR2SFOS = new std::ostringstream();
    CREE    = new std::ostringstream();
    CREM    = new std::ostringstream();
    CRMM    = new std::ostringstream();
    CR0SFOS = new std::ostringstream();
    CR1SFOS = new std::ostringstream();
    CR2SFOS = new std::ostringstream();
  }
  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int>    hbins;      hbins.clear();
  vector<float>  hlow;       hlow.clear();
  vector<float>  hup;        hup.clear();

  histonames.push_back("NVtx_xs69p2_observed");             hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs80_observed");               hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs63");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs64");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs65");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs66");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs67");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs68");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs69p2");                      hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs70");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs71");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs72");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs73");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs74");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs75");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_xs80");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);
  histonames.push_back("NVtx_raw");                        hbins.push_back(80); hlow.push_back(0); hup.push_back(80);

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

      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_") !=string::npos && gen_ht()>100.) continue;
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos && gen_ht()>100.) continue;
      //if(string(currentFile->GetTitle()).find("www_2l_mia")        !=string::npos) weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")   !=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(isData()) weight = 1.;
      //double rawweight = weight;
      bool checkevent = false;
            
      if(isData()){
        if(!passFilters())                      continue;
        if(checkevent) cout << "pass filter"    << endl;
        duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
        if( is_duplicate(id) )                  continue; 
        if(checkevent) cout << "pass duplicate" << endl;
        if( !goodrun(tas::run(), tas::lumi()) ) continue;
        if(checkevent) cout << "pass goodrun"   << endl;
      } 
      if(!passTriggers()) continue;//pass trigger for data, and offline lepton kinematic cuts for data/simulation
      if(checkevent) cout << "pass online/offline triggers" << endl;
      
      string sample   = skimFilePrefix;
      if(splitVH(fname)){ sample = "WHtoWWW"; }
      string sn = string(bkgtype().Data());
      if(vetophoton()) continue;
      if(isData()){
        fillhisto(  histos, "NVtx_raw",             sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs63",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs64",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs65",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs66",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs67",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs68",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs69p2",          sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs70",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs71",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs72",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs73",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs74",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs75",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs80",            sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs69p2_observed", sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs80_observed",   sample, sn, nVert(),     weight);
      } else {
        //nvertex = 0 starts at bin 1 --> bin = nVert()+1
        int bin = hPU["xs69p2"]->FindBin(nTrueInt());
        int binX = hPU["xs69p2"]->FindBin(nVert()-1);
        fillhisto(  histos, "NVtx_raw",             sample, sn, nVert(),     weight);
        fillhisto(  histos, "NVtx_xs63",            sample, sn, nVert(),     weight*hPU["xs63"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs64",            sample, sn, nVert(),     weight*hPU["xs64"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs65",            sample, sn, nVert(),     weight*hPU["xs65"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs66",            sample, sn, nVert(),     weight*hPU["xs66"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs67",            sample, sn, nVert(),     weight*hPU["xs67"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs68",            sample, sn, nVert(),     weight*hPU["xs68"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs69p2",          sample, sn, nVert(),     weight*hPU["xs69p2"]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs70",            sample, sn, nVert(),     weight*hPU["xs70"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs71",            sample, sn, nVert(),     weight*hPU["xs71"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs72",            sample, sn, nVert(),     weight*hPU["xs72"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs73",            sample, sn, nVert(),     weight*hPU["xs73"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs74",            sample, sn, nVert(),     weight*hPU["xs74"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs75",            sample, sn, nVert(),     weight*hPU["xs75"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs80",            sample, sn, nVert(),     weight*hPU["xs80"  ]->GetBinContent(bin));
        fillhisto(  histos, "NVtx_xs69p2_observed", sample, sn, nVert(),     weight*hPU["xs69p2_observed"]->GetBinContent(binX));
        fillhisto(  histos, "NVtx_xs80_observed",   sample, sn, nVert(),     weight*hPU["xs80_observed"  ]->GetBinContent(binX));
      }

      
      if(checkevent) cout << endl;
      
    }//event loop
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  SaveHistosToFile("rootfiles/PUTester.root",histos,true,true,(chainnumber==0));
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
