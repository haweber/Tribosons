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
//#define USE_CMS3_WWW100 

//#include "Functions112.h"
//#include "CMS3_WWW112.cc"
#include "Functions.h"
#include "CMS3_WWW121.cc"
/*
#include "Functions112.h" 
#ifdef USE_CMS3_WWW100
#include "CMS3_WWW112.cc"
#else
#include "CMS3_WWW0118.cc"
#endif
*/
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1, int year=-1) {

  TFile *fphot;
  TFile *fjet;
  TFile *femjet;
  if(year==2016) fphot  = TFile::Open("/home/users/haweber/L1prefiring_photonpt_2016BtoH.root");
  if(year==2016) fjet   = TFile::Open("/home/users/haweber/L1prefiring_jetpt_2016BtoH.root");
  if(year==2016) femjet = TFile::Open("/home/users/haweber/L1prefiring_jetempt_2016BtoH.root");
  TH2F* hphot  = (TH2F*)fphot ->Get("L1prefiring_photonpt_2016BtoH");
  TH2F* hjet   = (TH2F*)fjet  ->Get("L1prefiring_jetpt_2016BtoH");
  TH2F* hemjet = (TH2F*)femjet->Get("L1prefiring_jetempt_2016BtoH");
  

  vector<myevt> e;
  addeventtocheck(e, 274971,  241, 368840964);
  addeventtocheck(e, 283877, 1272, 2229337525);
  addeventtocheck(e, 283059,   43, 76595225);
  addeventtocheck(e, 283308,  384, 738858987);
  addeventtocheck(e, 283934,  970, 1680056919);
  addeventtocheck(e, 283478,  688, 1052704922);


  bool blindSR         = false;
  bool btagreweighting = false;
  bool applylepSF      = false;
  bool applytrigSF     = false;
  bool applyPUrewgt    = true;
  bool getJECunc       = false;
  
  const char* json_file2016 = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  const char* json_file2017 = "data/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt";
  if(year==2016) set_goodrun_file_json(json_file2016);
  if(year==2017) set_goodrun_file_json(json_file2017);

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  bool storeeventnumbers = true;
  std::ostringstream*  SREE    ;
  std::ostringstream*  SREM    ;
  std::ostringstream*  SRMM    ;
  std::ostringstream*  SREEsb  ;
  std::ostringstream*  SREMsb  ;
  std::ostringstream*  SRMMsb  ;
  std::ostringstream*  SR0SFOS ;
  std::ostringstream*  SR1SFOS ;
  std::ostringstream*  SR2SFOS ;
  std::ostringstream*  AREE    ;
  std::ostringstream*  AREM    ;
  std::ostringstream*  ARMM    ;
  std::ostringstream*  AREEsb  ;
  std::ostringstream*  AREMsb  ;
  std::ostringstream*  ARMMsb  ;
  std::ostringstream*  AR0SFOS ;
  std::ostringstream*  AR1SFOS ;
  std::ostringstream*  AR2SFOS ;
  std::ostringstream*  CREE    ;
  std::ostringstream*  CREM    ;
  std::ostringstream*  CRMM    ;
  std::ostringstream*  CR0SFOS ;
  std::ostringstream*  CR1SFOS ;
  std::ostringstream*  CR2SFOS ;
  std::ostringstream*  dummy1  ;
  std::ostringstream*  dummy2  ;
  std::ostringstream*  dummy3  ;
  std::ostringstream*  dummy4  ;
  std::ostringstream*  dummy5  ;
  std::ostringstream*  dummy6  ;
  if(storeeventnumbers){
    SREE    = new std::ostringstream();
    SREM    = new std::ostringstream();
    SRMM    = new std::ostringstream();
    SREEsb  = new std::ostringstream();
    SREMsb  = new std::ostringstream();
    SRMMsb  = new std::ostringstream();
    SR0SFOS = new std::ostringstream();
    SR1SFOS = new std::ostringstream();
    SR2SFOS = new std::ostringstream();
    AREEsb  = new std::ostringstream();
    AREMsb  = new std::ostringstream();
    ARMMsb  = new std::ostringstream();
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
    dummy1  = new std::ostringstream();
    dummy2  = new std::ostringstream();
    dummy3  = new std::ostringstream();
    dummy4  = new std::ostringstream();
    dummy5  = new std::ostringstream();
    dummy6  = new std::ostringstream();
  }
  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int>    hbins;      hbins.clear();
  vector<float>  hlow;       hlow.clear();
  vector<float>  hup;        hup.clear();

  histonames.push_back("NewSignalRegion");                        hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewSignalRegionPreselection");            hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewLLControlRegion");                     hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewLLControlRegionPreselection");         hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewApplicationRegion");                   hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewApplicationRegionPreselection");       hbins.push_back(9); hlow.push_back(0); hup.push_back(9);

  histonames.push_back("NewSignalRegion_Prefire");                hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewSignalRegion_PrefireUp");              hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewSignalRegion_PrefireDn");              hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewLLControlRegion_Prefire");             hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewLLControlRegion_PrefireUp");           hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewLLControlRegion_PrefireDn");           hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewApplicationRegion_Prefire");           hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewApplicationRegion_PrefireUp");         hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("NewApplicationRegion_PrefireDn");         hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  
  histonames.push_back("SignalRegionTightIso3lv1");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//0SFOS all tight
  histonames.push_back("SignalRegionTightIso3lv2");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//0SFOS only diflavor tight
  histonames.push_back("SignalRegionTightIso3lv3");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//0SFOS only pt<30 tight
  histonames.push_back("SignalRegionTightIso3lv4");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//0SFOS only pt<30 diflavor
  histonames.push_back("SignalRegionTightIso3lv5");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//0SFOS only sublead diflavor tight
  histonames.push_back("SignalRegionTightIso3lv6");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//0SFOS only sublead diflavor tight if pt<30
  histonames.push_back("SignalRegion");                        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("ApplicationRegion");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("WZControlRegion");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegionPresel");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("ApplicationRegionPresel");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("WZControlRegionPresel");               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawSignalRegion");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawApplicationRegion");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawWZControlRegion");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawSignalRegionPresel");               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawApplicationRegionPresel");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawWZControlRegionPresel");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegionv2");                      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("ApplicationRegionv2");                 hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  histonames.push_back("SignalRegion_Mjjside");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("ApplicationRegion_Mjjside");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("WZControlRegion_Mjjside");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegionPresel_Mjjside");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("ApplicationRegionPresel_Mjjside");     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("WZControlRegionPresel_Mjjside");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  histonames.push_back("SignalRegion_JECup");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_JECdn");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_lepSFup");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_lepSFdn");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_bHFSFup");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_bHFSFdn");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_bLFSFup");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_bLFSFdn");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_PUup");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_PUdn");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_trgSFup");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_trgSFdn");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Q2up");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Q2dn");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_aSup");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_aSdn");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_PDFup");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_PDFdn");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  histonames.push_back("SignalRegion_Mjjside_JECup");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_JECdn");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_lepSFup");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_lepSFdn");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_bHFSFup");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_bHFSFdn");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_bLFSFup");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_bLFSFdn");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_PUup");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_PUdn");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_trgSFup");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_trgSFdn");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_Q2up");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_Q2dn");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_aSup");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_aSdn");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_PDFup");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("SignalRegion_Mjjside_PDFdn");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

 
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
  TH1D* h_c;


  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    string fname = currentFile->GetTitle();
    h_c = (TH1D*)file->Get("h_neventsinfile");
    h_c->SetDirectory(0);

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

      double lumi = 1;
      if(year==2016) lumi = 35.9;
      if(year==2017) lumi = 41.9;
      double weight = evt_scale1fb()*lumi;
      //if(nEventsTotal==0) cout << weight << endl; 
      if(firstgoodvertex()!=0)   continue;
      if(nVert()<0)              continue;
      if(nLlep()<2)              continue;
      if(nTlepSS()<1&&nTlep()<2) continue;//preselection can be done already here
      if(nb()!=0)                continue;//preselection can be done already here

      if(string(currentFile->GetTitle()).find("wjets_incl_") !=string::npos && gen_ht()>100.) continue;
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos && gen_ht()>100. && year==2016) continue;
      //if(string(currentFile->GetTitle()).find("www_2l_mia")        !=string::npos) weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")   !=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(isData()) weight = 1.;
      //double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      if(!isData()&&applyPUrewgt)    weight *= purewgt();
      if(!isData()&&applylepSF)      weight *= lepsf();
      if(!isData()&&applytrigSF)     weight *= trigsf();
      bool checkevent = false;
      for(unsigned int i = 0; i<e.size();++i){
        if(e[i].run!=tas::run() ) continue;
        if(e[i].ls !=tas::lumi()) continue;
        if(e[i].evt!=tas::evt() ) continue;
        checkevent = true;
        cout << "Check event " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
        break;
      }
            
      if(isData()){
        if(!passFilters())                      continue;
        if(checkevent) cout << "pass filter"    << endl;
        duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
        if( is_duplicate(id) )                  continue; 
        if(checkevent) cout << "pass duplicate" << endl;
        if( !goodrun(tas::run(), tas::lumi()) ) continue;
        if(checkevent) cout << "pass goodrun"   << endl;
        weight = 1.;
      } 
      if(!passTriggers(true,true)) continue;//pass trigger for data, and offline lepton kinematic cuts for data/simulation
      if(checkevent) cout << "pass online/offline triggers" << endl;
      
      string sample   = skimFilePrefix;
      if(splitVH(fname)){ sample = "WHtoWWW"; }
      string sn = string(bkgtype().Data());
      if(vetophoton()) continue;
      
      int SRSS[30]; 
      int SR3l[30];
      int SRf(-1), ARf(-1), CRf(-1);
      int SRp(-1), ARp(-1), CRp(-1);
      for(int i = 0; i<30; ++i) { SRSS[i] = -1; SR3l[i] = -1;  }

      passAny(SRf,ARf,CRf);
      passAny(SRp,ARp,CRp,true);
      //SS
      passAnySS(SRSS[0],SRSS[2],SRSS[4]);      //full:         0: SR, 2: AR, 4: CR
      passAnySS(SRSS[1],SRSS[3],SRSS[5],true); //preselection: 1: SR, 3: AR, 5: CR
      if(getJECunc) SRSS[6] = isSRSS(false, 1);//6: SR JEC up
      if(getJECunc) SRSS[7] = isSRSS(false,-1);//7: SR JEC dn
      //3l
      passAny3l(SR3l[0],SR3l[2],SR3l[4]);      //full:         0: SR, 2: AR, 4: CR
      passAny3l(SR3l[1],SR3l[3],SR3l[5],true); //preselection: 1: SR, 3: AR, 5: CR
      if(getJECunc) SR3l[6] = isSR3l(false, 1);//6: SR JEC up
      if(getJECunc) SR3l[7] = isSR3l(false,-1);//7: SR JEC dn
      //Mjj side SS
      passAnySS(SRSS[10],SRSS[12],SRSS[14],false,0,false,false,true);      //full:         0: SR, 2: AR, 4: CR
      passAnySS(SRSS[11],SRSS[13],SRSS[15],true, 0,false,false,true); //preselection: 1: SR, 3: AR, 5: CR
      passAnySS(SRSS[16],SRSS[18],SRSS[18],false, 1,false,false,true);      //full, JESup:         0: SR, 2: AR, 4: CR
      passAnySS(SRSS[17],SRSS[18],SRSS[18],false,-1,false,false,true);      //full, JEDdn:         0: SR, 2: AR, 4: CR

      float maxMT = -1;
      LorentzVector METlv; METlv.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      if(nVlep()==3){
        for(int i = 0; i<nVlep();++i){
          if((lep_pass_VVV_cutbased_3l_fo()[i] || lep_pass_VVV_cutbased_3l_tight()[i]) && lep_p4()[i].Pt()>20.) {
            if(mT(lep_p4()[i],METlv)>maxMT) maxMT = mT(lep_p4()[i],METlv);
          }
        }
      }
      if(maxMT<=90.){
        if(SR3l[0]==0) SR3l[0] = -1;
        if(SR3l[1]==0) SR3l[1] = -1;
        if(SR3l[2]==0) SR3l[2] = -1;
        if(SR3l[3]==0) SR3l[3] = -1;
        if(SR3l[4]==0) SR3l[4] = -1;
        if(SR3l[5]==0) SR3l[5] = -1;
        if(SR3l[6]==0) SR3l[6] = -1;
        if(SR3l[7]==0) SR3l[7] = -1;
        if(SRf==6) SRf = -1;
        if(ARf==6) ARf = -1;
        if(CRf==6) CRf = -1;
      }
      if(met_pt()<=60.){
        if(SRf==5) SRf = -1;
        if(ARf==5) ARf = -1;
        if(CRf==5) CRf = -1;
        if(SRSS[10]==2) SRSS[10] = -1;
        if(SRSS[11]==2) SRSS[11] = -1;
        if(SRSS[12]==2) SRSS[12] = -1;
        if(SRSS[13]==2) SRSS[13] = -1;
        //if(SRSS[14]==2) SRSS[14] = -1;
        //if(SRSS[15]==2) SRSS[15] = -1;
        if(SRSS[16]==2) SRSS[16] = -1;
        if(SRSS[17]==2) SRSS[17] = -1;
        //if(SRSS[18]==2) SRSS[18] = -1;
      }
      SR3l[21] = SR3l[0];
      SR3l[22] = SR3l[0];
      SR3l[23] = SR3l[0];
      SR3l[24] = SR3l[0];
      SR3l[25] = SR3l[0];
      SR3l[26] = SR3l[0];
      if(SR3l[0]==0){//start new
        int e1(-1),e2(-1),m1(-1),m2(-1);
        for(int i =0; i<nVlep();++i){
          if(     abs(lep_pdgId()[i])==11&&e1<0) e1 = i;
          else if(abs(lep_pdgId()[i])==11&&e2<0) e2 = i;
          else if(abs(lep_pdgId()[i])==11) cout << "WTF " << __LINE__ << endl;
          if(     abs(lep_pdgId()[i])==13&&m1<0) m1 = i;
          else if(abs(lep_pdgId()[i])==13&&m2<0) m2 = i;
          else if(abs(lep_pdgId()[i])==13) cout << "WTF " << __LINE__ << endl;
        }
        if(e2==-1){//this is a e+mu-mu- event (or qs flipped)
          if(m2<0) cout << "WTF " << __LINE__ << endl;
          if(lep_p4()[m1].Pt()<lep_p4()[m2].Pt()) cout << "WTF " << __LINE__ << endl;
          if(!lep_pass_VVV_cutbased_tight()[e1]) {
            SR3l[21] = -1;
            if(lep_p4()[e1].Pt()<30) SR3l[23] = -1;
          }
          if(!lep_pass_VVV_cutbased_tight()[m1]) {
            SR3l[21] = -1; SR3l[22] = -1;
            if(lep_p4()[m1].Pt()<30){
              SR3l[23] = -1; SR3l[24] = -1;
            }
          }
          if(!lep_pass_VVV_cutbased_tight()[m2]) {
            SR3l[21] = -1; SR3l[22] = -1; SR3l[25] = -1;
            if(lep_p4()[m2].Pt()<30){
              SR3l[23] = -1; SR3l[24] = -1; SR3l[26] = -1;
            }
          }
        }
        else if(m2==-1){ //this is a e+e+mu- event (or qs flipped)
          if(e2<0) cout << "WTF " << __LINE__ << endl;
          if(lep_p4()[e1].Pt()<lep_p4()[e2].Pt()) cout << "WTF " << __LINE__ << endl;
          if(!lep_pass_VVV_cutbased_tight()[m1]) {
            SR3l[21] = -1;
            if(lep_p4()[m1].Pt()<30) SR3l[23] = -1;
          }
          if(!lep_pass_VVV_cutbased_tight()[e1]) {
            SR3l[21] = -1; SR3l[22] = -1;
            if(lep_p4()[e1].Pt()<30){
              SR3l[23] = -1; SR3l[24] = -1;
            }
          }
          if(!lep_pass_VVV_cutbased_tight()[e2]) {
            SR3l[21] = -1; SR3l[22] = -1; SR3l[25] = -1;
            if(lep_p4()[e2].Pt()<30){
              SR3l[23] = -1; SR3l[24] = -1; SR3l[26] = -1;
            }
          }
        }
        else cout << "WTF " << __LINE__ << endl;
      }//end new
      if(weight>50.) cout << __LINE__ << " " << fname << " " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(checkevent) {
        cout << "weight " << weight << " btag  " << weight_btagsf() << " PU " << purewgt() << " trig " << trigeff() << " lep " << lepsf() << endl;
        cout << "nj30 "   << nj30() << " nj " << nj() << " nb " << nb() << " Mjj " << Mjj() << " MjjL " << MjjL() << " Detajj " << DetajjL() << endl;
        cout << "MTmax " << MTmax() << " MTmax3L " << MTmax3L() << " myMATmax " << maxMT << " met " << met_pt() << " metphi " << met_phi() << " DPhi3lMET " << DPhi3lMET() << " Pt3l " << Pt3l() << endl;
        cout << "isSFOS01 " << isSFOS01() << " M01 " << M01() << " isSFOS02 " << isSFOS02() << " M02 " << M02() << " isSFOS12 " << isSFOS12() << " M12 " << M12() << " MeeSS " << MeeSS() << " MllSS " << MllSS() << " Mll3L " << Mll3L() << " Mll3L1 " << Mll3L1() << " Mee3L " << Mee3L() << " Mlll " << M3l() << endl;
        cout << " HLT_DoubleMu " << HLT_DoubleMu() << " HLT_MuEG " << HLT_MuEG() << " HLT_DoubleEl " << HLT_DoubleEl() << " HLT_DoubleEl_DZ " << HLT_DoubleEl_DZ() << " HLT_DoubleEl_DZ_2 " << HLT_DoubleEl_DZ_2() << endl;
        cout << " mc_HLT_DoubleMu " << mc_HLT_DoubleMu() << " mc_HLT_MuEG " << mc_HLT_MuEG() << " mc_HLT_DoubleEl " << mc_HLT_DoubleEl() << " mc_HLT_DoubleEl_DZ " << mc_HLT_DoubleEl_DZ() << " mc_HLT_DoubleEl_DZ_2 " << mc_HLT_DoubleEl_DZ_2() << endl;
        for(unsigned int i = 0; i<jets_p4().size(); ++i){
          cout << "jet pT " << jets_p4()[i].Pt() << " eta " << jets_p4()[i].Eta() << " phi " << jets_p4()[i].Phi() << " CSV " << jets_csv()[i];
          for(unsigned int j = i+1; j<jets_p4().size(); ++j) cout << " M"<<i<<j<< " " << (jets_p4()[i]+jets_p4()[j]).M() << " (dR " << dR(jets_p4()[i],jets_p4()[j]) << ")";
          cout << endl;
        }
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          cout << "lep " << lep_pdgId()[i] << " Pt " << lep_p4()[i].Pt() << " eta " << lep_p4()[i].Eta() << " phi " << lep_p4()[i].Phi();
          cout << " ID SSt/SSl/3lt/3ll " << lep_pass_VVV_cutbased_tight()[i] << "/" << lep_pass_VVV_cutbased_fo()[i] << "/" << lep_pass_VVV_cutbased_3l_tight()[i] << "/" << lep_pass_VVV_cutbased_3l_fo()[i];
          cout << " iso " << lep_relIso03EAv2()[i] << " ip3d " << lep_ip3d()[i] << " t.q " << lep_tightCharge()[i];
          for(unsigned int j = i+1; j<lep_pdgId().size();++j) { cout << " M" << i << j << " " << (lep_p4()[i]+lep_p4()[j]).M();
            for(unsigned int k = j+1; k<lep_pdgId().size();++k) cout << " M" << i << j << k << " " << (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).M() << " Pt " <<  (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).Pt(); }
          cout << endl;
        }
        cout << " weight_fr_r1_f1 " << weight_fr_r1_f1() << " weight_fr_r1_f0p5 " << weight_fr_r1_f0p5() << " weight_fr_r2_f1 " << weight_fr_r2_f1() << " weight_fr_r2_f2 " << weight_fr_r2_f2() << " weight_fr_r2_f0p5 " << weight_fr_r2_f0p5() << " weight_fr_r0p5_f1 " << weight_fr_r0p5_f1() << " weight_fr_r0p5_f2 " << weight_fr_r0p5_f2() << " weight_fr_r0p5_f0p5 " << weight_fr_r0p5_f0p5() << " weight_pdf_up " << weight_pdf_up() << " weight_pdf_down " << weight_pdf_down() << " weight_alphas_down " << weight_alphas_down() << " weight_alphas_up " << weight_alphas_up() << endl;

      }
      double prefire   = 1;
      double prefireup = 1;
      double prefiredn = 1;
      for(unsigned int i = 0; i<jets_p4().size(); ++i){
        if(fabs(jets_p4()[i].Eta())<2.0) continue;
        if(fabs(jets_p4()[i].Eta())>3.0) continue;
        LorentzVector photonsum(0,0,0,0);
        bool isphotonorelectron = false;
        bool print = false;
        for(unsigned int j = 0; j<genPart_pdgId().size(); ++j){
          if(fabs(genPart_p4()[j].Eta())<2.0) continue;
          if(fabs(genPart_p4()[j].Eta())>3.0) continue;
          if(fabs(genPart_pdgId()[j])!=21&&fabs(genPart_pdgId()[j])!=11) continue;
          //if(fabs(genPart_pdgId()[j])==11&&fabs(genPart_p4()[j].Eta())>2.4) continue;
          if(fabs(genPart_pdgId()[j])==21) continue;
          if(genPart_status()[j]!=71&&genPart_status()[j]!=1) continue;
          if(dR(jets_p4()[i],genPart_p4()[j])>0.1) continue;
          if(genPart_p4()[j].Pt()/jets_p4()[i].Pt()<0.85) continue;
          if(genPart_p4()[j].Pt()/jets_p4()[i].Pt()>1.15) continue;
          bool checkrecoele = false;
          if(fabs(genPart_pdgId()[j])==11&&fabs(genPart_p4()[j].Eta())<2.5){
            for(unsigned int k = 0; k<lep_pdgId().size();++k){
              if(abs(lep_pdgId()[k])!=11) continue;
              if(dR(lep_p4()[k],genPart_p4()[j])>0.1) continue;
              if(genPart_p4()[j].Pt()/lep_p4()[k].Pt()<0.85) continue;
              if(genPart_p4()[j].Pt()/lep_p4()[k].Pt()>1.15) continue;
              checkrecoele = true;
              break;
            }
          }
          if(checkrecoele) isphotonorelectron = false;//we reconstruct the lepton itself - use that one 
          else             isphotonorelectron = true;//use photon file
        }
        if(isphotonorelectron) {
          prefire   *= (1.-hphot->GetBinContent(hphot->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt())) );
          prefireup *= (1.-hphot->GetBinContent(hphot->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt()))+hphot->GetBinError(hphot->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt())) );
          prefiredn *= (1.-hphot->GetBinContent(hphot->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt()))-hphot->GetBinError(hphot->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt())) );
        }
        else {
          prefire   *= (1.-hjet ->GetBinContent(hjet ->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt())) );
          prefireup *= (1.-hjet ->GetBinContent(hjet ->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt()))+hjet ->GetBinError(hjet ->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt())) );
          prefiredn *= (1.-hjet ->GetBinContent(hjet ->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt()))-hjet ->GetBinError(hjet ->FindBin(jets_p4()[i].Eta(),jets_p4()[i].Pt())) );
        }
      }
      for(unsigned int i = 0; i<lep_p4().size(); ++i){
        if(fabs(lep_p4()[i].Eta())<2.0) continue;
        if(fabs(lep_p4()[i].Eta())>3.0) continue;
        if(abs(lep_pdgId()[i])!=11) continue;
        prefire   *= (1.-hphot->GetBinContent(hphot->FindBin(lep_p4()[i].Eta(),lep_p4()[i].Pt())) );
        prefireup *= (1.-hphot->GetBinContent(hphot->FindBin(lep_p4()[i].Eta(),lep_p4()[i].Pt()))+hphot->GetBinError(hphot->FindBin(lep_p4()[i].Eta(),lep_p4()[i].Pt())) );
        prefiredn *= (1.-hphot->GetBinContent(hphot->FindBin(lep_p4()[i].Eta(),lep_p4()[i].Pt()))-hphot->GetBinError(hphot->FindBin(lep_p4()[i].Eta(),lep_p4()[i].Pt())) );
      }

      if(!isData()||!blindSR){//SR is blinded
        if(SRf>=0) fillhisto(histos, "NewSignalRegion",               sample, sn, SRf, weight);
        if(SRp>=0) fillhisto(histos, "NewSignalRegionPreselection",   sample, sn, SRp, weight);

        if(SRf>=0) fillhisto(histos, "NewSignalRegion_Prefire",       sample, sn, SRf, weight*prefire);
        if(SRf>=0) fillhisto(histos, "NewSignalRegion_PrefireUp",     sample, sn, SRf, weight*prefireup);
        if(SRf>=0) fillhisto(histos, "NewSignalRegion_PrefireDn",     sample, sn, SRf, weight*prefiredn);
        int SR3lv2S = SR3l[0];
        if(SR3l[0]==0 && maxMT<90.) SR3lv2S = -1;
        fillSRhisto(histos, "SignalRegion",               sample, sn, SRSS[0], SR3l[0], weight);
        fillSRhisto(histos, "SignalRegionTightIso3lv1",   sample, sn, SRSS[0], SR3l[21],weight);
        fillSRhisto(histos, "SignalRegionTightIso3lv2",   sample, sn, SRSS[0], SR3l[22],weight);
        fillSRhisto(histos, "SignalRegionTightIso3lv3",   sample, sn, SRSS[0], SR3l[23],weight);
        fillSRhisto(histos, "SignalRegionTightIso3lv4",   sample, sn, SRSS[0], SR3l[24],weight);
        fillSRhisto(histos, "SignalRegionTightIso3lv5",   sample, sn, SRSS[0], SR3l[25],weight);
        fillSRhisto(histos, "SignalRegionTightIso3lv6",   sample, sn, SRSS[0], SR3l[26],weight);
        fillSRhisto(histos, "SignalRegionv2",             sample, sn, SRSS[0], SR3lv2S, weight);
        fillSRhisto(histos, "SignalRegionPresel",         sample, sn, SRSS[1], SR3l[1], weight);
        fillSRhisto(histos, "RawSignalRegion",            sample, sn, SRSS[0], SR3l[0], 1.);
        fillSRhisto(histos, "RawSignalRegionPresel",      sample, sn, SRSS[1], SR3l[1], 1.);
        fillSRhisto(histos, "SignalRegion_Mjjside",       sample, sn, SRSS[10], -1,     weight);
        fillSRhisto(histos, "SignalRegionPresel_Mjjside", sample, sn, SRSS[11], -1,     weight);
        if(!isData()&&getJECunc) {
          fillSRhisto(histos, "SignalRegion_JECup",       sample, sn, SRSS[6], SR3l[6],     weight);
          fillSRhisto(histos, "SignalRegion_JECdn",       sample, sn, SRSS[7], SR3l[7],     weight);
          fillSRhisto(histos, "SignalRegion_Mjjside_JECup",       sample, sn, SRSS[16], -1, weight);
          fillSRhisto(histos, "SignalRegion_Mjjside_JECdn",       sample, sn, SRSS[17], -1, weight);
        }
        if(!isData()&&applylepSF&&lepsf()!=0){
          fillSRhisto(histos, "SignalRegion_lepSFup",     sample, sn, SRSS[0], SR3l[0],     weight*lepsf_up()/lepsf());
          fillSRhisto(histos, "SignalRegion_lepSFdn",     sample, sn, SRSS[0], SR3l[0],     weight*lepsf_dn()/lepsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_lepSFup",     sample, sn, SRSS[10], -1, weight*lepsf_up()/lepsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_lepSFdn",     sample, sn, SRSS[10], -1, weight*lepsf_dn()/lepsf());
        }
        else if(!isData()&&applylepSF) if(weight!=0) cout << "WTF " << weight << endl;
        if(!isData()&&btagreweighting&&weight_btagsf()!=0){
          fillSRhisto(histos, "SignalRegion_bHFSFup",     sample, sn, SRSS[0], SR3l[0], weight*weight_btagsf_heavy_UP()/weight_btagsf());
          fillSRhisto(histos, "SignalRegion_bHFSFdn",     sample, sn, SRSS[0], SR3l[0], weight*weight_btagsf_heavy_DN()/weight_btagsf());
          fillSRhisto(histos, "SignalRegion_bLFSFup",     sample, sn, SRSS[0], SR3l[0], weight*weight_btagsf_light_UP()/weight_btagsf());
          fillSRhisto(histos, "SignalRegion_bLFSFdn",     sample, sn, SRSS[0], SR3l[0], weight*weight_btagsf_light_DN()/weight_btagsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_bHFSFup",     sample, sn, SRSS[10], -1, weight*weight_btagsf_heavy_UP()/weight_btagsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_bHFSFdn",     sample, sn, SRSS[10], -1, weight*weight_btagsf_heavy_DN()/weight_btagsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_bLFSFup",     sample, sn, SRSS[10], -1, weight*weight_btagsf_light_UP()/weight_btagsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_bLFSFdn",     sample, sn, SRSS[10], -1, weight*weight_btagsf_light_DN()/weight_btagsf());
        }
        if(!isData()&&applyPUrewgt&&!isData()&&purewgt()!=0){
          fillSRhisto(histos, "SignalRegion_PUup",        sample, sn, SRSS[0], SR3l[0],     weight*purewgt_up()/purewgt());
          fillSRhisto(histos, "SignalRegion_PUdn",        sample, sn, SRSS[0], SR3l[0],     weight*purewgt_dn()/purewgt());
          fillSRhisto(histos, "SignalRegion_Mjjside_PUup",        sample, sn, SRSS[10], -1, weight*purewgt_up()/purewgt());
          fillSRhisto(histos, "SignalRegion_Mjjside_PUdn",        sample, sn, SRSS[10], -1, weight*purewgt_dn()/purewgt());
        }
        if(!isData()&&applytrigSF&&!isData()&&trigsf()!=0){
          fillSRhisto(histos, "SignalRegion_trgSFup",        sample, sn, SRSS[0], SR3l[0],     weight*trigsf_up()/trigsf());
          fillSRhisto(histos, "SignalRegion_trgSFdn",        sample, sn, SRSS[0], SR3l[0],     weight*trigsf_dn()/trigsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_trgSFup",        sample, sn, SRSS[10], -1, weight*trigsf_up()/trigsf());
          fillSRhisto(histos, "SignalRegion_Mjjside_trgSFdn",        sample, sn, SRSS[10], -1, weight*trigsf_dn()/trigsf());
        }
        if(weight_fr_r1_f1()!=0){
          fillSRhisto(histos, "SignalRegion_Q2up",          sample, sn, SRSS[0], SR3l[0], weight*weight_fr_r2_f2()/weight_fr_r1_f1()    *h_c->GetBinContent(2)/h_c->GetBinContent(6));
          fillSRhisto(histos, "SignalRegion_Q2dn",          sample, sn, SRSS[0], SR3l[0], weight*weight_fr_r0p5_f0p5()/weight_fr_r1_f1()*h_c->GetBinContent(2)/h_c->GetBinContent(10));
          fillSRhisto(histos, "SignalRegion_Mjjside_Q2up",  sample, sn, SRSS[10], -1,     weight*weight_fr_r2_f2()/weight_fr_r1_f1()    *h_c->GetBinContent(2)/h_c->GetBinContent(6));
          fillSRhisto(histos, "SignalRegion_Mjjside_Q2dn",  sample, sn, SRSS[10], -1,     weight*weight_fr_r0p5_f0p5()/weight_fr_r1_f1()*h_c->GetBinContent(2)/h_c->GetBinContent(10));
          fillSRhisto(histos, "SignalRegion_PDFup",         sample, sn, SRSS[0], SR3l[0], weight*weight_pdf_up()/weight_fr_r1_f1()      *h_c->GetBinContent(2)/h_c->GetBinContent(11));
          fillSRhisto(histos, "SignalRegion_PDFdn",         sample, sn, SRSS[0], SR3l[0], weight*weight_pdf_down()/weight_fr_r1_f1()    *h_c->GetBinContent(2)/h_c->GetBinContent(12));
          fillSRhisto(histos, "SignalRegion_Mjjside_PDFup", sample, sn, SRSS[10], -1,     weight*weight_pdf_up()/weight_fr_r1_f1()      *h_c->GetBinContent(2)/h_c->GetBinContent(11));
          fillSRhisto(histos, "SignalRegion_Mjjside_PDFdn", sample, sn, SRSS[10], -1,     weight*weight_pdf_down()/weight_fr_r1_f1()    *h_c->GetBinContent(2)/h_c->GetBinContent(12));
          fillSRhisto(histos, "SignalRegion_aSup",          sample, sn, SRSS[0], SR3l[0], weight*weight_alphas_up()/weight_fr_r1_f1()   *h_c->GetBinContent(2)/h_c->GetBinContent(14));
          fillSRhisto(histos, "SignalRegion_aSdn",          sample, sn, SRSS[0], SR3l[0], weight*weight_alphas_down()/weight_fr_r1_f1() *h_c->GetBinContent(2)/h_c->GetBinContent(13));
          fillSRhisto(histos, "SignalRegion_Mjjside_aSup",  sample, sn, SRSS[10], -1,     weight*weight_alphas_up()/weight_fr_r1_f1()   *h_c->GetBinContent(2)/h_c->GetBinContent(14));
          fillSRhisto(histos, "SignalRegion_Mjjside_aSdn",  sample, sn, SRSS[10], -1,     weight*weight_alphas_down()/weight_fr_r1_f1() *h_c->GetBinContent(2)/h_c->GetBinContent(13));
        }
      }
      if(CRf>=0) fillhisto(histos, "NewLLControlRegion",               sample, sn, CRf, weight);
      if(CRf>=0) fillhisto(histos, "NewLLControlRegion_Prefire",       sample, sn, CRf, weight*prefire);
      if(CRf>=0) fillhisto(histos, "NewLLControlRegion_PrefireUp",     sample, sn, CRf, weight*prefireup);
      if(CRf>=0) fillhisto(histos, "NewLLControlRegion_PrefireDn",     sample, sn, CRf, weight*prefiredn);
      if(CRp>=0) fillhisto(histos, "NewLLControlRegionPreselection",   sample, sn, CRp, weight);
      if(ARf>=0) fillhisto(histos, "NewApplicationRegion",             sample, sn, ARf, weight);
      if(ARf>=0) fillhisto(histos, "NewApplicationRegion_Prefire",     sample, sn, ARf, weight*prefire);
      if(ARf>=0) fillhisto(histos, "NewApplicationRegion_PrefireUp",   sample, sn, ARf, weight*prefireup);
      if(ARf>=0) fillhisto(histos, "NewApplicationRegion_PrefireDn",   sample, sn, ARf, weight*prefiredn);
      if(ARp>=0) fillhisto(histos, "NewApplicationRegionPreselection", sample, sn, ARp, weight);
      int SR3lv2A = SR3l[2];
      if(SR3l[2]==0 && maxMT<90.) SR3lv2A = -1;
      fillSRhisto(  histos, "ApplicationRegion",          sample, sn, SRSS[2], SR3l[2],     weight);
      fillSRhisto(  histos, "ApplicationRegionv2",        sample, sn, SRSS[2], SR3lv2A,     weight);
      fillSRhisto(  histos, "ApplicationRegionPresel",    sample, sn, SRSS[3], SR3l[3],     weight);
      fillSRhisto(  histos, "ApplicationRegion_Mjjside",          sample, sn, SRSS[12], -1, weight);
      fillSRhisto(  histos, "ApplicationRegionPresel_Mjjside",    sample, sn, SRSS[12], -1, weight);
      fillSRhisto(  histos, "WZControlRegion",            sample, sn, SRSS[4], SR3l[4],     weight);
      fillSRhisto(  histos, "WZControlRegionPresel",      sample, sn, SRSS[5], SR3l[5],     weight);
      fillSRhisto(  histos, "WZControlRegion_Mjjside",            sample, sn, SRSS[14], -1, weight);
      fillSRhisto(  histos, "WZControlRegionPresel_Mjjside",      sample, sn, SRSS[15], -1, weight);
      fillSRhisto(  histos, "RawApplicationRegion",       sample, sn, SRSS[2], SR3l[2], 1.);
      fillSRhisto(  histos, "RawApplicationRegionPresel", sample, sn, SRSS[3], SR3l[3], 1.);
      fillSRhisto(  histos, "RawWZControlRegion",         sample, sn, SRSS[4], SR3l[4], 1.);
      fillSRhisto(  histos, "RawWZControlRegionPresel",   sample, sn, SRSS[5], SR3l[5], 1.);

      if(storeeventnumbers){
        if(fname.find("vh_nonbb_amcnlo")==string::npos || (splitVH(fname))){//don't store other if not WHtoWWW
          addeventtolist(SRSS[0], SR3l[0], SREE, SREM, SRMM, SR0SFOS, SR1SFOS, SR2SFOS);
          addeventtolist(SRSS[2], SR3l[2], AREE, AREM, ARMM, AR0SFOS, AR1SFOS, AR2SFOS);
          addeventtolist(SRSS[4], SR3l[4], CREE, CREM, CRMM, CR0SFOS, CR1SFOS, CR2SFOS);
          //addeventtolist(SRSS[1], SR3l[1], SREE, SREM, SRMM, SR0SFOS, SR1SFOS, SR2SFOS);
          //addeventtolist(SRSS[3], SR3l[3], AREE, AREM, ARMM, AR0SFOS, AR1SFOS, AR2SFOS);
          //addeventtolist(SRSS[5], SR3l[5], CREE, CREM, CRMM, CR0SFOS, CR1SFOS, CR2SFOS);
          addeventtolist(SRSS[10], -1, SREEsb, SREMsb, SRMMsb, dummy1, dummy2, dummy3);
          addeventtolist(SRSS[12], -1, AREEsb, AREMsb, ARMMsb, dummy1, dummy2, dummy3);
        }
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
  if(storeeventnumbers){
    storeeventlist("data/SR", skimFilePrefix, SREE, SREM, SRMM, SR0SFOS, SR1SFOS, SR2SFOS);
    storeeventlist("data/AR", skimFilePrefix, AREE, AREM, ARMM, AR0SFOS, AR1SFOS, AR2SFOS);
    storeeventlist("data/CR", skimFilePrefix, CREE, CREM, CRMM, CR0SFOS, CR1SFOS, CR2SFOS);
    storeeventlist("data/SRsb", skimFilePrefix, SREEsb, SREMsb, SRMMsb, dummy1, dummy2, dummy3);
    storeeventlist("data/ARsb", skimFilePrefix, AREEsb, AREMsb, ARMMsb, dummy4, dummy5, dummy6);
  }

  if(year==2016) SaveHistosToFile("rootfiles/SRLooper_prefire_2016.root",histos,true,true,(chainnumber==0));
  if(year==2017) SaveHistosToFile("rootfiles/SRLooper_prefire_2017.root",histos,true,true,(chainnumber==0));
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
