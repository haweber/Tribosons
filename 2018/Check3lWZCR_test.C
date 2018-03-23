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
#include "Functions_test.h"
#include "CMS3.cc"
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1) {

  bool blindSR = true;
  bool btagreweighting = true;
  bool applylepSF      = true;
  bool applytrigSF     = true;
  bool applyPUrewgt    = true;
  
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

  //SR
  histonames.push_back("YieldsSR");                                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_raw");                               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_rawweight");                         hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR");                                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_noMSFOSsel_SS");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  //SS JEC uncertainty
  histonames.push_back("YieldsSR_jesup");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_jesdn");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj_jesup");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj_jesdn");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_jesup");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_jesdn");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_cutonMjj");                          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_cutonMjj_jesup");                    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_cutonMjj_jesdn");                    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Mjj_CRlike_allSS_jesup");                     hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_allSS_jesdn");                     hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_allSS");                           hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  //PU unc
  histonames.push_back("YieldsSR_PUup");                              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_PUdn");                              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj_PUup");                      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj_PUdn");                      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_PUup");                              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_PUdn");                              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  //MSFOS/lep SF uncertainty
  histonames.push_back("YieldsCR_lepSFup");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_lepSFdn");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lepSFup");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lepSFdn");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("MSFOS_CRlike_allSS");                         hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_CRlike_allSS_lepSFup");                 hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_CRlike_allSS_lepSFdn");                 hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_CRlike_allSS");                      hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_CRlike_allSS_lepSFup");              hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_CRlike_allSS_lepSFdn");              hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  //validation SS
  histonames.push_back("YieldsSR_Mjjsideband");                       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_Mjjsideband_lowMET");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_Mjjsideband_lowMTmax");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  //validation 3l
  histonames.push_back("YieldsSR_inverteitherMETdPhiPt");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertMETdPhiPt");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_inverteitherMETdPhiPt");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertMETdPhiPt");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertMT3rd");                       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//new
  histonames.push_back("YieldsSR_inverteitherMT3rdMETdPhiPt");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//new
  histonames.push_back("YieldsSR_invertMT3rdMETdPhiPt");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//new
  histonames.push_back("YieldsCR_invertMT3rd");                       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//new
  histonames.push_back("YieldsCR_inverteitherMT3rdMETdPhiPt");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//new
  histonames.push_back("YieldsCR_invertMT3rdMETdPhiPt");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//new
  histonames.push_back("MSFOS_all3l");                                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_lepSFup");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_lepSFdn");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_inverteitherMETdPhiPt");          hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_invertMETdPhiPt");                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l");                             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_lepSFup");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_lepSFdn");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_inverteitherMETdPhiPt");       hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_invertMETdPhiPt");             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_1SFOS");                                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_1SFOS_lepSFup");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_1SFOS_lepSFdn");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_1SFOS_inverteitherMETdPhiPt");          hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_1SFOS_invertMETdPhiPt");                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_1SFOS_invertMT3rd");                    hbins.push_back(12);hlow.push_back(30);hup.push_back(150);//new
  histonames.push_back("MSFOS_1SFOS_inverteitherMT3rdMETdPhiPt");     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);//new
  histonames.push_back("MSFOS_1SFOS_invertMT3rdMETdPhiPt");           hbins.push_back(12);hlow.push_back(30);hup.push_back(150);//new
  histonames.push_back("MSFOS_2SFOS");                                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_lepSFup");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_lepSFdn");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_inverteitherMETdPhiPt");          hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_invertMETdPhiPt");                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS");                             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_lepSFup");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_lepSFdn");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_inverteitherMETdPhiPt");       hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_invertMETdPhiPt");             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);

  //CHECK OLD
  histonames.push_back("OldYieldsSR_dropMjj");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("OldYieldsSR_dropMjj_jesup");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("OldYieldsSR_dropMjj_jesdn");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("OldYieldsSR_Mjjsideband");                       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("OldYieldsSR_Mjjsideband_lowMET");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("OldYieldsSR_Mjjsideband_lowMTmax");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  
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

      if(nVert()<0)              continue;
      if(firstgoodvertex()!=0)   continue;
      if(nLlep()<2)              continue;
      
      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_") !=string::npos && gen_ht()>100.) continue;
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos && gen_ht()>100.) continue;
      //if(string(currentFile->GetTitle()).find("www_2l_mia")        !=string::npos) weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")   !=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(isData()) weight = 1.;
      double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      if(!isData()&&applyPUrewgt)    weight *= purewgt();
      if(!isData()&&applylepSF)      weight *= lepsf();
      if(!isData()&&applytrigSF)     weight *= trigeff();

      //need this for test only - old
      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());//for old
      LorentzVector MET_up; MET_up.SetPxPyPzE(met_up_pt()*TMath::Cos(met_up_phi()),met_up_pt()*TMath::Sin(met_up_phi()),0,met_up_pt());//for old
      LorentzVector MET_dn; MET_dn.SetPxPyPzE(met_dn_pt()*TMath::Cos(met_dn_phi()),met_dn_pt()*TMath::Sin(met_dn_phi()),0,met_dn_pt());//for old
      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      getleptonindices_v2(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l,25,20);//to be deleted old

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
      
      int SRSS[40];
      int SR3l[40];
      int dummy = -1;
      int SRSSold[40]; bool selects3l[40];
      int SR3lold[40];
      for(int i = 0; i<40; ++i) { SRSS[i] = -1; SR3l[i] = -1; }
      for(int i = 0; i<40; ++i) { SRSSold[i] = -1; SR3lold[i] = -1; selects3l[i] = false; }
      //SS
      passAnySS(SRSS[ 0],dummy,SRSS[ 2]);      //full:         0: SR, 2: CR
      passAnySS(SRSS[ 1],dummy,SRSS[ 3],true); //preselection: 1: SR, 3: CR
      SRSS[4] = isCRSS(false,0, true);//full CR - noZ
      SRSS[5] = isCRSS(true ,0, true);//preselection CR - noZ
      passAnySS(SRSS[ 6],dummy,SRSS[ 8],false, 1);//JEC UP
      passAnySS(SRSS[ 7],dummy,SRSS[ 9],true,  1);
      passAnySS(SRSS[10],dummy,SRSS[12],false,-1);//JEC DOWN
      passAnySS(SRSS[11],dummy,SRSS[13],true, -1);
      //4: CR - noZ - old
      SRSSold[ 4] = isCRSS_v0(iSS,i3l,  v3l,false,MTmax(),nj30(),nb(),Mjj(),MjjL(),DetajjL(),MET,0,true, false,false,1);
      //5: CR - noZ preselect - old
      SRSSold[ 5] = isCRSS_v0(iSS,i3l,  v3l,true ,MTmax(),nj30(),nb(),Mjj(),MjjL(),DetajjL(),MET,0,true, false,false,1);
      selects3l[4] = true; selects3l[5] = true;// - old
      bool fail = false;
      if(nVlep()<4){// - old
	if(SRSS[4]!=SRSSold[4]) { cout << "CR-noZ SS full   " << SRSS[4] << " " << SRSSold[4] << endl; fail = true; }
	if(SRSS[5]!=SRSSold[5]) { cout << "CR-noZ SS presel " << SRSS[5] << " " << SRSSold[5] << endl; fail = true; }
      }
      if(vetophoton()) continue;

      //Mjj sideband - like SR id + 20.
      //bool   passAnySS(int &SR, int &AR, int &CR, bool preselect=false, int jec=0, bool noZ=false, bool btag=false, bool Mjjside=false, int version=1);
      passAnySS(SRSS[20],dummy,SRSS[22],false, 0,false,false,true,1);
      passAnySS(SRSS[26],dummy,SRSS[28],false, 1,false,false,true,1);
      passAnySS(SRSS[30],dummy,SRSS[32],false,-1,false,false,true,1);
  
      //3l
      passAny3l(SR3l[ 0],dummy,SR3l[ 2]);      //full:         0: SR, 2 CR
      passAny3l(SR3l[ 1],dummy,SR3l[ 3],true); //preselection: 1: SR, 3 CR
      passAny3l(SR3l[ 6],dummy,SR3l[ 8],false, 1);//JEC UP
      passAny3l(SR3l[ 7],dummy,SR3l[ 9],true,  1);
      passAny3l(SR3l[10],dummy,SR3l[12],false,-1);//JEC DOWN
      passAny3l(SR3l[11],dummy,SR3l[13],true, -1);
      if(vetophoton()) continue;

      //old
      vector <float> MSFOSvecold;
      if(SRSS[5]>=0||SR3l[3]>=0||SR3l[1]>=0) MSFOSvecold = allMSFOS(i3l);//rewrite using MSFOS
      vector <float> MSFOSvec = allMSFOS();
      bool checkMSFOS = false;
      if(SRSS[5]>=0||SR3l[3]>=0||SR3l[1]>=0){
	if(MSFOSvecold.size()!=MSFOSvec.size()) checkMSFOS = true;
	else {
	  for(unsigned int i = 0; i<MSFOSvec.size(); ++i){
	    if(MSFOSvecold[i]==MSFOSvec[i]) continue;
	    checkMSFOS = true;
	    break;
	  }
	}
      }
      if(checkMSFOS){//to be deleted
	cout << "New MSFOSs: "; for(unsigned int i = 0; i<MSFOSvec   .size(); ++i) cout << MSFOSvec[i]    << " "; cout << endl;
	cout << "Old MSFOSs: "; for(unsigned int i = 0; i<MSFOSvecold.size(); ++i) cout << MSFOSvecold[i] << " "; cout << endl;
      }
      //cout << sn << endl;
      if(!(blindSR&&isData())){
	fillSRhisto(histos, "YieldsSR",           sample, sn, sn, SRSS[ 0], SR3l[ 0], weight, weight);
	fillSRhisto(histos, "YieldsSR_raw",       sample, sn, sn, SRSS[ 0], SR3l[ 0], 1., 1.);
	fillSRhisto(histos, "YieldsSR_rawweight", sample, sn, sn, SRSS[ 0], SR3l[ 0], rawweight, rawweight);
	fillSRhisto(histos, "YieldsSR_jesup",     sample, sn, sn, SRSS[ 6], SR3l[ 6], weight, weight);
	fillSRhisto(histos, "YieldsSR_jesdn",     sample, sn, sn, SRSS[10], SR3l[10], weight, weight);
	fillSRhisto(histos, "YieldsSR_PUup",      sample, sn, sn, SRSS[ 0], SR3l[ 0], weight*purewgt_up()/purewgt(), weight*purewgt_up()/purewgt());
	fillSRhisto(histos, "YieldsSR_PUdn",      sample, sn, sn, SRSS[ 0], SR3l[ 0], weight*purewgt_dn()/purewgt(), weight*purewgt_dn()/purewgt());
	if(SR3l[ 0]>=0||SRSS[ 1]>=0){
	  int t = SRSS[ 1];
	  if(MjjL()>400.||DetajjL()>1.5) t = -1;
	  if(SRSS[ 1]==0&&(met_pt()<=60.||MllSS()<=40.)              ) t = -1;
	  if(SRSS[ 1]==1&&(met_pt()<=60.||MllSS()<=30.||MTmax()<=90.)) t = -1;
	  if(SRSS[ 1]==2&&(               MllSS()<=40.)              ) t = -1;
	  fillSRhisto(histos, "OldYieldsSR_dropMjj",        sample, sn, sn, t, SR3l[ 0], weight, weight);
	}
	fillSRhisto(histos, "YieldsSR_dropMjj",        sample, sn, sn, (SRSS[0]>=0 ? SRSS[0] : SRSS[20]), SR3l[ 0], weight, weight);
	fillSRhisto(histos, "YieldsSR_dropMjj_PUup",   sample, sn, sn, (SRSS[0]>=0 ? SRSS[0] : SRSS[20]), SR3l[ 0], weight*purewgt_up()/purewgt(), weight*purewgt_up()/purewgt());
	fillSRhisto(histos, "YieldsSR_dropMjj_PUdn",   sample, sn, sn, (SRSS[0]>=0 ? SRSS[0] : SRSS[20]), SR3l[ 0], weight*purewgt_dn()/purewgt(), weight*purewgt_dn()/purewgt());
	if(SR3l[ 6]>=0||SRSS[ 7]>=0){
	  int t = SRSS[ 7];
	  if(MjjL_up()>400.||DetajjL_up()>1.5) t = -1;
	  if(SRSS[ 7]==0&&(met_up_pt()<=60.||MllSS()<=40.)                 ) t = -1;
	  if(SRSS[ 7]==1&&(met_up_pt()<=60.||MllSS()<=30.||MTmax_up()<=90.)) t = -1;
	  if(SRSS[ 7]==2&&(                  MllSS()<=40.)                 ) t = -1;
	  fillSRhisto(histos, "OldYieldsSR_dropMjj_jesup",   sample, sn, sn, t, SR3l[ 6], weight, weight);
	}
	fillSRhisto(histos, "YieldsSR_dropMjj_jesup",        sample, sn, sn, (SRSS[6]>=0 ? SRSS[6] : SRSS[26]), SR3l[ 6], weight, weight);
	if(SR3l[10]>=0||SRSS[11]>=0){
	  int t = SRSS[11];
	  if(MjjL_dn()>400.||DetajjL_dn()>1.5) t = -1;
	  if(SRSS[11]==0&&(met_dn_pt()<=60.||MllSS()<=40.)                 ) t = -1;
	  if(SRSS[11]==1&&(met_dn_pt()<=60.||MllSS()<=30.||MTmax_dn()<=90.)) t = -1;
	  if(SRSS[11]==2&&(                  MllSS()<=40.)                 ) t = -1;
	  fillSRhisto(histos, "OldYieldsSR_dropMjj_jesdn",   sample, sn, sn, t, SR3l[10], weight, weight);
	}
	fillSRhisto(histos, "YieldsSR_dropMjj_jesdn",        sample, sn, sn, (SRSS[10]>=0 ? SRSS[10] : SRSS[30]), SR3l[10], weight, weight);
	
	fillSRhisto(histos, "YieldsSR_lepSFup",     sample, sn, sn, SRSS[ 0], SR3l[ 0], weight*lepsf_up()/lepsf(), weight*lepsf_up()/lepsf());
	fillSRhisto(histos, "YieldsSR_lepSFdn",     sample, sn, sn, SRSS[ 0], SR3l[ 0], weight*lepsf_dn()/lepsf(), weight*lepsf_dn()/lepsf());
	//MSFOS - SR blinded
	if(SR3l[ 0]>=0&&MSFOSvec.size()>0){
	  fillhisto(histos, "MSFOS_all3l",         sample, sn, MSFOSvec[0], weight);
	  fillhisto(histos, "MSFOS_all3l_lepSFup", sample, sn, MSFOSvec[0], weight*lepsf_up()/lepsf());
	  fillhisto(histos, "MSFOS_all3l_lepSFdn", sample, sn, MSFOSvec[0], weight*lepsf_dn()/lepsf());
	  if(SR3l[ 0]==1){
	    fillhisto(histos, "MSFOS_1SFOS",         sample, sn, MSFOSvec[0], weight);
	    fillhisto(histos, "MSFOS_1SFOS_lepSFup", sample, sn, MSFOSvec[0], weight*lepsf_up()/lepsf());
	    fillhisto(histos, "MSFOS_1SFOS_lepSFdn", sample, sn, MSFOSvec[0], weight*lepsf_dn()/lepsf());
	  }
	  if(SR3l[ 0]==2){
	    fillhisto(histos, "MSFOS_2SFOS",         sample, sn, MSFOSvec[0], weight);
	    fillhisto(histos, "MSFOS_2SFOS_lepSFup", sample, sn, MSFOSvec[0], weight*lepsf_up()/lepsf());
	    fillhisto(histos, "MSFOS_2SFOS_lepSFdn", sample, sn, MSFOSvec[0], weight*lepsf_dn()/lepsf());
	  }
	  for(unsigned int i = 0; i<MSFOSvec.size();++i){
	    fillhisto(histos, "MSFOSall_all3l",         sample, sn, MSFOSvec[i], weight);
	    fillhisto(histos, "MSFOSall_all3l_lepSFup", sample, sn, MSFOSvec[i], weight*lepsf_up()/lepsf());
	    fillhisto(histos, "MSFOSall_all3l_lepSFdn", sample, sn, MSFOSvec[i], weight*lepsf_dn()/lepsf());
	    if(SR3l[ 0]==2){
	      fillhisto(histos, "MSFOSall_2SFOS",         sample, sn, MSFOSvec[i], weight);
	      fillhisto(histos, "MSFOSall_2SFOS_lepSFup", sample, sn, MSFOSvec[i], weight*lepsf_up()/lepsf());
	      fillhisto(histos, "MSFOSall_2SFOS_lepSFdn", sample, sn, MSFOSvec[i], weight*lepsf_dn()/lepsf());
	    }
	  }
	}

      }
      //MSFOS - CR not blinded
      if(SR3l[ 2]>=0&&MSFOSvec.size()>0){
	fillhisto(histos, "MSFOS_all3l",         sample, sn, MSFOSvec[0], weight);
	fillhisto(histos, "MSFOS_all3l_lepSFup", sample, sn, MSFOSvec[0], weight*lepsf_up()/lepsf());
	fillhisto(histos, "MSFOS_all3l_lepSFdn", sample, sn, MSFOSvec[0], weight*lepsf_dn()/lepsf());
	if(SR3l[ 2]==1){
	  fillhisto(histos, "MSFOS_1SFOS",         sample, sn, MSFOSvec[0], weight);
	  fillhisto(histos, "MSFOS_1SFOS_lepSFup", sample, sn, MSFOSvec[0], weight*lepsf_up()/lepsf());
	  fillhisto(histos, "MSFOS_1SFOS_lepSFdn", sample, sn, MSFOSvec[0], weight*lepsf_dn()/lepsf());
	}
	if(SR3l[ 2]==2){
	  fillhisto(histos, "MSFOS_2SFOS",         sample, sn, MSFOSvec[0], weight);
	  fillhisto(histos, "MSFOS_2SFOS_lepSFup", sample, sn, MSFOSvec[0], weight*lepsf_up()/lepsf());
	  fillhisto(histos, "MSFOS_2SFOS_lepSFdn", sample, sn, MSFOSvec[0], weight*lepsf_dn()/lepsf());
	}
	for(unsigned int i = 0; i<MSFOSvec.size();++i){
	  fillhisto(histos, "MSFOSall_all3l",         sample, sn, MSFOSvec[i], weight);
	  fillhisto(histos, "MSFOSall_all3l_lepSFup",  sample, sn, MSFOSvec[i], weight*lepsf_up()/lepsf());
	  fillhisto(histos, "MSFOSall_all3l_lepSFdn", sample, sn, MSFOSvec[i], weight*lepsf_dn()/lepsf());
	  if(SR3l[ 2]==2){
	    fillhisto(histos, "MSFOSall_2SFOS",         sample, sn, MSFOSvec[i], weight);
	    fillhisto(histos, "MSFOSall_2SFOS_lepSFup", sample, sn, MSFOSvec[i], weight*lepsf_up()/lepsf());
	    fillhisto(histos, "MSFOSall_2SFOS_lepSFdn", sample, sn, MSFOSvec[i], weight*lepsf_dn()/lepsf());
	  }
	}
      }
      //Mjj sideband
      int t1(-1),t2(-1);
      if(SRSS[ 1]>=0){//this if to be deleted - OLD
	int t = SRSS[ 1];
	if(MjjL()>400.||DetajjL()>1.5) t = -1;
	if(fabs(Mjj()     -80.)<=15.)  t = -1;
	if(SRSS[ 1]==0&&(met_pt()<=60.||MllSS()<=40.)              ) t = -1;
	if(SRSS[ 1]==1&&(met_pt()<=60.||MllSS()<=30.||MTmax()<=90.)) t = -1;
	if(SRSS[ 1]==2&&(               MllSS()<=40.)              ) t = -1;
	fillSRhisto(histos, "OldYieldsSR_Mjjsideband",   sample, sn, sn, t, -1, weight, weight);
	t = SRSS[ 1];
	if(MjjL()>400.||DetajjL()>1.5) t = -1;
	if(fabs(Mjj()     -80.)<=15.)  t = -1;
	if(SRSS[ 1]==0&&(met_pt()>40.||MllSS()<=40.)              ) t = -1;
	if(SRSS[ 1]==1&&(met_pt()>40.||MllSS()<=30.||MTmax()<=90.)) t = -1;
	if(SRSS[ 1]==2&&(met_pt()>40.||MllSS()<=40.)              ) t = -1;
	t1 = t;
	fillSRhisto(histos, "OldYieldsSR_Mjjsideband_lowMET",   sample, sn, sn, t, -1, weight, weight,true);
	t = SRSS[ 1];
	if(SRSS[ 1]!=1) t = -1;
	if(MjjL()>400.||DetajjL()>1.5) t = -1;
	if(fabs(Mjj()     -80.)<15.)  t = -1;
	if(SRSS[ 1]==1&&(met_pt()<=60.||MllSS()<=30.||MTmax()>90.)) t = -1;
	fillSRhisto(histos, "OldYieldsSR_Mjjsideband_lowMTmax",   sample, sn, sn, t, -1, weight, weight);
      }
      fillSRhisto(histos, "YieldsSR_Mjjsideband",   sample, sn, sn, SRSS[20], -1, weight, weight);
      if(SRSS[ 1]>=0){
	int t = SRSS[ 1];
	if(MjjL()>400.||DetajjL()>1.5) t = -1;
	if(fabs(Mjj()     -80.)<=15.)  t = -1;
	if(SRSS[ 1]==0&&(met_pt()>40.||MllSS()<=40.)              ) t = -1;
	if(SRSS[ 1]==1&&(met_pt()>40.||MllSS()<=30.||MTmax()<=90.)) t = -1;
	if(SRSS[ 1]==2&&(met_pt()>40.||MllSS()<=40.)              ) t = -1;
	t2 = t;
	fillSRhisto(histos, "YieldsSR_Mjjsideband_lowMET",   sample, sn, sn, t, -1, weight, weight,true);
	t = SRSS[ 1];
	if(SRSS[ 1]!=1) t = -1;
	if(MjjL()>400.||DetajjL()>1.5) t = -1;
	if(fabs(Mjj()     -80.)<=15.)  t = -1;
	if(SRSS[ 1]==1&&(met_pt()<=60.||MllSS()<=30.||MTmax()>90.)) t = -1;
	fillSRhisto(histos, "YieldsSR_Mjjsideband_lowMTmax",   sample, sn, sn, t, -1, weight, weight);
      }
      if(t1!=t2){
	cout << "Old low MET Mjj sideband " << t1 << " new " << t2 << " SR " << SRSS[1] << " " << SRSS[21] << " Mjj " << Mjj() << endl;
      }
      
      fillSRhisto(histos, "YieldsCR",           sample, sn,sn, SRSS[ 2], SR3l[ 2], weight, weight);
      fillSRhisto(histos, "YieldsCR_jesup",     sample, sn,sn, SRSS[ 8], SR3l[ 8], weight, weight);
      fillSRhisto(histos, "YieldsCR_jesdn",     sample, sn,sn, SRSS[12], SR3l[12], weight, weight);
      fillSRhisto(histos, "YieldsCR_PUup",      sample, sn,sn, SRSS[ 8], SR3l[ 8], weight*purewgt_up()/purewgt(), weight*purewgt_up()/purewgt());
      fillSRhisto(histos, "YieldsCR_PUdn",      sample, sn,sn, SRSS[12], SR3l[12], weight*purewgt_dn()/purewgt(), weight*purewgt_dn()/purewgt());
      if(fabs(Mjj()   -80.)<20.||SR3l[ 2]>=0) fillSRhisto(histos, "YieldsCR_cutonMjj",       sample, sn,sn, SRSS[ 2], SR3l[ 2], weight, weight);
      if(fabs(Mjj_up()-80.)<20.||SR3l[ 8]>=0) fillSRhisto(histos, "YieldsCR_cutonMjj_jesup", sample, sn,sn, SRSS[ 8], SR3l[ 8], weight, weight);
      if(fabs(Mjj_dn()-80.)<20.||SR3l[12]>=0) fillSRhisto(histos, "YieldsCR_cutonMjj_jesdn", sample, sn,sn, SRSS[12], SR3l[12], weight, weight);
      if(SRSS[ 2]>=0) fillhisto(histos, "Mjj_CRlike_allSS",       sample, sn, Mjj(),    weight);
      if(SRSS[ 8]>=0) fillhisto(histos, "Mjj_CRlike_allSS_jesup", sample, sn, Mjj_up(), weight);
      if(SRSS[12]>=0) fillhisto(histos, "Mjj_CRlike_allSS_jesdn", sample, sn, Mjj_dn(), weight);
      fillSRhisto(histos, "YieldsCR_lepSFup",     sample, sn,sn, SRSS[ 2], SR3l[ 2], weight*lepsf_up()/lepsf(), weight*lepsf_up()/lepsf());
      fillSRhisto(histos, "YieldsCR_lepSFdn",     sample, sn,sn, SRSS[ 2], SR3l[ 2], weight*lepsf_dn()/lepsf(), weight*lepsf_dn()/lepsf());
      if(SRSS[ 4]>=0&&MSFOSvec.size()>0){
	fillSRhisto(histos, "YieldsCR_noMSFOSsel_SS",     sample, sn,sn, SRSS[ 4], -1, weight, weight);
	fillhisto(  histos, "MSFOS_CRlike_allSS",         sample, sn, MSFOSvec[0],     weight);
	fillhisto(  histos, "MSFOS_CRlike_allSS_lepSFup", sample, sn, MSFOSvec[0],     weight*lepsf_up()/lepsf());
	fillhisto(  histos, "MSFOS_CRlike_allSS_lepSFdn", sample, sn, MSFOSvec[0],     weight*lepsf_dn()/lepsf());
	for(unsigned int i = 0; i<MSFOSvec.size();++i){
	  fillhisto(histos, "MSFOSall_CRlike_allSS",         sample, sn, MSFOSvec[i], weight);
	  fillhisto(histos, "MSFOSall_CRlike_allSS_lepSFup", sample, sn, MSFOSvec[i], weight*lepsf_up()/lepsf());
	  fillhisto(histos, "MSFOSall_CRlike_allSS_lepSFdn", sample, sn, MSFOSvec[i], weight*lepsf_dn()/lepsf());
	}
      }
      if((SR3l[ 1]>=0||SR3l[ 3]>=0)&&MSFOSvec.size()>0&&fabs(M3l()-90.)>10.){
	bool passMET = true;
	if((SR3l[ 1]==0||SR3l[ 3]==0)&&met_pt()<30.) passMET = false;
	if((SR3l[ 1]==1||SR3l[ 3]==1)&&met_pt()<40.) passMET = false;
	if((SR3l[ 1]==2||SR3l[ 3]==2)&&met_pt()<55.) passMET = false;
	bool passPTlll             = ((SR3l[ 1]==0||SR3l[ 3]==0) ? (true) : (Pt3l()>=60.));
	bool passDPhilllMET        = (DPhi3lMET()>=2.5);
	bool passMT3rd             = ((SR3l[ 1]==1||SR3l[ 3]==1)&&MT3rd()>90);
	bool passneither           =  !passMET&&!passPTlll&&!passDPhilllMET&& passMT3rd;
	bool passnotall            = !(passMET&& passPTlll&& passDPhilllMET)&&passMT3rd;
	bool passneitherInclMT3rd  =  !passMET&&!passPTlll&&!passDPhilllMET&&!passMT3rd;
	bool passnotallInclMT3rd   = !(passMET&& passPTlll&& passDPhilllMET&& passMT3rd);
	bool invertMT3rd           =   passMET&& passPTlll&& passDPhilllMET&&!passMT3rd;
	
	if(passneither){
	  fillSRhisto(histos, "YieldsSR_invertMETdPhiPt",         sample, sn, sn, -1, SR3l[ 1], weight, weight);
	  fillSRhisto(histos, "YieldsCR_invertMETdPhiPt",         sample, sn, sn, -1, SR3l[ 3], weight, weight);
	  fillhisto(  histos, "MSFOS_all3l_invertMETdPhiPt",      sample, sn, MSFOSvec[0], weight);
	  if((SR3l[ 1]==1)||(SR3l[ 3]==1)){
	    fillhisto(histos, "MSFOS_1SFOS_invertMETdPhiPt",      sample, sn, MSFOSvec[0], weight);
	  }
	  if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	    fillhisto(histos, "MSFOS_2SFOS_invertMETdPhiPt",      sample, sn, MSFOSvec[0], weight);
	  }
	  for(unsigned int i = 0; i<MSFOSvec.size();++i){
	    fillhisto(histos, "MSFOSall_all3l_invertMETdPhiPt",   sample, sn, MSFOSvec[i], weight);
	    if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	      fillhisto(histos, "MSFOSall_2SFOS_invertMETdPhiPt", sample, sn, MSFOSvec[i], weight);
	    }
	  }
	}
	if(passnotall){
	  fillSRhisto(histos, "YieldsSR_inverteitherMETdPhiPt",            sample, sn, sn, -1, SR3l[ 1], weight, weight);
	  fillSRhisto(histos, "YieldsCR_inverteitherMETdPhiPt",            sample, sn, sn, -1, SR3l[ 3], weight, weight);
	  fillhisto(  histos, "MSFOS_all3l_inverteitherMETdPhiPt",         sample, sn, MSFOSvec[0], weight);
	  if((SR3l[ 1]==1)||(SR3l[ 3]==1)){
	    fillhisto(histos, "MSFOS_1SFOS_inverteitherMETdPhiPt",      sample, sn, MSFOSvec[0], weight);
	  }
	  if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	    fillhisto(histos, "MSFOS_2SFOS_inverteitherMETdPhiPt",      sample, sn, MSFOSvec[0], weight);
	  }
	  for(unsigned int i = 0; i<MSFOSvec.size();++i){
	    fillhisto(histos, "MSFOSall_all3l_inverteitherMETdPhiPt",   sample, sn, MSFOSvec[i], weight);
	    if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	      fillhisto(histos, "MSFOSall_2SFOS_inverteitherMETdPhiPt", sample, sn, MSFOSvec[i], weight);
	    }
	  }
	}
	//new start
	if((SR3l[ 1]==1)||(SR3l[ 3]==1)){
	  if(MSFOSvec.size()!=1) cout << "WTF " << __LINE__ << " " << MSFOSvec.size() << endl;
	  if(invertMT3rd){
	    fillSRhisto(histos, "YieldsSR_invertMT3rd",                   sample, sn, sn, -1, SR3l[ 1], weight, weight);
	    fillSRhisto(histos, "YieldsCR_invertMT3rd",                   sample, sn, sn, -1, SR3l[ 3], weight, weight);
	    fillhisto(  histos, "MSFOS_1SFOS_invertMT3rd",                sample, sn, MSFOSvec[0],      weight);
	  }
	  if(passneitherInclMT3rd){
	    fillSRhisto(histos, "YieldsSR_invertMT3rdMETdPhiPt",          sample, sn, sn, -1, SR3l[ 1], weight, weight);
	    fillSRhisto(histos, "YieldsCR_invertMT3rdMETdPhiPt",          sample, sn, sn, -1, SR3l[ 3], weight, weight);
	    fillhisto(  histos, "MSFOS_1SFOS_invertMT3rdMETdPhiPt",       sample, sn, MSFOSvec[0],      weight);
	  }
	  if(passnotallInclMT3rd){
	    fillSRhisto(histos, "YieldsSR_inverteitherMT3rdMETdPhiPt",    sample, sn, sn, -1, SR3l[ 1], weight, weight);
	    fillSRhisto(histos, "YieldsCR_inverteitherMT3rdMETdPhiPt",    sample, sn, sn, -1, SR3l[ 3], weight, weight);
	    fillhisto(  histos, "MSFOS_1SFOS_inverteitherMT3rdMETdPhiPt", sample, sn, MSFOSvec[0],      weight);
	  }
	}
	//new end
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

  SaveHistosToFile("rootfiles/Check3lWZCR_test.root",histos,true,true,(chainnumber==0));
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
