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
#include "MyFunctions.h"
#include "CMS3_WWW0116.cc" //--> is this loaded by MyFunctions.C ?
#include "/home/users/haweber/CORE/Tools/dorky/dorky.h"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.cc"
#include "/home/users/haweber/CORE/Tools/goodrun.h"
#include "/home/users/haweber/CORE/Tools/goodrun.cc"
//MT2 variants

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  vector<myevt> e;
  addeventtocheck(e,1,36545,5920176);

  bool blindSR = true;
  bool btagreweighting = true;
  bool applylepSF = true;
  float muptmin = 20.1;                           float muptmax = 199.9; float muetamin =  0.01; float muetamax = 2.49;
  float elptmin = 10.1; float elptminReco = 25.1; float elptmax = 499.9; float eletamin = -2.49; float eletamax = 2.49;
  TFile *fSF;
  TH2F *hMu, *hX, *hElReco, *hElID;
  int SFloaded = -1;
  if(applylepSF) SFloaded = loadlepSFfile(fSF,hMu,hX,hElReco,hElID,"rootfiles/SF_TnP.root","muSF","","elSF_reco","elSF_ID");

  
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();

  histonames.push_back("SRyield");             hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);

  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir, true);
  cout << "Loaded histograms" << endl;

  const char* json_file = "/home/users/haweber/StopAnalysisMoriond2017/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);

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
      
      // bool checkevent = false;
      //for(unsigned int i = 0; i<e.size();++i){
      //if(e[i].run!=tas::run() ) continue;
      //if(e[i].ls !=tas::lumi()) continue;
      //if(e[i].evt!=tas::evt() ) continue;
      //checkevent = true;
      //cout << "Check event " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      //break;
      //}
      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);

      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l);
      float lepSFSS(1.), lepSFerrSS(1.), lepSF3l(0.), lepSFerr3l(0.);
      if(applylepSF){
	lepSFSS = getlepSFWeightandError(lepSFerrSS,iSS,iaSS,hMu,hX,hElReco,hElID, muptmin,muptmax,muetamin,muetamax, muptmin,muptmax,muetamin,muetamax, elptminReco,elptmax,eletamin,eletamax, elptmin,elptmax,eletamin,eletamax);
	lepSF3l = getlepSFWeightandError(lepSFerr3l,i3l,ia3l,hMu,hX,hElReco,hElID, muptmin,muptmax,muetamin,muetamax, muptmin,muptmax,muetamin,muetamax, elptminReco,elptmax,eletamin,eletamax, elptmin,elptmax,eletamin,eletamax);
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
	if( is_duplicate(id)        ) { /*cout << "Event " << tas::run() << ":" << tas::lumi() << ":" << tas::run() << " is duplicated." << endl;*/ continue; }
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = passonlineTriggers(i3l, ia3l);//currently applied only to data
	if(!passonlineTrigger) continue;
      }

      string sn = skimFilePrefix;
      string sn2 = skimFilePrefix;
      string temp = process(fname,true ,iSS,iaSS);
      string temp2 = process(fname,false,i3l,ia3l);
      bool isphotonSS = (temp =="photonfakes");
      bool isphoton3l = (temp2=="photonfakes");
      bool trueSS = (temp =="trueSS");
      bool true3l = (temp2=="true3l");
      if(sn.find("Other")!=string::npos&&splitVH(fname)){ sn = "WHtoWWW"; sn2 = sn; }
      else if(sn.find("Background")!=string::npos&&splitVH(fname)) continue;
      else if(sn.find("Background")!=string::npos&&!splitVH(fname)){
	sn  = temp;
	sn2 = temp2;
      }


      float MTmax = -1;
      if(iSS.size()==2) MTmax = calcMTmax(iSS,MET);
      else if(iSS.size()==1&&iaSS.size()>=1){
	vector<int> temp; temp.push_back(iSS[0]); temp.push_back(iaSS[0]);
	MTmax = calcMTmax(temp,MET);
      }
      float MTmax3l = calcMTmax(i3l,MET);

      int SRSS[6]; bool selects3l[6];
      int SR3l[6];
      for(int i = 0; i<6; ++i) { SRSS[i] = -1; SR3l[i] = -1; selects3l[i] = false; }
      //SS
      //0: SR
      SRSS[0] = isSRSS(iSS,      vSS,false,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //1: SR preselect
      //SRSS[1] = isSRSS(iSS,      vSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //2: AR
      //SRSS[2] = isARSS(iSS,iaSS,vaSS,false,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //3: AR preselect
      //SRSS[3] = isARSS(iSS,iaSS,vaSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //4: CR
      //SRSS[4] = isCRSS(iSS,i3l,  v3l,false,MTmax3l,nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //5: CR preselect
      //SRSS[5] = isCRSS(iSS,i3l,  v3l,true ,MTmax3l,nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      selects3l[4] = true; selects3l[5] = true;
      //3l
      //0: SR, 4: CR
      checkbothSRCR3l(SR3l[0],SR3l[4],i3l,false,nj,nb);
      //1: SR preselect, 5: CR preselect
      checkbothSRCR3l(SR3l[1],SR3l[5],i3l,true ,nj,nb);
      //2: AR
      //SR3l[2] = isAR3l(i3l,ia3l,false,nj,nb);
      //3: AR preselect
      //SR3l[3] = isAR3l(i3l,ia3l,true ,nj,nb);

      //if(SRSS[2]==0) cout << __LINE__ << endl;
      for(int i = 0; i<6; ++i) {
	if(!selects3l[i]){
	  if(vetophotonprocess(fname,isphotonSS))    { SRSS[i] = -1; }
	}
	else if(vetophotonprocess(fname,isphoton3l)){ SRSS[i] = -1; }
	if(vetophotonprocess(fname,isphoton3l))     { SR3l[i] = -1; }
      }

      bool anySS = (iSS.size()>=2);
      bool any3l = (i3l.size()>=3);
      if(     anySS&(skimFilePrefix=="WW")&&(fname.find("wpwpjj")!=string::npos)  &&trueSS){ sn = "WWVBS";     sn2 = "WWRest"; }
      else if(anySS&(skimFilePrefix=="WW")&&(fname.find("dbl_scat")!=string::npos)&&trueSS){ sn = "WWDPS";     sn2 = "WWRest"; }
      else if(anySS&(skimFilePrefix=="WW")                                        &&trueSS){ sn = "otherWWSS"; sn2 = "WWRest"; }
      else if(       skimFilePrefix=="WW")                                                 { sn = "WWRest";    sn2 = "WWRest"; }
      if(     anySS&&(skimFilePrefix=="ttV")&&trueSS){ sn  = "ttVSS";                 sn2 = "OtherttV"; }
      else if(any3l&&(skimFilePrefix=="ttV")&&true3l){ sn2 = "ttV3l"; if(sn!="ttVSS") sn  = "OtherttV"; }
      else if(       (skimFilePrefix=="ttV")        ){ sn  = "OtherttV";              sn2 = "OtherttV"; }


      if(!(blindSR&&isData())){
	fillSRhisto(histos, "SRyield",       sn, sn2, SRSS[0], SR3l[0], weight*lepSFSS, weight*lepSF3l);
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

  //SaveHistosToFile("testlogs/newLooper.root",histos,true,true);
  SaveHistosToFile("rootfiles/NewWWTTWSplitted.root",histos,true,true);

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
