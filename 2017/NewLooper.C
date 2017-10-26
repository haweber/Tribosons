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
  //addeventtocheck(e,276361,248, 475159525);

  bool btagreweighting = false;
  bool applylepSF = false;
  float muptmin = 20.1;                           float muptmax = 199.9; float muetamin =  0.01; float muetamax = 2.49;
  float elptmin = 10.1; float elptminReco = 25.1; float elptmax = 499.9; float eletamin = -2.49; float eletamax = 2.49;
  TFile *fSF;
  TH2F *hMu, *hX, *hElReco, *hElID;
  int SFloaded = -1;
  if(applylepSF) SFloaded = loadlepSFfile(fSF,hMu,hX,hElReco,hElID,"rootfiles/SF_TnP.root","muSF","","elSF_reco","elSF_ID");

  
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  bool storeeventnumbers = false;
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
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();

  histonames.push_back("RawSignalRegion");                     hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawApplicationRegion");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawWZControlRegion");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawSignalRegionPresel");               hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawApplicationRegionPresel");          hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawWZControlRegionPresel");            hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);

  histonames.push_back("SignalRegionPrecleaning");             hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion");                        hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ApplicationRegion");                   hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("WZControlRegion");                     hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ApplicationRegionPrecleaning");        hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("WZControlRegionPrecleaning");          hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegionPresel");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegionPreselPrecleaning");       hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ApplicationRegionPresel");             hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("WZControlRegionPresel");               hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ApplicationRegionPreselPrecleaning");  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("WZControlRegionPreselPrecleaning");    hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);

  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir);
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
      
      bool checkevent = false;
      for(unsigned int i = 0; i<e.size();++i){
	if(e[i].run!=tas::run() ) continue;
	if(e[i].ls !=tas::lumi()) continue;
	if(e[i].evt!=tas::evt() ) continue;
	checkevent = true;
	cout << "Check event " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	break;
      }
      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);
      if(checkevent) cout << "nj30 " << nj30 << " nj " << nj << " nb " << nb << " Mjj " << Mjj << " MjjL " << MjjL << " Detajj " << Detajj << endl;
      if(checkevent){
	for(unsigned int i = 0; i<jets_p4().size(); ++i){
	  cout << "jet pT " << jets_p4()[i].Pt() << " eta " << jets_p4()[i].Eta() << " CSV " << jets_csv()[i];// << endl;
	  for(unsigned int j = i+1; j<jets_p4().size(); ++j) cout << " M"<<i<<j<< " " << (jets_p4()[i]+jets_p4()[j]).M() << " (dR " << dR(jets_p4()[i],jets_p4()[j]) << ")";
	  cout << endl;
	}
      }

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
      if(checkevent) cout << "nSS " << nSS << " n3l " << n3l << " naSS " << naSS << " na3l " << na3l << " nvetoaSS " << nvetoaSS << " nvetoa3l " << nvetoa3l << " ntracks " << nisoTrack_mt2_cleaned_VVV_cutbased_veto() << endl;
      if(checkevent){
	for(unsigned int i = 0; i<lep_pdgId().size();++i){
	  cout << "lep " << lep_pdgId()[i] << " Pt " << lep_p4()[i].Pt() << " eta " << lep_p4()[i].Eta() << " ID t/l/v/trig " << lep_pass_VVV_cutbased_tight_noiso()[i] << "/" << lep_pass_VVV_cutbased_fo_noiso()[i] << "/" << lep_pass_VVV_cutbased_veto_noiso()[i] << "/" << lep_isTriggerSafe_v1()[i] << " iso " << lep_relIso03EAv2()[i] << " ip3d " << lep_ip3d()[i] << " losthits " << lep_lostHits()[i] << " t.q " << lep_tightCharge()[i];
	  for(unsigned int j = i+1; j<lep_pdgId().size();++j) { cout << " M" << i << j << " " << (lep_p4()[i]+lep_p4()[j]).M();
	    for(unsigned int k = j+1; k<lep_pdgId().size();++k) cout << " M" << i << j << k << " " << (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).M() << " Pt " <<  (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).Pt() << " DPhiMET " << dPhi( (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]),MET); }
	  cout << endl;
	}
      }
      
      if((n3l+na3l)<2) continue;
      bool passofflineforTrigger = passofflineTriggers(i3l, ia3l);
      if(!passofflineforTrigger) continue;
      if(checkevent) cout << "pass offline" << endl;
      
      if(isData()){
	if(!passFilters()) continue;
	if(checkevent) cout << "pass filter" << endl;
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) { continue; }
	if(checkevent) cout << "pass duplicate" << endl;
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	if(checkevent) cout << "pass goodrun" << endl;
	bool passonlineTrigger = passonlineTriggers(i3l, ia3l);//currently applied only to data
	if(!passonlineTrigger) continue;
      }
      if(checkevent) cout << "pass online" << endl;

      string sn = skimFilePrefix;
      string sn2 = skimFilePrefix;
      string temp = process(fname,true ,iSS,iaSS);
      string temp2 = process(fname,false,i3l,ia3l);
      bool isphotonSS = (temp =="photonfakes");
      bool isphoton3l = (temp2=="photonfakes");
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
      float MTmax3l = calcMTmax(i3l,MET,true);
      if(checkevent) cout << "MET " << MET.Pt() << " MTmax " << MTmax << " MTmax3l " << MTmax3l << endl;

      int SRSS[6]; bool selects3l[6];
      int SR3l[6];
      for(int i = 0; i<6; ++i) { SRSS[i] = -1; SR3l[i] = -1; selects3l[i] = false; }
      //SS
      //0: SR
      SRSS[0] = isSRSS(iSS,      vSS,false,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //1: SR preselect
      SRSS[1] = isSRSS(iSS,      vSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //2: AR
      SRSS[2] = isARSS(iSS,iaSS,vaSS,false,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //3: AR preselect
      SRSS[3] = isARSS(iSS,iaSS,vaSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //4: CR
      SRSS[4] = isCRSS(iSS,i3l,  v3l,false,MTmax3l,nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //5: CR preselect
      SRSS[5] = isCRSS(iSS,i3l,  v3l,true ,MTmax3l,nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      selects3l[4] = true; selects3l[5] = true;
      //3l
      //0: SR, 4: CR
      checkbothSRCR3l(SR3l[0],SR3l[4],i3l,false,nj,nb);
      //1: SR preselect, 5: CR preselect
      checkbothSRCR3l(SR3l[1],SR3l[5],i3l,true ,nj,nb);
      //2: AR
      SR3l[2] = isAR3l(i3l,ia3l,false,nj,nb);
      //3: AR preselect
      SR3l[3] = isAR3l(i3l,ia3l,true ,nj,nb);

      if(!isData()){//SR is blinded
	fillSRhisto(histos, "SignalRegionPrecleaning",       sn, sn2, SRSS[0], SR3l[0], weight*lepSFSS, weight*lepSF3l);
	fillSRhisto(histos, "SignalRegionPreselPrecleaning", sn, sn2, SRSS[1], SR3l[1], weight*lepSFSS, weight*lepSF3l);
      }
      fillSRhisto(histos, "ApplicationRegionPrecleaning",       sn, sn2, SRSS[2], SR3l[2], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "ApplicationRegionPreselPrecleaning", sn, sn2, SRSS[3], SR3l[3], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "WZControlRegionPrecleaning",         sn2,sn2, SRSS[4], SR3l[4], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "WZControlRegionPreselPrecleaning",   sn2,sn2, SRSS[5], SR3l[5], weight*lepSFSS, weight*lepSF3l);


      if(checkevent) cout << "passed       SRSS " << SRSS[0] << " SR3l " << SR3l[0] << " ARSS " << SRSS[2] << " AR3l " << SR3l[2] << " CRSS " << SRSS[4] << " CR3l " << SR3l[4] << endl;
      //if(SRSS[2]==0) cout << __LINE__ << endl;
      for(int i = 0; i<6; ++i) {
	if(!selects3l[i]){
	  if(vetophotonprocess(fname,isphotonSS))    { SRSS[i] = -1; }
	}
	else if(vetophotonprocess(fname,isphoton3l)){ SRSS[i] = -1; }
	if(vetophotonprocess(fname,isphoton3l))     { SR3l[i] = -1; }
      }
      if(checkevent) cout << "photonpassed SRSS " << SRSS[0] << " SR3l " << SR3l[0] << " ARSS " << SRSS[2] << " AR3l " << SR3l[2] << " CRSS " << SRSS[4] << " CR3l " << SR3l[4] << endl;


      if(!isData()){//SR is blinded
	fillSRhisto(histos, "SignalRegion",       sn, sn2, SRSS[0], SR3l[0], weight*lepSFSS, weight*lepSF3l);
	fillSRhisto(histos, "SignalRegionPresel", sn, sn2, SRSS[1], SR3l[1], weight*lepSFSS, weight*lepSF3l);
      }
      fillSRhisto(histos, "ApplicationRegion",       sn, sn2, SRSS[2], SR3l[2], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "ApplicationRegionPresel", sn, sn2, SRSS[3], SR3l[3], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "WZControlRegion",         sn2,sn2, SRSS[4], SR3l[4], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "WZControlRegionPresel",   sn2,sn2, SRSS[5], SR3l[5], weight*lepSFSS, weight*lepSF3l);

      if(!isData()){//SR is blinded
	fillSRhisto(histos, "RawSignalRegion",       sn, sn2, SRSS[0], SR3l[0], 1., 1.);
	fillSRhisto(histos, "RawSignalRegionPresel", sn, sn2, SRSS[1], SR3l[1], 1., 1.);
      }
      fillSRhisto(histos, "RawApplicationRegion",       sn, sn2, SRSS[2], SR3l[2], 1., 1.);
      fillSRhisto(histos, "RawApplicationRegionPresel", sn, sn2, SRSS[3], SR3l[3], 1., 1.);
      fillSRhisto(histos, "RawWZControlRegion",         sn2,sn2, SRSS[4], SR3l[4], 1., 1.);
      fillSRhisto(histos, "RawWZControlRegionPresel",   sn2,sn2, SRSS[5], SR3l[5], 1., 1.);

      if(storeeventnumbers){
	addeventtolist(SRSS[0], SR3l[0], SREE, SREM, SRMM, SR0SFOS, SR1SFOS, SR2SFOS);
	addeventtolist(SRSS[2], SR3l[2], AREE, AREM, ARMM, AR0SFOS, AR1SFOS, AR2SFOS);
	addeventtolist(SRSS[4], SR3l[4], CREE, CREM, CRMM, CR0SFOS, CR1SFOS, CR2SFOS);
	//addeventtolist(SRSS[1], SR3l[1], SREE, SREM, SRMM, SR0SFOS, SR1SFOS, SR2SFOS);
	//addeventtolist(SRSS[3], SR3l[3], AREE, AREM, ARMM, AR0SFOS, AR1SFOS, AR2SFOS);
	//addeventtolist(SRSS[5], SR3l[5], CREE, CREM, CRMM, CR0SFOS, CR1SFOS, CR2SFOS);
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
  }
  
  SaveHistosToFile("rootfiles/newLooper.root",histos,true,true);
  //SaveHistosToFile("testlogs/testV2.root",histos,true,true);

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
