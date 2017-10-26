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

  histonames.push_back("Mjj_SRSS_presel");             hbins.push_back(15);  hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("Mjj_SRSS_NmO");                hbins.push_back(15);  hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("Mll_SRSS_presel");             hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);
  histonames.push_back("Mll_SRSS_NmO");                hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);
  histonames.push_back("pTll_SRSS_presel");            hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);
  histonames.push_back("pTll_SRSS_NmO");               hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);
  histonames.push_back("MET_SRSS_presel");             hbins.push_back(15);  hlow.push_back( 0); hup.push_back(150);
  histonames.push_back("MET_SRSS_NmO");                hbins.push_back(15);  hlow.push_back( 0); hup.push_back(150);
  histonames.push_back("MTmax_SRSS_presel");           hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);
  histonames.push_back("MTmax_SRSS_NmO");              hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);
  histonames.push_back("MjjL_SRSS_presel");            hbins.push_back(10);  hlow.push_back( 0); hup.push_back(500);
  histonames.push_back("MjjL_SRSS_NmO");               hbins.push_back(10);  hlow.push_back( 0); hup.push_back(500);
  histonames.push_back("DetajjL_SRSS_presel");         hbins.push_back(12);  hlow.push_back( 0); hup.push_back(3.0);
  histonames.push_back("DetajjL_SRSS_NmO");            hbins.push_back(12);  hlow.push_back( 0); hup.push_back(3.0);

  histonames.push_back("DPhilllMET_SR3l_presel");      hbins.push_back(16);  hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("DPhilllMET_SR3l_NmO");         hbins.push_back(16);  hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("PTlll_SR3l_presel");           hbins.push_back(15);  hlow.push_back( 0); hup.push_back(150);
  histonames.push_back("PTlll_SR3l_NmO");              hbins.push_back(15);  hlow.push_back( 0); hup.push_back(150);
  histonames.push_back("Mlll_SR3l_presel");            hbins.push_back(10);  hlow.push_back(50); hup.push_back(150);
  histonames.push_back("Mlll_SR3l_NmO");               hbins.push_back(10);  hlow.push_back(50); hup.push_back(150);
  histonames.push_back("MET_SR3l_presel");             hbins.push_back(15);  hlow.push_back( 0); hup.push_back(150);
  histonames.push_back("MET_SR3l_NmO");                hbins.push_back(15);  hlow.push_back( 0); hup.push_back(150);
  histonames.push_back("MSFOSZlike_SR3l_presel");      hbins.push_back(16);  hlow.push_back(20); hup.push_back(180);
  histonames.push_back("MSFOSZlike_SR3l_NmO");         hbins.push_back(16);  hlow.push_back(20); hup.push_back(180);
  
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

      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);

      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l);
      float lepSFSS(1.), lepSFerrSS(0.), lepSF3l(1.), lepSFerr3l(0.);
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
	if( is_duplicate(id)        ) { continue; }
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

      int SRSSpresel = -1;
      int SR3lpresel = -1;
      int CR3lpresel = -1;
      //SS
      SRSSpresel = isSRSS(iSS,      vSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //3l
      SR3lpresel = isSR3l(i3l,true ,nj,nb);
      CR3lpresel = isCR3l(i3l,true ,nj,nb);

      if(vetophotonprocess(fname,isphotonSS))    { SRSSpresel = -1; }
      if(vetophotonprocess(fname,isphoton3l))    { SR3lpresel = -1; }

      if(SRSSpresel>=0){
	float Mll  = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();
	float pTll = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt();
	histos["Mjj_SRSS_presel_"+     sn]->Fill(Mjj,      weight*lepSFSS);
	histos["Mll_SRSS_presel_"+     sn]->Fill(Mll,      weight*lepSFSS);
	histos["pTll_SRSS_presel_"+    sn]->Fill(pTll,     weight*lepSFSS);
	histos["MET_SRSS_presel_"+     sn]->Fill(met_pt(), weight*lepSFSS);
	histos["MjjL_SRSS_presel_"+    sn]->Fill(MjjL,     weight*lepSFSS);
	histos["DetajjL_SRSS_presel_"+ sn]->Fill(Detajj,   weight*lepSFSS);
	if(SRSSpresel==1){
	  histos["MTmax_SRSS_presel_"+ sn]->Fill(MTmax,    weight*lepSFSS);
	}
	bool passMjj  = fabs(Mjj-80.)<20.;
	bool passMjjL = MjjL<400.;
	bool passDeta = Detajj<1.5;
	bool passMll  = ((SRSSpresel==1) ? ((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.) : ((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.));
	bool passMET  = ((SRSSpresel==2) ? (true) : (met_pt()>40.));
	bool passMT   = ((SRSSpresel==1) ? (MTmax>90.) : (true));
	if(         passMjjL&&passDeta&&passMll&&passMET&&passMT) histos["Mjj_SRSS_NmO_"+     sn]->Fill(Mjj,      weight*lepSFSS);
	if(passMjj&&passMjjL&&passDeta         &&passMET&&passMT) histos["Mll_SRSS_NmO_"+     sn]->Fill(Mll,      weight*lepSFSS);
	if(passMjj&&passMjjL&&passDeta&&passMll&&passMET&&passMT) histos["pTll_SRSS_NmO_"+    sn]->Fill(pTll,     weight*lepSFSS);
	if(passMjj&&passMjjL&&passDeta&&passMll         &&passMT) histos["MET_SRSS_NmO_"+     sn]->Fill(met_pt(), weight*lepSFSS);
	if(passMjj          &&passDeta&&passMll&&passMET&&passMT) histos["MjjL_SRSS_NmO_"+    sn]->Fill(MjjL,     weight*lepSFSS);
	if(passMjj&&passMjjL          &&passMll&&passMET&&passMT) histos["DetajjL_SRSS_NmO_"+ sn]->Fill(Detajj,   weight*lepSFSS);
	if(SRSSpresel==1){
	  if(passMjj&&passMjjL&&passDeta&&passMll&&passMET      ) histos["MTmax_SRSS_NmO_"+   sn]->Fill(MTmax,      weight*lepSFSS);
	}
      }
      if(SR3lpresel>=0){
	bool passZero  = true;
	if(SR3lpresel==0){
	  bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	  bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	  bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passZero = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) passZero = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) passZero = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-MZ)<15.) passZero = false;
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-MZ)<15.) passZero = false;
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-MZ)<15.) passZero = false;
	}
	if(passZero){
	  float DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	  float pTlll      = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	  float Mlll       = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M();
	  histos["DPhilllMET_SR3l_presel_"+   sn2]->Fill(DPhilllMET, weight*lepSF3l);
	  histos["PTlll_SR3l_presel_"+        sn2]->Fill(pTlll,      weight*lepSF3l);
	  histos["Mlll_SR3l_presel_"+         sn2]->Fill(Mlll,       weight*lepSF3l);
	  histos["MET_SR3l_presel_"+          sn2]->Fill(met_pt(),   weight*lepSF3l);
	   
	  bool passMET   = true; if(SR3lpresel==1) passMET = (met_pt()>45.); if(SR3lpresel==2) passMET = (met_pt()>55.);
	  bool passDPhi  = ((SR3lpresel==0) ? (DPhilllMET>2.7) : (DPhilllMET>2.5));
	  bool passpTlll = (pTlll>60.);
	  bool passMlll  = ((SR3lpresel>=1) ? (fabs(Mlll-MZ)>10.) : (true));
	  if(passMET          &&passpTlll&&passMlll) histos["DPhilllMET_SR3l_NmO_"+   sn2]->Fill(DPhilllMET, weight*lepSF3l);
	  if(passMET&&passDPhi           &&passMlll) histos["PTlll_SR3l_NmO_"+        sn2]->Fill(pTlll,      weight*lepSF3l);
	  if(passMET&&passDPhi&&passpTlll          ) histos["Mlll_SR3l_NmO_"+         sn2]->Fill(Mlll,       weight*lepSF3l);
	  if(         passDPhi&&passpTlll&&passMlll) histos["MET_SR3l_NmO_"+          sn2]->Fill(met_pt(),   weight*lepSF3l);
	}
      }
      if(SR3lpresel>=1||CR3lpresel>=1){
	vector<float> MSFOSvec = allMSFOS(i3l);
	float MSFOSZ = MSFOSvec[0];
	float DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	float pTlll      = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	float Mlll       = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M();
	bool passMET   = true; if(SR3lpresel==1) passMET = (met_pt()>45.); if(SR3lpresel==2) passMET = (met_pt()>55.);
	bool passDPhi  = ((SR3lpresel==0) ? (DPhilllMET>2.7) : (DPhilllMET>2.5));
	bool passpTlll = (pTlll>60.);
	bool passMlll  = ((SR3lpresel>=1) ? (fabs(Mlll-MZ)>10.) : (true));
	histos["MSFOSZlike_SR3l_presel_"+ sn2]->Fill(MSFOSZ,     weight*lepSF3l);
	if(passMET&&passDPhi&&passpTlll&&passMlll) histos["MSFOSZlike_SR3l_NmO_"+ sn2]->Fill(MSFOSZ,     weight*lepSF3l);
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
  
  SaveHistosToFile("rootfiles/newNminusOnePlots.root",histos,true,true);

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
