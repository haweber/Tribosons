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
//#include "CMS3_WWW0117.cc"
#ifdef USE_CMS3_WWW100
#include "CMS3_WWW100.cc"
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

  histonames.push_back("NJ30_SRpreselect_ge0j");            hbins.push_back(10); hlow.push_back(    0); hup.push_back(10);
  histonames.push_back("NJsoft_SRpreselect_ge0j");          hbins.push_back(10); hlow.push_back(    0); hup.push_back(10);
  histonames.push_back("NJsoft_SRpreselect_eq1j");          hbins.push_back(10); hlow.push_back(    0); hup.push_back(10);

  histonames.push_back("NJ30_SR_ge0j");                     hbins.push_back(10); hlow.push_back(    0); hup.push_back(10);
  histonames.push_back("NJsoft_SR_ge0j");                   hbins.push_back(10); hlow.push_back(    0); hup.push_back(10);
  histonames.push_back("NJsoft_SR_eq1j");                   hbins.push_back(10); hlow.push_back(    0); hup.push_back(10);

  histonames.push_back("Mother_NJ30_SR_ge2j");              hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Mother_NJ30_SR_Mjjsideband");       hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Mother_NJ30_SR_Mjjsidebandfull");   hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Mother_NJ30_SR_singlej");           hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Mother_NJ30_SR_eq1j");              hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Mother_NJsoft_SR_eq1j");            hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);

  histonames.push_back("Motherv2_NJ30_SR_ge2j");            hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Motherv2_NJ30_SR_Mjjsideband");     hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Motherv2_NJ30_SR_Mjjsidebandfull"); hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Motherv2_NJ30_SR_singlej");         hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Motherv2_NJ30_SR_eq1j");            hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);
  histonames.push_back("Motherv2_NJsoft_SR_eq1j");          hbins.push_back(30); hlow.push_back(    0); hup.push_back(30);

  histonames.push_back("SignalRegionsAllSS");               hbins.push_back(18); hlow.push_back(    0); hup.push_back(18);

  int nhistos = histonames.size();
  for(int i = 0; i<nhistos; ++i){
    for(int j = 0; j<3; ++j){
      string x;
      if(j==0) x = "_SSee";
      if(j==1) x = "_SSem";
      if(j==2) x = "_SSmm";
      histonames.push_back(histonames[i]+x);  hbins.push_back(hbins[i]); hlow.push_back(hlow[i]); hup.push_back(hup[i]);
    }
  }
  
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

      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);
      int njsoft = 0;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5) continue;
	if(   jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=3.0) ++njsoft;
      }

      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      //getleptonindices_BDT(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l, 1, 0.85);
      //getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l,1,25,20);
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
      
      int SRSS[20]; bool selects3l[20];
      int SR3l[20];
      for(int i = 0; i<20; ++i) { SRSS[i] = -1; SR3l[i] = -1; selects3l[i] = false; }
      
      if(i3l.size()==3){
	float pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0); bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0); bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0); bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	int SFOScounter = 0;
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	if(nj <=1)        SFOScounter =  -1;
	if(nb!=0)         SFOScounter =  -1;
	bool pass0(false), pass1(false), pass2(false);   //SR
	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-MZ)<10.) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) { pass0 = false; }
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) { pass0 = false; }
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) { pass0 = false; }
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-MZ)<15.) { pass0 = false; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-MZ)<15.) { pass0 = false; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-MZ)<15.) { pass0 = false; }
	  if(MET.Pt()<30) pass0 = false;
	}
	if(SFOScounter==1){
	  pass1 = true;
	  float metcut = 40.;
	  int thirdlep = -1;
	  if(OS01&&SF01) thirdlep = i3l[2];
	  if(OS02&&SF02) thirdlep = i3l[1];
	  if(OS12&&SF12) thirdlep = i3l[0];
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-MZ)<20.) { pass1 = false; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-MZ)<20.) { pass1 = false; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-MZ)<20.) { pass1 = false; }
	  if(thirdlep>=0&&mT(lep_p4()[thirdlep],MET)<90.) pass1 = false;
	  if(pTlll<60.)  pass1 = false;
	  if(MET.Pt()<40)  pass1 = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-MZ)<20.) { pass2 = false; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-MZ)<20.) { pass2 = false; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-MZ)<20.) { pass2 = false; }
	  if(pTlll<60.) pass2 = false;
	  if(MET.Pt()<55) pass2 = false;
	}
	if(pass0) SR3l[12] =  0;
	if(pass1) SR3l[12] =  1;
	if(pass2) SR3l[12] =  2;
      }

	  

      float minDRlj(999.);
      float minMT(999999.);
      for(unsigned int j1 = 0; j1<jets_p4().size();++j1){
	if(fabs(jets_p4()[j1].Eta())>2.4) continue;
	if(     jets_p4()[j1].Pt()  <30.) continue;
	for(unsigned int i = 0; i<iSS.size(); ++i){
	  if(minDRlj>dR(lep_p4()[iSS[i] ],jets_p4()[j1])) minDRlj = dR(lep_p4()[iSS[i] ],jets_p4()[j1]);
	}
      }
      for(unsigned int i = 0; i<iSS.size(); ++i){
	if(minMT  >mT(lep_p4()[iSS[i] ],MET)) minMT   = mT(lep_p4()[iSS[i] ],MET);
      }
      bool passSSsoftjets = true;
      float Massjj  = -9999;
      float MassjjL = -9999;
      float Deta    = -9999;
      float minDR=999.;
      int jDR1(-1), jDR2(-1);
      float minDRljsoft(999.);
      for(unsigned int j1 = 0; j1<jets_p4().size();++j1){
	if(fabs(jets_p4()[j1].Eta())>3.0) continue;
	if(     jets_p4()[j1].Pt()  <20.) continue;
	for(unsigned int i = 0; i<iSS.size(); ++i){
	  if(minDRljsoft>dR(lep_p4()[iSS[i] ],jets_p4()[j1])) minDRljsoft = dR(lep_p4()[iSS[i] ],jets_p4()[j1]);
	}
	for(unsigned int j2 = j1+1; j2<jets_p4().size();++j2){
	  if(fabs(jets_p4()[j2].Eta())>3.0) continue;
	  if(     jets_p4()[j2].Pt()  <20.) continue;
	  if(MassjjL<0) {
	    MassjjL    =     (jets_p4()[j1]+jets_p4()[j2]).M();
	    Deta       = dEta(jets_p4()[j1],jets_p4()[j2]);
	  }
	  if(       dR(jets_p4()[j1], jets_p4()[j2])<minDR){
	    minDR = dR(jets_p4()[j1], jets_p4()[j2]);
	    jDR1 = j1; jDR2 = j2;
	  }
	}
      }
      if(njsoft<1) passSSsoftjets = false;
      if(nj30!=1)  passSSsoftjets = false;
      if(nb!=0)    passSSsoftjets = false;
      if(jDR1>=0&&jDR2>=0) Massjj = (jets_p4()[jDR1]+jets_p4()[jDR2]).M();
      if(fabs(Deta)>1.5)       passSSsoftjets = false;
      if(fabs(MassjjL)>400.)   passSSsoftjets = false;
      //if(fabs(Massjj-80.)>=15) passSSsoftjets = false;
      if(fabs(Massjj)>=120)    passSSsoftjets = false;
      bool passSS = true;
      if(!passSSsoftjets) passSS = false;
      if(iSS.size()!=2)   passSS = false;
      if(vSS.size()!=0)   passSS = false;
      if(minDRljsoft>1.5) passSS = false;
      if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)          passSS = false;
      if(passSS&&(lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])<0) passSS = false;
      bool ee = false; bool em = false; bool mm = false;
      if(passSS&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==121&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-MZ)>10.) ee = true;
      if(passSS&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==143)  em = true;
      if(passSS&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==169)  mm = true;

      if(passSS&&ee&&MET.Pt()>60.&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.           ) SRSS[12] = 0;
      if(passSS&&em&&MET.Pt()>60.&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.) SRSS[12] = 1;
      if(passSS&&mm&&              (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.           ) SRSS[12] = 2;

      passSS = true;
      if(nj30!=1)         passSS = false;
      if(nb!=0)           passSS = false;
      if(iSS.size()!=2)   passSS = false;
      if(vSS.size()!=0)   passSS = false;
      if(minDRlj>1.5)     passSS = false;
      if(minMT<=50.||minMT>999990.)                             passSS = false;
      if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)          passSS = false;
      if(passSS&&(lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])<0) passSS = false;
      ee = false; em = false; mm = false;
      if(passSS&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==121&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-MZ)>10.) ee = true;
      if(passSS&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==143)  em = true;
      if(passSS&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==169)  mm = true;

      if(passSS&&ee&&MET.Pt()>60.&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.           ) SRSS[15] = 0;
      if(passSS&&em&&MET.Pt()>60.&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.) SRSS[15] = 1;
      if(passSS&&mm&&              (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.           ) SRSS[15] = 2;

    
      //1: SR preselect
      SRSS[1] = isSRSS(iSS,      vSS,true ,MTmax,  2,nb,Mjj,MjjL,Detajj);//enter variables for quicker calculation
      //11: SR preselect
      SRSS[10] = isSRSS(iSS,      vSS,false,MTmax,nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      SRSS[11] = isSRSS(iSS,      vSS,true ,MTmax,   2,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation

      SRSS[13] = isSRSS(iSS,      vSS,false,MTmax,nj30,nb,Mjj,MjjL,Detajj,MET,0,false,true,1);//enter variables for quicker calculation
      SRSS[14] = SRSS[13];
      if(Mjj>120) SRSS[13] = -1;

      SR3l[10] = isSR3l(i3l,false,nj,nb,MET,0,false,1);

      for(int i = 0; i<20; ++i) {
	if(!selects3l[i]){
	  if(vetophotonprocess(fname,isphotonSS))    { SRSS[i] = -1; }
	}
	else if(vetophotonprocess(fname,isphoton3l)){ SRSS[i] = -1; }
	if(vetophotonprocess(fname,isphoton3l))     { SR3l[i] = -1; }
      }
      vector<int> m2j,   m1j,   m2jhard,   m2jsoft,   mSB,   mSBfull;
      vector<int> m2jv2, m1jv2, m2jhardv2, m2jsoftv2, mSBv2, mSBfullv2;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	double minDeltaR (999.); int currentidx(-1);
	bool matchW = false; int Wmatchidx(-1);
	for(unsigned int i = 0; i<genPart_pdgId().size();++i){
	  if(abs(genPart_pdgId()[i])>5) continue;
	  if(dR(genPart_p4()[i],jets_p4()[n])>0.4) continue;
	  if(    genPart_status()[i]!=23) continue;
	  if(abs(genPart_motherId()[i])==24 || abs(genPart_motherId()[i])==25) { matchW = true; Wmatchidx = i; }
	  if(minDeltaR>dR(genPart_p4()[i],jets_p4()[n])){ minDeltaR = dR(genPart_p4()[i],jets_p4()[n]); currentidx = i; }
	  //cout << "jet " << jets_p4()[n].Pt() << " gen " << genPart_p4()[i].Pt() << " dR " << dR(genPart_p4()[i],jets_p4()[n]) << " genid " << genPart_pdgId()[i] << " genstatus " << genPart_status()[i] << " mother " << genPart_motherId()[i] << endl;
	}
	if(matchW == false) Wmatchidx = currentidx;
	if(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5) {
	  if(SRSS[10]>=0&&currentidx>=0) m2j      .push_back(abs(genPart_motherId()[currentidx]));
	  if(SRSS[10]>=0&&Wmatchidx >=0) m2jv2    .push_back(abs(genPart_motherId()[Wmatchidx ]));
	  if(SRSS[12]>=0&&currentidx>=0) m2jhard  .push_back(abs(genPart_motherId()[currentidx]));
	  if(SRSS[12]>=0&&Wmatchidx >=0) m2jhardv2.push_back(abs(genPart_motherId()[Wmatchidx ]));
	  if(SRSS[15]>=0&&currentidx>=0) m1j      .push_back(abs(genPart_motherId()[currentidx]));
	  if(SRSS[15]>=0&&Wmatchidx >=0) m1jv2    .push_back(abs(genPart_motherId()[Wmatchidx ]));
	  if(SRSS[13]>=0&&currentidx>=0) mSB      .push_back(abs(genPart_motherId()[currentidx]));
	  if(SRSS[13]>=0&&Wmatchidx >=0) mSBv2    .push_back(abs(genPart_motherId()[Wmatchidx ]));
	  if(SRSS[14]>=0&&currentidx>=0) mSBfull  .push_back(abs(genPart_motherId()[currentidx]));
	  if(SRSS[14]>=0&&Wmatchidx >=0) mSBfullv2.push_back(abs(genPart_motherId()[Wmatchidx ]));
	  continue;
	}
	if(   jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=3.0) {
	  if(SRSS[12]>=0&&currentidx>=0) m2jsoft  .push_back(abs(genPart_motherId()[currentidx]));
	  if(SRSS[12]>=0&&Wmatchidx >=0) m2jsoftv2.push_back(abs(genPart_motherId()[Wmatchidx ]));
	}
      }
      if(!isData()||!blindSR){//SR is blinded
	if(weight<80.){
	  if(SRSS[10]>=0&nj30>=2)  fillhisto(histos, "SignalRegionsAllSS",     sample, sn,  SRSS[11],      weight);
	  if(SRSS[12]>=0&&nj30==1) fillhisto(histos, "SignalRegionsAllSS",     sample, sn,  SRSS[12]+3,    weight);
	  if(SRSS[15]>=0&&nj30==1) fillhisto(histos, "SignalRegionsAllSS",     sample, sn,  SRSS[15]+6,    weight);
	  if(SRSS[13]>=0&&nj30>=2) fillhisto(histos, "SignalRegionsAllSS",     sample, sn,  SRSS[13]+9,    weight);
	  if(SR3l[10]>=0         ) fillhisto(histos, "SignalRegionsAllSS",     sample, sn2, SR3l[10]+12,    weight);
	  if(SR3l[12]>=0         ) fillhisto(histos, "SignalRegionsAllSS",     sample, sn2, SR3l[12]+15,    weight);
	}
	//fillSRhisto(histos, "SignalRegion",               sample, sn, sn2, SRSS[0], SR3l[0], weight, weight);
	if(SRSS[11]>=0)          fillhisto(histos, "NJ30_SRpreselect_ge0j",  sample, sn, nj30,      weight);
	if(SRSS[11]>=0)          fillhisto(histos, "NJsoft_SRpreselect_ge0j",sample, sn, njsoft,    weight);
	if(SRSS[11]>=0&&nj30==1) fillhisto(histos, "NJsoft_SRpreselect_eq1j",sample, sn, njsoft,    weight);
	if(SRSS[10]>=0&&nj30>=2) fillhisto(histos, "NJ30_SR_ge0j",     sample, sn, nj30,      weight);
	if(SRSS[15]>=0&&nj30==1) fillhisto(histos, "NJ30_SR_ge0j",     sample, sn, nj30,      weight);
	if(SRSS[12]>=0)          fillhisto(histos, "NJsoft_SR_ge0j",   sample, sn, njsoft,    weight);
	if(SRSS[12]>=0&&nj30==1) fillhisto(histos, "NJsoft_SR_eq1j",   sample, sn, njsoft,    weight);
	if(SRSS[10]>=0){
	  for(unsigned int i = 0; i<m2j  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_ge2j",   sample, sn, m2j[i],    weight);
	  for(unsigned int i = 0; i<m2jv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_ge2j", sample, sn, m2jv2[i],  weight);
	}
	if(SRSS[12]>=0&&nj30==1){
	  for(unsigned int i = 0; i<m2jhard  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_eq1j",     sample, sn, m2jhard[i],    weight);
	  for(unsigned int i = 0; i<m2jhardv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_eq1j",   sample, sn, m2jhardv2[i],  weight);
	  for(unsigned int i = 0; i<m2jsoft  .size(); ++i) fillhisto(histos, "Mother_NJsoft_SR_eq1j",   sample, sn, m2jsoft[i],    weight);
	  for(unsigned int i = 0; i<m2jsoftv2.size(); ++i) fillhisto(histos, "Motherv2_NJsoft_SR_eq1j", sample, sn, m2jsoftv2[i],  weight);
	}
	if(SRSS[15]>=0&&nj30==1){
	  for(unsigned int i = 0; i<m1j  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_singlej",   sample, sn, m1j[i],    weight);
	  for(unsigned int i = 0; i<m1jv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_singlej", sample, sn, m1jv2[i],  weight);
	}
	if(SRSS[13]>=0){
	  for(unsigned int i = 0; i<mSB  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_Mjjsideband",   sample, sn, mSB[i],    weight);
	  for(unsigned int i = 0; i<mSBv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_Mjjsideband", sample, sn, mSBv2[i],  weight);
	}
	if(SRSS[14]>=0){
	  for(unsigned int i = 0; i<mSBfull  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_Mjjsidebandfull",   sample, sn, mSBfull[i],    weight);
	  for(unsigned int i = 0; i<mSBfullv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_Mjjsidebandfull", sample, sn, mSBfullv2[i],  weight);
	}
	string x = "";
	if(iSS.size()>=2&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==121&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-MZ)>10.) x = "_SSee";
	if(iSS.size()>=2&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==143)  x = "_SSem";
	if(iSS.size()>=2&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iSS[1] ]))==169)  x = "_SSmm";
	if(x!=""){
	  if(SRSS[11]>=0)          fillhisto(histos, "NJ30_SRpreselect_ge0j"+x,     sample, sn, nj30,      weight);
	  if(SRSS[11]>=0)          fillhisto(histos, "NJsoft_SRpreselect_ge0j"+x,   sample, sn, njsoft,    weight);
	  if(SRSS[11]>=0&&nj30==1) fillhisto(histos, "NJsoft_SRpreselect_eq1j"+x,   sample, sn, njsoft,    weight);
	  if(SRSS[10]>=0&&nj30>=2) fillhisto(histos, "NJ30_SR_ge0j"+x,     sample, sn, nj30,      weight);
	  if(SRSS[15]>=0&&nj30<=1) fillhisto(histos, "NJ30_SR_ge0j"+x,     sample, sn, nj30,      weight);
	  if(SRSS[12]>=0)          fillhisto(histos, "NJsoft_SR_ge0j"+x,   sample, sn, njsoft,    weight);
	  if(SRSS[12]>=0&&nj30==1) fillhisto(histos, "NJsoft_SR_eq1j"+x,   sample, sn, njsoft,    weight);
	  if(SRSS[10]>=0){
	    for(unsigned int i = 0; i<m2j  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_ge2j"+x,   sample, sn, m2j[i],    weight);
	    for(unsigned int i = 0; i<m2jv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_ge2j"+x, sample, sn, m2jv2[i],  weight);
	  }
	  if(SRSS[12]>=0){
	    for(unsigned int i = 0; i<m2jhard  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_eq1j"+x,     sample, sn, m2jhard[i],    weight);
	    for(unsigned int i = 0; i<m2jhardv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_eq1j"+x,   sample, sn, m2jhardv2[i],  weight);
	    for(unsigned int i = 0; i<m2jsoft  .size(); ++i) fillhisto(histos, "Mother_NJsoft_SR_eq1j"+x,   sample, sn, m2jsoft[i],    weight);
	    for(unsigned int i = 0; i<m2jsoftv2.size(); ++i) fillhisto(histos, "Motherv2_NJsoft_SR_eq1j"+x, sample, sn, m2jsoftv2[i],  weight);
	  }
	  if(SRSS[15]>=0){
	    for(unsigned int i = 0; i<m1j  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_singlej"+x,   sample, sn, m1j[i],    weight);
	    for(unsigned int i = 0; i<m1jv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_singlej"+x, sample, sn, m1jv2[i],  weight);
	  }
	  if(SRSS[13]>=0){
	    for(unsigned int i = 0; i<mSB  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_Mjjsideband"+x,   sample, sn, mSB[i],    weight);
	    for(unsigned int i = 0; i<mSBv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_Mjjsideband"+x, sample, sn, mSBv2[i],  weight);
	  }
	  if(SRSS[14]>=0){
	    for(unsigned int i = 0; i<mSBfull  .size(); ++i) fillhisto(histos, "Mother_NJ30_SR_Mjjsidebandfull"+x,   sample, sn, mSBfull[i],    weight);
	    for(unsigned int i = 0; i<mSBfullv2.size(); ++i) fillhisto(histos, "Motherv2_NJ30_SR_Mjjsidebandfull"+x, sample, sn, mSBfullv2[i],  weight);
	  }
	}//x

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

  SaveHistosToFile("rootfiles/NJetsForSSSR.root",histos,true,true,(chainnumber==0));
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
