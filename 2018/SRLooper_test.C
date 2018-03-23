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

#include "Functions_test.h"
//#include "CMS3_WWW0117.cc"
#ifdef USE_CMS3_WWW100
#include "CMS3.cc"
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

  vector<myevt> e;
  addeventtocheck(e,1, 2842, 1443084);
  addeventtocheck(e,1,17341, 8807395);
  addeventtocheck(e,1,15849, 8049420);
  addeventtocheck(e,1,21682,11012507);

  bool blindSR = false;
  bool btagreweighting = true;
  bool applylepSF      = true;
  bool applytrigSF     = true;
  bool applyPUrewgt    = true;
  bool getJECunc       = true;
  
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
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();

  histonames.push_back("SignalRegion");                        hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ApplicationRegion");                   hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("WZControlRegion");                     hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegionPresel");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ApplicationRegionPresel");             hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("WZControlRegionPresel");               hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawSignalRegion");                     hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawApplicationRegion");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawWZControlRegion");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawSignalRegionPresel");               hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawApplicationRegionPresel");          hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("RawWZControlRegionPresel");            hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  
  histonames.push_back("SignalRegion_JECup");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_JECdn");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_lepSFup");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_lepSFdn");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_bHFSFup");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_bHFSFdn");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_bLFSFup");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_bLFSFdn");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_PUup");                   hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SignalRegion_PUdn");                   hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);

  
  histonames.push_back("OldSignalRegion");                        hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldApplicationRegion");                   hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldWZControlRegion");                     hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegionPresel");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldApplicationRegionPresel");             hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldWZControlRegionPresel");               hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldRawSignalRegion");                     hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldRawApplicationRegion");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldRawWZControlRegion");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldRawSignalRegionPresel");               hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldRawApplicationRegionPresel");          hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldRawWZControlRegionPresel");            hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  
  histonames.push_back("OldSignalRegion_JECup");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_JECdn");                  hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_lepSFup");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_lepSFdn");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_bHFSFup");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_bHFSFdn");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_bLFSFup");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_bLFSFdn");                hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_PUup");                   hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("OldSignalRegion_PUdn");                   hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);


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

      bool continueold(false), continuenew(false);
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
      //if(string(currentFile->GetTitle()).find("www_2l_mia")!=string::npos)      weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")!=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(weight>100) cout << weight << " " << currentFile->GetTitle() << endl;
      if(isData()) weight = 1.;
      double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      float PUweight(1.), PUweightup(1.), PUweightdn(1.);
      if(applyPUrewgt&&!isData()){
	PUweight = getPUWeightAndError(PUweightdn,PUweightup);
	//weight *= PUweight;
      }
      
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
      LorentzVector MET_up; MET_up.SetPxPyPzE(met_up_pt()*TMath::Cos(met_up_phi()),met_up_pt()*TMath::Sin(met_up_phi()),0,met_up_pt());
      LorentzVector MET_dn; MET_dn.SetPxPyPzE(met_dn_pt()*TMath::Cos(met_dn_phi()),met_dn_pt()*TMath::Sin(met_dn_phi()),0,met_dn_pt());

      int njold(0), nbold(0), nj30old(0);
      getalljetnumbers(njold,nj30old,nbold);
      float Mjjold = -1;
      float MjjLold = -1; float Detajjold = -1;
      getMjjAndDeta(Mjjold,MjjLold,Detajjold);
      if(checkevent) cout << "nj30old " << nj30old << " njold " << njold << " nbold " << nbold << " Mjjold " << Mjjold << " MjjLold " << MjjLold << " Detajjold " << Detajjold << endl;
      if(checkevent) cout << "nj30 " << nj30() << " nj " << nj() << " nb " << nb() << " Mjj " << Mjj() << " MjjL " << MjjL() << " Detajj " << DetajjL() << endl;//new
      if(checkevent){
	for(unsigned int i = 0; i<jets_p4().size(); ++i){
	  cout << "jet pT " << jets_p4()[i].Pt() << " eta " << jets_p4()[i].Eta() << " CSV " << jets_csv()[i];// << endl;
	  for(unsigned int j = i+1; j<jets_p4().size(); ++j) cout << " M"<<i<<j<< " " << (jets_p4()[i]+jets_p4()[j]).M() << " (dR " << dR(jets_p4()[i],jets_p4()[j]) << ")";
	  cout << endl;
	}
      }

      int njold_up(0), nbold_up(0), nj30old_up(0);
      if(getJECunc) getalljetnumbers(njold_up,nj30old_up,nbold_up,1);
      float Mjjold_up = -1;
      float MjjLold_up = -1; float Detajjold_up = -1;
      if(getJECunc) getMjjAndDeta(Mjjold_up,MjjLold_up,Detajjold_up,1);

      int njold_dn(0), nbold_dn(0), nj30old_dn(0);
      if(getJECunc) getalljetnumbers(njold_dn,nj30old_dn,nbold_dn,-1);
      float Mjjold_dn = -1;
      float MjjLold_dn = -1; float Detajjold_dn = -1;
      if(getJECunc) getMjjAndDeta(Mjjold_dn,MjjLold_dn,Detajjold_dn,-1);
      
      if(nj30old!=nj30()&&nj30old>=0) cout << __LINE__ << " " << nj30old << " " << nj30() << endl;
      if(njold!=nj()&&njold>=0) cout << __LINE__ << " " << njold << " " << nj() << endl;
      if(nbold!=nb()&&nbold>=0) cout << __LINE__ << " " << nbold << " " << nb() << endl;
      if(Mjjold!=Mjj()&&Mjjold>0) cout << __LINE__ << " " << Mjjold << " " << Mjj() << endl;
      if(MjjLold!=MjjL()&&MjjLold>0) cout << __LINE__ << " " << MjjLold << " " << MjjL() << endl;
      if(Detajjold!=DetajjL()&&Detajjold>0) cout << __LINE__ << " " << Detajjold << " " << DetajjL() << endl;
      if(nj30old_up!=nj30_up()&&nj30old_up>=0) cout << __LINE__ << " " << nj30old_up << " " << nj30_up() << endl;
      if(njold_up!=nj_up()&&njold_up>=0) cout << __LINE__ << " " << njold_up << " " << nj_up() << endl;
      if(nbold_up!=nb_up()&&nbold_up>=0) cout << __LINE__ << " " << nbold_up << " " << nb_up() << endl;
      if(Mjjold_up!=Mjj_up()&&Mjjold_up>=0) cout << __LINE__ << " " << Mjjold_up << " " << Mjj_up() << endl;
      if(MjjLold_up!=MjjL_up()&&MjjLold_up>=0) cout << __LINE__ << " " << MjjLold_up << " " << MjjL_up() << endl;
      if(Detajjold_up!=DetajjL_up()&&Detajjold_up>=0) cout << __LINE__ << " " << Detajjold_up << " " << DetajjL_up() << endl;
      if(nj30old_dn!=nj30_dn()&&nj30old_dn>=0) cout << __LINE__ << " " << nj30old_dn << " " << nj30_dn() << endl;
      if(njold_dn!=nj_dn()&&njold_dn>=0) cout << __LINE__ << " " << njold_dn << " " << nj_dn() << endl;
      if(nbold_dn!=nb_dn()&&nbold_dn>=0) cout << __LINE__ << " " << nbold_dn << " " << nb_dn() << endl;
      if(Mjjold_dn!=Mjj_dn()&&Mjjold_dn>=0) cout << __LINE__ << " " << Mjjold_dn << " " << Mjj_dn() << endl;
      if(MjjLold_dn!=MjjL_dn()&&MjjLold_dn>=0) cout << __LINE__ << " " << MjjLold_dn << " " << MjjL_dn() << endl;
      if(Detajjold_dn!=DetajjL_dn()&&Detajjold_dn>=0) cout << __LINE__ << " " << Detajjold_dn << " " << DetajjL_dn() << endl;
      
      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      getleptonindices_v2(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l,25,20);
      float lepSFold(1.), lepSFolderr(0.);//i3l and iSS have same ID
      if(applylepSF&&!isData()){
	lepSFold = getlepSFWeightandError(lepSFolderr,i3l,ia3l);
	//weight *= lepSF;
      }
      float trigSFold(1.), trigSFolderr(1.);
      if(applytrigSF&&!isData()){
	trigSFold    = getTriggerWeightandError(trigSFolderr, i3l,ia3l);
	//weight *= trigSF;
      }
      float weight_oldlepSFup = weight;
      float weight_oldlepSFdn = weight;
      float weight_oldPUup    = weight;
      float weight_oldPUdn    = weight;
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
	weight_oldPUup *= PUweightup/PUweight;
	weight_oldPUdn *= PUweightdn/PUweight;
      }
      if(applylepSF&&!isData()&&lepSFold!=0){
	weight_oldlepSFup *= (lepSFold+lepSFolderr)/lepSFold;
	weight_oldlepSFdn *= (lepSFold-lepSFolderr)/lepSFold;
      }
      if(checkevent) cout << "old weight " << weight << " btag  " << weight_btagsf() << " PU " << PUweight << " trigold " << trigSFold << " lepold " << lepSFold << endl;
      if(checkevent) cout << "new weight " << weight << " btag  " << weight_btagsf() << " PU " << purewgt() << " trig " << trigeff() << " lep " << lepsf() << endl;

      if(checkevent){
	for(unsigned int i = 0; i<lep_pdgId().size();++i){
	  cout << "lep " << lep_pdgId()[i] << " Pt " << lep_p4()[i].Pt() << " eta " << lep_p4()[i].Eta() << " ID t/l/v/trig " << lep_pass_VVV_cutbased_tight_noiso()[i] << "/" << lep_pass_VVV_cutbased_fo_noiso()[i] << "/" << lep_pass_VVV_cutbased_veto_noiso()[i] << "/" << lep_isTriggerSafe_v1()[i] << " iso " << lep_relIso03EAv2()[i] << " ip3d " << lep_ip3d()[i] << " losthits " << lep_lostHits()[i] << " t.q " << lep_tightCharge()[i];
	  for(unsigned int j = i+1; j<lep_pdgId().size();++j) { cout << " M" << i << j << " " << (lep_p4()[i]+lep_p4()[j]).M();
	    for(unsigned int k = j+1; k<lep_pdgId().size();++k) cout << " M" << i << j << k << " " << (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).M() << " Pt " <<  (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).Pt() << " DPhiMET " << dPhi( (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]),MET); }
	  cout << endl;
	}
      }
      
      if((i3l.size()+ia3l.size())<2) continue;
      bool passofflineforTrigger = passofflineTriggers(i3l, ia3l);
      if(!passofflineforTrigger) continueold = true;
      if(checkevent&&!continueold) cout << "pass offline" << endl;
      
      if(isData()){
	if(!passFilters_v0()) continueold = true;
	if(!passFilters()) continuenew = true;
	if(checkevent) cout << "pass filter" << endl;
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) { continue; }
	if(checkevent) cout << "pass duplicate" << endl;
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	if(checkevent) cout << "pass goodrun" << endl;
	bool passonlineTrigger = passonlineTriggers(i3l, ia3l);//currently applied only to data
	if(!passonlineTrigger) continueold = true;
      } 
      if(!passTriggers()) continuenew = true;
      if(checkevent&&!continuenew) cout << "pass online" << endl;

      
      string sample   = skimFilePrefix;
      string oldsn    = ((iSS.size()+iaSS.size())>=2) ? process(fname,true ,iSS,iaSS) : string("not2l");
      string oldsn2   = ((i3l.size()+ia3l.size())>=3) ? process(fname,false,i3l,ia3l) : string("not3l");
      bool isphotonSS = (oldsn =="photonfakes");
      bool isphoton3l = (oldsn2=="photonfakes");
      if(splitVH(fname)){ sample = "WHtoWWW"; }
      string sn = string(bkgtype().Data());
      bool isphoton = (sn =="photonfakes");
      if((iSS.size()+iaSS.size())==2&&vaSS.size()==0&&oldsn !=sn) cout << __LINE__ << " - " << oldsn << " " << oldsn2 << " " << sn << " " << iSS.size() << " " << iaSS.size() << " " << vaSS.size() << " " << i3l.size() << " " << ia3l.size() << endl;
      if((i3l.size()+ia3l.size())>=3&&oldsn2!=sn) { cout << __LINE__ << ":  " << oldsn << " " << oldsn2 << " " << sn << " " << iSS.size() << " " << iaSS.size() << " " << vaSS.size() << " " << i3l.size() << " " << ia3l.size() << endl;
      }

      float MTmaxold = -1;
      if(iSS.size()==2) MTmaxold = calcMTmax(iSS,MET);
      else if(iSS.size()==1&&iaSS.size()>=1){
	vector<int> temp; temp.push_back(iSS[0]); temp.push_back(iaSS[0]);
	MTmaxold = calcMTmax(temp,MET);
      }
      float MTmax3lold = calcMTmax(i3l,MET,true);
      if(checkevent) cout << "MET " << MET.Pt() << " MTmaxold " << MTmaxold << " MTmax3lold " << MTmax3lold << " MTmax " << MTmax() << endl;
      float MTmaxold_up(-1), MTmaxold_dn(-1), MTmax3lold_up(-1), MTmax3lold_dn(-1);
      if(getJECunc) {
	MTmaxold_up   = calcMTmax(iSS,MET_up);
	//MTmax3lold_up = calcMTmax(i3l,MET_up,true);
	MTmaxold_dn   = calcMTmax(iSS,MET_dn);
	//MTmax3lold_dn = calcMTmax(i3l,MET_dn,true);
      }
      /*
      bool checkthis  = false;
      if(iSS.size()+iaSS.size()==2&&&&MTmaxold!=MTmax()) { cout << __LINE__ << " " << MTmaxold << " " << MTmax() << " " << iSS.size() << " " << iaSS.size() << endl; checkthis = true; }
      if(iSS.size()==2&&MTmaxold_up!=MTmax_up()) { cout << __LINE__ << " " << MTmaxold_up << " " << MTmax_up() << endl; checkthis = true; }
      if(iSS.size()==2&&MTmaxold_dn!=MTmax_dn()) { cout << __LINE__ << " " << MTmaxold_dn << " " << MTmax_dn() << endl; checkthis = true; }
      if(i3l.size()>=2&&MTmax3lold!=MTmax()) { cout << __LINE__ << " " << MTmax3lold << " " << MTmax() << endl;  checkthis = true; }
      if(checkthis){
	for(unsigned int i = 0; i<lep_pdgId().size();++i){
	  cout << "lep " << lep_pdgId()[i] << " Pt " << lep_p4()[i].Pt() << " eta " << lep_p4()[i].Eta() << " ID SS t/l " << lep_pass_VVV_cutbased_tight()[i] << "/" << lep_pass_VVV_cutbased_fo()[i] << " ID 3l t/l " << lep_pass_VVV_cutbased_3l_tight()[i] << "/" << lep_pass_VVV_cutbased_3l_fo()[i] << " met " << met_pt() << " phi " << met_phi() << " MT (up/down) " << mT(lep_p4()[i],MET) << " (" << mT(lep_p4()[i],MET_up) << "/" << mT(lep_p4()[i],MET_dn) << ") MTmax() (up/down) " << MTmax() << " (" << MTmax_up() << "/" << MTmax_dn() << ")" << endl;
	}
      }
      */
      
      int SRSS[20]; bool selects3l[20];
      int SR3l[20];
      int SRSStest[20]; bool selects3ltest[20];
      int SR3ltest[20];
      for(int i = 0; i<20; ++i) { SRSS[i] = -1; SR3l[i] = -1; selects3l[i] = false; }
      for(int i = 0; i<20; ++i) { SRSStest[i] = -1; SR3ltest[i] = -1; selects3ltest[i] = false; }

      //THIS is new version
      //SS
      //0: SR
      SRSS[0] = isSRSS();//enter variables for quicker calculation
      //1: SR preselect
      SRSS[1] = isSRSS(true);//enter variables for quicker calculation
      //2: AR
      SRSS[2] = isARSS();//enter variables for quicker calculation
      //3: AR preselect
      SRSS[3] = isARSS(true);//enter variables for quicker calculation
      //4: CR
      SRSS[4] = isCRSS();//enter variables for quicker calculation
      //5: CR preselect
      SRSS[5] = isCRSS(true);//enter variables for quicker calculation
      //6: SR JEC up
      if(getJECunc) SRSS[6] = isSRSS(false,1);
      //7: SR JEC dn
      if(getJECunc) SRSS[7] = isSRSS(false,-1);
      passAnySS(SRSStest[0],SRSStest[2],SRSStest[4]);
      passAnySS(SRSStest[1],SRSStest[3],SRSStest[5],true);
      selects3l[4] = true; selects3l[5] = true;
      selects3l[14] = true; selects3l[15] = true;
      bool inconsistent = false;
      if(SRSS[0]!=SRSStest[0]) { cout << __LINE__ << " SRf " << SRSS[0] << " " << SRSStest[0] << endl; inconsistent = true; }
      if(SRSS[1]!=SRSStest[1]) { cout << __LINE__ << " SRp " << SRSS[1] << " " << SRSStest[1] << endl; inconsistent = true; }
      if(SRSS[2]!=SRSStest[2]) { cout << __LINE__ << " ARf " << SRSS[2] << " " << SRSStest[2] << endl; inconsistent = true; }
      if(SRSS[3]!=SRSStest[3]) { cout << __LINE__ << " ARp " << SRSS[3] << " " << SRSStest[3] << endl; inconsistent = true; }
      if(SRSS[4]!=SRSStest[4]) { cout << __LINE__ << " CRf " << SRSS[4] << " " << SRSStest[4] << endl; inconsistent = true; }
      if(SRSS[5]!=SRSStest[5]) { cout << __LINE__ << " CRp " << SRSS[5] << " " << SRSStest[5] << endl; inconsistent = true; }
      //3l
      //0: SR
      SR3l[0] = isSR3l();
      //1: SR preselect
      SR3l[1] = isSR3l(true);
      //2: AR
      SR3l[2] = isAR3l();
      //3: AR preselect
      SR3l[3] = isAR3l(true);
      //4: CR
      SR3l[4] = isCR3l();
      //5: CR preselect
      SR3l[5] = isCR3l(true);
      //6: SR JEC up
      if(getJECunc) SR3l[6] = isSR3l(false,1);
      //7: SR JEC dn
      if(getJECunc) SR3l[7] = isSR3l(false,-1);
      passAny3l(SR3ltest[0],SR3ltest[2],SR3ltest[4]);
      passAny3l(SR3ltest[1],SR3ltest[3],SR3ltest[5],true);
      if(SR3l[0]!=SR3ltest[0]) { cout << __LINE__ << " " << SR3l[0] << " " << SR3ltest[0] << endl; inconsistent = true; }
      if(SR3l[1]!=SR3ltest[1]) { cout << __LINE__ << " " << SR3l[1] << " " << SR3ltest[1] << endl; inconsistent = true; }
      if(SR3l[2]!=SR3ltest[2]) { cout << __LINE__ << " " << SR3l[2] << " " << SR3ltest[2] << endl; inconsistent = true; }
      if(SR3l[3]!=SR3ltest[3]) { cout << __LINE__ << " " << SR3l[3] << " " << SR3ltest[3] << endl; inconsistent = true; }
      if(SR3l[4]!=SR3ltest[4]) { cout << __LINE__ << " " << SR3l[4] << " " << SR3ltest[4] << endl; inconsistent = true; }
      if(SR3l[5]!=SR3ltest[5]) { cout << __LINE__ << " " << SR3l[5] << " " << SR3ltest[5] << endl; inconsistent = true; }
      
      //THIS IS THE OLD VERSION
      //10: SR
      SRSS[10] = isSRSS_v0(iSS,      vSS,false,MTmaxold,  nj30old,nbold,Mjjold,MjjLold,Detajjold,MET,0,false,false,1);//enter variables for quicker calculation
      //11: SR preselect
      SRSS[11] = isSRSS_v0(iSS,      vSS,true ,MTmaxold,  nj30old,nbold,Mjjold,MjjLold,Detajjold,MET,0,false,false,1);//enter variables for quicker calculation
      //12: AR
      SRSS[12] = isARSS_v0(iSS,iaSS,vaSS,false,MTmaxold,  nj30old,nbold,Mjjold,MjjLold,Detajjold,MET,0,false,false,1);//enter variables for quicker calculation
      //13: AR preselect
      SRSS[13] = isARSS_v0(iSS,iaSS,vaSS,true ,MTmaxold,  nj30old,nbold,Mjjold,MjjLold,Detajjold,MET,0,false,false,1);//enter variables for quicker calculation
      //14: CR
      SRSS[14] = isCRSS_v0(iSS,i3l,  v3l,false,-1.*MTmax3lold,nj30old,nbold,Mjjold,MjjLold,Detajjold,MET,0,false,false,false,1);//enter variables for quicker calculation
      //15: CR preselect
      SRSS[15] = isCRSS_v0(iSS,i3l,  v3l,true ,-1.*MTmax3lold,nj30old,nbold,Mjjold,MjjLold,Detajjold,MET,0,false,false,false,1);//enter variables for quicker calculation
      selects3l[14] = true; selects3l[15] = true;
      //16: SR JEC up
      if(getJECunc) SRSS[16] = isSRSS_v0(iSS, vSS,false,MTmaxold_up, nj30old_up,nbold_up,Mjjold_up,MjjLold_up,Detajjold_up,MET_up,1,false,false,1);
      //17: SR JEC dn
      if(getJECunc) SRSS[17] = isSRSS_v0(iSS, vSS,false,MTmaxold_dn, nj30old_dn,nbold_dn,Mjjold_dn,MjjLold_dn,Detajjold_dn,MET_dn,-1,false,false,1);
      //10: SR, 14: CR
      checkbothSRCR3l_v0(SR3l[10],SR3l[14],i3l,false,njold,nbold,MET,0,false,1);
      //11: SR preselect, 15: CR preselect
      checkbothSRCR3l_v0(SR3l[11],SR3l[15],i3l,true ,njold,nbold,MET,0,false,1);
      //12: AR
      SR3l[12] = isAR3l_v0(i3l,ia3l,false,njold,nbold,MET,0,false,1);
      //13: AR preselect
      SR3l[13] = isAR3l_v0(i3l,ia3l,true ,njold,nbold,MET,0,false,1);
      //16: SR JEC up
      if(getJECunc) SR3l[16] = isSR3l_v0(i3l,false,njold_up,nbold_up,MET_up,1,false,1);
      //17: SR JEC dn
      if(getJECunc) SR3l[17] = isSR3l_v0(i3l,false,njold_dn,nbold_dn,MET_dn,-1,false,1);

      if(checkevent) cout << "passed          SRSS " << SRSS[ 0] << " SR3l " << SR3l[ 0] << " ARSS " << SRSS[ 2] << " AR3l " << SR3l[ 2] << " CRSS " << SRSS[ 4] << " CR3l " << SR3l[ 4] << endl;
      if(checkevent) cout << "passedold       SRSS " << SRSS[10] << " SR3l " << SR3l[10] << " ARSS " << SRSS[12] << " AR3l " << SR3l[12] << " CRSS " << SRSS[14] << " CR3l " << SR3l[14] << endl;
      for(int i = 0; i<10; ++i) {
	if(vetophoton())    { SRSS[i] = -1; SR3l[i] = -1; }
      }
      for(int i = 10; i<20; ++i) {
	if(!selects3l[i]){
	  if(vetophotonprocess(fname,isphotonSS))    { SRSS[i] = -1; }
	}
	else if(vetophotonprocess(fname,isphoton3l)){ SRSS[i] = -1; }
	if(vetophotonprocess(fname,isphoton3l))     { SR3l[i] = -1; }
      }
      //this is a very bad test - I guess I run into many mistakes
      bool fail = false;
      if(nVlep()<4){
	if(SRSS[0]!=SRSS[10]) { cout << "SR SS full   " << SRSS[0] << " " << SRSS[10] << endl; fail = true; }
	if(SRSS[1]!=SRSS[11]) { cout << "SR SS presel " << SRSS[1] << " " << SRSS[11] << endl; fail = true; }
	if(SRSS[2]!=SRSS[12]) { cout << "AR SS full   " << SRSS[2] << " " << SRSS[12] << endl; fail = true; }
	if(SRSS[3]!=SRSS[13]) { cout << "AR SS presel " << SRSS[3] << " " << SRSS[13] << endl; fail = true; }
	if(SRSS[4]!=SRSS[14]) { cout << "CR SS full   " << SRSS[4] << " " << SRSS[14] << endl; fail = true; }
	if(SRSS[5]!=SRSS[15]) { cout << "CR SS presel " << SRSS[5] << " " << SRSS[15] << endl; fail = true; }
	if(SRSS[6]!=SRSS[16]) { cout << "JU SS full   " << SRSS[6] << " " << SRSS[16] << endl; fail = true; }
	if(SRSS[7]!=SRSS[17]) { cout << "JD SS full   " << SRSS[7] << " " << SRSS[17] << endl; fail = true; }
	if(SR3l[0]!=SR3l[10]) { cout << "SR 3l full   " << SR3l[0] << " " << SR3l[10] << endl; fail = true; }
	if(SR3l[1]!=SR3l[11]) { cout << "SR 3l presel " << SR3l[1] << " " << SR3l[11] << endl; fail = true; }
	if(SR3l[2]!=SR3l[12]) { cout << "AR 3l full   " << SR3l[2] << " " << SR3l[12] << endl; fail = true; }
	if(SR3l[3]!=SR3l[13]) { cout << "AR 3l presel " << SR3l[3] << " " << SR3l[13] << endl; fail = true; }
	if(SR3l[4]!=SR3l[14]) { cout << "CR 3l full   " << SR3l[4] << " " << SR3l[14] << endl; fail = true; }
	if(SR3l[5]!=SR3l[15]) { cout << "CR 3l presel " << SR3l[5] << " " << SR3l[15] << endl; fail = true; }
	if(SR3l[6]!=SR3l[16]) { cout << "JU 3l full   " << SR3l[6] << " " << SR3l[16] << endl; fail = true; }
	if(SR3l[7]!=SR3l[17]) { cout << "JD 3l full   " << SR3l[7] << " " << SR3l[17] << endl; fail = true; }
      }

      if(checkevent) cout << "photonpassed    SRSS " << SRSS[ 0] << " SR3l " << SR3l[ 0] << " ARSS " << SRSS[ 2] << " AR3l " << SR3l[ 2] << " CRSS " << SRSS[ 4] << " CR3l " << SR3l[ 4] << endl;
      if(checkevent) cout << "photonpassedold SRSS " << SRSS[10] << " SR3l " << SR3l[10] << " ARSS " << SRSS[12] << " AR3l " << SR3l[12] << " CRSS " << SRSS[14] << " CR3l " << SR3l[14] << endl;
      if(!continueold){
	if(!isData()||!blindSR){//SR is blinded
	  fillSRhisto(histos, "OldSignalRegion",               sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight*lepSFold*trigSFold*PUweight,            weight*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegionPresel",         sample, oldsn, oldsn2, SRSS[11], SR3l[11], weight*lepSFold*trigSFold*PUweight,            weight*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_JECup",         sample, oldsn, oldsn2, SRSS[16], SR3l[16], weight*lepSFold*trigSFold*PUweight,            weight*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_JECdn",         sample, oldsn, oldsn2, SRSS[17], SR3l[17], weight*lepSFold*trigSFold*PUweight,            weight*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_lepSFup",       sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_oldlepSFup*lepSFold*trigSFold*PUweight, weight_oldlepSFup*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_lepSFdn",       sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_oldlepSFdn*lepSFold*trigSFold*PUweight, weight_oldlepSFdn*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_bHFSFup",       sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_bHFSFup*lepSFold*trigSFold*PUweight,    weight_bHFSFup*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_bHFSFdn",       sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_bHFSFdn*lepSFold*trigSFold*PUweight,    weight_bHFSFdn*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_bLFSFup",       sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_bLFSFup*lepSFold*trigSFold*PUweight,    weight_bLFSFup*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_bLFSFdn",       sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_bLFSFdn*lepSFold*trigSFold*PUweight,    weight_bLFSFdn*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_PUup",          sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_oldPUup*lepSFold*trigSFold*PUweight,    weight_oldPUup*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldSignalRegion_PUdn",          sample, oldsn, oldsn2, SRSS[10], SR3l[10], weight_oldPUdn*lepSFold*trigSFold*PUweight,    weight_oldPUdn*lepSFold*trigSFold*PUweight);
	  fillSRhisto(histos, "OldRawSignalRegion",            sample, oldsn, oldsn2, SRSS[10], SR3l[10], 1., 1.);
	  fillSRhisto(histos, "OldRawSignalRegionPresel",      sample, oldsn, oldsn2, SRSS[11], SR3l[11], 1., 1.);	
	}
	fillSRhisto(  histos, "OldApplicationRegion",          sample, oldsn, oldsn2, SRSS[12], SR3l[12], weight*lepSFold*trigSFold*PUweight, weight*lepSFold*trigSFold*PUweight);
	fillSRhisto(  histos, "OldApplicationRegionPresel",    sample, oldsn, oldsn2, SRSS[13], SR3l[13], weight*lepSFold*trigSFold*PUweight, weight*lepSFold*trigSFold*PUweight);
	fillSRhisto(  histos, "OldWZControlRegion",            sample, oldsn2,oldsn2, SRSS[14], SR3l[14], weight*lepSFold*trigSFold*PUweight, weight*lepSFold*trigSFold*PUweight);
	fillSRhisto(  histos, "OldWZControlRegionPresel",      sample, oldsn2,oldsn2, SRSS[15], SR3l[15], weight*lepSFold*trigSFold*PUweight, weight*lepSFold*trigSFold*PUweight);
	fillSRhisto(  histos, "OldRawApplicationRegion",       sample, oldsn, oldsn2, SRSS[12], SR3l[12], 1., 1.);
	fillSRhisto(  histos, "OldRawApplicationRegionPresel", sample, oldsn, oldsn2, SRSS[13], SR3l[13], 1., 1.);
	fillSRhisto(  histos, "OldRawWZControlRegion",         sample, oldsn2,oldsn2, SRSS[14], SR3l[14], 1., 1.);
	fillSRhisto(  histos, "OldRawWZControlRegionPresel",   sample, oldsn2,oldsn2, SRSS[15], SR3l[15], 1., 1.);

      }
      if(!continuenew){
	float addfac = 1.;
	if(applylepSF)   addfac = lepsf();
	if(applytrigSF)  addfac = trigeff();
	if(applyPUrewgt) addfac = purewgt();
	if(!isData()||!blindSR){//SR is blinded
	  fillSRhisto(histos, "SignalRegion",               sample, sn, sn, SRSS[0], SR3l[0], weight*addfac,            weight*addfac);
	  fillSRhisto(histos, "SignalRegionPresel",         sample, sn, sn, SRSS[1], SR3l[1], weight*addfac,            weight*addfac);
	  fillSRhisto(histos, "SignalRegion_JECup",         sample, sn, sn, SRSS[6], SR3l[6], weight*addfac,            weight*addfac);
	  fillSRhisto(histos, "SignalRegion_JECdn",         sample, sn, sn, SRSS[7], SR3l[7], weight*addfac,            weight*addfac);
	  fillSRhisto(histos, "SignalRegion_lepSFup",       sample, sn, sn, SRSS[0], SR3l[0], weight*lepsf_up()*addfac/lepsf(), weight*lepsf_up()*addfac/lepsf());
	  fillSRhisto(histos, "SignalRegion_lepSFdn",       sample, sn, sn, SRSS[0], SR3l[0], weight*lepsf_dn()*addfac/lepsf(), weight*lepsf_dn()*addfac/lepsf());
	  fillSRhisto(histos, "SignalRegion_bHFSFup",       sample, sn, sn, SRSS[0], SR3l[0], weight_bHFSFup*addfac,    weight_bHFSFup*addfac);
	  fillSRhisto(histos, "SignalRegion_bHFSFdn",       sample, sn, sn, SRSS[0], SR3l[0], weight_bHFSFdn*addfac,    weight_bHFSFdn*addfac);
	  fillSRhisto(histos, "SignalRegion_bLFSFup",       sample, sn, sn, SRSS[0], SR3l[0], weight_bLFSFup*addfac,    weight_bLFSFup*addfac);
	  fillSRhisto(histos, "SignalRegion_bLFSFdn",       sample, sn, sn, SRSS[0], SR3l[0], weight_bLFSFdn*addfac,    weight_bLFSFdn*addfac);
	  fillSRhisto(histos, "SignalRegion_PUup",          sample, sn, sn, SRSS[0], SR3l[0], weight*purewgt_up()*addfac/purewgt(),    weight*purewgt_up()*addfac/purewgt());
	  fillSRhisto(histos, "SignalRegion_PUdn",          sample, sn, sn, SRSS[0], SR3l[0], weight*purewgt_dn()*addfac/purewgt(),    weight*purewgt_dn()*addfac/purewgt());
	  fillSRhisto(histos, "RawSignalRegion",            sample, sn, sn, SRSS[0], SR3l[0], 1., 1.);
	  fillSRhisto(histos, "RawSignalRegionPresel",      sample, sn, sn, SRSS[1], SR3l[1], 1., 1.);	
	}
	fillSRhisto(  histos, "ApplicationRegion",          sample, sn, sn, SRSS[2], SR3l[2], weight*addfac, weight*addfac);
	fillSRhisto(  histos, "ApplicationRegionPresel",    sample, sn, sn, SRSS[3], SR3l[3], weight*addfac, weight*addfac);
	fillSRhisto(  histos, "WZControlRegion",            sample, sn, sn, SRSS[4], SR3l[4], weight*addfac, weight*addfac);
	fillSRhisto(  histos, "WZControlRegionPresel",      sample, sn, sn, SRSS[5], SR3l[5], weight*addfac, weight*addfac);
	fillSRhisto(  histos, "RawApplicationRegion",       sample, sn, sn, SRSS[2], SR3l[2], 1., 1.);
	fillSRhisto(  histos, "RawApplicationRegionPresel", sample, sn, sn, SRSS[3], SR3l[3], 1., 1.);
	fillSRhisto(  histos, "RawWZControlRegion",         sample, sn, sn, SRSS[4], SR3l[4], 1., 1.);
	fillSRhisto(  histos, "RawWZControlRegionPresel",   sample, sn, sn, SRSS[5], SR3l[5], 1., 1.);

      }
      if(storeeventnumbers){
	addeventtolist(SRSS[0], SR3l[0], SREE, SREM, SRMM, SR0SFOS, SR1SFOS, SR2SFOS);
	addeventtolist(SRSS[2], SR3l[2], AREE, AREM, ARMM, AR0SFOS, AR1SFOS, AR2SFOS);
	addeventtolist(SRSS[4], SR3l[4], CREE, CREM, CRMM, CR0SFOS, CR1SFOS, CR2SFOS);
	//addeventtolist(SRSS[], SR3l[], SREE, SREM, SRMM, SR0SFOS, SR1SFOS, SR2SFOS);
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
  
  SaveHistosToFile("rootfiles/SRLooper_Test.root",histos,true,true,(chainnumber==0));
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
