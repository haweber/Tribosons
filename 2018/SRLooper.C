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

  bool blindSR         = false;
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
  vector<int>    hbins;      hbins.clear();
  vector<float>  hlow;       hlow.clear();
  vector<float>  hup;        hup.clear();

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
      bool checkevent = false;
      for(unsigned int i = 0; i<e.size();++i){
	if(e[i].run!=tas::run() ) continue;
	if(e[i].ls !=tas::lumi()) continue;
	if(e[i].evt!=tas::evt() ) continue;
	checkevent = true;
	cout << "Check event " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	break;
      }
      if(checkevent) {
	cout << "weight " << weight << " btag  " << weight_btagsf() << " PU " << purewgt() << " trig " << trigeff() << " lep " << lepsf() << endl;
	cout << "nj30 "   << nj30() << " nj " << nj() << " nb " << nb() << " Mjj " << Mjj() << " MjjL " << MjjL() << " Detajj " << DetajjL() << endl;
	for(unsigned int i = 0; i<jets_p4().size(); ++i){
	  cout << "jet pT " << jets_p4()[i].Pt() << " eta " << jets_p4()[i].Eta() << " CSV " << jets_csv()[i];
	  for(unsigned int j = i+1; j<jets_p4().size(); ++j) cout << " M"<<i<<j<< " " << (jets_p4()[i]+jets_p4()[j]).M() << " (dR " << dR(jets_p4()[i],jets_p4()[j]) << ")";
	  cout << endl;
	}
	for(unsigned int i = 0; i<lep_pdgId().size();++i){
	  cout << "lep " << lep_pdgId()[i] << " Pt " << lep_p4()[i].Pt() << " eta " << lep_p4()[i].Eta() << " IDnoiso t/l/v/trig " << lep_pass_VVV_cutbased_tight_noiso()[i] << "/" << lep_pass_VVV_cutbased_fo_noiso()[i] << "/" << lep_pass_VVV_cutbased_veto_noiso()[i] << "/" << lep_isTriggerSafe_v1()[i];
	  cout << " ID SSt/SSl/3lt/3ll " << lep_pass_VVV_cutbased_tight()[i] << "/" << lep_pass_VVV_cutbased_fo()[i] << "/" << lep_pass_VVV_cutbased_3l_tight()[i] << "/" << lep_pass_VVV_cutbased_3l_fo()[i];
	  cout << " iso " << lep_relIso03EAv2()[i] << " ip3d " << lep_ip3d()[i] << " losthits " << lep_lostHits()[i] << " t.q " << lep_tightCharge()[i];
	  for(unsigned int j = i+1; j<lep_pdgId().size();++j) { cout << " M" << i << j << " " << (lep_p4()[i]+lep_p4()[j]).M();
	    for(unsigned int k = j+1; k<lep_pdgId().size();++k) cout << " M" << i << j << k << " " << (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).M() << " Pt " <<  (lep_p4()[i]+lep_p4()[j]+lep_p4()[k]).Pt(); }
	  cout << endl;
	}
      }
            
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
      
      int SRSS[20]; 
      int SR3l[20];
      for(int i = 0; i<20; ++i) { SRSS[i] = -1; SR3l[i] = -1;  }

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
      if(vetophoton()) continue;

      if(!isData()||!blindSR){//SR is blinded
	fillSRhisto(histos, "SignalRegion",               sample, sn, sn, SRSS[0], SR3l[0], weight, weight);
	fillSRhisto(histos, "SignalRegionPresel",         sample, sn, sn, SRSS[1], SR3l[1], weight, weight);
	fillSRhisto(histos, "RawSignalRegion",            sample, sn, sn, SRSS[0], SR3l[0], 1.,     1.);
	fillSRhisto(histos, "RawSignalRegionPresel",      sample, sn, sn, SRSS[1], SR3l[1], 1.,     1.);
	if(!isData()&&getJECunc) {
	  fillSRhisto(histos, "SignalRegion_JECup",       sample, sn, sn, SRSS[6], SR3l[6], weight, weight);
	  fillSRhisto(histos, "SignalRegion_JECdn",       sample, sn, sn, SRSS[7], SR3l[7], weight, weight);
	}
	if(!isData()&&applylepSF&&lepsf()!=0){
	  fillSRhisto(histos, "SignalRegion_lepSFup",     sample, sn, sn, SRSS[0], SR3l[0], weight*lepsf_up()/lepsf(),                       weight*lepsf_up()/lepsf());
	  fillSRhisto(histos, "SignalRegion_lepSFdn",     sample, sn, sn, SRSS[0], SR3l[0], weight*lepsf_dn()/lepsf(),                       weight*lepsf_dn()/lepsf());
	}
	if(!isData()&&btagreweighting&&weight_btagsf()!=0){
	  fillSRhisto(histos, "SignalRegion_bHFSFup",     sample, sn, sn, SRSS[0], SR3l[0], weight*weight_btagsf_heavy_UP()/weight_btagsf(), weight*weight_btagsf_heavy_UP()/weight_btagsf());
	  fillSRhisto(histos, "SignalRegion_bHFSFdn",     sample, sn, sn, SRSS[0], SR3l[0], weight*weight_btagsf_heavy_DN()/weight_btagsf(), weight*weight_btagsf_heavy_DN()/weight_btagsf());
	  fillSRhisto(histos, "SignalRegion_bLFSFup",     sample, sn, sn, SRSS[0], SR3l[0], weight*weight_btagsf_light_UP()/weight_btagsf(), weight*weight_btagsf_light_UP()/weight_btagsf());
	  fillSRhisto(histos, "SignalRegion_bLFSFdn",     sample, sn, sn, SRSS[0], SR3l[0], weight*weight_btagsf_light_DN()/weight_btagsf(), weight*weight_btagsf_light_DN()/weight_btagsf());
	}
	if(!isData()&&applyPUrewgt&&!isData()&&purewgt()!=0){
	  fillSRhisto(histos, "SignalRegion_PUup",        sample, sn, sn, SRSS[0], SR3l[0], weight*purewgt_up()/purewgt(),                   weight*purewgt_up()/purewgt());
	  fillSRhisto(histos, "SignalRegion_PUdn",        sample, sn, sn, SRSS[0], SR3l[0], weight*purewgt_dn()/purewgt(),                   weight*purewgt_dn()/purewgt());
	}	
      }
      fillSRhisto(  histos, "ApplicationRegion",          sample, sn, sn, SRSS[2], SR3l[2], weight, weight);
      fillSRhisto(  histos, "ApplicationRegionPresel",    sample, sn, sn, SRSS[3], SR3l[3], weight, weight);
      fillSRhisto(  histos, "WZControlRegion",            sample, sn, sn, SRSS[4], SR3l[4], weight, weight);
      fillSRhisto(  histos, "WZControlRegionPresel",      sample, sn, sn, SRSS[5], SR3l[5], weight, weight);
      fillSRhisto(  histos, "RawApplicationRegion",       sample, sn, sn, SRSS[2], SR3l[2], 1.,     1.);
      fillSRhisto(  histos, "RawApplicationRegionPresel", sample, sn, sn, SRSS[3], SR3l[3], 1.,     1.);
      fillSRhisto(  histos, "RawWZControlRegion",         sample, sn, sn, SRSS[4], SR3l[4], 1.,     1.);
      fillSRhisto(  histos, "RawWZControlRegionPresel",   sample, sn, sn, SRSS[5], SR3l[5], 1.,     1.);

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
  
  SaveHistosToFile("rootfiles/SRLooper.root",histos,true,true,(chainnumber==0));
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
