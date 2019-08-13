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

  bool blindSR         = true;
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
  histonames.push_back("MTmax_emu_presel");            hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);//new
  histonames.push_back("MTmax_emu_NmO");               hbins.push_back(15);  hlow.push_back( 0); hup.push_back(225);//new
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
  histonames.push_back("MT3rd_1SFOS_presel");          hbins.push_back(12);  hlow.push_back( 0); hup.push_back(180);//new
  histonames.push_back("MT3rd_1SFOS_NmO");             hbins.push_back(12);  hlow.push_back( 0); hup.push_back(180);//new
  
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

      int SRSSpresel = isSRSS(true); 
      int SR3lpresel = isSR3l(true);
      int CR3lpresel = isCR3l(true);
      
      if(vetophoton()) continue;

      if(SRSSpresel>=0){
	fillhisto(histos, "Mjj_SRSS_presel",     sample, sn, Mjj(),                          weight);
	fillhisto(histos, "Mll_SRSS_presel",     sample, sn, MllSS(),                        weight);
	fillhisto(histos, "pTll_SRSS_presel",    sample, sn, (lep_p4()[0]+lep_p4()[1]).Pt(), weight);
	fillhisto(histos, "MET_SRSS_presel",     sample, sn, met_pt(),                       weight);
	fillhisto(histos, "MjjL_SRSS_presel",    sample, sn, MjjL(),                         weight);
	fillhisto(histos, "DetajjL_SRSS_presel", sample, sn, DetajjL(),                      weight);
	fillhisto(histos, "MTmax_SRSS_presel",   sample, sn, MTmax(),                        weight);
	if(SRSSpresel==1){
	  fillhisto(histos, "MTmax_emu_presel",  sample, sn, MTmax(),                        weight);
	}
	if(!isData()||!blindSR){
	  bool passMjj  = fabs(Mjj()-80.)<15.;
	  bool passMjjL = MjjL()<400.;
	  bool passDeta = DetajjL()<1.5;
	  bool passMll  = ((SRSSpresel==1) ? (MllSS()>30.) : (MllSS()>40.));
	  bool passMET  = ((SRSSpresel==2) ? (true)        : (met_pt()>60.));
	  bool passMT   = ((SRSSpresel==1) ? (MTmax()>30.) : (true));
	  if(         passMjjL&&passDeta&&passMll&&passMET&&passMT) fillhisto(histos, "Mjj_SRSS_NmO",     sample, sn, Mjj(),                          weight);
	  if(passMjj&&passMjjL&&passDeta         &&passMET&&passMT) fillhisto(histos, "Mll_SRSS_NmO",     sample, sn, MllSS(),                        weight);
	  if(passMjj&&passMjjL&&passDeta&&passMll&&passMET&&passMT) fillhisto(histos, "pTll_SRSS_NmO",    sample, sn, (lep_p4()[0]+lep_p4()[1]).Pt(), weight);
	  if(passMjj&&passMjjL&&passDeta&&passMll         &&passMT) fillhisto(histos, "MET_SRSS_NmO",     sample, sn, met_pt(),                       weight);
	  if(passMjj          &&passDeta&&passMll&&passMET&&passMT) fillhisto(histos, "MjjL_SRSS_NmO",    sample, sn, MjjL(),                         weight);
	  if(passMjj&&passMjjL          &&passMll&&passMET&&passMT) fillhisto(histos, "DetajjL_SRSS_NmO", sample, sn, DetajjL(),                      weight);
	  if(passMjj&&passMjjL&&passDeta&&passMll&&passMET        ) fillhisto(histos, "MTmax_SRSS_NmO",   sample, sn, MTmax(),                        weight);
	  if(SRSSpresel==1){
	    if(passMjj&&passMjjL&&passDeta&&passMll&&passMET      ) fillhisto(histos, "MTmax_emu_NmO",    sample, sn, MTmax(),                        weight);
	  }
	}
      }
      if(SR3lpresel>=0){
	bool passMlll       = (fabs(M3l()-MZ)>=10.);
	bool passMET        = ((SR3lpresel==0) ? (met_pt()>30.) : ((SRSSpresel==1) ? (met_pt()>40.) : (met_pt()>55.)));
	bool passMT3rd      = ((SR3lpresel==1) ? (MT3rd()>90.) : (true));
	bool passDPhilllMET = (DPhi3lMET()>2.5);
	bool passPtlll      = ((SR3lpresel==0) ? (true) : (Pt3l()>60.));
	fillhisto(  histos, "DPhilllMET_SR3l_presel", sample, sn, DPhi3lMET(), weight);
	fillhisto(  histos, "PTlll_SR3l_presel",      sample, sn, Pt3l(),      weight);
	fillhisto(  histos, "Mlll_SR3l_presel",       sample, sn, M3l(),       weight);
	fillhisto(  histos, "MET_SR3l_presel",        sample, sn, met_pt(),    weight);
	if(SR3lpresel==1)
	  fillhisto(histos, "MT3rd_1SFOS_presel",     sample, sn, MT3rd(),     weight);
	if(SR3lpresel>=1)
	  fillhisto(histos, "MSFOSZlike_SR3l_presel", sample, sn, Mll3L(),     weight);
	if(               passMlll&&passMET&&passMT3rd&&                passPtlll) fillhisto(histos, "DPhilllMET_SR3l_NmO", sample, sn, DPhi3lMET(), weight);
	if(               passMlll&&passMET&&passMT3rd&&passDPhilllMET           ) fillhisto(histos, "PTlll_SR3l_NmO",      sample, sn, Pt3l(),      weight);
	if(                         passMET&&passMT3rd&&passDPhilllMET&&passPtlll) fillhisto(histos, "Mlll_SR3l_NmO",       sample, sn, M3l(),       weight);
	if(               passMlll&&         passMT3rd&&passDPhilllMET&&passPtlll) fillhisto(histos, "MET_SR3l_NmO",        sample, sn, met_pt(),    weight);
	if(SR3lpresel>=1&&passMlll&&passMET&&           passDPhilllMET&&passPtlll) fillhisto(histos, "MT3rd_1SFOS_NmO",     sample, sn, MT3rd(),     weight);
	if(SR3lpresel>=1&&passMlll&&passMET&&passMT3rd&&passDPhilllMET&&passPtlll) fillhisto(histos, "MSFOSZlike_SR3l_NmO", sample, sn, Mll3L(),     weight);
      }
      if(CR3lpresel>=1){
	bool passMlll       = (fabs(M3l()-MZ)>=10.);
	bool passMET        = ((CR3lpresel==0) ? (met_pt()>30.) : ((CR3lpresel==1) ? (met_pt()>40.) : (met_pt()>55.)));
	bool passMT3rd      = ((CR3lpresel==1) ? (MT3rd()>90.) : (true));
	bool passDPhilllMET = (DPhi3lMET()>2.5);
	bool passPtlll      = ((CR3lpresel==0) ? (true) : (Pt3l()>60.));
	if(CR3lpresel==1) fillhisto(histos, "MSFOSZlike_SR3l_presel", sample, sn, Mll3L(),     weight);
	if(passMlll&&passMET&&passMT3rd&&passDPhilllMET&&passPtlll) fillhisto(histos, "MSFOSZlike_SR3l_NmO", sample, sn, Mll3L(),     weight);
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
  
  SaveHistosToFile("rootfiles/NminusOnePlots.root",histos,true,true,(chainnumber==0));
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
