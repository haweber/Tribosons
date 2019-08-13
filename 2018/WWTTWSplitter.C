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

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int>    hbins;      hbins.clear();
  vector<float>  hlow;       hlow.clear();
  vector<float>  hup;        hup.clear();

  histonames.push_back("SRyield");                        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir,true);
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
      if(sample!="ttV"&&sample!="WW") continue;
      
      int SRSS = isSRSS(); 
      int SR3l = isSR3l();

      if(vetophoton()) continue;
      
      bool trueSS = (sn =="trueSS");
      bool true3l = (sn=="true3l"||sn=="trueWWW");      
      bool isWWnotVBSDPS = (fname.find("ww_2l2nu_powheg")  !=string::npos)||(fname.find("ww_2l2nu_powheg")  !=string::npos);
      bool isTTX = (fname.find("ttg_")  !=string::npos)||(fname.find("tth_")  !=string::npos)||(fname.find("ttw_")  !=string::npos)||(fname.find("ttz_")  !=string::npos);
      string sn2 = sn;
      
      if(     NtightSS()==2&&NVlep()==2&&(fname.find("wpwpjj_ewk-")      !=string::npos)&&trueSS){ sn = "WWVBS";     sn2 = "WWRest"; }
      else if(NtightSS()==2&&NVlep()==2&&(fname.find("ww_2l2nu_dbl_scat")!=string::npos)&&trueSS){ sn = "WWDPS";     sn2 = "WWRest"; }
      else if(NtightSS()==2&&NVlep()==2&&isWWnotVBSDPS&&                                  trueSS){ sn = "otherWWSS"; sn2 = "WWRest"; }
      else if(                           isWWnotVBSDPS                                          ){ sn = "WWRest";    sn2 = "WWRest"; }
      if(     NtightSS()==2&&NVlep()==2&&isTTX&&trueSS){ sn  = "ttVSS";                 sn2 = "OtherttV"; }
      else if(Ntight3l()==3&&NVlep()==3&&isTTX&&true3l){ sn2 = "ttV3l"; if(sn!="ttVSS") sn  = "OtherttV"; }
      else if(                           isTTX        ){ sn  = "OtherttV";              sn2 = "OtherttV"; }

      if(!(blindSR&&isData())) fillSRhisto(histos, "SRyield",       sample, sn, sn2, SRSS, SR3l, weight, weight);      
      
    }//event loop
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  SaveHistosToFile("rootfiles/WWTTWSplitted.root",histos,true,true,(chainnumber==0));
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
