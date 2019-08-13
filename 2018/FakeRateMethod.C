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
//using namespace mt2_bisect;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1) {

  bool blindSR         = true;
  bool btagreweighting = true;
  bool applylepSF      = true;
  bool applytrigSF     = true;
  bool applyPUrewgt    = true;
  bool closuretest     = false;//if true use QCD fake rate
  
  bool isBtagSB        = false;
  
  bool checkMjjSB      = false;
  bool checkMjjlowMETSB= false;
  
  bool skippreselection= false;
  bool fillshapes      = false;
  
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
  histonames.push_back("FakeEstimationEvsMu");                             hbins.push_back( 2); hlow.push_back(0); hup.push_back(2);
  histonames.push_back("SRyieldEvsMu");                                    hbins.push_back( 2); hlow.push_back(0); hup.push_back(2);
  histonames.push_back("ARyieldEvsMu");                                    hbins.push_back( 2); hlow.push_back(0); hup.push_back(2);
  histonames.push_back("PreselFakeEstimationEvsMu");                       hbins.push_back( 2); hlow.push_back(0); hup.push_back(2);
  histonames.push_back("PreselSRyieldEvsMu");                              hbins.push_back( 2); hlow.push_back(0); hup.push_back(2);
  histonames.push_back("PreselARyieldEvsMu");                              hbins.push_back( 2); hlow.push_back(0); hup.push_back(2);

  histonames.push_back("SRyield");                                         hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselSRyield");                                   hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  //use cone correction only to extract fake rate
  histonames.push_back("ARyield");                                         hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimation");                                  hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRup");                              hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRdn");                              hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationClosureup");                         hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationClosuredn");                         hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  //use cone correction only to extract fake rate - preselection
  histonames.push_back("PreselARyield");                                   hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimation");                            hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRup");                        hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRdn");                        hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationClosureup");                   hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationClosuredn");                   hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);

  histonames.push_back("njets_SRSS");                                      hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("nbjets_SRSS");                                     hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("pTl1_SRSS");                                       hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("pTl2_SRSS");                                       hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("Mjj_SRSS");                                        hbins.push_back(15); hlow.push_back(0); hup.push_back(300);
  histonames.push_back("Mll_SRSS");                                        hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("MET_SRSS");                                        hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("MTmax_SRSS");                                      hbins.push_back(15); hlow.push_back(0); hup.push_back(225);

  histonames.push_back("SRyield_Mjjsideband");                             hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselSRyield_Mjjsideband");                       hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("ARyield_Mjjsideband");                             hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimation_Mjjsideband");                      hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationUp_Mjjsideband");                    hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationDn_Mjjsideband");                    hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselARyield_Mjjsideband");                       hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimation_Mjjsideband");                hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationUp_Mjjsideband");              hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationDn_Mjjsideband");              hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  
  histonames.push_back("njets_SRSS_Mjjsideband");                          hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("nbjets_SRSS_Mjjsideband");                         hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("pTl1_SRSS_Mjjsideband");                           hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("pTl2_SRSS_Mjjsideband");                           hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("Mjj_SRSS_Mjjsideband");                            hbins.push_back(15); hlow.push_back(0); hup.push_back(300);
  histonames.push_back("Mll_SRSS_Mjjsideband");                            hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("MET_SRSS_Mjjsideband");                            hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("MTmax_SRSS_Mjjsideband");                          hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  
  histonames.push_back("SRyield_MjjsidebandlowMET");                       hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselSRyield_MjjsidebandlowMET");                 hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("ARyield_MjjsidebandlowMET");                       hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimation_MjjsidebandlowMET");                hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationUp_MjjsidebandlowMET");              hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("FakeEstimationDn_MjjsidebandlowMET");              hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselARyield_MjjsidebandlowMET");                 hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimation_MjjsidebandlowMET");          hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationUp_MjjsidebandlowMET");        hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationDn_MjjsidebandlowMET");        hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);

  histonames.push_back("njets_SRSS_MjjsidebandlowMET");                    hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("nbjets_SRSS_MjjsidebandlowMET");                   hbins.push_back( 6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("pTl1_SRSS_MjjsidebandlowMET");                     hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("pTl2_SRSS_MjjsidebandlowMET");                     hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("Mjj_SRSS_MjjsidebandlowMET");                      hbins.push_back(15); hlow.push_back(0); hup.push_back(300);
  histonames.push_back("Mll_SRSS_MjjsidebandlowMET");                      hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  histonames.push_back("MET_SRSS_MjjsidebandlowMET");                      hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("MTmax_SRSS_MjjsidebandlowMET");                    hbins.push_back(15); hlow.push_back(0); hup.push_back(225);
  
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

      if(firstgoodvertex()!=0)                 continue;
      if(nVert()<0)                            continue;
      if(nTlepSS()<1&&nTlep()<2)               continue;
      if(nLlep()<2)                            continue;
      if(!isBtagSB&&nb() >0)                   continue;
      if( isBtagSB&&nb()==0)                   continue;
      
      //weight = 1;
      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      //if(string(currentFile->GetTitle()).find("www_2l_mia")!=string::npos) weight *= 0.066805*91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")!=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(weight>100) cout << weight << " " << currentFile->GetTitle() << endl;
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
      bool isfake     = (sn =="fakes");
      float subtract = 1;
      if( !isData()&&(!isfake)) subtract = -1;
      if(vetophoton()) continue;
      
      int SRSS[12];
      int SR3l[12];
      int dummy = -1;
      for(int i = 0; i<4; ++i) { SRSS[i] = -1; SR3l[i] = -1; }

      passAnySS(                      SRSS[ 0],SRSS[ 2],dummy,false,0,isBtagSB,false,1); //full:         0: SR,  2: AR
      if(!skippreselection) passAnySS(SRSS[ 1],SRSS[ 3],dummy,true, 0,isBtagSB,false,1); //preselection: 1: SR,  3: AR
      //3l
      passAny3l(                      SR3l[ 0],SR3l[ 2],dummy,false,0,isBtagSB,1);       //full:         0: SR,  2: AR
      if(!skippreselection) passAny3l(SR3l[ 1],SR3l[ 3],dummy,true, 0,isBtagSB,1);       //preselection: 1: SR,  3: AR
      //Mjj sideband
      if(checkMjjSB||checkMjjlowMETSB){
	//SS
	passAnySS(                      SRSS[ 4],SRSS[ 6],dummy,false,0,isBtagSB,true,1);   //full:         4: SR,  6: AR
	if(!skippreselection) passAnySS(SRSS[ 5],SRSS[ 7],dummy,true, 0,isBtagSB,true,1);   //preselection: 5: SR,  7: AR
      }

      float SF(0.), SFerr(0.), closureerr(0.);
      int idSS(-1), id3l(-1);
      for(int i = 0; i< nLlep(); ++i){
	if(lep_pass_VVV_cutbased_fo()[i]   &&!lep_pass_VVV_cutbased_tight()[i]   ) idSS = i;
	if(lep_pass_VVV_cutbased_3l_fo()[i]&&!lep_pass_VVV_cutbased_3l_tight()[i]) id3l = i;
      }
      int idboth = (nLlep()==2) ? idSS : id3l;
      if(!closuretest){
	float FR = ffwgt();
	float FRerr = fabs(ffwgt_up()-ffwgt());
	SF = FR/(1.-FR) * subtract;
	SFerr = FRerr/pow(1-FR,2) * subtract;
      } else {
	float FR = ffwgtqcd();
	float FRerr = fabs(ffwgtqcd_up()-ffwgtqcd());
	SF = FR/(1.-FR) * subtract;
	SFerr = FRerr/pow(1-FR,2) * subtract;
      }
      closureerr = (idboth>=0&&abs(lep_pdgId()[idboth])==11) ? 0.17 * SF : 0.31 * SF;
      float toterr = sqrt(pow(SFerr,2)+pow(closureerr,2));
      
      float fillid = (idboth>=0&&abs(lep_pdgId()[idboth])==11) ? 0. : 1.;
      if(SRSS[ 2]>=0&&fillid>=0){
	fillhisto(histos, "FakeEstimationEvsMu",       sample, sn, fillid, weight*SF);
	fillhisto(histos, "ARyieldEvsMu",              sample, sn, fillid, weight);
      }
      if(SRSS[ 3]>=0&&fillid>=0){
	fillhisto(histos, "PreselFakeEstimationEvsMu", sample, sn, fillid, weight*SF);
	fillhisto(histos, "PreselARyieldEvsMu",        sample, sn, fillid, weight);
      }
      if(isfake){
	int lepfakeidx = -1;
	for (int lepindex = 0;lepindex<nLlep();++lepindex){
	  if(lep_motherIdSS().at(lepindex) > 0) continue;
	  if(!lep_pass_VVV_cutbased_tight().at(lepindex)) continue;
	  else if(lep_motherIdSS().at(lepindex) == -3) continue;
	  lepfakeidx = lepindex;
	  break;
	}
	fillid = (lepfakeidx>=0&&abs(lep_pdgId()[lepfakeidx])==11) ? 0. : 1.;
	if(SRSS[ 0]>=0&&fillid>=0) histos["SRyieldEvsMu_"         +sn]->Fill(fillid,    weight);
	if(SRSS[ 1]>=0&&fillid>=0) histos["PreselSRyieldEvsMu_"   +sn]->Fill(fillid,    weight);
      }
      
      if(!(blindSR&&isData())){
	fillSRhisto(histos, "SRyield",                         sample, sn, sn, SRSS[ 0], SR3l[ 0], weight,                 weight);
	fillSRhisto(histos, "PreselSRyield",                   sample, sn, sn, SRSS[ 1], SR3l[ 1], weight,                 weight);
      }
      fillSRhisto(histos, "ARyield",                           sample, sn, sn, SRSS[ 2], SR3l[ 2], weight, weight);
      fillSRhisto(histos, "FakeEstimation",                    sample, sn, sn, SRSS[ 2], SR3l[ 2], weight* SF,             weight* SF);
      fillSRhisto(histos, "FakeEstimationFRup",                sample, sn, sn, SRSS[ 2], SR3l[ 2], weight*(SF+SFerr),      weight*(SF+SFerr));
      fillSRhisto(histos, "FakeEstimationFRdn",                sample, sn, sn, SRSS[ 2], SR3l[ 2], weight*(SF-SFerr),      weight*(SF-SFerr));
      fillSRhisto(histos, "FakeEstimationClosureup",           sample, sn, sn, SRSS[ 2], SR3l[ 2], weight*(SF+closureerr), weight*(SF+closureerr));
      fillSRhisto(histos, "FakeEstimationClosuredn",           sample, sn, sn, SRSS[ 2], SR3l[ 2], weight*(SF-closureerr), weight*(SF-closureerr));
      fillSRhisto(histos, "PreselARyield",                     sample, sn, sn, SRSS[ 3], SR3l[ 3], weight,                 weight);
      fillSRhisto(histos, "PreselFakeEstimation",              sample, sn, sn, SRSS[ 3], SR3l[ 3], weight* SF,             weight* SF);
      fillSRhisto(histos, "PreselFakeEstimationFRup",          sample, sn, sn, SRSS[ 3], SR3l[ 3], weight*(SF+SFerr),      weight*(SF+SFerr));
      fillSRhisto(histos, "PreselFakeEstimationFRdn",          sample, sn, sn, SRSS[ 3], SR3l[ 3], weight*(SF-SFerr),      weight*(SF-SFerr));
      fillSRhisto(histos, "PreselFakeEstimationClosureup",     sample, sn, sn, SRSS[ 3], SR3l[ 3], weight*(SF+closureerr), weight*(SF+closureerr));
      fillSRhisto(histos, "PreselFakeEstimationClosuredn",     sample, sn, sn, SRSS[ 3], SR3l[ 3], weight*(SF-closureerr), weight*(SF-closureerr));

      if(checkMjjSB){
	fillSRhisto(histos, "SRyield_Mjjsideband",                           sample, sn, sn, SRSS[ 4], SR3l[ 4], weight,                weight);
	fillSRhisto(histos, "PreselSRyield_Mjjsideband",                     sample, sn, sn, SRSS[ 5], SR3l[ 5], weight,                weight);
	fillSRhisto(histos, "ARyield_Mjjsideband",                           sample, sn, sn, SRSS[ 6], SR3l[ 6], weight,                weight);
	fillSRhisto(histos, "FakeEstimation_Mjjsideband",                    sample, sn, sn, SRSS[ 6], SR3l[ 6], weight* SF,            weight* SF);
	fillSRhisto(histos, "FakeEstimationUp_Mjjsideband",                  sample, sn, sn, SRSS[ 6], SR3l[ 6], weight*(SF+toterr),    weight*(SF+toterr));
	fillSRhisto(histos, "FakeEstimationDn_Mjjsideband",                  sample, sn, sn, SRSS[ 6], SR3l[ 6], weight*(SF-toterr),    weight*(SF-toterr));
	fillSRhisto(histos, "PreselARyield_Mjjsideband",                     sample, sn, sn, SRSS[ 7], SR3l[ 7], weight,                weight);
	fillSRhisto(histos, "PreselFakeEstimation_Mjjsideband",              sample, sn, sn, SRSS[ 7], SR3l[ 7], weight* SF,            weight* SF);
	fillSRhisto(histos, "PreselFakeEstimationUp_Mjjsideband",            sample, sn, sn, SRSS[ 7], SR3l[ 7], weight*(SF+toterr),    weight*(SF+toterr));
	fillSRhisto(histos, "PreselFakeEstimationDn_Mjjsideband",            sample, sn, sn, SRSS[ 7], SR3l[ 7], weight*(SF-toterr),    weight*(SF-toterr));
      }
      bool passMETsbloose = (fabs(Mjj()-80.)>=15.) && (met_pt()<40.);
      bool passMETsbtight = passMETsbloose && (MjjL()<400.) && (DetajjL()<1.5) && ( passSSem() ? (MllSS()>30.) : (MllSS()>40.)) && ( passSSem() ? (MTmax()>90.) : (true));
      if(checkMjjlowMETSB&&passMETsbtight){
	fillSRhisto(histos, "SRyield_MjjsidebandlowMET",                           sample, sn, sn, SRSS[ 5], SR3l[ 5], weight,                 weight);
	fillSRhisto(histos, "ARyield_MjjsidebandlowMET",                           sample, sn, sn, SRSS[ 7], SR3l[ 7], weight,                 weight);
	fillSRhisto(histos, "FakeEstimation_MjjsidebandlowMET",                    sample, sn, sn, SRSS[ 7], SR3l[ 7], weight* SF,             weight* SF);
	fillSRhisto(histos, "FakeEstimationUp_MjjsidebandlowMET",                  sample, sn, sn, SRSS[ 7], SR3l[ 7], weight*(SF+toterr),     weight*(SF+toterr));
	fillSRhisto(histos, "FakeEstimationDn_MjjsidebandlowMET",                  sample, sn, sn, SRSS[ 7], SR3l[ 7], weight*(SF-toterr),     weight*(SF-toterr));
      }
      if(checkMjjlowMETSB&&passMETsbloose){
	fillSRhisto(histos, "PreselSRyield_MjjsidebandlowMET",                     sample, sn, sn, SRSS[ 5], SR3l[ 5], weight,                 weight);
	fillSRhisto(histos, "PreselARyield_MjjsidebandlowMET",                     sample, sn, sn, SRSS[ 7], SR3l[ 7], weight,                 weight);
	fillSRhisto(histos, "PreselFakeEstimation_MjjsidebandlowMET",              sample, sn, sn, SRSS[ 7], SR3l[ 7], weight* SF,             weight* SF);
	fillSRhisto(histos, "PreselFakeEstimationUp_MjjsidebandlowMET",            sample, sn, sn, SRSS[ 7], SR3l[ 7], weight*(SF+toterr),     weight*(SF+toterr));
	fillSRhisto(histos, "PreselFakeEstimationDn_MjjsidebandlowMET",            sample, sn, sn, SRSS[ 7], SR3l[ 7], weight*(SF-toterr),     weight*(SF-toterr));
      }

      if(!isfake&&sample.find("WWW")==string::npos){ //skip fakes when filling in fake estimations from data

	fillSRhisto(histos, "SRyield",                                   sample, "fakesPred", "fakesPred", SRSS[ 2], SR3l[ 2], weight* SF,             weight* SF,             false);
	fillSRhisto(histos, "FakeEstimationFRup",                        sample, "fakesPred", "fakesPred", SRSS[ 2], SR3l[ 2], weight*(SF+SFerr),      weight*(SF+SFerr),      false);
	fillSRhisto(histos, "FakeEstimationFRdn",                        sample, "fakesPred", "fakesPred", SRSS[ 2], SR3l[ 2], weight*(SF-SFerr),      weight*(SF-SFerr),      false);
	fillSRhisto(histos, "FakeEstimationClosureup",                   sample, "fakesPred", "fakesPred", SRSS[ 2], SR3l[ 2], weight*(SF+closureerr), weight*(SF+closureerr), false);
	fillSRhisto(histos, "FakeEstimationClosuredn",                   sample, "fakesPred", "fakesPred", SRSS[ 2], SR3l[ 2], weight*(SF-closureerr), weight*(SF-closureerr), false);
	if(checkMjjSB){
	  fillSRhisto(histos, "SRyield_Mjjsideband",                       sample, "fakesPred", "fakesPred", SRSS[ 6], SR3l[ 6], weight* SF,             weight* SF,             false);
	  fillSRhisto(histos, "FakeEstimationUp_Mjjsideband",              sample, "fakesPred", "fakesPred", SRSS[ 6], SR3l[ 6], weight*(SF+toterr),     weight*(SF+toterr),     false);
	  fillSRhisto(histos, "FakeEstimationDn_Mjjsideband",              sample, "fakesPred", "fakesPred", SRSS[ 6], SR3l[ 6], weight*(SF-toterr),     weight*(SF-toterr),     false);
	}
	if(checkMjjlowMETSB&&passMETsbtight){
	  fillSRhisto(histos, "SRyield_MjjsidebandlowMET",               sample, "fakesPred", "fakesPred", SRSS[ 7], SR3l[ 7], weight* SF,             weight* SF,             false);
	  fillSRhisto(histos, "FakeEstimationUp_MjjsidebandlowMET",      sample, "fakesPred", "fakesPred", SRSS[ 7], SR3l[ 7], weight*(SF+toterr),     weight*(SF+toterr),     false);
	  fillSRhisto(histos, "FakeEstimationDn_MjjsidebandlowMET",      sample, "fakesPred", "fakesPred", SRSS[ 7], SR3l[ 7], weight*(SF-toterr),     weight*(SF-toterr),     false);
	}
	if(fillshapes){
	  if(SRSS[2]>=0) {
	    //fillhisto(histos, "Mjj_SRSS",                      sample, "fakesPred", Mjj(),            weight*SF,true);
	    fillhisto(histos, "njets_SRSS",                    sample, "fakesPred", nj30(),           weight*SF,true);
	    //fillhisto(histos, "nbjets_SRSS",                   sample, "fakesPred", nb(),             weight*SF,true);
	    fillhisto(histos, "pTl1_SRSS",                     sample, "fakesPred", lep_p4()[0].pt(), weight*SF,true);
	    fillhisto(histos, "pTl2_SRSS",                     sample, "fakesPred", lep_p4()[1].pt(), weight*SF,true);
	    fillhisto(histos, "Mll_SRSS",                      sample, "fakesPred", MllSS(),          weight*SF,true);
	    fillhisto(histos, "MTmax_SRSS",                    sample, "fakesPred", MTmax(),          weight*SF,true);
	    fillhisto(histos, "MET_SRSS",                      sample, "fakesPred", met_pt(),         weight*SF,true);
	  }
	  if(checkMjjSB&&SRSS[6]>=0) {
	    //fillhisto(histos, "Mjj_SRSS_Mjjsideband",          sample, "fakesPred", Mjj(),            weight*SF,true);
	    fillhisto(histos, "njets_SRSS_Mjjsideband",        sample, "fakesPred", nj30(),           weight*SF,true);
	    //fillhisto(histos, "nbjets_SRSS_Mjjsideband",       sample, "fakesPred", nb(),             weight*SF,true);
	    fillhisto(histos, "pTl1_SRSS_Mjjsideband",         sample, "fakesPred", lep_p4()[0].pt(), weight*SF,true);
	    fillhisto(histos, "pTl2_SRSS_Mjjsideband",         sample, "fakesPred", lep_p4()[1].pt(), weight*SF,true);
	    fillhisto(histos, "Mll_SRSS_Mjjsideband",          sample, "fakesPred", MllSS(),          weight*SF,true);
	    fillhisto(histos, "MTmax_SRSS_Mjjsideband",        sample, "fakesPred", MTmax(),          weight*SF,true);
	    fillhisto(histos, "MET_SRSS_Mjjsideband",          sample, "fakesPred", met_pt(),         weight*SF,true);
	  }
	  if(checkMjjlowMETSB&&SRSS[7]>=0&&passMETsbtight) {
	    //fillhisto(histos, "Mjj_SRSS_MjjsidebandlowMET",    sample, "fakesPred", Mjj(),            weight*SF,true);
	    fillhisto(histos, "njets_SRSS_MjjsidebandlowMET",  sample, "fakesPred", nj30(),           weight*SF,true);
	    //fillhisto(histos, "nbjets_SRSS_MjjsidebandlowMET", sample, "fakesPred", nb(),             weight*SF,true);
	    fillhisto(histos, "pTl1_SRSS_MjjsidebandlowMET",   sample, "fakesPred", lep_p4()[0].pt(), weight*SF,true);
	    fillhisto(histos, "pTl2_SRSS_MjjsidebandlowMET",   sample, "fakesPred", lep_p4()[1].pt(), weight*SF,true);
	    fillhisto(histos, "Mll_SRSS_MjjsidebandlowMET",    sample, "fakesPred", MllSS(),          weight*SF,true);
	    fillhisto(histos, "MTmax_SRSS_MjjsidebandlowMET",  sample, "fakesPred", MTmax(),          weight*SF,true);
	    //fillhisto(histos, "MET_SRSS_MjjsidebandlowMET",    sample, "fakesPred", met_pt(),         weight*SF,true);
	  }
	}
      }//note l.449: if(isfake||sample.find("WWW")!=string::npos)  continue; //skip fakes when filling in fake estimations from data
    }//event loop

    // Clean Up
    delete tree;
    file->Close();
    delete file;

  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  string savename = "rootfiles/FakeRateHistograms";
  if(isBtagSB)    savename = savename + "BTagVal";
  if(closuretest) savename = savename + "_closure.root";
  else            savename = savename + ".root";
  SaveHistosToFile(savename,histos,true,true,(chainnumber==0));
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
