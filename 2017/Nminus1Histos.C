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
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"
#include <sstream>
#include <iostream>
#include <fstream>

// CMS3
//#include "CMS3_old20150505.cc"
//#include "CMS3_fuckingsync.cc"
//#include "CMS3_Moriond17.cc"
#include "CMS3_WWW0116.cc"
#include "/home/users/haweber/CORE/Tools/MT2/MT2Utility.h"
#include "/home/users/haweber/CORE/Tools/MT2/MT2Utility.cc"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.h"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.cc"
#include "/home/users/haweber/CORE/Tools/goodrun.h"
#include "/home/users/haweber/CORE/Tools/goodrun.cc"
//MT2 variants

using namespace std;
using namespace tas;
//using namespace mt2_bisect;

int gentype_v2(unsigned lep1_index=0,unsigned lep2_index=1, int lep3_index=-1){
  bool gammafake = false;
  bool jetfake   = false;
  unsigned int ngenlep = ngenLepFromTau()+ ngenLep();
  unsigned int nW(0), nZ(0);
  bool lep1_real = lep_motherIdSS().at(lep1_index) > 0;
  bool lep2_real = lep_motherIdSS().at(lep2_index) > 0;
  bool lep3_real = false;
  if(lep3_index>0) lep3_real = lep_motherIdSS().at(lep3_index) > 0;
  vector<int> reallepindex;

  for (unsigned int lepindex = 0;lepindex<lep_p4().size();++lepindex){
      if(lep_motherIdSS().at(lepindex) > 0) reallepindex.push_back(lepindex);
      else if(lep_motherIdSS().at(lepindex) == -3) gammafake = true;
      else                                           jetfake = true;
      if(lep_isFromW().at(lepindex)) nW++;
      if(lep_isFromZ().at(lepindex)) nZ++;
  }
  //found two real leptons
  if(lep3_index<0){
    bool ischargeflip = false;
    bool isSS = false;
    if(lep1_real&&lep2_real) {
      int ilep1 =   lep_genPart_index().at(lep1_index);
      int ilep2 =   lep_genPart_index().at(lep2_index);
      bool lep1_chargeflip  =genPart_charge().at(ilep1)!= lep_charge().at(lep1_index);
      bool lep2_chargeflip  =genPart_charge().at(ilep2)!= lep_charge().at(lep2_index);
      if (!lep1_chargeflip&&!lep2_chargeflip&&nW==2) return 0; // true SS
      else if (!lep1_chargeflip&&!lep2_chargeflip) isSS = true; // true SS - but could be still lost lepton WZ
      if (lep1_chargeflip||lep2_chargeflip)   ischargeflip = true; // charge flip
    }
    
    if(ngenlep>2 || reallepindex.size()>2 || (nW>0 && nZ>0)) return 3; // lostlep
    if((ngenlep<2 ||!lep1_real||!lep2_real)&&    jetfake) return 4; // jetfake - if double fake with one jet fake and one gamma fake call it jet fake
    if((ngenlep<2 ||!lep1_real||!lep2_real)&&  gammafake) return 5; // gammafake
    if((ngenlep<2 ||!lep1_real||!lep2_real)&& !gammafake) return 4; // call all without gamma fake jetfake - safety cut
    if(isSS) return 0;
    if(ischargeflip) return 2;
    
    cout << "This event was not classified - 2 lepton event - v2" << endl;
    return 1;
  } else {
    //found three real leptons
    bool ischargeflip = false;
    bool isthreelep = false;
    if(lep1_real&&lep2_real&&lep3_real) {
      int ilep1 =   lep_genPart_index().at(lep1_index);
      int ilep2 =   lep_genPart_index().at(lep2_index);
      int ilep3 =   lep_genPart_index().at(lep3_index);
      bool lep1_chargeflip  =genPart_charge().at(ilep1)!= lep_charge().at(lep1_index);
      bool lep2_chargeflip  =genPart_charge().at(ilep2)!= lep_charge().at(lep2_index);
      bool lep3_chargeflip  =genPart_charge().at(ilep3)!= lep_charge().at(lep3_index);
      if (!lep1_chargeflip&&!lep2_chargeflip&&!lep3_chargeflip&&nW==3) return 0; // true WWW
      else if (!lep1_chargeflip&&!lep2_chargeflip&&!lep3_chargeflip) isthreelep = true; // true 3l, but could be lost lepton ZZ
      if (lep1_chargeflip||lep2_chargeflip||lep3_chargeflip)   ischargeflip = true; // charge flip
    }
    if(ngenlep>3 || reallepindex.size()>3 || (nW>=2 && nZ>=1) || (nZ>=3)) return 3; // lostlep (2 lep from W and 2 from Z, or 4 from Z)
    //there is the case of having WZZ with two lost leptons --> ngenlep>3 - correctly put has lostlep
    if((ngenlep<3 ||!lep1_real||!lep2_real||!lep3_real)&&    jetfake) return 4; // jetfake
    if((ngenlep<3 ||!lep1_real||!lep2_real||!lep3_real)&&  gammafake) return 5; // gammafake
    if((ngenlep<3 ||!lep1_real||!lep2_real||!lep3_real)&& !gammafake) return 4; // jetfake
    if(isthreelep) return 1;
    if(ischargeflip) return 2;
    
    cout << "This event was not classified - 3 lepton event - v2" << endl;
    return 0;
  }
}


float dR(LorentzVector vec1,LorentzVector vec2 ){                                                                                                              
  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);
}

float dEta(LorentzVector vec1,LorentzVector vec2 ){                                                                                                             
  return fabs(vec1.Eta() - vec2.Eta());
}

float dPhi(LorentzVector vec1,LorentzVector vec2 ){                                                                                                             
  return fabs(std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi())));
}

float deltaPhi(float phi1,float phi2 ){                                                                                                              
  return fabs(std::min(float(fabs(phi1-phi2)),float(2*M_PI-fabs(phi1-phi2))));
}

float mT(LorentzVector p4, LorentzVector met){
  float phi1 = p4.Phi();
  float phi2 = met.Phi();
  float Et1  = p4.Et();
  float Et2  = met.Et();

  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  vector<string> eventstocheck;

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");


  map<string, TH1D*> histos;
  vector<string> histonames; histonames.clear();
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();

  histonames.push_back("Detajj_SRlike_allSS");               hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Detajj_SRlike_SSee");                hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Detajj_SRlike_SSem");                hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Detajj_SRlike_SSmm");                hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Mjj_SRlike_allSS");                  hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("Mjj_SRlike_SSee");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("Mjj_SRlike_SSem");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("Mjj_SRlike_SSmm");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("MET_SRlike_allSS");                  hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("MET_SRlike_SSee");                   hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("MET_SRlike_SSem");                   hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("MET_SRlike_SSmm");                   hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("Mll_SRlike_allSS");                  hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("Mll_SRlike_SSee");                   hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("Mll_SRlike_SSem");                   hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("Mll_SRlike_SSmm");                   hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("MTmax_SRlike_SSem");                 hbins.push_back(14);hlow.push_back( 0);hup.push_back(210);//blind the data

  histonames.push_back("MET_SRlike_allSFOS");                 hbins.push_back( 8);hlow.push_back( 0);hup.push_back(160);//blind the data
  histonames.push_back("MET_SRlike_0SFOS");                   hbins.push_back( 8);hlow.push_back( 0);hup.push_back(160);//blind the data
  histonames.push_back("MET_SRlike_1SFOS");                   hbins.push_back(10);hlow.push_back( 0);hup.push_back(150);//blind the data
  histonames.push_back("MET_SRlike_2SFOS");                   hbins.push_back(10);hlow.push_back(10);hup.push_back(160);//blind the data
  histonames.push_back("MSFOS_SRlike_allSFOS");               hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("MSF_SRlike_0SFOS");                   hbins.push_back(15);hlow.push_back( 0);hup.push_back(150);//blind the data
  histonames.push_back("MSFOS_SRlike_1SFOS");                 hbins.push_back(14);hlow.push_back(22);hup.push_back(176);//blind the data
  histonames.push_back("MSFOS_SRlike_2SFOS");                 hbins.push_back(14);hlow.push_back(30);hup.push_back(170);//blind the data
  histonames.push_back("pTlll_SRlike_allSFOS");               hbins.push_back(12);hlow.push_back( 0);hup.push_back(120);//blind the data
  histonames.push_back("pTlll_SRlike_0SFOS");                 hbins.push_back(12);hlow.push_back( 0);hup.push_back(120);//blind the data
  histonames.push_back("pTlll_SRlike_1SFOS");                 hbins.push_back(12);hlow.push_back( 0);hup.push_back(120);//blind the data
  histonames.push_back("pTlll_SRlike_2SFOS");                 hbins.push_back(12);hlow.push_back(10);hup.push_back(120);//blind the data
  histonames.push_back("dPhiMETlll_SRlike_allSFOS");          hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data
  histonames.push_back("dPhiMETlll_SRlike_0SFOS");            hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data
  histonames.push_back("dPhiMETlll_SRlike_1SFOS");            hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data
  histonames.push_back("dPhiMETlll_SRlike_2SFOS");            hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data

  histonames.push_back("Detajj_SRloose_allSS");               hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Detajj_SRloose_SSee");                hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Detajj_SRloose_SSem");                hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Detajj_SRloose_SSmm");                hbins.push_back(15);hlow.push_back( 0);hup.push_back(4.5);//blind the data
  histonames.push_back("Mjj_SRloose_allSS");                  hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("Mjj_SRloose_SSee");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("Mjj_SRloose_SSem");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("Mjj_SRloose_SSmm");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("MET_SRloose_allSS");                  hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("MET_SRloose_SSee");                   hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("MET_SRloose_SSem");                   hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("MET_SRloose_SSmm");                   hbins.push_back(14);hlow.push_back( 0);hup.push_back(140);//blind the data
  histonames.push_back("Mll_SRloose_allSS");                  hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("Mll_SRloose_SSee");                   hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("Mll_SRloose_SSem");                   hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("Mll_SRloose_SSmm");                   hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("MTmax_SRloose_SSem");                 hbins.push_back(14);hlow.push_back( 0);hup.push_back(210);//blind the data

  histonames.push_back("MET_SRloose_allSFOS");                 hbins.push_back( 8);hlow.push_back( 0);hup.push_back(160);//blind the data
  histonames.push_back("MET_SRloose_0SFOS");                   hbins.push_back( 8);hlow.push_back( 0);hup.push_back(160);//blind the data
  histonames.push_back("MET_SRloose_1SFOS");                   hbins.push_back(10);hlow.push_back( 0);hup.push_back(150);//blind the data
  histonames.push_back("MET_SRloose_2SFOS");                   hbins.push_back(10);hlow.push_back(10);hup.push_back(160);//blind the data
  histonames.push_back("MSFOS_SRloose_allSFOS");               hbins.push_back(16);hlow.push_back(10);hup.push_back(170);//blind the data
  histonames.push_back("MSF_SRloose_0SFOS");                   hbins.push_back(15);hlow.push_back( 0);hup.push_back(150);//blind the data
  histonames.push_back("MSFOS_SRloose_1SFOS");                 hbins.push_back(14);hlow.push_back(22);hup.push_back(176);//blind the data
  histonames.push_back("MSFOS_SRloose_2SFOS");                 hbins.push_back(14);hlow.push_back(30);hup.push_back(170);//blind the data
  histonames.push_back("pTlll_SRloose_allSFOS");               hbins.push_back(12);hlow.push_back( 0);hup.push_back(120);//blind the data
  histonames.push_back("pTlll_SRloose_0SFOS");                 hbins.push_back(12);hlow.push_back( 0);hup.push_back(120);//blind the data
  histonames.push_back("pTlll_SRloose_1SFOS");                 hbins.push_back(12);hlow.push_back( 0);hup.push_back(120);//blind the data
  histonames.push_back("pTlll_SRloose_2SFOS");                 hbins.push_back(12);hlow.push_back(10);hup.push_back(120);//blind the data
  histonames.push_back("dPhiMETlll_SRloose_allSFOS");          hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data
  histonames.push_back("dPhiMETlll_SRloose_0SFOS");            hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data
  histonames.push_back("dPhiMETlll_SRloose_1SFOS");            hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data
  histonames.push_back("dPhiMETlll_SRloose_2SFOS");            hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);//blind the data

  //after finding good Mll region test with 2l OS region the in/out ratio
  
  //ee,em,mm,0SFOS,1SFOS,2SFOS

  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    string mapname = histonames[i];
    //cout << histonames[i] << endl;
    if(skimFilePrefix.find("Other")!=string::npos){
	mapname = histonames[i] + "_Other";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_WHtoWWW";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
    else if(skimFilePrefix.find("Background")!=string::npos){
      //cout << skimFilePrefix << endl;
      mapname = histonames[i] + "_trueSS";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_chargeflips";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_SSLL";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fakes";//jetfakes
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_photonfakes";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_others";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_trueWWW";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_3lLL";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_true3L";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
    else {
      mapname = histonames[i] + "_"+skimFilePrefix;
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      //histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
    }
  }
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->Sumw2(); h->second->SetDirectory(rootdir);
  }

  cout << "done" << endl;

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
      if(nlep()<2)               continue;
      if(firstgoodvertex()!=0)   continue;
      //if(nTaus20()!=0)           continue;
      //if(nBJetMedium()!=0)       continue;
      //if(nisoTrack_mt2_cleaned_VVV_cutbased_veto()!=0)  continue;

      //weight = 1;
      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("www_2l_mia")!=string::npos) weight *= 0.066805*91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")!=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(weight>100) cout << weight << " " << currentFile->GetTitle() << endl;
      if(isData()) weight = 1.;
      
      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());

      int nj(0),nb(0);
      int nj20(0),nj30(0);
      vector<int> i2p5;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	if(fabs(jets_p4()[n].Eta())<2.5) i2p5.push_back(n);
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<5) ++nj;
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.5) ++nj30;
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.5) ++nj20;
	if(jets_csv()[n]>0.5426&&jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.4) ++nb;
      }
      float minDR=999.;
      int jDR1(-1), jDR2(-1);
      for(unsigned int j1 = 0; j1<i2p5.size();++j1){
	if(jets_p4()[i2p5[j1] ].Pt()<30) continue;
	for(unsigned int j2 = j1+1; j2<i2p5.size();++j2){
	  if(jets_p4()[i2p5[j2] ].Pt()<30) continue;
	  if(dR(jets_p4()[i2p5[j1] ], jets_p4()[i2p5[j2] ])<minDR){
	    minDR = dR(jets_p4()[i2p5[j1] ], jets_p4()[i2p5[j2] ]);
	    jDR1 = i2p5[j1]; jDR2 = i2p5[j2];
	  }
	}
      }
      float Mjj = -1;
      if(jDR1>=0&&jDR2>=0) Mjj = (jets_p4()[jDR1]+jets_p4()[jDR2]).M();
      float MjjL = -1; float Detajj = -1;
      if(i2p5.size()>1&&jets_p4()[i2p5[0] ].Pt()>30&&jets_p4()[i2p5[1] ].Pt()>30) {
	MjjL = (jets_p4()[i2p5[0] ]+jets_p4()[i2p5[1] ]).M();
	Detajj = dEta(jets_p4()[i2p5[0] ],jets_p4()[i2p5[1] ]);
      }
      //if(jets_p4().size()>1&&njets()>1&&Mjj>0&&fabs(mjj()-Mjj)>0.01) cout << "MJJ " << Mjj << " " << mjj() << endl;
      //if(jets_p4().size()>1&&njets()>1&&Detajj>0&&fabs(deta_jj()-Detajj)>0.01) cout << "Deta " << Detajj << " " << deta_jj() << endl;
      bool passMDetajj = true;
      if(nj30<2) passMDetajj = false;
      if(fabs(Detajj)>1.5)  passMDetajj = false;
      if(fabs(MjjL)>400.)   passMDetajj = false;
      if(fabs(Mjj-80.)>20.) passMDetajj = false;
      bool passMjj = true;
      if(nj30<2) passMjj = false;
      if(fabs(Mjj-80.)>20.) passMjj = false;
      bool passDetajj = true;
      if(nj30<2) passDetajj = false;
      if(fabs(Detajj)>1.5)  passDetajj = false;
      if(fabs(MjjL)>400.)   passDetajj = false;
      bool passDetajj2 = true;
      if(nj30<2) passDetajj2 = false;
      if(fabs(Detajj)>1.5)  passDetajj2 = false;
      vector<int> vSS, v3l, v, iSS, i3l;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool isSS = false; bool is3l = false;
	bool isaSS = false; bool isa3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v2()[i]){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06&&lep_isTriggerSafe_v2()[i]){
	    if(lep_p4()[i ].Pt()>20) { i3l.push_back(i); is3l = true; }
	    if(lep_p4()[i ].Pt()>30) { iSS.push_back(i); isSS = true; }
	  }
	}
      }

      int nvetoSS = vSS.size();
      int nveto3l = v3l.size();
      int nveto = v.size();
      int nSS = iSS.size();
      int n3l = i3l.size();

      bool passofflineforTrigger = false;
      int nel25 = 0;
      int nel = 0;
      int nmu25 = 0;
      int nmu = 0;
      for(int i = 0; i<n3l; ++i){
	if(abs(lep_pdgId()[i3l[i] ])==11){
	  ++nel;
	  if(lep_p4()[i3l[i] ].Pt()>25) ++nel25;
	} else if(abs(lep_pdgId()[i3l[i] ])==13){
	  ++nmu;
	  if(lep_p4()[i3l[i] ].Pt()>25) ++nmu25;
	}
      }
      if(nmu>=2)           passofflineforTrigger = true;
      if(nmu25>=1&&nel>=1) passofflineforTrigger = true;
      if(nel25>=1&&nel>=2) passofflineforTrigger = true;

      if(n3l<2) continue;
      if(!passofflineforTrigger) continue;
      
      if(isData()){
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) continue;
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = false;
	if(nmu>=2&&HLT_DoubleMu() )                               passonlineTrigger = true;
	if(nmu25>=1&&nel>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) passonlineTrigger = true;
	if(!passonlineTrigger) continue;
      }


      string sn = skimFilePrefix;
      string sn2 = skimFilePrefix;
      if(sn.find("Other")!=string::npos){
	bool isHtoWW = false;
	bool isWnotFromH = false;
	for(unsigned int i = 0; i<genPart_pdgId().size();++i){
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])==25&&genPart_status()[i]==22) { isHtoWW = true; }
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])!=25&&abs(genPart_grandmaId()[i])!=25&&genPart_status()[i]==22) { isWnotFromH = true; }
	  if(isHtoWW&&isWnotFromH) break;
	}
	bool isWHtoWWW = isHtoWW&&isWnotFromH;
	if(isWHtoWWW) {
	  sn = "WHtoWWW";
	  sn2 = "WHtoWWW";
	}
      }
      if(sn.find("Background")!=string::npos){
	bool isHtoWW = false;
	bool isWnotFromH = false;
	for(unsigned int i = 0; i<genPart_pdgId().size();++i){
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])==25&&genPart_status()[i]==22) { isHtoWW = true; }
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])!=25&&abs(genPart_grandmaId()[i])!=25&&genPart_status()[i]==22) { isWnotFromH = true; }
	  if(isHtoWW&&isWnotFromH) break;
	}
	if(isHtoWW&&isWnotFromH) continue;
	if((iSS.size())>=2){
	  int l1(-1),l2(-1);
	  if(iSS.size()>=2) { l1 = iSS[0]; l2 = iSS[1]; }
	  //else if(iSS.size()==1&&iaSS.size()>=1) { l1 = iSS[0]; l2 = iaSS[0]; }
	  //else if(iaSS.size()>=2) { l1 = iaSS[0]; l2 = iaSS[1]; }
	  int gentype = gentype_v2(l1,l2,-1);
	  if(     gentype==0) sn = "trueSS";
	  else if(gentype==2) sn = "chargeflips";
	  else if(gentype==3) sn = "SSLL";
	  else if(gentype==4) sn = "fakes";
	  else if(gentype==5) sn = "photonfakes";
	  else                sn = "others";
	}
	if(i3l.size()>=3){
	  int l1(-1), l2(-1), l3(-1);
	  if(i3l.size()>=3) { l1 = i3l[0]; l2 = i3l[1]; l3 = i3l[2]; }
	  int gentype = gentype_v2(l1,l2,l3);
	  if(     gentype==0) sn2 = "trueWWW";
	  else if(gentype==1) sn2 = "true3L";
	  else if(gentype==2) sn2 = "chargeflips";
	  else if(gentype==3) sn2 = "3lLL";
	  else if(gentype==4) sn2 = "fakes";
	  else if(gentype==5) sn2 = "photonfakes";
	  else                sn2 = "others";
	}
      }
      
      float MTmax = -1;
      if(iSS.size()>=2){
	if(mT(lep_p4()[iSS[1] ],MET)>mT(lep_p4()[iSS[0] ],MET)) MTmax = mT(lep_p4()[iSS[1] ],MET);
	else MTmax = mT(lep_p4()[iSS[0] ],MET);
      }

      double MllSS = -1;
      bool ee[20],mm[20],em[20];
      for(int i = 0; i<20; ++i) { ee[i] = false; em[i] = false; mm[i] = false; }
      //0: SR, 5: SR but no Mjj cut
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	MllSS = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11){
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.&&passMDetajj) ee[0] = true; //all
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.&&passMDetajj)               ee[1] = true; //all but MET
	  if(met_pt()>40.&&passMDetajj)                                                                                                         ee[2] = true; //all but Mll
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.&&passDetajj)  ee[3] = true; //all but Mjj
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.&&passMjj)     ee[5] = true; //all but Detajj
	  ee[6] = true; //none
	  if(fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[7] = true; //none but MllZ

	}
	if((abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11)||(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13)){
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&met_pt()>40.&&MTmax>90.&&passMDetajj) em[0] = true;//all
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.&&passMDetajj)               em[1] = true;//all but MET
	  if(met_pt()>40.&&MTmax>90.&&passMDetajj)                                                em[2] = true;//all but Mll
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&met_pt()>40.&&passDetajj&&MTmax>90.)  em[3] = true;//all but Mjj
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&met_pt()>40.&&passMDetajj)            em[4] = true;//all but MTmax
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&met_pt()>40.&&passMjj&&MTmax>90.)     em[5] = true;//all but Detajj
	  em[6] = true;//none
	  em[7] = true;//none

	}
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13){
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj) mm[0] = true;//all
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj) mm[1] = true;//all 'but MET'
	  if(passMDetajj)                                                mm[2] = true;//all but Mll
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passDetajj)  mm[3] = true;//all but Mjj
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMjj)     mm[5] = true;//all but Detajj
	  mm[6] = true;//none
	  mm[7] = true;//none

	}
      }

      int SFOS[20];
      for(int i = 0; i<20; ++i) { SFOS[i] = -1; }
      //0 : SR, 1: inverted Mll-Z, 2: SR but inverted DPhi,Pt, 3: inverted Mll-Z and inverted DPhi,Pt
      double DPhi(-1), pT(-1);
      if(nj<2&&nb==0/*&&nveto3l==0*/&&n3l==3){
	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	int SFOScounter = 0;
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	bool pass0X(false), pass1X(false), pass2X(false);
	bool pass0SFcut = true;

	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90)<10.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0SFcut = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0SFcut = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0SFcut = false;
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  //if(!atleastoneSFOSZ) pass0X = false;
	}
	if(SFOScounter==1){
	  pass1 = true;
	  pass1X = true;
	  //if(met_pt()<45) pass1=false;
	  //if(met_pt()<45) pass1X=false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  //if(met_pt()<55) pass2=false;
	  //if(met_pt()<55) pass2X=false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	DPhi = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	pT = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	if(pass0&& pass0SFcut&&DPhi >2.7&&pT >60.) SFOS[0] = 0;//all
	if(pass0&& pass0SFcut&&DPhi >2.7&&pT >60.) SFOS[1] = 0;//all 'but MET'
	if(pass0&&!pass0SFcut&&DPhi >2.7&&pT >60.) SFOS[2] = 0;//all but inverted MSFOS cut (MSF for 0SFOS)
	if(pass0&& pass0SFcut&&DPhi >2.7         ) SFOS[3] = 0;//all but pT cut
	if(pass0&& pass0SFcut           &&pT >60.) SFOS[4] = 0;//all but dPhi cut
	if(pass0 ) SFOS[5] = 0;//none
	if(pass0X) SFOS[6] = 0;//XXX

	if(pass1 &&met_pt()>45.&&DPhi >2.5&&pT >60.) SFOS[0] = 1;//all
	if(pass1               &&DPhi >2.5&&pT >60.) SFOS[1] = 1;//all but MET cut
	if(pass1X&&met_pt()>45.&&DPhi >2.5&&pT >60.) SFOS[2] = 1;//all but invertedMSFOS cut (MSF for 0SFOS)
	if(pass1 &&met_pt()>45.&&DPhi >2.5         ) SFOS[3] = 1;//all but pT cut
	if(pass1 &&met_pt()>45.           &&pT >60.) SFOS[4] = 1;//all but dPhi cut
	if(pass1 ) SFOS[5] = 1;//none but Mll
	if(pass1X) SFOS[6] = 1;//none and inverted Mll

	if(pass2 &&met_pt()>55.&&DPhi >2.5&&pT >60.) SFOS[0] = 2;//all
	if(pass2               &&DPhi >2.5&&pT >60.) SFOS[1] = 2;//all but MET cut
	if(pass2X&&met_pt()>55.&&DPhi >2.5&&pT >60.) SFOS[2] = 2;//all but invertedMSFOS cut (MSF for 0SFOS)
	if(pass2 &&met_pt()>55.&&DPhi >2.5         ) SFOS[3] = 2;//all but pT cut
	if(pass2 &&met_pt()>55.           &&pT >60.) SFOS[4] = 2;//all but dPhi cut
	if(pass2 ) SFOS[5] = 2;//none but Mll
	if(pass2X) SFOS[6] = 2;//none and inverted Mll
      }

      if(!isData()){
	if(ee[5]&&ee[0]){ histos["Detajj_SRlike_SSee_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRlike_allSS_"+ sn]->Fill(Detajj,  weight); }
	if(em[5]&&em[0]){ histos["Detajj_SRlike_SSem_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRlike_allSS_"+ sn]->Fill(Detajj,  weight); }
	if(mm[5]&&mm[0]){ histos["Detajj_SRlike_SSmm_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRlike_allSS_"+ sn]->Fill(Detajj,  weight); }
	if(ee[3]&&ee[0]){ histos["Mjj_SRlike_SSee_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRlike_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(em[3]&&em[0]){ histos["Mjj_SRlike_SSem_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRlike_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(mm[3]&&mm[0]){ histos["Mjj_SRlike_SSmm_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRlike_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(ee[1]&&ee[0]){ histos["MET_SRlike_SSee_"     + sn]->Fill(met_pt(),weight); histos["MET_SRlike_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(em[1]&&em[0]){ histos["MET_SRlike_SSem_"     + sn]->Fill(met_pt(),weight); histos["MET_SRlike_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(mm[1]&&mm[0]){ histos["MET_SRlike_SSmm_"     + sn]->Fill(met_pt(),weight); histos["MET_SRlike_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(ee[2]&&ee[0]){ histos["Mll_SRlike_SSee_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRlike_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(em[2]&&em[0]){ histos["Mll_SRlike_SSem_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRlike_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(mm[2]&&mm[0]){ histos["Mll_SRlike_SSmm_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRlike_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(em[4]&&em[0]){ histos["MTmax_SRlike_SSem_"   + sn]->Fill(MTmax,   weight);  }
      }
	if(ee[5]&&!ee[0]){ histos["Detajj_SRlike_SSee_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRlike_allSS_"+ sn]->Fill(Detajj,  weight); }
	if(em[5]&&!em[0]){ histos["Detajj_SRlike_SSem_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRlike_allSS_"+ sn]->Fill(Detajj,  weight); }
	if(mm[5]&&!mm[0]){ histos["Detajj_SRlike_SSmm_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRlike_allSS_"+ sn]->Fill(Detajj,  weight); }
      	if(ee[3]&&!ee[0]){ histos["Mjj_SRlike_SSee_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRlike_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(em[3]&&!em[0]){ histos["Mjj_SRlike_SSem_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRlike_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(mm[3]&&!mm[0]){ histos["Mjj_SRlike_SSmm_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRlike_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(ee[1]&&!ee[0]){ histos["MET_SRlike_SSee_"     + sn]->Fill(met_pt(),weight); histos["MET_SRlike_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(em[1]&&!em[0]){ histos["MET_SRlike_SSem_"     + sn]->Fill(met_pt(),weight); histos["MET_SRlike_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(mm[1]&&!mm[0]){ histos["MET_SRlike_SSmm_"     + sn]->Fill(met_pt(),weight); histos["MET_SRlike_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(ee[2]&&!ee[0]){ histos["Mll_SRlike_SSee_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRlike_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(em[2]&&!em[0]){ histos["Mll_SRlike_SSem_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRlike_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(mm[2]&&!mm[0]){ histos["Mll_SRlike_SSmm_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRlike_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(em[4]&&!em[0]){ histos["MTmax_SRlike_SSem_"   + sn]->Fill(MTmax,   weight);  }

	if(ee[7]){ histos["Detajj_SRloose_SSee_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRloose_allSS_"+ sn]->Fill(Detajj,  weight); }
	if(em[6]){ histos["Detajj_SRloose_SSem_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRloose_allSS_"+ sn]->Fill(Detajj,  weight); }
	if(mm[6]){ histos["Detajj_SRloose_SSmm_"  + sn]->Fill(Detajj,  weight); histos["Detajj_SRloose_allSS_"+ sn]->Fill(Detajj,  weight); }
      	if(ee[7]){ histos["Mjj_SRloose_SSee_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRloose_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(em[6]){ histos["Mjj_SRloose_SSem_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRloose_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(mm[6]){ histos["Mjj_SRloose_SSmm_"     + sn]->Fill(Mjj,     weight); histos["Mjj_SRloose_allSS_"   + sn]->Fill(Mjj,     weight); }
	if(ee[7]){ histos["MET_SRloose_SSee_"     + sn]->Fill(met_pt(),weight); histos["MET_SRloose_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(em[6]){ histos["MET_SRloose_SSem_"     + sn]->Fill(met_pt(),weight); histos["MET_SRloose_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(mm[6]){ histos["MET_SRloose_SSmm_"     + sn]->Fill(met_pt(),weight); histos["MET_SRloose_allSS_"   + sn]->Fill(met_pt(),weight); }
	if(ee[6]){ histos["Mll_SRloose_SSee_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRloose_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(em[6]){ histos["Mll_SRloose_SSem_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRloose_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(mm[6]){ histos["Mll_SRloose_SSmm_"     + sn]->Fill(MllSS,   weight); histos["Mll_SRloose_allSS_"   + sn]->Fill(MllSS,   weight); }
	if(em[6]){ histos["MTmax_SRloose_SSem_"   + sn]->Fill(MTmax,   weight);  }

      if(n3l>=3){
      	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	if(!isData()){
	  if(SFOS[0]==0){
	    if(SF01) histos["MSF_SRlike_0SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	    if(SF02) histos["MSF_SRlike_0SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	    if(SF12) histos["MSF_SRlike_0SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	  }
	  if(SFOS[0]==1){
	    if(OS01&&SF01) histos["MSFOS_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	    if(OS02&&SF02) histos["MSFOS_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	    if(OS12&&SF12) histos["MSFOS_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	  }
	  if(SFOS[0]==2){
	    if(OS01&&SF01) histos["MSFOS_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	    if(OS02&&SF02) histos["MSFOS_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	    if(OS12&&SF12) histos["MSFOS_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	  }
	  if(SFOS[0]>=1){
	    if(OS01&&SF01) histos["MSFOS_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	    if(OS02&&SF02) histos["MSFOS_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	    if(OS12&&SF12) histos["MSFOS_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	  }
	  if(SFOS[1]==0&&SFOS[0]==0){ histos["MET_SRlike_0SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRlike_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	  if(SFOS[3]==0&&SFOS[0]==0){ histos["pTlll_SRlike_0SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRlike_allSFOS_"     +sn2]->Fill(pT,      weight); }
	  if(SFOS[4]==0&&SFOS[0]==0){ histos["dPhiMETlll_SRlike_0SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRlike_allSFOS_"+sn2]->Fill(DPhi,    weight); }
	  if(SFOS[1]==1&&SFOS[0]==1){ histos["MET_SRlike_1SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRlike_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	  if(SFOS[3]==1&&SFOS[0]==1){ histos["pTlll_SRlike_1SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRlike_allSFOS_"     +sn2]->Fill(pT,      weight); }
	  if(SFOS[4]==1&&SFOS[0]==1){ histos["dPhiMETlll_SRlike_1SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRlike_allSFOS_"+sn2]->Fill(DPhi,    weight); }
	  if(SFOS[1]==2&&SFOS[0]==2){ histos["MET_SRlike_2SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRlike_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	  if(SFOS[3]==2&&SFOS[0]==2){ histos["pTlll_SRlike_2SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRlike_allSFOS_"     +sn2]->Fill(pT,      weight); }
	  if(SFOS[4]==2&&SFOS[0]==2){ histos["dPhiMETlll_SRlike_2SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRlike_allSFOS_"+sn2]->Fill(DPhi,    weight); }
	}//no data
	if(SFOS[2]==0){
	  if(SF01) histos["MSF_SRlike_0SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(SF02) histos["MSF_SRlike_0SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(SF12) histos["MSF_SRlike_0SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	if(SFOS[2]==1){
	  if(OS01&&SF01) histos["MSFOS_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(OS02&&SF02) histos["MSFOS_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(OS12&&SF12) histos["MSFOS_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	if(SFOS[2]==2){
	  if(OS01&&SF01) histos["MSFOS_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(OS02&&SF02) histos["MSFOS_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(OS12&&SF12) histos["MSFOS_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	if(SFOS[2]>=1){
	  if(OS01&&SF01) histos["MSFOS_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(OS02&&SF02) histos["MSFOS_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(OS12&&SF12) histos["MSFOS_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	if(SFOS[1]==0&&SFOS[0]<0){ histos["MET_SRlike_0SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRlike_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	if(SFOS[3]==0&&SFOS[0]<0){ histos["pTlll_SRlike_0SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRlike_allSFOS_"     +sn2]->Fill(pT,      weight); }
	if(SFOS[4]==0&&SFOS[0]<0){ histos["dPhiMETlll_SRlike_0SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRlike_allSFOS_"+sn2]->Fill(DPhi,    weight); }
	if(SFOS[1]==1&&SFOS[0]<0){ histos["MET_SRlike_1SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRlike_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	if(SFOS[3]==1&&SFOS[0]<0){ histos["pTlll_SRlike_1SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRlike_allSFOS_"     +sn2]->Fill(pT,      weight); }
	if(SFOS[4]==1&&SFOS[0]<0){ histos["dPhiMETlll_SRlike_1SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRlike_allSFOS_"+sn2]->Fill(DPhi,    weight); }
	if(SFOS[1]==2&&SFOS[0]<0){ histos["MET_SRlike_2SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRlike_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	if(SFOS[3]==2&&SFOS[0]<0){ histos["pTlll_SRlike_2SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRlike_allSFOS_"     +sn2]->Fill(pT,      weight); }
	if(SFOS[4]==2&&SFOS[0]<0){ histos["dPhiMETlll_SRlike_2SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRlike_allSFOS_"+sn2]->Fill(DPhi,    weight); }

	if(SFOS[5]==0||SFOS[6]==0){
	  if(SF01) histos["MSF_SRloose_0SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(SF02) histos["MSF_SRloose_0SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(SF12) histos["MSF_SRloose_0SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	if(SFOS[5]==1||SFOS[6]==1){
	  if(OS01&&SF01) histos["MSFOS_SRloose_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(OS02&&SF02) histos["MSFOS_SRloose_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(OS12&&SF12) histos["MSFOS_SRloose_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	if(SFOS[5]==2||SFOS[6]==2){
	  if(OS01&&SF01) histos["MSFOS_SRloose_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(OS02&&SF02) histos["MSFOS_SRloose_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(OS12&&SF12) histos["MSFOS_SRloose_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	  if((SFOS[5]==1||SFOS[6]==1)||(SFOS[5]==2||SFOS[6]==2)){
	  if(OS01&&SF01) histos["MSFOS_SRloose_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(OS02&&SF02) histos["MSFOS_SRloose_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(OS12&&SF12) histos["MSFOS_SRloose_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
	if(SFOS[5]==0){ histos["MET_SRloose_0SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRloose_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	if(SFOS[5]==0){ histos["pTlll_SRloose_0SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRloose_allSFOS_"     +sn2]->Fill(pT,      weight); }
	if(SFOS[5]==0){ histos["dPhiMETlll_SRloose_0SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRloose_allSFOS_"+sn2]->Fill(DPhi,    weight); }
	if(SFOS[5]==1){ histos["MET_SRloose_1SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRloose_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	if(SFOS[5]==1){ histos["pTlll_SRloose_1SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRloose_allSFOS_"     +sn2]->Fill(pT,      weight); }
	if(SFOS[5]==1){ histos["dPhiMETlll_SRloose_1SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRloose_allSFOS_"+sn2]->Fill(DPhi,    weight); }
	if(SFOS[5]==2){ histos["MET_SRloose_2SFOS_"       +sn2]->Fill(met_pt(),weight); histos["MET_SRloose_allSFOS_"       +sn2]->Fill(met_pt(),weight); }
	if(SFOS[5]==2){ histos["pTlll_SRloose_2SFOS_"     +sn2]->Fill(pT,      weight); histos["pTlll_SRloose_allSFOS_"     +sn2]->Fill(pT,      weight); }
	if(SFOS[5]==2){ histos["dPhiMETlll_SRloose_2SFOS_"+sn2]->Fill(DPhi,    weight); histos["dPhiMETlll_SRloose_allSFOS_"+sn2]->Fill(DPhi,    weight); }

      }//n3l>2
    }//event loop
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //add overflow
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
    //add underflow
    h->second->SetBinContent(1, h->second->GetBinContent(1)+ h->second->GetBinContent(0) );
    h->second->SetBinError(1, sqrt(pow(h->second->GetBinError(1),2)+pow(h->second->GetBinError(0),2) ) );
  }
  
  string filename = "rootfiles/Nminus1Histos.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //string mapname = histonames[i]+"_"+skimFilePrefix;
    //string writename = histonames[i];
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;

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
