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
//#include "CMS3_old20150505.cc"
//#include "CMS3_fuckingsync.cc"
//#include "CMS3_Moriond17.cc"
#include "CMS3_WWW0111.cc"
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


int gentype(unsigned lep1_index,unsigned lep2_index){
  bool gammafake = false;
  unsigned int ngenlep = ngenLepFromTau()+ ngenLep();
  unsigned int nW(0), nZ(0);
  bool lep1_real = lep_motherIdSS().at(lep1_index) > 0;
  bool lep2_real = lep_motherIdSS().at(lep2_index) > 0;
  vector<int> reallepindex;

  for (unsigned int lepindex = 0;lepindex<lep_p4().size();++lepindex){
      if(lep_motherIdSS().at(lepindex) > 0) reallepindex.push_back(lepindex);
      if(lep_motherIdSS().at(lepindex) == -3) gammafake = true;
      if(lep_isFromW().at(lepindex)) nW++;
      if(lep_isFromZ().at(lepindex)) nZ++;
  }

  if(ngenlep>2 || reallepindex.size()>2 || (nW>0 && nZ>0)) return 3; // lostlep
  if((ngenlep<2 ||!lep1_real||!lep2_real)&& !gammafake) return 4; // jetfake
  if((ngenlep<2 ||!lep1_real||!lep2_real)&& gammafake) return 5; // gammafake

//found two real leptons
  if(lep1_real&&lep2_real) {
      int ilep1 =   lep_genPart_index().at(lep1_index);
      int ilep2 =   lep_genPart_index().at(lep2_index);
      bool lep1_chargeflip  =genPart_charge().at(ilep1)!= lep_charge().at(lep1_index);
      bool lep2_chargeflip  =genPart_charge().at(ilep2)!= lep_charge().at(lep2_index);
      if (!lep1_chargeflip&&!lep2_chargeflip) return 1; // true SS
      if (lep1_chargeflip||lep2_chargeflip)   return 2; // charge flip
   }

  cout << "This event was not classified - v1" << endl;
  return 1;
 }

int gentype_v2(int lep1_index=0,int lep2_index=1, int lep3_index=-1){
  //cout << "Calling " << lep1_index << " " << lep2_index << " " << lep3_index << endl;
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

float coneCorrPt(int lepi){//apply only to loose but not tight leptons
        float coneptcorr = 0;
        float relIso = lep_relIso03EAv2().at(lepi);
        //bool passId = isGoodLepton(lepi, "ss");
        if ( abs(lep_pdgId().at(lepi)) == 11 )
        {
            if ( fabs(lep_eta().at(lepi)) <= 1.479 )
                coneptcorr = std::max( 0., relIso - 0.0588 );
            else
                coneptcorr = std::max( 0., relIso - 0.0571 );
        }

        if ( abs(lep_pdgId().at(lepi)) == 13 )
            coneptcorr = std::max( 0., relIso - 0.06 );

        // If the lepton passed the tight ID no correction is needed.
        //if ( passId )  coneptcorr = 0;
 return coneptcorr;
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
  std::ostringstream*  eventstocheckEE = 0;
  std::ostringstream*  eventstocheckEM = 0;
  std::ostringstream*  eventstocheckMM = 0;
  std::ostringstream*  eventstocheck0SFOS = 0;
  std::ostringstream*  eventstocheck1SFOS = 0;
  std::ostringstream*  eventstocheck2SFOS = 0;
  eventstocheckEE   = new std::ostringstream();
  eventstocheckEM   = new std::ostringstream();
  eventstocheckMM   = new std::ostringstream();
  eventstocheck0SFOS = new std::ostringstream();
  eventstocheck1SFOS = new std::ostringstream();
  eventstocheck2SFOS = new std::ostringstream();
  vector<unsigned int> r,l,e;
r.push_back(276502);l.push_back( 371);e.push_back( 565320149);
  vector<unsigned int> rr,ll,ee;
  rr.push_back(0);ll.push_back(0);ee.push_back(0);

  
  float wgtsum1(0), wgtsum2(0), wgtsum3(0);
  int counter = 0;
  bool applylepSF = false;

  TFile *fSF = new TFile("rootfiles/SF_TnP.root","read");
  TH2F *hMu = (TH2F*)fSF->Get("muSF");
  TH2F *hElReco = (TH2F*)fSF->Get("elSF_reco");
  TH2F *hElID = (TH2F*)fSF->Get("elSF_ID");
  float muptmin = 20.1;                           float muptmax = 199.9; float muetamin =  0.01; float muetamax = 2.49;
  float elptmin = 10.1; float elptminReco = 25.1; float elptmax = 499.9; float eletamin = -2.49; float eletamax = 2.49;

    
  TFile *fFR = new TFile("rootfiles/fakerate_pt_v_eta_20170911.root","read");
  TH2D *hMuFR = (TH2D*)fFR->Get("muon_fakerate_conecorrpt_v_eta");
  TH2D *hElFR = (TH2D*)fFR->Get("elec_fakerate_conecorrpt_v_eta");
  double muptminFR = 10.1; double muptmaxFR = 119.9; double muetaminFR =  0.01; double muetamaxFR = 2.39;
  double elptminFR = 10.1; double elptmaxFR = 119.9; double eletaminFR =  0.01; double eletamaxFR = 2.49;
   
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
  histonames.push_back("SRyield");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselSRyield");              hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  //use cone correction only to extract fake rate
  histonames.push_back("ARyield");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimation");             hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRup");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRdn");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  //do not use cone corrections to extract fake rate
  histonames.push_back("NoConeCorrARyield");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("NoConeCorrFakeEstimation");             hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("NoConeCorrFakeEstimationFRup");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("NoConeCorrFakeEstimationFRdn");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  //use cone correction to also check lepton pT thresholds
  histonames.push_back("AllConeCorrARyield");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("AllConeCorrFakeEstimation");             hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("AllConeCorrFakeEstimationFRup");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("AllConeCorrFakeEstimationFRdn");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  //use cone correction only to extract fake rate - preselection
  histonames.push_back("PreselARyield");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimation");             hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRup");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRdn");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);

  histonames.push_back("SRyield_Mjjsideband");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselSRyield_Mjjsideband");              hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ARyield_Mjjsideband");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimation_Mjjsideband");             hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRup_Mjjsideband");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRdn_Mjjsideband");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselARyield_Mjjsideband");              hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimation_Mjjsideband");       hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRup_Mjjsideband");   hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRdn_Mjjsideband");   hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("SRyield_MjjsidebandlowMET");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselSRyield_MjjsidebandlowMET");              hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("ARyield_MjjsidebandlowMET");                    hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimation_MjjsidebandlowMET");             hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRup_MjjsidebandlowMET");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("FakeEstimationFRdn_MjjsidebandlowMET");         hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselARyield_MjjsidebandlowMET");              hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimation_MjjsidebandlowMET");       hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRup_MjjsidebandlowMET");   hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  histonames.push_back("PreselFakeEstimationFRdn_MjjsidebandlowMET");   hbins.push_back( 6); hlow.push_back(    0); hup.push_back(6);
  
  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    string mapname = histonames[i];
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
    } else {
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

      bool testevent = false;
      for(unsigned int n = 0; n<r.size(); ++n){
	if(tas::run()==r[n]&&tas::lumi()==l[n]&&tas::evt()==e[n]){
	  testevent = true;
	  break;
	}
      }
      if(testevent) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " test this event" << endl;
      bool checkevent = false;
      for(unsigned int n = 0; n<rr.size(); ++n){
	if(tas::run()==rr[n]&&tas::lumi()==ll[n]&&tas::evt()==ee[n]){
	  checkevent = true;
	  break;
	}
      }
    
      // Progress
      CMS3::progress( nEventsTotal, nEventsChain );

      
      double weight = evt_scale1fb()*35.9;
      //if(nEventsTotal==0) cout << weight << endl; 

      if(nVert()<0)              continue;
      if(nlep()<2)               continue;

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
      bool passMDetajj = true;
      if(nj30<2)            passMDetajj = false;
      if(fabs(Detajj)>1.5)  passMDetajj = false;
      if(fabs(MjjL)>400.)   passMDetajj = false;
      if(fabs(Mjj-80.)>20.) passMDetajj = false;
      bool passMjj = true;
      if(nj30<2)            passMjj = false;
      if(fabs(Mjj-80.)>20.) passMjj = false;
      bool passDetajj = true;
      if(nj30<2)            passDetajj = false;
      if(fabs(Detajj)>1.5)  passDetajj = false;
      if(testevent) cout << "nj " << nj30 << "(" << nj << ") nb " << nb << " Mjj " << Mjj << " MjjL " << MjjL << " Detajj " << Detajj << endl;

      vector<int> vSS, v3l, v, iSS, i3l;
      vector<int> vaSS, va3l, va, iaSS, ia3l;
      vector<int> vlSS, vl3l, vl, ilSS, il3l;
      if(testevent) cout << "nlep_VVV_cutbased_veto " << nlep_VVV_cutbased_veto() << " nlep_VVV_cutbased_veto_noiso " << nlep_VVV_cutbased_veto_noiso() << " nlep_VVV_cutbased_fo " << nlep_VVV_cutbased_fo() << " nlep_VVV_cutbased_tight " << nlep_VVV_cutbased_tight() << " nlep_VVV_baseline " << nlep_VVV_baseline() << endl;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool isSS = false; bool is3l = false;
	bool isaSS = false; bool isa3l = false;
	bool islSS = false; bool isl3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015){
	  if(abs(lep_pdgId()[i])==11&&lep_lostHits()[i]==0){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06){
	    if(lep_p4()[i ].Pt()>20) { i3l.push_back(i); is3l = true; }
	    if(lep_p4()[i ].Pt()>30) { iSS.push_back(i); isSS = true; }
	  }
	}

	if(fabs(lep_p4()[i].Eta())<2.4){
	  if(lep_pass_VVV_cutbased_fo_noiso()[i]&&abs(lep_pdgId()[i])==11&&lep_lostHits()[i]==0/*&&fabs(lep_ip3d()[i])<0.015*/){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.2){
	      if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	      if(((1.+coneCorrPt(i))*lep_p4()[i ].Pt())>20&&!is3l) { il3l.push_back(i); isl3l = true; }
	      if(((1.+coneCorrPt(i))*lep_p4()[i ].Pt())>30&&lep_tightCharge()[i]==2&&!isSS) { ilSS.push_back(i); islSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.2){
	      if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	      if(((1.+coneCorrPt(i))*lep_p4()[i ].Pt())>20&&!is3l) { il3l.push_back(i); isl3l = true; }
	      if(((1.+coneCorrPt(i))*lep_p4()[i ].Pt())>30&&lep_tightCharge()[i]==2&&!isSS) { ilSS.push_back(i); islSS = true; }
	    }
	  }
	  if(lep_pass_VVV_cutbased_fo_noiso()[i]&&abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.4&&fabs(lep_ip3d()[i])<0.015){
	    if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	    if(lep_p4()[i ].Pt()>30&&!isSS) { iaSS.push_back(i); isaSS = true; }
	  }
	  if(lep_pass_VVV_cutbased_veto_noiso()[i]&&abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.4){
	    if(((1.+coneCorrPt(i))*lep_p4()[i ].Pt())>20&&!is3l) { il3l.push_back(i); isl3l = true; }
	    if(((1.+coneCorrPt(i))*lep_p4()[i ].Pt())>30&&!isSS) { ilSS.push_back(i); islSS = true; }
	  }
	}
	bool bla = false;
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_p4()[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4) {
	  //if(lep_pass_VVV_cutbased_veto()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_p4()[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4) {
	  v.push_back(i);
	  if(!isSS) vSS.push_back(i);
	  if(!is3l) {
	    v3l.push_back(i);
	  }
	  va.push_back(i);
	  if(!isSS&&!isaSS) { vaSS.push_back(i); bla=true;}
	  if(!is3l&&!isa3l) va3l.push_back(i);
	  if(!isSS&&!islSS) vlSS.push_back(i);
	  if(!is3l&&!isl3l) vl3l.push_back(i);
	}
	if(testevent)
	  cout << "lep " << i << " pdgid " << lep_pdgId()[i] << " pt " << lep_p4()[i ].Pt() << " eta " << fabs(lep_p4()[i].Eta()) << "(" << fabs(lep_etaSC()[i]) << ") iso " << lep_relIso03EA()[i] << "("<<lep_relIso03EAv2()[i] <<") IP3D " << fabs(lep_ip3d()[i]) << " losthits " << lep_lostHits()[i] << " tightcharge " << lep_tightCharge()[i] << " ID " << lep_pass_VVV_cutbased_tight_noiso()[i] << "/" << lep_pass_VVV_cutbased_fo_noiso()[i] << "/" << lep_pass_VVV_cutbased_fo_noiso()[i] << " lep_pass_VVV_cutbased_veto " << lep_pass_VVV_cutbased_veto()[i] << " conecor " << coneCorrPt(i) << " iSS " << isSS << " iaSS " << isaSS << " isveto " << bla << endl;

      }

      int nvetoSS = vSS.size();
      int nveto3l = v3l.size();
      int nvetoaSS = vaSS.size();
      int nvetoa3l = va3l.size();
      int nveto = v.size();
      int nSS = iSS.size();
      int n3l = i3l.size();
      int naSS = iaSS.size();
      int na3l = ia3l.size();

      int nlSS = ilSS.size();
      int nl3l = il3l.size();
      int nvetolSS = vlSS.size();
      int nvetol3l = vl3l.size();
      if(testevent) cout << "nSS " << nSS << " naSS " << naSS << " n3l " << n3l << " na3l " << na3l << " nvetoaSS " << nvetoaSS << endl;

	    
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
      for(int i = 0; i<na3l; ++i){
	if(abs(lep_pdgId()[ia3l[i] ])==11){
	  ++nel;
	  if(lep_p4()[ia3l[i] ].Pt()>25) ++nel25;
	} else if(abs(lep_pdgId()[ia3l[i] ])==13){
	  ++nmu;
	  if(lep_p4()[ia3l[i] ].Pt()>25) ++nmu25;
	}
      }
   
      if(nmu>=2)           passofflineforTrigger = true;
      if(nmu25>=1&&nel>=1) passofflineforTrigger = true;
      if(nel25>=1&&nel>=2) passofflineforTrigger = true;
      if((nSS+naSS)>=2)    passofflineforTrigger = true;

      if(n3l<2&&(nSS+naSS)<2) continue;
      //if(n3l<2&&(nSS+nlSS)<2) continue;
      //if(nj30<2) continue;
      //if(nb!=0) continue;
      if(!passofflineforTrigger) continue;
      if((nSS+naSS)==0&&lep_p4()[i3l[0] ].Pt()<25) continue;
      //if((nSS+nlSS)==0&&lep_p4()[i3l[0] ].Pt()<25) continue;

      if(testevent) cout << "pass up to trigger" << endl;
      if(isData()){
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) { /*cout << "Event " << tas::run() << ":" << tas::lumi() << ":" << tas::run() << " is duplicated." << endl;*/ continue; }
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = false;
	if(nmu>=2&&(HLT_DoubleMu()) )                             passonlineTrigger = true;
	if(nmu25>=1&&nel>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) passonlineTrigger = true;
	//if(nmu>=2&&(HLT_DoubleMu_noiso()) )                             passonlineTrigger = true;
	//if(nmu25>=1&&nel>=1&&(HLT_MuEG_noiso()||HLT_MuEG_noiso_2()))    passonlineTrigger = true;
	//if(nel25>=1&&nel>=2&&(HLT_DoubleEl_noiso()))                    passonlineTrigger = true;
	if(!passonlineTrigger) continue;
      }
      if(testevent) cout << "pass beyond trigger" << endl;

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
      else if(sn.find("Background")!=string::npos){
	bool isHtoWW = false;
	bool isWnotFromH = false;
	for(unsigned int i = 0; i<genPart_pdgId().size();++i){
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])==25&&genPart_status()[i]==22) { isHtoWW = true; }
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])!=25&&abs(genPart_grandmaId()[i])!=25&&genPart_status()[i]==22) { isWnotFromH = true; }
	  if(isHtoWW&&isWnotFromH) break;
	}
	if(isHtoWW&&isWnotFromH) continue;
	if(((iSS.size()+iaSS.size())>=2)||((iSS.size()+ilSS.size())>=2)){
	  int l1(-1),l2(-1);
	  if(iSS.size()>=2) { l1 = iSS[0]; l2 = iSS[1]; }
	  else if(iSS.size()==1&&iaSS.size()>=1) { l1 = iSS[0]; l2 = iaSS[0]; }
	  else if(iSS.size()==1&&ilSS.size()>=1) { l1 = iSS[0]; l2 = ilSS[0]; }
	  else if(iaSS.size()>=2) { l1 = iaSS[0]; l2 = iaSS[1]; }
	  else if(ilSS.size()>=2) { l1 = ilSS[0]; l2 = ilSS[1]; }
	  //cout << "Put into gentype SS " <<  l1 << " " << l2 << endl;
	  int gentype = gentype_v2(l1,l2,-1);
	  //cout << "Coming out of gentype SS " <<  l1 << " " << l2 << endl;
	  if(     gentype==0) sn = "trueSS";
	  else if(gentype==2) sn = "chargeflips";
	  else if(gentype==3) sn = "SSLL";
	  else if(gentype==4) sn = "fakes";
	  else if(gentype==5) sn = "photonfakes";
	  else                sn = "others";
	}
	if((((i3l.size()+ia3l.size())>=3)||((i3l.size()+il3l.size())>=3))&&i3l.size()>=2){
	  int l1(-1), l2(-1), l3(-1);
	  if(i3l.size()>=3) { l1 = i3l[0]; l2 = i3l[1]; l3 = i3l[2]; }
	  if(i3l.size()==2&&ia3l.size()>=1) { l1 = i3l[0]; l2 = i3l[1]; l3 = ia3l[0]; }
	  if(i3l.size()==2&&il3l.size()>=1) { l1 = i3l[0]; l2 = i3l[1]; l3 = il3l[0]; }
	  //cout << "Put into gentype 3l " <<  l1 << " " << l2 << endl;
	  int gentype = gentype_v2(l1,l2,l3);
	  //cout << "Coming out of gentype 3l " <<  l1 << " " << l2 << endl;
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
      float aMTmax = -1;
      if(iSS.size()==1&&iaSS.size()>=1){
	//cout << __LINE__<<endl;
	if(mT(lep_p4()[iaSS[0] ],MET)>mT(lep_p4()[iSS[0] ],MET)) aMTmax = mT(lep_p4()[iaSS[0] ],MET);
	else aMTmax = mT(lep_p4()[iSS[0] ],MET);
	//cout << aMTmax << endl;
      }
      else if(iSS.size()==0&&iaSS.size()>=2){
	//cout << __LINE__<<endl;
	if(mT(lep_p4()[iaSS[0] ],MET)>mT(lep_p4()[iaSS[1] ],MET)) aMTmax = mT(lep_p4()[iaSS[0] ],MET);
	else aMTmax = mT(lep_p4()[iaSS[1] ],MET);
	//cout << aMTmax << endl;
      }
      else if(iSS.size()==0&&ilSS.size()>=2){
	//cout << __LINE__<<endl;
	if(mT(lep_p4()[ilSS[0] ],MET)>mT(lep_p4()[ilSS[1] ],MET)) aMTmax = mT(lep_p4()[ilSS[0] ],MET);
	else aMTmax = mT(lep_p4()[ilSS[1] ],MET);
	//cout << aMTmax << endl;
      }
      bool ee[50],mm[50],em[50];
      for(int i = 0; i<50; ++i) { ee[i] = false; em[i] = false; mm[i] = false; }
      //0: SR
      //cout << __LINE__ << endl;
      if(testevent&&nSS==1&&naSS==1){
	cout << "isotrack " << nisoTrack_mt2_cleaned_VVV_cutbased_veto() << " Mll " << (lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M() << " met " << met_pt() << " MTmax " << MTmax << endl;
      }
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	//cout << __LINE__ << endl;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                        mm[0] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[0] = false; mm[0] = false; }
	//cout << __LINE__ << endl;
      }
      //1: SRpreselect
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) ee[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11) em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13) em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) mm[1] = true;
      }
      //cout << __LINE__ << endl;
      //2: AR
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.){
	//cout << __LINE__ << endl;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) ee[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&aMTmax>90.)                                               em[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13&&met_pt()>40.&&aMTmax>90.)                                               em[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13)                                                                         mm[2] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()<40){ ee[2] = false; mm[2] = false; }
	//cout << __LINE__ << endl;
      }
      //3: ARpreselect
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11) ee[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11) em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13) em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13) mm[3] = true;
      }
      //4: AR-cc
      //cout << __LINE__ << endl;
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetolSS==0&&nSS==1&&nlSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[ilSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[ilSS[0] ]).M()>30.){
	//cout << __LINE__ << endl;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[ilSS[0] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[ilSS[0] ]).M()-90.)>10.) ee[4] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[ilSS[0] ])==11&&met_pt()>40.&&aMTmax>90.)                                               em[4] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[ilSS[0] ])==13&&met_pt()>40.&&aMTmax>90.)                                               em[4] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[ilSS[0] ])==13)                                                                         mm[4] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[ilSS[0] ]).M()<40){ ee[4] = false; mm[4] = false; }
	//cout << __LINE__ << endl;
      }
      //5: ARpreselect-cc
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetolSS==0&&nSS==1&&nlSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[ilSS[0] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[ilSS[0] ])==11) ee[5] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[ilSS[0] ])==11) em[5] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[ilSS[0] ])==13) em[5] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[ilSS[0] ])==13) mm[5] = true;
      }
      //6: SR: Mjjsideband, no MET
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	//cout << __LINE__ << endl;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[6] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&MTmax>90.)                                               em[6] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&MTmax>90.)                                               em[6] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                          mm[6] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[6] = false; mm[6] = false; }
	if(fabs(Mjj-80.)<20.){ ee[6] = false; em[6] = false; mm[6] = false; }
	//cout << __LINE__ << endl;
      }
      //7: SRpresel: Mjjsideband, no MET
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[7] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11)                                                          em[7] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13)                                                          em[7] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                          mm[7] = true;
	if(fabs(Mjj-80.)<20.){ ee[7] = false; em[7] = false; mm[7] = false; }
      }
      //8: AR: Mjjsideband, no MET
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) ee[8] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11&&aMTmax>90.)                                               em[8] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13&&aMTmax>90.)                                               em[8] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13)                                                           mm[8] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()<40){ ee[8] = false; mm[8] = false; }
	if(fabs(Mjj-80.)<20.){ ee[8] = false; em[8] = false; mm[8] = false; }
      }
      //9: ARpresel: Mjjsideband, no MET
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) ee[9] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11)                                                            em[9] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13)                                                            em[9] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13)                                                            mm[9] = true;
	if(fabs(Mjj-80.)<20.){ ee[9] = false; em[9] = false; mm[9] = false; }
      }
      for(int i = 0; i<50; ++i) {
	//cout << __LINE__ << endl;
	if((ee[i]||em[i]||mm[i]) && nSS==2){
	  if(fname.find("wjets")!=string::npos&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")!=string::npos&&lep_motherIdSS()[iSS[0] ]==(-3))   { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")!=string::npos&&lep_motherIdSS()[iSS[1] ]==(-3))   { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if((ee[i]||em[i]||mm[i]) && nSS==1&&naSS==1){
	  if(fname.find("wjets")!=string::npos&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&lep_motherIdSS()[iaSS[0] ]==(-3)){ ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")!=string::npos&&lep_motherIdSS()[iSS[0] ]==(-3))   { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")!=string::npos&&lep_motherIdSS()[iaSS[0] ]==(-3))  { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if((ee[i]||em[i]||mm[i]) && nSS==1&&nlSS==1){
	  if(fname.find("wjets")!=string::npos&&lep_motherIdSS()[iaSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&lep_motherIdSS()[iaSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")!=string::npos&&lep_motherIdSS()[iaSS[0] ]==(-3))   { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")!=string::npos&&lep_motherIdSS()[iaSS[1] ]==(-3))   { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	//cout << __LINE__ << endl;
      }
      int SFOS[50];
      for(int i = 0; i<50; ++i) { SFOS[i] = -1; }
      double pTlll(-1), DPhilllMET(-1); double Mmumu(-1), Mmumue(-1);
      //cout << __LINE__ << endl;
      if(i3l.size()>=3){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
      } else if(i3l.size()==2&&ia3l.size()>=1){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[ia3l[0]  ]).Pt();
      } else if(i3l.size()==2&&il3l.size()>=1){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[il3l[0]  ]).Pt();
      }
      //cout << __LINE__ << endl;
      bool passlowMSFOS = true;
      bool passlowMlll = true;
      //1,2: SR
      if(nj<2&&nb==0&&(nveto3l==0&&n3l==3)){
	//cout << __LINE__ << endl;
	int SFOScounter = 0;
	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	bool pass0X(false), pass1X(false), pass2X(false);
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) SFOScounter = -1;
	//upper three lines require an mu+mu- pair
	  
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) { pass0 = false; pass0X = false; }
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) { pass0 = false; pass0X = false; }
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) { pass0 = false; pass0X = false; }
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass0X = false;
	  else cout<<OS01<<" "<<SF01<<" ("<<(OS01&&SF01)<<") "<<OS02<<" "<<SF02<<" ("<<(OS02&&SF02)<<") "<<OS12<<" "<<SF12<<" ("<<(OS12&&SF12)<<") "<<endl;	  
	}
	if(SFOScounter==1){
	  pass1 = true;
	  pass1X = true;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passlowMSFOS = false;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) passlowMSFOS = false;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) passlowMSFOS = false;
	  if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<10.) passlowMlll = false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passlowMSFOS = false;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) passlowMSFOS = false;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) passlowMSFOS = false;
	  if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<10.) passlowMlll = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	if(pass0){
	  if((DPhilllMET>2.7&&pTlll>60.)&&passlowMSFOS)               SFOS[0] = 0;
	}
	if(pass1){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)&&passlowMSFOS&&passlowMlll) SFOS[0] = 1;
	}
	if(pass2){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)&&passlowMSFOS&&passlowMlll) SFOS[0] = 2;
	}
	if(pass0X){
	  if((DPhilllMET>2.7&&pTlll>60.))               SFOS[1] = 0;
	}
	if(pass1X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[1] = 1;
	}
	if(pass2X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[1] = 2;
	}
	//cout << __LINE__ << endl;
      }
      //3-4:AR
      if(nj<2&&nb==0&&(nvetoa3l==0&&n3l==2&&na3l==1)){
	//cout << __LINE__ << endl;
	int SFOScounter = 0;
	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[ia3l[0] ]<0);
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[ia3l[0] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[ia3l[0] ]));
	bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[ia3l[0] ]));
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	bool pass0X(false), pass1X(false), pass2X(false);
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[ia3l[0] ])==3) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[ia3l[0] ])==11&&lep_motherIdSS()[ia3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[ia3l[0] ])==11&&lep_motherIdSS()[ia3l[0] ]==(-3)) SFOScounter = -1;
	//upper three lines require an mu+mu- pair
	  
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) { pass0 = false; pass0X = false; }
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[ia3l[0] ]).M()<20.) { pass0 = false; pass0X = false; }
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()<20.) { pass0 = false; pass0X = false; }
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[ia3l[0] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass0X = false;
	  else cout<<OS01<<" "<<SF01<<" ("<<(OS01&&SF01)<<") "<<OS02<<" "<<SF02<<" ("<<(OS02&&SF02)<<") "<<OS12<<" "<<SF12<<" ("<<(OS12&&SF12)<<") "<<endl;	  
	}
	if(SFOScounter==1){
	  pass1 = true;
	  pass1X = true;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passlowMSFOS = false;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[ia3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()-90.)<10.) passlowMlll = false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[ia3l[0] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[ia3l[0] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()>55.&&(lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passlowMSFOS = false;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[ia3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()-90.)<10.) passlowMlll = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[ia3l[0] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[ia3l[0] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	if(pass0){
	  if((DPhilllMET>2.7&&pTlll>60.)&&passlowMSFOS)               SFOS[2] = 0;
	}
	if(pass1){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)&&passlowMSFOS&&passlowMlll) SFOS[2] = 1;
	}
	if(pass2){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)&&passlowMSFOS&&passlowMlll) SFOS[2] = 2;
	}
	if(pass0X){
	  if((DPhilllMET>2.7&&pTlll>60.))               SFOS[3] = 0;
	}
	if(pass1X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[3] = 1;
	}
	if(pass2X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[3] = 2;
	}
	//cout << __LINE__ << endl;
      }
      //5-6:AR
      if(nj<2&&nb==0&&(nvetol3l==0&&n3l==2&&nl3l==1)){
	//cout << __LINE__ << endl;
	int SFOScounter = 0;
	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[il3l[0] ]<0);
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[il3l[0] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[il3l[0] ]));
	bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[il3l[0] ]));
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	bool pass0X(false), pass1X(false), pass2X(false);
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[il3l[0] ])==3) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[il3l[0] ])==11&&lep_motherIdSS()[il3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&abs(lep_pdgId()[il3l[0] ])==11&&lep_motherIdSS()[il3l[0] ]==(-3)) SFOScounter = -1;
	//upper three lines require an mu+mu- pair
	  
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) { pass0 = false; pass0X = false; }
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[il3l[0] ]).M()<20.) { pass0 = false; pass0X = false; }
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()<20.) { pass0 = false; pass0X = false; }
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[il3l[0] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass0X = false;
	  else cout<<OS01<<" "<<SF01<<" ("<<(OS01&&SF01)<<") "<<OS02<<" "<<SF02<<" ("<<(OS02&&SF02)<<") "<<OS12<<" "<<SF12<<" ("<<(OS12&&SF12)<<") "<<endl;	  
	}
	if(SFOScounter==1){
	  pass1 = true;
	  pass1X = true;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passlowMSFOS = false;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[il3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()-90.)<10.) passlowMlll = false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[il3l[0] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[il3l[0] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()>55.&&(lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passlowMSFOS = false;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[il3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()<20.) passlowMSFOS = false;
	  if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()-90.)<10.) passlowMlll = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[il3l[0] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[il3l[0] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	if(pass0){
	  if((DPhilllMET>2.7&&pTlll>60.)&&passlowMSFOS)               SFOS[4] = 0;
	}
	if(pass1){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)&&passlowMSFOS&&passlowMlll) SFOS[4] = 1;
	}
	if(pass2){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)&&passlowMSFOS&&passlowMlll) SFOS[4] = 2;
	}
	if(pass0X){
	  if((DPhilllMET>2.7&&pTlll>60.))               SFOS[5] = 0;
	}
	if(pass1X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[5] = 1;
	}
	if(pass2X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[5] = 2;
	}
	//cout << __LINE__ << endl;
      }
      /*
        if( lep_pdgId().size() ==2 && (nSS+naSS)==2&&naSS==1&&nvetoaSS==0){
	cout << sn << " " << endl;
	cout << ee[0] << " " << em[0] << " " << mm[0] << "  " << ee[2] << " " << em[2] << " " << mm[2] << "  " << ee[4] << " " << em[4] << " " << mm[4] << endl;
	cout << ee[1] << " " << em[1] << " " << mm[1] << "  " << ee[3] << " " << em[3] << " " << mm[3] << "  " << ee[5] << " " << em[5] << " " << mm[5] << endl;
	cout << lep_pdgId().size() << " " << nSS << " " << n3l << " a " << naSS << " " << na3l << " v " << nvetoSS << " " << nveto3l << " va " << nvetoaSS << " " << nvetoa3l << endl;
	for(unsigned int i = 0; i<lep_pdgId().size();++i){
	  cout << i << " " << lep_pdgId()[i] << " " << lep_p4()[i].Pt() << " " << fabs(lep_p4()[i].Eta()) << " " << lep_relIso03EAv2()[i] << " " << fabs(lep_ip3d()[i]) << " " << lep_pass_VVV_cutbased_tight_noiso()[i] << " " << lep_pass_VVV_cutbased_fo_noiso()[i] <<  " " << lep_tightCharge()[i] << " " << lep_lostHits()[i] << endl;
	}
      }
      */
      //cout << __LINE__ << " " << sn << " " << sn2 << endl;
      if(ee[1]     ) histos["PreselSRyield_"+sn ]->Fill(0.,weight);
      if(em[1]     ) histos["PreselSRyield_"+sn ]->Fill(1.,weight);
      if(mm[1]     ) histos["PreselSRyield_"+sn ]->Fill(2.,weight);
      if(ee[0]     ) histos["SRyield_"+sn ]->Fill(0.,weight);
      if(em[0]     ) histos["SRyield_"+sn ]->Fill(1.,weight);
      if(mm[0]     ) histos["SRyield_"+sn ]->Fill(2.,weight);
      if(SFOS[0]==0) histos["SRyield_"+sn2]->Fill(3.,weight);
      if(SFOS[0]==1) histos["SRyield_"+sn2]->Fill(4.,weight);
      if(SFOS[0]==2) histos["SRyield_"+sn2]->Fill(5.,weight);
       //cout << __LINE__ << endl;
      if(ee[6]&&met_pt()>40.     ) histos["SRyield_Mjjsideband_"+sn ]->Fill(0.,weight);
      if(em[6]&&met_pt()>40.     ) histos["SRyield_Mjjsideband_"+sn ]->Fill(1.,weight);
      if(mm[6]                   ) histos["SRyield_Mjjsideband_"+sn ]->Fill(2.,weight);
      if(ee[6]&&met_pt()<40.     ) histos["SRyield_MjjsidebandlowMET_"+sn ]->Fill(0.,weight);
      if(em[6]&&met_pt()<40.     ) histos["SRyield_MjjsidebandlowMET_"+sn ]->Fill(1.,weight);
      if(mm[6]&&met_pt()<40.     ) histos["SRyield_MjjsidebandlowMET_"+sn ]->Fill(2.,weight);
      if(ee[7]                   ) histos["PreselSRyield_Mjjsideband_"+sn ]->Fill(0.,weight);
      if(em[7]                   ) histos["PreselSRyield_Mjjsideband_"+sn ]->Fill(1.,weight);
      if(mm[7]                   ) histos["PreselSRyield_Mjjsideband_"+sn ]->Fill(2.,weight);
      if(ee[7]&&met_pt()<40.     ) histos["PreselSRyield_MjjsidebandlowMET_"+sn ]->Fill(0.,weight);
      if(em[7]&&met_pt()<40.     ) histos["PreselSRyield_MjjsidebandlowMET_"+sn ]->Fill(1.,weight);
      if(mm[7]&&met_pt()<40.     ) histos["PreselSRyield_MjjsidebandlowMET_"+sn ]->Fill(2.,weight);

      
      if(ee[3]||em[3]||mm[3]){
      //cout << __LINE__ << endl;
	double myFRT = -1;
	double myFRTerr = -1;
	if(abs(lep_pdgId()[iaSS[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	double myFR = myFRT/(1.-myFRT);
	double myFRerr = myFRTerr/pow(1.-myFRT,2);
	//cout << __LINE__ << endl;
	if(ee[3]){
	  histos["PreselARyield_"+sn]->Fill(0.,weight);
	  histos["PreselFakeEstimation_"+sn]->Fill(0.,weight*myFR);
	  histos["PreselFakeEstimationFRup_"+sn]->Fill(0.,weight*(myFR+myFRerr));
	  histos["PreselFakeEstimationFRdn_"+sn]->Fill(0.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(em[3]){
	  histos["PreselARyield_"+sn]->Fill(1.,weight);
	  histos["PreselFakeEstimation_"+sn]->Fill(1.,weight*myFR);
	  histos["PreselFakeEstimationFRup_"+sn]->Fill(1.,weight*(myFR+myFRerr));
	  histos["PreselFakeEstimationFRdn_"+sn]->Fill(1.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(mm[3]){
	  histos["PreselARyield_"+sn]->Fill(2.,weight);
	  histos["PreselFakeEstimation_"+sn]->Fill(2.,weight*myFR);
	  histos["PreselFakeEstimationFRup_"+sn]->Fill(2.,weight*(myFR+myFRerr));
	  histos["PreselFakeEstimationFRdn_"+sn]->Fill(2.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
      }
      if(isData()){
	if(ee[2])      *eventstocheckEE    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	if(em[2])      *eventstocheckEM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	if(mm[2])      *eventstocheckMM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	if(SFOS[2]==0) *eventstocheck0SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	if(SFOS[2]==1) *eventstocheck1SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	if(SFOS[2]==2) *eventstocheck2SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      }
     if(ee[2]||em[2]||mm[2]){
      //cout << __LINE__ << endl;
	double myFRT = -1;
	double myFRTerr = -1;
	if(abs(lep_pdgId()[iaSS[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	double myFR = myFRT/(1.-myFRT);
	double myFRerr = myFRTerr/pow(1.-myFRT,2);
	//cout << __LINE__ << endl;
	if(ee[2]){
	  histos["ARyield_"+sn]->Fill(0.,weight);
	  histos["FakeEstimation_"+sn]->Fill(0.,weight*myFR);
	  histos["FakeEstimationFRup_"+sn]->Fill(0.,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_"+sn]->Fill(0.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(em[2]){
	  histos["ARyield_"+sn]->Fill(1.,weight);
	  histos["FakeEstimation_"+sn]->Fill(1.,weight*myFR);
	  histos["FakeEstimationFRup_"+sn]->Fill(1.,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_"+sn]->Fill(1.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(mm[2]){
	  histos["ARyield_"+sn]->Fill(2.,weight);
	  histos["FakeEstimation_"+sn]->Fill(2.,weight*myFR);
	  histos["FakeEstimationFRup_"+sn]->Fill(2.,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_"+sn]->Fill(2.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
	if(abs(lep_pdgId()[iaSS[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.)*lep_p4()[iaSS[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.)*lep_p4()[iaSS[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	myFR = myFRT/(1.-myFRT);
	myFRerr = myFRTerr/pow(1.-myFRT,2);
	if(ee[2]){
	  histos["NoConeCorrARyield_"+sn]->Fill(0.,weight);
	  histos["NoConeCorrFakeEstimation_"+sn]->Fill(0.,weight*myFR);
	  histos["NoConeCorrFakeEstimationFRup_"+sn]->Fill(0.,weight*(myFR+myFRerr));
	  histos["NoConeCorrFakeEstimationFRdn_"+sn]->Fill(0.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(em[2]){
	  histos["NoConeCorrARyield_"+sn]->Fill(1.,weight);
	  histos["NoConeCorrFakeEstimation_"+sn]->Fill(1.,weight*myFR);
	  histos["NoConeCorrFakeEstimationFRup_"+sn]->Fill(1.,weight*(myFR+myFRerr));
	  histos["NoConeCorrFakeEstimationFRdn_"+sn]->Fill(1.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(mm[2]){
	  histos["NoConeCorrARyield_"+sn]->Fill(2.,weight);
	  histos["NoConeCorrFakeEstimation_"+sn]->Fill(2.,weight*myFR);
	  histos["NoConeCorrFakeEstimationFRup_"+sn]->Fill(2.,weight*(myFR+myFRerr));
	  histos["NoConeCorrFakeEstimationFRdn_"+sn]->Fill(2.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
      }
      if(SFOS[2]>=0){
	//cout << __LINE__ <<  " " << ia3l.size() << endl;
	double myFRT = -1;
	double myFRTerr = -1;
	if(abs(lep_pdgId()[ia3l[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(ia3l[0]))*lep_p4()[ia3l[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[ia3l[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(ia3l[0]))*lep_p4()[ia3l[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[ia3l[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	double myFR = myFRT/(1.-myFRT);
	double myFRerr = myFRTerr/pow(1.-myFRT,2);
	//cout << __LINE__ << endl;
	if(SFOS[2]==0){
	  histos["ARyield_"+sn2]->Fill(3.,weight);
	  histos["FakeEstimation_"+sn2]->Fill(3.,weight*myFR);
	  histos["FakeEstimationFRup_"+sn2]->Fill(3.,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_"+sn2]->Fill(3.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
	if(SFOS[2]==1){
	  histos["ARyield_"+sn2]->Fill(4.,weight);
	  histos["FakeEstimation_"+sn2]->Fill(4.,weight*myFR);
	  histos["FakeEstimationFRup_"+sn2]->Fill(4.,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_"+sn2]->Fill(4.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
	if(SFOS[2]==2){
	  histos["ARyield_"+sn2]->Fill(5.,weight);
	  histos["FakeEstimation_"+sn2]->Fill(5.,weight*myFR);
	  histos["FakeEstimationFRup_"+sn2]->Fill(5.,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_"+sn2]->Fill(5.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
	if(abs(lep_pdgId()[ia3l[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.)*lep_p4()[ia3l[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[ia3l[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.)*lep_p4()[ia3l[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[ia3l[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	myFR = myFRT/(1.-myFRT);
	myFRerr = myFRTerr/pow(1.-myFRT,2);
	if(SFOS[2]==0){
	  histos["NoConeCorrARyield_"+sn2]->Fill(3.,weight);
	  histos["NoConeCorrFakeEstimation_"+sn2]->Fill(3.,weight*myFR);
	  histos["NoConeCorrFakeEstimationFRup_"+sn2]->Fill(3.,weight*(myFR+myFRerr));
	  histos["NoConeCorrFakeEstimationFRdn_"+sn2]->Fill(3.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(SFOS[2]==1){
	  histos["NoConeCorrARyield_"+sn2]->Fill(4.,weight);
	  histos["NoConeCorrFakeEstimation_"+sn2]->Fill(4.,weight*myFR);
	  histos["NoConeCorrFakeEstimationFRup_"+sn2]->Fill(4.,weight*(myFR+myFRerr));
	  histos["NoConeCorrFakeEstimationFRdn_"+sn2]->Fill(4.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(SFOS[2]==2){
	  histos["NoConeCorrARyield_"+sn2]->Fill(5.,weight);
	  histos["NoConeCorrFakeEstimation_"+sn2]->Fill(5.,weight*myFR);
	  histos["NoConeCorrFakeEstimationFRup_"+sn2]->Fill(5.,weight*(myFR+myFRerr));
	  histos["NoConeCorrFakeEstimationFRdn_"+sn2]->Fill(5.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
      }
      if(ee[4]||em[4]||mm[4]){
	//cout << __LINE__ << endl;
	double myFRT = -1;
	double myFRTerr = -1;
	if(abs(lep_pdgId()[ilSS[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(ilSS[0]))*lep_p4()[ilSS[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[ilSS[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(ilSS[0]))*lep_p4()[ilSS[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[ilSS[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	double myFR = myFRT/(1.-myFRT);
	double myFRerr = myFRTerr/pow(1.-myFRT,2);
	if(ee[4]){
	  histos["AllConeCorrARyield_"+sn]->Fill(0.,weight);
	  histos["AllConeCorrFakeEstimation_"+sn]->Fill(0.,weight*myFR);
	  histos["AllConeCorrFakeEstimationFRup_"+sn]->Fill(0.,weight*(myFR+myFRerr));
	  histos["AllConeCorrFakeEstimationFRdn_"+sn]->Fill(0.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(em[4]){
	  histos["AllConeCorrARyield_"+sn]->Fill(1.,weight);
	  histos["AllConeCorrFakeEstimation_"+sn]->Fill(1.,weight*myFR);
	  histos["AllConeCorrFakeEstimationFRup_"+sn]->Fill(1.,weight*(myFR+myFRerr));
	  histos["AllConeCorrFakeEstimationFRdn_"+sn]->Fill(1.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(mm[4]){
	  histos["AllConeCorrARyield_"+sn]->Fill(2.,weight);
	  histos["AllConeCorrFakeEstimation_"+sn]->Fill(2.,weight*myFR);
	  histos["AllConeCorrFakeEstimationFRup_"+sn]->Fill(2.,weight*(myFR+myFRerr));
	  histos["AllConeCorrFakeEstimationFRdn_"+sn]->Fill(2.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
      }
      if(SFOS[4]>=0){
	//cout << __LINE__ << endl;
	double myFRT = -1;
	double myFRTerr = -1;
	if(abs(lep_pdgId()[il3l[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(il3l[0]))*lep_p4()[il3l[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[il3l[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(il3l[0]))*lep_p4()[il3l[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[il3l[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	double myFR = myFRT/(1.-myFRT);
	double myFRerr = myFRTerr/pow(1.-myFRT,2);
	if(SFOS[4]==0){
	  histos["AllConeCorrARyield_"+sn2]->Fill(3.,weight);
	  histos["AllConeCorrFakeEstimation_"+sn2]->Fill(3.,weight*myFR);
	  histos["AllConeCorrFakeEstimationFRup_"+sn2]->Fill(3.,weight*(myFR+myFRerr));
	  histos["AllConeCorrFakeEstimationFRdn_"+sn2]->Fill(3.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(SFOS[4]==1){
	  histos["AllConeCorrARyield_"+sn2]->Fill(4.,weight);
	  histos["AllConeCorrFakeEstimation_"+sn2]->Fill(4.,weight*myFR);
	  histos["AllConeCorrFakeEstimationFRup_"+sn2]->Fill(4.,weight*(myFR+myFRerr));
	  histos["AllConeCorrFakeEstimationFRdn_"+sn2]->Fill(4.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if(SFOS[4]==2){
	  histos["AllConeCorrARyield_"+sn2]->Fill(5.,weight);
	  histos["AllConeCorrFakeEstimation_"+sn2]->Fill(5.,weight*myFR);
	  histos["AllConeCorrFakeEstimationFRup_"+sn2]->Fill(5.,weight*(myFR+myFRerr));
	  histos["AllConeCorrFakeEstimationFRdn_"+sn2]->Fill(5.,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	//cout << __LINE__ << endl;
      }

      if(ee[9]||em[9]||mm[9]){
      //cout << __LINE__ << endl;
	float entry = -1;
	if(ee[9]) entry = 0.;
	if(em[9]) entry = 1.;
	if(mm[9]) entry = 2.;
	double myFRT = -1;
	double myFRTerr = -1;
	if(abs(lep_pdgId()[iaSS[0] ])==11) {
	  int bin = hElFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	  myFRT = hElFR->GetBinContent(bin);
	  myFRTerr = hElFR->GetBinError(bin);
	} else {
	  int bin = hMuFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				   std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	  myFRT = hMuFR->GetBinContent(bin);
	  myFRTerr = hMuFR->GetBinError(bin);
	}
	double myFR = myFRT/(1.-myFRT);
	double myFRerr = myFRTerr/pow(1.-myFRT,2);
	if(ee[9]||em[9]||mm[9]){
	  histos["PreselARyield_Mjjsideband_"+sn]->Fill(entry,weight);
	  histos["PreselFakeEstimation_Mjjsideband_"+sn]->Fill(entry,weight*myFR);
	  histos["PreselFakeEstimationFRup_Mjjsideband_"+sn]->Fill(entry,weight*(myFR+myFRerr));
	  histos["PreselFakeEstimationFRdn_Mjjsideband_"+sn]->Fill(entry,weight*std::max((double)myFR-myFRerr,(double)0.));
	  if(met_pt()<40.){
	    histos["PreselARyield_MjjsidebandlowMET_"+sn]->Fill(entry,weight);
	    histos["PreselFakeEstimation_MjjsidebandlowMET_"+sn]->Fill(entry,weight*myFR);
	    histos["PreselFakeEstimationFRup_MjjsidebandlowMET_"+sn]->Fill(entry,weight*(myFR+myFRerr));
	    histos["PreselFakeEstimationFRdn_MjjsidebandlowMET_"+sn]->Fill(entry,weight*std::max((double)myFR-myFRerr,(double)0.));
	  }
	}
	if(((ee[8]||em[8])&&met_pt()>40.)||mm[8]){
	  histos["ARyield_Mjjsideband_"+sn]->Fill(entry,weight);
	  histos["FakeEstimation_Mjjsideband_"+sn]->Fill(entry,weight*myFR);
	  histos["FakeEstimationFRup_Mjjsideband_"+sn]->Fill(entry,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_Mjjsideband_"+sn]->Fill(entry,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
	if((ee[8]||em[8]||mm[8])&&met_pt()<40.){
	  histos["ARyield_MjjsidebandlowMET_"+sn]->Fill(entry,weight);
	  histos["FakeEstimation_MjjsidebandlowMET_"+sn]->Fill(entry,weight*myFR);
	  histos["FakeEstimationFRup_MjjsidebandlowMET_"+sn]->Fill(entry,weight*(myFR+myFRerr));
	  histos["FakeEstimationFRdn_MjjsidebandlowMET_"+sn]->Fill(entry,weight*std::max((double)myFR-myFRerr,(double)0.));
	}
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

  cout << "Are these too much? " << counter << endl;
  cout << "ee " << wgtsum1 << " em " << wgtsum2 << " mm " << wgtsum3 << endl;
  
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //add overflow
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
    //add underflow
    h->second->SetBinContent(1, h->second->GetBinContent(1)+ h->second->GetBinContent(0) );
    h->second->SetBinError(1, sqrt(pow(h->second->GetBinError(1),2)+pow(h->second->GetBinError(0),2) ) );
  }
  
  string filename = "rootfiles/FakeRatePrediction.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //string mapname = histonames[i]+"_"+skimFilePrefix;
    //string writename = histonames[i];
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;

  cout << "ee AR events: " << endl;
  cout    << eventstocheckEE   ->str();
  cout << endl << "emu AR events: " << endl;
  cout    << eventstocheckEM   ->str();
  cout << endl << "mumu AR events: " << endl;
  cout    << eventstocheckMM   ->str();
  cout << endl << "0SFOS AR events: " << endl;
  cout << eventstocheck0SFOS->str();
  cout << endl << "1SFOS AR events: " << endl;
  cout << eventstocheck1SFOS->str();
  cout << endl << "2SFOS AR events: " << endl;
  cout << eventstocheck2SFOS->str();
  cout << endl;
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  fSF->Close();
  delete fSF;
  fFR->Close();
  delete fFR;
  return 0;
}
