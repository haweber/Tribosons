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

float MT2(LorentzVector vec1,LorentzVector vec2, LorentzVector met, bool massive ){

  double pa[3];
  double pb[3];
  double pmiss[3];
  
  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (met.Px());
  pmiss[2] = static_cast<double> (met.Py());
  
  pa[0] = static_cast<double> (massive ? vec1.M() : 0);
  pa[1] = static_cast<double> (vec1.Px());
  pa[2] = static_cast<double> (vec1.Py());
  
  pb[0] = static_cast<double> (massive ? vec2.M() : 0);
  pb[1] = static_cast<double> (vec2.Px());
  pb[2] = static_cast<double> (vec2.Py());
  
  mt2_bisect::mt2 *mymt2 = new mt2_bisect::mt2();
  mymt2->set_momenta(pa, pb, pmiss);
  mymt2->set_mn(0);
  Float_t MT2=mymt2->get_mt2();
  delete mymt2;
  return MT2;
  
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

  bool applylepSF = true;

  int cCRee(0), cCR2ee(0), cCR3ee(0);
  int cCRem(0), cCR2em(0), cCR3em(0);
  int cCRmm(0), cCR2mm(0), cCR3mm(0);
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
  TFile *fSF = new TFile("rootfiles/SF_TnP.root","read");
  TH2F *hMu = (TH2F*)fSF->Get("muSF");
  TH2F *hElReco = (TH2F*)fSF->Get("elSF_reco");
  TH2F *hElID = (TH2F*)fSF->Get("elSF_ID");
  float muptmin = 20.1;                           float muptmax = 199.9; float muetamin =  0.01; float muetamax = 2.49;
  float elptmin = 10.1; float elptminReco = 25.1; float elptmax = 499.9; float eletamin = -2.49; float eletamax = 2.49;
   
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

  histonames.push_back("YieldsSR");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsCR");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsAR");                           hbins.push_back(3); hlow.push_back(0); hup.push_back(3);//application region - 1 loose but not tight lepton
  histonames.push_back("YieldsSRpresel");                     hbins.push_back(3); hlow.push_back(0); hup.push_back(3);//blind the data
  histonames.push_back("YieldsARpresel");                     hbins.push_back(3); hlow.push_back(0); hup.push_back(3);//application region - 1 loose but not tight lepton
  histonames.push_back("SRee_cuts");                          hbins.push_back(13);hlow.push_back(0); hup.push_back(13);
  histonames.push_back("SRem_cuts");                          hbins.push_back(17);hlow.push_back(0); hup.push_back(17);
  histonames.push_back("ARee_cuts");                          hbins.push_back(13);hlow.push_back(0); hup.push_back(13);
  histonames.push_back("ARem_cuts");                          hbins.push_back(17);hlow.push_back(0); hup.push_back(17);
  
  histonames.push_back("elSR_motherIDSS");                    hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);
  histonames.push_back("elAR_motherIDSS");                    hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);
  histonames.push_back("elSRpresel_motherIDSS");              hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);
  histonames.push_back("elARpresel_motherIDSS");              hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);
  histonames.push_back("muSR_motherIDSS");                    hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);
  histonames.push_back("muAR_motherIDSS");                    hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);
  histonames.push_back("muSRpresel_motherIDSS");              hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);
  histonames.push_back("muARpresel_motherIDSS");              hbins.push_back(9); hlow.push_back(-4); hup.push_back(5);

  histonames.push_back("elSR_isFromFlags");                   hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);
  histonames.push_back("elAR_isFromFlags");                   hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);
  histonames.push_back("elSRpresel_isFromFlags");             hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);
  histonames.push_back("elARpresel_isFromFlags");             hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);
  histonames.push_back("muSR_isFromFlags");                   hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);
  histonames.push_back("muAR_isFromFlags");                   hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);
  histonames.push_back("muSRpresel_isFromFlags");             hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);
  histonames.push_back("muARpresel_isFromFlags");             hbins.push_back(8); hlow.push_back(-1); hup.push_back(7);

  histonames.push_back("elSR_myMatch");                       hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);
  histonames.push_back("elAR_myMatch");                       hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);
  histonames.push_back("elSRpresel_myMatch");                 hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);
  histonames.push_back("elARpresel_myMatch");                 hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);
  histonames.push_back("muSR_myMatch");                       hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);
  histonames.push_back("muAR_myMatch");                       hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);
  histonames.push_back("muSRpresel_myMatch");                 hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);
  histonames.push_back("muARpresel_myMatch");                 hbins.push_back(9); hlow.push_back(-1); hup.push_back(8);

  histonames.push_back("elAR_tight_pT");                      hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("elAR_tight_eta");                     hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("elAR_tight_phi");                     hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("elAR_tight_reliso");                  hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);
  histonames.push_back("elAR_loose_pT");                      hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("elAR_loose_eta");                     hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("elAR_loose_phi");                     hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("elAR_loose_reliso");                  hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);
  histonames.push_back("muAR_tight_pT");                      hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("muAR_tight_eta");                     hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("muAR_tight_phi");                     hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("muAR_tight_reliso");                  hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);
  histonames.push_back("muAR_loose_pT");                      hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("muAR_loose_eta");                     hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("muAR_loose_phi");                     hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("muAR_loose_reliso");                  hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);

  histonames.push_back("elARpresel_tight_pT");                hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("elARpresel_tight_eta");               hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("elARpresel_tight_phi");               hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("elARpresel_tight_reliso");            hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);
  histonames.push_back("elARpresel_loose_pT");                hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("elARpresel_loose_eta");               hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("elARpresel_loose_phi");               hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("elARpresel_loose_reliso");            hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);
  histonames.push_back("muARpresel_tight_pT");                hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("muARpresel_tight_eta");               hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("muARpresel_tight_phi");               hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("muARpresel_tight_reliso");            hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);
  histonames.push_back("muARpresel_loose_pT");                hbins.push_back(19); hlow.push_back(20); hup.push_back(200);
  histonames.push_back("muARpresel_loose_eta");               hbins.push_back(12); hlow.push_back( 0); hup.push_back(2.4);
  histonames.push_back("muARpresel_loose_phi");               hbins.push_back(16); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("muARpresel_loose_reliso");            hbins.push_back(16); hlow.push_back( 0); hup.push_back(0.4);

  histonames.push_back("YieldsSR_tightcharge");                        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsCR_tightcharge");                        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_nisotrack");                          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsCR_nisotrack");                          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_nisotrack_tightcharge");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsCR_nisotrack_tightcharge");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSRmumue_1mhits_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCRmumue_1mhits_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSRmumue_0mhits_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCRmumue_0mhits_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSRmumue_2mhits_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCRmumue_2mhits_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("eMHits_CRmumue_invertMETDPhiOrPt");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("eMHits_SRmumue_invertMETDPhiOrPt");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("eIso_CRmumue_1mhits_invertMETDPhiOrPt");       hbins.push_back(10);hlow.push_back(0); hup.push_back(0.2);
  histonames.push_back("eIso_SRmumue_1mhits_invertMETDPhiOrPt");       hbins.push_back(10);hlow.push_back(0); hup.push_back(0.2);
  histonames.push_back("eIso_CRmumue_2mhits_invertMETDPhiOrPt");       hbins.push_back(10);hlow.push_back(0); hup.push_back(0.2);
  histonames.push_back("eIso_SRmumue_2mhits_invertMETDPhiOrPt");       hbins.push_back(10);hlow.push_back(0); hup.push_back(0.2);
  histonames.push_back("eIso_CRmumue_0mhits_invertMETDPhiOrPt");       hbins.push_back(10);hlow.push_back(0); hup.push_back(0.2);
  histonames.push_back("eIso_SRmumue_0mhits_invertMETDPhiOrPt");       hbins.push_back(10);hlow.push_back(0); hup.push_back(0.2);
  histonames.push_back("epT_CRmumue_1mhits_invertMETDPhiOrPt");        hbins.push_back(18);hlow.push_back(20);hup.push_back(200);
  histonames.push_back("epT_SRmumue_1mhits_invertMETDPhiOrPt");        hbins.push_back(18);hlow.push_back(20);hup.push_back(200);
  histonames.push_back("epT_CRmumue_0mhits_invertMETDPhiOrPt");        hbins.push_back(18);hlow.push_back(20);hup.push_back(200);
  histonames.push_back("epT_SRmumue_0mhits_invertMETDPhiOrPt");        hbins.push_back(18);hlow.push_back(20);hup.push_back(200);
  histonames.push_back("epT_CRmumue_2mhits_invertMETDPhiOrPt");        hbins.push_back(18);hlow.push_back(20);hup.push_back(200);
  histonames.push_back("epT_SRmumue_2mhits_invertMETDPhiOrPt");        hbins.push_back(18);hlow.push_back(20);hup.push_back(200);
  histonames.push_back("eIP_CRmumue_1mhits_invertMETDPhiOrPt");        hbins.push_back(15);hlow.push_back(0); hup.push_back(0.015);
  histonames.push_back("eIP_SRmumue_1mhits_invertMETDPhiOrPt");        hbins.push_back(15);hlow.push_back(0); hup.push_back(0.015);
  histonames.push_back("eIP_CRmumue_0mhits_invertMETDPhiOrPt");        hbins.push_back(15);hlow.push_back(0); hup.push_back(0.015);
  histonames.push_back("eIP_SRmumue_0mhits_invertMETDPhiOrPt");        hbins.push_back(15);hlow.push_back(0); hup.push_back(0.015);
  histonames.push_back("eIP_CRmumue_2mhits_invertMETDPhiOrPt");        hbins.push_back(15);hlow.push_back(0); hup.push_back(0.015);
  histonames.push_back("eIP_SRmumue_2mhits_invertMETDPhiOrPt");        hbins.push_back(15);hlow.push_back(0); hup.push_back(0.015);
  histonames.push_back("Mmumu_mumue_0mhits_invertMETDPhiOrPt");        hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mmumu_mumue_1mhits_invertMETDPhiOrPt");        hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mmumu_mumue_2mhits_invertMETDPhiOrPt");        hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mmumue_mumue_0mhits_invertMETDPhiOrPt");       hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mmumue_mumue_1mhits_invertMETDPhiOrPt");       hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mmumue_mumue_2mhits_invertMETDPhiOrPt");       hbins.push_back(14);hlow.push_back(22);hup.push_back(176);

    
  histonames.push_back("Mlll_SRlike_allSFOS");               hbins.push_back(15);hlow.push_back(30);hup.push_back(180);
  histonames.push_back("MSFOS_SRlike_allSFOS");              hbins.push_back(10);hlow.push_back(0); hup.push_back(200);
  histonames.push_back("minDRll_SRlike_allSFOS");            hbins.push_back(25);hlow.push_back(0.);hup.push_back(2.5);
  histonames.push_back("YieldsSR_lowMSFOS");                 hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_lowMSFOS_Mlll");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_dRllmin");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data

  histonames.push_back("Mlll_SRlike_allSFOS_invertMETDPhiOrPt");               hbins.push_back(15);hlow.push_back(30);hup.push_back(180);
  histonames.push_back("MSFOS_SRlike_allSFOS_invertMETDPhiOrPt");              hbins.push_back(10);hlow.push_back(0); hup.push_back(200);
  histonames.push_back("minDRll_SRlike_allSFOS_invertMETDPhiOrPt");            hbins.push_back(25);hlow.push_back(0.);hup.push_back(2.5);
  histonames.push_back("YieldsSR_lowMSFOS_invertMETDPhiOrPt");                 hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_lowMSFOS_Mlll_invertMETDPhiOrPt");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_dRllmin_invertMETDPhiOrPt");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_invertMETDPhiOrPt");                          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  
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
      mapname = histonames[i] + "_fakes";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_doublefakes";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_photonfakes";//new
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_photondoublefakes";//new
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_photontriplefakes";//new
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_otherphotonfakes";//new
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fakesphotonfakes";//new
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_others";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_trueWWW";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_3lLL";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_true3L";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fromWZ";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fromBC";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fromL";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fromLF";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fromG";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_fromO";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_tightFake_looseNotFake";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      mapname = histonames[i] + "_tightFake_looseFake";
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    } else {
      mapname = histonames[i] + "_"+skimFilePrefix;
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
      //histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
    }
  }
  TH2D *h2D;
  if(skimFilePrefix.find("Background")!=string::npos){
    h2D = new TH2D("GenFlags_vs_motherSSID","",8,-1.,7.,9,-4.,5.); h2D->Sumw2(); h2D->SetDirectory(rootdir);
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
    
      // Progress
      CMS3::progress( nEventsTotal, nEventsChain );

      
      double weight = evt_scale1fb()*35.9;
      //if(nEventsTotal==0) cout << weight << endl; 

      if(nVert()<0)              continue;
      if(nlep()<2)               continue;
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
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<5) ++nj;
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
      vector<int> vSS, v3l, v, iSS, i3l;
      vector<int> vaSS, va3l, va, iaSS, ia3l;
      int looseEle = -1; int veton3lspec = 0;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool islooseEle = false;
	bool isSS = false; bool is3l = false;
	bool isaSS = false; bool isa3l = false;
	if(looseEle<0&&abs(lep_pdgId()[i])==11&&lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015){
	  //don't have triggersafebla, or ID bla
	  if(fabs(lep_dxy()[i])<0.05&&fabs(lep_dz()[i])<0.1&&(fabs(lep_ip3d()[i]/lep_ip3derr()[i])<4.)&&lep_p4()[i ].Pt()>20.){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EA()[i]<0.0588)      { looseEle = i; islooseEle = true; }
	    else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EA()[i]<0.0571) { looseEle = i; islooseEle = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EA()[i]<0.0588){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EA()[i]<0.0571){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EA()[i]<0.06){
	    if(lep_p4()[i ].Pt()>20) { i3l.push_back(i); is3l = true; }
	    if(lep_p4()[i ].Pt()>30) { iSS.push_back(i); isSS = true; }
	  }
	}

	if(fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015){
	  if(lep_pass_VVV_cutbased_veto_noiso()[i]&&abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EA()[i]<0.175){
	      if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EA()[i]<0.159){
	      if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	    }
	  } else if(lep_pass_VVV_cutbased_fo_noiso()[i]&&abs(lep_pdgId()[i])==13&&lep_relIso03EA()[i]<0.4){
	    if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	    if(lep_p4()[i ].Pt()>30&&!isSS) { iaSS.push_back(i); isaSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_p4()[i ].Pt()>10&&lep_relIso03EA()[i]<=0.4) {
	  v.push_back(i);
	  if(!isSS) vSS.push_back(i);
	  if(!is3l) {
	    v3l.push_back(i);
	    if(!islooseEle) ++veton3lspec;
	  }
	  va.push_back(i);
	  if(!isSS&&!isaSS) vaSS.push_back(i);
	  if(!is3l&&!isa3l) va3l.push_back(i);
	}
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
      if(looseEle>=0) ++nel;
      if(looseEle>=0&&lep_p4()[looseEle].Pt()>25.) ++nel25;
      if(nmu>=2)           passofflineforTrigger = true;
      if(nmu25>=1&&nel>=1) passofflineforTrigger = true;
      if(nel25>=1&&nel>=2) passofflineforTrigger = true;
      if((nSS+naSS)>=2)    passofflineforTrigger = true;

      if(n3l<2&&(nSS+naSS)<2) continue;
      //if(nj30<2) continue;
      if(nb!=0) continue;
      if(!passofflineforTrigger) continue;
      if((nSS+naSS)==0&&lep_p4()[i3l[0] ].Pt()<25) continue;
      
      if(isData()){
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) continue;
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = false;
	if(nmu>=2&&(HLT_DoubleMu()||HLT_singleMu()) )             passonlineTrigger = true;
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
      else if(sn.find("Background")!=string::npos){
	bool isHtoWW = false;
	bool isWnotFromH = false;
	for(unsigned int i = 0; i<genPart_pdgId().size();++i){
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])==25) { isHtoWW = true; }
	  if(abs(genPart_pdgId()[i])==24&&abs(genPart_motherId()[i])!=25&&abs(genPart_grandmaId()[i])!=25&&genPart_status()[i]==22) { isWnotFromH = true; }
	  if(isHtoWW&&isWnotFromH) break;
	}
	if(isHtoWW&&isWnotFromH) continue;
	if(iSS.size()>=2){
	  int nW(0), nZ(0), nG(0), nF(0);
	  if(lep_isFromW()[iSS[0] ]) ++nW;
	  else if(lep_isFromZ()[iSS[0] ]) ++nZ;
	  else if(lep_isFromB()[iSS[0] ]||lep_isFromC()[iSS[0] ]||lep_isFromL()[iSS[0] ]||lep_isFromLF()[iSS[0] ]) ++nF;
	  else if(lep_motherIdSS()[iSS[0] ]) ++nG;
	  if(lep_isFromW()[iSS[1] ]) ++nW;
	  else if(lep_isFromZ()[iSS[1] ]) ++nZ;
	  else if(lep_isFromB()[iSS[1] ]||lep_isFromC()[iSS[1] ]||lep_isFromL()[iSS[1] ]||lep_isFromLF()[iSS[1] ]) ++nF;
	  else if(lep_motherIdSS()[iSS[1] ]) ++nG;
	  if(nW==2&&lep_mc_Id()[iSS[0] ]*lep_mc_Id()[iSS[1] ]>0) sn = "trueSS";//W+W+
	  else if(nW==2) sn = "chargeflips";//W+W-
	  else if(nZ==2&&lep_mc_Id()[iSS[0] ]*lep_mc_Id()[iSS[1] ]<=0) sn = "chargeflips";//Z
	  else if(nZ==2) sn = "SSLL";//ZZ both with a lost lepton
	  else if(nW==1&&nZ==1) sn = "SSLL";//WZ
	  else if((nW+nZ)==1&&nG==1) sn = "photonfakes";
	  else if((nW+nZ)==1) sn = "fakes";
	  else if((nW+nZ)==0&&nG==2) sn = "photondoublefakes";
	  else if((nW+nZ)==0&&nG==1) sn = "fakesphotonfakes";
	  else if((nW+nZ)==0) sn = "doublefakes";
	  else { 
	    if(nG>=1) sn = "otherphotonfakes";
	    else      sn = "others";
	  }
	} else if(iSS.size()==1&&iaSS.size()>=1){
	  int nW(0), nZ(0), nG(0), nF(0);
	  if(lep_isFromW()[iSS[0] ]) ++nW;
	  else if(lep_isFromZ()[iSS[0] ]) ++nZ;
	  else if(lep_isFromB()[iSS[0] ]||lep_isFromC()[iSS[0] ]||lep_isFromL()[iSS[0] ]||lep_isFromLF()[iSS[0] ]) ++nF;
	  else if(lep_motherIdSS()[iSS[0] ]) ++nG;
	  if(lep_isFromW()[iaSS[0] ]) ++nW;
	  else if(lep_isFromZ()[iaSS[0] ]) ++nZ;
	  else if(lep_isFromB()[iaSS[0] ]||lep_isFromC()[iaSS[0] ]||lep_isFromL()[iaSS[0] ]||lep_isFromLF()[iaSS[0] ]) ++nF;
	  else if(lep_motherIdSS()[iaSS[0] ]) ++nG;
	  if(nW==2&&lep_mc_Id()[iSS[0] ]*lep_mc_Id()[iaSS[0] ]>0) sn = "trueSS";//W+W+
	  else if(nW==2) sn = "chargeflips";//W+W-
	  else if(nZ==2&&lep_mc_Id()[iSS[0] ]*lep_mc_Id()[iaSS[0] ]<=0) sn = "chargeflips";//Z
	  else if(nZ==2) sn = "SSLL";//ZZ both with a lost lepton
	  else if(nW==1&&nZ==1) sn = "SSLL";//WZ
	  else if((nW+nZ)==1&&nG==1) sn = "photonfakes";
	  else if((nW+nZ)==1) sn = "fakes";
	  else if((nW+nZ)==0&&nG==2) sn = "photondoublefakes";
	  else if((nW+nZ)==0&&nG==1) sn = "fakesphotonfakes";
	  else if((nW+nZ)==0) sn = "doublefakes";
	  else { 
	    if(nG>=1) sn = "otherphotonfakes";
	    else      sn = "others";
	  }
	}
	if(i3l.size()>=3){
	  int nW(0), nZ(0), nG(0), nF(0);
	  if(lep_isFromW()[i3l[0] ]) ++nW;
	  else if(lep_isFromZ()[i3l[0] ]) ++nZ;
	  else if(lep_isFromB()[i3l[0] ]||lep_isFromC()[i3l[0] ]||lep_isFromL()[i3l[0] ]||lep_isFromLF()[i3l[0] ]) ++nF;
	  else if(lep_motherIdSS()[i3l[0] ]) ++nG;
	  if(lep_isFromW()[i3l[1] ]) ++nW;
	  else if(lep_isFromZ()[i3l[1] ]) ++nZ;
	  else if(lep_isFromB()[i3l[1] ]||lep_isFromC()[i3l[1] ]||lep_isFromL()[i3l[1] ]||lep_isFromLF()[i3l[1] ]) ++nF;
	  else if(lep_motherIdSS()[i3l[1] ]) ++nG;
	  if(lep_isFromW()[i3l[2] ]) ++nW;
	  else if(lep_isFromZ()[i3l[2] ]) ++nZ;
	  else if(lep_isFromB()[i3l[2] ]||lep_isFromC()[i3l[2] ]||lep_isFromL()[i3l[2] ]||lep_isFromLF()[i3l[2] ]) ++nF;
	  else if(lep_motherIdSS()[i3l[2] ]) ++nG;
	  if(nW==3&&(lep_mc_Id()[i3l[0] ]>0&&lep_mc_Id()[i3l[1] ]>0&&lep_mc_Id()[i3l[2] ]>0)) sn2 = "chargeflips";//W+W+W+ - it could be +++ final state, but at the end this final state will be vetoed, so if reco is ++- (e.g.), then this is a chargeflip
	  else if(nW==3&&(lep_mc_Id()[i3l[0] ]<0&&lep_mc_Id()[i3l[1] ]<0&&lep_mc_Id()[i3l[2] ]<0)) sn2 = "chargeflips";//W+W+W+ - it could be +++ final state, but at the end this final state will be vetoed, so if reco is ++- (e.g.), then this is a chargeflip
	  else if(nW==3) sn2 = "trueWWW";
	  else if(nW==2&&nZ==1) sn2 = "3lLL";//ttZ w/ LL
	  else if(nW==1&&nZ==2) sn2 = "true3L";//WZ, neglect WZZ as LL
	  else if(nZ==3) sn2 = "3lLL";//ZZ
	  else if((nW+nZ)==2) {
	    if(nG==1) sn2 = "photonfakes";
	    else      sn2 = "fakes";
	  }
	  else if((nW+nZ)==1) {
	    if(nG==2) sn2 = "photondoublefakes";
	    else if(nG==1) sn2 = "fakesphotonfakes";
	    else sn2 = "doublefakes";
	  }
	  else {
	    if(nG==3) sn2 = "photontriplefakes";
	    else if(nG>=1) sn2 = "otherphotonfakes";
	    else sn2 = "others";//could be triple fakes
	  }
	}
	else if(i3l.size()==2&&looseEle>=0){
	  int nW(0), nZ(0), nG(0), nF(0);
	  if(lep_isFromW()[i3l[0] ]) ++nW;
	  else if(lep_isFromZ()[i3l[0] ]) ++nZ;
	  else if(lep_isFromB()[i3l[0] ]||lep_isFromC()[i3l[0] ]||lep_isFromL()[i3l[0] ]||lep_isFromLF()[i3l[0] ]) ++nF;
	  else if(lep_motherIdSS()[i3l[0] ]) ++nG;
	  if(lep_isFromW()[i3l[1] ]) ++nW;
	  else if(lep_isFromZ()[i3l[1] ]) ++nZ;
	  else if(lep_isFromB()[i3l[1] ]||lep_isFromC()[i3l[1] ]||lep_isFromL()[i3l[1] ]||lep_isFromLF()[i3l[1] ]) ++nF;
	  else if(lep_motherIdSS()[i3l[1] ]) ++nG;
	  if(lep_isFromW()[looseEle]) ++nW;
	  else if(lep_isFromZ()[looseEle]) ++nZ;
	  else if(lep_isFromB()[looseEle]||lep_isFromC()[looseEle]||lep_isFromL()[looseEle]||lep_isFromLF()[looseEle]) ++nF;
	  else if(lep_motherIdSS()[looseEle]) ++nG;
	  if(nW==3&&(lep_mc_Id()[i3l[0] ]>0&&lep_mc_Id()[i3l[1] ]>0&&lep_mc_Id()[looseEle]>0)) sn2 = "chargeflips";//W+W+W+ - it could be +++ final state, but at the end this final state will be vetoed, so if reco is ++- (e.g.), then this is a chargeflip
	  else if(nW==3&&(lep_mc_Id()[i3l[0] ]<0&&lep_mc_Id()[i3l[1] ]<0&&lep_mc_Id()[looseEle]<0)) sn2 = "chargeflips";//W+W+W+ - it could be +++ final state, but at the end this final state will be vetoed, so if reco is ++- (e.g.), then this is a chargeflip
	  else if(nW==3) sn2 = "trueWWW";
	  else if(nW==2&&nZ==1) sn2 = "3lLL";//ttZ w/ LL
	  else if(nW==1&&nZ==2) sn2 = "true3L";//WZ, neglect WZZ as LL
	  else if(nZ==3) sn2 = "3lLL";//ZZ
	  else if((nW+nZ)==2) {
	    if(nG==1) sn2 = "photonfakes";
	    else      sn2 = "fakes";
	  }
	  else if((nW+nZ)==1) {
	    if(nG==2) sn2 = "photondoublefakes";
	    else if(nG==1) sn2 = "fakesphotonfakes";
	    else sn2 = "doublefakes";
	  }
	  else {
	    if(nG==3) sn2 = "photontriplefakes";
	    else if(nG>=1) sn2 = "otherphotonfakes";
	    else sn2 = "others";//could be triple fakes
	  }
	}
      }
      //cout << sn << " " << sn2 << endl;
      string lepsn1=""; string lepsn2="";
      if(iSS.size()>=1){
	if(abs(lep_pdgId()[iSS[0] ])==11) lepsn1="el";
	if(abs(lep_pdgId()[iSS[0] ])==13) lepsn1="mu";
      }
      if(iSS.size()>=2){
	if(abs(lep_pdgId()[iSS[1] ])==11) lepsn2="el";
	if(abs(lep_pdgId()[iSS[1] ])==13) lepsn2="mu";
      }
      if(iaSS.size()>=1){
	if(abs(lep_pdgId()[iaSS[0] ])==11) lepsn2="el";
	if(abs(lep_pdgId()[iaSS[0] ])==13) lepsn2="mu";
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

      bool ee[50],mm[50],em[50];
      for(int i = 0; i<50; ++i) { ee[i] = false; em[i] = false; mm[i] = false; }
      //0: SR
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                        mm[0] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[0] = false; mm[0] = false; }
      }
      //1: SRpreselect
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) ee[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11) em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13) em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) mm[1] = true;
      }
      //2: AR
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) ee[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&aMTmax>90.)                                               em[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13&&met_pt()>40.&&aMTmax>90.)                                               em[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13)                                                                         mm[2] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()<40){ ee[2] = false; mm[2] = false; }
	//if((abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11)||(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13)) cout<<met_pt()<<" "<<aMTmax<<endl;
      }
      //3: ARpreselect
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11) ee[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11) em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13) em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13) mm[3] = true;
      }
 
      bool inif = false;
      bool printout = false;
      if(ee[3]||em[3]||mm[3]||ee[1]||em[1]||mm[1]){
	inif = true;
	int second = -1;
	if(ee[3]||em[3]||mm[3]) second = iaSS[0];
	else second = iSS[1];
	int idxclose(-1), idxmatch(-1), idxphot(-1);
	int adxclose(-1), adxmatch(-1), adxphot(-1);
	for(unsigned int i = 0; i<genPart_pdgId().size(); ++i){
	  if(abs(lep_pdgId()[iSS[0] ])==abs(genPart_pdgId()[i])&&dR(lep_p4()[iSS[0] ],genPart_p4()[i])<0.3){
	    if(idxclose<0) idxclose = i;
	    else if(dR(lep_p4()[iSS[0] ],genPart_p4()[i])<dR(lep_p4()[iSS[0] ],genPart_p4()[idxclose])) idxclose = i;
	  }
	  if(abs(genPart_pdgId()[i])==22&&dR(lep_p4()[iSS[0] ],genPart_p4()[i])<0.3){
	    if(idxphot<0) idxphot = i;
	    else if(dR(lep_p4()[iSS[0] ],genPart_p4()[i])<dR(lep_p4()[iSS[0] ],genPart_p4()[idxphot])) idxphot = i;
	  }
	  if(dR(lep_p4()[iSS[0] ],genPart_p4()[i])<0.3){
	    if(printout) cout << "signal lep0 " << lep_pdgId()[iSS[0] ] << " has dR " << dR(lep_p4()[iSS[0] ],genPart_p4()[i]) << " to genpart " << genPart_pdgId()[i] << " that has status " << genPart_status()[i] << " (" << genPart_isp6status3()[i] << ") and (grand)mother " << genPart_motherId()[i] << " " << genPart_grandmaId()[i] << endl;
	    if(idxmatch<0) idxmatch = i;
	    else if(dR(lep_p4()[iSS[0] ],genPart_p4()[i])<dR(lep_p4()[iSS[0] ],genPart_p4()[idxmatch])) idxmatch = i;
	  }

	  if(abs(lep_pdgId()[second])==abs(genPart_pdgId()[i])&&dR(lep_p4()[second],genPart_p4()[i])<0.3){
	    if(adxclose<0) adxclose = i;
	    else if(dR(lep_p4()[second],genPart_p4()[i])<dR(lep_p4()[second],genPart_p4()[adxclose])) adxclose = i;
	  }
	  if(abs(genPart_pdgId()[i])==22&&dR(lep_p4()[second],genPart_p4()[i])<0.3){
	    if(adxphot<0) adxphot = i;
	    else if(dR(lep_p4()[second],genPart_p4()[i])<dR(lep_p4()[second],genPart_p4()[adxphot])) adxphot = i;
	  }
	  if(dR(lep_p4()[second],genPart_p4()[i])<0.3){
	    if(printout) cout << "signal lep1 " << lep_pdgId()[second] << " has dR " << dR(lep_p4()[second],genPart_p4()[i]) << " to genpart " << genPart_pdgId()[i] << " that has status " << genPart_status()[i] << " (" << genPart_isp6status3()[i] << ") and (grand)mother " << genPart_motherId()[i] << " " << genPart_grandmaId()[i] << endl;
	    if(adxmatch<0) adxmatch = i;
	    else if(dR(lep_p4()[second],genPart_p4()[i])<dR(lep_p4()[second],genPart_p4()[adxmatch])) adxmatch = i;
	  }
	}

	if(printout) {
	  if(abs(lep_pdgId()[iSS[0] ])==11) cout << "tight lepton is electron" << endl;
	  else cout << "tight lepton is muon" << endl;
	  cout << "motherIdSS " << lep_motherIdSS()[iSS[0] ] << " mcMatchId " << lep_mcMatchId()[iSS[0] ] << " isfromW/Z/B/C/L/LF " << lep_isFromW()[iSS[0] ] << "/"  << lep_isFromZ()[iSS[0] ] << "/" << lep_isFromB()[iSS[0] ] << "/" << lep_isFromC()[iSS[0] ] << "/" << lep_isFromL()[iSS[0] ] << "/" << lep_isFromLF()[iSS[0] ] << endl;
	  if(idxclose<0) cout << "no genlepton matched" << endl;
	  else cout << "closest genlepton (dR=" << dR(lep_p4()[iSS[0] ],genPart_p4()[idxclose]) << ") has status " << genPart_status()[idxclose] << " (" << genPart_isp6status3()[idxclose] << ") and (grand)mother " << genPart_motherId()[idxclose] << " " << genPart_grandmaId()[idxclose] << endl;
	  if(idxphot>=0) cout << "lepton has matched genphoton: (dR=" << dR(lep_p4()[iSS[0] ],genPart_p4()[idxphot]) << ") has status " << genPart_status()[idxphot] << " (" << genPart_isp6status3()[idxphot] << ") and (grand)mother " << genPart_motherId()[idxphot] << " " << genPart_grandmaId()[idxphot] << endl;
	  if(idxmatch>=0) cout << "closest arbitrarily matched genparticle (dR=" << dR(lep_p4()[iSS[0] ],genPart_p4()[idxmatch]) << ") is " << genPart_pdgId()[idxmatch] << " that has status " << genPart_status()[idxmatch] << " (" << genPart_isp6status3()[idxmatch] << ") and (grand)mother " << genPart_motherId()[idxmatch] << " " << genPart_grandmaId()[idxmatch] << endl;
	  if(ee[3]||em[3]||mm[3]) {
	    if(abs(lep_pdgId()[second])==11) cout << "loose lepton is electron" << endl;
	    else cout << "loose lepton is muon" << endl;
	  } else{
	    if(abs(lep_pdgId()[second])==11) cout << "other tight lepton is electron" << endl;
	    else cout << "other tight lepton is muon" << endl;
	  } 
	  cout << "motherIdSS " << lep_motherIdSS()[second] << " mcMatchId " << lep_mcMatchId()[second] << " isfromW/Z/B/C/L/LF " << lep_isFromW()[second] << "/"  << lep_isFromZ()[second] << "/" << lep_isFromB()[second] << "/" << lep_isFromC()[second] << "/" << lep_isFromL()[second] << "/" << lep_isFromLF()[second] << endl;
	  if(adxclose<0) cout << "no genlepton matched" << endl;
	  else cout << "closest genlepton (dR=" << dR(lep_p4()[second],genPart_p4()[adxclose]) << ") has status " << genPart_status()[adxclose] << " (" << genPart_isp6status3()[adxclose] << ") and (grand)mother " << genPart_motherId()[adxclose] << " " << genPart_grandmaId()[adxclose] << endl;
	  if(adxphot>=0) cout << "lepton has matched genphoton: (dR=" << dR(lep_p4()[second],genPart_p4()[adxphot]) << ") has status " << genPart_status()[adxphot] << " (" << genPart_isp6status3()[adxphot] << ") and (grand)mother " << genPart_motherId()[adxphot] << " " << genPart_grandmaId()[adxphot] << endl;
	  if(adxmatch>=0) cout << "closest arbitrarily matched genparticle (dR=" << dR(lep_p4()[second],genPart_p4()[adxmatch]) << ") is " << genPart_pdgId()[adxmatch] << " that has status " << genPart_status()[adxmatch] << " (" << genPart_isp6status3()[adxmatch] << ") and (grand)mother " << genPart_motherId()[adxmatch] << " " << genPart_grandmaId()[adxmatch] << endl;
	}
      }
      if(printout&&inif) cout << endl;

      int SFOS[50];
      for(int i = 0; i<50; ++i) { SFOS[i] = -1; }
      double pTlll(-1), DPhilllMET(-1); double Mmumu(-1), Mmumue(-1);
      if(i3l.size()>=3){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
      } else if(i3l.size()==2&&looseEle>=0){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[looseEle],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[looseEle ]).Pt();
      }
      bool passlowMSFOS = true;
      bool passlowMlll = true;
      if(nj<2&&nb==0&&(nveto3l==0&&n3l==3)){
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
	  if((DPhilllMET>2.7&&pTlll>60.))               SFOS[0] = 0;
	  if(!(DPhilllMET>2.7&&pTlll>60.))               SFOS[2] = 0;
	}
	if(pass1){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[0] = 1;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[2] = 1;
	}
	if(pass2){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[0] = 2;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[2] = 2;
	}
	if(pass0X){
	  if((DPhilllMET>2.7&&pTlll>60.))               SFOS[1] = 0;
	  if(!(DPhilllMET>2.7&&pTlll>60.))               SFOS[3] = 0;
	}
	if(pass1X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[1] = 1;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[3] = 1;
	}
	if(pass2X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[1] = 2;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[3] = 2;
	}
	//if(SFOScounter==1) cout << pass1 << " " << pass1X << " " << DPhilllMET << " " << pTlll << " " << met_pt() << " " << SFOScounter << endl;
      }
      int imu0(-1), imu1(-1);
      if(nj<2&&nb==0&&((veton3lspec==0&&n3l==2&&looseEle>=0)||(nveto3l==0&&n3l==3))){
	int SFOScounter = 0;
	int i0(-1),i1(-1),i2(-1);
	if(veton3lspec==0&&n3l==2&&looseEle>=0){
	  i0 = i3l[0]; i1 = i3l[1]; i2 = looseEle;
	} else {
	  if(fabs(lep_pdgId()[i3l[0] ])==13&&fabs(lep_pdgId()[i3l[1] ])==13) { i0 = i3l[0]; i1 = i3l[1]; i2 = i3l[2]; }
	  else if(fabs(lep_pdgId()[i3l[0] ])==13&&fabs(lep_pdgId()[i3l[2] ])==13) { i0 = i3l[0]; i1 = i3l[2]; i2 = i3l[1]; }
	  else if(fabs(lep_pdgId()[i3l[1] ])==13&&fabs(lep_pdgId()[i3l[2] ])==13) { i0 = i3l[1]; i1 = i3l[2]; i2 = i3l[0]; }
	  else { i0 = i3l[0]; i1 = i3l[1]; i2 = i3l[2]; }
	}
	imu0 = i0; imu1 = i1;

	bool OS01 = (lep_pdgId()[i0]*lep_pdgId()[i1]<0);
	bool OS02 = (lep_pdgId()[i0]*lep_pdgId()[i2]<0);
	bool OS12 = (lep_pdgId()[i1]*lep_pdgId()[i2]<0);
	bool SF01 = (abs(lep_pdgId()[i0])==abs(lep_pdgId()[i1]));
	bool SF02 = (abs(lep_pdgId()[i0])==abs(lep_pdgId()[i2]));
	bool SF12 = (abs(lep_pdgId()[i1])==abs(lep_pdgId()[i2]));
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	bool pass0X(false), pass1X(false), pass2X(false);
	if(abs(lep_charge()[i0]+lep_charge()[i1]+lep_charge()[i2])==3) SFOScounter = -1;
	if(lep_charge()[i0]==lep_charge()[i1]) SFOScounter = -1;
	if(fabs(lep_pdgId()[i0])!=13) SFOScounter = -1;
	if(fabs(lep_pdgId()[i1])!=13) SFOScounter = -1;
	if(fabs(lep_pdgId()[i2])!=11) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&lep_motherIdSS()[i2]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")!=string::npos&&lep_motherIdSS()[i2]==(-3)) SFOScounter = -1;
	//cout << "lep1ID " << lep_pdgId()[i0] << " lep2ID " << lep_pdgId()[i1] << " looseEle " << i2<< " SFOScounter " << SFOScounter << endl;

	//upper three lines require an mu+mu- pair
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i0]+lep_p4()[i1]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i0]+lep_p4()[i2]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i1]+lep_p4()[i2]).M()<20.) pass0 = false;
	  if(SF01&&(lep_p4()[i0]+lep_p4()[i1]).M()<20.) pass0X = false;
	  if(SF02&&(lep_p4()[i0]+lep_p4()[i2]).M()<20.) pass0X = false;
	  if(SF12&&(lep_p4()[i1]+lep_p4()[i2]).M()<20.) pass0X = false;
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i0])==11&&fabs((lep_p4()[i0]+lep_p4()[i1]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i0])==11&&fabs((lep_p4()[i0]+lep_p4()[i2]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i1])==11&&fabs((lep_p4()[i1]+lep_p4()[i2]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass0X = false;
	  else cout<<OS01<<" "<<SF01<<" ("<<(OS01&&SF01)<<") "<<OS02<<" "<<SF02<<" ("<<(OS02&&SF02)<<") "<<OS12<<" "<<SF12<<" ("<<(OS12&&SF12)<<") "<<endl;	  
	}
	if(SFOScounter==1){
	  pass1 = true;
	  pass1X = true;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i0]+lep_p4()[i1]).M()>55.&&(lep_p4()[i0]+lep_p4()[i1]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i0]+lep_p4()[i2]).M()>55.&&(lep_p4()[i0]+lep_p4()[i2]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i1]+lep_p4()[i2]).M()>55.&&(lep_p4()[i1]+lep_p4()[i2]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	  Mmumu = (lep_p4()[i0]+lep_p4()[i1]).M();
	  Mmumue = (lep_p4()[i0]+lep_p4()[i1]+lep_p4()[i2]).M();
	  if(looseEle<0) looseEle = i2;//overwrite looseEle to contain also tight ele
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i0]+lep_p4()[i1]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i0]+lep_p4()[i2]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i1]+lep_p4()[i2]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	//SR inverted DPhi,Pt, OR MET
	if(pass0){
	  if(!(DPhilllMET>2.7&&pTlll>60.))               SFOS[10] = 0;//invert either of the two
	}
	if(pass1){
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[10] = 1;//invert either of the three
	}
	if(pass2){
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[10] = 2;//invert either of the three
	}
	//SR inverted DPhi,Pt, OR MET
	if(pass0X){
	  if(!(DPhilllMET>2.7&&pTlll>60.))               SFOS[11] = 0;//invert either of the two
	}
	if(pass1X){
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[11] = 1;//invert either of the three
	}
	if(pass2X){
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[11] = 2;//invert either of the three
	}
	//if(SFOScounter==1) cout << pass1 << " " << pass1X << " " << DPhilllMET << " " << pTlll << " " << met_pt() << " " << SFOScounter << endl;
      }
      //cout << __LINE__ << endl;
      if(!isData()){
	//	if(SFOS[0]>=0) cout << "YieldsSR_"      +sn2 << endl;
	if(SFOS[0]==0)  histos["YieldsSR_"      +sn2]->Fill(3.,weight);
	if(SFOS[0]==1)  histos["YieldsSR_"      +sn2]->Fill(4.,weight);
	if(SFOS[0]==2)  histos["YieldsSR_"      +sn2]->Fill(5.,weight);
	if(SFOS[0]==0              )  histos["YieldsSR_lowMSFOS_"      +sn2]->Fill(3.,weight);
	if(SFOS[0]==1&&passlowMSFOS)  histos["YieldsSR_lowMSFOS_"      +sn2]->Fill(4.,weight);
	if(SFOS[0]==2&&passlowMSFOS)  histos["YieldsSR_lowMSFOS_"      +sn2]->Fill(5.,weight);
	if(SFOS[0]==0                           )  histos["YieldsSR_lowMSFOS_Mlll_"      +sn2]->Fill(3.,weight);
	if(SFOS[0]==1&&passlowMSFOS&&passlowMlll)  histos["YieldsSR_lowMSFOS_Mlll_"      +sn2]->Fill(4.,weight);
	if(SFOS[0]==2&&passlowMSFOS&&passlowMlll)  histos["YieldsSR_lowMSFOS_Mlll_"      +sn2]->Fill(5.,weight);
	if(SFOS[0]>=0){
	  float DRmin = 10.;
	  if(dR(lep_p4()[i3l[0] ],lep_p4()[i3l[1] ])<dR(lep_p4()[i3l[0] ],lep_p4()[i3l[2] ])&&dR(lep_p4()[i3l[0] ],lep_p4()[i3l[1] ])<dR(lep_p4()[i3l[1] ],lep_p4()[i3l[2] ]))
	    DRmin = dR(lep_p4()[i3l[0] ],lep_p4()[i3l[1] ]);
	  else if(dR(lep_p4()[i3l[0] ],lep_p4()[i3l[2] ])<dR(lep_p4()[i3l[1] ],lep_p4()[i3l[2] ]))
	    DRmin = dR(lep_p4()[i3l[0] ],lep_p4()[i3l[2] ]);
	  else
	    DRmin = dR(lep_p4()[i3l[1] ],lep_p4()[i3l[2] ]);
	  histos["minDRll_SRlike_allSFOS_"      +sn2]->Fill(DRmin,weight);
	  if(DRmin>0.2){
	    if(SFOS[0]==0)  histos["YieldsSR_dRllmin_"      +sn2]->Fill(3.,weight);
	    if(SFOS[0]==1)  histos["YieldsSR_dRllmin_"      +sn2]->Fill(4.,weight);
	    if(SFOS[0]==2)  histos["YieldsSR_dRllmin_"      +sn2]->Fill(5.,weight);
	  }
	  if(SFOS[0]==0||(SFOS[0]>=1&&passlowMSFOS)) histos["Mlll_SRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(SFOS[0]>=1){
	    bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	    bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	    bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	    bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	    bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	    bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	    if(OS01&&SF01)  histos["MSFOS_SRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	    if(OS02&&SF02)  histos["MSFOS_SRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	    if(OS12&&SF12)  histos["MSFOS_SRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	  }
	}
	if((SFOS[0]>=0||ee[0]||em[0]||mm[0])&&(skimFilePrefix.find("WG")!=string::npos||skimFilePrefix.find("ZG")!=string::npos)){
	  cout << "sample " << skimFilePrefix << " run:ls:evt " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
	  cout << "Ml1l2 " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M() << " Ml1l3 " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M() << " Ml2l3 " << (lep_p4()[i3l[2] ]+lep_p4()[i3l[1] ]).M() << " Mlll " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M() << endl;
	  for(unsigned int i = 0; i<i3l.size(); ++i){
	    cout << "lep " << i3l[i] << " id " << lep_pdgId()[i3l[i] ] << " pt " << lep_p4()[i3l[i] ].Pt() << " eta " << lep_p4()[i3l[i] ].Eta() << " phi " << lep_p4()[i3l[i] ].Phi() << " fromW/Z/B/C/L/LF " << lep_isFromW()[i3l[i] ] << "/" <<  lep_isFromZ()[i3l[i] ] << "/" <<  lep_isFromB()[i3l[i] ] << "/" <<  lep_isFromC()[i3l[i] ] << "/" <<  lep_isFromL()[i3l[i] ] << "/" <<  lep_isFromLF()[i3l[i] ] << " ssid " << lep_motherIdSS()[i3l[i] ] << " and mass " <<  lep_p4()[i3l[i] ].M() << endl;
	  }
	  for(unsigned int i = 0; i<genPart_status().size();++i){
	    if(abs(genPart_pdgId()[i])==11||abs(genPart_pdgId()[i])==13||abs(genPart_pdgId()[i])==15||abs(genPart_pdgId()[i])==22){
	      cout << "genPart " << i << " id " << genPart_pdgId()[i] << " status " << genPart_status()[i] << " pt " << genPart_p4()[i].Pt() << " eta " << genPart_p4()[i].Eta() << " phi " << genPart_p4()[i].Phi() << " motherid " << genPart_motherId()[i] << " grandmaid " << genPart_grandmaId()[i] << " and mass " << genPart_p4()[i].M() << endl;
	    }
	  }
	  cout << endl;
	}
	
	if(SFOS[0]==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsSR_nisotrack_"      +sn2]->Fill(3.,weight);
	if(SFOS[0]==1&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsSR_nisotrack_"      +sn2]->Fill(4.,weight);
	if(SFOS[0]==2&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsSR_nisotrack_"      +sn2]->Fill(5.,weight);
	if(n3l>=3){
	  bool elepasstightcharge = true;
	  if(fabs(lep_pdgId()[i3l[0] ])==11&&lep_tightCharge()[i3l[0] ]!=2) elepasstightcharge = false;
	  if(fabs(lep_pdgId()[i3l[1] ])==11&&lep_tightCharge()[i3l[1] ]!=2) elepasstightcharge = false;
	  if(fabs(lep_pdgId()[i3l[2] ])==11&&lep_tightCharge()[i3l[2] ]!=2) elepasstightcharge = false;
	  if(SFOS[0]==0&&elepasstightcharge) histos["YieldsSR_tightcharge_"      +sn2]->Fill(3.,weight);
	  if(SFOS[0]==1&&elepasstightcharge) histos["YieldsSR_tightcharge_"      +sn2]->Fill(4.,weight);
	  if(SFOS[0]==2&&elepasstightcharge) histos["YieldsSR_tightcharge_"      +sn2]->Fill(5.,weight);
	  if(SFOS[0]==0&&elepasstightcharge&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsSR_nisotrack_tightcharge_"      +sn2]->Fill(3.,weight);
	  if(SFOS[0]==1&&elepasstightcharge&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsSR_nisotrack_tightcharge_"      +sn2]->Fill(4.,weight);
	  if(SFOS[0]==2&&elepasstightcharge&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsSR_nisotrack_tightcharge_"      +sn2]->Fill(5.,weight);
	}
      }


      if(SFOS[2]==0)  histos["YieldsSR_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[2]==1)  histos["YieldsSR_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[2]==2)  histos["YieldsSR_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[2]==0              )  histos["YieldsSR_lowMSFOS_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[2]==1&&passlowMSFOS)  histos["YieldsSR_lowMSFOS_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[2]==2&&passlowMSFOS)  histos["YieldsSR_lowMSFOS_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[2]==0                           )  histos["YieldsSR_lowMSFOS_Mlll_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[2]==1&&passlowMSFOS&&passlowMlll)  histos["YieldsSR_lowMSFOS_Mlll_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[2]==2&&passlowMSFOS&&passlowMlll)  histos["YieldsSR_lowMSFOS_Mlll_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[2]>=0){
	float DRmin = 10.;
	if(dR(lep_p4()[i3l[0] ],lep_p4()[i3l[1] ])<dR(lep_p4()[i3l[0] ],lep_p4()[i3l[2] ])&&dR(lep_p4()[i3l[0] ],lep_p4()[i3l[1] ])<dR(lep_p4()[i3l[1] ],lep_p4()[i3l[2] ]))
	  DRmin = dR(lep_p4()[i3l[0] ],lep_p4()[i3l[1] ]);
	else if(dR(lep_p4()[i3l[0] ],lep_p4()[i3l[2] ])<dR(lep_p4()[i3l[1] ],lep_p4()[i3l[2] ]))
	  DRmin = dR(lep_p4()[i3l[0] ],lep_p4()[i3l[2] ]);
	else
	  DRmin = dR(lep_p4()[i3l[1] ],lep_p4()[i3l[2] ]);
	histos["minDRll_SRlike_allSFOS_invertMETDPhiOrPt_"      +sn2]->Fill(DRmin,weight);
	if(DRmin>0.2){
	  if(SFOS[2]==0)  histos["YieldsSR_dRllmin_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
	  if(SFOS[2]==1)  histos["YieldsSR_dRllmin_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
	  if(SFOS[2]==2)  histos["YieldsSR_dRllmin_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
	}
	if(SFOS[2]==0||(SFOS[2]>=1&&passlowMSFOS)) histos["Mlll_SRlike_allSFOS_invertMETDPhiOrPt_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	if(SFOS[2]>=1){
	  bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	  bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	  bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	  bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	  bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	  bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	  if(OS01&&SF01)  histos["MSFOS_SRlike_allSFOS_invertMETDPhiOrPt_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight);
	  if(OS02&&SF02)  histos["MSFOS_SRlike_allSFOS_invertMETDPhiOrPt_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight);
	  if(OS12&&SF12)  histos["MSFOS_SRlike_allSFOS_invertMETDPhiOrPt_"      +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight);
	}
      }




      
      if(SFOS[1]==0)  histos["YieldsCR_"      +sn2]->Fill(3.,weight);
      if(SFOS[1]==1)  histos["YieldsCR_"      +sn2]->Fill(4.,weight);
      if(SFOS[1]==2)  histos["YieldsCR_"      +sn2]->Fill(5.,weight);
      if(SFOS[1]==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsCR_nisotrack_"      +sn2]->Fill(3.,weight);
      if(SFOS[1]==1&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsCR_nisotrack_"      +sn2]->Fill(4.,weight);
      if(SFOS[1]==2&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsCR_nisotrack_"      +sn2]->Fill(5.,weight);
      if(n3l>=3){
	bool elepasstightcharge = true;
	if(fabs(lep_pdgId()[i3l[0] ])==11&&lep_tightCharge()[i3l[0] ]!=2) elepasstightcharge = false;
	if(fabs(lep_pdgId()[i3l[1] ])==11&&lep_tightCharge()[i3l[1] ]!=2) elepasstightcharge = false;
	if(fabs(lep_pdgId()[i3l[2] ])==11&&lep_tightCharge()[i3l[2] ]!=2) elepasstightcharge = false;
	if(SFOS[1]==0&&elepasstightcharge) histos["YieldsCR_tightcharge_"      +sn2]->Fill(3.,weight);
	if(SFOS[1]==1&&elepasstightcharge) histos["YieldsCR_tightcharge_"      +sn2]->Fill(4.,weight);
	if(SFOS[1]==2&&elepasstightcharge) histos["YieldsCR_tightcharge_"      +sn2]->Fill(5.,weight);
	if(SFOS[1]==0&&elepasstightcharge&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsCR_nisotrack_tightcharge_"      +sn2]->Fill(3.,weight);
	if(SFOS[1]==1&&elepasstightcharge&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsCR_nisotrack_tightcharge_"      +sn2]->Fill(4.,weight);
	if(SFOS[1]==2&&elepasstightcharge&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0) histos["YieldsCR_nisotrack_tightcharge_"      +sn2]->Fill(5.,weight);
      }

      if(SFOS[10]==1&&lep_lostHits()[looseEle]==0) histos["Mmumu_mumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(Mmumu ,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==1) histos["Mmumu_mumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(Mmumu ,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==2) histos["Mmumu_mumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(Mmumu ,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==0) histos["Mmumu_mumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(Mmumu ,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==1) histos["Mmumu_mumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(Mmumu ,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==2) histos["Mmumu_mumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(Mmumu ,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==0) histos["Mmumue_mumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(Mmumue,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==1) histos["Mmumue_mumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(Mmumue,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==2) histos["Mmumue_mumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(Mmumue,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==0) histos["Mmumue_mumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(Mmumue,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==1) histos["Mmumue_mumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(Mmumue,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==2) histos["Mmumue_mumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(Mmumue,weight);
      
      if(SFOS[10]==0&&lep_lostHits()[looseEle]==0) histos["YieldsSRmumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[10]==0&&lep_lostHits()[looseEle]==1) histos["YieldsSRmumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[10]==0&&lep_lostHits()[looseEle]==2) histos["YieldsSRmumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==0) histos["YieldsSRmumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==1) histos["YieldsSRmumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[10]==1&&lep_lostHits()[looseEle]==2) histos["YieldsSRmumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[10]==2&&lep_lostHits()[looseEle]==0) histos["YieldsSRmumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[10]==2&&lep_lostHits()[looseEle]==1) histos["YieldsSRmumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[10]==2&&lep_lostHits()[looseEle]==2) histos["YieldsSRmumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      //cout << __LINE__ << endl;
      if(SFOS[11]==0&&lep_lostHits()[looseEle]==0) histos["YieldsCRmumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[11]==0&&lep_lostHits()[looseEle]==1) histos["YieldsCRmumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[11]==0&&lep_lostHits()[looseEle]==2) histos["YieldsCRmumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(3.,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==0) histos["YieldsCRmumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==1) histos["YieldsCRmumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[11]==1&&lep_lostHits()[looseEle]==2) histos["YieldsCRmumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(4.,weight);
      if(SFOS[11]==2&&lep_lostHits()[looseEle]==0) histos["YieldsCRmumue_0mhits_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[11]==2&&lep_lostHits()[looseEle]==1) histos["YieldsCRmumue_1mhits_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[11]==2&&lep_lostHits()[looseEle]==2) histos["YieldsCRmumue_2mhits_invertMETDPhiOrPt_"      +sn2]->Fill(5.,weight);
      if(SFOS[10]==1){
	//cout << __LINE__ << endl;
	histos["eMHits_SRmumue_invertMETDPhiOrPt_"      +sn2]->Fill(lep_lostHits()[looseEle],weight);
	if(lep_lostHits()[looseEle]==0){
	  //cout << __LINE__ << endl;
	  histos["eIso_SRmumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_relIso03EA()[looseEle],weight);
	  histos[ "epT_SRmumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_p4()[looseEle].Pt(),   weight);
	  histos[ "eIP_SRmumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_ip3d()[looseEle],      weight);
	}
	if(lep_lostHits()[looseEle]==1){
	  //cout << __LINE__ << endl;
	  histos["eIso_SRmumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_relIso03EA()[looseEle],weight);
	  histos[ "epT_SRmumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_p4()[looseEle].Pt(),   weight);
	  histos[ "eIP_SRmumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_ip3d()[looseEle],      weight);
	}
	if(lep_lostHits()[looseEle]==2){
	  //cout << __LINE__ << endl;
	  histos["eIso_SRmumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_relIso03EA()[looseEle],weight);
	  histos[ "epT_SRmumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_p4()[looseEle].Pt(),   weight);
	  histos[ "eIP_SRmumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_ip3d()[looseEle],      weight);
	}
	//cout << __LINE__ << endl;
     }
      if(SFOS[11]==1){
	//cout << __LINE__ << endl;
	histos["eMHits_CRmumue_invertMETDPhiOrPt_"      +sn2]->Fill(lep_lostHits()[looseEle],weight);
	if(lep_lostHits()[looseEle]==0){
	  histos["eIso_CRmumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_relIso03EA()[looseEle],weight);
	  histos[ "epT_CRmumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_p4()[looseEle].Pt(),   weight);
	  histos[ "eIP_CRmumue_0mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_ip3d()[looseEle],      weight);
	}
	if(lep_lostHits()[looseEle]==1){
	  histos["eIso_CRmumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_relIso03EA()[looseEle],weight);
	  histos[ "epT_CRmumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_p4()[looseEle].Pt(),   weight);
	  histos[ "eIP_CRmumue_1mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_ip3d()[looseEle],      weight);
	}
	if(lep_lostHits()[looseEle]==2){
	  histos["eIso_CRmumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_relIso03EA()[looseEle],weight);
	  histos[ "epT_CRmumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_p4()[looseEle].Pt(),   weight);
	  histos[ "eIP_CRmumue_2mhits_invertMETDPhiOrPt_"     +sn2]->Fill(lep_ip3d()[looseEle],      weight);
	}
	//cout << __LINE__ << endl;
      }

      
      if(skimFilePrefix.find("Background")!=string::npos){
	string tempsn1=""; string tempsn2="";
	if(iSS.size()>=1){
	  if(lep_isFromW()[iSS[0] ]||lep_isFromZ()[iSS[0] ])      tempsn1="fromWZ";
	  else if(lep_isFromB()[iSS[0] ]||lep_isFromC()[iSS[0] ]) tempsn1="fromBC";
	  else if(lep_isFromL()[iSS[0] ])                         tempsn1="fromL";
	  else if(lep_isFromLF()[iSS[0] ])                        tempsn1="fromLF";
	  else if(lep_motherIdSS()[iSS[0] ]==(-3))                tempsn1="fromG";
	  else                                                    tempsn1="fromO";
	} if(iSS.size()>=2){
	  if(lep_isFromW()[iSS[1] ]||lep_isFromZ()[iSS[1] ])      tempsn2="fromWZ";
	  else if(lep_isFromB()[iSS[1] ]||lep_isFromC()[iSS[1] ]) tempsn2="fromBC";
	  else if(lep_isFromL()[iSS[1] ])                         tempsn2="fromL";
	  else if(lep_isFromLF()[iSS[1] ])                        tempsn2="fromLF";
	  else if(lep_motherIdSS()[iSS[1] ]==(-3))                tempsn2="fromG";
	  else                                                    tempsn2="fromO";
	} else if(iaSS.size()>=1){
	  if(lep_isFromW()[iaSS[0] ]||lep_isFromZ()[iaSS[0] ])      tempsn2="fromWZ";
	  else if(lep_isFromB()[iaSS[0] ]||lep_isFromC()[iaSS[0] ]) tempsn2="fromBC";
	  else if(lep_isFromL()[iaSS[0] ])                          tempsn2="fromL";
	  else if(lep_isFromLF()[iaSS[0] ])                         tempsn2="fromLF";
	  else if(lep_motherIdSS()[iaSS[0] ]==(-3))                 tempsn2="fromG";
	  else                                                      tempsn2="fromO";
	}
	if(SFOS[10]==1||SFOS[11]==1){
	  string tmp = "";
	  if(lep_isFromW()[looseEle]||lep_isFromZ()[looseEle])      tmp="fromWZ";
	  else if(lep_isFromB()[looseEle]||lep_isFromC()[looseEle]) tmp="fromBC";
	  else if(lep_isFromL()[looseEle])                          tmp="fromL";
	  else if(lep_isFromLF()[looseEle])                         tmp="fromLF";
	  else if(lep_motherIdSS()[looseEle]==(-3))                 tmp="fromG";
	  else                                                      tmp="fromO";
	  if(SFOS[10]==1){
	    histos["eMHits_SRmumue_invertMETDPhiOrPt_"      +tmp]->Fill(lep_lostHits()[looseEle],weight);
	    if(lep_lostHits()[looseEle]==0){
	      histos["eIso_SRmumue_0mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_relIso03EA()[looseEle],weight);
	      histos[ "epT_SRmumue_0mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_p4()[looseEle].Pt(),   weight);
	      histos[ "eIP_SRmumue_0mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_ip3d()[looseEle],      weight);
	    }
	    if(lep_lostHits()[looseEle]==1){
	      histos["eIso_SRmumue_1mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_relIso03EA()[looseEle],weight);
	      histos[ "epT_SRmumue_1mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_p4()[looseEle].Pt(),   weight);
	      histos[ "eIP_SRmumue_1mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_ip3d()[looseEle],      weight);
	    }
	    if(lep_lostHits()[looseEle]==2){
	      histos["eIso_SRmumue_2mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_relIso03EA()[looseEle],weight);
	      histos[ "epT_SRmumue_2mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_p4()[looseEle].Pt(),   weight);
	      histos[ "eIP_SRmumue_2mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_ip3d()[looseEle],      weight);
	    }
	  }
	  if(SFOS[11]==1){
	    histos["eMHits_CRmumue_invertMETDPhiOrPt_"      +tmp]->Fill(lep_lostHits()[looseEle],weight);
	    if(lep_lostHits()[looseEle]==0){
	      histos["eIso_CRmumue_0mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_relIso03EA()[looseEle],weight);
	      histos[ "epT_CRmumue_0mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_p4()[looseEle].Pt(),   weight);
	      histos[ "eIP_CRmumue_0mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_ip3d()[looseEle],      weight);
	    }
	    if(lep_lostHits()[looseEle]==1){
	      histos["eIso_CRmumue_1mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_relIso03EA()[looseEle],weight);
	      histos[ "epT_CRmumue_1mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_p4()[looseEle].Pt(),   weight);
	      histos[ "eIP_CRmumue_1mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_ip3d()[looseEle],      weight);
	    }
	    if(lep_lostHits()[looseEle]==2){
	      histos["eIso_CRmumue_2mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_relIso03EA()[looseEle],weight);
	      histos[ "epT_CRmumue_2mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_p4()[looseEle].Pt(),   weight);
	      histos[ "eIP_CRmumue_2mhits_invertMETDPhiOrPt_"     +tmp]->Fill(lep_ip3d()[looseEle],      weight);
	    }
	  }
	}

	
	if((ee[0]||em[0]||mm[0])&&!isData()){
	  histos[lepsn1+"SR_motherIDSS_"+tempsn1]->Fill(lep_motherIdSS()[iSS[0] ],weight);
	  histos[lepsn2+"SR_motherIDSS_"+tempsn2]->Fill(lep_motherIdSS()[iSS[1] ],weight);
	  if(tempsn1=="fromWZ"){
	    if(ee[0]) histos["YieldsSR_"+tempsn2]->Fill(0.,weight);
	    if(em[0]) histos["YieldsSR_"+tempsn2]->Fill(1.,weight);
	    if(mm[0]) histos["YieldsSR_"+tempsn2]->Fill(2.,weight);
	  } else {
	    if(tempsn2=="fromWZ"){
	      if(ee[0]) histos["YieldsSR_"+tempsn1]->Fill(0.,weight);
	      if(em[0]) histos["YieldsSR_"+tempsn1]->Fill(1.,weight);
	      if(mm[0]) histos["YieldsSR_"+tempsn1]->Fill(2.,weight);
	    } else {
	      if(ee[0]) histos["YieldsSR_tightFake_looseFake"]->Fill(0.,weight);
	      if(em[0]) histos["YieldsSR_tightFake_looseFake"]->Fill(1.,weight);
	      if(mm[0]) histos["YieldsSR_tightFake_looseFake"]->Fill(2.,weight);   
	    }
	  } 	  
	}
	if((ee[1]||em[1]||mm[1])&&!isData()){
	  histos[lepsn1+"SRpresel_motherIDSS_"+tempsn1]->Fill(lep_motherIdSS()[iSS[0] ],weight);
	  histos[lepsn2+"SRpresel_motherIDSS_"+tempsn2]->Fill(lep_motherIdSS()[iSS[1] ],weight);
	  if(tempsn1=="fromWZ"){
	    if(ee[1]) histos["YieldsSRpresel_"+tempsn2]->Fill(0.,weight);
	    if(em[1]) histos["YieldsSRpresel_"+tempsn2]->Fill(1.,weight);
	    if(mm[1]) histos["YieldsSRpresel_"+tempsn2]->Fill(2.,weight);
	  } else {
	    if(tempsn2=="fromWZ"){
	      if(ee[1]) histos["YieldsSRpresel_"+tempsn1]->Fill(0.,weight);
	      if(em[1]) histos["YieldsSRpresel_"+tempsn1]->Fill(1.,weight);
	      if(mm[1]) histos["YieldsSRpresel_"+tempsn1]->Fill(2.,weight);
	    } else {
	      if(ee[1]) histos["YieldsSRpresel_tightFake_looseFake"]->Fill(0.,weight);
	      if(em[1]) histos["YieldsSRpresel_tightFake_looseFake"]->Fill(1.,weight);
	      if(mm[1]) histos["YieldsSRpresel_tightFake_looseFake"]->Fill(2.,weight);   
	    }
	  }
	}
	if(ee[2]||em[2]||mm[2]){
	  histos[lepsn1+"AR_motherIDSS_"+tempsn1]->Fill(lep_motherIdSS()[iSS[0] ], weight);
	  histos[lepsn2+"AR_motherIDSS_"+tempsn2]->Fill(lep_motherIdSS()[iaSS[0] ],weight);
	  histos[lepsn1+"AR_tight_pT_"    +tempsn1]->Fill(lep_p4()[iSS[0] ].Pt(),        weight);
	  histos[lepsn1+"AR_tight_eta_"   +tempsn1]->Fill(fabs(lep_p4()[iSS[0] ].Eta()), weight);
	  histos[lepsn1+"AR_tight_phi_"   +tempsn1]->Fill(fabs(lep_p4()[iSS[0] ].Phi()), weight);
	  histos[lepsn1+"AR_tight_reliso_"+tempsn1]->Fill(lep_relIso03EA()[iSS[0] ],     weight);
	  histos[lepsn2+"AR_loose_pT_"    +tempsn2]->Fill(lep_p4()[iaSS[0] ].Pt(),       weight);
	  histos[lepsn2+"AR_loose_eta_"   +tempsn2]->Fill(fabs(lep_p4()[iaSS[0] ].Eta()),weight);
	  histos[lepsn2+"AR_loose_phi_"   +tempsn2]->Fill(fabs(lep_p4()[iaSS[0] ].Phi()),weight);
	  histos[lepsn2+"AR_loose_reliso_"+tempsn2]->Fill(lep_relIso03EA()[iaSS[0] ],    weight);
	  if(tempsn1.find("fromWZ")!=string::npos){
	    if(ee[2]) histos["YieldsAR_"+tempsn2]->Fill(0.,weight);
	    if(em[2]) histos["YieldsAR_"+tempsn2]->Fill(1.,weight);
	    if(mm[2]) histos["YieldsAR_"+tempsn2]->Fill(2.,weight);
	  } else {
	    if(tempsn2.find("fromWZ")!=string::npos){
	      if(ee[2]) histos["YieldsAR_tightFake_looseNotFake"]->Fill(0.,weight);
	      if(em[2]) histos["YieldsAR_tightFake_looseNotFake"]->Fill(1.,weight);
	      if(mm[2]) histos["YieldsAR_tightFake_looseNotFake"]->Fill(2.,weight);
	    } else {
	      if(ee[2]) histos["YieldsAR_tightFake_looseFake"]->Fill(0.,weight);
	      if(em[2]) histos["YieldsAR_tightFake_looseFake"]->Fill(1.,weight);
	      if(mm[2]) histos["YieldsAR_tightFake_looseFake"]->Fill(2.,weight);   
	    }
	  } 
	}
	if(ee[3]||em[3]||mm[3]){
	  if(lep_isFromW()[iSS[0] ])       h2D->Fill( 1.,lep_motherIdSS()[iSS[0] ], weight);
	  else if(lep_isFromZ()[iSS[0] ])  h2D->Fill( 2.,lep_motherIdSS()[iSS[0] ], weight);
	  else if(lep_isFromB()[iSS[0] ])  h2D->Fill( 3.,lep_motherIdSS()[iSS[0] ], weight);
	  else if(lep_isFromC()[iSS[0] ])  h2D->Fill( 4.,lep_motherIdSS()[iSS[0] ], weight);
	  else if(lep_isFromL()[iSS[0] ])  h2D->Fill( 5.,lep_motherIdSS()[iSS[0] ], weight);
	  else if(lep_isFromLF()[iSS[0] ]) h2D->Fill( 6.,lep_motherIdSS()[iSS[0] ], weight);
	  else                             h2D->Fill(-1.,lep_motherIdSS()[iSS[0] ], weight);
	  if(lep_isFromW()[iaSS[0] ])       h2D->Fill( 1.,lep_motherIdSS()[iaSS[0] ], weight);
	  else if(lep_isFromZ()[iaSS[0] ])  h2D->Fill( 2.,lep_motherIdSS()[iaSS[0] ], weight);
	  else if(lep_isFromB()[iaSS[0] ])  h2D->Fill( 3.,lep_motherIdSS()[iaSS[0] ], weight);
	  else if(lep_isFromC()[iaSS[0] ])  h2D->Fill( 4.,lep_motherIdSS()[iaSS[0] ], weight);
	  else if(lep_isFromL()[iaSS[0] ])  h2D->Fill( 5.,lep_motherIdSS()[iaSS[0] ], weight);
	  else if(lep_isFromLF()[iaSS[0] ]) h2D->Fill( 6.,lep_motherIdSS()[iaSS[0] ], weight);
	  else                              h2D->Fill(-1.,lep_motherIdSS()[iaSS[0] ], weight);
	  histos[lepsn1+"ARpresel_motherIDSS_"+tempsn1]->Fill(lep_motherIdSS()[iSS[0] ], weight);
	  histos[lepsn2+"ARpresel_motherIDSS_"+tempsn2]->Fill(lep_motherIdSS()[iaSS[0] ],weight);
	  histos[lepsn1+"ARpresel_tight_pT_"    +tempsn1]->Fill(lep_p4()[iSS[0] ].Pt(),        weight);
	  histos[lepsn1+"ARpresel_tight_eta_"   +tempsn1]->Fill(fabs(lep_p4()[iSS[0] ].Eta()), weight);
	  histos[lepsn1+"ARpresel_tight_phi_"   +tempsn1]->Fill(fabs(lep_p4()[iSS[0] ].Phi()), weight);
	  histos[lepsn1+"ARpresel_tight_reliso_"+tempsn1]->Fill(lep_relIso03EA()[iSS[0] ],     weight);
	  histos[lepsn2+"ARpresel_loose_pT_"    +tempsn2]->Fill(lep_p4()[iaSS[0] ].Pt(),       weight);
	  histos[lepsn2+"ARpresel_loose_eta_"   +tempsn2]->Fill(fabs(lep_p4()[iaSS[0] ].Eta()),weight);
	  histos[lepsn2+"ARpresel_loose_phi_"   +tempsn2]->Fill(fabs(lep_p4()[iaSS[0] ].Phi()),weight);
	  histos[lepsn2+"ARpresel_loose_reliso_"+tempsn2]->Fill(lep_relIso03EA()[iaSS[0] ],    weight);
	  if(tempsn1.find("fromWZ")!=string::npos){
	    if(ee[3]) histos["YieldsARpresel_"+tempsn2]->Fill(0.,weight);
	    if(em[3]) histos["YieldsARpresel_"+tempsn2]->Fill(1.,weight);
	    if(mm[3]) histos["YieldsARpresel_"+tempsn2]->Fill(2.,weight);
	  } else {
	    if(tempsn2.find("fromWZ")!=string::npos){
	      if(ee[3]) histos["YieldsARpresel_tightFake_looseNotFake"]->Fill(0.,weight);
	      if(em[3]) histos["YieldsARpresel_tightFake_looseNotFake"]->Fill(1.,weight);
	      if(mm[3]) histos["YieldsARpresel_tightFake_looseNotFake"]->Fill(2.,weight);
	    } else {
	      if(ee[3]) histos["YieldsARpresel_tightFake_looseFake"]->Fill(0.,weight);
	      if(em[3]) histos["YieldsARpresel_tightFake_looseFake"]->Fill(1.,weight);
	      if(mm[3]) histos["YieldsARpresel_tightFake_looseFake"]->Fill(2.,weight);   
	    }
	  }
	}

	if(ee[1]&&!isData()){
	  histos["SRee_cuts_"+tempsn1]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.)                                                          histos["SRee_cuts_"+tempsn1]->Fill( 1.,weight*0.5);
	  if(fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.)                                                histos["SRee_cuts_"+tempsn1]->Fill( 2.,weight*0.5);
	  if(met_pt()>40.)                                                                                           histos["SRee_cuts_"+tempsn1]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                             histos["SRee_cuts_"+tempsn1]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                histos["SRee_cuts_"+tempsn1]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) histos["SRee_cuts_"+tempsn1]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&met_pt()>40.)                                            histos["SRee_cuts_"+tempsn1]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                            histos["SRee_cuts_"+tempsn1]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj)                                             histos["SRee_cuts_"+tempsn1]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                               histos["SRee_cuts_"+tempsn1]->Fill(10.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj&&met_pt()>40.)                               histos["SRee_cuts_"+tempsn1]->Fill(11.,weight*0.5);
	  if(ee[0]) histos["SRee_cuts_"+tempsn1]->Fill(12.,weight*0.5);
	  histos["SRee_cuts_"+tempsn2]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.)                                                          histos["SRee_cuts_"+tempsn2]->Fill( 1.,weight*0.5);
	  if(fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.)                                                histos["SRee_cuts_"+tempsn2]->Fill( 2.,weight*0.5);
	  if(met_pt()>40.)                                                                                           histos["SRee_cuts_"+tempsn2]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                             histos["SRee_cuts_"+tempsn2]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                histos["SRee_cuts_"+tempsn2]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) histos["SRee_cuts_"+tempsn2]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&met_pt()>40.)                                            histos["SRee_cuts_"+tempsn2]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                            histos["SRee_cuts_"+tempsn2]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj)                                             histos["SRee_cuts_"+tempsn2]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                               histos["SRee_cuts_"+tempsn2]->Fill(10.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj&&met_pt()>40.)                               histos["SRee_cuts_"+tempsn2]->Fill(11.,weight*0.5);
	  if(ee[0]) histos["SRee_cuts_"+tempsn2]->Fill(12.,weight*0.5);
	}
	if(em[1]&&!isData()){
	  histos["SRem_cuts_"+tempsn1]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.)                                                          histos["SRem_cuts_"+tempsn1]->Fill( 1.,weight*0.5);
	  if(met_pt()>40.)                                                                                           histos["SRem_cuts_"+tempsn1]->Fill( 2.,weight*0.5);
	  if(MTmax>90.)                                                                                              histos["SRem_cuts_"+tempsn1]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                             histos["SRem_cuts_"+tempsn1]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                histos["SRem_cuts_"+tempsn1]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&met_pt()>40.)                                            histos["SRem_cuts_"+tempsn1]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.)                                               histos["SRem_cuts_"+tempsn1]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                            histos["SRem_cuts_"+tempsn1]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj)                                             histos["SRem_cuts_"+tempsn1]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                               histos["SRem_cuts_"+tempsn1]->Fill(10.,weight*0.5);
	  if(MTmax>90.&&passMDetajj)                                                                                 histos["SRem_cuts_"+tempsn1]->Fill(11.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj&&met_pt()>40.)                               histos["SRem_cuts_"+tempsn1]->Fill(12.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj&&MTmax>90.)                                  histos["SRem_cuts_"+tempsn1]->Fill(13.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.&&met_pt()>40.)                                 histos["SRem_cuts_"+tempsn1]->Fill(14.,weight*0.5);
	  if(MTmax>90.&&passMDetajj&&met_pt()>40.)                                                                   histos["SRem_cuts_"+tempsn1]->Fill(15.,weight*0.5);
	  if(em[0]) histos["SRem_cuts_"+tempsn1]->Fill(16.,weight*0.5);
	  histos["SRem_cuts_"+tempsn2]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.)                                                          histos["SRem_cuts_"+tempsn2]->Fill( 1.,weight*0.5);
	  if(met_pt()>40.)                                                                                           histos["SRem_cuts_"+tempsn2]->Fill( 2.,weight*0.5);
	  if(MTmax>90.)                                                                                              histos["SRem_cuts_"+tempsn2]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                             histos["SRem_cuts_"+tempsn2]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                histos["SRem_cuts_"+tempsn2]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&met_pt()>40.)                                            histos["SRem_cuts_"+tempsn2]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.)                                               histos["SRem_cuts_"+tempsn2]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                            histos["SRem_cuts_"+tempsn2]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj)                                             histos["SRem_cuts_"+tempsn2]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                               histos["SRem_cuts_"+tempsn2]->Fill(10.,weight*0.5);
	  if(MTmax>90.&&passMDetajj)                                                                                 histos["SRem_cuts_"+tempsn2]->Fill(11.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj&&met_pt()>40.)                               histos["SRem_cuts_"+tempsn2]->Fill(12.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj&&MTmax>90.)                                  histos["SRem_cuts_"+tempsn2]->Fill(13.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.&&met_pt()>40.)                                 histos["SRem_cuts_"+tempsn2]->Fill(14.,weight*0.5);
	  if(MTmax>90.&&passMDetajj&&met_pt()>40.)                                                                   histos["SRem_cuts_"+tempsn2]->Fill(15.,weight*0.5);
	  if(em[0]) histos["SRem_cuts_"+tempsn2]->Fill(16.,weight*0.5);
	}
	if(ee[3]){
	  histos["ARee_cuts_"+tempsn1]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.)                                                           histos["ARee_cuts_"+tempsn1]->Fill( 1.,weight*0.5);
	  if(fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.)                                                 histos["ARee_cuts_"+tempsn1]->Fill( 2.,weight*0.5);
	  if(met_pt()>40.)                                                                                             histos["ARee_cuts_"+tempsn1]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                               histos["ARee_cuts_"+tempsn1]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                  histos["ARee_cuts_"+tempsn1]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) histos["ARee_cuts_"+tempsn1]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&met_pt()>40.)                                             histos["ARee_cuts_"+tempsn1]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                              histos["ARee_cuts_"+tempsn1]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&passMDetajj)                                              histos["ARee_cuts_"+tempsn1]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                                 histos["ARee_cuts_"+tempsn1]->Fill(10.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&passMDetajj&&met_pt()>40.)                                histos["ARee_cuts_"+tempsn1]->Fill(11.,weight*0.5);
	  if(ee[2]) histos["ARee_cuts_"+tempsn1]->Fill(12.,weight*0.5);
	  histos["ARee_cuts_"+tempsn2]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.)                                                           histos["ARee_cuts_"+tempsn2]->Fill( 1.,weight*0.5);
	  if(fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.)                                                 histos["ARee_cuts_"+tempsn2]->Fill( 2.,weight*0.5);
	  if(met_pt()>40.)                                                                                             histos["ARee_cuts_"+tempsn2]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                               histos["ARee_cuts_"+tempsn2]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                  histos["ARee_cuts_"+tempsn2]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) histos["ARee_cuts_"+tempsn2]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&met_pt()>40.)                                             histos["ARee_cuts_"+tempsn2]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                              histos["ARee_cuts_"+tempsn2]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&passMDetajj)                                              histos["ARee_cuts_"+tempsn2]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                                 histos["ARee_cuts_"+tempsn2]->Fill(10.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&passMDetajj&&met_pt()>40.)                                histos["ARee_cuts_"+tempsn2]->Fill(11.,weight*0.5);
	  if(ee[2]) histos["ARee_cuts_"+tempsn2]->Fill(12.,weight*0.5);
	}
	if(em[3]){
	  histos["ARem_cuts_"+tempsn1]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.)                                                         histos["ARem_cuts_"+tempsn1]->Fill( 1.,weight*0.5);
	  if(met_pt()>40.)                                                                                           histos["ARem_cuts_"+tempsn1]->Fill( 2.,weight*0.5);
	  if(aMTmax>90.)                                                                                             histos["ARem_cuts_"+tempsn1]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                             histos["ARem_cuts_"+tempsn1]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                histos["ARem_cuts_"+tempsn1]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&met_pt()>40.)                                           histos["ARem_cuts_"+tempsn1]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&aMTmax>90.)                                             histos["ARem_cuts_"+tempsn1]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                            histos["ARem_cuts_"+tempsn1]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj)                                            histos["ARem_cuts_"+tempsn1]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                               histos["ARem_cuts_"+tempsn1]->Fill(10.,weight*0.5);
	  if(aMTmax>90.&&passMDetajj)                                                                                histos["ARem_cuts_"+tempsn1]->Fill(11.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj&&met_pt()>40.)                              histos["ARem_cuts_"+tempsn1]->Fill(12.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj&&aMTmax>90.)                                histos["ARem_cuts_"+tempsn1]->Fill(13.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&aMTmax>90.&&met_pt()>40.)                               histos["ARem_cuts_"+tempsn1]->Fill(14.,weight*0.5);
	  if(aMTmax>90.&&passMDetajj&&met_pt()>40.)                                                                  histos["ARem_cuts_"+tempsn1]->Fill(15.,weight*0.5);
	  if(em[2]) histos["ARem_cuts_"+tempsn1]->Fill(16.,weight*0.5);
	  histos["ARem_cuts_"+tempsn2]->Fill(0.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.)                                                         histos["ARem_cuts_"+tempsn2]->Fill( 1.,weight*0.5);
	  if(met_pt()>40.)                                                                                           histos["ARem_cuts_"+tempsn2]->Fill( 2.,weight*0.5);
	  if(aMTmax>90.)                                                                                             histos["ARem_cuts_"+tempsn2]->Fill( 3.,weight*0.5);
	  if(passDetajj)                                                                                             histos["ARem_cuts_"+tempsn2]->Fill( 4.,weight*0.5);
	  if(passMjj)                                                                                                histos["ARem_cuts_"+tempsn2]->Fill( 5.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&met_pt()>40.)                                           histos["ARem_cuts_"+tempsn2]->Fill( 6.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&aMTmax>90.)                                             histos["ARem_cuts_"+tempsn2]->Fill( 7.,weight*0.5);
	  if(passMDetajj)                                                                                            histos["ARem_cuts_"+tempsn2]->Fill( 8.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj)                                            histos["ARem_cuts_"+tempsn2]->Fill( 9.,weight*0.5);
	  if(met_pt()>40&&passMDetajj)                                                                               histos["ARem_cuts_"+tempsn2]->Fill(10.,weight*0.5);
	  if(aMTmax>90.&&passMDetajj)                                                                                histos["ARem_cuts_"+tempsn2]->Fill(11.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj&&met_pt()>40.)                              histos["ARem_cuts_"+tempsn2]->Fill(12.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj&&aMTmax>90.)                                histos["ARem_cuts_"+tempsn2]->Fill(13.,weight*0.5);
	  if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&aMTmax>90.&&met_pt()>40.)                               histos["ARem_cuts_"+tempsn2]->Fill(14.,weight*0.5);
	  if(aMTmax>90.&&passMDetajj&&met_pt()>40.)                                                                  histos["ARem_cuts_"+tempsn2]->Fill(15.,weight*0.5);
	  if(em[2]) histos["ARem_cuts_"+tempsn2]->Fill(16.,weight*0.5);
	}
      }//background

      //if(ee[0]||em[0]||mm[0]) cout << __LINE__ <<  " " << "YieldsSR_"      +sn << endl;
      if(ee[0]&&!isData()) histos["YieldsSR_"      +sn]->Fill(0.,weight);
      if(em[0]&&!isData()) histos["YieldsSR_"      +sn]->Fill(1.,weight);
      if(mm[0]&&!isData()) histos["YieldsSR_"      +sn]->Fill(2.,weight);
      if(ee[1]&&!isData()) histos["YieldsSRpresel_"+sn]->Fill(0.,weight);
      if(em[1]&&!isData()) histos["YieldsSRpresel_"+sn]->Fill(1.,weight);
      if(mm[1]&&!isData()) histos["YieldsSRpresel_"+sn]->Fill(2.,weight);
      if(ee[2])            histos["YieldsAR_"      +sn]->Fill(0.,weight);
      if(em[2])            histos["YieldsAR_"      +sn]->Fill(1.,weight);
      if(mm[2])            histos["YieldsAR_"      +sn]->Fill(2.,weight);
      if(ee[3])            histos["YieldsARpresel_"+sn]->Fill(0.,weight);
      if(em[3])            histos["YieldsARpresel_"+sn]->Fill(1.,weight);
      if(mm[3])            histos["YieldsARpresel_"+sn]->Fill(2.,weight);

      if(ee[1]&&!isData()){
	histos["SRee_cuts_"+sn]->Fill(0.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.)                                                          histos["SRee_cuts_"+sn]->Fill( 1.,weight);
	if(fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.)                                                histos["SRee_cuts_"+sn]->Fill( 2.,weight);
	if(met_pt()>40.)                                                                                           histos["SRee_cuts_"+sn]->Fill( 3.,weight);
	if(passDetajj)                                                                                             histos["SRee_cuts_"+sn]->Fill( 4.,weight);
	if(passMjj)                                                                                                histos["SRee_cuts_"+sn]->Fill( 5.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) histos["SRee_cuts_"+sn]->Fill( 6.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&met_pt()>40.)                                            histos["SRee_cuts_"+sn]->Fill( 7.,weight);
	if(passMDetajj)                                                                                            histos["SRee_cuts_"+sn]->Fill( 8.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj)                                             histos["SRee_cuts_"+sn]->Fill( 9.,weight);
	if(met_pt()>40&&passMDetajj)                                                                               histos["SRee_cuts_"+sn]->Fill(10.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.&&passMDetajj&&met_pt()>40.)                               histos["SRee_cuts_"+sn]->Fill(11.,weight);
	if(ee[0]) histos["SRee_cuts_"+sn]->Fill(12.,weight);
      }
      if(em[1]&&!isData()){
	histos["SRem_cuts_"+sn]->Fill(0.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.)                                                          histos["SRem_cuts_"+sn]->Fill( 1.,weight);
	if(met_pt()>40.)                                                                                           histos["SRem_cuts_"+sn]->Fill( 2.,weight);
	if(MTmax>90.)                                                                                              histos["SRem_cuts_"+sn]->Fill( 3.,weight);
	if(passDetajj)                                                                                             histos["SRem_cuts_"+sn]->Fill( 4.,weight);
	if(passMjj)                                                                                                histos["SRem_cuts_"+sn]->Fill( 5.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&met_pt()>40.)                                            histos["SRem_cuts_"+sn]->Fill( 6.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.)                                               histos["SRem_cuts_"+sn]->Fill( 7.,weight);
	if(passMDetajj)                                                                                            histos["SRem_cuts_"+sn]->Fill( 8.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj)                                             histos["SRem_cuts_"+sn]->Fill( 9.,weight);
	if(met_pt()>40&&passMDetajj)                                                                               histos["SRem_cuts_"+sn]->Fill(10.,weight);
	if(MTmax>90.&&passMDetajj)                                                                                 histos["SRem_cuts_"+sn]->Fill(11.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj&&met_pt()>40.)                               histos["SRem_cuts_"+sn]->Fill(12.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&passMDetajj&&MTmax>90.)                                  histos["SRem_cuts_"+sn]->Fill(13.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.&&MTmax>90.&&met_pt()>40.)                                 histos["SRem_cuts_"+sn]->Fill(14.,weight);
	if(MTmax>90.&&passMDetajj&&met_pt()>40.)                                                                   histos["SRem_cuts_"+sn]->Fill(15.,weight);
	if(em[0]) histos["SRem_cuts_"+sn]->Fill(16.,weight);
      }
      if(ee[3]){
	histos["ARee_cuts_"+sn]->Fill(0.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.)                                                           histos["ARee_cuts_"+sn]->Fill( 1.,weight);
	if(fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.)                                                 histos["ARee_cuts_"+sn]->Fill( 2.,weight);
	if(met_pt()>40.)                                                                                             histos["ARee_cuts_"+sn]->Fill( 3.,weight);
	if(passDetajj)                                                                                               histos["ARee_cuts_"+sn]->Fill( 4.,weight);
	if(passMjj)                                                                                                  histos["ARee_cuts_"+sn]->Fill( 5.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) histos["ARee_cuts_"+sn]->Fill( 6.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&met_pt()>40.)                                             histos["ARee_cuts_"+sn]->Fill( 7.,weight);
	if(passMDetajj)                                                                                              histos["ARee_cuts_"+sn]->Fill( 8.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&passMDetajj)                                              histos["ARee_cuts_"+sn]->Fill( 9.,weight);
	if(met_pt()>40&&passMDetajj)                                                                                 histos["ARee_cuts_"+sn]->Fill(10.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>40.&&passMDetajj&&met_pt()>40.)                                histos["ARee_cuts_"+sn]->Fill(11.,weight);
	if(ee[2]) histos["ARee_cuts_"+sn]->Fill(12.,weight);
      }
      if(em[3]){
	histos["ARem_cuts_"+sn]->Fill(0.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.)                                                         histos["ARem_cuts_"+sn]->Fill( 1.,weight);
	if(met_pt()>40.)                                                                                           histos["ARem_cuts_"+sn]->Fill( 2.,weight);
	if(aMTmax>90.)                                                                                             histos["ARem_cuts_"+sn]->Fill( 3.,weight);
	if(passDetajj)                                                                                             histos["ARem_cuts_"+sn]->Fill( 4.,weight);
	if(passMjj)                                                                                                histos["ARem_cuts_"+sn]->Fill( 5.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&met_pt()>40.)                                           histos["ARem_cuts_"+sn]->Fill( 6.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&aMTmax>90.)                                             histos["ARem_cuts_"+sn]->Fill( 7.,weight);
	if(passMDetajj)                                                                                            histos["ARem_cuts_"+sn]->Fill( 8.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj)                                            histos["ARem_cuts_"+sn]->Fill( 9.,weight);
	if(met_pt()>40&&passMDetajj)                                                                               histos["ARem_cuts_"+sn]->Fill(10.,weight);
	if(aMTmax>90.&&passMDetajj)                                                                                histos["ARem_cuts_"+sn]->Fill(11.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj&&met_pt()>40.)                              histos["ARem_cuts_"+sn]->Fill(12.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&passMDetajj&&aMTmax>90.)                                histos["ARem_cuts_"+sn]->Fill(13.,weight);
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.&&aMTmax>90.&&met_pt()>40.)                               histos["ARem_cuts_"+sn]->Fill(14.,weight);
	if(aMTmax>90.&&passMDetajj&&met_pt()>40.)                                                                  histos["ARem_cuts_"+sn]->Fill(15.,weight);
	if(em[2]) histos["ARem_cuts_"+sn]->Fill(16.,weight);
      }

      if((ee[0]||em[0]||mm[0])&&!isData()){
	histos[lepsn1+"SR_motherIDSS_"+sn]->Fill(lep_motherIdSS()[iSS[0] ],weight);
	histos[lepsn2+"SR_motherIDSS_"+sn]->Fill(lep_motherIdSS()[iSS[1] ],weight);
	if(lep_isFromW()[iSS[0] ])       histos[lepsn1+"SR_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[ iSS[0] ]) histos[lepsn1+"SR_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[ iSS[0] ]) histos[lepsn1+"SR_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[ iSS[0] ]) histos[lepsn1+"SR_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[ iSS[0] ]) histos[lepsn1+"SR_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[iSS[0] ]) histos[lepsn1+"SR_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iSS[0] ]==(-3)) histos[lepsn1+"SR_isFromFlags_"+sn]->Fill( 7.,weight);
	else                             histos[lepsn1+"SR_isFromFlags_"+sn]->Fill(-1.,weight);
	if(lep_isFromW()[iSS[1] ])       histos[lepsn2+"SR_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[ iSS[1] ]) histos[lepsn2+"SR_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[ iSS[1] ]) histos[lepsn2+"SR_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[ iSS[1] ]) histos[lepsn2+"SR_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[ iSS[1] ]) histos[lepsn2+"SR_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[iSS[1] ]) histos[lepsn2+"SR_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iSS[1] ]==(-3)) histos[lepsn2+"SR_isFromFlags_"+sn]->Fill( 7.,weight);
	else                             histos[lepsn2+"SR_isFromFlags_"+sn]->Fill(-1.,weight);
      }
      if((ee[1]||em[1]||mm[1])&&!isData()){
	histos[lepsn1+"SRpresel_motherIDSS_"+sn]->Fill(lep_motherIdSS()[iSS[0] ],weight);
	histos[lepsn2+"SRpresel_motherIDSS_"+sn]->Fill(lep_motherIdSS()[iSS[1] ],weight);
	if(lep_isFromW()[iSS[0] ])       histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[ iSS[0] ]) histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[ iSS[0] ]) histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[ iSS[0] ]) histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[ iSS[0] ]) histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[iSS[0] ]) histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iSS[0] ]==(-3)) histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill( 7.,weight);
	else                             histos[lepsn1+"SRpresel_isFromFlags_"+sn]->Fill(-1.,weight);
	if(lep_isFromW()[iSS[1] ])       histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[ iSS[1] ]) histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[ iSS[1] ]) histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[ iSS[1] ]) histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[ iSS[1] ]) histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[iSS[1] ]) histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iSS[1] ]==(-3)) histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill( 7.,weight);
	else                             histos[lepsn2+"SRpresel_isFromFlags_"+sn]->Fill(-1.,weight);
      }

      if(ee[2]||em[2]||mm[2]){
	histos[lepsn1+"AR_motherIDSS_"+sn]->Fill(lep_motherIdSS()[ iSS[0] ],weight);
	histos[lepsn2+"AR_motherIDSS_"+sn]->Fill(lep_motherIdSS()[iaSS[0] ],weight);
	if(lep_isFromW()[ iSS[0] ])       histos[lepsn1+"AR_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[  iSS[0] ]) histos[lepsn1+"AR_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[  iSS[0] ]) histos[lepsn1+"AR_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[  iSS[0] ]) histos[lepsn1+"AR_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[  iSS[0] ]) histos[lepsn1+"AR_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[ iSS[0] ]) histos[lepsn1+"AR_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iSS[0] ]==(-3)) histos[lepsn1+"AR_isFromFlags_"+sn]->Fill( 7.,weight);
	else                              histos[lepsn1+"AR_isFromFlags_"+sn]->Fill(-1.,weight);
	if(lep_isFromW()[iaSS[0] ])       histos[lepsn2+"AR_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[ iaSS[0] ]) histos[lepsn2+"AR_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[ iaSS[0] ]) histos[lepsn2+"AR_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[ iaSS[0] ]) histos[lepsn2+"AR_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[ iaSS[0] ]) histos[lepsn2+"AR_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[iaSS[0] ]) histos[lepsn2+"AR_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iaSS[0] ]==(-3)) histos[lepsn2+"AR_isFromFlags_"+sn]->Fill( 7.,weight);
	else                              histos[lepsn2+"AR_isFromFlags_"+sn]->Fill(-1.,weight);
	histos[lepsn1+"AR_tight_pT_"    +sn]->Fill(lep_p4()[iSS[0] ].Pt(),       weight);
	histos[lepsn1+"AR_tight_eta_"   +sn]->Fill(fabs(lep_p4()[iSS[0] ].Eta()),weight);
	histos[lepsn1+"AR_tight_phi_"   +sn]->Fill(fabs(lep_p4()[iSS[0] ].Phi()),weight);
	histos[lepsn1+"AR_tight_reliso_"+sn]->Fill(lep_relIso03EA()[iSS[0] ],    weight);
	histos[lepsn2+"AR_loose_pT_"    +sn]->Fill(lep_p4()[iaSS[0] ].Pt(),       weight);
	histos[lepsn2+"AR_loose_eta_"   +sn]->Fill(fabs(lep_p4()[iaSS[0] ].Eta()),weight);
	histos[lepsn2+"AR_loose_phi_"   +sn]->Fill(fabs(lep_p4()[iaSS[0] ].Phi()),weight);
	histos[lepsn2+"AR_loose_reliso_"+sn]->Fill(lep_relIso03EA()[iaSS[0] ],    weight);
      }
      if(ee[3]||em[3]||mm[3]){
	histos[lepsn1+"ARpresel_motherIDSS_"+sn]->Fill(lep_motherIdSS()[ iSS[0] ],weight);
	histos[lepsn2+"ARpresel_motherIDSS_"+sn]->Fill(lep_motherIdSS()[iaSS[0] ],weight);
	if(lep_isFromW()[ iSS[0] ])       histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[  iSS[0] ]) histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[  iSS[0] ]) histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[  iSS[0] ]) histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[  iSS[0] ]) histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[ iSS[0] ]) histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iSS[0] ]==(-3)) histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill( 7.,weight);
	else                              histos[lepsn1+"ARpresel_isFromFlags_"+sn]->Fill(-1.,weight);
	if(lep_isFromW()[iaSS[0] ])       histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill( 1.,weight);
	else if(lep_isFromZ()[ iaSS[0] ]) histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill( 2.,weight);
	else if(lep_isFromB()[ iaSS[0] ]) histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill( 3.,weight);
	else if(lep_isFromC()[ iaSS[0] ]) histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill( 4.,weight);
	else if(lep_isFromL()[ iaSS[0] ]) histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill( 5.,weight);
	else if(lep_isFromLF()[iaSS[0] ]) histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill( 6.,weight);
	else if(lep_motherIdSS()[iaSS[0] ]==(-3)) histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill( 7.,weight);
	else                              histos[lepsn2+"ARpresel_isFromFlags_"+sn]->Fill(-1.,weight);
	histos[lepsn1+"ARpresel_tight_pT_"    +sn]->Fill(lep_p4()[iSS[0] ].Pt(),       weight);
	histos[lepsn1+"ARpresel_tight_eta_"   +sn]->Fill(fabs(lep_p4()[iSS[0] ].Eta()),weight);
	histos[lepsn1+"ARpresel_tight_phi_"   +sn]->Fill(fabs(lep_p4()[iSS[0] ].Phi()),weight);
	histos[lepsn1+"ARpresel_tight_reliso_"+sn]->Fill(lep_relIso03EA()[iSS[0] ],    weight);
	histos[lepsn2+"ARpresel_loose_pT_"    +sn]->Fill(lep_p4()[iaSS[0] ].Pt(),       weight);
	histos[lepsn2+"ARpresel_loose_eta_"   +sn]->Fill(fabs(lep_p4()[iaSS[0] ].Eta()),weight);
	histos[lepsn2+"ARpresel_loose_phi_"   +sn]->Fill(fabs(lep_p4()[iaSS[0] ].Phi()),weight);
	histos[lepsn2+"ARpresel_loose_reliso_"+sn]->Fill(lep_relIso03EA()[iaSS[0] ],    weight);
      }
      //cout << __LINE__ << endl;
	
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
  
  string filename = "rootfiles/CheckWGrelatedstuff.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //string mapname = histonames[i]+"_"+skimFilePrefix;
    //string writename = histonames[i];
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }
  if(skimFilePrefix.find("Background")!=string::npos) h2D->Write(h2D->GetName(),TObject::kOverwrite);
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;

  /*
  ofstream eventstocheckEElog    ((string("logs/eventListEE_"   +skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheckEMlog    ((string("logs/eventListEM_"   +skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheckMMlog    ((string("logs/eventListMM_"   +skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheck0SFOSlog ((string("logs/eventList0SFOS_"+skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheck1SFOSlog ((string("logs/eventList1SFOS_"+skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheck2SFOSlog ((string("logs/eventList2SFOS_"+skimFilePrefix+".log")).c_str(), ios::trunc);
  eventstocheckEElog    << eventstocheckEE   ->str();
  eventstocheckEMlog    << eventstocheckEM   ->str();
  eventstocheckMMlog    << eventstocheckMM   ->str();
  eventstocheck0SFOSlog << eventstocheck0SFOS->str();
  eventstocheck1SFOSlog << eventstocheck1SFOS->str();
  eventstocheck2SFOSlog << eventstocheck2SFOS->str();
  */
  cout << "ee " << cCRee << " " << cCR2ee << " " << cCR3ee << endl;
  cout << "em " << cCRem << " " << cCR2em << " " << cCR3em << endl;
  cout << "mm " << cCRmm << " " << cCR2mm << " " << cCR3mm << endl;
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
  return 0;
}
