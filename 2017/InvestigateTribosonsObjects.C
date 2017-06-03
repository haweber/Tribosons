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
#include "CMS3.cc"
#include "/home/users/haweber/CORE/Tools/MT2/MT2Utility.h"
#include "/home/users/haweber/CORE/Tools/MT2/MT2Utility.cc"
//#include "/home/users/haweber/CORE/Tools/dorky/dorky.h"
//#include "/home/users/haweber/CORE/Tools/dorky/dorky.cc"
//#include "/home/users/haweber/CORE/Tools/goodrun.h"
//#include "/home/users/haweber/CORE/Tools/goodrun.cc"

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

  vector<string> histonames2; histonames2.clear();
  vector<int> hbins2; hbins2.clear();
  vector<float> hlow2; hlow2.clear();
  vector<float> hup2; hup2.clear();
  histonames.push_back("e_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("e_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightID_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightID_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightID_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightID_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightID_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("e_tightID_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightIP_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightIP_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightIP_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightIP_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightIP_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("e_tightIP_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightID_tightIP_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightID_tightIP_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightID_tightIP_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightID_tightIP_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightID_tightIP_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("e_tightID_tightIP_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightID_tightIP_tightreliso_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightID_tightIP_tightreliso_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightID_tightIP_tightreliso_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("e_tightID_tightIP_tightreliso_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_tightreliso_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightreliso_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("e_tightreliso_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("e_tightreliso_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("e_ID");                     hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_tightIP_tightreliso_ID"); hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_tightreliso_ID");         hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_tightIP_ID");             hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_IDv2");                     hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_tightIP_tightreliso_IDv2"); hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_tightreliso_IDv2");         hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_tightIP_IDv2");             hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("e_tq");                             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_tightIP_tightreliso_tq"); hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_tightIP_tq");             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_tightreliso_tq");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_tq");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightIP_tightreliso_tq");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightreliso_tq");                 hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightIP_tq");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_3q");                             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_tightIP_tightreliso_3q"); hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_tightIP_3q");             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_tightreliso_3q");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightID_3q");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightIP_tightreliso_3q");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightreliso_3q");                 hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_tightIP_3q");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("e_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("e_tightID_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_tightID_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_tightID_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("e_tightID_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_tightID_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_tightID_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("e_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("e_tightreliso_tightID_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_tightreliso_tightID_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightreliso_tightID_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_tightreliso_tightID_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightreliso_tightID_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("e_tightID_tightIP_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_tightID_tightIP_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_tightIP_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_tightID_tightIP_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_tightIP_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("e_tightIP_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_tightIP_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightIP_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_tightIP_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightIP_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("e_tightID_tightIP_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("e_tightID_tightIP_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_tightIP_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("e_tightID_tightIP_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("e_tightID_tightIP_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);

  histonames.push_back("mu_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("mu_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightID_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightID_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightID_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightID_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightID_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("mu_tightID_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightIP_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightIP_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightIP_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightIP_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightIP_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("mu_tightIP_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightID_tightIP_reliso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightID_tightIP_reliso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightID_tightIP_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightID_tightIP_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightID_tightIP_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("mu_tightID_tightIP_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightID_tightIP_tightreliso_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightID_tightIP_tightreliso_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightID_tightIP_tightreliso_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("mu_tightID_tightIP_tightreliso_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_tightreliso_absiso03");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightreliso_absiso04");       hbins.push_back(30); hlow.push_back(0); hup.push_back(15);
  histonames.push_back("mu_tightreliso_ptrel");          hbins.push_back(20); hlow.push_back(0); hup.push_back(10);
  histonames.push_back("mu_tightreliso_ptratio");        hbins.push_back(30); hlow.push_back(0); hup.push_back(1.5);
  histonames.push_back("mu_ID");                     hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_tightIP_tightreliso_ID"); hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_tightreliso_ID");         hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_tightIP_ID");             hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_IDv2");                     hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_tightIP_tightreliso_IDv2"); hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_tightreliso_IDv2");         hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_tightIP_IDv2");             hbins.push_back(4); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("mu_tq");                             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_tightIP_tightreliso_tq"); hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_tightIP_tq");             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_tightreliso_tq");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_tq");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightIP_tightreliso_tq");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightreliso_tq");                 hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightIP_tq");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_3q");                             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_tightIP_tightreliso_3q"); hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_tightIP_3q");             hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_tightreliso_3q");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightID_3q");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightIP_tightreliso_3q");         hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightreliso_3q");                 hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_tightIP_3q");                     hbins.push_back(4); hlow.push_back(-1); hup.push_back(3);
  histonames.push_back("mu_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("mu_tightID_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_tightID_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_tightID_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("mu_tightID_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_tightID_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_tightID_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("mu_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("mu_tightreliso_tightID_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_tightreliso_tightID_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightreliso_tightID_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_tightreliso_tightID_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightreliso_tightID_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("mu_tightID_tightIP_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_tightID_tightIP_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_tightIP_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_tightID_tightIP_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_tightIP_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("mu_tightIP_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_tightIP_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightIP_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_tightIP_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightIP_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);
  histonames.push_back("mu_tightID_tightIP_tightreliso_dxy");            hbins.push_back(20); hlow.push_back(-0.05); hup.push_back(0.05);
  histonames.push_back("mu_tightID_tightIP_tightreliso_dz");             hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_tightIP_tightreliso_dzOverdzerr");    hbins.push_back(30); hlow.push_back(-15);   hup.push_back(15);
  histonames.push_back("mu_tightID_tightIP_tightreliso_IP3");            hbins.push_back(20); hlow.push_back(-0.1);  hup.push_back(0.1);
  histonames.push_back("mu_tightID_tightIP_tightreliso_IP3overIP3err");  hbins.push_back(20); hlow.push_back(-10);   hup.push_back(10);

  histonames2.push_back("itrl_reliso03");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);
  histonames2.push_back("itrl_reliso04");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);
  histonames2.push_back("itrl_absiso03");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrl_absiso04");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrl_ptrel");          hbins2.push_back(20); hlow2.push_back(0); hup2.push_back(10);
  histonames2.push_back("itrl_ptratio");        hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);
  histonames2.push_back("itrh_reliso03");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);
  histonames2.push_back("itrh_reliso04");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);
  histonames2.push_back("itrh_absiso03");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrh_absiso04");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrh_ptrel");          hbins2.push_back(20); hlow2.push_back(0); hup2.push_back(10);
  histonames2.push_back("itrh_ptratio");        hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);
  histonames2.push_back("itrl_tightreliso_absiso03");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrl_tightreliso_absiso04");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrl_tightreliso_ptrel");          hbins2.push_back(20); hlow2.push_back(0); hup2.push_back(10);
  histonames2.push_back("itrl_tightreliso_ptratio");        hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);
  histonames2.push_back("itrh_tightreliso_absiso03");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrh_tightreliso_absiso04");       hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(15);
  histonames2.push_back("itrh_tightreliso_ptrel");          hbins2.push_back(20); hlow2.push_back(0); hup2.push_back(10);
  histonames2.push_back("itrh_tightreliso_ptratio");        hbins2.push_back(30); hlow2.push_back(0); hup2.push_back(1.5);

  cout << "booking histograms" << endl;
  string prefix[3] = {"","SS","l3"};
  string postfix[4] = {"_fromWZ","_fromCB","_fromL","_fromOth"};
  for(unsigned int i = 0; i<histonames.size(); ++i){
    for(unsigned int j = 0; j<3; ++j){
      for(unsigned int k = 0; k<4; ++k){
	string postfixx = postfix[k];
	if( histonames[i].find("itrl_")!=string::npos || histonames[i].find("itrh_")!=string::npos ){
	  if(k!=0) continue;
	  postfixx = "";
	}
	string mapname = prefix[j] + histonames[i] + "_"+skimFilePrefix + postfixx;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
	//cout << mapname << endl;
      }
    }
  }
  string prefix2[3] = {"","SS","l3"};
  string postfix2[4] = {"_matchedLepWZ","_matchedLep","_matchedTau","_unmatched"};//DR=0.1 for e/mu, 0.2 for tau
  for(unsigned int i = 0; i<histonames2.size(); ++i){
    for(unsigned int j = 0; j<3; ++j){
      for(unsigned int k = 0; k<4; ++k){
	string postfixx = postfix2[k];
	string mapname = prefix2[j] + histonames2[i] + "_"+skimFilePrefix + postfix2[k];
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
	//cout << i << " " << mapname << endl;
      }
    }
  }
  cout << "done" << endl;

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
  int WZcount(0), Othcount(0);

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


      string sn = skimFilePrefix;
      
      double weight = evt_scale1fb()*70.;

      if(nVert()<0)               continue;
      if(nlep()<2)               continue;

      bool SS = false;
      if(nlep()==2&&((lep_pdgId()[0]*lep_pdgId()[1])>0)) SS = true;
      if(false){
	for(unsigned int n = 0; n<lep_pdgId().size();++n){
	  cout << "lep ID " << lep_pdgId()[n] << " px " << lep_p4()[n].Px() << " py " << lep_p4()[n].Py() << " py " << lep_p4()[n].Py() << " fromW " << lep_isfromW()[n] << " Z " << lep_isfromZ()[n] << " B " << lep_isfromB()[n] << " C " << lep_isfromC()[n] << " L " << lep_isfromL()[n] << " mcID " << lep_mcMatchId()[n] << endl;
	}
	for(unsigned int ng = 0; ng<genPart_pdgId().size();++ng){
	  if(abs(genPart_pdgId()[ng])!=11&&abs(genPart_pdgId()[ng])!=13) continue;
	  cout << "gen ID " << genPart_pdgId()[ng] << " px " << genPart_p4()[ng].Px() << " py " << genPart_p4()[ng].Py() << " pz " << genPart_p4()[ng].Pz() << " mom " << genPart_motherId()[ng] << " grandma " << genPart_grandmaId()[ng] << endl;
	}
      }

      bool l3Z = false;
      bool l3 = false;
      if(nlep()>=3){
	l3 = true;
	if((lep_pdgId()[0]==(-lep_pdgId()[1]))&&fabs((lep_p4()[0]+lep_p4()[1]).M()-91.)<15) { l3Z = true; }
	if((lep_pdgId()[0]==(-lep_pdgId()[2]))&&fabs((lep_p4()[0]+lep_p4()[2]).M()-91.)<15) { l3Z = true; }
	if((lep_pdgId()[3]==(-lep_pdgId()[1]))&&fabs((lep_p4()[2]+lep_p4()[1]).M()-91.)<15) { l3Z = true; }
      }
      //  string prefix[3] = {"","SS","l3"};
      //  string postfix[4] = {"_fromWZ","_fromCB","_fromL","_fromOth"};
      for(unsigned int n = 0; n<lep_pdgId().size();++n){
	if(lep_p4()[n].Pt()<10) continue;
	if(fabs(lep_p4()[n].Eta())>2.4) continue;
	string id = "";
	if(abs(lep_pdgId()[n])==11) id = "e";
	else if(abs(lep_pdgId()[n])==13) id = "mu";
	else {cout << "WTF lep id " << lep_pdgId()[n] << endl; continue; }
	bool tightID = lep_tightId()[n];
	bool tightreliso = lep_relIso04()[n]<0.1;
	bool tightIPdz = fabs(lep_dz()[n])<0.05;
	bool tightIPdxy = fabs(lep_dxy()[n])<0.02;
	for(int pr = 0; pr<3; ++pr){
	  if(pr==1&&!SS) continue;
	  if(pr==2&&!l3Z) continue;
	  bool mycontinue = true;
	  for(int po = 0; po<4; ++po){
	    if(po==0) {
	      if(lep_isfromW()[n]) mycontinue = false;
	      if(lep_isfromZ()[n]) mycontinue = false;
	      //if(mycontinue==false&&(lep_isfromB()[n]||lep_isfromC()[n]||lep_isfromL()[n])) cout << "So there is double counting: " <<lep_isfromB()[n]<<" "<<lep_isfromC()[n]<<" "<<lep_isfromL()[n]<<" "<<lep_isfromW()[n]<<" "<<lep_isfromZ()[n]<<endl;
	      //if(!(lep_isfromW()[n]||lep_isfromZ()[n])) continue;
	    } else if(po==1){
	      if(lep_isfromB()[n]) mycontinue = false;
	      if(lep_isfromC()[n]) mycontinue = false;
	      if(lep_isfromW()[n]) mycontinue = true;//avoid double counting
	      if(lep_isfromZ()[n]) mycontinue = true;//avoid double counting
	      //if(!(lep_isfromB()[n]||lep_isfromC()[n])) continue;
	      //if((lep_isfromW()[n]||lep_isfromZ()[n])) {cout<<"WTF "<<__LINE__<<" "<<lep_isfromB()[n]<<" "<<lep_isfromC()[n]<<" "<<lep_isfromL()[n]<<" "<<lep_isfromW()[n]<<" "<<lep_isfromZ()[n]<<endl;continue;}
	    } else if(po==2){
	      if(lep_isfromL()[n]) mycontinue = false;
	      if(lep_isfromB()[n]) mycontinue = true;//avoid double counting
	      if(lep_isfromC()[n]) mycontinue = true;//avoid double counting
	      if(lep_isfromW()[n]) mycontinue = true;//avoid double counting
	      if(lep_isfromZ()[n]) mycontinue = true;//avoid double counting
	      //if(!(lep_isfromL()[n])) continue;
	      //if((lep_isfromB()[n]||lep_isfromC()[n])) {cout<<"WTF "<<__LINE__<<" "<<lep_isfromB()[n]<<" "<<lep_isfromC()[n]<<" "<<lep_isfromL()[n]<<" "<<lep_isfromW()[n]<<" "<<lep_isfromZ()[n]<<endl;continue;}
	      //if((lep_isfromW()[n]||lep_isfromZ()[n])) {cout<<"WTF "<<__LINE__<<" "<<lep_isfromB()[n]<<" "<<lep_isfromC()[n]<<" "<<lep_isfromL()[n]<<" "<<lep_isfromW()[n]<<" "<<lep_isfromZ()[n]<<endl;continue;}
	      
	    } else if(po==3){
	      mycontinue = false;
	      if(lep_isfromW()[n]) mycontinue = true;
	      if(lep_isfromZ()[n]) mycontinue = true;
	      if(lep_isfromB()[n]) mycontinue = true;
	      if(lep_isfromC()[n]) mycontinue = true;
	      if(lep_isfromL()[n]) mycontinue = true;
	    }
	    if(mycontinue) continue;
	    string mypr = prefix[pr]+id +"_";
	    string mypo = "_" + skimFilePrefix + postfix[po];
	    if(pr==1&&po==0){ ++WZcount; /*cout << n << " This is a WZ lepton, right " << lep_isfromW()[n] << " " << lep_isfromZ()[n] << " and it is lepID " << lep_pdgId()[n] << " " << id << " " << lep_p4()[n].Pt()  << " " << mypr+mypo << endl;*/ }
	    else if (pr==1) { ++Othcount; /*cout << n << " This is NOT a WZ lepton, right " << lep_isfromW()[n] << " " << lep_isfromZ()[n] << " and it is lepID " << lep_pdgId()[n] << " " << id << " " << lep_p4()[n].Pt() << " " << po << " " << mypr+mypo << endl;*/ }
	    int partid = 0;
	    if(mypo.find("WZ")!=string::npos&&((lep_isfromB()[n]||lep_isfromC()[n]||lep_isfromL()[n])||!(lep_isfromW()[n]||lep_isfromZ()[n]))){
		cout<<"WTF "<<__LINE__<<" "<<mypo<<" " << po << " " <<lep_isfromB()[n]<<" "<<lep_isfromC()[n]<<" "<<lep_isfromL()[n]<<" "<<lep_isfromW()[n]<<" "<<lep_isfromZ()[n]<<endl;continue;
	      }
	    if(lep_looseId()[n])  partid = 1;
	    if(lep_mediumId()[n]) partid = 2;
	    if(lep_tightId()[n])  partid = 3;
	    //cout << __LINE__ << " " << mypr+"reliso03"+mypo << endl;
	    histos[mypr+"reliso03"+mypo]->Fill(lep_relIso03()[n],weight);
	    histos[mypr+"reliso04"+mypo]->Fill(lep_relIso04()[n],weight);
	    histos[mypr+"absiso03"+mypo]->Fill(lep_relIso03()[n]*lep_p4()[n].Pt(),weight);
	    histos[mypr+"absiso04"+mypo]->Fill(lep_relIso04()[n]*lep_p4()[n].Pt(),weight);
	    histos[mypr+"ptrel"   +mypo]->Fill(lep_ptrel()[n],weight);
	    histos[mypr+"ptratio" +mypo]->Fill(lep_ptratio()[n],weight);
	    histos[mypr+"ID" +mypo]->Fill(partid,weight);
	    histos[mypr+"IDv2" +mypo]->Fill(0.,weight);
	    if(lep_looseId()[n])  histos[mypr+"IDv2" +mypo]->Fill(1.,weight);
	    if(lep_mediumId()[n]) histos[mypr+"IDv2" +mypo]->Fill(2.,weight);
	    if(lep_tightId()[n])  histos[mypr+"IDv2" +mypo]->Fill(3.,weight);
	    histos[mypr+"dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	    histos[mypr+"dz"           +mypo]->Fill(lep_dz()[n],weight);
	    histos[mypr+"dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	    histos[mypr+"IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	    histos[mypr+"IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    histos[mypr+"3q" +mypo]->Fill(lep_threecharge()[n],weight);
	    histos[mypr+"tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	    //cout << __LINE__ << endl;
	    if(tightID){
	      histos[mypr+"tightID_reliso03"+mypo]->Fill(lep_relIso03()[n],weight);
	      histos[mypr+"tightID_reliso04"+mypo]->Fill(lep_relIso04()[n],weight);
	      histos[mypr+"tightID_absiso03"+mypo]->Fill(lep_relIso03()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightID_absiso04"+mypo]->Fill(lep_relIso04()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightID_ptrel"   +mypo]->Fill(lep_ptrel()[n],weight);
	      histos[mypr+"tightID_ptratio" +mypo]->Fill(lep_ptratio()[n],weight);
	      histos[mypr+"tightID_dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	      histos[mypr+"tightID_dz"           +mypo]->Fill(lep_dz()[n],weight);
	      histos[mypr+"tightID_dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	      histos[mypr+"tightID_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	      histos[mypr+"tightID_3q" +mypo]->Fill(lep_threecharge()[n],weight);
	      histos[mypr+"tightID_tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	    }
	    //cout << __LINE__ << endl;
	    if(tightIPdz&&tightIPdxy){
	      histos[mypr+"tightIP_reliso03"+mypo]->Fill(lep_relIso03()[n],weight);
	      histos[mypr+"tightIP_reliso04"+mypo]->Fill(lep_relIso04()[n],weight);
	      histos[mypr+"tightIP_absiso03"+mypo]->Fill(lep_relIso03()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightIP_absiso04"+mypo]->Fill(lep_relIso04()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightIP_ptrel"   +mypo]->Fill(lep_ptrel()[n],weight);
	      histos[mypr+"tightIP_ptratio" +mypo]->Fill(lep_ptratio()[n],weight);
	      histos[mypr+"tightIP_ID"     +mypo]->Fill(partid,weight);
	      histos[mypr+"tightIP_3q" +mypo]->Fill(lep_threecharge()[n],weight);
	      histos[mypr+"tightIP_tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	      histos[mypr+"tightIP_IDv2" +mypo]->Fill(0.,weight);
	      if(lep_looseId()[n])  histos[mypr+"tightIP_IDv2" +mypo]->Fill(1.,weight);
	      if(lep_mediumId()[n]) histos[mypr+"tightIP_IDv2" +mypo]->Fill(2.,weight);
	      if(lep_tightId()[n])  histos[mypr+"tightIP_IDv2" +mypo]->Fill(3.,weight);
	    }
	    //cout << __LINE__ << endl;
	    if(tightreliso){
	      histos[mypr+"tightreliso_ID"           +mypo]->Fill(partid,weight);
	      histos[mypr+"tightreliso_dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	      histos[mypr+"tightreliso_dz"           +mypo]->Fill(lep_dz()[n],weight);
	      histos[mypr+"tightreliso_dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	      histos[mypr+"tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	      histos[mypr+"tightreliso_3q" +mypo]->Fill(lep_threecharge()[n],weight);
	      histos[mypr+"tightreliso_tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	      histos[mypr+"tightreliso_absiso03"+mypo]->Fill(lep_relIso03()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightreliso_absiso04"+mypo]->Fill(lep_relIso04()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightreliso_ptrel"   +mypo]->Fill(lep_ptrel()[n],weight);
	      histos[mypr+"tightreliso_ptratio" +mypo]->Fill(lep_ptratio()[n],weight);
	      histos[mypr+"tightreliso_IDv2" +mypo]->Fill(0.,weight);
	      if(lep_looseId()[n])  histos[mypr+"tightreliso_IDv2" +mypo]->Fill(1.,weight);
	      if(lep_mediumId()[n]) histos[mypr+"tightreliso_IDv2" +mypo]->Fill(2.,weight);
	      if(lep_tightId()[n])  histos[mypr+"tightreliso_IDv2" +mypo]->Fill(3.,weight);
	    }
	    //cout << __LINE__ << endl;
	    if(tightID&&tightIPdz&&tightIPdxy){
	      histos[mypr+"tightID_tightIP_reliso03"+mypo]->Fill(lep_relIso03()[n],weight);
	      histos[mypr+"tightID_tightIP_reliso04"+mypo]->Fill(lep_relIso04()[n],weight);
	      histos[mypr+"tightID_tightIP_absiso03"+mypo]->Fill(lep_relIso03()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightID_tightIP_absiso04"+mypo]->Fill(lep_relIso04()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightID_tightIP_ptrel"   +mypo]->Fill(lep_ptrel()[n],weight);
	      histos[mypr+"tightID_tightIP_ptratio" +mypo]->Fill(lep_ptratio()[n],weight);
	      histos[mypr+"tightID_tightIP_3q" +mypo]->Fill(lep_threecharge()[n],weight);
	      histos[mypr+"tightID_tightIP_tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	    }
	    //cout << __LINE__ << endl;
	    if(tightID&&tightreliso){
	      histos[mypr+"tightID_tightreliso_dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	      histos[mypr+"tightID_tightreliso_dz"           +mypo]->Fill(lep_dz()[n],weight);
	      histos[mypr+"tightID_tightreliso_dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	      histos[mypr+"tightID_tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	      histos[mypr+"tightID_tightreliso_3q" +mypo]->Fill(lep_threecharge()[n],weight);
	      histos[mypr+"tightID_tightreliso_tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	    }
	    //cout << __LINE__ << endl;
	    if(tightIPdz&&tightIPdxy&&tightreliso){
	      histos[mypr+"tightIP_tightreliso_ID" +mypo]->Fill(partid,weight);
	      histos[mypr+"tightIP_tightreliso_3q" +mypo]->Fill(lep_threecharge()[n],weight);
	      histos[mypr+"tightIP_tightreliso_tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	      histos[mypr+"tightIP_tightreliso_IDv2" +mypo]->Fill(0.,weight);
	      if(lep_looseId()[n])  histos[mypr+"tightIP_tightreliso_IDv2" +mypo]->Fill(1.,weight);
	      if(lep_mediumId()[n]) histos[mypr+"tightIP_tightreliso_IDv2" +mypo]->Fill(2.,weight);
	      if(lep_tightId()[n])  histos[mypr+"tightIP_tightreliso_IDv2" +mypo]->Fill(3.,weight);
	    }
	    if(tightID&&tightIPdz&&tightIPdxy&&tightreliso){
	      histos[mypr+"tightID_tightIP_tightreliso_3q" +mypo]->Fill(lep_threecharge()[n],weight);
	      histos[mypr+"tightID_tightIP_tightreliso_tq" +mypo]->Fill(lep_tightcharge()[n],weight);
	      //cout << mypr+"tightID_tightIP_tightreliso_absiso03"+mypo << endl;
	      histos[mypr+"tightID_tightIP_tightreliso_absiso03"+mypo]->Fill(lep_relIso03()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightID_tightIP_tightreliso_absiso04"+mypo]->Fill(lep_relIso04()[n]*lep_p4()[n].Pt(),weight);
	      histos[mypr+"tightID_tightIP_tightreliso_ptrel"   +mypo]->Fill(lep_ptrel()[n],weight);
	      histos[mypr+"tightID_tightIP_tightreliso_ptratio" +mypo]->Fill(lep_ptratio()[n],weight);
	    }
	    //special for dz
	    //cout << __LINE__ << endl;
	    if(tightID&&tightIPdxy){
	      histos[mypr+"tightID_tightIP_dz"           +mypo]->Fill(lep_dz()[n],weight);
	      histos[mypr+"tightID_tightIP_dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	      histos[mypr+"tightID_tightIP_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_tightIP_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    if(tightIPdxy&&tightreliso){
	      histos[mypr+"tightIP_tightreliso_dz"           +mypo]->Fill(lep_dz()[n],weight);
	      histos[mypr+"tightIP_tightreliso_dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	      histos[mypr+"tightIP_tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightIP_tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    if(tightID&&tightreliso){
	      histos[mypr+"tightID_tightreliso_dz"           +mypo]->Fill(lep_dz()[n],weight);
	      histos[mypr+"tightID_tightreliso_dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	      histos[mypr+"tightID_tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    if(tightID&&tightIPdxy&&tightreliso){
	      histos[mypr+"tightID_tightIP_tightreliso_dz"           +mypo]->Fill(lep_dz()[n],weight);
	      histos[mypr+"tightID_tightIP_tightreliso_dzOverdzerr"  +mypo]->Fill(lep_dz()[n]/lep_dzerr()[n],weight);
	      histos[mypr+"tightID_tightIP_tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_tightIP_tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    //special for dxy
	    //cout << __LINE__ << endl;
	    if(tightID&&tightIPdz){
	      histos[mypr+"tightID_tightIP_dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	      histos[mypr+"tightID_tightIP_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_tightIP_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    if(tightIPdz&&tightreliso){
	      histos[mypr+"tightIP_tightreliso_dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	      histos[mypr+"tightIP_tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightIP_tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    if(tightID&&tightreliso){
	      histos[mypr+"tightID_tightreliso_dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	      histos[mypr+"tightID_tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    if(tightID&&tightIPdz&&tightreliso){
	      histos[mypr+"tightID_tightIP_tightreliso_dxy"          +mypo]->Fill(lep_dxy()[n],weight);
	      histos[mypr+"tightID_tightIP_tightreliso_IP3"          +mypo]->Fill(lep_ip3d()[n],weight);
	      histos[mypr+"tightID_tightIP_tightreliso_IP3overIP3err"+mypo]->Fill(lep_ip3d()[n]/lep_ip3derr()[n],weight);
	    }
	    //cout << __LINE__ << endl;
	  }
	}
      }

      vector<int> gentauid; gentauid.clear();
      for(unsigned int ng = 0; ng<genPart_pdgId().size();++ng){
	if(abs(genPart_pdgId()[ng])!=11&&abs(genPart_pdgId()[ng])!=13) continue;
	if(abs(genPart_motherId()[ng])==15&&abs(genPart_grandmaId()[ng])==24) gentauid.push_back(genPart_motherId()[ng]);
	if(abs(genPart_motherId()[ng])==15&&abs(genPart_grandmaId()[ng])==23) gentauid.push_back(genPart_motherId()[ng]);
	//if(abs(genPart_motherId()[ng])==15&&abs(genPart_grandmaId()[ng])==22) gentauid.push_back(genPart_motherId()[ng]);
      }
      for(unsigned int n = 0; n<isotr_pdgId().size();++n){
	if(isotr_p4()[n].Pt()<5) continue;
	if(fabs(isotr_p4()[n].Eta())>2.4) continue;
	bool idlep = false;
	if(abs(isotr_pdgId()[n])==11) idlep = true;
	else if(abs(isotr_pdgId()[n])==13) idlep = true;
	if(!idlep){
	if(isotr_p4()[n].Pt()<10) continue;
	}
	bool matchedrecolep = false;
	for(unsigned int nl = 0; nl<lep_pdgId().size();++nl){
	  if(lep_p4()[nl].Pt()<10) continue;
	  if(fabs(lep_p4()[nl].Eta())>2.4) continue;
	  if(dR(lep_p4()[nl],isotr_p4()[n])>0.1) continue;
	  matchedrecolep = true;
	  break;
	}
	if(matchedrecolep) continue;
	bool matchgenlepWZ = false;
	bool matchgenlep = false;
	bool matchgentau = false;//damn - cannot test if this is hadronic - drop any potential lep-tau candidate
	for(unsigned int ng = 0; ng<genPart_pdgId().size();++ng){
	  if(abs(genPart_pdgId()[ng])!=11&&abs(genPart_pdgId()[ng])!=13&&abs(genPart_pdgId()[ng])!=15) continue;
	  if(abs(genPart_pdgId()[ng])==11||abs(genPart_pdgId()[ng])==13){
	    if(abs(genPart_motherId()[ng])<=24&&abs(genPart_motherId()[ng])>=23&&dR(genPart_p4()[ng],isotr_p4()[n])<0.15) {
	      matchgenlepWZ = true;
	      break;
	    }
	    if(abs(genPart_motherId()[ng])==15&&abs(genPart_grandmaId()[ng])<=24&&abs(genPart_grandmaId()[ng])>=23&&dR(genPart_p4()[ng],isotr_p4()[n])<0.15) {
	      matchgenlepWZ = true;
	      break;
	    }
	  }
	  //test if it fits a tau only if it does not fit a lepton
	  if(abs(genPart_pdgId()[ng])==15){
	    bool dropthis = false;
	    for(unsigned int nt = 0; nt<gentauid.size(); ++nt){
	      if(gentauid[nt]==genPart_pdgId()[ng]) { dropthis = true; break; }
	    }
	    if(dropthis) continue;
	    if(abs(genPart_motherId()[ng])<=24&&abs(genPart_motherId()[ng])>=23&&dR(genPart_p4()[ng],isotr_p4()[n])<0.3) {
	      matchgentau = true;
	      break;
	    }
	  }
	  if((abs(genPart_pdgId()[ng])==11||abs(genPart_pdgId()[ng])==13)&&dR(genPart_p4()[ng],isotr_p4()[n])<0.15) {
	    matchgenlep = true;
	    break;
	  }
	}
	    
	for(int pr = 0; pr<3; ++pr){
	  if(pr==1&&!SS) continue;
	    if(pr==2&&!l3Z) continue;
	    string id = "itrh";
	    if(idlep) id = "itrl";
	    string mypr = prefix[pr]+id +"_";
	    for(int po = 0; po<4; ++po){
	      string mypo = "_" + skimFilePrefix + postfix2[po];
	      if(po==0&&!matchgenlepWZ) continue;
	      if(po==1&&!matchgenlep)   continue;
	      if(po==1&&matchgenlepWZ)  continue;
	      if(po==2&&!matchgentau)   continue;
	      if(po==2&&matchgenlepWZ)  continue;
	      if(po==2&&matchgenlep)    continue;
	      if(po==3&&matchgenlepWZ)  continue;
	      if(po==3&&matchgenlep)    continue;
	      if(po==3&&matchgentau)    continue;
	      bool tightreliso = isotr_relIso04()[n]<0.1;

	      //cout << __LINE__ << endl;
	      //cout << mypr+"reliso03"+mypo << endl;
	      histos[mypr+"reliso03"+mypo]->Fill(isotr_relIso03()[n],weight);
	      histos[mypr+"reliso04"+mypo]->Fill(isotr_relIso04()[n],weight);
	      histos[mypr+"absiso03"+mypo]->Fill(isotr_relIso03()[n]*isotr_p4()[n].Pt(),weight);
	      histos[mypr+"absiso04"+mypo]->Fill(isotr_relIso04()[n]*isotr_p4()[n].Pt(),weight);
	      histos[mypr+"ptrel"   +mypo]->Fill(isotr_ptrel()[n],weight);
	      histos[mypr+"ptratio" +mypo]->Fill(isotr_ptratio()[n],weight);
	      if(tightreliso){
		//cout << mypr+"tightreliso_absiso03"+mypo << endl;
		histos[mypr+"tightreliso_absiso03"+mypo]->Fill(isotr_relIso03()[n]*isotr_p4()[n].Pt(),weight);
		histos[mypr+"tightreliso_absiso04"+mypo]->Fill(isotr_relIso04()[n]*isotr_p4()[n].Pt(),weight);
		histos[mypr+"tightreliso_ptrel"   +mypo]->Fill(isotr_ptrel()[n],weight);
		histos[mypr+"tightreliso_ptratio" +mypo]->Fill(isotr_ptratio()[n],weight);
	      }
	      //cout << __LINE__ << endl;
	    }
	}
      }      
       
    }//event loop
    cout << "File as " << WZcount << " WZ SS leptons and " << Othcount << " other SS leptons." << endl;
    //float integralWZ = 0; float integraloth = 0; bool ct1=true; bool ct2=true;
    for(unsigned int i = 0; i<histonames.size(); ++i){
      string mapname = prefix[1] + histonames[i] + "_"+skimFilePrefix;
      if(mapname.find("tight")!=string::npos) continue;
      int np1 = histos[mapname+postfix[0] ]->GetNbinsX()+1;
      cout << mapname << " " << postfix[0] << " " << histos[mapname+postfix[0] ]->Integral(0,np1)  << " " << postfix[1] << " " << histos[mapname+postfix[1] ]->Integral(0,np1) << " " << postfix[2] << " " << histos[mapname+postfix[2] ]->Integral(0,np1) << " " << postfix[3] << " " << histos[mapname+postfix[3] ]->Integral(0,np1) << endl;
    }
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
  
  string filename = "rootfiles/InvestigateTribosonsObjects_onlyTT1lW1lWW2lWZ3lWWW.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  float integralWZ = 0; float integraloth = 0; bool ct1=true; bool ct2=true;
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    if(h->first.find("fromWZ")!=string::npos){
      if(ct1) { cout << h->first << " " << h->second->Integral() << endl; ct1=false;}
      integralWZ += h->second->Integral();}
    else {
      if(ct2) { cout << h->first << " " << h->second->Integral() << endl; ct2=false;}
      integraloth += h->second->Integral();}
    //string mapname = histonames[i]+"_"+skimFilePrefix;
    //string writename = histonames[i];
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }
  f->Close();
  cout << "Total WZ " << integralWZ << " total Other " << integraloth << endl;
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
