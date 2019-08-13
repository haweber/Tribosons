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
#include "BobaksTuples.cc"
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

  //0: no jet requirement, 1: >=1j, 2: >=2j, 3: Mjj (ATLAS) cut, 4: Mjj+DEta cut
  histonames.push_back("Yield_jreq");                  hbins.push_back(5); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_Mjj_largest_MjjWDR");             hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_Mjj_largest_MjjWDR_Deta2lead");   hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_Mjj_largest_MjjWDR_DetaLarge");   hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_Mjj_largest_MjjWDR_DetaDR");      hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_Mjj_2leading_MjjWDR");            hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_Mjj_2leading_MjjWDR_Deta2lead");  hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_Mjj_2leading_MjjWDR_DetaLarge");  hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_Mjj_2leading_MjjWDR_DetaDR");     hbins.push_back(25); hlow.push_back(0); hup.push_back(500);//x
  histonames.push_back("Yield_DEta_2leading_MjjWDR");  hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_largest");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_largest_MjjWDR");   hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_Mjj_closestDR_Deta1p5"); hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_Mjj_2leading_Deta1p5");  hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_Mjj_closestDR");         hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_Mjj_closestDEta");       hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_Mjj_closestMW");         hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_Mjj_mostcentral");       hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_Mjj_2leading");          hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_Mjj_balanced");          hbins.push_back(22); hlow.push_back(30); hup.push_back(140);
  histonames.push_back("Yield_DR_closestDR");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_closestDEta");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_closestMW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_mostcentral");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_2leading");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_balanced");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_closestDR_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_closestDEta_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_closestMW_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_mostcentral_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_2leading_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DR_balanced_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_closestDR");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_closestDEta");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_closestMW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_mostcentral");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_2leading");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_balanced");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_closestDR_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_closestDEta_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_closestMW_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_mostcentral_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_2leading_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_balanced_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DPhi_closestDR");         hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_closestDEta");       hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_closestMW");         hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_mostcentral");       hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_2leading");          hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_balanced");          hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_closestDR_MjjW");         hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_closestDEta_MjjW");       hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_closestMW_MjjW");         hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_mostcentral_MjjW");       hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_2leading_MjjW");          hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_balanced_MjjW");          hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_centrality_closestDR");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_closestDEta");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_closestMW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_mostcentral");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_2leading");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_balanced");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_closestDR_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_closestDEta_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_closestMW_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_mostcentral_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_2leading_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_centrality_balanced_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_Pt_closestDR");         hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_closestDEta");       hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_closestMW");         hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_mostcentral");       hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_2leading");          hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_balanced");          hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_closestDR_MjjW");         hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_closestDEta_MjjW");       hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_closestMW_MjjW");         hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_mostcentral_MjjW");       hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_2leading_MjjW");          hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_Pt_balanced_MjjW");          hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_PtRatio_closestDR");         hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_closestDEta");       hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_closestMW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_mostcentral");       hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_2leading");          hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_balanced");          hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_closestDR_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_closestDEta_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_closestMW_MjjW");         hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_mostcentral_MjjW");       hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_2leading_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_PtRatio_balanced_MjjW");          hbins.push_back(20); hlow.push_back(0); hup.push_back(1);

  histonames.push_back("Yield_DR_closestDR_noMjjW");                hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DEta_closestDR_noMjjW");              hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_DPhi_closestDR_noMjjW");              hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_centrality_closestDR_noMjjW");        hbins.push_back(20); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Yield_Pt_closestDR_noMjjW");                hbins.push_back(25); hlow.push_back(0); hup.push_back(250);
  histonames.push_back("Yield_PtRatio_closestDR_noMjjW");           hbins.push_back(20); hlow.push_back(0); hup.push_back(1);
  histonames.push_back("Yield_DPhi_jj_ll_closestDR_noMjjW");        hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_jj_MET_closestDR_noMjjW");       hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_jj_llMET_closestDR_noMjjW");     hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_PtRatio_jj_ll_closestDR_noMjjW");     hbins.push_back(30); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("Yield_PtRatio_jj_MET_closestDR_noMjjW");    hbins.push_back(30); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("Yield_PtRatio_jj_llMET_closestDR_noMjjW");  hbins.push_back(30); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("Yield_MT_jj_ll_closestDR_noMjjW");          hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT_jj_MET_closestDR_noMjjW");         hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT_jj_llMET_closestDR_noMjjW");       hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT2_jj_ll_MET_closestDR_noMjjW");     hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT2ml_jj_ll_MET_closestDR_noMjjW");   hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_DPhi_jj_ll_closestDR_MjjW");        hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_jj_MET_closestDR_MjjW");       hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_DPhi_jj_llMET_closestDR_MjjW");     hbins.push_back(32); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("Yield_PtRatio_jj_ll_closestDR_MjjW");     hbins.push_back(30); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("Yield_PtRatio_jj_MET_closestDR_MjjW");    hbins.push_back(30); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("Yield_PtRatio_jj_llMET_closestDR_MjjW");  hbins.push_back(30); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("Yield_MT_jj_ll_closestDR_MjjW");          hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT_jj_MET_closestDR_MjjW");         hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT_jj_llMET_closestDR_MjjW");       hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT2_jj_ll_MET_closestDR_MjjW");     hbins.push_back(25); hlow.push_back(0); hup.push_back(500);
  histonames.push_back("Yield_MT2ml_jj_ll_MET_closestDR_MjjW");   hbins.push_back(25); hlow.push_back(0); hup.push_back(300);


  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    for(int b = 0; b < 3; ++b){
	string mapname = histonames[i] + "_"+skimFilePrefix;
	if(b==0) mapname = "ee"+mapname;
	if(b==1) mapname = "em"+mapname;
	if(b==2) mapname = "mm"+mapname;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
	//cout << mapname << endl;
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
      
      double weight = evt_scale1fb()*35.9;
      //if(nEventsTotal==0) cout << weight << endl; 

      if(nVert()<0)              continue;
      if(nlep()<2)               continue;
      //if(nTaus20()!=0)           continue;
      //if(nBJetMedium()!=0)       continue;
      if(nisoTrack_mt2()!=0)  continue;

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
      
      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());


      //cout << __LINE__ << endl;
      int nj(0),nb(0);
      int nj20(0),nj30(0);
      vector<int> ij, ij25, ij20, ij30;
      vector<int>  bm30, bl20;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.5) ij.push_back(n);
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.5) ij30.push_back(n);
	if(jets_p4()[n].Pt()>25&&fabs(jets_p4()[n].Eta())<2.5) ij25.push_back(n);
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.5) ij20.push_back(n);
	if(jets_csv()[n]>0.8484&&jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.4) bm30.push_back(n);
	if(jets_csv()[n]>0.5426&&jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.4) bl20.push_back(n);
      }
      
      //cout << __LINE__ << endl;
      float Mjj_2p5 = -1; float Detajj_2p5 = -1;
      if(ij.size()>1){
	Mjj_2p5 = (jets_p4()[ij[0] ]+jets_p4()[ij[1] ]).M();
	Detajj_2p5 = dEta(jets_p4()[ij[0] ],jets_p4()[ij[1] ]);
      }
      bool passMDetajj_2p5 = false;
      bool passMDetajj_3p0 = false;
      bool passMDetajj_5p0 = false;
      if(ij.size()>=2&&fabs(Detajj_2p5)<1.5&&fabs(Mjj_2p5-85.)<20.) passMDetajj_2p5 = true;
      
      //cout << __LINE__ << endl;
     vector<int> i0, v, iSS, i3l;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	if(abs(lep_pdgId()[i])==11&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSS.push_back(i);
	if(abs(lep_pdgId()[i])==13&&lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSS.push_back(i);
	if(abs(lep_pdgId()[i])==11&&/*lep_tightCharge()[i]==2&&*/lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3l.push_back(i);
	if(abs(lep_pdgId()[i])==13&&lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3l.push_back(i);
	if(lep_relIso03EA()[i]<=0.1 && fabs(lep_ip3d()[i])<=0.015&&lep_p4()[i].Pt()>30) i0.push_back(i);
	//if(lep_relIso03EA()[i]<=0.1 && fabs(lep_ip3d()[i])<=0.015&&lep_p4()[i].Pt()>25) i0.push_back(i);//loosened
	if(lep_relIso03EA()[i]<=0.3) v.push_back(i);
      }
      //cout << __LINE__ << endl;

      bool preselectSS = false;//all cuts that are not ee/mumu/emu specific and unrelated to jets
      if(bl20.size()==0&&nisoTrack_mt2()==0&&nlep()==2&&iSS.size()>=2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&lep_p4()[iSS[0] ].Pt()>30&&lep_p4()[iSS[1] ].Pt()>30&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) preselectSS = true;
      //if(bm30.size()==0&&nisoTrack_mt2()==0&&v.size()==2&&iSS.size()>=2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&lep_p4()[iSS[0] ].Pt()>25&&lep_p4()[iSS[1] ].Pt()>25&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) preselectSS = true;//loosened
      if(!preselectSS) continue;
      bool preselectee(false), preselectemu(false), preselectmumu(false);
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[0] ]==2&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>55.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) preselectee = true;
      //if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[0] ]==2&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>45.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) preselectee = true;//loosened
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>55.) preselectemu = true;
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&lep_tightCharge()[iSS[0] ]==2&&met_pt()>55.) preselectemu = true;
      //if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>45.) preselectemu = true;//loosened
      //if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&lep_tightCharge()[iSS[0] ]==2&&met_pt()>45.) preselectemu = true;//loosened
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) preselectmumu = true;
      bool passoneSS = false;
      if(preselectee||preselectemu||preselectmumu) passoneSS = true;
      if(!passoneSS) continue;
      string SS = "";
      if(preselectee)   SS = "ee";
      if(preselectemu)  SS = "em";
      if(preselectmumu) SS = "mm";
      
      /*if(ij.size()>=0)*/                                          histos[SS+"Yield_jreq_"+sn]->Fill(0.,weight);
      if(ij.size()>=1)                                              histos[SS+"Yield_jreq_"+sn]->Fill(1.,weight);
      if(ij.size()>=2)                                              histos[SS+"Yield_jreq_"+sn]->Fill(2.,weight);
      if(ij.size()>=2&&fabs(Mjj_2p5-85.)<20.)                       histos[SS+"Yield_jreq_"+sn]->Fill(3.,weight);
      if(ij.size()>=2&&fabs(Mjj_2p5-85.)<20.&&fabs(Detajj_2p5)<1.5) histos[SS+"Yield_jreq_"+sn]->Fill(4.,weight);
	
      if(ij.size()<2) continue;
      //if(ij30.size()<1) continue;

      int jmE1(-1), jmE2(-1);
      int jDR1(-1), jDR2(-1);
      int jDE1(-1), jDE2(-1);
      int jMW1(-1), jMW2(-1);
      int jct1(-1), jct2(-1);
      int jpT1(-1), jpT2(-1);
      int jbl1(-1), jbl2(-1);
      float minDR(999.), minDEta(999.), minDMW(999), minCent(999.), minBal(999.), maxDEta(-999.);
      jpT1 = ij[0]; jpT2 = ij[1];//vector was pT ordered.
      for(unsigned int j1 = 0; j1<ij.size();++j1){
	for(unsigned int j2 = j1+1; j2<ij.size();++j2){
	  if(dEta(jets_p4()[ij[j1] ], jets_p4()[ij[j2] ])>maxDEta){
	    maxDEta = dEta(jets_p4()[ij[j1] ], jets_p4()[ij[j2] ]);
	    jmE1 = ij[j1]; jmE2 = ij[j2];
	  }
	  if(dR(jets_p4()[ij[j1] ], jets_p4()[ij[j2] ])<minDR){
	    minDR = dR(jets_p4()[ij[j1] ], jets_p4()[ij[j2] ]);
	    jDR1 = ij[j1]; jDR2 = ij[j2];
	  }
	  if(dEta(jets_p4()[ij[j1] ], jets_p4()[ij[j2] ])<minDEta){
	    minDEta = dEta(jets_p4()[ij[j1] ], jets_p4()[ij[j2] ]);
	    jDE1 = ij[j1]; jDE2 = ij[j2];
	  }
	  if(fabs((jets_p4()[ij[j1] ]+jets_p4()[ij[j2] ]).M()-80.)<minDMW){
	    minDMW = fabs((jets_p4()[ij[j1] ]+jets_p4()[ij[j2] ]).M()-80.);
	    jMW1 = ij[j1]; jMW2 = ij[j2];
	  }
	  if(sqrt(pow(jets_p4()[ij[j1] ].Eta(),2)+pow(jets_p4()[ij[j2] ].Eta(),2))<minCent){
	    minCent = sqrt(pow(jets_p4()[ij[j1] ].Eta(),2)+pow(jets_p4()[ij[j2] ].Eta(),2));
	    jct1 = ij[j1]; jct2 = ij[j2];
	  }
	  if(fabs((jets_p4()[ij[j2] ].Pt()/jets_p4()[ij[j1] ].Pt())-1.)<minBal){
	    minBal = fabs((jets_p4()[ij[j2] ].Pt()/jets_p4()[ij[j1] ].Pt())-1.);
	    jbl1 = ij[j1]; jbl2 = ij[j2];
	  }
	}
      }
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) {
	histos[SS+"Yield_Mjj_largest_MjjWDR_"  +sn]->Fill((jets_p4()[jmE1]+jets_p4()[jmE2]).M(),weight);
	if(dEta(jets_p4()[jmE1],jets_p4()[jmE2])<2.5) histos[SS+"Yield_Mjj_largest_MjjWDR_DetaLarge_"  +sn]->Fill((jets_p4()[jmE1]+jets_p4()[jmE2]).M(),weight);
	if(dEta(jets_p4()[jpT1],jets_p4()[jpT2])<1.5) histos[SS+"Yield_Mjj_largest_MjjWDR_Deta2lead_"  +sn]->Fill((jets_p4()[jmE1]+jets_p4()[jmE2]).M(),weight);
	if(dEta(jets_p4()[jDR1],jets_p4()[jDR2])<1.5) histos[SS+"Yield_Mjj_largest_MjjWDR_DetaDR_"  +sn]->Fill((jets_p4()[jmE1]+jets_p4()[jmE2]).M(),weight);
	histos[SS+"Yield_Mjj_2leading_MjjWDR_"  +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).M(),weight);
	if(dEta(jets_p4()[jmE1],jets_p4()[jmE2])<2.5) histos[SS+"Yield_Mjj_2leading_MjjWDR_DetaLarge_"  +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).M(),weight);
	if(dEta(jets_p4()[jpT1],jets_p4()[jpT2])<1.5) histos[SS+"Yield_Mjj_2leading_MjjWDR_Deta2lead_"  +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).M(),weight);
	if(dEta(jets_p4()[jDR1],jets_p4()[jDR2])<1.5) histos[SS+"Yield_Mjj_2leading_MjjWDR_DetaDR_"  +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).M(),weight);
      }
      if(dEta(jets_p4()[jDR1],jets_p4()[jDR2])<1.5) histos[SS+"Yield_Mjj_closestDR_Deta1p5_"  +sn]->Fill((jets_p4()[jDR1]+jets_p4()[jDR2]).M(),weight);
      if(dEta(jets_p4()[jpT1],jets_p4()[jpT2])<1.5) histos[SS+"Yield_Mjj_2leading_Deta1p5_"   +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).M(),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_DEta_2leading_MjjWDR_"  +sn]->Fill(dEta(jets_p4()[jpT1],jets_p4()[jpT2]),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_DEta_largest_MjjWDR_"  +sn]->Fill(dEta(jets_p4()[jmE1],jets_p4()[jmE2]),weight);
      histos[SS+"Yield_DEta_largest_"  +sn]->Fill(dR(jets_p4()[jmE1],jets_p4()[jmE2]),weight);

      histos[SS+"Yield_Mjj_closestDR_"  +sn]->Fill((jets_p4()[jDR1]+jets_p4()[jDR2]).M(),weight);
      histos[SS+"Yield_Mjj_closestDEta_"+sn]->Fill((jets_p4()[jDE1]+jets_p4()[jDE2]).M(),weight);
      histos[SS+"Yield_Mjj_closestMW_"  +sn]->Fill((jets_p4()[jMW1]+jets_p4()[jMW2]).M(),weight);
      histos[SS+"Yield_Mjj_mostcentral_"+sn]->Fill((jets_p4()[jct1]+jets_p4()[jct2]).M(),weight);
      histos[SS+"Yield_Mjj_2leading_"   +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).M(),weight);
      histos[SS+"Yield_Mjj_balanced_"   +sn]->Fill((jets_p4()[jbl1]+jets_p4()[jbl2]).M(),weight);

      histos[SS+"Yield_DR_closestDR_"  +sn]->Fill(dR(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
      histos[SS+"Yield_DR_closestDEta_"+sn]->Fill(dR(jets_p4()[jDE1],jets_p4()[jDE2]),weight);
      histos[SS+"Yield_DR_closestMW_"  +sn]->Fill(dR(jets_p4()[jMW1],jets_p4()[jMW2]),weight);
      histos[SS+"Yield_DR_mostcentral_"+sn]->Fill(dR(jets_p4()[jct1],jets_p4()[jct2]),weight);
      histos[SS+"Yield_DR_2leading_"   +sn]->Fill(dR(jets_p4()[jpT1],jets_p4()[jpT2]),weight);
      histos[SS+"Yield_DR_balanced_"   +sn]->Fill(dR(jets_p4()[jbl1],jets_p4()[jbl2]),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_DR_closestDR_MjjW_"  +sn]->Fill(dR(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
      if(fabs((jets_p4()[jDE1]+jets_p4()[jDE2]).M()-80.)<20.) histos[SS+"Yield_DR_closestDEta_MjjW_"+sn]->Fill(dR(jets_p4()[jDE1],jets_p4()[jDE2]),weight);
      if(fabs((jets_p4()[jMW1]+jets_p4()[jMW2]).M()-80.)<20.) histos[SS+"Yield_DR_closestMW_MjjW_"  +sn]->Fill(dR(jets_p4()[jMW1],jets_p4()[jMW2]),weight);
      if(fabs((jets_p4()[jct1]+jets_p4()[jct2]).M()-80.)<20.) histos[SS+"Yield_DR_mostcentral_MjjW_"+sn]->Fill(dR(jets_p4()[jct1],jets_p4()[jct2]),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_DR_2leading_MjjW_"   +sn]->Fill(dR(jets_p4()[jpT1],jets_p4()[jpT2]),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_DR_balanced_MjjW_"   +sn]->Fill(dR(jets_p4()[jbl1],jets_p4()[jbl2]),weight);
      histos[SS+"Yield_DEta_closestDR_"  +sn]->Fill(dEta(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
      histos[SS+"Yield_DEta_closestDEta_"+sn]->Fill(dEta(jets_p4()[jDE1],jets_p4()[jDE2]),weight);
      histos[SS+"Yield_DEta_closestMW_"  +sn]->Fill(dEta(jets_p4()[jMW1],jets_p4()[jMW2]),weight);
      histos[SS+"Yield_DEta_mostcentral_"+sn]->Fill(dEta(jets_p4()[jct1],jets_p4()[jct2]),weight);
      histos[SS+"Yield_DEta_2leading_"   +sn]->Fill(dEta(jets_p4()[jpT1],jets_p4()[jpT2]),weight);
      histos[SS+"Yield_DEta_balanced_"   +sn]->Fill(dEta(jets_p4()[jbl1],jets_p4()[jbl2]),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_DEta_closestDR_MjjW_"  +sn]->Fill(dEta(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
      if(fabs((jets_p4()[jDE1]+jets_p4()[jDE2]).M()-80.)<20.) histos[SS+"Yield_DEta_closestDEta_MjjW_"+sn]->Fill(dEta(jets_p4()[jDE1],jets_p4()[jDE2]),weight);
      if(fabs((jets_p4()[jMW1]+jets_p4()[jMW2]).M()-80.)<20.) histos[SS+"Yield_DEta_closestMW_MjjW_"  +sn]->Fill(dEta(jets_p4()[jMW1],jets_p4()[jMW2]),weight);
      if(fabs((jets_p4()[jct1]+jets_p4()[jct2]).M()-80.)<20.) histos[SS+"Yield_DEta_mostcentral_MjjW_"+sn]->Fill(dEta(jets_p4()[jct1],jets_p4()[jct2]),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_DEta_2leading_MjjW_"   +sn]->Fill(dEta(jets_p4()[jpT1],jets_p4()[jpT2]),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_DEta_balanced_MjjW_"   +sn]->Fill(dEta(jets_p4()[jbl1],jets_p4()[jbl2]),weight);
      histos[SS+"Yield_DPhi_closestDR_"  +sn]->Fill(dPhi(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
      histos[SS+"Yield_DPhi_closestDEta_"+sn]->Fill(dPhi(jets_p4()[jDE1],jets_p4()[jDE2]),weight);
      histos[SS+"Yield_DPhi_closestMW_"  +sn]->Fill(dPhi(jets_p4()[jMW1],jets_p4()[jMW2]),weight);
      histos[SS+"Yield_DPhi_mostcentral_"+sn]->Fill(dPhi(jets_p4()[jct1],jets_p4()[jct2]),weight);
      histos[SS+"Yield_DPhi_2leading_"   +sn]->Fill(dPhi(jets_p4()[jpT1],jets_p4()[jpT2]),weight);
      histos[SS+"Yield_DPhi_balanced_"   +sn]->Fill(dPhi(jets_p4()[jbl1],jets_p4()[jbl2]),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_DPhi_closestDR_MjjW_"  +sn]->Fill(dPhi(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
      if(fabs((jets_p4()[jDE1]+jets_p4()[jDE2]).M()-80.)<20.) histos[SS+"Yield_DPhi_closestDEta_MjjW_"+sn]->Fill(dPhi(jets_p4()[jDE1],jets_p4()[jDE2]),weight);
      if(fabs((jets_p4()[jMW1]+jets_p4()[jMW2]).M()-80.)<20.) histos[SS+"Yield_DPhi_closestMW_MjjW_"  +sn]->Fill(dPhi(jets_p4()[jMW1],jets_p4()[jMW2]),weight);
      if(fabs((jets_p4()[jct1]+jets_p4()[jct2]).M()-80.)<20.) histos[SS+"Yield_DPhi_mostcentral_MjjW_"+sn]->Fill(dPhi(jets_p4()[jct1],jets_p4()[jct2]),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_DPhi_2leading_MjjW_"   +sn]->Fill(dPhi(jets_p4()[jpT1],jets_p4()[jpT2]),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_DPhi_balanced_MjjW_"   +sn]->Fill(dPhi(jets_p4()[jbl1],jets_p4()[jbl2]),weight);
      histos[SS+"Yield_centrality_closestDR_"  +sn]->Fill(sqrt(pow(jets_p4()[jDR1].Eta(),2)+pow(jets_p4()[jDR2].Eta(),2)),weight);
      histos[SS+"Yield_centrality_closestDEta_"+sn]->Fill(sqrt(pow(jets_p4()[jDE1].Eta(),2)+pow(jets_p4()[jDE2].Eta(),2)),weight);
      histos[SS+"Yield_centrality_closestMW_"  +sn]->Fill(sqrt(pow(jets_p4()[jMW1].Eta(),2)+pow(jets_p4()[jMW2].Eta(),2)),weight);
      histos[SS+"Yield_centrality_mostcentral_"+sn]->Fill(sqrt(pow(jets_p4()[jct1].Eta(),2)+pow(jets_p4()[jct2].Eta(),2)),weight);
      histos[SS+"Yield_centrality_2leading_"   +sn]->Fill(sqrt(pow(jets_p4()[jpT1].Eta(),2)+pow(jets_p4()[jpT2].Eta(),2)),weight);
      histos[SS+"Yield_centrality_balanced_"   +sn]->Fill(sqrt(pow(jets_p4()[jbl1].Eta(),2)+pow(jets_p4()[jbl2].Eta(),2)),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_centrality_closestDR_MjjW_"  +sn]->Fill(sqrt(pow(jets_p4()[jDR1].Eta(),2)+pow(jets_p4()[jDR2].Eta(),2)),weight);
      if(fabs((jets_p4()[jDE1]+jets_p4()[jDE2]).M()-80.)<20.) histos[SS+"Yield_centrality_closestDEta_MjjW_"+sn]->Fill(sqrt(pow(jets_p4()[jDE1].Eta(),2)+pow(jets_p4()[jDE2].Eta(),2)),weight);
      if(fabs((jets_p4()[jMW1]+jets_p4()[jMW2]).M()-80.)<20.) histos[SS+"Yield_centrality_closestMW_MjjW_"  +sn]->Fill(sqrt(pow(jets_p4()[jMW1].Eta(),2)+pow(jets_p4()[jMW2].Eta(),2)),weight);
      if(fabs((jets_p4()[jct1]+jets_p4()[jct2]).M()-80.)<20.) histos[SS+"Yield_centrality_mostcentral_MjjW_"+sn]->Fill(sqrt(pow(jets_p4()[jct1].Eta(),2)+pow(jets_p4()[jct2].Eta(),2)),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_centrality_2leading_MjjW_"   +sn]->Fill(sqrt(pow(jets_p4()[jpT1].Eta(),2)+pow(jets_p4()[jpT2].Eta(),2)),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_centrality_balanced_MjjW_"   +sn]->Fill(sqrt(pow(jets_p4()[jbl1].Eta(),2)+pow(jets_p4()[jbl2].Eta(),2)),weight);
      histos[SS+"Yield_Pt_closestDR_"  +sn]->Fill((jets_p4()[jDR1]+jets_p4()[jDR2]).Pt(),weight);
      histos[SS+"Yield_Pt_closestDEta_"+sn]->Fill((jets_p4()[jDE1]+jets_p4()[jDE2]).Pt(),weight);
      histos[SS+"Yield_Pt_closestMW_"  +sn]->Fill((jets_p4()[jMW1]+jets_p4()[jMW2]).Pt(),weight);
      histos[SS+"Yield_Pt_mostcentral_"+sn]->Fill((jets_p4()[jct1]+jets_p4()[jct2]).Pt(),weight);
      histos[SS+"Yield_Pt_2leading_"   +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).Pt(),weight);
      histos[SS+"Yield_Pt_balanced_"   +sn]->Fill((jets_p4()[jbl1]+jets_p4()[jbl2]).Pt(),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_Pt_closestDR_MjjW_"  +sn]->Fill((jets_p4()[jDR1]+jets_p4()[jDR2]).Pt(),weight);
      if(fabs((jets_p4()[jDE1]+jets_p4()[jDE2]).M()-80.)<20.) histos[SS+"Yield_Pt_closestDEta_MjjW_"+sn]->Fill((jets_p4()[jDE1]+jets_p4()[jDE2]).Pt(),weight);
      if(fabs((jets_p4()[jMW1]+jets_p4()[jMW2]).M()-80.)<20.) histos[SS+"Yield_Pt_closestMW_MjjW_"  +sn]->Fill((jets_p4()[jMW1]+jets_p4()[jMW2]).Pt(),weight);
      if(fabs((jets_p4()[jct1]+jets_p4()[jct2]).M()-80.)<20.) histos[SS+"Yield_Pt_mostcentral_MjjW_"+sn]->Fill((jets_p4()[jct1]+jets_p4()[jct2]).Pt(),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_Pt_2leading_MjjW_"   +sn]->Fill((jets_p4()[jpT1]+jets_p4()[jpT2]).Pt(),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_Pt_balanced_MjjW_"   +sn]->Fill((jets_p4()[jbl1]+jets_p4()[jbl2]).Pt(),weight);
      histos[SS+"Yield_PtRatio_closestDR_"  +sn]->Fill(jets_p4()[jDR2].Pt()/jets_p4()[jDR1].Pt(),weight);
      histos[SS+"Yield_PtRatio_closestDEta_"+sn]->Fill(jets_p4()[jDE2].Pt()/jets_p4()[jDE1].Pt(),weight);
      histos[SS+"Yield_PtRatio_closestMW_"  +sn]->Fill(jets_p4()[jMW2].Pt()/jets_p4()[jMW1].Pt(),weight);
      histos[SS+"Yield_PtRatio_mostcentral_"+sn]->Fill(jets_p4()[jct2].Pt()/jets_p4()[jct1].Pt(),weight);
      histos[SS+"Yield_PtRatio_2leading_"   +sn]->Fill(jets_p4()[jpT2].Pt()/jets_p4()[jpT1].Pt(),weight);
      histos[SS+"Yield_PtRatio_balanced_"   +sn]->Fill(jets_p4()[jbl1].Pt()/jets_p4()[jbl2].Pt(),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.) histos[SS+"Yield_PtRatio_closestDR_MjjW_"  +sn]->Fill(jets_p4()[jDR2].Pt()/jets_p4()[jDR1].Pt(),weight);
      if(fabs((jets_p4()[jDE1]+jets_p4()[jDE2]).M()-80.)<20.) histos[SS+"Yield_PtRatio_closestDEta_MjjW_"+sn]->Fill(jets_p4()[jDE2].Pt()/jets_p4()[jDE1].Pt(),weight);
      if(fabs((jets_p4()[jMW1]+jets_p4()[jMW2]).M()-80.)<20.) histos[SS+"Yield_PtRatio_closestMW_MjjW_"  +sn]->Fill(jets_p4()[jMW2].Pt()/jets_p4()[jMW1].Pt(),weight);
      if(fabs((jets_p4()[jct1]+jets_p4()[jct2]).M()-80.)<20.) histos[SS+"Yield_PtRatio_mostcentral_MjjW_"+sn]->Fill(jets_p4()[jct2].Pt()/jets_p4()[jct1].Pt(),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_PtRatio_2leading_MjjW_"   +sn]->Fill(jets_p4()[jpT2].Pt()/jets_p4()[jpT1].Pt(),weight);
      if(fabs((jets_p4()[jpT1]+jets_p4()[jpT2]).M()-80.)<20.) histos[SS+"Yield_PtRatio_balanced_MjjW_"   +sn]->Fill(jets_p4()[jbl2].Pt()/jets_p4()[jbl1].Pt(),weight);
      if(fabs((jets_p4()[jDR1]+jets_p4()[jDR2]).M()-80.)<20.){
	histos[SS+"Yield_DPhi_jj_ll_closestDR_MjjW_"       +sn]->Fill(dPhi(jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]),          weight);
	histos[SS+"Yield_DPhi_jj_MET_closestDR_MjjW_"      +sn]->Fill(dPhi(jets_p4()[jDR1]+jets_p4()[jDR2],                                         MET),      weight);
	histos[SS+"Yield_DPhi_jj_llMET_closestDR_MjjW_"    +sn]->Fill(dPhi(jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]+MET),      weight);
	histos[SS+"Yield_PtRatio_jj_ll_closestDR_MjjW_"    +sn]->Fill(    (jets_p4()[jDR1]+jets_p4()[jDR2]).Pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt(),     weight);
	histos[SS+"Yield_PtRatio_jj_MET_closestDR_MjjW_"   +sn]->Fill(    (jets_p4()[jDR1]+jets_p4()[jDR2]).Pt()/(                                  MET).Pt(), weight);
	histos[SS+"Yield_PtRatio_jj_llMET_closestDR_MjjW_" +sn]->Fill(    (jets_p4()[jDR1]+jets_p4()[jDR2]).Pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]+MET).Pt(), weight);
	histos[SS+"Yield_MT_jj_ll_closestDR_MjjW_"         +sn]->Fill(mT(  jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]),          weight);
	histos[SS+"Yield_MT_jj_MET_closestDR_MjjW_"        +sn]->Fill(mT(  jets_p4()[jDR1]+jets_p4()[jDR2],                                         MET),      weight);
	histos[SS+"Yield_MT_jj_llMET_closestDR_MjjW_"      +sn]->Fill(mT(  jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]+MET),      weight);
	histos[SS+"Yield_MT2_jj_ll_MET_closestDR_MjjW_"    +sn]->Fill(MT2( jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ],MET,true), weight);
	histos[SS+"Yield_MT2ml_jj_ll_MET_closestDR_MjjW_"  +sn]->Fill(MT2( jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ],MET,false),weight);
      } else {
	histos[SS+"Yield_DPhi_jj_ll_closestDR_noMjjW_"       +sn]->Fill(dPhi(jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]),          weight);
	histos[SS+"Yield_DPhi_jj_MET_closestDR_noMjjW_"      +sn]->Fill(dPhi(jets_p4()[jDR1]+jets_p4()[jDR2],                                         MET),      weight);
	histos[SS+"Yield_DPhi_jj_llMET_closestDR_noMjjW_"    +sn]->Fill(dPhi(jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]+MET),      weight);
	histos[SS+"Yield_PtRatio_jj_ll_closestDR_noMjjW_"    +sn]->Fill(    (jets_p4()[jDR1]+jets_p4()[jDR2]).Pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt(),     weight);
	histos[SS+"Yield_PtRatio_jj_MET_closestDR_noMjjW_"   +sn]->Fill(    (jets_p4()[jDR1]+jets_p4()[jDR2]).Pt()/(                                  MET).Pt(), weight);
	histos[SS+"Yield_PtRatio_jj_llMET_closestDR_noMjjW_" +sn]->Fill(    (jets_p4()[jDR1]+jets_p4()[jDR2]).Pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]+MET).Pt(), weight);
	histos[SS+"Yield_MT_jj_ll_closestDR_noMjjW_"         +sn]->Fill(mT(  jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]),          weight);
	histos[SS+"Yield_MT_jj_MET_closestDR_noMjjW_"        +sn]->Fill(mT(  jets_p4()[jDR1]+jets_p4()[jDR2],                                         MET),      weight);
	histos[SS+"Yield_MT_jj_llMET_closestDR_noMjjW_"      +sn]->Fill(mT(  jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]+MET),      weight);
	histos[SS+"Yield_MT2_jj_ll_MET_closestDR_noMjjW_"    +sn]->Fill(MT2( jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ],MET,true), weight);
	histos[SS+"Yield_MT2ml_jj_ll_MET_closestDR_noMjjW_"  +sn]->Fill(MT2( jets_p4()[jDR1]+jets_p4()[jDR2],       lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ],MET,false),weight);
	histos[SS+"Yield_PtRatio_closestDR_noMjjW_"  +sn]->Fill(jets_p4()[jDR2].Pt()/jets_p4()[jDR1].Pt(),weight);
	histos[SS+"Yield_Pt_closestDR_noMjjW_"  +sn]->Fill((jets_p4()[jDR1]+jets_p4()[jDR2]).Pt(),weight);
	histos[SS+"Yield_centrality_closestDR_noMjjW_"  +sn]->Fill(sqrt(pow(jets_p4()[jDR1].Eta(),2)+pow(jets_p4()[jDR2].Eta(),2)),weight);
	histos[SS+"Yield_DPhi_closestDR_noMjjW_"  +sn]->Fill(dPhi(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
	histos[SS+"Yield_DEta_closestDR_noMjjW_"  +sn]->Fill(dEta(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
	histos[SS+"Yield_DR_closestDR_noMjjW_"  +sn]->Fill(dR(jets_p4()[jDR1],jets_p4()[jDR2]),weight);
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
  
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //add overflow
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
    //add underflow
    h->second->SetBinContent(1, h->second->GetBinContent(1)+ h->second->GetBinContent(0) );
    h->second->SetBinError(1, sqrt(pow(h->second->GetBinError(1),2)+pow(h->second->GetBinError(0),2) ) );
  }

  //loosened up lepton-pT to pT > 25 GeV
  //loosened b-veto to medium 30 GeV
  //relax MET cut to 45 GeV
  //string filename = "rootfiles/StupidStuff_Mjj_loosened.root";
  string filename = "rootfiles/StupidStuff_Mjj.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //string mapname = histonames[i]+"_"+skimFilePrefix;
    //string writename = histonames[i];
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }  f->Close();
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
