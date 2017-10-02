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

  histonames.push_back("YieldsSR_jesup");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_jesdown");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsCR_dropMjj_jesup");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_dropMjj_jesdown");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_using3lorMll_jesup");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_using3lorMll_jesdown");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_dropcutsbutNJ_jesup");               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_dropcutsbutNJ_jesdown");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Mjj_CRlike_allSS_jesup");                     hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_allSS_jesdown");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRloose_allSS_jesup");                    hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRloose_allSS_jesdown");                  hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_SRlike_allSS_jesup");                     hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("Mjj_SRlike_allSS_jesdown");                   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data
  histonames.push_back("YieldsCR_SSany_dropMjj_jesup");               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_dropMjj_jesdown");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_dropMjj");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_dropMjj_test");                 hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_dropMjj_butnoMll");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_dropMjj_rawweight");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_jesup");                       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_jesdown");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Mjj_CRlike_SSany_allSS_jesup");               hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_SSany_allSS_jesdown");             hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_SSany_allSS");                     hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRloose_butnotlike_SSany_allSS");         hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRloose_butnotlike_SSany_allSS_wMjjL");   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRloose_butnotlike_SSany_allSS_wMjjLDeta");hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("YieldsCR_SSany_dropMjj_lepSFup");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_SSany_dropMjj_lepSFdn");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lepSFup");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lepSFdn");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  
  histonames.push_back("YieldsSR_lowMSFOS");                 hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_lowMSFOS_Mlll");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_lowMSFOS_dRllmin");         hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR");                          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_raw");                      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_rawweight");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_preselection");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_invertdPhiPtfor3Mjjfor2l"); hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertdPhiPtfor3Mjjfor2l"); hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_using3lorMll");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_dropcutsbutNJisotr");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_dropcutsbutNJvetolep");     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_dropcutsbutNJ");            hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCRnoMll_dropcutsbutNJisotr");  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCRnoMll_dropcutsbutNJvetolep");hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCRnoMll_dropcutsbutNJ");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_dropMjj");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lowMET");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lowMETinvertdPhiPtMjj");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_lowMET");                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_lowMETinvertdPhiPtMjj");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  
  histonames.push_back("YieldsSR_invertDPhi");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertDPhi");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertPt");               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertPt");               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertMET");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertMET");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertDPhiPt");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertDPhiPt");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertDPhiOrPt");         hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertDPhiOrPt");         hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertMETDPhiPt");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertMETDPhiPt");        hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_invertMETDPhiOrPt");      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  
  histonames.push_back("Mlll_SRlike_allSFOS");               hbins.push_back(15);hlow.push_back(30);hup.push_back(180);
  histonames.push_back("minDRll_SRlike_allSFOS");            hbins.push_back(25);hlow.push_back(0.);hup.push_back(2.5);
  histonames.push_back("Mjj_lowMET_SRlike_allSS");           hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_lowMET_CRlike_allSS");           hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_nolowMll_allSS");         hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mll_invertdPhiPt_1SFOS");            hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_invertdPhiOrPt_1SFOS");          hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_invertdPhi_1SFOS");              hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_invertPt_1SFOS");                hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_invertMET_1SFOS");               hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_invertMETdPhiPt_1SFOS");         hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_invertdPhiPt_2SFOS");            hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertdPhiOrPt_2SFOS");          hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertdPhi_2SFOS");              hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertPt_2SFOS");                hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertMET_2SFOS");               hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertMETdPhiPt_2SFOS");         hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_inverteitherMETdPhiPt_1SFOS");   hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_inverteitherMETdPhiPt_2SFOS");   hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_SRlike_1SFOS");                  hbins.push_back(14);hlow.push_back(22);hup.push_back(176);//blind the data
  histonames.push_back("Mll_SRlike_2SFOS");                  hbins.push_back(14);hlow.push_back(30);hup.push_back(170);//blind the data
  histonames.push_back("Mll_ge2j_1SFOS");                    hbins.push_back(14);hlow.push_back(22);hup.push_back(176);
  histonames.push_back("Mll_ge2j_2SFOS");                    hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_SRlike_MllclosestZ_2SFOS");      hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_SRlike_MllfurthestZ_2SFOS");     hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mjj_CRlike_allSS_SSpairnotlead2");   hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_allSS_SSpairnotlead2v2"); hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_SRlike_allSS");                  hbins.push_back(12);hlow.push_back(20);hup.push_back(260);//blind the data

  histonames.push_back("Mll_invertdPhiPt_allSFOS");            hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertdPhiOrPt_allSFOS");          hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertdPhi_allSFOS");              hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertPt_allSFOS");                hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertMET_allSFOS");               hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_invertMETdPhiPt_allSFOS");         hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_inverteitherMETdPhiPt_allSFOS");   hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_SRlike_allSFOS");                  hbins.push_back(14);hlow.push_back(30);hup.push_back(170);//blind the data
  histonames.push_back("Mll_ge2j_allSFOS");                    hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  
  histonames.push_back("Mll_CRlike_allSS");                  hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_CRlike_allSS_SSpairnotlead2");   hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_CRlike_allSS_SSpairnotlead2v2"); hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_CRnoMjj_allSS");                 hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_CRnoMjj_allSS_SSpairnotlead2");  hbins.push_back(14);hlow.push_back(30);hup.push_back(170);
  histonames.push_back("Mll_CRnoMjj_allSS_SSpairnotlead2v2");hbins.push_back(14);hlow.push_back(30);hup.push_back(170);

  histonames.push_back("Mjj_CRlike_allSS");                  hbins.push_back(12);hlow.push_back( 20);hup.push_back(260);
  histonames.push_back("Mjj_CRloose_allSS");                 hbins.push_back(12);hlow.push_back( 20);hup.push_back(260);
  histonames.push_back("Mll_CRlike_allSS");                  hbins.push_back(16);hlow.push_back( 10);hup.push_back(170);
  histonames.push_back("Mll_CRloose_allSS");                 hbins.push_back(16);hlow.push_back( 10);hup.push_back(170);
  histonames.push_back("MET_CRlike_allSS");                  hbins.push_back(14);hlow.push_back(  0);hup.push_back(140);
  histonames.push_back("MET_CRloose_allSS");                 hbins.push_back(14);hlow.push_back(  0);hup.push_back(140);
  histonames.push_back("Detajj_CRlike_allSS");               hbins.push_back(15);hlow.push_back(  0);hup.push_back(4.5);
  histonames.push_back("Detajj_CRloose_allSS");              hbins.push_back(15);hlow.push_back(  0);hup.push_back(4.5);
  histonames.push_back("MTmax_CRlike_allSS");                hbins.push_back(14);hlow.push_back(  0);hup.push_back(210);
  histonames.push_back("MTmax_CRloose_allSS");               hbins.push_back(14);hlow.push_back(  0);hup.push_back(210);
  histonames.push_back("MET_CRlike_allSFOS");                hbins.push_back(14);hlow.push_back(  0);hup.push_back(140);
  histonames.push_back("MET_CRloose_allSFOS");               hbins.push_back(14);hlow.push_back(  0);hup.push_back(140);
  histonames.push_back("MSFOS_CRlike_allSFOS");              hbins.push_back(16);hlow.push_back( 10);hup.push_back(170);
  histonames.push_back("MSFOS_CRloose_allSFOS");             hbins.push_back(16);hlow.push_back( 10);hup.push_back(170);
  histonames.push_back("pTlll_CRlike_allSFOS");              hbins.push_back(12);hlow.push_back( 10);hup.push_back(120);
  histonames.push_back("pTlll_CRloose_allSFOS");             hbins.push_back(12);hlow.push_back( 10);hup.push_back(120);
  histonames.push_back("dPhiMETlll_CRlike_allSFOS");         hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);
  histonames.push_back("dPhiMETlll_CRloose_allSFOS");        hbins.push_back(16);hlow.push_back(0.1);hup.push_back(3.3);

  histonames.push_back("NJ_mumumu_CRlike_allSFOS");     hbins.push_back( 7);hlow.push_back(0.1);hup.push_back(7);
  histonames.push_back("NJ_eee_CRlike_allSFOS");        hbins.push_back( 7);hlow.push_back(0.1);hup.push_back(7);
  histonames.push_back("NJ_noteeemmm_CRlike_allSFOS");  hbins.push_back( 7);hlow.push_back(0.1);hup.push_back(7);

  //after finding good Mll region test with 2l OS region the in/out ratio
  
  //ee,em,mm,0SFOS,1SFOS,2SFOS

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
      LorentzVector METx;   METx  .SetPxPyPzE(met_T1CHS_miniAOD_CORE_pt()   *TMath::Cos(met_T1CHS_miniAOD_CORE_phi()   ),met_T1CHS_miniAOD_CORE_pt()   *TMath::Sin(met_T1CHS_miniAOD_CORE_phi()   ),0,met_T1CHS_miniAOD_CORE_pt()   );
      LorentzVector MET_up; MET_up.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*TMath::Cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*TMath::Sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
      LorentzVector MET_dn; MET_dn.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*TMath::Cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*TMath::Sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
      //if(fabs(MET.Pt()-METx.Pt())>0.01*MET.Pt()) cout << "Damn " << MET.Pt() << " " << METx.Pt() << endl;

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

      int nj_up(0),nb_up(0);
      int nj20_up(0),nj30_up(0);
      vector<int> i2p5_up;
      for(unsigned int n = 0; n<jets_up_csv().size();++n){
	if(fabs(jets_up_p4()[n].Eta())<2.5) i2p5_up.push_back(n);
	if(jets_up_p4()[n].Pt()>20&&fabs(jets_up_p4()[n].Eta())<5) ++nj_up;
	if(jets_up_p4()[n].Pt()>30&&fabs(jets_up_p4()[n].Eta())<2.5) ++nj30_up;
	if(jets_up_p4()[n].Pt()>20&&fabs(jets_up_p4()[n].Eta())<2.5) ++nj20_up;
	if(jets_up_csv()[n]>0.5426&&jets_up_p4()[n].Pt()>20&&fabs(jets_up_p4()[n].Eta())<2.4) ++nb_up;
      }
      int nj_dn(0),nb_dn(0);
      int nj20_dn(0),nj30_dn(0);
      vector<int> i2p5_dn;
      for(unsigned int n = 0; n<jets_dn_csv().size();++n){
	if(fabs(jets_dn_p4()[n].Eta())<2.5) i2p5_dn.push_back(n);
	if(jets_dn_p4()[n].Pt()>20&&fabs(jets_dn_p4()[n].Eta())<5) ++nj_dn;
	if(jets_dn_p4()[n].Pt()>30&&fabs(jets_dn_p4()[n].Eta())<2.5) ++nj30_dn;
	if(jets_dn_p4()[n].Pt()>20&&fabs(jets_dn_p4()[n].Eta())<2.5) ++nj20_dn;
	if(jets_dn_csv()[n]>0.5426&&jets_dn_p4()[n].Pt()>20&&fabs(jets_dn_p4()[n].Eta())<2.4) ++nb_dn;
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

      minDR=999.;
      jDR1 = -1; jDR2 = -1;
      for(unsigned int j1 = 0; j1<i2p5_up.size();++j1){
	if(jets_up_p4()[i2p5_up[j1] ].Pt()<30) continue;
	for(unsigned int j2 = j1+1; j2<i2p5_up.size();++j2){
	  if(jets_up_p4()[i2p5_up[j2] ].Pt()<30) continue;
	  if(dR(jets_up_p4()[i2p5_up[j1] ], jets_up_p4()[i2p5_up[j2] ])<minDR){
	    minDR = dR(jets_up_p4()[i2p5_up[j1] ], jets_up_p4()[i2p5_up[j2] ]);
	    jDR1 = i2p5_up[j1]; jDR2 = i2p5_up[j2];
	  }
	}
      }
      float Mjj_up = -1;
      if(jDR1>=0&&jDR2>=0) Mjj_up = (jets_up_p4()[jDR1]+jets_up_p4()[jDR2]).M();
      float MjjL_up = -1; float Detajj_up = -1;
      if(i2p5_up.size()>1&&jets_up_p4()[i2p5_up[0] ].Pt()>30&&jets_up_p4()[i2p5_up[1] ].Pt()>30) {
	MjjL_up = (jets_up_p4()[i2p5_up[0] ]+jets_up_p4()[i2p5_up[1] ]).M();
	Detajj_up = dEta(jets_up_p4()[i2p5_up[0] ],jets_up_p4()[i2p5_up[1] ]);
      }
      minDR=999.;
      jDR1 = -1; jDR2 = -1;
      for(unsigned int j1 = 0; j1<i2p5_dn.size();++j1){
	if(jets_dn_p4()[i2p5_dn[j1] ].Pt()<30) continue;
	for(unsigned int j2 = j1+1; j2<i2p5_dn.size();++j2){
	  if(jets_dn_p4()[i2p5_dn[j2] ].Pt()<30) continue;
	  if(dR(jets_dn_p4()[i2p5_dn[j1] ], jets_dn_p4()[i2p5_dn[j2] ])<minDR){
	    minDR = dR(jets_dn_p4()[i2p5_dn[j1] ], jets_dn_p4()[i2p5_dn[j2] ]);
	    jDR1 = i2p5_dn[j1]; jDR2 = i2p5_dn[j2];
	  }
	}
      }
      float Mjj_dn = -1;
      if(jDR1>=0&&jDR2>=0) Mjj_dn = (jets_dn_p4()[jDR1]+jets_dn_p4()[jDR2]).M();
      float MjjL_dn = -1; float Detajj_dn = -1;
      if(i2p5_dn.size()>1&&jets_dn_p4()[i2p5_dn[0] ].Pt()>30&&jets_dn_p4()[i2p5_dn[1] ].Pt()>30) {
	MjjL_dn = (jets_dn_p4()[i2p5_dn[0] ]+jets_dn_p4()[i2p5_dn[1] ]).M();
	Detajj_dn = dEta(jets_dn_p4()[i2p5_dn[0] ],jets_dn_p4()[i2p5_dn[1] ]);
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
      if(fabs(MjjL)>400.)   passDetajj = false;
      bool passDetajj2 = true;
      if(nj30<2)            passDetajj2 = false;
      if(fabs(Detajj)>1.5)  passDetajj2 = false;

      bool passMDetajj_up = true;
      if(nj30_up<2)            passMDetajj_up = false;
      if(fabs(Detajj_up)>1.5)  passMDetajj_up = false;
      if(fabs(MjjL_up)>400.)   passMDetajj_up = false;
      if(fabs(Mjj_up-80.)>20.) passMDetajj_up = false;
      bool passMjj_up = true;
      if(nj30_up<2)            passMjj_up = false;
      if(fabs(Mjj_up-80.)>20.) passMjj_up = false;
      bool passDetajj_up = true;
      if(nj30_up<2)            passDetajj_up = false;
      if(fabs(Detajj_up)>1.5)  passDetajj_up = false;
      if(fabs(MjjL_up)>400.)   passDetajj_up = false;
      bool passDetajj2_up = true;
      if(nj30_up<2)            passDetajj2_up = false;
      if(fabs(Detajj_up)>1.5)  passDetajj2_up = false;
      bool passMDetajj_dn = true;
      if(nj30_dn<2)            passMDetajj_dn = false;
      if(fabs(Detajj_dn)>1.5)  passMDetajj_dn = false;
      if(fabs(MjjL_dn)>400.)   passMDetajj_dn = false;
      if(fabs(Mjj_dn-80.)>20.) passMDetajj_dn = false;
      bool passMjj_dn = true;
      if(nj30_dn<2)            passMjj_dn = false;
      if(fabs(Mjj_dn-80.)>20.) passMjj_dn = false;
      bool passDetajj_dn = true;
      if(nj30_dn<2)            passDetajj_dn = false;
      if(fabs(Detajj_dn)>1.5)  passDetajj_dn = false;
      if(fabs(MjjL_dn)>400.)   passDetajj_dn = false;
      bool passDetajj2_dn = true;
      if(nj30_dn<2)            passDetajj2_dn = false;
      if(fabs(Detajj_dn)>1.5)  passDetajj2_dn = false;
      vector<int> vSS, v3l, v, iSS, i3l;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool isSS = false; bool is3l = false;
	bool isaSS = false; bool isa3l = false;
	bool islSS = false; bool isl3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015) {
	  if(abs(lep_pdgId()[i])==11){
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
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_p4()[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4) {
	  v.push_back(i);
	  if(!isSS) vSS.push_back(i);
	  if(!is3l) v3l.push_back(i);
	}
      }
      float lepSFSS = 1.; float lepSFerrSS = 0.;
      float lepSF3l = 1.; float lepSFerr3l = 0.;
      if(applylepSF&&!isData()&&iSS.size()>=2){
	vector<float> eff1, err1, eff2, err2;
	int bin;
	for(unsigned i = 0; i<iSS.size(); ++i){
	  if(abs(lep_pdgId()[iSS[i] ])==11){
	    bin = hElReco->FindBin(std::max(eletamin,std::min(eletamax,lep_etaSC()[iSS[i] ])), std::max(elptminReco,std::min(elptmax,lep_p4()[iSS[i] ].Pt())));
	    eff2.push_back(hElReco->GetBinContent(bin));
	    err2.push_back(hElReco->GetBinError(bin));
	    bin = hElID->FindBin(std::max(eletamin,std::min(eletamax,lep_etaSC()[iSS[i] ])), std::max(elptmin,std::min(elptmax,lep_p4()[iSS[i] ].Pt())));
	    eff1.push_back(hElID->GetBinContent(bin));
	    err1.push_back(hElID->GetBinError(bin));
	  } else if(abs(lep_pdgId()[iSS[i] ])==13){
	    bin = hMu->FindBin(std::max(muetamin,std::min(muetamax,lep_etaSC()[iSS[i] ])), std::max(muptmin,std::min(muptmax,lep_p4()[iSS[i] ].Pt())));
	    eff1.push_back(hMu->GetBinContent(bin));
	    err1.push_back(hMu->GetBinError(bin));
	  }
	}
	for(unsigned int i = 0; i<eff2.size();++i) lepSFSS *= eff2[i];
	for(unsigned int i = 0; i<eff1.size();++i) lepSFSS *= eff1[i];
	for(unsigned int i = 0; i<eff2.size();++i) lepSFerrSS += pow(lepSFSS*err2[i]/eff2[i],2);
	for(unsigned int i = 0; i<eff1.size();++i) lepSFerrSS += pow(lepSFSS*err1[i]/eff1[i],2);
	lepSFerrSS = sqrt(lepSFerrSS);
	eff1.clear(); err1.clear(); eff2.clear(); err2.clear();
      }
      if(applylepSF&&!isData()&&i3l.size()>=3){
	vector<float> eff1, err1, eff2, err2;
	int bin;
	for(unsigned i = 0; i<i3l.size(); ++i){
	  if(abs(lep_pdgId()[i3l[i] ])==11){
	    bin = hElReco->FindBin(std::max(eletamin,std::min(eletamax,lep_etaSC()[i3l[i] ])), std::max(elptminReco,std::min(elptmax,lep_p4()[i3l[i] ].Pt())));
	    eff2.push_back(hElReco->GetBinContent(bin));
	    err2.push_back(hElReco->GetBinError(bin));
	    bin = hElID->FindBin(std::max(eletamin,std::min(eletamax,lep_etaSC()[i3l[i] ])), std::max(elptmin,std::min(elptmax,lep_p4()[i3l[i] ].Pt())));
	    eff1.push_back(hElID->GetBinContent(bin));
	    err1.push_back(hElID->GetBinError(bin));
	  } else if(abs(lep_pdgId()[i3l[i] ])==13){
	    bin = hMu->FindBin(std::max(muetamin,std::min(muetamax,lep_etaSC()[i3l[i] ])), std::max(muptmin,std::min(muptmax,lep_p4()[i3l[i] ].Pt())));
	    eff1.push_back(hMu->GetBinContent(bin));
	    err1.push_back(hMu->GetBinError(bin));
	  }
	}
	for(unsigned int i = 0; i<eff2.size();++i) lepSF3l *= eff2[i];
	for(unsigned int i = 0; i<eff1.size();++i) lepSF3l *= eff1[i];
	for(unsigned int i = 0; i<eff2.size();++i) lepSFerr3l += pow(lepSF3l*err2[i]/eff2[i],2);
	for(unsigned int i = 0; i<eff1.size();++i) lepSFerr3l += pow(lepSF3l*err1[i]/eff1[i],2);
	lepSFerr3l = sqrt(lepSFerr3l);
	eff1.clear(); err1.clear(); eff2.clear(); err2.clear();
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
      if(lep_p4()[i3l[0] ].Pt()<25) continue;
      if(!passofflineforTrigger) continue;
      
      if(isData()){
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) continue;
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = false;
	if(nmu>=2&&HLT_DoubleMu() )                               passonlineTrigger = true;
	if(nmu25>=1&&nel>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nmu>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) passonlineTrigger = true;
	//if(nmu>=2&&(HLT_DoubleMu_noiso()) )                             passonlineTrigger = true;
	//if(nmu25>=1&&nel>=1&&(HLT_MuEG_noiso()||HLT_MuEG_noiso_2()))    passonlineTrigger = true;
	//if(nel25>=1&&nel>=2&&(HLT_DoubleEl_noiso()))                    passonlineTrigger = true;
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
	if(iSS.size()>=2){
	  int l1(-1),l2(-1);
	  if(iSS.size()>=2) { l1 = iSS[0]; l2 = iSS[1]; }
	  int gentype = gentype_v2(l1,l2-1);
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
	  
	  if(sn2=="others"){
	    int nW(0), nZ(0), nG(0), nF(0);
	    if(lep_isFromW()[l1]) ++nW;
	    else if(lep_isFromZ()[l1]) ++nZ;
	    else if(lep_isFromB()[l1]||lep_isFromC()[l1]||lep_isFromL()[l1]||lep_isFromLF()[l1]) ++nF;
	    else if(lep_motherIdSS()[l1]==(-3)) ++nG;
	    if(lep_isFromW()[l2]) ++nW;
	    else if(lep_isFromZ()[l2]) ++nZ;
	    else if(lep_isFromB()[l2]||lep_isFromC()[l2]||lep_isFromL()[l2]||lep_isFromLF()[l2]) ++nF;
	    else if(lep_motherIdSS()[l2]==(-3)) ++nG;
	    if(lep_isFromW()[l3]) ++nW;
	    else if(lep_isFromZ()[l3]) ++nZ;
	    else if(lep_isFromB()[l3]||lep_isFromC()[l3]||lep_isFromL()[l3]||lep_isFromLF()[l3]) ++nF;
	    else if(lep_motherIdSS()[l3]==(-3)) ++nG;
	    if(nW==3&&(lep_mc_Id()[l1]>0&&lep_mc_Id()[l2]>0&&lep_mc_Id()[l3]>0)) sn2 = "chargeflips";
	    else if(nW==3&&(lep_mc_Id()[l1]<0&&lep_mc_Id()[l2]<0&&lep_mc_Id()[l3]<0)) sn2 = "chargeflips";
	    else if(nW==3) sn2 = "trueWWW";
	    else if(nW==2&&nZ==1) sn2 = "3lLL";//ttZ w/ LL
	    else if(nW==1&&nZ==2) sn2 = "true3L";//WZ, neglect WZZ as LL
	    else if(nZ==3) sn2 = "3lLL";//ZZ
	    else if((nW+nZ)<3&&nF>0) sn2 = "fakes";
	    else if((nW+nZ)<3&&nG>0) sn2 = "photonfakes";
	    else if((nW+nZ)<3&&(nW+nZ)>=1) sn2 = "fakes";//count those as fakes
	    else sn2 = "others";
	  
	    cout << " nW/nZ/nG/nF " << nW << "/" << nZ << "/" << nF << "/" << nG << " fname " << fname << endl;
	  }
	}
      }
      
      float MTmax = -1;
      if(iSS.size()>=2){
	if(mT(lep_p4()[iSS[1] ],MET)>mT(lep_p4()[iSS[0] ],MET)) MTmax = mT(lep_p4()[iSS[1] ],MET);
	else MTmax = mT(lep_p4()[iSS[0] ],MET);
      }
      float MTmax3l = -1;
      if(i3l.size()>=3){
	if((mT(lep_p4()[i3l[1] ],MET)>mT(lep_p4()[i3l[0] ],MET))&&(mT(lep_p4()[i3l[1] ],MET)>mT(lep_p4()[i3l[2] ],MET))) MTmax3l = mT(lep_p4()[i3l[1] ],MET);
	else if(mT(lep_p4()[i3l[2] ],MET)>mT(lep_p4()[i3l[0] ],MET)) MTmax3l = mT(lep_p4()[i3l[2] ],MET);
	else MTmax3l = mT(lep_p4()[i3l[0] ],MET);
      }

      float MTmax_up = -1;
      if(iSS.size()>=2){
	if(mT(lep_p4()[iSS[1] ],MET_up)>mT(lep_p4()[iSS[0] ],MET_up)) MTmax_up = mT(lep_p4()[iSS[1] ],MET_up);
	else MTmax_up = mT(lep_p4()[iSS[0] ],MET_up);
      }
      float MTmax3l_up = -1;
      if(i3l.size()>=3){
	if((mT(lep_p4()[i3l[1] ],MET_up)>mT(lep_p4()[i3l[0] ],MET_up))&&(mT(lep_p4()[i3l[1] ],MET_up)>mT(lep_p4()[i3l[2] ],MET_up))) MTmax3l_up = mT(lep_p4()[i3l[1] ],MET_up);
	else if(mT(lep_p4()[i3l[2] ],MET_up)>mT(lep_p4()[i3l[0] ],MET_up)) MTmax3l_up = mT(lep_p4()[i3l[2] ],MET_up);
	else MTmax3l_up = mT(lep_p4()[i3l[0] ],MET_up);
      }
      float MTmax_dn = -1;
      if(iSS.size()>=2){
	if(mT(lep_p4()[iSS[1] ],MET_dn)>mT(lep_p4()[iSS[0] ],MET_dn)) MTmax_dn = mT(lep_p4()[iSS[1] ],MET_dn);
	else MTmax_dn = mT(lep_p4()[iSS[0] ],MET_dn);
      }
      float MTmax3l_dn = -1;
      if(i3l.size()>=3){
	if((mT(lep_p4()[i3l[1] ],MET_dn)>mT(lep_p4()[i3l[0] ],MET_dn))&&(mT(lep_p4()[i3l[1] ],MET_dn)>mT(lep_p4()[i3l[2] ],MET_dn))) MTmax3l_dn = mT(lep_p4()[i3l[1] ],MET_dn);
	else if(mT(lep_p4()[i3l[2] ],MET_dn)>mT(lep_p4()[i3l[0] ],MET_dn)) MTmax3l_dn = mT(lep_p4()[i3l[2] ],MET_dn);
	else MTmax3l_dn = mT(lep_p4()[i3l[0] ],MET_dn);
      }

      float MllCR(-1), MllCRb(-1), MllCRbv2(-1);
      bool ee[50],mm[50],em[50];
      for(int i = 0; i<50; ++i) { ee[i] = false; em[i] = false; mm[i] = false; }
      //0: SR, 5: SR but no Mjj cut
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                        mm[0] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[0] = false; mm[0] = false; }
	ee[5] = ee[0]; em[5] = em[0]; mm[5] = mm[0];
	if(!passMDetajj){ ee[0] = false; em[0] = false; mm[0] = false; }
      }
      //40: SR_jesup, 41: SR_jesup but no Mjj cut, 42: SR_jesup but no Mjj but jesup only for jets
      if(nj30_up>=2&&nb_up==0&&passDetajj_up&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&MET_up.Pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[30] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&MET_up.Pt()>40.&&MTmax_up>90.)                                            em[30] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&MET_up.Pt()>40.&&MTmax_up>90.)                                            em[30] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                           mm[30] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[30] = false; mm[30] = false; }
	ee[31] = ee[30]; em[31] = em[30]; mm[31] = mm[30];
	if(!passMDetajj_up){ ee[30] = false; em[30] = false; mm[30] = false; }
      }
      //32: SR_jesdn, 33: SR_jesdn but no Mjj cut
      if(nj30_dn>=2&&nb_dn==0&&passDetajj_dn&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&MET_dn.Pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[32] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&MET_dn.Pt()>40.&&MTmax_dn>90.)                                            em[32] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&MET_dn.Pt()>40.&&MTmax_dn>90.)                                            em[32] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                           mm[32] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[32] = false; mm[32] = false;  }
	ee[33] = ee[32]; em[33] = em[32]; mm[33] = mm[32];
	if(!passMDetajj_dn){ ee[32] = false; em[32] = false; mm[32] = false; }
      }
      //1: 3l, 6: 3l but no Mjj cut
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	int lep3 = -1;
	if(i3l[0]!=iSS[0]&&i3l[0]!=iSS[1]) lep3 = i3l[0];
	else if(i3l[1]!=iSS[0]&&i3l[1]!=iSS[1]) lep3 = i3l[1];
	else if(i3l[2]!=iSS[0]&&i3l[2]!=iSS[1]) lep3 = i3l[2];
	bool hasSFOSZ = false;
	if((lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[0] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if((lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[1] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	bool hasSFOS = false;
	double t1(-1), t2(-1);
	if(lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3])) { hasSFOS = true; t1 = (lep_p4()[iSS[0] ]+lep_p4()[lep3]).M(); }
	if(lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3])) { hasSFOS = true; t2 = (lep_p4()[iSS[1] ]+lep_p4()[lep3]).M(); }
	if(t1>0&&t2>0) MllCR = fabs(t1-90.)<fabs(t2-90) ? t1 : t2;
	else if(t1>0) MllCR = t1; else if(t2>0) MllCR = t2;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[1] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&MTmax3l>90.)                                             em[1] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()>40.&&MTmax3l>90.)                                             em[1] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                        mm[1] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[1] = false; mm[1] = false; }
	ee[6] = ee[1]; em[6] = em[1]; mm[6] = mm[1];
	if(!passMDetajj){ ee[1] = false; em[1] = false; mm[1] = false; }
	//14: no Mll cut, 15: no Mll cut no Mjj cut
	if(hasSFOS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[14] = true;
	if(hasSFOS&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&MTmax3l>90.)                                             em[14] = true;
	if(hasSFOS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()>40.&&MTmax3l>90.)                                             em[14] = true;
	if(hasSFOS&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                        mm[14] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[14] = false; mm[14] = false; }
	ee[15] = ee[14]; em[15] = em[14]; mm[15] = mm[14];
	if(!passMDetajj){ ee[14] = false; em[14] = false; mm[14] = false; }
      }
      if(ee[6]) ++cCRee;
      if(em[6]) ++cCRem;
      if(mm[6]) ++cCRmm;
      //34: 3l jesup, 35: 3l but no Mjj cut jesup
      if(nj30_up>=2&&nb_up==0&&passDetajj_up&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	int lep3 = -1;
	if(i3l[0]!=iSS[0]&&i3l[0]!=iSS[1]) lep3 = i3l[0];
	else if(i3l[1]!=iSS[0]&&i3l[1]!=iSS[1]) lep3 = i3l[1];
	else if(i3l[2]!=iSS[0]&&i3l[2]!=iSS[1]) lep3 = i3l[2];
	bool hasSFOSZ = false;
	if((lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[0] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if((lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[1] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&MET_up.Pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[34] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&MET_up.Pt()>40.&&MTmax3l_up>90.)                                          em[34] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&MET_up.Pt()>40.&&MTmax3l_up>90.)                                          em[34] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                           mm[34] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[34] = false; mm[34] = false;}
	ee[35] = ee[34]; em[35] = em[34]; mm[35] = mm[34];
	if(!passMDetajj_up){ ee[34] = false; em[34] = false; mm[34] = false; }
      }
      //36: 3l jesdn, 37: 3l but no Mjj cut jesdn
      if(nj30_dn>=2&&nb_dn==0&&passDetajj_dn&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	int lep3 = -1;
	if(i3l[0]!=iSS[0]&&i3l[0]!=iSS[1]) lep3 = i3l[0];
	else if(i3l[1]!=iSS[0]&&i3l[1]!=iSS[1]) lep3 = i3l[1];
	else if(i3l[2]!=iSS[0]&&i3l[2]!=iSS[1]) lep3 = i3l[2];
	bool hasSFOSZ = false;
	if((lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[0] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if((lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[1] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&MET_dn.Pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[36] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&MET_dn.Pt()>40.&&MTmax3l_dn>90.)                                          em[36] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&MET_dn.Pt()>40.&&MTmax3l_dn>90.)                                          em[36] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                           mm[36] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[36] = false; mm[36] = false; }
	ee[37] = ee[36]; em[37] = em[36]; mm[37] = mm[36];
	if(!passMDetajj_dn){ ee[36] = false; em[36] = false; mm[36] = false; }
      }
      //3: SR-like, low MET (i.e. no mumu), 8: SR-like, low MET (i.e. no mumu) but no Mjj cut
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()<40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()<40.&&MTmax3l>90.)                                             em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()<40.&&MTmax3l>90.)                                             em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()<0)                                                            mm[3] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[3] = false; mm[3] = false; }
	ee[8] = ee[3]; em[8] = em[3]; mm[8] = mm[3];
	if(!passMDetajj){ ee[3] = false; em[3] = false; mm[3] = false; }
      }
      //4: 3l, low MET (i.e. no mumu), 9: 3l, low MET (i.e. no mumu) but no Mjj cut
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	int lep3 = -1;
	if(i3l[0]!=iSS[0]&&i3l[0]!=iSS[1]) lep3 = i3l[0];
	else if(i3l[1]!=iSS[0]&&i3l[1]!=iSS[1]) lep3 = i3l[1];
	else if(i3l[2]!=iSS[0]&&i3l[2]!=iSS[1]) lep3 = i3l[2];
	bool hasSFOSZ = false;
	if((lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[0] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if((lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[1] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()<30.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[4] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()<30.&&MTmax3l>90.)                                             em[4] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()<30.&&MTmax3l>90.)                                             em[4] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()<0)                                                            mm[4] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[4] = false; mm[4] = false; }
	ee[9] = ee[4]; em[9] = em[4]; mm[9] = mm[4];
	if(!passMDetajj){ ee[4] = false; em[4] = false; mm[4] = false; }
      }
      //10: 3l but SS not leading leps, 11: as 10 but Mll cut on leading leptons, 12: as 10 but no Mjj cut, 13: as 11 but no Mjj cut
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3&&((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ])<0)/*&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.*/){
	bool isSS = false;
	int lep3 = -1; int SS1 = -1; int SS2 = -1;
	if((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[1]; SS1 = i3l[0]; SS2 = i3l[2]; isSS = true; }
	if((lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[0]; SS1 = i3l[1]; SS2 = i3l[2]; isSS = true; }
	if(isSS){
	  bool hasSFOSZ = false;
	  if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  bool hasSFOS = false;
	  double t1(-1), t2(-1);
	  if(lep_pdgId()[SS1]==(-lep_pdgId()[lep3])) { hasSFOS = true; t1 = (lep_p4()[SS1]+lep_p4()[lep3]).M(); }
	  if(lep_pdgId()[SS2]==(-lep_pdgId()[lep3])) { hasSFOS = true; t2 = (lep_p4()[SS2]+lep_p4()[lep3]).M(); }
	  if(t1>0&&t2>0) MllCRb = fabs(t1-90.)<fabs(t2-90) ? t1 : t2;
	  else if(t1>0) MllCRb = t1; else if(t2>0) MllCRb = t2;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==11&&met_pt()>40.&&fabs((lep_p4()[SS1]+lep_p4()[SS2]).M()-90.)>10.) ee[10] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==11&&met_pt()>40.&&MTmax3l>90.)                                     em[10] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==13&&met_pt()>40.&&MTmax3l>90.)                                     em[10] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==13)                                                                mm[10] = true;
	  ee[11] = ee[10]; em[11] = em[10]; mm[11] = mm[10];
	  if((lep_p4()[SS1]+lep_p4()[SS2]).M()<40){ ee[10] = false; mm[10] = false; }
	  if((lep_p4()[SS1]+lep_p4()[SS2]).M()<30){ em[10] = false; }
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<40){ ee[11] = false; mm[11] = false; }
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<30){ em[11] = false; }
	  ee[12] = ee[10]; em[12] = em[10]; mm[12] = mm[10];
	  ee[13] = ee[11]; em[13] = em[11]; mm[13] = mm[11];
	  if(!passMDetajj){ ee[10] = false; em[10] = false; mm[10] = false; ee[11] = false; em[11] = false; mm[11] = false; }
	  //16: no Mll cut, 17: no Mll cut no Mjj cut - SSpairnotlead2
	  //18: no Mll cut, 18: no Mll cut no Mjj cut - SSpairnotlead2v2
	  if(hasSFOS&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==11&&met_pt()>40.&&fabs((lep_p4()[SS1]+lep_p4()[SS2]).M()-90.)>10.) ee[16] = true;
	  if(hasSFOS&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==11&&met_pt()>40.&&MTmax3l>90.)                                     em[16] = true;
	  if(hasSFOS&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==13&&met_pt()>40.&&MTmax3l>90.)                                     em[16] = true;
	  if(hasSFOS&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==13)                                                                mm[16] = true;
	  ee[18] = ee[16]; em[18] = em[16]; mm[18] = mm[16];
	  if((lep_p4()[SS1]+lep_p4()[SS2]).M()<40){ ee[16] = false; mm[16] = false; }
	  if((lep_p4()[SS1]+lep_p4()[SS2]).M()<30){ em[16] = false; }
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<40){ ee[18] = false; mm[18] = false; }
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<30){ em[18] = false; }
	  ee[17] = ee[16]; em[17] = em[16]; mm[17] = mm[16];
	  ee[19] = ee[18]; em[19] = em[18]; mm[19] = mm[18];
	  if(!passMDetajj){ ee[16] = false; em[16] = false; mm[16] = false; ee[18] = false; em[18] = false; mm[18] = false; }
	}
      }
      if(ee[12]) ++cCR2ee;
      if(em[12]) ++cCR2em;
      if(mm[12]) ++cCR2mm;
      if(ee[13]) ++cCR3ee;
      if(em[13]) ++cCR3em;
      if(mm[13]) ++cCR3mm;
      //44: CR loose, but 2 leading leps are not SS
      if(nj30>=2&&nb==0&&nSS>=2&&n3l>=3&&((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ])<0)){
	bool isSS = false;
	int lep3 = -1; int SS1 = -1; int SS2 = -1;
	if((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[1]; SS1 = i3l[0]; SS2 = i3l[2]; isSS = true; }
	if((lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[0]; SS1 = i3l[1]; SS2 = i3l[2]; isSS = true; }
	if(isSS){
	  bool hasSFOSZ = false;
	  if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==11) ee[44] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==11) em[44] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==13) em[44] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==13) mm[44] = true;
	}
      }
      //20: SS+2j+eq3l(incl. isotrack) with MZ, 22: SS+2j+eq3l(excl. isotrack) with MZ, 24: SS+2j+ge3l with MZ
      //21: SS+2j+eq3l(incl. isotrack)  w/o MZ, 23: SS+2j+eq3l(excl. isotrack)  w/o MZ, 25: SS+2j+ge3l  w/o MZ
      bool passMll = false;//little helper
      if(nj30>=2&&nb==0&&nSS>=2&&n3l>=3&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	int lep3 = -1;
	if(i3l[0]!=iSS[0]&&i3l[0]!=iSS[1]) lep3 = i3l[0];
	else if(i3l[1]!=iSS[0]&&i3l[1]!=iSS[1]) lep3 = i3l[1];
	else if(i3l[2]!=iSS[0]&&i3l[2]!=iSS[1]) lep3 = i3l[2];
	bool hasSFOSZ = false;
	if((lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[0] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if((lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[1] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	bool hasSFOS = false;
	if(lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3])) hasSFOS = true;
	if(lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3])) hasSFOS = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) ee[20] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) ee[21] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&                    nveto3l==0&&n3l==3) ee[22] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&                    nveto3l==0&&n3l==3) ee[23] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11                                        ) ee[24] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11                                        ) ee[25] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) em[20] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) em[21] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&                    nveto3l==0&&n3l==3) em[22] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&                    nveto3l==0&&n3l==3) em[23] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13                                        ) em[24] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13                                        ) em[25] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) em[20] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) em[21] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&                    nveto3l==0&&n3l==3) em[22] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&                    nveto3l==0&&n3l==3) em[23] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11                                        ) em[24] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11                                        ) em[25] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) mm[20] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&n3l==3) mm[21] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13&&                    nveto3l==0&&n3l==3) mm[22] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13&&                    nveto3l==0&&n3l==3) mm[23] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13                                        ) mm[24] = true;
	if(hasSFOS &&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13                                        ) mm[25] = true;
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11){
	  if(((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40)&&(fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.)) passMll  = true;
	}
	if(hasSFOSZ&&((abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11)||(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13))){
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30) passMll  = true;
	}
	if(hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13){
	  if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40) passMll  = true;
	}
      }
      //38: CR loose jesup, 39: CR loose jesdn
      if(nSS>=2&&n3l>=3&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	int lep3 = -1;
	if(i3l[0]!=iSS[0]&&i3l[0]!=iSS[1]) lep3 = i3l[0];
	else if(i3l[1]!=iSS[0]&&i3l[1]!=iSS[1]) lep3 = i3l[1];
	else if(i3l[2]!=iSS[0]&&i3l[2]!=iSS[1]) lep3 = i3l[2];
	bool hasSFOSZ = false;
	if((lep_pdgId()[iSS[0] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[0] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if((lep_pdgId()[iSS[1] ]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[iSS[1] ]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	if(nj30_up>=2&&nb_up==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) ee[38] = true;
	if(nj30_up>=2&&nb_up==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13) em[38] = true;
	if(nj30_up>=2&&nb_up==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11) em[38] = true;
	if(nj30_up>=2&&nb_up==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) mm[38] = true;
	if(nj30_dn>=2&&nb_dn==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) ee[39] = true;
	if(nj30_dn>=2&&nb_dn==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13) em[39] = true;
	if(nj30_dn>=2&&nb_dn==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11) em[39] = true;
	if(nj30_dn>=2&&nb_dn==0&&hasSFOSZ&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) mm[39] = true;
      }
      //40: 3l CR but SS not 2 leading leps, jesdown; 41: as 40 but no Mjj
      if(nj30_dn>=2&&nb_dn==0&&passDetajj_dn&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3&&((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ])<0)){
	bool isSS = false;
	int lep3 = -1; int SS1 = -1; int SS2 = -1;
	if((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[1]; SS1 = i3l[0]; SS2 = i3l[2]; isSS = true; }
	if((lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[0]; SS1 = i3l[1]; SS2 = i3l[2]; isSS = true; }
	if(isSS){
	  bool hasSFOSZ = false;
	  if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==11&&MET_dn.Pt()>40.&&fabs((lep_p4()[SS1]+lep_p4()[SS2]).M()-90.)>10.) ee[40] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==11&&MET_dn.Pt()>40.&&MTmax3l_dn>90.)                                  em[40] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==13&&MET_dn.Pt()>40.&&MTmax3l_dn>90.)                                  em[40] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==13)                                                                   mm[40] = true;
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<40){ ee[40] = false; mm[40] = false; }
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<30){ em[40] = false; }
	  ee[41] = ee[40]; em[41] = em[40]; mm[41] = mm[40];
	  if(!passMDetajj_dn){ ee[40] = false; em[40] = false; mm[40] = false; }
	}
      }
      //42: 3l CR but SS not 2 leading leps, jesup; 43: as 42 but no Mjj
      if(nj30_up>=2&&nb_up==0&&passDetajj_up&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3&&((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ])<0)){
	bool isSS = false;
	int lep3 = -1; int SS1 = -1; int SS2 = -1;
	if((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[1]; SS1 = i3l[0]; SS2 = i3l[2]; isSS = true; }
	if((lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[0]; SS1 = i3l[1]; SS2 = i3l[2]; isSS = true; }
	if(isSS){
	  bool hasSFOSZ = false;
	  if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==11&&MET_up.Pt()>40.&&fabs((lep_p4()[SS1]+lep_p4()[SS2]).M()-90.)>10.) ee[42] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==11&&MET_up.Pt()>40.&&MTmax3l_up>90.)                                  em[42] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==13&&MET_up.Pt()>40.&&MTmax3l_up>90.)                                  em[42] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==13)                                                                   mm[42] = true;
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<40){ ee[42] = false; mm[42] = false; }
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<30){ em[42] = false; }
	  ee[43] = ee[42]; em[43] = em[42]; mm[43] = mm[42];
	  if(!passMDetajj_up){ ee[42] = false; em[42] = false; mm[42] = false; }
	}
      }
      //45: SR preselection
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) { ee[45] = true; }
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) { ee[45] = true; }
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11) { em[45] = true; }
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13) { em[45] = true; }
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) { mm[45] = true; }
      }
      if(nj30>=2&&nb==0&&passDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nveto3l==0&&nSS>=2&&n3l==3 ){
	//SS is just tighter pt, so can loop over i3l
	bool isSS = false;
	int lep3 = -1; int SS1 = -1; int SS2 = -1;
	if((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ])>0) { lep3 = i3l[2]; SS1 = i3l[0]; SS2 = i3l[1]; isSS = true; }
	if((lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[1]; SS1 = i3l[0]; SS2 = i3l[2]; isSS = true; }
	if((lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ])>0) { lep3 = i3l[0]; SS1 = i3l[1]; SS2 = i3l[2]; isSS = true; }
	if(isSS){
	  bool hasSFOSZ = false;
	  if((lep_pdgId()[SS1]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS1]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if((lep_pdgId()[SS2]==(-lep_pdgId()[lep3]))&&(fabs((lep_p4()[SS2]+lep_p4()[lep3]).M()-90.)<=10)) hasSFOSZ = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==11&&met_pt()>40.&&fabs((lep_p4()[SS1]+lep_p4()[SS2]).M()-90.)>10.) ee[46] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==11&&met_pt()>40.&&MTmax3l>90.)                                     em[46] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==11&&abs(lep_pdgId()[SS2])==13&&met_pt()>40.&&MTmax3l>90.)                                     em[46] = true;
	  if(hasSFOSZ&&abs(lep_pdgId()[SS1])==13&&abs(lep_pdgId()[SS2])==13)                                                                mm[46] = true;
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<40){ ee[46] = false; mm[46] = false; }
	  if((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<30){ em[46] = false; }
	}
      }

      int SFOS[50];
      double pTlll(-1), DPhilllMET(-1);
      if(i3l.size()>=3){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
      }
      for(int i = 0; i<50; ++i) { SFOS[i] = -1; }
      //0 : SR, 1: inverted Mll-Z, 2: SR but inverted DPhi,Pt, 3: inverted Mll-Z and inverted DPhi,Pt
      if(nj<2&&nb==0/*&&nveto3l==0*/&&n3l==3&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)>10.){
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
	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0X = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
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
	  if(met_pt()<45) pass1=false;
	  if(met_pt()<45) pass1X=false;

	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  if(met_pt()<55) pass2=false;
	  if(met_pt()<55) pass2X=false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	//SR
	if(DPhilllMET>2.7&&pTlll>60&&pass0) SFOS[0] = 0;
	if(DPhilllMET>2.5&&pTlll>60&&pass1) SFOS[0] = 1;
	if(DPhilllMET>2.5&&pTlll>60&&pass2) SFOS[0] = 2;
	//CR
	if(DPhilllMET>2.7&&pTlll>60&&pass0X) SFOS[1] = 0;
	if(DPhilllMET>2.5&&pTlll>60&&pass1X) SFOS[1] = 1;
	if(DPhilllMET>2.5&&pTlll>60&&pass2X) SFOS[1] = 2;
	//SR inverted DPhi,Pt
	if(DPhilllMET<2.7&&pTlll<60&&pass0) SFOS[2] = 0;
	if(DPhilllMET<2.5&&pTlll<60&&pass1) SFOS[2] = 1;
	if(DPhilllMET<2.5&&pTlll<60&&pass2) SFOS[2] = 2;
	//CR inverted DPhi,Pt
	if(DPhilllMET<2.7&&pTlll<60&&pass0X) SFOS[3] = 0;
	if(DPhilllMET<2.5&&pTlll<60&&pass1X) SFOS[3] = 1;
	if(DPhilllMET<2.5&&pTlll<60&&pass2X) SFOS[3] = 2;
      }
      //4-7: as 0-3 but low MET
      if(nj<2&&nb==0/*&&nveto3l==0*/&&n3l==3&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)>10.){
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
	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(met_pt()>0) pass0=false;
	  if(met_pt()>0) pass0X=false;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0X = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
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
	  if(met_pt()>45) pass1=false;
	  if(met_pt()>45) pass1X=false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  if(met_pt()>55) pass2=false;
	  if(met_pt()>55) pass2X=false;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	//SR
	if(DPhilllMET>2.7&&pTlll>60&&pass0) SFOS[4] = 0;
	if(DPhilllMET>2.5&&pTlll>60&&pass1) SFOS[4] = 1;
	if(DPhilllMET>2.5&&pTlll>60&&pass2) SFOS[4] = 2;
	//CR
	if(DPhilllMET>2.7&&pTlll>60&&pass0X) SFOS[5] = 0;
	if(DPhilllMET>2.5&&pTlll>60&&pass1X) SFOS[5] = 1;
	if(DPhilllMET>2.5&&pTlll>60&&pass2X) SFOS[5] = 2;
	//SR inverted DPhi,Pt
	if(DPhilllMET<2.7&&pTlll<60&&pass0) SFOS[6] = 0;
	if(DPhilllMET<2.5&&pTlll<60&&pass1) SFOS[6] = 1;
	if(DPhilllMET<2.5&&pTlll<60&&pass2) SFOS[6] = 2;
	//CR inverted DPhi,Pt
	if(DPhilllMET<2.7&&pTlll<60&&pass0X) SFOS[7] = 0;
	if(DPhilllMET<2.5&&pTlll<60&&pass1X) SFOS[7] = 1;
	if(DPhilllMET<2.5&&pTlll<60&&pass2X) SFOS[7] = 2;
      }
      //8: SR but >=2j, 9: CR but >=2j, 10: SR but invert either Pt,DPhi,MET, 11: CR but invert either Pt,DPhi,MET
      if(/*nj<2&&*/nb==0/*&&nveto3l==0*/&&n3l==3&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)>10.){
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
	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0X = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass0X = false;
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
	//SR but ge2j
	if(nj>=2&&DPhilllMET>2.7&&pTlll>60&&pass0) SFOS[8] = 0;
	if(nj>=2&&DPhilllMET>2.5&&pTlll>60&&pass1) SFOS[8] = 1;
	if(nj>=2&&DPhilllMET>2.5&&pTlll>60&&pass2) SFOS[8] = 2;
	//SR inverted DPhi,Pt, OR MET
	if(nj<2&&pass0){
	  if(DPhilllMET<2.7&&pTlll>60) SFOS[10] = 0;//invert DPhi
	  if(DPhilllMET>2.7&&pTlll<60) SFOS[10] = 0;//invert Pt
	  if(DPhilllMET<2.7&&pTlll<60) SFOS[10] = 0;//invert DPhi/Pt
	}
	if(nj<2&&pass1){
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()>45) SFOS[10] = 1;//invert DPhi
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()>45) SFOS[10] = 1;//invert Pt
	  if(DPhilllMET>2.5&&pTlll>60&&met_pt()<45) SFOS[10] = 1;//invert MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()>45) SFOS[10] = 1;//invert DPhi/Pt
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()<45) SFOS[10] = 1;//invert DPhi/MET
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()<45) SFOS[10] = 1;//invert Pt/MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()<45) SFOS[10] = 1;//invert all3
	}
	if(nj<2&&pass2){
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()>55) SFOS[10] = 2;//invert DPhi
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()>55) SFOS[10] = 2;//invert Pt
	  if(DPhilllMET>2.5&&pTlll>60&&met_pt()<55) SFOS[10] = 2;//invert MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()>55) SFOS[10] = 2;//invert DPhi/Pt
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()<55) SFOS[10] = 2;//invert DPhi/MET
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()<55) SFOS[10] = 2;//invert Pt/MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()<55) SFOS[10] = 2;//invert all3
	}
	//CR but ge2j
	if(nj>=2&&DPhilllMET>2.7&&pTlll>60&&pass0X) SFOS[9] = 0;
	if(nj>=2&&DPhilllMET>2.5&&pTlll>60&&pass1X) SFOS[9] = 1;
	if(nj>=2&&DPhilllMET>2.5&&pTlll>60&&pass2X) SFOS[9] = 2;
	//CR inverted DPhi,Pt, OR MET
	if(nj<2&&pass0X){
	  if(DPhilllMET<2.7&&pTlll>60) SFOS[11] = 0;//invert DPhi
	  if(DPhilllMET>2.7&&pTlll<60) SFOS[11] = 0;//invert Pt
	  if(DPhilllMET<2.7&&pTlll<60) SFOS[11] = 0;//invert DPhi/Pt
	}
	if(nj<2&&pass1X){
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()>45) SFOS[11] = 1;//invert DPhi
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()>45) SFOS[11] = 1;//invert Pt
	  if(DPhilllMET>2.5&&pTlll>60&&met_pt()<45) SFOS[11] = 1;//invert MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()>45) SFOS[11] = 1;//invert DPhi/Pt
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()<45) SFOS[11] = 1;//invert DPhi/MET
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()<45) SFOS[11] = 1;//invert Pt/MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()<45) SFOS[11] = 1;//invert all3
	}
	if(nj<2&&pass2X){
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()>55) SFOS[11] = 2;//invert DPhi
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()>55) SFOS[11] = 2;//invert Pt
	  if(DPhilllMET>2.5&&pTlll>60&&met_pt()<55) SFOS[11] = 2;//invert MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()>55) SFOS[11] = 2;//invert DPhi/Pt
	  if(DPhilllMET<2.5&&pTlll>60&&met_pt()<55) SFOS[11] = 2;//invert DPhi/MET
	  if(DPhilllMET>2.5&&pTlll<60&&met_pt()<55) SFOS[11] = 2;//invert Pt/MET
	  if(DPhilllMET<2.5&&pTlll<60&&met_pt()<55) SFOS[11] = 2;//invert all3
	}
	if(nj<2&&pass0&&DPhilllMET<2.7&&pTlll>60)              SFOS[13] = 0;//invert DPhi
	if(nj<2&&pass1&&DPhilllMET<2.5&&pTlll>60&&met_pt()>45) SFOS[13] = 1;//invert DPhi
	if(nj<2&&pass2&&DPhilllMET<2.5&&pTlll>60&&met_pt()>55) SFOS[13] = 2;//invert DPhi

	if(nj<2&&pass0&&DPhilllMET>2.7&&pTlll<60)              SFOS[15] = 0;//invert Pt
	if(nj<2&&pass1&&DPhilllMET>2.5&&pTlll<60&&met_pt()>45) SFOS[15] = 1;//invert Pt
	if(nj<2&&pass2&&DPhilllMET>2.5&&pTlll<60&&met_pt()>55) SFOS[15] = 2;//invert Pt
	
	if(nj<2&&pass0X&&DPhilllMET<2.7&&pTlll>60)              SFOS[14] = 0;//invert DPhi
	if(nj<2&&pass1X&&DPhilllMET<2.5&&pTlll>60&&met_pt()>45) SFOS[14] = 1;//invert DPhi
	if(nj<2&&pass2X&&DPhilllMET<2.5&&pTlll>60&&met_pt()>55) SFOS[14] = 2;//invert DPhi

	if(nj<2&&pass0X&&DPhilllMET>2.7&&pTlll<60)              SFOS[16] = 0;//invert Pt
	if(nj<2&&pass1X&&DPhilllMET>2.5&&pTlll<60&&met_pt()>45) SFOS[16] = 1;//invert Pt
	if(nj<2&&pass2X&&DPhilllMET>2.5&&pTlll<60&&met_pt()>55) SFOS[16] = 2;//invert Pt
      }
      //17: 3lSR+1j+eq3l with MZ, 18: 3lCR+1j+eq3l  w/o MZ, 19: 3lSR+2j+ge3l with MZ, 20: 3lCR+2j+ge3l  w/o MZ
      if(nj<2&&nb==0&&n3l>=3&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)>10.){
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
	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass0X = false;	  
	}
	if(SFOScounter==1){
	  pass1 = true;
	  pass1X = true;
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
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	if(pass0 &&nveto3l==0&&n3l==3) SFOS[17] = 0;
	if(pass1 &&nveto3l==0&&n3l==3) SFOS[17] = 1;
	if(pass2 &&nveto3l==0&&n3l==3) SFOS[17] = 2;
	if(pass0X&&nveto3l==0&&n3l==3) SFOS[18] = 0;
	if(pass1X&&nveto3l==0&&n3l==3) SFOS[18] = 1;
	if(pass2X&&nveto3l==0&&n3l==3) SFOS[18] = 2;

	if(pass0 &&            n3l>=3) SFOS[19] = 0;
	if(pass1 &&            n3l>=3) SFOS[19] = 1;
	if(pass2 &&            n3l>=3) SFOS[19] = 2;
	if(pass0X&&            n3l>=3) SFOS[20] = 0;
	if(pass1X&&            n3l>=3) SFOS[20] = 1;
	if(pass2X&&            n3l>=3) SFOS[20] = 2;
      }

      //SRpresel but no NJ>=2
      if(/*nj<2&&*/nb==0/*&&nveto3l==0*/&&n3l==3&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)>10.){
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
	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  pass0X = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0X = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0X = false;
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass0X = false;
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
	//CR but ge2j
	if(/*nj>=2&&DPhilllMET>2.7&&pTlll>60&&*/pass0X) SFOS[21] = 0;
	if(/*nj>=2&&DPhilllMET>2.5&&pTlll>60&&*/pass1X) SFOS[21] = 1;
	if(/*nj>=2&&DPhilllMET>2.5&&pTlll>60&&*/pass2X) SFOS[21] = 2;
      }

      for(int i = 0; i<50; ++i){
	if((ee[i]||em[i]||mm[i])&&n3l>=3){
	  bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	  bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	  bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	  bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	  bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	  bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<10.) { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if(nSS>=1){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[0] ])==13&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[0] ])==13&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[iSS[0] ])==13&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if(nSS>=2){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[1] ])==13&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[1] ])==13&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[iSS[1] ])==13&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if(n3l>=1){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==13&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[0] ])==13&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[i3l[0] ])==13&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i]=-1; }
	}
	if(n3l>=2){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==13&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[1] ])==13&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[i3l[1] ])==13&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i]=-1; }
	}
	if(n3l>=3){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[2] ])==13&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[2] ])==13&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i]=-1; }
	  if(fname.find("ttbar_")  !=string::npos&&abs(lep_pdgId()[i3l[2] ])==13&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i]=-1; }
	}
      }
      
      if(ee[0])      *eventstocheckEE    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(em[0])      *eventstocheckEM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(mm[0])      *eventstocheckMM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(SFOS[0]==0) *eventstocheck0SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(SFOS[0]==1) *eventstocheck1SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(SFOS[0]==2) *eventstocheck2SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;

      if(     SFOS[21]>=1&&abs(lep_pdgId()[i3l[0] ])==11&&abs(lep_pdgId()[i3l[1] ])==11&&abs(lep_pdgId()[i3l[2] ])==11) histos["NJ_eee_CRlike_allSFOS_"      +sn2]->Fill(nj,weight);
      else if(SFOS[21]>=1&&abs(lep_pdgId()[i3l[0] ])==13&&abs(lep_pdgId()[i3l[1] ])==13&&abs(lep_pdgId()[i3l[2] ])==13) histos["NJ_mumumu_CRlike_allSFOS_"   +sn2]->Fill(nj,weight);
      else if(SFOS[21]>=1                                                                                             ) histos["NJ_noteeemmm_CRlike_allSFOS_"+sn2]->Fill(nj,weight);

      //cout << sn << " " << sn2 << endl;
      if(!isData()){
	if(ee[0])       histos["YieldsSR_raw_"     +sn ]->Fill(0.,1.);
	if(em[0])       histos["YieldsSR_raw_"     +sn ]->Fill(1.,1.);
	if(mm[0])       histos["YieldsSR_raw_"     +sn ]->Fill(2.,1.);
	if(SFOS[0]==0)  histos["YieldsSR_raw_"     +sn2]->Fill(3.,1.);
	if(SFOS[0]==1)  histos["YieldsSR_raw_"     +sn2]->Fill(4.,1.);
	if(SFOS[0]==2)  histos["YieldsSR_raw_"     +sn2]->Fill(5.,1.);
	if(ee[0])       histos["YieldsSR_rawweight_"     +sn ]->Fill(0.,weight);
	if(em[0])       histos["YieldsSR_rawweight_"     +sn ]->Fill(1.,weight);
	if(mm[0])       histos["YieldsSR_rawweight_"     +sn ]->Fill(2.,weight);
	if(SFOS[0]==0)  histos["YieldsSR_rawweight_"     +sn2]->Fill(3.,weight);
	if(SFOS[0]==1)  histos["YieldsSR_rawweight_"     +sn2]->Fill(4.,weight);
	if(SFOS[0]==2)  histos["YieldsSR_rawweight_"     +sn2]->Fill(5.,weight);
	if(ee[0])       histos["YieldsSR_"     +sn ]->Fill(0.,weight*lepSFSS);
	if(em[0])       histos["YieldsSR_"     +sn ]->Fill(1.,weight*lepSFSS);
	if(mm[0])       histos["YieldsSR_"     +sn ]->Fill(2.,weight*lepSFSS);
	if(SFOS[0]==0)  histos["YieldsSR_"     +sn2]->Fill(3.,weight*lepSF3l);
	if(SFOS[0]==1)  histos["YieldsSR_"     +sn2]->Fill(4.,weight*lepSF3l);
	if(SFOS[0]==2)  histos["YieldsSR_"     +sn2]->Fill(5.,weight*lepSF3l);
	if(ee[0])       histos["YieldsSR_preselection_"     +sn ]->Fill(0.,weight*lepSFSS);
	if(em[0])       histos["YieldsSR_preselection_"     +sn ]->Fill(1.,weight*lepSFSS);
	if(mm[0])       histos["YieldsSR_preselection_"     +sn ]->Fill(2.,weight*lepSFSS);
	if(SFOS[0]==0)  histos["YieldsSR_preselection_"     +sn2]->Fill(3.,weight*lepSF3l);
	if(SFOS[0]==1)  histos["YieldsSR_preselection_"     +sn2]->Fill(4.,weight*lepSF3l);
	if(SFOS[0]==2)  histos["YieldsSR_preselection_"     +sn2]->Fill(5.,weight*lepSF3l);
	if(ee[30])       histos["YieldsSR_jesup_"  +sn ]->Fill(0.,weight*lepSFSS);
	if(em[30])       histos["YieldsSR_jesup_"  +sn ]->Fill(1.,weight*lepSFSS);
	if(mm[30])       histos["YieldsSR_jesup_"  +sn ]->Fill(2.,weight*lepSFSS);
	if(ee[32])       histos["YieldsSR_jesdown_"+sn ]->Fill(0.,weight*lepSFSS);
	if(em[32])       histos["YieldsSR_jesdown_"+sn ]->Fill(1.,weight*lepSFSS);
	if(mm[32])       histos["YieldsSR_jesdown_"+sn ]->Fill(2.,weight*lepSFSS);
	if(ee[0])       histos["YieldsSR_lepSFup_"     +sn ]->Fill(0.,weight*(lepSFSS+lepSFerrSS));
	if(em[0])       histos["YieldsSR_lepSFup_"     +sn ]->Fill(1.,weight*(lepSFSS+lepSFerrSS));
	if(mm[0])       histos["YieldsSR_lepSFup_"     +sn ]->Fill(2.,weight*(lepSFSS+lepSFerrSS));
	if(SFOS[0]==0)  histos["YieldsSR_lepSFup_"     +sn2]->Fill(3.,weight*(lepSF3l+lepSFerr3l));
	if(SFOS[0]==1)  histos["YieldsSR_lepSFup_"     +sn2]->Fill(4.,weight*(lepSF3l+lepSFerr3l));
	if(SFOS[0]==2)  histos["YieldsSR_lepSFup_"     +sn2]->Fill(5.,weight*(lepSF3l+lepSFerr3l));
	if(ee[0])       histos["YieldsSR_lepSFdn_"     +sn ]->Fill(0.,weight*(lepSFSS-lepSFerrSS));
	if(em[0])       histos["YieldsSR_lepSFdn_"     +sn ]->Fill(1.,weight*(lepSFSS-lepSFerrSS));
	if(mm[0])       histos["YieldsSR_lepSFdn_"     +sn ]->Fill(2.,weight*(lepSFSS-lepSFerrSS));
	if(SFOS[0]==0)  histos["YieldsSR_lepSFdn_"     +sn2]->Fill(3.,weight*(lepSF3l-lepSFerr3l));
	if(SFOS[0]==1)  histos["YieldsSR_lepSFdn_"     +sn2]->Fill(4.,weight*(lepSF3l-lepSFerr3l));
	if(SFOS[0]==2)  histos["YieldsSR_lepSFdn_"     +sn2]->Fill(5.,weight*(lepSF3l-lepSFerr3l));
      }
      if(!isData()){
	if(ee[5]&&ee[0]) histos["Mjj_SRlike_allSS_"+sn]->Fill(Mjj,weight*lepSFSS);
	if(em[5]&&em[0]) histos["Mjj_SRlike_allSS_"+sn]->Fill(Mjj,weight*lepSFSS);
	if(mm[5]&&mm[0]) histos["Mjj_SRlike_allSS_"+sn]->Fill(Mjj,weight*lepSFSS);
	if(ee[31]&&ee[30]) histos["Mjj_SRlike_allSS_jesup_"  +sn]->Fill(Mjj,weight*lepSFSS);
	if(em[31]&&em[30]) histos["Mjj_SRlike_allSS_jesup_"  +sn]->Fill(Mjj,weight*lepSFSS);
	if(mm[31]&&mm[30]) histos["Mjj_SRlike_allSS_jesup_"  +sn]->Fill(Mjj,weight*lepSFSS);
	if(ee[33]&&ee[32]) histos["Mjj_SRlike_allSS_jesdown_"+sn]->Fill(Mjj,weight*lepSFSS);
	if(em[33]&&em[32]) histos["Mjj_SRlike_allSS_jesdown_"+sn]->Fill(Mjj,weight*lepSFSS);
	if(mm[33]&&mm[32]) histos["Mjj_SRlike_allSS_jesdown_"+sn]->Fill(Mjj,weight*lepSFSS);
      }
      if(ee[45]&&!ee[0])       histos["YieldsSR_preselection_"     +sn ]->Fill(0.,weight*lepSFSS);
      if(em[45]&&!em[0])       histos["YieldsSR_preselection_"     +sn ]->Fill(1.,weight*lepSFSS);
      if(mm[45]&&!mm[0])       histos["YieldsSR_preselection_"     +sn ]->Fill(2.,weight*lepSFSS);
      if(SFOS[17]==0&&!(SFOS[0]==0))  histos["YieldsSR_preselection_"     +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[17]==1&&!(SFOS[0]==1))  histos["YieldsSR_preselection_"     +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[17]==2&&!(SFOS[0]==2))  histos["YieldsSR_preselection_"     +sn2]->Fill(5.,weight*lepSFSS);
      if(ee[5]&&!ee[0]) histos["Mjj_SRlike_allSS_"+sn]->Fill(Mjj,weight*lepSFSS);
      if(em[5]&&!em[0]) histos["Mjj_SRlike_allSS_"+sn]->Fill(Mjj,weight*lepSFSS);
      if(mm[5]&&!mm[0]) histos["Mjj_SRlike_allSS_"+sn]->Fill(Mjj,weight*lepSFSS);
      if(ee[31]&&!ee[30]) histos["Mjj_SRlike_allSS_jesup_"  +sn]->Fill(Mjj,weight*lepSFSS);
      if(em[31]&&!em[30]) histos["Mjj_SRlike_allSS_jesup_"  +sn]->Fill(Mjj,weight*lepSFSS);
      if(mm[31]&&!mm[30]) histos["Mjj_SRlike_allSS_jesup_"  +sn]->Fill(Mjj,weight*lepSFSS);
      if(ee[33]&&!ee[32]) histos["Mjj_SRlike_allSS_jesdown_"+sn]->Fill(Mjj,weight*lepSFSS);
      if(em[33]&&!em[32]) histos["Mjj_SRlike_allSS_jesdown_"+sn]->Fill(Mjj,weight*lepSFSS);
      if(mm[33]&&!mm[32]) histos["Mjj_SRlike_allSS_jesdown_"+sn]->Fill(Mjj,weight*lepSFSS);
      //Mll SR is below	
      if(ee[5]&&fabs(Mjj-80.)>20.) histos["YieldsSR_invertdPhiPtfor3Mjjfor2l_"+sn ]->Fill(0.,weight*lepSFSS);
      if(em[5]&&fabs(Mjj-80.)>20.) histos["YieldsSR_invertdPhiPtfor3Mjjfor2l_"+sn ]->Fill(1.,weight*lepSFSS);
      if(mm[5]&&fabs(Mjj-80.)>20.) histos["YieldsSR_invertdPhiPtfor3Mjjfor2l_"+sn ]->Fill(2.,weight*lepSFSS);
      if(SFOS[2]==0)  histos["YieldsSR_invertdPhiPtfor3Mjjfor2l_"             +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[2]==1)  histos["YieldsSR_invertdPhiPtfor3Mjjfor2l_"             +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[2]==2)  histos["YieldsSR_invertdPhiPtfor3Mjjfor2l_"             +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[6]&&fabs(Mjj-80.)>20.) histos["YieldsCR_invertdPhiPtfor3Mjjfor2l_"+sn2]->Fill(0.,weight*lepSF3l);
      if(em[6]&&fabs(Mjj-80.)>20.) histos["YieldsCR_invertdPhiPtfor3Mjjfor2l_"+sn2]->Fill(1.,weight*lepSF3l);
      if(mm[6]&&fabs(Mjj-80.)>20.) histos["YieldsCR_invertdPhiPtfor3Mjjfor2l_"+sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[3]==0)  histos["YieldsCR_invertdPhiPtfor3Mjjfor2l_"             +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[3]==1)  histos["YieldsCR_invertdPhiPtfor3Mjjfor2l_"             +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[3]==2)  histos["YieldsCR_invertdPhiPtfor3Mjjfor2l_"             +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[1])       histos["YieldsCR_using3lorMll_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(em[1])       histos["YieldsCR_using3lorMll_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[1])       histos["YieldsCR_using3lorMll_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[1]==0)  histos["YieldsCR_using3lorMll_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[1]==1)  histos["YieldsCR_using3lorMll_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[1]==2)  histos["YieldsCR_using3lorMll_"     +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[34])      histos["YieldsCR_using3lorMll_jesup_"  +sn2]->Fill(0.,weight*lepSF3l);
      if(em[34])      histos["YieldsCR_using3lorMll_jesup_"  +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[34])      histos["YieldsCR_using3lorMll_jesup_"  +sn2]->Fill(2.,weight*lepSF3l);
      if(ee[36])      histos["YieldsCR_using3lorMll_jesdown_"+sn2]->Fill(0.,weight*lepSF3l);
      if(em[36])      histos["YieldsCR_using3lorMll_jesdown_"+sn2]->Fill(1.,weight*lepSF3l);
      if(mm[36])      histos["YieldsCR_using3lorMll_jesdown_"+sn2]->Fill(2.,weight*lepSF3l);
      if(ee[6])       histos["YieldsCR_dropMjj_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(em[6])       histos["YieldsCR_dropMjj_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[6])       histos["YieldsCR_dropMjj_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[1]==0)  histos["YieldsCR_dropMjj_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[1]==1)  histos["YieldsCR_dropMjj_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[1]==2)  histos["YieldsCR_dropMjj_"     +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[35])      histos["YieldsCR_dropMjj_jesup_"  +sn2]->Fill(0.,weight*lepSF3l);
      if(em[35])      histos["YieldsCR_dropMjj_jesup_"  +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[35])      histos["YieldsCR_dropMjj_jesup_"  +sn2]->Fill(2.,weight*lepSF3l);
      if(ee[37])      histos["YieldsCR_dropMjj_jesdown_"+sn2]->Fill(0.,weight*lepSF3l);
      if(em[37])      histos["YieldsCR_dropMjj_jesdown_"+sn2]->Fill(1.,weight*lepSF3l);
      if(mm[37])      histos["YieldsCR_dropMjj_jesdown_"+sn2]->Fill(2.,weight*lepSF3l);
      if(ee[3])       histos["YieldsSR_lowMET_"     +sn ]->Fill(0.,weight*lepSFSS);
      if(em[3])       histos["YieldsSR_lowMET_"     +sn ]->Fill(1.,weight*lepSFSS);
      if(mm[3])       histos["YieldsSR_lowMET_"     +sn ]->Fill(2.,weight*lepSFSS);
      if(SFOS[4]==0)  histos["YieldsSR_lowMET_"     +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[4]==1)  histos["YieldsSR_lowMET_"     +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[4]==2)  histos["YieldsSR_lowMET_"     +sn2]->Fill(5.,weight*lepSFSS);
      if(ee[8]&&fabs(Mjj-80.)>20.) histos["YieldsSR_lowMETinvertdPhiPtMjj_"     +sn ]->Fill(0.,weight*lepSFSS);
      if(em[8]&&fabs(Mjj-80.)>20.) histos["YieldsSR_lowMETinvertdPhiPtMjj_"     +sn ]->Fill(1.,weight*lepSFSS);
      if(mm[8]&&fabs(Mjj-80.)>20.) histos["YieldsSR_lowMETinvertdPhiPtMjj_"     +sn ]->Fill(2.,weight*lepSFSS);
      if(SFOS[6]==0)               histos["YieldsSR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[6]==1)               histos["YieldsSR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[6]==2)               histos["YieldsSR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(5.,weight*lepSFSS);
      if(ee[4])       histos["YieldsCR_lowMET_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(em[4])       histos["YieldsCR_lowMET_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[4])       histos["YieldsCR_lowMET_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[5]==0)  histos["YieldsCR_lowMET_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[5]==1)  histos["YieldsCR_lowMET_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[5]==2)  histos["YieldsCR_lowMET_"     +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[9]&&fabs(Mjj-80.)>20.) histos["YieldsCR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(em[9]&&fabs(Mjj-80.)>20.) histos["YieldsCR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[9]&&fabs(Mjj-80.)>20.) histos["YieldsCR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[7]==0)               histos["YieldsCR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[7]==1)               histos["YieldsCR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[7]==2)               histos["YieldsCR_lowMETinvertdPhiPtMjj_"     +sn2]->Fill(5.,weight*lepSF3l);

      if(SFOS[ 4]==0)  histos["YieldsSR_invertMET_"        +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[ 4]==1)  histos["YieldsSR_invertMET_"        +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[ 4]==2)  histos["YieldsSR_invertMET_"        +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[ 5]==0)  histos["YieldsCR_invertMET_"        +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[ 5]==1)  histos["YieldsCR_invertMET_"        +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[ 5]==2)  histos["YieldsCR_invertMET_"        +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[13]==0)  histos["YieldsSR_invertDPhi_"       +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[13]==1)  histos["YieldsSR_invertDPhi_"       +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[13]==2)  histos["YieldsSR_invertDPhi_"       +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[14]==0)  histos["YieldsCR_invertDPhi_"       +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[14]==1)  histos["YieldsCR_invertDPhi_"       +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[14]==2)  histos["YieldsCR_invertDPhi_"       +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[15]==0)  histos["YieldsSR_invertPt_"         +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[15]==1)  histos["YieldsSR_invertPt_"         +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[15]==2)  histos["YieldsSR_invertPt_"         +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[16]==0)  histos["YieldsCR_invertPt_"         +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[16]==1)  histos["YieldsCR_invertPt_"         +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[16]==2)  histos["YieldsCR_invertPt_"         +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[ 2]==0)  histos["YieldsSR_invertDPhiPt_"     +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[ 2]==1)  histos["YieldsSR_invertDPhiPt_"     +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[ 2]==2)  histos["YieldsSR_invertDPhiPt_"     +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[ 3]==0)  histos["YieldsCR_invertDPhiPt_"     +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[ 3]==1)  histos["YieldsCR_invertDPhiPt_"     +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[ 3]==2)  histos["YieldsCR_invertDPhiPt_"     +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[ 6]==0)  histos["YieldsSR_invertMETDPhiPt_"  +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[ 6]==1)  histos["YieldsSR_invertMETDPhiPt_"  +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[ 6]==2)  histos["YieldsSR_invertMETDPhiPt_"  +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[ 7]==0)  histos["YieldsCR_invertMETDPhiPt_"  +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[ 7]==1)  histos["YieldsCR_invertMETDPhiPt_"  +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[ 7]==2)  histos["YieldsCR_invertMETDPhiPt_"  +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[10]==0)  histos["YieldsSR_invertMETDPhiOrPt_"+sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[10]==1)  histos["YieldsSR_invertMETDPhiOrPt_"+sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[10]==2)  histos["YieldsSR_invertMETDPhiOrPt_"+sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[11]==0)  histos["YieldsCR_invertMETDPhiOrPt_"+sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[11]==1)  histos["YieldsCR_invertMETDPhiOrPt_"+sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[11]==2)  histos["YieldsCR_invertMETDPhiOrPt_"+sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[10]==0             )  histos["YieldsSR_invertDPhiOrPt_"     +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[10]==1&&met_pt()>45)  histos["YieldsSR_invertDPhiOrPt_"     +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[10]==2&&met_pt()>55)  histos["YieldsSR_invertDPhiOrPt_"     +sn2]->Fill(5.,weight*lepSFSS);
      if(SFOS[11]==0             )  histos["YieldsCR_invertDPhiOrPt_"     +sn2]->Fill(3.,weight*lepSFSS);
      if(SFOS[11]==1&&met_pt()>45)  histos["YieldsCR_invertDPhiOrPt_"     +sn2]->Fill(4.,weight*lepSFSS);
      if(SFOS[11]==2&&met_pt()>55)  histos["YieldsCR_invertDPhiOrPt_"     +sn2]->Fill(5.,weight*lepSFSS);

      if(ee[ 6]||ee[13]) histos["YieldsCR_SSany_dropMjj_rawweight_"     +sn2]->Fill(0.,weight);
      if(em[ 6]||em[13]) histos["YieldsCR_SSany_dropMjj_rawweight_"     +sn2]->Fill(1.,weight);
      if(mm[ 6]||mm[13]) histos["YieldsCR_SSany_dropMjj_rawweight_"     +sn2]->Fill(2.,weight);
      if(SFOS[1]==0)     histos["YieldsCR_SSany_dropMjj_rawweight_"     +sn2]->Fill(3.,weight);
      if(SFOS[1]==1)     histos["YieldsCR_SSany_dropMjj_rawweight_"     +sn2]->Fill(4.,weight);
      if(SFOS[1]==2)     histos["YieldsCR_SSany_dropMjj_rawweight_"     +sn2]->Fill(5.,weight);
      
      if(ee[46])         histos["YieldsCR_SSany_dropMjj_test_"+sn2]->Fill(0.,weight*lepSF3l);
      if(em[46])         histos["YieldsCR_SSany_dropMjj_test_"+sn2]->Fill(1.,weight*lepSF3l);
      if(mm[46])         histos["YieldsCR_SSany_dropMjj_test_"+sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[1]==0)     histos["YieldsCR_SSany_dropMjj_test_"+sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[1]==1)     histos["YieldsCR_SSany_dropMjj_test_"+sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[1]==2)     histos["YieldsCR_SSany_dropMjj_test_"+sn2]->Fill(5.,weight*lepSF3l);
      
      if(ee[ 6]||ee[13]) histos["YieldsCR_SSany_dropMjj_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(em[ 6]||em[13]) histos["YieldsCR_SSany_dropMjj_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[ 6]||mm[13]) histos["YieldsCR_SSany_dropMjj_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[1]==0)     histos["YieldsCR_SSany_dropMjj_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[1]==1)     histos["YieldsCR_SSany_dropMjj_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[1]==2)     histos["YieldsCR_SSany_dropMjj_"     +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[35]||ee[43]) histos["YieldsCR_SSany_dropMjj_jesup_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(em[35]||em[43]) histos["YieldsCR_SSany_dropMjj_jesup_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[35]||mm[43]) histos["YieldsCR_SSany_dropMjj_jesup_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(ee[37]||ee[41]) histos["YieldsCR_SSany_dropMjj_jesdown_"   +sn2]->Fill(0.,weight*lepSF3l);
      if(em[37]||em[41]) histos["YieldsCR_SSany_dropMjj_jesdown_"   +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[37]||mm[41]) histos["YieldsCR_SSany_dropMjj_jesdown_"   +sn2]->Fill(2.,weight*lepSF3l);
      if(ee[14]||ee[19]) histos["YieldsCR_SSany_dropMjj_butnoMll_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(em[14]||em[19]) histos["YieldsCR_SSany_dropMjj_butnoMll_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(mm[14]||mm[19]) histos["YieldsCR_SSany_dropMjj_butnoMll_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(SFOS[1]==0)     histos["YieldsCR_SSany_dropMjj_butnoMll_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[1]==1)     histos["YieldsCR_SSany_dropMjj_butnoMll_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[1]==2)     histos["YieldsCR_SSany_dropMjj_butnoMll_"     +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[ 6]||ee[13]||em[ 6]||em[13]||mm[ 6]||mm[13]) histos["Mjj_CRlike_SSany_allSS_"+sn2]->Fill(Mjj,weight*lepSF3l);
      if(ee[35]||ee[43]||em[35]||em[43]||mm[35]||mm[43]) histos["Mjj_CRlike_SSany_allSS_jesup_"  +sn2]->Fill(Mjj_up,weight*lepSF3l);
      if(ee[37]||ee[41]||em[37]||em[41]||mm[37]||mm[41]) histos["Mjj_CRlike_SSany_allSS_jesdown_"+sn2]->Fill(Mjj_dn,weight*lepSF3l);
      if((ee[24]||ee[44]||em[24]||em[44]||mm[24]||mm[44])&&!(ee[ 6]||ee[13]||em[ 6]||em[13]||mm[ 6]||mm[13])) histos["Mjj_CRloose_butnotlike_SSany_allSS_"+sn2]->Fill(Mjj,weight*lepSF3l);
      if(MjjL<400.&&((ee[24]||ee[44]||em[24]||em[44]||mm[24]||mm[44])&&!(ee[ 6]||ee[13]||em[ 6]||em[13]||mm[ 6]||mm[13]))) histos["Mjj_CRloose_butnotlike_SSany_allSS_wMjjL_"+sn2]->Fill(Mjj,weight*lepSF3l);
      if(MjjL<400.&&fabs(Detajj)<1.5&&((ee[24]||ee[44]||em[24]||em[44]||mm[24]||mm[44])&&!(ee[ 6]||ee[13]||em[ 6]||em[13]||mm[ 6]||mm[13]))) histos["Mjj_CRloose_butnotlike_SSany_allSS_wMjjLDeta_"+sn2]->Fill(Mjj,weight*lepSF3l);
      if(fabs(Mjj   -80.)<=20.&&(ee[ 6]||ee[13])) histos["YieldsCR_SSany_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(fabs(Mjj   -80.)<=20.&&(em[ 6]||em[13])) histos["YieldsCR_SSany_"     +sn2]->Fill(1.,weight*lepSF3l); 
      if(fabs(Mjj   -80.)<=20.&&(mm[ 6]||mm[13])) histos["YieldsCR_SSany_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(fabs(Mjj   -80.)<=20.&&(SFOS[1]==0))     histos["YieldsCR_SSany_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(fabs(Mjj   -80.)<=20.&&(SFOS[1]==1))     histos["YieldsCR_SSany_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(fabs(Mjj   -80.)<=20.&&(SFOS[1]==2))     histos["YieldsCR_SSany_"     +sn2]->Fill(5.,weight*lepSF3l);
      if(fabs(Mjj_up-80.)<=20.&&(ee[35]||ee[43])) histos["YieldsCR_SSany_jesup_"     +sn2]->Fill(0.,weight*lepSF3l);
      if(fabs(Mjj_up-80.)<=20.&&(em[35]||em[43])) histos["YieldsCR_SSany_jesup_"     +sn2]->Fill(1.,weight*lepSF3l);
      if(fabs(Mjj_up-80.)<=20.&&(mm[35]||mm[43])) histos["YieldsCR_SSany_jesup_"     +sn2]->Fill(2.,weight*lepSF3l);
      if(fabs(Mjj_dn-80.)<=20.&&(ee[37]||ee[41])) histos["YieldsCR_SSany_jesdown_"   +sn2]->Fill(0.,weight*lepSF3l);
      if(fabs(Mjj_dn-80.)<=20.&&(em[37]||em[41])) histos["YieldsCR_SSany_jesdown_"   +sn2]->Fill(1.,weight*lepSF3l);
      if(fabs(Mjj_dn-80.)<=20.&&(mm[37]||mm[41])) histos["YieldsCR_SSany_jesdown_"   +sn2]->Fill(2.,weight*lepSF3l);

      if(ee[ 6]||ee[13]) histos["YieldsCR_SSany_dropMjj_lepSFup_"     +sn2]->Fill(0.,weight*(lepSF3l+lepSFerr3l));
      if(em[ 6]||em[13]) histos["YieldsCR_SSany_dropMjj_lepSFup_"     +sn2]->Fill(1.,weight*(lepSF3l+lepSFerr3l));
      if(mm[ 6]||mm[13]) histos["YieldsCR_SSany_dropMjj_lepSFup_"     +sn2]->Fill(2.,weight*(lepSF3l+lepSFerr3l));
      if(SFOS[1]==0)     histos["YieldsCR_SSany_dropMjj_lepSFup_"     +sn2]->Fill(3.,weight*(lepSF3l+lepSFerr3l));
      if(SFOS[1]==1)     histos["YieldsCR_SSany_dropMjj_lepSFup_"     +sn2]->Fill(4.,weight*(lepSF3l+lepSFerr3l));
      if(SFOS[1]==2)     histos["YieldsCR_SSany_dropMjj_lepSFup_"     +sn2]->Fill(5.,weight*(lepSF3l+lepSFerr3l));
      if(ee[ 6]||ee[13]) histos["YieldsCR_SSany_dropMjj_lepSFdn_"     +sn2]->Fill(0.,weight*(lepSF3l-lepSFerr3l));
      if(em[ 6]||em[13]) histos["YieldsCR_SSany_dropMjj_lepSFdn_"     +sn2]->Fill(1.,weight*(lepSF3l-lepSFerr3l));
      if(mm[ 6]||mm[13]) histos["YieldsCR_SSany_dropMjj_lepSFdn_"     +sn2]->Fill(2.,weight*(lepSF3l-lepSFerr3l));
      if(SFOS[1]==0)     histos["YieldsCR_SSany_dropMjj_lepSFdn_"     +sn2]->Fill(3.,weight*(lepSF3l-lepSFerr3l));
      if(SFOS[1]==1)     histos["YieldsCR_SSany_dropMjj_lepSFdn_"     +sn2]->Fill(4.,weight*(lepSF3l-lepSFerr3l));
      if(SFOS[1]==2)     histos["YieldsCR_SSany_dropMjj_lepSFdn_"     +sn2]->Fill(5.,weight*(lepSF3l-lepSFerr3l));
      
      if(ee[20])      histos["YieldsCR_dropcutsbutNJisotr_"       +sn ]->Fill(0.,weight*lepSF3l);
      if(em[20])      histos["YieldsCR_dropcutsbutNJisotr_"       +sn ]->Fill(1.,weight*lepSF3l);
      if(mm[20])      histos["YieldsCR_dropcutsbutNJisotr_"       +sn ]->Fill(2.,weight*lepSF3l);
      if(ee[21])      histos["YieldsCRnoMll_dropcutsbutNJisotr_"  +sn ]->Fill(0.,weight*lepSF3l);
      if(em[21])      histos["YieldsCRnoMll_dropcutsbutNJisotr_"  +sn ]->Fill(1.,weight*lepSF3l);
      if(mm[21])      histos["YieldsCRnoMll_dropcutsbutNJisotr_"  +sn ]->Fill(2.,weight*lepSF3l);
      if(ee[22])      histos["YieldsCR_dropcutsbutNJvetolep_"     +sn ]->Fill(0.,weight*lepSF3l);
      if(em[22])      histos["YieldsCR_dropcutsbutNJvetolep_"     +sn ]->Fill(1.,weight*lepSF3l);
      if(mm[22])      histos["YieldsCR_dropcutsbutNJvetolep_"     +sn ]->Fill(2.,weight*lepSF3l);
      if(SFOS[18]==0) histos["YieldsCR_dropcutsbutNJvetolep_"     +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[18]==1) histos["YieldsCR_dropcutsbutNJvetolep_"     +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[18]==2) histos["YieldsCR_dropcutsbutNJvetolep_"     +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[23])      histos["YieldsCRnoMll_dropcutsbutNJvetolep_"+sn ]->Fill(0.,weight*lepSF3l);
      if(em[23])      histos["YieldsCRnoMll_dropcutsbutNJvetolep_"+sn ]->Fill(1.,weight*lepSF3l);
      if(mm[23])      histos["YieldsCRnoMll_dropcutsbutNJvetolep_"+sn ]->Fill(2.,weight*lepSF3l);
      if(ee[24])      histos["YieldsCR_dropcutsbutNJ_"            +sn ]->Fill(0.,weight*lepSF3l);
      if(em[24])      histos["YieldsCR_dropcutsbutNJ_"            +sn ]->Fill(1.,weight*lepSF3l);
      if(mm[24])      histos["YieldsCR_dropcutsbutNJ_"            +sn ]->Fill(2.,weight*lepSF3l);
      if(SFOS[20]==0) histos["YieldsCR_dropcutsbutNJ_"            +sn2]->Fill(3.,weight*lepSF3l);
      if(SFOS[20]==1) histos["YieldsCR_dropcutsbutNJ_"            +sn2]->Fill(4.,weight*lepSF3l);
      if(SFOS[20]==2) histos["YieldsCR_dropcutsbutNJ_"            +sn2]->Fill(5.,weight*lepSF3l);
      if(ee[38])      histos["YieldsCR_dropcutsbutNJ_jesup_"  +sn ]->Fill(0.,weight*lepSF3l);
      if(em[38])      histos["YieldsCR_dropcutsbutNJ_jesup_"  +sn ]->Fill(1.,weight*lepSF3l);
      if(mm[38])      histos["YieldsCR_dropcutsbutNJ_jesup_"  +sn ]->Fill(2.,weight*lepSF3l);
      if(ee[39])      histos["YieldsCR_dropcutsbutNJ_jesdown_"+sn ]->Fill(0.,weight*lepSF3l);
      if(em[39])      histos["YieldsCR_dropcutsbutNJ_jesdown_"+sn ]->Fill(1.,weight*lepSF3l);
      if(mm[39])      histos["YieldsCR_dropcutsbutNJ_jesdown_"+sn ]->Fill(2.,weight*lepSF3l);
      if(ee[25])      histos["YieldsCRnoMll_dropcutsbutNJ_"       +sn ]->Fill(0.,weight*lepSF3l);
      if(em[25])      histos["YieldsCRnoMll_dropcutsbutNJ_"       +sn ]->Fill(1.,weight*lepSF3l);
      if(mm[25])      histos["YieldsCRnoMll_dropcutsbutNJ_"       +sn ]->Fill(2.,weight*lepSF3l);
      if(ee[14]||em[14]||mm[14]) histos["Mll_CRlike_allSS_"                  +sn2]->Fill(MllCR ,weight*lepSF3l);
      if(ee[16]||em[16]||mm[16]) histos["Mll_CRlike_allSS_SSpairnotlead2_"   +sn2]->Fill(MllCRb,weight*lepSF3l);
      if(ee[18]||em[18]||mm[18]) histos["Mll_CRlike_allSS_SSpairnotlead2v2_" +sn2]->Fill(MllCRb,weight*lepSF3l);
      if(ee[15]||em[15]||mm[15]) histos["Mll_CRnoMjj_allSS_"                 +sn2]->Fill(MllCR ,weight*lepSF3l);
      if(ee[17]||em[17]||mm[17]) histos["Mll_CRnoMjj_allSS_SSpairnotlead2_"  +sn2]->Fill(MllCRb,weight*lepSF3l);
      if(ee[19]||em[19]||mm[19]) histos["Mll_CRnoMjj_allSS_SSpairnotlead2v2_"+sn2]->Fill(MllCRb,weight*lepSF3l);
      if(ee[ 8]||em[ 8]||mm[ 8]) histos["Mjj_lowMET_SRlike_allSS_"          +sn ]->Fill(Mjj,weight*lepSFSS);
      if(ee[ 9]||em[ 9]||mm[ 9]) histos["Mjj_lowMET_CRlike_allSS_"          +sn2]->Fill(Mjj,weight*lepSF3l);
      if(ee[ 7]||em[ 7]||mm[ 7]) histos["Mjj_CRlike_nolowMll_allSS_"        +sn2]->Fill(Mjj,weight*lepSF3l);
      if(ee[12]||em[12]||mm[12]) histos["Mjj_CRlike_allSS_SSpairnotlead2_"  +sn2]->Fill(Mjj,weight*lepSF3l);
      if(ee[13]||em[13]||mm[13]) histos["Mjj_CRlike_allSS_SSpairnotlead2v2_"+sn2]->Fill(Mjj,weight*lepSF3l);
      if(ee[ 6]||em[ 6]||mm[ 6]) histos["Mjj_CRlike_allSS_"                 +sn2]->Fill(Mjj,weight*lepSF3l);
      if(ee[6]!=(ee[20]&&passMll&&passDetajj&&met_pt()>40.))               cout << __LINE__ << " " << nj30 << " " << nb << " " << Mjj << " " << MjjL << " " << Detajj << " " << nisoTrack_mt2_cleaned_VVV_cutbased_veto() << " " << nveto3l << " " << nSS << " " << n3l << " " << lep_pdgId()[iSS[0] ]<< ":" << lep_pdgId()[iSS[1] ] << " " << (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M() << " " << met_pt() << endl;
      if(ee[6]!=(ee[20]&&passMll&&passDetajj&&met_pt()>40.))               cout <<ee[6] << " " << ee[20] << " " << passMll << " " << passDetajj << " " << (met_pt()>40.) << endl;
      if(em[6]!=(em[20]&&passMll&&passDetajj&&met_pt()>40.&&MTmax3l>=90.)) cout <<__LINE__<<endl;
      if(mm[6]!=(mm[20]&&passMll&&passDetajj))                             cout <<__LINE__<<endl;
      if(ee[24]||em[24]||mm[24]) histos["Mjj_CRloose_allSS_"   +sn2]->Fill(Mjj,weight);
      if(ee[35]||em[35]||mm[35]) histos["Mjj_CRlike_allSS_jesup_"                 +sn2]->Fill(Mjj_up,weight*lepSF3l);
      if(ee[37]||em[37]||mm[37]) histos["Mjj_CRlike_allSS_jesdown_"               +sn2]->Fill(Mjj_dn,weight*lepSF3l);
      if(ee[38]||em[38]||mm[38]) histos["Mjj_CRloose_allSS_jesup_"                +sn2]->Fill(Mjj_up,weight*lepSF3l);
      if(ee[39]||em[39]||mm[39]) histos["Mjj_CRloose_allSS_jesdown_"              +sn2]->Fill(Mjj_dn,weight*lepSF3l);

      if(ee[24]||em[24]||mm[24]) histos["Mll_CRloose_allSS_"   +sn2]->Fill((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M(),weight);
      if(ee[24]||em[24]||mm[24]) histos["MET_CRloose_allSS_"   +sn2]->Fill(met_pt(),weight*lepSF3l);
      if(ee[24]||em[24]||mm[24]) histos["Detajj_CRloose_allSS_"+sn2]->Fill(Detajj,weight*lepSF3l);
      if(        em[24]        ) histos["MTmax_CRloose_allSS_" +sn2]->Fill(MTmax3l,weight*lepSF3l);
      if(ee[20]&&passDetajj&&met_pt()>40.)               histos["Mll_CRlike_allSS_"   +sn2]->Fill((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M(),weight*lepSF3l);
      if(em[20]&&passDetajj&&met_pt()>40.&&MTmax3l>=90.) histos["Mll_CRlike_allSS_"   +sn2]->Fill((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M(),weight*lepSF3l);
      if(mm[20]&&passDetajj)                             histos["Mll_CRlike_allSS_"   +sn2]->Fill((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M(),weight*lepSF3l);
      if(ee[20]&&passMll&&passDetajj)                    histos["MET_CRlike_allSS_"   +sn2]->Fill(met_pt(),weight*lepSF3l);
      if(em[20]&&passMll&&passDetajj&&MTmax3l>=90.)      histos["MET_CRlike_allSS_"   +sn2]->Fill(met_pt(),weight*lepSF3l);
      if(mm[20]&&passMll&&passDetajj)                    histos["MET_CRlike_allSS_"   +sn2]->Fill(met_pt(),weight*lepSF3l);
      if(ee[20]&&passMll)                                histos["Detajj_CRlike_allSS_"+sn2]->Fill(Detajj,weight*lepSF3l);
      if(em[20]&&passMll&&MTmax3l>=90.)                  histos["Detajj_CRlike_allSS_"+sn2]->Fill(Detajj,weight*lepSF3l);
      if(mm[20]&&passMll)                                histos["Detajj_CRlike_allSS_"+sn2]->Fill(Detajj,weight*lepSF3l);
      if(em[20]&&passMll&&passDetajj&&met_pt()>40.)      histos["MTmax_CRlike_allSS_" +sn2]->Fill(MTmax3l,weight*lepSF3l);

      if(n3l>=3){
      	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	if(!isData()){
	  if(SFOS[0]==1){
	    if(OS01&&SF01) histos["Mll_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	    if(OS02&&SF02) histos["Mll_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	    if(OS12&&SF12) histos["Mll_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	    if(OS01&&SF01) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	    if(OS02&&SF02) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	    if(OS12&&SF12) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  }
	  if(SFOS[0]==2){
	    if(OS01&&SF01) histos["Mll_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	    if(OS02&&SF02) histos["Mll_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	    if(OS12&&SF12) histos["Mll_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	    if(OS01&&SF01) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	    if(OS02&&SF02) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	    if(OS12&&SF12) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	    if(OS01&&SF01&&OS02&&SF02){
	      if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)>fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)){
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	      }
	      else {
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	    }
	    if(OS01&&SF01&&OS12&&SF12){
	      if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)>fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)){
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	      }
	      else {
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	    }
	    if(OS02&&SF02&&OS12&&SF12){
	      if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)>fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)){
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	      else {
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	    }
	  }
	}
	if(SFOS[1]==1){
	  if(OS01&&SF01) histos["Mll_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_SRlike_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(SFOS[1]==2){
	  if(OS01&&SF01) histos["Mll_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_SRlike_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_SRlike_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01&&OS02&&SF02){
	      if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)>fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)){
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	      }
	      else {
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	    }
	    if(OS01&&SF01&&OS12&&SF12){
	      if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)>fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)){
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	      }
	      else {
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	    }
	    if(OS02&&SF02&&OS12&&SF12){
	      if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)>fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)){
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	      else {
		histos["Mll_SRlike_MllclosestZ_2SFOS_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
		histos["Mll_SRlike_MllfurthestZ_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	      }
	    }
	}
	if((SFOS[2]==1)||(SFOS[3]==1)){
	  if(OS01&&SF01) histos["Mll_invertdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(((SFOS[10]==1)||(SFOS[11]==1))&&met_pt()>45){
	  if(OS01&&SF01) histos["Mll_invertdPhiOrPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiOrPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiOrPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertdPhiOrPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiOrPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiOrPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[4]==1)||(SFOS[5]==1)){
	  if(OS01&&SF01) histos["Mll_invertMET_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMET_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMET_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertMET_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMET_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMET_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[6]==1)||(SFOS[7]==1)){
	  if(OS01&&SF01) histos["Mll_invertMETdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMETdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMETdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[2]==2)||(SFOS[3]==2)){
	  if(OS01&&SF01) histos["Mll_invertdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(((SFOS[10]==2)||(SFOS[11]==2))&&met_pt()>55){
	  if(OS01&&SF01) histos["Mll_invertdPhiOrPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiOrPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiOrPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertdPhiOrPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhiOrPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhiOrPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[4]==2)||(SFOS[5]==2)){
	  if(OS01&&SF01) histos["Mll_invertMET_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMET_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMET_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertMET_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMET_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMET_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[6]==2)||(SFOS[7]==2)){
	  if(OS01&&SF01) histos["Mll_invertMETdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMETdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMETdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[8]==1)||(SFOS[9]==1)){
	  if(OS01&&SF01) histos["Mll_ge2j_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_ge2j_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_ge2j_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_ge2j_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_ge2j_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_ge2j_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[8]==2)||(SFOS[9]==2)){
	  if(OS01&&SF01) histos["Mll_ge2j_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_ge2j_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_ge2j_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_ge2j_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_ge2j_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_ge2j_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[10]==1)||(SFOS[11]==1)){
	  if(OS01&&SF01) histos["Mll_inverteitherMETdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_inverteitherMETdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_inverteitherMETdPhiPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_inverteitherMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_inverteitherMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_inverteitherMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[10]==2)||(SFOS[11]==2)){
	  if(OS01&&SF01) histos["Mll_inverteitherMETdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_inverteitherMETdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_inverteitherMETdPhiPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_inverteitherMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_inverteitherMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_inverteitherMETdPhiPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[13]==1)||(SFOS[14]==1)){
	  if(OS01&&SF01) histos["Mll_invertdPhi_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhi_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhi_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertdPhi_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhi_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhi_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[13]==2)||(SFOS[14]==2)){
	  if(OS01&&SF01) histos["Mll_invertdPhi_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhi_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhi_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertdPhi_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertdPhi_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertdPhi_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[15]==1)||(SFOS[16]==1)){
	  if(OS01&&SF01) histos["Mll_invertPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertPt_1SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if((SFOS[15]==2)||(SFOS[16]==2)){
	  if(OS01&&SF01) histos["Mll_invertPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertPt_2SFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS01&&SF01) histos["Mll_invertPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02) histos["Mll_invertPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12) histos["Mll_invertPt_allSFOS_"+sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(SFOS[20]>=0&&SFOS[20]<=2) histos["MET_CRloose_allSFOS_"        +sn2]->Fill(met_pt(),weight*lepSF3l);
	if(SFOS[20]>=0&&SFOS[20]<=2) {
	  if(OS01&&SF01)             histos["MSFOS_CRloose_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02)             histos["MSFOS_CRloose_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12)             histos["MSFOS_CRloose_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(SFOS[20]>=0&&SFOS[20]<=2) histos["pTlll_CRloose_allSFOS_"      +sn2]->Fill(pTlll,weight*lepSF3l);
	if(SFOS[20]>=0&&SFOS[20]<=2) histos["dPhiMETlll_CRloose_allSFOS_" +sn2]->Fill(DPhilllMET,weight*lepSF3l);
	if(SFOS[18]==0&&pTlll>60.&&DPhilllMET<2.7)    histos["MET_CRlike_allSFOS_"        +sn2]->Fill(met_pt(),weight*lepSF3l);
	if(SFOS[18]==1&&pTlll>60.&&DPhilllMET<2.7)    histos["MET_CRlike_allSFOS_"        +sn2]->Fill(met_pt(),weight*lepSF3l);
	if(SFOS[18]==2&&pTlll>60.&&DPhilllMET<2.7)    histos["MET_CRlike_allSFOS_"        +sn2]->Fill(met_pt(),weight*lepSF3l);
	if(SFOS[18]==0&&pTlll>60.&&DPhilllMET<2.7) {
	  if(OS01&&SF01)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(SFOS[18]==1&&pTlll>60.&&DPhilllMET<2.7&&met_pt()>45.) {
	  if(OS01&&SF01)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(SFOS[18]==2&&pTlll>60.&&DPhilllMET<2.7&&met_pt()>55.) {
	  if(OS01&&SF01)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(),weight*lepSF3l);
	  if(OS02&&SF02)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	  if(OS12&&SF12)                              histos["MSFOS_CRlike_allSFOS_"      +sn2]->Fill((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(),weight*lepSF3l);
	}
	if(SFOS[18]==0&&DPhilllMET<2.7              ) histos["pTlll_CRlike_allSFOS_"      +sn2]->Fill(pTlll,weight*lepSF3l);
	if(SFOS[18]==1&&DPhilllMET<2.7&&met_pt()>45.) histos["pTlll_CRlike_allSFOS_"      +sn2]->Fill(pTlll,weight*lepSF3l);
	if(SFOS[18]==2&&DPhilllMET<2.7&&met_pt()>55.) histos["pTlll_CRlike_allSFOS_"      +sn2]->Fill(pTlll,weight*lepSF3l);
	if(SFOS[18]==0&&pTlll>60.              )      histos["dPhiMETlll_CRlike_allSFOS_" +sn2]->Fill(DPhilllMET,weight*lepSF3l);
	if(SFOS[18]==1&&pTlll>60.&&met_pt()>45.)      histos["dPhiMETlll_CRlike_allSFOS_" +sn2]->Fill(DPhilllMET,weight*lepSF3l);
	if(SFOS[18]==2&&pTlll>60.&&met_pt()>55.)      histos["dPhiMETlll_CRlike_allSFOS_" +sn2]->Fill(DPhilllMET,weight*lepSF3l);
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
  
  string filename = "rootfiles/Check3lCRv2.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //string mapname = histonames[i]+"_"+skimFilePrefix;
    //string writename = histonames[i];
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;

  
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
