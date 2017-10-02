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
  //this one has no mm
  histonames.push_back("SRpresel_MjjW_lMET");           hbins.push_back( 3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("SR_MjjW_lMET");                 hbins.push_back( 3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("Mll_SRpresel_ee_MjjW_lMET");    hbins.push_back(12); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_MjjW_em_lMET");             hbins.push_back(12); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("DEtall_SRpresel_ee_MjjW_lMET"); hbins.push_back(12); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("DEtall_MjjW_em_lMET");          hbins.push_back(12); hlow.push_back(0); hup.push_back(4);
  histonames.push_back("EPt_SRpresel_ee_MjjW_lMET");    hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("EPt_MjjW_em_lMET");             hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("NJ_SRpresel_MjjW_lMET");        hbins.push_back( 5); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("NJ_SR_MjjW_lMET");              hbins.push_back( 5); hlow.push_back(0); hup.push_back(5);
  //do three plots all in Mjj sideband, and with Mll>10 GeV - low MET, high MET, incl

  histonames.push_back("SR_lMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("SR_hMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("SR_iMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("SRpresel_lMET");   hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("SRpresel_hMET");   hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("SRpresel_iMET");   hbins.push_back(3); hlow.push_back(0); hup.push_back(3);

  histonames.push_back("SR_ee");         hbins.push_back(7); hlow.push_back(0); hup.push_back(7);
  histonames.push_back("SR_em");         hbins.push_back(7); hlow.push_back(0); hup.push_back(7);
  histonames.push_back("SR_mm");         hbins.push_back(7); hlow.push_back(0); hup.push_back(7);
  histonames.push_back("SRpresel_ee");   hbins.push_back(7); hlow.push_back(0); hup.push_back(7);
  histonames.push_back("SRpresel_em");   hbins.push_back(7); hlow.push_back(0); hup.push_back(7);
  histonames.push_back("SRpresel_mm");   hbins.push_back(7); hlow.push_back(0); hup.push_back(7);
   
  histonames.push_back("NB_ee_lMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_ee_hMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_ee_iMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_em_lMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_em_hMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_em_iMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_mm_lMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_mm_hMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);
  histonames.push_back("NB_mm_iMET");         hbins.push_back(3); hlow.push_back(0); hup.push_back(3);

  histonames.push_back("MET_ee_0b");    hbins.push_back(12); hlow.push_back(0); hup.push_back(120);
  histonames.push_back("MET_em_0b");    hbins.push_back(12); hlow.push_back(0); hup.push_back(120);
  histonames.push_back("MET_mm_0b");    hbins.push_back(12); hlow.push_back(0); hup.push_back(120);
  histonames.push_back("MTmax_ee_0b");  hbins.push_back(20); hlow.push_back(0); hup.push_back(200);
  histonames.push_back("MTmax_em_0b");  hbins.push_back(20); hlow.push_back(0); hup.push_back(200);
  histonames.push_back("MTmax_mm_0b");  hbins.push_back(20); hlow.push_back(0); hup.push_back(200);

  histonames.push_back("q_ee_0b_lMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_ee_0b_hMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_ee_0b_iMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_em_0b_lMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_em_0b_hMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_em_0b_iMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_mm_0b_lMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_mm_0b_hMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);
  histonames.push_back("q_mm_0b_iMET");         hbins.push_back(3); hlow.push_back(-1); hup.push_back(2);

  histonames.push_back("NJ_ee_0b_lMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_ee_0b_hMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_ee_0b_iMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_em_0b_lMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_em_0b_hMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_em_0b_iMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_mm_0b_lMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_mm_0b_hMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  histonames.push_back("NJ_mm_0b_iMET");         hbins.push_back(8); hlow.push_back(0); hup.push_back(8);
  
  histonames.push_back("Mll_ee_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_ee_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_ee_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_em_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_em_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_em_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_mm_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_mm_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mll_mm_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);

  histonames.push_back("Mjj_ee_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_ee_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_ee_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_em_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_em_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_em_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_mm_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_mm_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("Mjj_mm_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);

  histonames.push_back("MjjL_ee_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_ee_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_ee_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_em_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_em_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_em_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_mm_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_mm_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);
  histonames.push_back("MjjL_mm_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(450);

  histonames.push_back("pTll_ee_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_ee_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_ee_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_em_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_em_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_em_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_mm_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_mm_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);
  histonames.push_back("pTll_mm_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(150);

  histonames.push_back("dPhillMET_ee_0b_lMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_ee_0b_hMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_ee_0b_iMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_em_0b_lMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_em_0b_hMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_em_0b_iMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_mm_0b_lMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_mm_0b_hMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhillMET_mm_0b_iMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);

  histonames.push_back("dPhill_ee_0b_lMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_ee_0b_hMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_ee_0b_iMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_em_0b_lMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_em_0b_hMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_em_0b_iMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_mm_0b_lMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_mm_0b_hMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);
  histonames.push_back("dPhill_mm_0b_iMET");         hbins.push_back(16); hlow.push_back(0); hup.push_back(3.2);

  histonames.push_back("dEtall_ee_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_ee_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_ee_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_em_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_em_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_em_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_mm_0b_lMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_mm_0b_hMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);
  histonames.push_back("dEtall_mm_0b_iMET");         hbins.push_back(15); hlow.push_back(0); hup.push_back(4.5);

  histonames.push_back("minDRlj_ee_0b_lMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_ee_0b_hMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_ee_0b_iMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_em_0b_lMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_em_0b_hMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_em_0b_iMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_mm_0b_lMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_mm_0b_hMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);
  histonames.push_back("minDRlj_mm_0b_iMET");         hbins.push_back(12); hlow.push_back(0); hup.push_back(2.4);

  histonames.push_back("MuPt_em_0b_lMET");          hbins.push_back(12); hlow.push_back(30); hup.push_back(150);
  histonames.push_back("MuPt_em_0b_hMET");          hbins.push_back(12); hlow.push_back(30); hup.push_back(150);
  histonames.push_back("MuPt_em_0b_iMET");          hbins.push_back(12); hlow.push_back(30); hup.push_back(150);
  histonames.push_back("MuPt_mm_0b_lMET");          hbins.push_back(12); hlow.push_back(30); hup.push_back(150);
  histonames.push_back("MuPt_mm_0b_hMET");          hbins.push_back(12); hlow.push_back(30); hup.push_back(150);
  histonames.push_back("MuPt_mm_0b_iMET");          hbins.push_back(12); hlow.push_back(30); hup.push_back(150);

  histonames.push_back("MuEta_em_0b_lMET");         hbins.push_back(13); hlow.push_back(-2.6); hup.push_back(2.6);
  histonames.push_back("MuEta_em_0b_hMET");         hbins.push_back(13); hlow.push_back(-2.6); hup.push_back(2.6);
  histonames.push_back("MuEta_em_0b_iMET");         hbins.push_back(13); hlow.push_back(-2.6); hup.push_back(2.6);
  histonames.push_back("MuEta_mm_0b_lMET");         hbins.push_back(13); hlow.push_back(-2.6); hup.push_back(2.6);
  histonames.push_back("MuEta_mm_0b_hMET");         hbins.push_back(13); hlow.push_back(-2.6); hup.push_back(2.6);
  histonames.push_back("MuEta_mm_0b_iMET");         hbins.push_back(13); hlow.push_back(-2.6); hup.push_back(2.6);

  histonames.push_back("MuPhi_em_0b_lMET");         hbins.push_back(16); hlow.push_back(-3.2); hup.push_back(3.2);
  histonames.push_back("MuPhi_em_0b_hMET");         hbins.push_back(16); hlow.push_back(-3.2); hup.push_back(3.2);
  histonames.push_back("MuPhi_em_0b_iMET");         hbins.push_back(16); hlow.push_back(-3.2); hup.push_back(3.2);
  histonames.push_back("MuPhi_mm_0b_lMET");         hbins.push_back(16); hlow.push_back(-3.2); hup.push_back(3.2);
  histonames.push_back("MuPhi_mm_0b_hMET");         hbins.push_back(16); hlow.push_back(-3.2); hup.push_back(3.2);
  histonames.push_back("MuPhi_mm_0b_iMET");         hbins.push_back(16); hlow.push_back(-3.2); hup.push_back(3.2);

  histonames.push_back("MuRelIso_em_0b_lMET");      hbins.push_back(10); hlow.push_back(0.); hup.push_back(0.1);
  histonames.push_back("MuRelIso_em_0b_hMET");      hbins.push_back(10); hlow.push_back(0.); hup.push_back(0.1);
  histonames.push_back("MuRelIso_em_0b_iMET");      hbins.push_back(10); hlow.push_back(0.); hup.push_back(0.1);
  histonames.push_back("MuRelIso_mm_0b_lMET");      hbins.push_back(10); hlow.push_back(0.); hup.push_back(0.1);
  histonames.push_back("MuRelIso_mm_0b_hMET");      hbins.push_back(10); hlow.push_back(0.); hup.push_back(0.1);
  histonames.push_back("MuRelIso_mm_0b_iMET");      hbins.push_back(10); hlow.push_back(0.); hup.push_back(0.1);

  histonames.push_back("MuIP3D_em_0b_lMET");        hbins.push_back(15); hlow.push_back(0.); hup.push_back(0.015);
  histonames.push_back("MuIP3D_em_0b_hMET");        hbins.push_back(15); hlow.push_back(0.); hup.push_back(0.015);
  histonames.push_back("MuIP3D_em_0b_iMET");        hbins.push_back(15); hlow.push_back(0.); hup.push_back(0.015);
  histonames.push_back("MuIP3D_mm_0b_lMET");        hbins.push_back(15); hlow.push_back(0.); hup.push_back(0.015);
  histonames.push_back("MuIP3D_mm_0b_hMET");        hbins.push_back(15); hlow.push_back(0.); hup.push_back(0.015);
  histonames.push_back("MuIP3D_mm_0b_iMET");        hbins.push_back(15); hlow.push_back(0.); hup.push_back(0.015);

  histonames.push_back("MuPtErrOvPt_em_0b_lMET");   hbins.push_back(15); hlow.push_back(-0.); hup.push_back(0.3);
  histonames.push_back("MuPtErrOvPt_em_0b_hMET");   hbins.push_back(15); hlow.push_back(-0.); hup.push_back(0.3);
  histonames.push_back("MuPtErrOvPt_em_0b_iMET");   hbins.push_back(15); hlow.push_back(-0.); hup.push_back(0.3);
  histonames.push_back("MuPtErrOvPt_mm_0b_lMET");   hbins.push_back(15); hlow.push_back(-0.); hup.push_back(0.3);
  histonames.push_back("MuPtErrOvPt_mm_0b_hMET");   hbins.push_back(15); hlow.push_back(-0.); hup.push_back(0.3);
  histonames.push_back("MuPtErrOvPt_mm_0b_iMET");   hbins.push_back(15); hlow.push_back(-0.); hup.push_back(0.3);

  histonames.push_back("MuValidFrac_em_0b_lMET");   hbins.push_back(18); hlow.push_back( 0.4); hup.push_back(1);
  histonames.push_back("MuValidFrac_em_0b_hMET");   hbins.push_back(18); hlow.push_back( 0.4); hup.push_back(1);
  histonames.push_back("MuValidFrac_em_0b_iMET");   hbins.push_back(18); hlow.push_back( 0.4); hup.push_back(1);
  histonames.push_back("MuValidFrac_mm_0b_lMET");   hbins.push_back(18); hlow.push_back( 0.4); hup.push_back(1);
  histonames.push_back("MuValidFrac_mm_0b_hMET");   hbins.push_back(18); hlow.push_back( 0.4); hup.push_back(1);
  histonames.push_back("MuValidFrac_mm_0b_iMET");   hbins.push_back(18); hlow.push_back( 0.4); hup.push_back(1);

  histonames.push_back("MuChi2N_em_0b_lMET");       hbins.push_back(20); hlow.push_back(   0); hup.push_back(6);
  histonames.push_back("MuChi2N_em_0b_hMET");       hbins.push_back(20); hlow.push_back(   0); hup.push_back(6);
  histonames.push_back("MuChi2N_em_0b_iMET");       hbins.push_back(20); hlow.push_back(   0); hup.push_back(6);
  histonames.push_back("MuChi2N_mm_0b_lMET");       hbins.push_back(20); hlow.push_back(   0); hup.push_back(6);
  histonames.push_back("MuChi2N_mm_0b_hMET");       hbins.push_back(20); hlow.push_back(   0); hup.push_back(6);
  histonames.push_back("MuChi2N_mm_0b_iMET");       hbins.push_back(20); hlow.push_back(   0); hup.push_back(6);

  histonames.push_back("MuLostHits_em_0b_lMET");    hbins.push_back( 3); hlow.push_back(   0); hup.push_back(3);
  histonames.push_back("MuLostHits_em_0b_hMET");    hbins.push_back( 3); hlow.push_back(   0); hup.push_back(3);
  histonames.push_back("MuLostHits_em_0b_iMET");    hbins.push_back( 3); hlow.push_back(   0); hup.push_back(3);
  histonames.push_back("MuLostHits_mm_0b_lMET");    hbins.push_back( 3); hlow.push_back(   0); hup.push_back(3);
  histonames.push_back("MuLostHits_mm_0b_hMET");    hbins.push_back( 3); hlow.push_back(   0); hup.push_back(3);
  histonames.push_back("MuLostHits_mm_0b_iMET");    hbins.push_back( 3); hlow.push_back(   0); hup.push_back(3);

    
  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    string mapname = histonames[i];
    if(skimFilePrefix.find("Data")!=string::npos){
	mapname = histonames[i] + "_Data";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_FakePred";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
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
      vector<int> vSS, v3l, v, iSS, i3l;
      vector<int> vaSS, va3l, va, iaSS, ia3l;
      int looseEle = -1; int veton3lspec = 0;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool isSS = false; bool is3l = false;
	bool isaSS = false; bool isa3l = false;
	bool islSS = false; bool isl3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015){
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

	if(fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015){
	  if(lep_pass_VVV_cutbased_fo_noiso()[i]&&abs(lep_pdgId()[i])==11&&lep_lostHits()[i]==0){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.2){
	      if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.2){
	      if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	    }
	  }
	  if(lep_pass_VVV_cutbased_fo_noiso()[i]&&abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.4){
	    if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	    if(lep_p4()[i ].Pt()>30&&!isSS) { iaSS.push_back(i); isaSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_p4()[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4) {
	  v.push_back(i);
	  if(!isSS) vSS.push_back(i);
	  if(!is3l) {
	    v3l.push_back(i);
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

      if(nmu>=2)           passofflineforTrigger = true;
      if(nmu25>=1&&nel>=1) passofflineforTrigger = true;
      if(nel25>=1&&nel>=2) passofflineforTrigger = true;
      if((nSS+naSS)>=2)    passofflineforTrigger = true;

      if(n3l<2&&(nSS+naSS)<2) continue;
      //if(nj30<2) continue;
      //if(nb!=0) continue;
      if(!passofflineforTrigger) continue;
      if((nSS+naSS)==0&&lep_p4()[i3l[0] ].Pt()<25) continue;
      
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
	if((iSS.size()+iaSS.size())>=2){
	  int l1(-1),l2(-1);
	  if(iSS.size()>=2) { l1 = iSS[0]; l2 = iSS[1]; }
	  else if(iSS.size()==1&&iaSS.size()>=1) { l1 = iSS[0]; l2 = iaSS[0]; }
	  else if(iaSS.size()>=2) { l1 = iaSS[0]; l2 = iaSS[1]; }
	  int gentype = gentype_v2(l1,l2,-1);
	  if(     gentype==0) sn = "trueSS";
	  else if(gentype==2) sn = "chargeflips";
	  else if(gentype==3) sn = "SSLL";
	  else if(gentype==4) sn = "fakes";
	  else if(gentype==5) sn = "photonfakes";
	  else                sn = "others";
	}
	if((i3l.size()>=3)||(i3l.size()==2&&looseEle>=0)){
	  int l1(-1), l2(-1), l3(-1);
	  if(i3l.size()>=3) { l1 = i3l[0]; l2 = i3l[1]; l3 = i3l[2]; }
	  else if(i3l.size()>=2&&looseEle>=0) { l1 = i3l[0]; l2 = i3l[1]; l3 = looseEle; }
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
      float Mjjmu1(-1), Mjjmu2(-2);
      if(iSS.size()>=2&&nj30>=2){
	Mjjmu1 = (jets_p4()[i2p5[0] ]+jets_p4()[i2p5[1] ]+lep_p4()[iSS[0] ]).M();
	Mjjmu2 = (jets_p4()[i2p5[0] ]+jets_p4()[i2p5[1] ]+lep_p4()[iSS[1] ]).M();
      }
      float Mll = -1;
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
      if(nj30>=2/*&&nb==0*/&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>10.){
	Mll = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) ee[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11) em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13) em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) mm[1] = true;
	//if(ee[1]==true&&fabs(Mll-90.)<10.) Mll = -1.;
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
      if(nj30>=2/*&&nb==0*/&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>10.){
	Mll = (lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M();
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11) ee[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11) em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13) em[3] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13) mm[3] = true;
	//if(ee[3]==true&&fabs(Mll-90.)<10.) Mll = -1.;
      }
      //4: AR-loose-loose
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&naSS==2&&((lep_pdgId()[iaSS[1] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iaSS[1] ]+lep_p4()[iaSS[0] ]).M()>30.){
	if(abs(lep_pdgId()[iaSS[1] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&fabs((lep_p4()[iaSS[1] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) ee[4] = true;
	if(abs(lep_pdgId()[iaSS[1] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&aMTmax>90.)                                                em[4] = true;
	if(abs(lep_pdgId()[iaSS[1] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13&&met_pt()>40.&&aMTmax>90.)                                                em[4] = true;
	if(abs(lep_pdgId()[iaSS[1] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13)                                                                          mm[4] = true;
	if((lep_p4()[iaSS[1] ]+lep_p4()[iaSS[0] ]).M()<40){ ee[4] = false; mm[4] = false; }
      }
      //5: ARpreselect-loose-loose
      if(nj30>=2/*&&nb==0*/&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&naSS==2&&((lep_pdgId()[iaSS[1] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iaSS[1] ]+lep_p4()[iaSS[0] ]).M()>10.){
	Mll = (lep_p4()[iaSS[0] ]+lep_p4()[iaSS[1] ]).M();
	if(abs(lep_pdgId()[iaSS[1] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11) ee[5] = true;
	if(abs(lep_pdgId()[iaSS[1] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11) em[5] = true;
	if(abs(lep_pdgId()[iaSS[1] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13) em[5] = true;
	if(abs(lep_pdgId()[iaSS[1] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13) mm[5] = true;
	//if(ee[5]==true&&fabs(Mll-90.)<10.) Mll = -1.;
      }
      //6: SRpreselect - w/o NJ
      if(/*nj30>=2&&nb==0&&*/nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>10.){
	Mll = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11) ee[6] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11) em[6] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13) em[6] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) mm[6] = true;
	//if(ee[1]==true&&fabs(Mll-90.)<10.) Mll = -1.;
      }
      //7: SR no Mjj and no NB
      if(nj30>=2/*&&nb==0*/&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[7] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&MTmax>90.)                                               em[7] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()>40.&&MTmax>90.)                                               em[7] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                        mm[7] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[7] = false; mm[7] = false; }
	if(fabs(Detajj)>1.5) { ee[7] = false; em[7] = false; mm[7] = false; }
	if(fabs(MjjL)>400.)  { ee[7] = false; em[7] = false; mm[7] = false; }
	//if(fabs(Mjj-80.)<20.){ ee[8] = false; em[8] = false; mm[8] = false; }
      }
      //8: AR no Mjj and no NB
      if(nj30>=2/*&&nb==0*/&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&naSS==1&&nSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) ee[8] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&aMTmax>90.)                                               em[8] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13&&met_pt()>40.&&aMTmax>90.)                                               em[8] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13)                                                                         mm[8] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()<40){ ee[8] = false; mm[8] = false; }
	if(fabs(Detajj)>1.5) { ee[8] = false; em[8] = false; mm[8] = false; }
	if(fabs(MjjL)>400.)  { ee[8] = false; em[8] = false; mm[8] = false; }
	//if(fabs(Mjj-80.)<20.){ ee[8] = false; em[8] = false; mm[8] = false; }
      }
      //9: ARpreselect - w/o NJ
      if(/*nj30>=2&&nb==0&&*/nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>10.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11) ee[9] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11) em[9] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13) em[9] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==13) mm[9] = true;
	//if(ee[1]==true&&fabs(Mll-90.)<10.) Mll = -1.;
      }
      for(int i = 0; i<50; ++i){
	if(nSS>=1){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[0] ])==13&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[0] ])==13&&lep_motherIdSS()[iSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if(nSS>=2){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iSS[1] ])==13&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iSS[1] ])==13&&lep_motherIdSS()[iSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if(naSS>=1){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iaSS[0] ])==11&&lep_motherIdSS()[iaSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iaSS[0] ])==13&&lep_motherIdSS()[iaSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iaSS[0] ])==11&&lep_motherIdSS()[iaSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iaSS[0] ])==13&&lep_motherIdSS()[iaSS[0] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	}
	if(naSS>=2){
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iaSS[1] ])==11&&lep_motherIdSS()[iaSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[iaSS[1] ])==13&&lep_motherIdSS()[iaSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iaSS[1] ])==11&&lep_motherIdSS()[iaSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[iaSS[1] ])==13&&lep_motherIdSS()[iaSS[1] ]==(-3)) { ee[i] = false; em[i] = false; mm[i] = false; }
	}
      }

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
	  if((DPhilllMET>2.7&&pTlll>60.)&&passlowMSFOS)               SFOS[0] = 0;
	  if(!(DPhilllMET>2.7&&pTlll>60.))               SFOS[2] = 0;
	}
	if(pass1){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)&&passlowMSFOS&&passlowMlll) SFOS[0] = 1;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[2] = 1;
	}
	if(pass2){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)&&passlowMSFOS&&passlowMlll) SFOS[0] = 2;
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

      if(nb==0&&sn=="photonfakes"){
	if(ee[7]){
	  if(fabs(Mjj-80.)<20.) cout << "SR ee photon fake from " << fname << endl;
	  else                  cout << "Mjj sideband ee photon fake from " << fname << endl;
	}
	if(em[7]){
	  if(fabs(Mjj-80.)<20.) cout << "SR em photon fake from " << fname << endl;
	  else                  cout << "Mjj sideband em photon fake from " << fname << endl;
	}
	if(mm[7]){
	  if(fabs(Mjj-80.)<20.) cout << "SR em photon fake from " << fname << endl;
	  else                  cout << "Mjj sideband em photon fake from " << fname << endl;
	}
	if(SFOS[0]==0) cout << "SR 0SFOS photon fake from " << fname << endl;
	if(SFOS[0]==1) cout << "SR 1SFOS photon fake from " << fname << endl;
	if(SFOS[0]==2) cout << "SR 2SFOS photon fake from " << fname << endl;
      }
	

      //SR:0
      //SR w/o NB, Mjj: 7
      //SRpresel: 1 (no nb)
      //SRpresel + no nj cut: 6
      if(fabs(Mjj-80.)<20.) continue;

      
      if(!isData()){
	//cout << __LINE__ << endl;
	if(ee[7]&&nb==0){
	  histos["SR_ee_"+sn ]->Fill(0.,weight*5.788348400057/35.867059982507);
	  histos["SR_ee_"+sn ]->Fill(1.,weight*2.573399420069/35.867059982507);
	  histos["SR_ee_"+sn ]->Fill(2.,weight*4.248383597366/35.867059982507);
	  histos["SR_ee_"+sn ]->Fill(3.,weight*4.009132485903/35.867059982507);
	  histos["SR_ee_"+sn ]->Fill(4.,weight*3.101618425944/35.867059982507);
	  histos["SR_ee_"+sn ]->Fill(5.,weight*7.540487746602/35.867059982507);
	  histos["SR_ee_"+sn ]->Fill(6.,weight*8.605689906566/35.867059982507);
	}
	if(em[7]&&nb==0){
	  histos["SR_em_"+sn ]->Fill(0.,weight*5.788348400057/35.867059982507);
	  histos["SR_em_"+sn ]->Fill(1.,weight*2.573399420069/35.867059982507);
	  histos["SR_em_"+sn ]->Fill(2.,weight*4.248383597366/35.867059982507);
	  histos["SR_em_"+sn ]->Fill(3.,weight*4.009132485903/35.867059982507);
	  histos["SR_em_"+sn ]->Fill(4.,weight*3.101618425944/35.867059982507);
	  histos["SR_em_"+sn ]->Fill(5.,weight*7.540487746602/35.867059982507);
	  histos["SR_em_"+sn ]->Fill(6.,weight*8.605689906566/35.867059982507);
	}
	if(mm[7]&&nb==0){
	  histos["SR_mm_"+sn ]->Fill(0.,weight*5.788348400057/35.867059982507);
	  histos["SR_mm_"+sn ]->Fill(1.,weight*2.573399420069/35.867059982507);
	  histos["SR_mm_"+sn ]->Fill(2.,weight*4.248383597366/35.867059982507);
	  histos["SR_mm_"+sn ]->Fill(3.,weight*4.009132485903/35.867059982507);
	  histos["SR_mm_"+sn ]->Fill(4.,weight*3.101618425944/35.867059982507);
	  histos["SR_mm_"+sn ]->Fill(5.,weight*7.540487746602/35.867059982507);
	  histos["SR_mm_"+sn ]->Fill(6.,weight*8.605689906566/35.867059982507);
	}
	//cout << __LINE__ << endl;
	if(ee[1]&&nb==0){
	  histos["SRpresel_ee_"+sn ]->Fill(0.,weight*5.788348400057/35.867059982507);
	  histos["SRpresel_ee_"+sn ]->Fill(1.,weight*2.573399420069/35.867059982507);
	  histos["SRpresel_ee_"+sn ]->Fill(2.,weight*4.248383597366/35.867059982507);
	  histos["SRpresel_ee_"+sn ]->Fill(3.,weight*4.009132485903/35.867059982507);
	  histos["SRpresel_ee_"+sn ]->Fill(4.,weight*3.101618425944/35.867059982507);
	  histos["SRpresel_ee_"+sn ]->Fill(5.,weight*7.540487746602/35.867059982507);
	  histos["SRpresel_ee_"+sn ]->Fill(6.,weight*8.605689906566/35.867059982507);
	}
	if(em[1]&&nb==0){
	  histos["SRpresel_em_"+sn ]->Fill(0.,weight*5.788348400057/35.867059982507);
	  histos["SRpresel_em_"+sn ]->Fill(1.,weight*2.573399420069/35.867059982507);
	  histos["SRpresel_em_"+sn ]->Fill(2.,weight*4.248383597366/35.867059982507);
	  histos["SRpresel_em_"+sn ]->Fill(3.,weight*4.009132485903/35.867059982507);
	  histos["SRpresel_em_"+sn ]->Fill(4.,weight*3.101618425944/35.867059982507);
	  histos["SRpresel_em_"+sn ]->Fill(5.,weight*7.540487746602/35.867059982507);
	  histos["SRpresel_em_"+sn ]->Fill(6.,weight*8.605689906566/35.867059982507);
	}
	if(mm[1]&&nb==0){
	  histos["SRpresel_mm_"+sn ]->Fill(0.,weight*5.788348400057/35.867059982507);
	  histos["SRpresel_mm_"+sn ]->Fill(1.,weight*2.573399420069/35.867059982507);
	  histos["SRpresel_mm_"+sn ]->Fill(2.,weight*4.248383597366/35.867059982507);
	  histos["SRpresel_mm_"+sn ]->Fill(3.,weight*4.009132485903/35.867059982507);
	  histos["SRpresel_mm_"+sn ]->Fill(4.,weight*3.101618425944/35.867059982507);
	  histos["SRpresel_mm_"+sn ]->Fill(5.,weight*7.540487746602/35.867059982507);
	  histos["SRpresel_mm_"+sn ]->Fill(6.,weight*8.605689906566/35.867059982507);
	}
	//cout << __LINE__ << endl;
      } else { //data
	if(fname.find("data_Run2016B")!=string::npos&&nb==0){
	  if(ee[7]     ) histos["SR_ee_" +sn ]->Fill(0.,1.);
	  if(em[7]     ) histos["SR_em_" +sn ]->Fill(0.,1.);
	  if(mm[7]     ) histos["SR_mm_" +sn ]->Fill(0.,1.);
	  if(ee[1]     ) histos["SRpresel_ee_" +sn ]->Fill(0.,1.);
	  if(em[1]     ) histos["SRpresel_em_" +sn ]->Fill(0.,1.);
	  if(mm[1]     ) histos["SRpresel_mm_" +sn ]->Fill(0.,1.);
	}
	if(fname.find("data_Run2016C")!=string::npos&&nb==0){
	  if(ee[7]     ) histos["SR_ee_" +sn ]->Fill(1.,1.);
	  if(em[7]     ) histos["SR_em_" +sn ]->Fill(1.,1.);
	  if(mm[7]     ) histos["SR_mm_" +sn ]->Fill(1.,1.);
	  if(ee[1]     ) histos["SRpresel_ee_" +sn ]->Fill(1.,1.);
	  if(em[1]     ) histos["SRpresel_em_" +sn ]->Fill(1.,1.);
	  if(mm[1]     ) histos["SRpresel_mm_" +sn ]->Fill(1.,1.);
	}
	if(fname.find("data_Run2016D")!=string::npos&&nb==0){
	  if(ee[7]     ) histos["SR_ee_" +sn ]->Fill(2.,1.);
	  if(em[7]     ) histos["SR_em_" +sn ]->Fill(2.,1.);
	  if(mm[7]     ) histos["SR_mm_" +sn ]->Fill(2.,1.);
	  if(ee[1]     ) histos["SRpresel_ee_" +sn ]->Fill(2.,1.);
	  if(em[1]     ) histos["SRpresel_em_" +sn ]->Fill(2.,1.);
	  if(mm[1]     ) histos["SRpresel_mm_" +sn ]->Fill(2.,1.);
	}
	if(fname.find("data_Run2016E")!=string::npos&&nb==0){
	  if(ee[7]     ) histos["SR_ee_" +sn ]->Fill(3.,1.);
	  if(em[7]     ) histos["SR_em_" +sn ]->Fill(3.,1.);
	  if(mm[7]     ) histos["SR_mm_" +sn ]->Fill(3.,1.);
	  if(ee[1]     ) histos["SRpresel_ee_" +sn ]->Fill(3.,1.);
	  if(em[1]     ) histos["SRpresel_em_" +sn ]->Fill(3.,1.);
	  if(mm[1]     ) histos["SRpresel_mm_" +sn ]->Fill(3.,1.);
	}
	if(fname.find("data_Run2016F")!=string::npos&&nb==0){
	  if(ee[7]     ) histos["SR_ee_" +sn ]->Fill(4.,1.);
	  if(em[7]     ) histos["SR_em_" +sn ]->Fill(4.,1.);
	  if(mm[7]     ) histos["SR_mm_" +sn ]->Fill(4.,1.);
	  if(ee[1]     ) histos["SRpresel_ee_" +sn ]->Fill(4.,1.);
	  if(em[1]     ) histos["SRpresel_em_" +sn ]->Fill(4.,1.);
	  if(mm[1]     ) histos["SRpresel_mm_" +sn ]->Fill(4.,1.);
	}
	if(fname.find("data_Run2016G")!=string::npos&&nb==0){
	  if(ee[7]     ) histos["SR_ee_" +sn ]->Fill(5.,1.);
	  if(em[7]     ) histos["SR_em_" +sn ]->Fill(5.,1.);
	  if(mm[7]     ) histos["SR_mm_" +sn ]->Fill(5.,1.);
	  if(ee[1]     ) histos["SRpresel_ee_" +sn ]->Fill(5.,1.);
	  if(em[1]     ) histos["SRpresel_em_" +sn ]->Fill(5.,1.);
	  if(mm[1]     ) histos["SRpresel_mm_" +sn ]->Fill(5.,1.);
	}
	if(fname.find("data_Run2016H")!=string::npos&&nb==0){
	  if(ee[7]     ) histos["SR_ee_" +sn ]->Fill(6.,1.);
	  if(em[7]     ) histos["SR_em_" +sn ]->Fill(6.,1.);
	  if(mm[7]     ) histos["SR_mm_" +sn ]->Fill(6.,1.);
	  if(ee[1]     ) histos["SRpresel_ee_" +sn ]->Fill(6.,1.);
	  if(em[1]     ) histos["SRpresel_em_" +sn ]->Fill(6.,1.);
	  if(mm[1]     ) histos["SRpresel_mm_" +sn ]->Fill(6.,1.);
	}
      }
      //cout << __LINE__ << endl;
      if(ee[7]&&nb==0     ) histos["SR_iMET_"+sn ]->Fill(0.,weight);
      if(em[7]&&nb==0     ) histos["SR_iMET_"+sn ]->Fill(1.,weight);
      if(mm[7]&&nb==0     ) histos["SR_iMET_"+sn ]->Fill(2.,weight);
      //cout << __LINE__ << endl;
	if(ee[7]&&nb==0&&met_pt()<40     ) histos["SR_lMET_"+sn ]->Fill(0.,weight);
      if(em[7]&&nb==0&&met_pt()<40     ) histos["SR_lMET_"+sn ]->Fill(1.,weight);
      if(mm[7]&&nb==0&&met_pt()<40     ) histos["SR_lMET_"+sn ]->Fill(2.,weight);
      //cout << __LINE__ << endl;
      if(ee[7]&&nb==0&&met_pt()>40     ) histos["SR_hMET_"+sn ]->Fill(0.,weight);
      if(em[7]&&nb==0&&met_pt()>40     ) histos["SR_hMET_"+sn ]->Fill(1.,weight);
      if(mm[7]&&nb==0&&met_pt()>40     ) histos["SR_hMET_"+sn ]->Fill(2.,weight);
      //cout << __LINE__ << endl;
      if(ee[1]&&nb==0     ) histos["SRpresel_iMET_"+sn ]->Fill(0.,weight);
      if(em[1]&&nb==0     ) histos["SRpresel_iMET_"+sn ]->Fill(1.,weight);
      if(mm[1]&&nb==0     ) histos["SRpresel_iMET_"+sn ]->Fill(2.,weight);
      //cout << __LINE__ << endl;
      if(ee[1]&&nb==0&&met_pt()<40     ) histos["SRpresel_lMET_"+sn ]->Fill(0.,weight);
      if(em[1]&&nb==0&&met_pt()<40     ) histos["SRpresel_lMET_"+sn ]->Fill(1.,weight);
      if(mm[1]&&nb==0&&met_pt()<40     ) histos["SRpresel_lMET_"+sn ]->Fill(2.,weight);
      //cout << __LINE__ << endl;
      if(ee[1]&&nb==0&&met_pt()>40     ) histos["SRpresel_hMET_"+sn ]->Fill(0.,weight);
      if(em[1]&&nb==0&&met_pt()>40     ) histos["SRpresel_hMET_"+sn ]->Fill(1.,weight);
      if(mm[1]&&nb==0&&met_pt()>40     ) histos["SRpresel_hMET_"+sn ]->Fill(2.,weight);

      //cout << __LINE__ << endl;
      if(ee[6]&&nb==0     ) histos["NJ_ee_0b_iMET_"+sn ]->Fill(nj30,weight);
      if(em[6]&&nb==0     ) histos["NJ_em_0b_iMET_"+sn ]->Fill(nj30,weight);
      if(mm[6]&&nb==0     ) histos["NJ_mm_0b_iMET_"+sn ]->Fill(nj30,weight);
      //cout << __LINE__ << endl;
      if(ee[6]&&nb==0&&met_pt()<40     ) histos["NJ_ee_0b_lMET_"+sn ]->Fill(nj30,weight);
      if(em[6]&&nb==0&&met_pt()<40     ) histos["NJ_em_0b_lMET_"+sn ]->Fill(nj30,weight);
      if(mm[6]&&nb==0&&met_pt()<40     ) histos["NJ_mm_0b_lMET_"+sn ]->Fill(nj30,weight);
      //cout << __LINE__ << endl;
      if(ee[6]&&nb==0&&met_pt()>40     ) histos["NJ_ee_0b_hMET_"+sn ]->Fill(nj30,weight);
      if(em[6]&&nb==0&&met_pt()>40     ) histos["NJ_em_0b_hMET_"+sn ]->Fill(nj30,weight);
      if(mm[6]&&nb==0&&met_pt()>40     ) histos["NJ_mm_0b_hMET_"+sn ]->Fill(nj30,weight);

      //cout << __LINE__ << endl;
      if(ee[1]     ) histos["NB_ee_iMET_"+sn ]->Fill(nb,weight);
      if(em[1]     ) histos["NB_em_iMET_"+sn ]->Fill(nb,weight);
      if(mm[1]     ) histos["NB_mm_iMET_"+sn ]->Fill(nb,weight);
      //cout << __LINE__ << endl;
      if(ee[1]&&met_pt()<40     ) histos["NB_ee_lMET_"+sn ]->Fill(nb,weight);
      if(em[1]&&met_pt()<40     ) histos["NB_em_lMET_"+sn ]->Fill(nb,weight);
      if(mm[1]&&met_pt()<40     ) histos["NB_mm_lMET_"+sn ]->Fill(nb,weight);
      //cout << __LINE__ << endl;
      if(ee[1]&&met_pt()>40     ) histos["NB_ee_hMET_"+sn ]->Fill(nb,weight);
      if(em[1]&&met_pt()>40     ) histos["NB_em_hMET_"+sn ]->Fill(nb,weight);
      if(mm[1]&&met_pt()>40     ) histos["NB_mm_hMET_"+sn ]->Fill(nb,weight);

      float dRminSR(99.);//any jet pT > 20 GeV
      if(nb==0){
	//if(!(ee[1]||em[1]||mm[1])) continue;
	for(unsigned int n = 0; n<jets_csv().size();++n){
	  if(jets_p4()[n].Pt()<20.) continue;
	  if(fabs(jets_p4()[n].Eta())>5.0) continue;
	  if(ee[1]||em[1]||mm[1]) {
	    if(dR(jets_p4()[n],lep_p4()[ iSS[0] ])<dRminSR) dRminSR = dR(jets_p4()[n],lep_p4()[ iSS[0] ]);
	    if(dR(jets_p4()[n],lep_p4()[ iSS[1] ])<dRminSR) dRminSR = dR(jets_p4()[n],lep_p4()[ iSS[1] ]);
	  }
	  if(ee[9]||em[9]||mm[9]) {
	    if(dR(jets_p4()[n],lep_p4()[ iSS[0] ])<dRminSR) dRminSR = dR(jets_p4()[n],lep_p4()[ iSS[0] ]);
	    if(dR(jets_p4()[n],lep_p4()[iaSS[0] ])<dRminSR) dRminSR = dR(jets_p4()[n],lep_p4()[iaSS[0] ]);
	  }
	}

	if(ee[1]){
	  histos["MET_ee_0b_"+sn ]->Fill(met_pt(),weight);
	  histos["MTmax_ee_0b_"+sn ]->Fill(MTmax,weight);
	  if(lep_pdgId()[iSS[0] ]>0) histos["q_ee_0b_iMET_"+sn ]->Fill(-1.,weight);
	  else                       histos["q_ee_0b_iMET_"+sn ]->Fill( 1.,weight);
	  histos["Mll_ee_0b_iMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	  histos["Mjj_ee_0b_iMET_"+sn ]->Fill(Mjj,weight);
	  histos["MjjL_ee_0b_iMET_"+sn ]->Fill(MjjL,weight);
	  histos["pTll_ee_0b_iMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	  histos["dPhillMET_ee_0b_iMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	  histos["dPhill_ee_0b_iMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	  histos["dEtall_ee_0b_iMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	  histos["minDRlj_ee_0b_iMET_"+sn ]->Fill(dRminSR,weight);
	  if(met_pt()<40){
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_ee_0b_lMET_"+sn ]->Fill(-1.,weight);
	    else                       histos["q_ee_0b_lMET_"+sn ]->Fill( 1.,weight);
	    histos["Mll_ee_0b_lMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	    histos["Mjj_ee_0b_lMET_"+sn ]->Fill(Mjj,weight);
	    histos["MjjL_ee_0b_lMET_"+sn ]->Fill(MjjL,weight);
	    histos["pTll_ee_0b_lMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	    histos["dPhillMET_ee_0b_lMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	    histos["dPhill_ee_0b_lMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["dEtall_ee_0b_lMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["minDRlj_ee_0b_lMET_"+sn ]->Fill(dRminSR,weight);
	  }
	  else {
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_ee_0b_hMET_"+sn ]->Fill(-1.,weight);
	    else                       histos["q_ee_0b_hMET_"+sn ]->Fill( 1.,weight);
	    histos["Mll_ee_0b_hMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	    histos["Mjj_ee_0b_hMET_"+sn ]->Fill(Mjj,weight);
	    histos["MjjL_ee_0b_hMET_"+sn ]->Fill(MjjL,weight);
	    histos["pTll_ee_0b_hMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	    histos["dPhillMET_ee_0b_hMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	    histos["dPhill_ee_0b_hMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["dEtall_ee_0b_hMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["minDRlj_ee_0b_hMET_"+sn ]->Fill(dRminSR,weight);
	  }
	}//ee[1]
	if(em[1]){
	  histos["MET_em_0b_"+sn ]->Fill(met_pt(),weight);
	  histos["MTmax_em_0b_"+sn ]->Fill(MTmax,weight);
	  if(lep_pdgId()[iSS[0] ]>0) histos["q_em_0b_iMET_"+sn ]->Fill(-1.,weight);
	  else                       histos["q_em_0b_iMET_"+sn ]->Fill( 1.,weight);
	  histos["Mll_em_0b_iMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	  histos["Mjj_em_0b_iMET_"+sn ]->Fill(Mjj,weight);
	  histos["MjjL_em_0b_iMET_"+sn ]->Fill(MjjL,weight);
	  histos["pTll_em_0b_iMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	  histos["dPhillMET_em_0b_iMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	  histos["dPhill_em_0b_iMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	  histos["dEtall_em_0b_iMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	  histos["minDRlj_em_0b_iMET_"+sn ]->Fill(dRminSR,weight);
	  int l = -1;
	  if(abs(lep_pdgId()[iSS[0] ])==13) l = iSS[0];
	  else if(abs(lep_pdgId()[iSS[1] ])==13) l = iSS[1];
	  else cout << "WTF em " << __LINE__ << endl;
	  //cout << __LINE__ << endl;
	  histos["MuPt_em_0b_iMET_"+sn ]->Fill(lep_p4()[l].Pt(),weight);
	  //cout << __LINE__ << endl;
	  histos["MuPhi_em_0b_iMET_"+sn ]->Fill(lep_p4()[l].Phi(),weight);
	  //cout << __LINE__ << endl;
	  histos["MuEta_em_0b_iMET_"+sn ]->Fill(lep_p4()[l].Eta(),weight);
	  //cout << __LINE__ << endl;
	  histos["MuRelIso_em_0b_iMET_"+sn ]->Fill(lep_relIso03EA()[l],weight);
	  //cout << __LINE__ << endl;
	  histos["MuIP3D_em_0b_iMET_"+sn ]->Fill(lep_ip3d()[l],weight);
	  //cout << __LINE__ << endl;
	  histos["MuPtErrOvPt_em_0b_iMET_"+sn ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weight);
	  //cout << __LINE__ << endl;
	  histos["MuValidFrac_em_0b_iMET_"+sn ]->Fill(lep_validfraction()[l],weight);
	  //cout << __LINE__ << endl;
	  histos["MuChi2N_em_0b_iMET_"+sn ]->Fill(lep_x2ondof()[l],weight);
	  //cout << __LINE__ << endl;
	  histos["MuLostHits_em_0b_iMET_"+sn ]->Fill(lep_lostHits()[l],weight);
	  //cout << __LINE__ << endl;
	  if(met_pt()<40){
	    //cout << __LINE__ << endl;
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_em_0b_lMET_"+sn ]->Fill(-1.,weight);
	    else                       histos["q_em_0b_lMET_"+sn ]->Fill( 1.,weight);
	    histos["Mll_em_0b_lMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	    histos["Mjj_em_0b_lMET_"+sn ]->Fill(Mjj,weight);
	    histos["MjjL_em_0b_lMET_"+sn ]->Fill(MjjL,weight);
	    histos["pTll_em_0b_lMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	    histos["dPhillMET_em_0b_lMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	    histos["dPhill_em_0b_lMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["dEtall_em_0b_lMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["minDRlj_em_0b_lMET_"+sn ]->Fill(dRminSR,weight);
	    histos["MuPt_em_0b_lMET_"+sn ]->Fill(lep_p4()[l].Pt(),weight);
	    histos["MuPhi_em_0b_lMET_"+sn ]->Fill(lep_p4()[l].Phi(),weight);
	    histos["MuEta_em_0b_lMET_"+sn ]->Fill(lep_p4()[l].Eta(),weight);
	    histos["MuRelIso_em_0b_lMET_"+sn ]->Fill(lep_relIso03EA()[l],weight);
	    histos["MuIP3D_em_0b_lMET_"+sn ]->Fill(lep_ip3d()[l],weight);
	    histos["MuPtErrOvPt_em_0b_lMET_"+sn ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weight);
	    histos["MuValidFrac_em_0b_lMET_"+sn ]->Fill(lep_validfraction()[l],weight);
	    histos["MuChi2N_em_0b_lMET_"+sn ]->Fill(lep_x2ondof()[l],weight);
	    histos["MuLostHits_em_0b_lMET_"+sn ]->Fill(lep_lostHits()[l],weight);
	    //cout << __LINE__ << endl;
	  }
	  else {
	    //cout << __LINE__ << endl;
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_em_0b_hMET_"+sn ]->Fill(-1.,weight);
	    else                       histos["q_em_0b_hMET_"+sn ]->Fill( 1.,weight);
	    histos["Mll_em_0b_hMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	    histos["Mjj_em_0b_hMET_"+sn ]->Fill(Mjj,weight);
	    histos["MjjL_em_0b_hMET_"+sn ]->Fill(MjjL,weight);
	    histos["pTll_em_0b_hMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	    histos["dPhillMET_em_0b_hMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	    histos["dPhill_em_0b_hMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["dEtall_em_0b_hMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["minDRlj_em_0b_hMET_"+sn ]->Fill(dRminSR,weight);
	    histos["MuPt_em_0b_hMET_"+sn ]->Fill(lep_p4()[l].Pt(),weight);
	    histos["MuPhi_em_0b_hMET_"+sn ]->Fill(lep_p4()[l].Phi(),weight);
	    histos["MuEta_em_0b_hMET_"+sn ]->Fill(lep_p4()[l].Eta(),weight);
	    histos["MuRelIso_em_0b_hMET_"+sn ]->Fill(lep_relIso03EA()[l],weight);
	    histos["MuIP3D_em_0b_hMET_"+sn ]->Fill(lep_ip3d()[l],weight);
	    histos["MuPtErrOvPt_em_0b_hMET_"+sn ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weight);
	    histos["MuValidFrac_em_0b_hMET_"+sn ]->Fill(lep_validfraction()[l],weight);
	    histos["MuChi2N_em_0b_hMET_"+sn ]->Fill(lep_x2ondof()[l],weight);
	    histos["MuLostHits_em_0b_hMET_"+sn ]->Fill(lep_lostHits()[l],weight);
	    //cout << __LINE__ << endl;
	  }
	}//em[1]
	if(mm[1]){
	  //cout << __LINE__ << endl;
	  histos["MET_mm_0b_"+sn ]->Fill(met_pt(),weight);
	  histos["MTmax_mm_0b_"+sn ]->Fill(MTmax,weight);
	  if(lep_pdgId()[iSS[0] ]>0) histos["q_mm_0b_iMET_"+sn ]->Fill(-1.,weight);
	  else                       histos["q_mm_0b_iMET_"+sn ]->Fill( 1.,weight);
	  histos["Mll_mm_0b_iMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	  histos["Mjj_mm_0b_iMET_"+sn ]->Fill(Mjj,weight);
	  histos["MjjL_mm_0b_iMET_"+sn ]->Fill(MjjL,weight);
	  histos["pTll_mm_0b_iMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	  histos["dPhillMET_mm_0b_iMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	  histos["dPhill_mm_0b_iMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	  histos["dEtall_mm_0b_iMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	  histos["minDRlj_mm_0b_iMET_"+sn ]->Fill(dRminSR,weight);
	  int l = iSS[0];
	  int m = iSS[1];
	  histos["MuPt_mm_0b_iMET_"+sn ]->Fill(lep_p4()[l].Pt(),weight);
	  histos["MuPhi_mm_0b_iMET_"+sn ]->Fill(lep_p4()[l].Phi(),weight);
	  histos["MuEta_mm_0b_iMET_"+sn ]->Fill(lep_p4()[l].Eta(),weight);
	  histos["MuRelIso_mm_0b_iMET_"+sn ]->Fill(lep_relIso03EA()[l],weight);
	  histos["MuIP3D_mm_0b_iMET_"+sn ]->Fill(lep_ip3d()[l],weight);
	  histos["MuPtErrOvPt_mm_0b_iMET_"+sn ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weight);
	  histos["MuValidFrac_mm_0b_iMET_"+sn ]->Fill(lep_validfraction()[l],weight);
	  histos["MuChi2N_mm_0b_iMET_"+sn ]->Fill(lep_x2ondof()[l],weight);
	  histos["MuLostHits_mm_0b_iMET_"+sn ]->Fill(lep_lostHits()[l],weight);
	  histos["MuPt_mm_0b_iMET_"+sn ]->Fill(lep_p4()[m].Pt(),weight);
	  histos["MuPhi_mm_0b_iMET_"+sn ]->Fill(lep_p4()[m].Phi(),weight);
	  histos["MuEta_mm_0b_iMET_"+sn ]->Fill(lep_p4()[m].Eta(),weight);
	  histos["MuRelIso_mm_0b_iMET_"+sn ]->Fill(lep_relIso03EA()[m],weight);
	  histos["MuIP3D_mm_0b_iMET_"+sn ]->Fill(lep_ip3d()[m],weight);
	  histos["MuPtErrOvPt_mm_0b_iMET_"+sn ]->Fill(lep_pterr()[m]/lep_p4()[m].Pt(),weight);
	  histos["MuValidFrac_mm_0b_iMET_"+sn ]->Fill(lep_validfraction()[m],weight);
	  histos["MuChi2N_mm_0b_iMET_"+sn ]->Fill(lep_x2ondof()[m],weight);
	  histos["MuLostHits_mm_0b_iMET_"+sn ]->Fill(lep_lostHits()[m],weight);
	  if(met_pt()<40){
	    //cout << __LINE__ << endl;
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_mm_0b_lMET_"+sn ]->Fill(-1.,weight);
	    else                       histos["q_mm_0b_lMET_"+sn ]->Fill( 1.,weight);
	    histos["Mll_mm_0b_lMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	    histos["Mjj_mm_0b_lMET_"+sn ]->Fill(Mjj,weight);
	    histos["MjjL_mm_0b_lMET_"+sn ]->Fill(MjjL,weight);
	    histos["pTll_mm_0b_lMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	    histos["dPhillMET_mm_0b_lMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	    histos["dPhill_mm_0b_lMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["dEtall_mm_0b_lMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["minDRlj_mm_0b_lMET_"+sn ]->Fill(dRminSR,weight);
	    histos["MuPt_mm_0b_lMET_"+sn ]->Fill(lep_p4()[l].Pt(),weight);
	    histos["MuPhi_mm_0b_lMET_"+sn ]->Fill(lep_p4()[l].Phi(),weight);
	    histos["MuEta_mm_0b_lMET_"+sn ]->Fill(lep_p4()[l].Eta(),weight);
	    histos["MuRelIso_mm_0b_lMET_"+sn ]->Fill(lep_relIso03EA()[l],weight);
	    histos["MuIP3D_mm_0b_lMET_"+sn ]->Fill(lep_ip3d()[l],weight);
	    histos["MuPtErrOvPt_mm_0b_lMET_"+sn ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weight);
	    histos["MuValidFrac_mm_0b_lMET_"+sn ]->Fill(lep_validfraction()[l],weight);
	    histos["MuChi2N_mm_0b_lMET_"+sn ]->Fill(lep_x2ondof()[l],weight);
	    histos["MuLostHits_mm_0b_lMET_"+sn ]->Fill(lep_lostHits()[l],weight);
	    histos["MuPt_mm_0b_lMET_"+sn ]->Fill(lep_p4()[m].Pt(),weight);
	    histos["MuPhi_mm_0b_lMET_"+sn ]->Fill(lep_p4()[m].Phi(),weight);
	    histos["MuEta_mm_0b_lMET_"+sn ]->Fill(lep_p4()[m].Eta(),weight);
	    histos["MuRelIso_mm_0b_lMET_"+sn ]->Fill(lep_relIso03EA()[m],weight);
	    histos["MuIP3D_mm_0b_lMET_"+sn ]->Fill(lep_ip3d()[m],weight);
	    histos["MuPtErrOvPt_mm_0b_lMET_"+sn ]->Fill(lep_pterr()[m]/lep_p4()[m].Pt(),weight);
	    histos["MuValidFrac_mm_0b_lMET_"+sn ]->Fill(lep_validfraction()[m],weight);
	    histos["MuChi2N_mm_0b_lMET_"+sn ]->Fill(lep_x2ondof()[m],weight);
	    histos["MuLostHits_mm_0b_lMET_"+sn ]->Fill(lep_lostHits()[m],weight);
	  }
	  else {
	    //cout << __LINE__ << endl;
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_mm_0b_hMET_"+sn ]->Fill(-1.,weight);
	    else                       histos["q_mm_0b_hMET_"+sn ]->Fill( 1.,weight);
	    histos["Mll_mm_0b_hMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).M(),weight);
	    histos["Mjj_mm_0b_hMET_"+sn ]->Fill(Mjj,weight);
	    histos["MjjL_mm_0b_hMET_"+sn ]->Fill(MjjL,weight);
	    histos["pTll_mm_0b_hMET_"+sn ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ]).Pt(),weight);
	    histos["dPhillMET_mm_0b_hMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iSS[1] ],MET),weight);
	    histos["dPhill_mm_0b_hMET_"+sn ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["dEtall_mm_0b_hMET_"+sn ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iSS[1] ]),weight);
	    histos["minDRlj_mm_0b_hMET_"+sn ]->Fill(dRminSR,weight);
	    histos["MuPt_mm_0b_hMET_"+sn ]->Fill(lep_p4()[l].Pt(),weight);
	    histos["MuPhi_mm_0b_hMET_"+sn ]->Fill(lep_p4()[l].Phi(),weight);
	    histos["MuEta_mm_0b_hMET_"+sn ]->Fill(lep_p4()[l].Eta(),weight);
	    histos["MuRelIso_mm_0b_hMET_"+sn ]->Fill(lep_relIso03EA()[l],weight);
	    histos["MuIP3D_mm_0b_hMET_"+sn ]->Fill(lep_ip3d()[l],weight);
	    histos["MuPtErrOvPt_mm_0b_hMET_"+sn ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weight);
	    histos["MuValidFrac_mm_0b_hMET_"+sn ]->Fill(lep_validfraction()[l],weight);
	    histos["MuChi2N_mm_0b_hMET_"+sn ]->Fill(lep_x2ondof()[l],weight);
	    histos["MuLostHits_mm_0b_hMET_"+sn ]->Fill(lep_lostHits()[l],weight);
	    histos["MuPt_mm_0b_hMET_"+sn ]->Fill(lep_p4()[m].Pt(),weight);
	    histos["MuPhi_mm_0b_hMET_"+sn ]->Fill(lep_p4()[m].Phi(),weight);
	    histos["MuEta_mm_0b_hMET_"+sn ]->Fill(lep_p4()[m].Eta(),weight);
	    histos["MuRelIso_mm_0b_hMET_"+sn ]->Fill(lep_relIso03EA()[m],weight);
	    histos["MuIP3D_mm_0b_hMET_"+sn ]->Fill(lep_ip3d()[m],weight);
	    histos["MuPtErrOvPt_mm_0b_hMET_"+sn ]->Fill(lep_pterr()[m]/lep_p4()[m].Pt(),weight);
	    histos["MuValidFrac_mm_0b_hMET_"+sn ]->Fill(lep_validfraction()[m],weight);
	    histos["MuChi2N_mm_0b_hMET_"+sn ]->Fill(lep_x2ondof()[m],weight);
	    histos["MuLostHits_mm_0b_hMET_"+sn ]->Fill(lep_lostHits()[m],weight);
	  }
	  //cout << __LINE__ << endl;
	}//mm[1]

      }//nb==0



      if(isData()){ //data
	string snFR = "FakePred";
	double weightFR = weight;
	if(naSS>=1){
	  if(abs(lep_pdgId()[iaSS[0] ])==11) {
	    int bin = hElFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)elptmaxFR),(double)elptminFR),
				     std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)eletamaxFR),(double)eletaminFR) );
	    weightFR *= hElFR->GetBinContent(bin)/(1.-hElFR->GetBinContent(bin));
	  } else {
	    int bin = hMuFR->FindBin(std::max(std::min((double)(1.+coneCorrPt(iaSS[0]))*lep_p4()[iaSS[0] ].Pt(),(double)muptmaxFR),(double)muptminFR),
				     std::max(std::min(fabs((double)lep_p4()[iaSS[0] ].Eta()),(double)muetamaxFR),(double)muetaminFR) );
	    weightFR *= hMuFR->GetBinContent(bin)/(1.-hmuFR->GetBinContent(bin));
	  }
	}
	if(ee[3]     ) histos["NB_ee_iMET_"+snFR ]->Fill(nb,weightFR);
	if(em[3]     ) histos["NB_em_iMET_"+snFR ]->Fill(nb,weightFR);
	if(mm[3]     ) histos["NB_mm_iMET_"+snFR ]->Fill(nb,weightFR);
	if(ee[3]&&met_pt()<40     ) histos["NB_ee_lMET_"+snFR ]->Fill(nb,weightFR);
	if(em[3]&&met_pt()<40     ) histos["NB_em_lMET_"+snFR ]->Fill(nb,weightFR);
	if(mm[3]&&met_pt()<40     ) histos["NB_mm_lMET_"+snFR ]->Fill(nb,weightFR);
	if(ee[3]&&met_pt()>40     ) histos["NB_ee_hMET_"+snFR ]->Fill(nb,weightFR);
	if(em[3]&&met_pt()>40     ) histos["NB_em_hMET_"+snFR ]->Fill(nb,weightFR);
	if(mm[3]&&met_pt()>40     ) histos["NB_mm_hMET_"+snFR ]->Fill(nb,weightFR);
	if(nb==0){
	  if(fname.find("data_Run2016B")!=string::npos&&nb==0){
	    if(ee[8]     ) histos["SR_ee_" +snFR ]->Fill(0.,weightFR);
	    if(em[8]     ) histos["SR_em_" +snFR ]->Fill(0.,weightFR);
	    if(mm[8]     ) histos["SR_mm_" +snFR ]->Fill(0.,weightFR);
	    if(ee[3]     ) histos["SRpresel_ee_" +snFR ]->Fill(0.,weightFR);
	    if(em[3]     ) histos["SRpresel_em_" +snFR ]->Fill(0.,weightFR);
	    if(mm[3]     ) histos["SRpresel_mm_" +snFR ]->Fill(0.,weightFR);
	  }
	  if(fname.find("data_Run2016C")!=string::npos&&nb==0){
	    if(ee[8]     ) histos["SR_ee_" +snFR ]->Fill(1.,weightFR);
	    if(em[8]     ) histos["SR_em_" +snFR ]->Fill(1.,weightFR);
	    if(mm[8]     ) histos["SR_mm_" +snFR ]->Fill(1.,weightFR);
	    if(ee[3]     ) histos["SRpresel_ee_" +snFR ]->Fill(1.,weightFR);
	    if(em[3]     ) histos["SRpresel_em_" +snFR ]->Fill(1.,weightFR);
	    if(mm[3]     ) histos["SRpresel_mm_" +snFR ]->Fill(1.,weightFR);
	  }
	  if(fname.find("data_Run2016D")!=string::npos&&nb==0){
	    if(ee[8]     ) histos["SR_ee_" +snFR ]->Fill(2.,weightFR);
	    if(em[8]     ) histos["SR_em_" +snFR ]->Fill(2.,weightFR);
	    if(mm[8]     ) histos["SR_mm_" +snFR ]->Fill(2.,weightFR);
	    if(ee[3]     ) histos["SRpresel_ee_" +snFR ]->Fill(2.,weightFR);
	    if(em[3]     ) histos["SRpresel_em_" +snFR ]->Fill(2.,weightFR);
	    if(mm[3]     ) histos["SRpresel_mm_" +snFR ]->Fill(2.,weightFR);
	  }
	  if(fname.find("data_Run2016E")!=string::npos&&nb==0){
	    if(ee[8]     ) histos["SR_ee_" +snFR ]->Fill(3.,weightFR);
	    if(em[8]     ) histos["SR_em_" +snFR ]->Fill(3.,weightFR);
	    if(mm[8]     ) histos["SR_mm_" +snFR ]->Fill(3.,weightFR);
	    if(ee[3]     ) histos["SRpresel_ee_" +snFR ]->Fill(3.,weightFR);
	    if(em[3]     ) histos["SRpresel_em_" +snFR ]->Fill(3.,weightFR);
	    if(mm[3]     ) histos["SRpresel_mm_" +snFR ]->Fill(3.,weightFR);
	  }
	  if(fname.find("data_Run2016F")!=string::npos&&nb==0){
	    if(ee[8]     ) histos["SR_ee_" +snFR ]->Fill(4.,weightFR);
	    if(em[8]     ) histos["SR_em_" +snFR ]->Fill(4.,weightFR);
	    if(mm[8]     ) histos["SR_mm_" +snFR ]->Fill(4.,weightFR);
	    if(ee[3]     ) histos["SRpresel_ee_" +snFR ]->Fill(4.,weightFR);
	    if(em[3]     ) histos["SRpresel_em_" +snFR ]->Fill(4.,weightFR);
	    if(mm[3]     ) histos["SRpresel_mm_" +snFR ]->Fill(4.,weightFR);
	  }
	  if(fname.find("data_Run2016G")!=string::npos&&nb==0){
	    if(ee[8]     ) histos["SR_ee_" +snFR ]->Fill(5.,weightFR);
	    if(em[8]     ) histos["SR_em_" +snFR ]->Fill(5.,weightFR);
	    if(mm[8]     ) histos["SR_mm_" +snFR ]->Fill(5.,weightFR);
	    if(ee[3]     ) histos["SRpresel_ee_" +snFR ]->Fill(5.,weightFR);
	    if(em[3]     ) histos["SRpresel_em_" +snFR ]->Fill(5.,weightFR);
	    if(mm[3]     ) histos["SRpresel_mm_" +snFR ]->Fill(5.,weightFR);
	  }
	  if(fname.find("data_Run2016H")!=string::npos&&nb==0){
	    if(ee[8]     ) histos["SR_ee_" +snFR ]->Fill(6.,weightFR);
	    if(em[8]     ) histos["SR_em_" +snFR ]->Fill(6.,weightFR);
	    if(mm[8]     ) histos["SR_mm_" +snFR ]->Fill(6.,weightFR);
	    if(ee[3]     ) histos["SRpresel_ee_" +snFR ]->Fill(6.,weightFR);
	    if(em[3]     ) histos["SRpresel_em_" +snFR ]->Fill(6.,weightFR);
	    if(mm[3]     ) histos["SRpresel_mm_" +snFR ]->Fill(6.,weightFR);
	  }
	  if(ee[8]&&nb==0     ) histos["SR_iMET_"+snFR ]->Fill(0.,weightFR);
	  if(em[8]&&nb==0     ) histos["SR_iMET_"+snFR ]->Fill(1.,weightFR);
	  if(mm[8]&&nb==0     ) histos["SR_iMET_"+snFR ]->Fill(2.,weightFR);
	  if(ee[8]&&nb==0&&met_pt()<40     ) histos["SR_lMET_"+snFR ]->Fill(0.,weightFR);
	  if(em[8]&&nb==0&&met_pt()<40     ) histos["SR_lMET_"+snFR ]->Fill(1.,weightFR);
	  if(mm[8]&&nb==0&&met_pt()<40     ) histos["SR_lMET_"+snFR ]->Fill(2.,weightFR);
	  if(ee[8]&&nb==0&&met_pt()>40     ) histos["SR_hMET_"+snFR ]->Fill(0.,weightFR);
	  if(em[8]&&nb==0&&met_pt()>40     ) histos["SR_hMET_"+snFR ]->Fill(1.,weightFR);
	  if(mm[8]&&nb==0&&met_pt()>40     ) histos["SR_hMET_"+snFR ]->Fill(2.,weightFR);
	  if(ee[3]&&nb==0     ) histos["SRpresel_iMET_"+snFR ]->Fill(0.,weightFR);
	  if(em[3]&&nb==0     ) histos["SRpresel_iMET_"+snFR ]->Fill(1.,weightFR);
	  if(mm[3]&&nb==0     ) histos["SRpresel_iMET_"+snFR ]->Fill(2.,weightFR);
	  if(ee[3]&&nb==0&&met_pt()<40     ) histos["SRpresel_lMET_"+snFR ]->Fill(0.,weightFR);
	  if(em[3]&&nb==0&&met_pt()<40     ) histos["SRpresel_lMET_"+snFR ]->Fill(1.,weightFR);
	  if(mm[3]&&nb==0&&met_pt()<40     ) histos["SRpresel_lMET_"+snFR ]->Fill(2.,weightFR);
	  if(ee[3]&&nb==0&&met_pt()>40     ) histos["SRpresel_hMET_"+snFR ]->Fill(0.,weightFR);
	  if(em[3]&&nb==0&&met_pt()>40     ) histos["SRpresel_hMET_"+snFR ]->Fill(1.,weightFR);
	  if(mm[3]&&nb==0&&met_pt()>40     ) histos["SRpresel_hMET_"+snFR ]->Fill(2.,weightFR);

	  if(ee[9]&&nb==0     ) histos["NJ_ee_0b_iMET_"+snFR ]->Fill(nj30,weightFR);
	  if(em[9]&&nb==0     ) histos["NJ_em_0b_iMET_"+snFR ]->Fill(nj30,weightFR);
	  if(mm[9]&&nb==0     ) histos["NJ_mm_0b_iMET_"+snFR ]->Fill(nj30,weightFR);
	  if(ee[9]&&nb==0&&met_pt()<40     ) histos["NJ_ee_0b_lMET_"+snFR ]->Fill(nj30,weightFR);
	  if(em[9]&&nb==0&&met_pt()<40     ) histos["NJ_em_0b_lMET_"+snFR ]->Fill(nj30,weightFR);
	  if(mm[9]&&nb==0&&met_pt()<40     ) histos["NJ_mm_0b_lMET_"+snFR ]->Fill(nj30,weightFR);
	  if(ee[9]&&nb==0&&met_pt()>40     ) histos["NJ_ee_0b_hMET_"+snFR ]->Fill(nj30,weightFR);
	  if(em[9]&&nb==0&&met_pt()>40     ) histos["NJ_em_0b_hMET_"+snFR ]->Fill(nj30,weightFR);
	  if(mm[9]&&nb==0&&met_pt()>40     ) histos["NJ_mm_0b_hMET_"+snFR ]->Fill(nj30,weightFR);

	  if(ee[3]){
	    histos["MET_ee_0b_"+snFR ]->Fill(met_pt(),weightFR);
	    histos["MTmax_ee_0b_"+snFR ]->Fill(MTmax,weightFR);
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_ee_0b_iMET_"+snFR ]->Fill(-1.,weightFR);
	    else                       histos["q_ee_0b_iMET_"+snFR ]->Fill( 1.,weightFR);
	    histos["Mll_ee_0b_iMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	    histos["Mjj_ee_0b_iMET_"+snFR ]->Fill(Mjj,weightFR);
	    histos["MjjL_ee_0b_iMET_"+snFR ]->Fill(MjjL,weightFR);
	    histos["pTll_ee_0b_iMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	    histos["dPhillMET_ee_0b_iMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	    histos["dPhill_ee_0b_iMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	    histos["dEtall_ee_0b_iMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	    histos["minDRlj_ee_0b_iMET_"+snFR ]->Fill(dRminSR,weightFR);
	    if(met_pt()<40){
	      if(lep_pdgId()[iSS[0] ]>0) histos["q_ee_0b_lMET_"+snFR ]->Fill(-1.,weightFR);
	      else                       histos["q_ee_0b_lMET_"+snFR ]->Fill( 1.,weightFR);
	      histos["Mll_ee_0b_lMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	      histos["Mjj_ee_0b_lMET_"+snFR ]->Fill(Mjj,weightFR);
	      histos["MjjL_ee_0b_lMET_"+snFR ]->Fill(MjjL,weightFR);
	      histos["pTll_ee_0b_lMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	      histos["dPhillMET_ee_0b_lMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	      histos["dPhill_ee_0b_lMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["dEtall_ee_0b_lMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["minDRlj_ee_0b_lMET_"+snFR ]->Fill(dRminSR,weightFR);
	    }
	    else {
	      if(lep_pdgId()[iSS[0] ]>0) histos["q_ee_0b_hMET_"+snFR ]->Fill(-1.,weightFR);
	      else                       histos["q_ee_0b_hMET_"+snFR ]->Fill( 1.,weightFR);
	      histos["Mll_ee_0b_hMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	      histos["Mjj_ee_0b_hMET_"+snFR ]->Fill(Mjj,weightFR);
	      histos["MjjL_ee_0b_hMET_"+snFR ]->Fill(MjjL,weightFR);
	      histos["pTll_ee_0b_hMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	      histos["dPhillMET_ee_0b_hMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	      histos["dPhill_ee_0b_hMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["dEtall_ee_0b_hMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["minDRlj_ee_0b_hMET_"+snFR ]->Fill(dRminSR,weightFR);
	    }
	  }//ee[3]
	  if(em[3]){
	    histos["MET_em_0b_"+snFR ]->Fill(met_pt(),weightFR);
	    histos["MTmax_em_0b_"+snFR ]->Fill(MTmax,weightFR);
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_em_0b_iMET_"+snFR ]->Fill(-1.,weightFR);
	    else                       histos["q_em_0b_iMET_"+snFR ]->Fill( 1.,weightFR);
	    histos["Mll_em_0b_iMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	    histos["Mjj_em_0b_iMET_"+snFR ]->Fill(Mjj,weightFR);
	    histos["MjjL_em_0b_iMET_"+snFR ]->Fill(MjjL,weightFR);
	    histos["pTll_em_0b_iMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	    histos["dPhillMET_em_0b_iMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	    histos["dPhill_em_0b_iMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	    histos["dEtall_em_0b_iMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	    histos["minDRlj_em_0b_iMET_"+snFR ]->Fill(dRminSR,weightFR);
	    int l = -1;
	    if(abs(lep_pdgId()[iSS[0] ])==13) l = iSS[0];
	    else if(abs(lep_pdgId()[iaSS[0] ])==13) l = iaSS[0];
	    else cout << "WTF em" << endl;
	    histos["MuPt_em_0b_iMET_"+snFR ]->Fill(lep_p4()[l].Pt(),weightFR);
	    histos["MuPhi_em_0b_iMET_"+snFR ]->Fill(lep_p4()[l].Phi(),weightFR);
	    histos["MuEta_em_0b_iMET_"+snFR ]->Fill(lep_p4()[l].Eta(),weightFR);
	    histos["MuRelIso_em_0b_iMET_"+snFR ]->Fill(lep_relIso03EA()[l],weightFR);
	    histos["MuIP3D_em_0b_iMET_"+snFR ]->Fill(lep_ip3d()[l],weightFR);
	    histos["MuPtErrOvPt_em_0b_iMET_"+snFR ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weightFR);
	    histos["MuValidFrac_em_0b_iMET_"+snFR ]->Fill(lep_validfraction()[l],weightFR);
	    histos["MuChi2N_em_0b_iMET_"+snFR ]->Fill(lep_x2ondof()[l],weightFR);
	    histos["MuLostHits_em_0b_iMET_"+snFR ]->Fill(lep_lostHits()[l],weightFR);
	    if(met_pt()<40){
	      if(lep_pdgId()[iSS[0] ]>0) histos["q_em_0b_lMET_"+snFR ]->Fill(-1.,weightFR);
	      else                       histos["q_em_0b_lMET_"+snFR ]->Fill( 1.,weightFR);
	      histos["Mll_em_0b_lMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	      histos["Mjj_em_0b_lMET_"+snFR ]->Fill(Mjj,weightFR);
	      histos["MjjL_em_0b_lMET_"+snFR ]->Fill(MjjL,weightFR);
	      histos["pTll_em_0b_lMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	      histos["dPhillMET_em_0b_lMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	      histos["dPhill_em_0b_lMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["dEtall_em_0b_lMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["minDRlj_em_0b_lMET_"+snFR ]->Fill(dRminSR,weightFR);
	      histos["MuPt_em_0b_lMET_"+snFR ]->Fill(lep_p4()[l].Pt(),weightFR);
	      histos["MuPhi_em_0b_lMET_"+snFR ]->Fill(lep_p4()[l].Phi(),weightFR);
	      histos["MuEta_em_0b_lMET_"+snFR ]->Fill(lep_p4()[l].Eta(),weightFR);
	      histos["MuRelIso_em_0b_lMET_"+snFR ]->Fill(lep_relIso03EA()[l],weightFR);
	      histos["MuIP3D_em_0b_lMET_"+snFR ]->Fill(lep_ip3d()[l],weightFR);
	      histos["MuPtErrOvPt_em_0b_lMET_"+snFR ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weightFR);
	      histos["MuValidFrac_em_0b_lMET_"+snFR ]->Fill(lep_validfraction()[l],weightFR);
	      histos["MuChi2N_em_0b_lMET_"+snFR ]->Fill(lep_x2ondof()[l],weightFR);
	      histos["MuLostHits_em_0b_lMET_"+snFR ]->Fill(lep_lostHits()[l],weightFR);
	    }
	    else {
	      if(lep_pdgId()[iSS[0] ]>0) histos["q_em_0b_hMET_"+snFR ]->Fill(-1.,weightFR);
	      else                       histos["q_em_0b_hMET_"+snFR ]->Fill( 1.,weightFR);
	      histos["Mll_em_0b_hMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	      histos["Mjj_em_0b_hMET_"+snFR ]->Fill(Mjj,weightFR);
	      histos["MjjL_em_0b_hMET_"+snFR ]->Fill(MjjL,weightFR);
	      histos["pTll_em_0b_hMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	      histos["dPhillMET_em_0b_hMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	      histos["dPhill_em_0b_hMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["dEtall_em_0b_hMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["minDRlj_em_0b_hMET_"+snFR ]->Fill(dRminSR,weightFR);
	      histos["MuPt_em_0b_hMET_"+snFR ]->Fill(lep_p4()[l].Pt(),weightFR);
	      histos["MuPhi_em_0b_hMET_"+snFR ]->Fill(lep_p4()[l].Phi(),weightFR);
	      histos["MuEta_em_0b_hMET_"+snFR ]->Fill(lep_p4()[l].Eta(),weightFR);
	      histos["MuRelIso_em_0b_hMET_"+snFR ]->Fill(lep_relIso03EA()[l],weightFR);
	      histos["MuIP3D_em_0b_hMET_"+snFR ]->Fill(lep_ip3d()[l],weightFR);
	      histos["MuPtErrOvPt_em_0b_hMET_"+snFR ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weightFR);
	      histos["MuValidFrac_em_0b_hMET_"+snFR ]->Fill(lep_validfraction()[l],weightFR);
	      histos["MuChi2N_em_0b_hMET_"+snFR ]->Fill(lep_x2ondof()[l],weightFR);
	      histos["MuLostHits_em_0b_hMET_"+snFR ]->Fill(lep_lostHits()[l],weightFR);
	    }
	  }//em[3]
	  if(mm[3]){
	    histos["MET_mm_0b_"+snFR ]->Fill(met_pt(),weightFR);
	    histos["MTmax_mm_0b_"+snFR ]->Fill(MTmax,weightFR);
	    if(lep_pdgId()[iSS[0] ]>0) histos["q_mm_0b_iMET_"+snFR ]->Fill(-1.,weightFR);
	    else                       histos["q_mm_0b_iMET_"+snFR ]->Fill( 1.,weightFR);
	    histos["Mll_mm_0b_iMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	    histos["Mjj_mm_0b_iMET_"+snFR ]->Fill(Mjj,weightFR);
	    histos["MjjL_mm_0b_iMET_"+snFR ]->Fill(MjjL,weightFR);
	    histos["pTll_mm_0b_iMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	    histos["dPhillMET_mm_0b_iMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	    histos["dPhill_mm_0b_iMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	    histos["dEtall_mm_0b_iMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	    histos["minDRlj_mm_0b_iMET_"+snFR ]->Fill(dRminSR,weightFR);
	    int l = iSS[0];
	    int m = iSS[1];
	    histos["MuPt_mm_0b_iMET_"+snFR ]->Fill(lep_p4()[l].Pt(),weightFR);
	    histos["MuPhi_mm_0b_iMET_"+snFR ]->Fill(lep_p4()[l].Phi(),weightFR);
	    histos["MuEta_mm_0b_iMET_"+snFR ]->Fill(lep_p4()[l].Eta(),weightFR);
	    histos["MuRelIso_mm_0b_iMET_"+snFR ]->Fill(lep_relIso03EA()[l],weightFR);
	    histos["MuIP3D_mm_0b_iMET_"+snFR ]->Fill(lep_ip3d()[l],weightFR);
	    histos["MuPtErrOvPt_mm_0b_iMET_"+snFR ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weightFR);
	    histos["MuValidFrac_mm_0b_iMET_"+snFR ]->Fill(lep_validfraction()[l],weightFR);
	    histos["MuChi2N_mm_0b_iMET_"+snFR ]->Fill(lep_x2ondof()[l],weightFR);
	    histos["MuLostHits_mm_0b_iMET_"+snFR ]->Fill(lep_lostHits()[l],weightFR);
	    histos["MuPt_mm_0b_iMET_"+snFR ]->Fill(lep_p4()[m].Pt(),weightFR);
	    histos["MuPhi_mm_0b_iMET_"+snFR ]->Fill(lep_p4()[m].Phi(),weightFR);
	    histos["MuEta_mm_0b_iMET_"+snFR ]->Fill(lep_p4()[m].Eta(),weightFR);
	    histos["MuRelIso_mm_0b_iMET_"+snFR ]->Fill(lep_relIso03EA()[m],weightFR);
	    histos["MuIP3D_mm_0b_iMET_"+snFR ]->Fill(lep_ip3d()[m],weightFR);
	    histos["MuPtErrOvPt_mm_0b_iMET_"+snFR ]->Fill(lep_pterr()[m]/lep_p4()[m].Pt(),weightFR);
	    histos["MuValidFrac_mm_0b_iMET_"+snFR ]->Fill(lep_validfraction()[m],weightFR);
	    histos["MuChi2N_mm_0b_iMET_"+snFR ]->Fill(lep_x2ondof()[m],weightFR);
	    histos["MuLostHits_mm_0b_iMET_"+snFR ]->Fill(lep_lostHits()[m],weightFR);
	    if(met_pt()<40){
	      if(lep_pdgId()[iSS[0] ]>0) histos["q_mm_0b_lMET_"+snFR ]->Fill(-1.,weightFR);
	      else                       histos["q_mm_0b_lMET_"+snFR ]->Fill( 1.,weightFR);
	      histos["Mll_mm_0b_lMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	      histos["Mjj_mm_0b_lMET_"+snFR ]->Fill(Mjj,weightFR);
	      histos["MjjL_mm_0b_lMET_"+snFR ]->Fill(MjjL,weightFR);
	      histos["pTll_mm_0b_lMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	      histos["dPhillMET_mm_0b_lMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	      histos["dPhill_mm_0b_lMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["dEtall_mm_0b_lMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["minDRlj_mm_0b_lMET_"+snFR ]->Fill(dRminSR,weightFR);
	      histos["MuPt_mm_0b_lMET_"+snFR ]->Fill(lep_p4()[l].Pt(),weightFR);
	      histos["MuPhi_mm_0b_lMET_"+snFR ]->Fill(lep_p4()[l].Phi(),weightFR);
	      histos["MuEta_mm_0b_lMET_"+snFR ]->Fill(lep_p4()[l].Eta(),weightFR);
	      histos["MuRelIso_mm_0b_lMET_"+snFR ]->Fill(lep_relIso03EA()[l],weightFR);
	      histos["MuIP3D_mm_0b_lMET_"+snFR ]->Fill(lep_ip3d()[l],weightFR);
	      histos["MuPtErrOvPt_mm_0b_lMET_"+snFR ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weightFR);
	      histos["MuValidFrac_mm_0b_lMET_"+snFR ]->Fill(lep_validfraction()[l],weightFR);
	      histos["MuChi2N_mm_0b_lMET_"+snFR ]->Fill(lep_x2ondof()[l],weightFR);
	      histos["MuLostHits_mm_0b_lMET_"+snFR ]->Fill(lep_lostHits()[l],weightFR);
	      histos["MuPt_mm_0b_lMET_"+snFR ]->Fill(lep_p4()[m].Pt(),weightFR);
	      histos["MuPhi_mm_0b_lMET_"+snFR ]->Fill(lep_p4()[m].Phi(),weightFR);
	      histos["MuEta_mm_0b_lMET_"+snFR ]->Fill(lep_p4()[m].Eta(),weightFR);
	      histos["MuRelIso_mm_0b_lMET_"+snFR ]->Fill(lep_relIso03EA()[m],weightFR);
	      histos["MuIP3D_mm_0b_lMET_"+snFR ]->Fill(lep_ip3d()[m],weightFR);
	      histos["MuPtErrOvPt_mm_0b_lMET_"+snFR ]->Fill(lep_pterr()[m]/lep_p4()[m].Pt(),weightFR);
	      histos["MuValidFrac_mm_0b_lMET_"+snFR ]->Fill(lep_validfraction()[m],weightFR);
	      histos["MuChi2N_mm_0b_lMET_"+snFR ]->Fill(lep_x2ondof()[m],weightFR);
	      histos["MuLostHits_mm_0b_lMET_"+snFR ]->Fill(lep_lostHits()[m],weightFR);
	    }
	    else {
	      if(lep_pdgId()[iSS[0] ]>0) histos["q_mm_0b_hMET_"+snFR ]->Fill(-1.,weightFR);
	      else                       histos["q_mm_0b_hMET_"+snFR ]->Fill( 1.,weightFR);
	      histos["Mll_mm_0b_hMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).M(),weightFR);
	      histos["Mjj_mm_0b_hMET_"+snFR ]->Fill(Mjj,weightFR);
	      histos["MjjL_mm_0b_hMET_"+snFR ]->Fill(MjjL,weightFR);
	      histos["pTll_mm_0b_hMET_"+snFR ]->Fill((lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ]).Pt(),weightFR);
	      histos["dPhillMET_mm_0b_hMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ]+lep_p4()[ iaSS[0] ],MET),weightFR);
	      histos["dPhill_mm_0b_hMET_"+snFR ]->Fill(dPhi(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["dEtall_mm_0b_hMET_"+snFR ]->Fill(dEta(lep_p4()[ iSS[0] ],lep_p4()[ iaSS[0] ]),weightFR);
	      histos["minDRlj_mm_0b_hMET_"+snFR ]->Fill(dRminSR,weightFR);
	      histos["MuPt_mm_0b_hMET_"+snFR ]->Fill(lep_p4()[l].Pt(),weightFR);
	      histos["MuPhi_mm_0b_hMET_"+snFR ]->Fill(lep_p4()[l].Phi(),weightFR);
	      histos["MuEta_mm_0b_hMET_"+snFR ]->Fill(lep_p4()[l].Eta(),weightFR);
	      histos["MuRelIso_mm_0b_hMET_"+snFR ]->Fill(lep_relIso03EA()[l],weightFR);
	      histos["MuIP3D_mm_0b_hMET_"+snFR ]->Fill(lep_ip3d()[l],weightFR);
	      histos["MuPtErrOvPt_mm_0b_hMET_"+snFR ]->Fill(lep_pterr()[l]/lep_p4()[l].Pt(),weightFR);
	      histos["MuValidFrac_mm_0b_hMET_"+snFR ]->Fill(lep_validfraction()[l],weightFR);
	      histos["MuChi2N_mm_0b_hMET_"+snFR ]->Fill(lep_x2ondof()[l],weightFR);
	      histos["MuLostHits_mm_0b_hMET_"+snFR ]->Fill(lep_lostHits()[l],weightFR);
	      histos["MuPt_mm_0b_hMET_"+snFR ]->Fill(lep_p4()[m].Pt(),weightFR);
	      histos["MuPhi_mm_0b_hMET_"+snFR ]->Fill(lep_p4()[m].Phi(),weightFR);
	      histos["MuEta_mm_0b_hMET_"+snFR ]->Fill(lep_p4()[m].Eta(),weightFR);
	      histos["MuRelIso_mm_0b_hMET_"+snFR ]->Fill(lep_relIso03EA()[m],weightFR);
	      histos["MuIP3D_mm_0b_hMET_"+snFR ]->Fill(lep_ip3d()[m],weightFR);
	      histos["MuPtErrOvPt_mm_0b_hMET_"+snFR ]->Fill(lep_pterr()[m]/lep_p4()[m].Pt(),weightFR);
	      histos["MuValidFrac_mm_0b_hMET_"+snFR ]->Fill(lep_validfraction()[m],weightFR);
	      histos["MuChi2N_mm_0b_hMET_"+snFR ]->Fill(lep_x2ondof()[m],weightFR);
	      histos["MuLostHits_mm_0b_hMET_"+snFR ]->Fill(lep_lostHits()[m],weightFR);
	    }
	  }//mm[3]
	}//nb==0
      }//isData for FR

    
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
  
  string filename = "rootfiles/ExcessPlots.root";
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
  fSF->Close();
  delete fSF;
  fFR->Close();
  delete fFR;
  
  return 0;
}
