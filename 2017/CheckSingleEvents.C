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

struct myevt{
  unsigned int run;
  unsigned int ls;
  long long evt;
};

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

  vector<myevt> e;
  myevt ee;
  /*
  ee.run = 274335; ee.ls =  929; ee.evt = 1713663464; e.push_back(ee);
  ee.run = 275068; ee.ls =  688; ee.evt = 1275633476; e.push_back(ee);
  ee.run = 275001; ee.ls = 1492; ee.evt = 2023852463; e.push_back(ee);
  ee.run = 275293; ee.ls =   35; ee.evt = 63013573;   e.push_back(ee);
  ee.run = 275890; ee.ls =  433; ee.evt = 819661965;  e.push_back(ee);
  ee.run = 277070; ee.ls =  920; ee.evt = 1719156153; e.push_back(ee);
  ee.run = 278018; ee.ls =   63; ee.evt = 114757971;  e.push_back(ee);
  ee.run = 280251; ee.ls =  370; ee.evt = 656924162;  e.push_back(ee);
  ee.run = 280018; ee.ls =  723; ee.evt = 1255769696; e.push_back(ee);
  ee.run = 279667; ee.ls =  398; ee.evt = 610845548;  e.push_back(ee);
  ee.run = 278873; ee.ls =  124; ee.evt = 46701135;   e.push_back(ee);
  ee.run = 279658; ee.ls =  363; ee.evt = 640869423;  e.push_back(ee);
  ee.run = 283885; ee.ls =  628; ee.evt = 1150203211; e.push_back(ee);
  */
  ee.run = 275068; ee.ls =   59; ee.evt = 121092033;  e.push_back(ee);
  ee.run = 275338; ee.ls =  205; ee.evt = 345522396;  e.push_back(ee);
  ee.run = 275310; ee.ls =  659; ee.evt = 1156930149; e.push_back(ee);
  ee.run = 275001; ee.ls = 2013; ee.evt = 2630023366; e.push_back(ee);
  ee.run = 273403; ee.ls =   13; ee.evt = 15292405;   e.push_back(ee);
  ee.run = 274971; ee.ls =  210; ee.evt = 321797534;  e.push_back(ee);
  ee.run = 274969; ee.ls =  525; ee.evt = 957161415;  e.push_back(ee);
  ee.run = 275847; ee.ls = 1825; ee.evt = 2237393680; e.push_back(ee);
  ee.run = 276097; ee.ls =  426; ee.evt = 818694347;  e.push_back(ee);
  ee.run = 276384; ee.ls =  743; ee.evt = 1081674149; e.push_back(ee);
  ee.run = 276582; ee.ls =  652; ee.evt = 1172445440; e.push_back(ee);
  ee.run = 276655; ee.ls =  353; ee.evt = 661756787;  e.push_back(ee);
  ee.run = 277194; ee.ls =  939; ee.evt = 1632932890; e.push_back(ee);
  ee.run = 277202; ee.ls =   42; ee.evt = 32078429;   e.push_back(ee);
  ee.run = 276948; ee.ls =  220; ee.evt = 431740861;  e.push_back(ee);
  ee.run = 278239; ee.ls =  467; ee.evt = 772807981;  e.push_back(ee);
  ee.run = 278018; ee.ls =  595; ee.evt = 1095478828; e.push_back(ee);
  ee.run = 275001; ee.ls = 1052; ee.evt = 1464075959; e.push_back(ee);
  ee.run = 280015; ee.ls =   57; ee.evt = 98802122;   e.push_back(ee);
  ee.run = 277220; ee.ls =  255; ee.evt = 475950652;  e.push_back(ee);
  ee.run = 276327; ee.ls =   45; ee.evt = 55844247;   e.push_back(ee);
  
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");


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

      bool checkevent = false;
      for(unsigned int i = 0; i<e.size();++i){
	if(e[i].run!=tas::run() ) continue;
	if(e[i].ls !=tas::lumi()) continue;
	if(e[i].evt!=tas::evt() ) continue;
	checkevent = true;
	break;
      }
      if(!checkevent) continue;
      //checkevent = false;

      
      double weight = evt_scale1fb()*35.9;
      //if(nEventsTotal==0) cout << weight << endl; 

      if(nVert()<0) {
	if(checkevent) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " fails NVertex cut " << nVert() << endl;
	continue;
      }
      if(nlep()<2) {
	if(checkevent) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " fails N(stored leps) " << nlep() << endl;
	continue;
      }
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
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool isSS = false;  bool is3l = false;
	bool isaSS = false; bool isa3l = false;
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
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EA()[i]<0.2){
	      if(lep_p4()[i ].Pt()>20&&!is3l) { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EA()[i]<0.2){
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


      if(isData()){
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) continue;
	if( !goodrun(tas::run(), tas::lumi()) ){
	  if(checkevent) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " is not part of the good runlist json" << endl;
	  continue;
	}
	bool passonlineTrigger = false;
	if(nmu>=2&&(HLT_DoubleMu()||HLT_singleMu()) )             passonlineTrigger = true;
	if(nmu25>=1&&nel>=1&&(HLT_MuEG()||HLT_singleMu()))        passonlineTrigger = true;
	if(nel25>=1&&nmu>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) passonlineTrigger = true;
	if(!passonlineTrigger) {
	  if(checkevent) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " fails trigger: Nmu25/Nmu20 = " << nmu25 << "/" << nmu << " Nel25/Nel = " << nel25 << "/" << nel << " HLT: DiMu/SingleMu " << HLT_DoubleMu() << "/" << HLT_singleMu() << " EMu " << HLT_MuEG() << " DiEl(dZ) " << HLT_DoubleEl() << "(" << HLT_DoubleEl_DZ() << ")" << endl;
	  continue;
	}
      }

      if((!passofflineforTrigger)||(n3l<2&&(nSS+naSS)<2)||((nSS+naSS)==0&&lep_p4()[i3l[0] ].Pt()<25)){
	if(checkevent){
	  if(n3l<2&&(nSS+naSS)<2) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " fails minimum lepton requirement: n3l/(nSS+naSS) = " << n3l << "/" << (nSS+naSS) << endl;
	  else if((nSS+naSS)==0&&lep_p4()[i3l[0] ].Pt()<25) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " has no lepton with pT > 25 GeV: " << lep_p4()[i3l[0] ].Pt() << "/" << (nSS+naSS) << endl;
	  else if(!passofflineforTrigger) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " fails offline trigger cuts: Nmu25/Nmu20 = " << nmu25 << "/" << nmu << " Nel25/Nel = " << nel25 << "/" << nel << " nSS+naSS = " << (nSS+naSS) << endl;
	  for(unsigned int i = 0; i<lep_pdgId().size();++i){
	    cout << "lep" << i << " id " << lep_pdgId()[i] << " pT " << lep_p4()[i ].Pt() << " eta(SC) " << lep_p4()[i].Eta() << "(" << lep_etaSC()[i] <<")" << " IP3D " << lep_ip3d()[i] << " relIso " << lep_relIso03EA()[i] << " tightcharge " << lep_tightCharge()[i] << " cutbased veto/fo/tight " << lep_pass_VVV_cutbased_veto_noiso()[i] << "/" << lep_pass_VVV_cutbased_fo_noiso()[i] << "/" << lep_pass_VVV_cutbased_tight_noiso()[i] << endl;
	  }
	  continue;
	}
      }
      if(nb!=0) {
	if(checkevent) {
	  cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " have >=1 b jets: " << nb << endl;
	  for(unsigned int n = 0; n<jets_csv().size();++n){
	    cout << "jet" << n << " pT " << jets_p4()[n].Pt() << "(>20) eta " << jets_p4()[n].Eta() << "(<2.4) CSVv2 " << jets_csv()[n] << "(>0.5426)" << endl;
	  }
	}	
	continue;
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
	if(iSS.size()>=2||(iSS.size()==1&&iaSS.size()>=1)){
	  int nW(0), nZ(0), nG(0), nF(0);
	  if(lep_isFromW()[iSS[0] ]) ++nW;
	  else if(lep_isFromZ()[iSS[0] ]) ++nZ;
	  else if(lep_isFromB()[iSS[0] ]||lep_isFromC()[iSS[0] ]||lep_isFromL()[iSS[0] ]||lep_isFromLF()[iSS[0] ]) ++nF;
	  else if(lep_motherIdSS()[iSS[0] ]) ++nG;
	  if(iSS.size()>=2){
	    if(lep_isFromW()[iSS[1] ]) ++nW;
	    else if(lep_isFromZ()[iSS[1] ]) ++nZ;
	    else if(lep_isFromB()[iSS[1] ]||lep_isFromC()[iSS[1] ]||lep_isFromL()[iSS[1] ]||lep_isFromLF()[iSS[1] ]) ++nF;
	    else if(lep_motherIdSS()[iSS[1] ]) ++nG;
	  } else if(iSS.size()==1&&iaSS.size()>=1){
	    if(lep_isFromW()[iaSS[0] ]) ++nW;
	    else if(lep_isFromZ()[iaSS[0] ]) ++nZ;
	    else if(lep_isFromB()[iaSS[0] ]||lep_isFromC()[iaSS[0] ]||lep_isFromL()[iaSS[0] ]||lep_isFromLF()[iaSS[0] ]) ++nF;
	    else if(lep_motherIdSS()[iaSS[0] ]) ++nG;
	  }
	  if(iSS.size()>=2&&nW==2&&lep_mc_Id()[iSS[0] ]*lep_mc_Id()[iSS[1] ]>0) sn = "trueSS";//W+W+
	  else if(nW==2&&lep_mc_Id()[iSS[0] ]*lep_mc_Id()[iaSS[0] ]>0) sn = "trueSS";//W+W+
	  else if(nW==2) sn = "chargeflips";//W+W-
	  else if(iSS.size()>=2&&nZ==2&&lep_mc_Id()[iSS[0] ]*lep_mc_Id()[iSS[1] ]<=0) sn = "chargeflips";//Z
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
      }
      float MTmax = -1;
      if(iSS.size()>=2){
	if(mT(lep_p4()[iSS[1] ],MET)>mT(lep_p4()[iSS[0] ],MET)) MTmax = mT(lep_p4()[iSS[1] ],MET);
	else MTmax = mT(lep_p4()[iSS[0] ],MET);
      }
      else if(iSS.size()==1&&iaSS.size()>=1){
	if(mT(lep_p4()[iaSS[0] ],MET)>mT(lep_p4()[iSS[0] ],MET)) MTmax = mT(lep_p4()[iaSS[0] ],MET);
	else MTmax = mT(lep_p4()[iSS[0] ],MET);
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
      //if(ee[1]) cout << "SRee " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      //if(em[1]) cout << "SRem " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      //if(mm[1]) cout << "SRmm " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      //2: AR
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoaSS==0&&nSS==1&&naSS==1&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iaSS[0] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M()-90.)>10.) ee[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iaSS[0] ])==11&&met_pt()>40.&&MTmax>90.)                                                em[2] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iaSS[0] ])==13&&met_pt()>40.&&MTmax>90.)                                                em[2] = true;
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
 
      int SFOS[50];
      for(int i = 0; i<50; ++i) { SFOS[i] = -1; }
      double pTlll(-1), DPhilllMET(-1); double Mmumu(-1), Mmumue(-1);
      if(i3l.size()>=3){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
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
	if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) SFOScounter = -1;
	if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) SFOScounter = -1;
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
	  if((DPhilllMET>2.7&&pTlll>60.))                SFOS[0] = 0;
	  if(!(DPhilllMET>2.7&&pTlll>60.))               SFOS[2] = 0;
	}
	if(pass1){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.))  SFOS[0] = 1;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[2] = 1;
	}
	if(pass2){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.))  SFOS[0] = 2;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[2] = 2;
	}
	if(pass0X){
	  if((DPhilllMET>2.7&&pTlll>60.))                SFOS[1] = 0;
	  if(!(DPhilllMET>2.7&&pTlll>60.))               SFOS[3] = 0;
	}
	if(pass1X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.))  SFOS[1] = 1;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>45.)) SFOS[3] = 1;
	}
	if(pass2X){
	  if((DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.))  SFOS[1] = 2;
	  if(!(DPhilllMET>2.5&&pTlll>60.&&met_pt()>55.)) SFOS[3] = 2;
	}
      }

      if(checkevent){
	cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " is " << (ee[0] ? " SRee " : "") << (em[0] ? " SRem " : "") << (mm[0] ? " SRmm " : "") << (ee[2] ? " ARee " : "") << (em[2] ? " ARem " : "") << (mm[2] ? " ARmm " : "") << " or has " << SFOS[0] << " SRSFOS pairs" << endl;
	cout << "NJ(SS) " << nj30 << " NJ(3l) " << nj << " NB " << nb << " nisotrack " << nisoTrack_mt2_cleaned_VVV_cutbased_veto() << " nSS(nARSS) " << nSS << "(" << naSS << ") nvetoSS(AR) " << nvetoSS << "(" << nvetoaSS << ")" << " n3l " << n3l << " nveto3l " << nveto3l << " MET " << met_pt() << endl;
	if(nSS==2) cout << "SR-SS id1/2 " << lep_pdgId()[iSS[0] ] << "/" << lep_pdgId()[iSS[1] ] << " Mll " << (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M() << " MTmax " << MTmax << " Mjj " << Mjj << " MjjL " << MjjL << " Detajj " << Detajj << endl;
	if(nSS==1&&naSS==1) cout << "AR-SS id1/2 " << lep_pdgId()[iSS[0] ] << "/" << lep_pdgId()[iaSS[0] ] << " Mll " << (lep_p4()[iSS[0] ]+lep_p4()[iaSS[0] ]).M() << " MTmax " << MTmax << " Mjj " << Mjj << " MjjL " << MjjL << " Detajj " << Detajj << endl;
	if(n3l>=3) {
	  cout << "SR-SFOS id1/2/3 " << lep_pdgId()[i3l[0] ] << "/" << lep_pdgId()[i3l[1] ] << "/" << lep_pdgId()[i3l[2] ] << " DPhilllMET " << DPhilllMET << " pTlll " << pTlll << " Mlll " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M() << " Ml0l1 " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M() << " Ml0l2 " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M() << " Ml1l2 " << (lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M() << endl;
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
	  cout << "NSFOS " << SFOScounter << endl;
	  if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) SFOScounter = -1;
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) SFOScounter = -1;
	  if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) SFOScounter = -1;
	  if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) SFOScounter = -1;
	  cout << "but event vetoed because of gamma fake in wjets/dy or 3charge agreement" << endl;
	}
	for(unsigned int i = 0; i<lep_pdgId().size();++i){
	  cout << "lep" << i << " id " << lep_pdgId()[i] << " pT " << lep_p4()[i ].Pt() << " eta(SC) " << lep_p4()[i].Eta() << "(" << lep_etaSC()[i] <<")" << " IP3D " << lep_ip3d()[i] << " relIso " << lep_relIso03EA()[i] << " tightcharge " << lep_tightCharge()[i] << " cutbased veto/fo/tight " << lep_pass_VVV_cutbased_veto_noiso()[i] << "/" << lep_pass_VVV_cutbased_fo_noiso()[i] << "/" << lep_pass_VVV_cutbased_tight_noiso()[i] << endl;
	}
	for(unsigned int n = 0; n<jets_csv().size();++n){
	  cout << "jet" << n << " pT " << jets_p4()[n].Pt() << "(>20) eta " << jets_p4()[n].Eta() << "(<2.4) CSVv2 " << jets_csv()[n] << "(>0.5426)" << endl;
	}
      }//checkevent
    
    }//event loop
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }


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
