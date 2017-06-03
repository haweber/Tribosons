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

  histonames.push_back("YieldsSyncTighttest");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSyncTightIso");     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync");             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync3025_252020");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync3020_252015");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync3015_252010");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync3030_251515");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync3030_251510");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync2525_202015");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync2520_202010");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync2515_201515");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync2525_201510");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSync2525_201010");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSyncTighttest"); hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSyncTightIso");  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync3025_252020");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync3020_252015");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync3015_252010");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync3030_251515");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync3030_251510");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync2525_202015");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync2520_202010");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync2515_201515");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync2525_201510");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("RawYieldsSync2525_201010");    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  //ee,em,mm,0SFOS,1SFOS,2SFOS

  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
	string mapname = histonames[i] + "_"+skimFilePrefix;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
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
      //if(nisoTrack_mt2()!=0)  continue;

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

      int nj(0),nb(0);
      int nj20(0),nj30(0);
      vector<int> i2p5;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	if(tas::run()==827&&tas::lumi()==2&&tas::evt()==83) cout << currentFile->GetTitle() << " " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " j" << n << " pT " << jets_p4()[n].Pt() << " eta " << jets_p4()[n].Eta() << endl;
	if(fabs(jets_p4()[n].Eta())<2.5) i2p5.push_back(n);
	if(jets_p4()[n].Pt()>25&&fabs(jets_p4()[n].Eta())<4.5) ++nj;
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
      vector<int> iSStest, i3ltest, v, iSS, i3l, iSSv2, i3lv2, iSS25, i3l25, iSS20, i3l15, iSS15, i3l10;
      if(tas::run()==827&&tas::lumi()==2&&tas::evt()==83) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " NJ " << nj30 << " NB " << nb << " Deta " << Detajj << " MjjL " <<MjjL << " Mjj " << Mjj << endl;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	if(tas::run()==827&&tas::lumi()==2&&tas::evt()==83) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " id " << lep_pdgId()[i] << " 3q " << lep_tightCharge()[i] << " reliso " << lep_relIso03EA()[i] << " ip3d " << lep_ip3d()[i] << " pT " << lep_p4()[i ].Pt() << endl;
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSS.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&                         lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSS.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&/*lep_tightCharge()[i]==2&&*/lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3l.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&                             lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3l.push_back(i);
	/*
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso04DB()[i]<0.15   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso04DB()[i]<0.15   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso04DB()[i]<0.15   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso04DB()[i]<0.15   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	*/
	/*
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.12   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.12   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.12   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.12   && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	*/
	/*
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.1    && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.1    && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.1    && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.1    && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	*/
	
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	
	//bool passMultiIsoCuts(float cutMiniIso, float cutPtRatio, float cutPtRel, float miniIsoValue, float ptRatioValue, float ptRelValue){
  //  MiniRelIsoCMS3_EA<0.16 && (ptRat>0.76 || ptRel > 7.2) //muons
  //  MiniRelIsoCMS3_EA<0.12 && (ptRat>0.80 || ptRel > 7.2) //electrons


	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) iSS25.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) iSS25.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) iSS25.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) iSS25.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) i3l25.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) i3l25.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) i3l25.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>25) i3l25.push_back(i);

	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) iSS20.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) iSS20.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) iSS20.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) iSS20.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) i3l15.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) i3l15.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) i3l15.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) i3l15.push_back(i);

	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) iSS15.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) iSS15.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) iSS15.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>15) iSS15.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>10) i3l10.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>10) i3l10.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>10) i3l10.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>10) i3l10.push_back(i);
	
	// (miniIsoValue < cutMiniIso && (ptRatioValue>cutPtRatio || ptRelValue > cutPtRel))
	  
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&lep_tightCharge()[i]==2&&lep_miniRelIsoCMS3_EA()[i]<0.12&&(lep_ptRatio()[i]>0.80||lep_ptRel()[i]>7.2) && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSStest.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&                         lep_miniRelIsoCMS3_EA()[i]<0.16&&(lep_ptRatio()[i]>0.76||lep_ptRel()[i]>7.2) && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSStest.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&/*lep_tightCharge()[i]==2&&*/lep_miniRelIsoCMS3_EA()[i]<0.12&&(lep_ptRatio()[i]>0.80||lep_ptRel()[i]>7.2) && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3ltest.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&                             lep_miniRelIsoCMS3_EA()[i]<0.16&&(lep_ptRatio()[i]>0.76||lep_ptRel()[i]>7.2) && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3ltest.push_back(i);
	
	if(lep_relIso03EA()[i]<=0.3) v.push_back(i);
      }
      bool Bobak1(false), Bobak2(false);

      bool ee[7],mm[7],em[7];
      for(int i = 0; i<7; ++i) { ee[i] = false; em[i] = false; mm[i] = false; }
      //0 - select 2 signal leptons i1 + no taus
      //1 - select 2 signal leptons i1 + no taus + no vetotracks + no veto leptons + no loose b
      //2 - select 2 signal leptons i2 + no taus + no vetotracks + no veto leptons + no loose b
      //3 - select 2 signal leptons i3 + no taus + no vetotracks + no veto leptons + no loose b
      if(tas::run()==827&&tas::lumi()==2&&tas::evt()==83) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " isotrack " << nisoTrack_mt2() << " met " << met_pt() << " M " << (iSS.size()>=2 ? (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M() : -1.) << endl;
      //if(tas::run()==827&&tas::lumi()==2&&tas::evt()==83) cout << " nj30 " << nj30 << " nj20 " << nj20 << " Mjjpass " << passMDetajj << " nisotr " << nisoTrack_mt2() << " nlep " << nlep() << " nL " << iSS.size() << " SS " << ((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0) << " pT1/2 " << lep_p4()[iSS[0] ].Pt() << "/" << lep_p4()[iSS[1] ].Pt() << " M " << (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M() << endl; 
      if(nj30>1&&nj20>1&&nb==0&&passMDetajj&&nisoTrack_mt2()==0&&nlep()==2&&iSS.size()==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&lep_p4()[iSS[0] ].Pt()>30&&lep_p4()[iSS[1] ].Pt()>30&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[0] ]==2&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>40.)                                                                                         em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&lep_tightCharge()[iSS[0] ]==2&&met_pt()>40.)                                                                                         em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                                                                                      mm[0] = true;
      }
      if(nj30>1&&nj20>1&&nb==0&&passMDetajj&&nisoTrack_mt2()==0&&nlep()==2&&iSSv2.size()==2&&((lep_pdgId()[iSSv2[0] ]*lep_pdgId()[iSSv2[1] ])>0)&&lep_p4()[iSSv2[0] ].Pt()>30&&lep_p4()[iSSv2[1] ].Pt()>30&&(lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ]).M()>40.){
	if(abs(lep_pdgId()[iSSv2[0] ])==11&&abs(lep_pdgId()[iSSv2[1] ])==11&&lep_tightCharge()[iSSv2[0] ]==2&&lep_tightCharge()[iSSv2[1] ]==2&&met_pt()>40.&&fabs((lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ]).M()-90.)>10.) ee[1] = true;
	if(abs(lep_pdgId()[iSSv2[0] ])==13&&abs(lep_pdgId()[iSSv2[1] ])==11&&lep_tightCharge()[iSSv2[1] ]==2&&met_pt()>40.)                                                                                               em[1] = true;
	if(abs(lep_pdgId()[iSSv2[0] ])==11&&abs(lep_pdgId()[iSSv2[1] ])==13&&lep_tightCharge()[iSSv2[0] ]==2&&met_pt()>40.)                                                                                               em[1] = true;
	if(abs(lep_pdgId()[iSSv2[0] ])==13&&abs(lep_pdgId()[iSSv2[1] ])==13)                                                                                                                                              mm[1] = true;
      }
      if(nj30>1&&nj20>1&&nb==0&&passMDetajj&&nisoTrack_mt2()==0&&nlep()==2&&iSStest.size()==2&&((lep_pdgId()[iSStest[0] ]*lep_pdgId()[iSStest[1] ])>0)&&lep_p4()[iSStest[0] ].Pt()>30&&lep_p4()[iSStest[1] ].Pt()>30&&(lep_p4()[iSStest[0] ]+lep_p4()[iSStest[1] ]).M()>40.){
	if(abs(lep_pdgId()[iSStest[0] ])==11&&abs(lep_pdgId()[iSStest[1] ])==11&&lep_tightCharge()[iSStest[0] ]==2&&lep_tightCharge()[iSStest[1] ]==2&&met_pt()>40.&&fabs((lep_p4()[iSStest[0] ]+lep_p4()[iSStest[1] ]).M()-90.)>10.) ee[2] = true;
	if(abs(lep_pdgId()[iSStest[0] ])==13&&abs(lep_pdgId()[iSStest[1] ])==11&&lep_tightCharge()[iSStest[1] ]==2&&met_pt()>40.)                                                                                                     em[2] = true;
	if(abs(lep_pdgId()[iSStest[0] ])==11&&abs(lep_pdgId()[iSStest[1] ])==13&&lep_tightCharge()[iSStest[0] ]==2&&met_pt()>40.)                                                                                                     em[2] = true;
	if(abs(lep_pdgId()[iSStest[0] ])==13&&abs(lep_pdgId()[iSStest[1] ])==13)                                                                                                                                                      mm[2] = true;
      }
      if(nj30>1&&nj20>1&&nb==0&&passMDetajj&&nisoTrack_mt2()==0&&nlep()==2&&iSS15.size()==2&&((lep_pdgId()[iSS15[0] ]*lep_pdgId()[iSS15[1] ])>0)&&lep_p4()[iSS15[0] ].Pt()>15&&lep_p4()[iSS15[1] ].Pt()>15&&(lep_p4()[iSS15[0] ]+lep_p4()[iSS15[1] ]).M()>40.){
	if(abs(lep_pdgId()[iSS15[0] ])==11&&abs(lep_pdgId()[iSS15[1] ])==11&&lep_tightCharge()[iSS15[0] ]==2&&lep_tightCharge()[iSS15[1] ]==2&&met_pt()>40.&&fabs((lep_p4()[iSS15[0] ]+lep_p4()[iSS15[1] ]).M()-90.)>10.) ee[3] = true;
	if(abs(lep_pdgId()[iSS15[0] ])==13&&abs(lep_pdgId()[iSS15[1] ])==11&&lep_tightCharge()[iSS15[1] ]==2&&met_pt()>40.)                                                                                               em[3] = true;
	if(abs(lep_pdgId()[iSS15[0] ])==11&&abs(lep_pdgId()[iSS15[1] ])==13&&lep_tightCharge()[iSS15[0] ]==2&&met_pt()>40.)                                                                                               em[3] = true;
	if(abs(lep_pdgId()[iSS15[0] ])==13&&abs(lep_pdgId()[iSS15[1] ])==13)                                                                                                                                              mm[3] = true;
      }

      int SFOS[7];
      for(int i = 0; i<7; ++i) { SFOS[i] = -1; }
      //0 - select 3 signal leptons i1 + no taus
      //1 - select 3 signal leptons i1 + no taus + no vetotracks + no veto leptons + no loose b
      //2 - select 3 signal leptons i2 + no taus + no vetotracks + no veto leptons + no loose b
      //3 - select 3 signal leptons i3 + no taus + no vetotracks + no veto leptons + no loose b
      if(nj<2&&nb==0/*&&nisoTrack_mt2()==0&&v.size()==3*/&&nlep()==3&&i3l.size()==3&&lep_p4()[i3l[0] ].Pt()>20&&lep_p4()[i3l[1] ].Pt()>20&&lep_p4()[i3l[2] ].Pt()>20&&dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET)>2.5){
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
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-91.)<15.) pass0 = false;
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-91.)<15.) pass0 = false;
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-91.)<15.) pass0 = false;
	}
	if(SFOScounter==1){
	  pass1 = true;
	  if(met_pt()<45) pass1=false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>56.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<111.) pass1 = false;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()>56.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<111.) pass1 = false;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()>56.&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<111.) pass1 = false;
	}
	if(SFOScounter==2){
	  //cout << __LINE__ << endl;
	  pass2 = true;
	  if(met_pt()<55) pass2=false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-91.)<20.) pass2 = false;
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-91.)<20.) pass2 = false;
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-91.)<20.) pass2 = false;
	}
	if(pass0) SFOS[0] = 0;
	if(pass1) SFOS[0] = 1;
	if(pass2) SFOS[0] = 2;
      }
      if(nj<2&&nb==0/*&&nisoTrack_mt2()==0&&v.size()==3*/&&nlep()==3&&i3lv2.size()==3&&lep_p4()[i3lv2[0] ].Pt()>20&&lep_p4()[i3lv2[1] ].Pt()>20&&lep_p4()[i3lv2[2] ].Pt()>20&&dPhi(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ],MET)>2.5){
	bool OS01 = (lep_pdgId()[i3lv2[0] ]*lep_pdgId()[i3lv2[1] ]<0);
	bool OS02 = (lep_pdgId()[i3lv2[0] ]*lep_pdgId()[i3lv2[2] ]<0);
	bool OS12 = (lep_pdgId()[i3lv2[1] ]*lep_pdgId()[i3lv2[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3lv2[0] ])==abs(lep_pdgId()[i3lv2[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3lv2[0] ])==abs(lep_pdgId()[i3lv2[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3lv2[1] ])==abs(lep_pdgId()[i3lv2[2] ]));
	int SFOScounter = 0;
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	if(abs(lep_charge()[i3lv2[0] ]+lep_charge()[i3lv2[1] ]+lep_charge()[i3lv2[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  if(SF01&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()<20.) pass0 = false;
	  if(SF01&&abs(lep_pdgId()[i3lv2[0] ])==11&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()-91.)<15.) pass0 = false;
	  if(SF02&&abs(lep_pdgId()[i3lv2[0] ])==11&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()-91.)<15.) pass0 = false;
	  if(SF12&&abs(lep_pdgId()[i3lv2[1] ])==11&&fabs((lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()-91.)<15.) pass0 = false;
	}
	if(SFOScounter==1){
	  pass1 = true;
	  if(met_pt()<45) pass1=false;
	  if(OS01&&SF01&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()>56.&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()<111.) pass1 = false;
	  if(OS02&&SF02&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()>56.&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()<111.) pass1 = false;
	  if(OS12&&SF12&&(lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()>56.&&(lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()<111.) pass1 = false;
	}
	if(SFOScounter==2){
	  //cout << __LINE__ << endl;
	  pass2 = true;
	  if(met_pt()<55) pass2=false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()-91.)<20.) pass2 = false;
	  if(OS02&&SF02&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()-91.)<20.) pass2 = false;
	  if(OS12&&SF12&&fabs((lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()-91.)<20.) pass2 = false;
	}
	if(pass0) SFOS[1] = 0;
	if(pass1) SFOS[1] = 1;
	if(pass2) SFOS[1] = 2;
      }
      if(nj<2&&nb==0/*&&nisoTrack_mt2()==0&&v.size()==3*/&&nlep()==3&&i3ltest.size()==3&&lep_p4()[i3ltest[0] ].Pt()>20&&lep_p4()[i3ltest[1] ].Pt()>20&&lep_p4()[i3ltest[2] ].Pt()>20&&dPhi(lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[1] ]+lep_p4()[i3ltest[2] ],MET)>2.5){
	bool OS01 = (lep_pdgId()[i3ltest[0] ]*lep_pdgId()[i3ltest[1] ]<0);
	bool OS02 = (lep_pdgId()[i3ltest[0] ]*lep_pdgId()[i3ltest[2] ]<0);
	bool OS12 = (lep_pdgId()[i3ltest[1] ]*lep_pdgId()[i3ltest[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3ltest[0] ])==abs(lep_pdgId()[i3ltest[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3ltest[0] ])==abs(lep_pdgId()[i3ltest[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3ltest[1] ])==abs(lep_pdgId()[i3ltest[2] ]));
	int SFOScounter = 0;
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	if(abs(lep_charge()[i3ltest[0] ]+lep_charge()[i3ltest[1] ]+lep_charge()[i3ltest[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  if(SF01&&(lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3ltest[1] ]+lep_p4()[i3ltest[2] ]).M()<20.) pass0 = false;
	  if(SF01&&abs(lep_pdgId()[i3ltest[0] ])==11&&fabs((lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[1] ]).M()-91.)<15.) pass0 = false;
	  if(SF02&&abs(lep_pdgId()[i3ltest[0] ])==11&&fabs((lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[2] ]).M()-91.)<15.) pass0 = false;
	  if(SF12&&abs(lep_pdgId()[i3ltest[1] ])==11&&fabs((lep_p4()[i3ltest[1] ]+lep_p4()[i3ltest[2] ]).M()-91.)<15.) pass0 = false;
	}
	if(SFOScounter==1){
	  pass1 = true;
	  if(met_pt()<45) pass1=false;
	  if(OS01&&SF01&&(lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[1] ]).M()>56.&&(lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[1] ]).M()<111.) pass1 = false;
	  if(OS02&&SF02&&(lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[2] ]).M()>56.&&(lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[2] ]).M()<111.) pass1 = false;
	  if(OS12&&SF12&&(lep_p4()[i3ltest[1] ]+lep_p4()[i3ltest[2] ]).M()>56.&&(lep_p4()[i3ltest[1] ]+lep_p4()[i3ltest[2] ]).M()<111.) pass1 = false;
	}
	if(SFOScounter==2){
	  //cout << __LINE__ << endl;
	  pass2 = true;
	  if(met_pt()<55) pass2=false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[1] ]).M()-91.)<20.) pass2 = false;
	  if(OS02&&SF02&&fabs((lep_p4()[i3ltest[0] ]+lep_p4()[i3ltest[2] ]).M()-91.)<20.) pass2 = false;
	  if(OS12&&SF12&&fabs((lep_p4()[i3ltest[1] ]+lep_p4()[i3ltest[2] ]).M()-91.)<20.) pass2 = false;
	}
	if(pass0) SFOS[2] = 0;
	if(pass1) SFOS[2] = 1;
	if(pass2) SFOS[2] = 2;
      }

      if(nj<2&&nb==0/*&&nisoTrack_mt2()==0&&v.size()==3*/&&nlep()==3&&i3l10.size()==3&&lep_p4()[i3l10[0] ].Pt()>20&&lep_p4()[i3l10[1] ].Pt()>10&&lep_p4()[i3l10[2] ].Pt()>10&&dPhi(lep_p4()[i3l10[0] ]+lep_p4()[i3l10[1] ]+lep_p4()[i3l10[2] ],MET)>2.5){
	bool OS01 = (lep_pdgId()[i3l10[0] ]*lep_pdgId()[i3l10[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l10[0] ]*lep_pdgId()[i3l10[2] ]<0);
	bool OS12 = (lep_pdgId()[i3l10[1] ]*lep_pdgId()[i3l10[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l10[0] ])==abs(lep_pdgId()[i3l10[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l10[0] ])==abs(lep_pdgId()[i3l10[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3l10[1] ])==abs(lep_pdgId()[i3l10[2] ]));
	int SFOScounter = 0;
	if(OS01&&SF01) ++SFOScounter;
	if(OS02&&SF02) ++SFOScounter;
	if(OS12&&SF12) ++SFOScounter;
	bool pass0(false), pass1(false), pass2(false);
	if(abs(lep_charge()[i3l10[0] ]+lep_charge()[i3l10[1] ]+lep_charge()[i3l10[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  if(SF01&&(lep_p4()[i3l10[0] ]+lep_p4()[i3l10[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3l10[0] ]+lep_p4()[i3l10[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3l10[1] ]+lep_p4()[i3l10[2] ]).M()<20.) pass0 = false;
	  if(SF01&&abs(lep_pdgId()[i3l10[0] ])==11&&fabs((lep_p4()[i3l10[0] ]+lep_p4()[i3l10[1] ]).M()-91.)<15.) pass0 = false;
	  if(SF02&&abs(lep_pdgId()[i3l10[0] ])==11&&fabs((lep_p4()[i3l10[0] ]+lep_p4()[i3l10[2] ]).M()-91.)<15.) pass0 = false;
	  if(SF12&&abs(lep_pdgId()[i3l10[1] ])==11&&fabs((lep_p4()[i3l10[1] ]+lep_p4()[i3l10[2] ]).M()-91.)<15.) pass0 = false;
	}
	if(SFOScounter==1){
	  pass1 = true;
	  if(met_pt()<45) pass1=false;
	  if(OS01&&SF01&&(lep_p4()[i3l10[0] ]+lep_p4()[i3l10[1] ]).M()>56.&&(lep_p4()[i3l10[0] ]+lep_p4()[i3l10[1] ]).M()<111.) pass1 = false;
	  if(OS02&&SF02&&(lep_p4()[i3l10[0] ]+lep_p4()[i3l10[2] ]).M()>56.&&(lep_p4()[i3l10[0] ]+lep_p4()[i3l10[2] ]).M()<111.) pass1 = false;
	  if(OS12&&SF12&&(lep_p4()[i3l10[1] ]+lep_p4()[i3l10[2] ]).M()>56.&&(lep_p4()[i3l10[1] ]+lep_p4()[i3l10[2] ]).M()<111.) pass1 = false;
	}
	if(SFOScounter==2){
	  //cout << __LINE__ << endl;
	  pass2 = true;
	  if(met_pt()<55) pass2=false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l10[0] ]+lep_p4()[i3l10[1] ]).M()-91.)<20.) pass2 = false;
	  if(OS02&&SF02&&fabs((lep_p4()[i3l10[0] ]+lep_p4()[i3l10[2] ]).M()-91.)<20.) pass2 = false;
	  if(OS12&&SF12&&fabs((lep_p4()[i3l10[1] ]+lep_p4()[i3l10[2] ]).M()-91.)<20.) pass2 = false;
	}
	if(pass0) SFOS[3] = 0;
	if(pass1) SFOS[3] = 1;
	if(pass2) SFOS[3] = 2;
      }

      
      if(ee[0])      *eventstocheckEE    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(em[0])      *eventstocheckEM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(mm[0])      *eventstocheckMM    << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(SFOS[0]==0) *eventstocheck0SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(SFOS[0]==1) *eventstocheck1SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
      if(SFOS[0]==2) *eventstocheck2SFOS << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;

      if(ee[0])       histos["YieldsSync_"     +skimFilePrefix]->Fill(0.,weight);
      if(em[0])       histos["YieldsSync_"     +skimFilePrefix]->Fill(1.,weight);
      if(mm[0])       histos["YieldsSync_"     +skimFilePrefix]->Fill(2.,weight);
      if(SFOS[0]==0)  histos["YieldsSync_"     +skimFilePrefix]->Fill(3.,weight);
      if(SFOS[0]==1)  histos["YieldsSync_"     +skimFilePrefix]->Fill(4.,weight);
      if(SFOS[0]==2)  histos["YieldsSync_"     +skimFilePrefix]->Fill(5.,weight);
      if(ee[1])       histos["YieldsSyncTightIso_"     +skimFilePrefix]->Fill(0.,weight);
      if(em[1])       histos["YieldsSyncTightIso_"     +skimFilePrefix]->Fill(1.,weight);
      if(mm[1])       histos["YieldsSyncTightIso_"     +skimFilePrefix]->Fill(2.,weight);
      if(SFOS[1]==0)  histos["YieldsSyncTightIso_"     +skimFilePrefix]->Fill(3.,weight);
      if(SFOS[1]==1)  histos["YieldsSyncTightIso_"     +skimFilePrefix]->Fill(4.,weight);
      if(SFOS[1]==2)  histos["YieldsSyncTightIso_"     +skimFilePrefix]->Fill(5.,weight);
      if(ee[2])       histos["YieldsSyncTighttest_"     +skimFilePrefix]->Fill(0.,weight);
      if(em[2])       histos["YieldsSyncTighttest_"     +skimFilePrefix]->Fill(1.,weight);
      if(mm[2])       histos["YieldsSyncTighttest_"     +skimFilePrefix]->Fill(2.,weight);
      if(SFOS[2]==0)  histos["YieldsSyncTighttest_"     +skimFilePrefix]->Fill(3.,weight);
      if(SFOS[2]==1)  histos["YieldsSyncTighttest_"     +skimFilePrefix]->Fill(4.,weight);
      if(SFOS[2]==2)  histos["YieldsSyncTighttest_"     +skimFilePrefix]->Fill(5.,weight);
      if(iSS.size()>=1&&iSS25.size()>=2){
	if(ee[3])       histos["YieldsSync3025_252020_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync3025_252020_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync3025_252020_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync3025_252020_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync3025_252020_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync3025_252020_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=1&&iSS20.size()>=2){
	if(ee[3])       histos["YieldsSync3020_252015_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync3020_252015_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync3020_252015_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l.size()>=2&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync3020_252015_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync3020_252015_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync3020_252015_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=1&&iSS15.size()>=2){
	if(ee[3])       histos["YieldsSync3015_252010_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync3015_252010_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync3015_252010_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync3015_252010_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync3015_252010_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync3015_252010_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=2){
	if(ee[3])       histos["YieldsSync3030_251515_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync3030_251515_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync3030_251515_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync3030_251515_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync3030_251515_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync3030_251515_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=2){
	if(ee[3])       histos["YieldsSync3030_251510_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync3030_251510_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync3030_251510_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l15.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync3030_251510_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync3030_251510_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync3030_251510_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=2){
	if(ee[3])       histos["YieldsSync2525_202015_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync2525_202015_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync2525_202015_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=2&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync2525_202015_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync2525_202015_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync2525_202015_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=1&&iSS20.size()>=2){
	if(ee[3])       histos["YieldsSync2520_202010_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync2520_202010_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync2520_202010_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync2520_202010_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync2520_202010_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync2520_202010_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=1&&iSS15.size()>=2){
	if(ee[3])       histos["YieldsSync2515_201515_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync2515_201515_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync2515_201515_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=1&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync2515_201515_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync2515_201515_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync2515_201515_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=2){
	if(ee[3])       histos["YieldsSync2525_201510_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync2525_201510_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync2525_201510_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=1&&i3l15.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync2525_201510_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync2525_201510_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync2525_201510_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=2){
	if(ee[3])       histos["YieldsSync2525_201010_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["YieldsSync2525_201010_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["YieldsSync2525_201010_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=1&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["YieldsSync2525_201010_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["YieldsSync2525_201010_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["YieldsSync2525_201010_"     +skimFilePrefix]->Fill(5.,weight);
      }


      
      if(Bobak1) histos["RawBobak_"+skimFilePrefix]->Fill(0.,1);
      if(Bobak2) histos["RawBobak_"+skimFilePrefix]->Fill(1.,1);
      if(ee[0])       histos["RawYieldsSync_"     +skimFilePrefix]->Fill(0.,1);
      if(em[0])       histos["RawYieldsSync_"     +skimFilePrefix]->Fill(1.,1);
      if(mm[0])       histos["RawYieldsSync_"     +skimFilePrefix]->Fill(2.,1);
      if(SFOS[0]==0)  histos["RawYieldsSync_"     +skimFilePrefix]->Fill(3.,1);
      if(SFOS[0]==1)  histos["RawYieldsSync_"     +skimFilePrefix]->Fill(4.,1);
      if(SFOS[0]==2)  histos["RawYieldsSync_"     +skimFilePrefix]->Fill(5.,1);
      if(ee[1])       histos["RawYieldsSyncTightIso_"     +skimFilePrefix]->Fill(0.,1);
      if(em[1])       histos["RawYieldsSyncTightIso_"     +skimFilePrefix]->Fill(1.,1);
      if(mm[1])       histos["RawYieldsSyncTightIso_"     +skimFilePrefix]->Fill(2.,1);
      if(SFOS[1]==0)  histos["RawYieldsSyncTightIso_"     +skimFilePrefix]->Fill(3.,1);
      if(SFOS[1]==1)  histos["RawYieldsSyncTightIso_"     +skimFilePrefix]->Fill(4.,1);
      if(SFOS[1]==2)  histos["RawYieldsSyncTightIso_"     +skimFilePrefix]->Fill(5.,1);
      if(ee[2])       histos["RawYieldsSyncTighttest_"     +skimFilePrefix]->Fill(0.,1);
      if(em[2])       histos["RawYieldsSyncTighttest_"     +skimFilePrefix]->Fill(1.,1);
      if(mm[2])       histos["RawYieldsSyncTighttest_"     +skimFilePrefix]->Fill(2.,1);
      if(SFOS[2]==0)  histos["RawYieldsSyncTighttest_"     +skimFilePrefix]->Fill(3.,1);
      if(SFOS[2]==1)  histos["RawYieldsSyncTighttest_"     +skimFilePrefix]->Fill(4.,1);
      if(SFOS[2]==2)  histos["RawYieldsSyncTighttest_"     +skimFilePrefix]->Fill(5.,1);
      if(iSS.size()>=1&&iSS25.size()>=2){
	if(ee[3])       histos["RawYieldsSync3025_252020_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync3025_252020_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync3025_252020_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync3025_252020_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync3025_252020_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync3025_252020_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=1&&iSS20.size()>=2){
	if(ee[3])       histos["RawYieldsSync3020_252015_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync3020_252015_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync3020_252015_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l.size()>=2&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync3020_252015_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync3020_252015_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync3020_252015_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=1&&iSS15.size()>=2){
	if(ee[3])       histos["RawYieldsSync3015_252010_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync3015_252010_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync3015_252010_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync3015_252010_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync3015_252010_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync3015_252010_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=2){
	if(ee[3])       histos["RawYieldsSync3030_251515_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync3030_251515_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync3030_251515_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync3030_251515_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync3030_251515_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync3030_251515_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS.size()>=2){
	if(ee[3])       histos["RawYieldsSync3030_251510_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync3030_251510_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync3030_251510_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l25.size()>=1&&i3l15.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync3030_251510_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync3030_251510_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync3030_251510_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=2){
	if(ee[3])       histos["RawYieldsSync2525_202015_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync2525_202015_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync2525_202015_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=2&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync2525_202015_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync2525_202015_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync2525_202015_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=1&&iSS20.size()>=2){
	if(ee[3])       histos["RawYieldsSync2520_202010_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync2520_202010_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync2520_202010_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync2520_202010_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync2520_202010_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync2520_202010_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=1&&iSS15.size()>=2){
	if(ee[3])       histos["RawYieldsSync2515_201515_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync2515_201515_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync2515_201515_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=1&&i3l15.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync2515_201515_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync2515_201515_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync2515_201515_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=2){
	if(ee[3])       histos["RawYieldsSync2525_201510_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync2525_201510_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync2525_201510_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=1&&i3l15.size()>=2&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync2525_201510_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync2525_201510_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync2525_201510_"     +skimFilePrefix]->Fill(5.,weight);
      }
      if(iSS25.size()>=2){
	if(ee[3])       histos["RawYieldsSync2525_201010_"     +skimFilePrefix]->Fill(0.,weight);
	if(em[3])       histos["RawYieldsSync2525_201010_"     +skimFilePrefix]->Fill(1.,weight);
	if(mm[3])       histos["RawYieldsSync2525_201010_"     +skimFilePrefix]->Fill(2.,weight);
      } if(i3l.size()>=1&&i3l10.size()>=3){
	if(SFOS[3]==0)  histos["RawYieldsSync2525_201010_"     +skimFilePrefix]->Fill(3.,weight);
	if(SFOS[3]==1)  histos["RawYieldsSync2525_201010_"     +skimFilePrefix]->Fill(4.,weight);
	if(SFOS[3]==2)  histos["RawYieldsSync2525_201010_"     +skimFilePrefix]->Fill(5.,weight);
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
  /*
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //add overflow
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
    //add underflow
    h->second->SetBinContent(1, h->second->GetBinContent(1)+ h->second->GetBinContent(0) );
    h->second->SetBinError(1, sqrt(pow(h->second->GetBinError(1),2)+pow(h->second->GetBinError(0),2) ) );
  }
  */
  string filename = "rootfiles/StupidSync.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //string mapname = histonames[i]+"_"+skimFilePrefix;
    //string writename = histonames[i];
    h->second->Write(h->first.c_str(),TObject::kOverwrite);
  }  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;
  for(unsigned int i = 0; i< eventstocheck.size(); ++i) cout << eventstocheck[i] << endl;
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
