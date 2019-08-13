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

  histonames.push_back("MET");           hbins.push_back(25); hlow.push_back( 0); hup.push_back(250);
  histonames.push_back("MTclosest");     hbins.push_back(25); hlow.push_back( 0); hup.push_back(250);
  histonames.push_back("MTsum");         hbins.push_back(25); hlow.push_back( 0); hup.push_back(500);
  histonames.push_back("MTmin");         hbins.push_back(35); hlow.push_back( 0); hup.push_back(175);
  histonames.push_back("MTmax");         hbins.push_back(25); hlow.push_back( 0); hup.push_back(250);
  histonames.push_back("MTlsys");        hbins.push_back(25); hlow.push_back( 0); hup.push_back(175);
  histonames.push_back("ptlsum");        hbins.push_back(25); hlow.push_back(60); hup.push_back(310);
  histonames.push_back("ptlsys");        hbins.push_back(25); hlow.push_back( 0); hup.push_back(250);
  histonames.push_back("Mlsys");         hbins.push_back(35); hlow.push_back( 0); hup.push_back(350);
  histonames.push_back("MllSFOS");       hbins.push_back(35); hlow.push_back( 0); hup.push_back(350);
  histonames.push_back("DPhillSFOS");    hbins.push_back(32); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("DRllSFOS");      hbins.push_back(25); hlow.push_back( 0); hup.push_back(5);
  histonames.push_back("DPhillmin");     hbins.push_back(32); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("DPhillmax");     hbins.push_back(32); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("Mljclosestmin"); hbins.push_back(25); hlow.push_back( 0); hup.push_back(250);
  histonames.push_back("Mljclosestmax"); hbins.push_back(25); hlow.push_back( 0); hup.push_back(250);
  histonames.push_back("DPhilsysMET");   hbins.push_back(32); hlow.push_back( 0); hup.push_back(3.2);
  histonames.push_back("HT");            hbins.push_back(30); hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("ST");            hbins.push_back(30); hlow.push_back( 0); hup.push_back(600);
  histonames.push_back("ptjsys");        hbins.push_back(25); hlow.push_back( 0); hup.push_back(250);
  histonames.push_back("MT2ljmin");      hbins.push_back(35); hlow.push_back( 0); hup.push_back(350);
  //ee,em,mm,0SFOS,1SFOS,2SFOS

  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    for(int b = 0; b<6; ++b){
      string bb;
      if(b==0)      bb = "SSee";
      else if(b==1) bb = "SSem";
      else if(b==2) bb = "SSmm";
      else if(b==3) bb = "0SFOS";
      else if(b==4) bb = "1SFOS";
      else          bb = "2SFOS";
	string mapname = histonames[i] + "_"+bb+"_"+skimFilePrefix;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
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

	
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSSv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i])<=1.479&&                         lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==11&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_etaSC()[i]) >1.479&&                         lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta())<=1.479&&                      lep_relIso03EA()[i]<0.0588 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	if(abs(lep_pdgId()[i])==13&&fabs(lep_p4()[i].Eta())<2.5&&fabs(lep_p4()[i].Eta()) >1.479&&                      lep_relIso03EA()[i]<0.0571 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3lv2.push_back(i);
	
	
	if(lep_relIso03EA()[i]<=0.3) v.push_back(i);
      }
      bool Bobak1(false), Bobak2(false);

      bool ee[7],mm[7],em[7];
      for(int i = 0; i<7; ++i) { ee[i] = false; em[i] = false; mm[i] = false; }
      //0 - select 2 signal leptons i1 + no taus
      //1 - select 2 signal leptons i1 + no taus + no vetotracks + no veto leptons + no loose b
      //2 - select 2 signal leptons i2 + no taus + no vetotracks + no veto leptons + no loose b
      //3 - select 2 signal leptons i3 + no taus + no vetotracks + no veto leptons + no loose b
      bool passSelectionMET = true;
      bool passSelectionMll = true;
      bool passSelectionDPhi = true;

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
      if(nj30>1&&nj20>1&&nb==0&&passMDetajj&&nisoTrack_mt2()==0&&nlep()==2&&iSSv2.size()==2&&((lep_pdgId()[iSSv2[0] ]*lep_pdgId()[iSSv2[1] ])>0)&&lep_p4()[iSSv2[0] ].Pt()>30&&lep_p4()[iSSv2[1] ].Pt()>30){
	if((lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ]).M()<=40) passSelectionMll = false;
	if(abs(lep_pdgId()[iSSv2[0] ])==11&&abs(lep_pdgId()[iSSv2[1] ])==11&&lep_tightCharge()[iSSv2[0] ]==2&&lep_tightCharge()[iSSv2[1] ]==2) ee[2] = true;
	if(ee[2]&&fabs((lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ]).M()-90.)<=10.) passSelectionMll = false;
	if(ee[2]&&met_pt()<=40) passSelectionMET = false;
	if(abs(lep_pdgId()[iSSv2[0] ])==13&&abs(lep_pdgId()[iSSv2[1] ])==11&&lep_tightCharge()[iSSv2[1] ]==2)                                  em[2] = true;
	if(abs(lep_pdgId()[iSSv2[0] ])==11&&abs(lep_pdgId()[iSSv2[1] ])==13&&lep_tightCharge()[iSSv2[0] ]==2)                                  em[2] = true;
	if(em[2]&&met_pt()<=40) passSelectionMET = false;
	if(abs(lep_pdgId()[iSSv2[0] ])==13&&abs(lep_pdgId()[iSSv2[1] ])==13)                                                                   mm[2] = true;
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

      if(nj<2&&nb==0/*&&nisoTrack_mt2()==0&&v.size()==3*/&&nlep()==3&&i3lv2.size()==3&&lep_p4()[i3lv2[0] ].Pt()>20&&lep_p4()[i3lv2[1] ].Pt()>20&&lep_p4()[i3lv2[2] ].Pt()>20){
	if(dPhi(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ],MET)<=2.5) passSelectionDPhi = false;
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
	  if(SF01&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()<20.) passSelectionMll = false;
	  if(SF02&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()<20.) passSelectionMll = false;
	  if(SF12&&(lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()<20.) passSelectionMll = false;
	  if(SF01&&abs(lep_pdgId()[i3lv2[0] ])==11&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()-91.)<15.) passSelectionMll = false;
	  if(SF02&&abs(lep_pdgId()[i3lv2[0] ])==11&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()-91.)<15.) passSelectionMll = false;
	  if(SF12&&abs(lep_pdgId()[i3lv2[1] ])==11&&fabs((lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()-91.)<15.) passSelectionMll = false;
	}
	if(SFOScounter==1){
	  pass1 = true;
	  if(met_pt()<45) passSelectionMET=false;
	  if(OS01&&SF01&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()>56.&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()<111.) passSelectionMll = false;
	  if(OS02&&SF02&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()>56.&&(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()<111.) passSelectionMll = false;
	  if(OS12&&SF12&&(lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()>56.&&(lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()<111.) passSelectionMll = false;
	}
	if(SFOScounter==2){
	  //cout << __LINE__ << endl;
	  pass2 = true;
	  if(met_pt()<55) passSelectionMET=false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M()-91.)<20.) passSelectionMll = false;
	  if(OS02&&SF02&&fabs((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M()-91.)<20.) passSelectionMll = false;
	  if(OS12&&SF12&&fabs((lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M()-91.)<20.) passSelectionMll = false;
	}
	if(pass0) SFOS[2] = 0;
	if(pass1) SFOS[2] = 1;
	if(pass2) SFOS[2] = 2;
      }
      
      if(ee[2]||em[2]||mm[2]){
	string bb = "";
	if(ee[2]) bb = "_SSee_";
	if(em[2]) bb = "_SSem_";
	if(mm[2]) bb = "_SSmm_";
	if(passSelectionMll) histos["MET" + bb +skimFilePrefix]->Fill(met_pt(),weight);
	if(passSelectionMET) histos["Mlsys" + bb +skimFilePrefix]->Fill((lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ]).M(),weight);
	if(passSelectionMET&&passSelectionMll){
	  if(dPhi(MET,(lep_p4()[iSSv2[0] ]))<dPhi(MET,(lep_p4()[iSSv2[1] ]))) histos["MTclosest" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[0] ]),weight);
	  else                                                                histos["MTclosest" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[1] ]),weight);
	  if(mT(MET,(lep_p4()[iSSv2[0] ]))<mT(MET,(lep_p4()[iSSv2[1] ]))) {
	    histos["MTmin" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[0] ]),weight);
	    histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[1] ]),weight);
	  } else {
	    histos["MTmin" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[1] ]),weight);
	    histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[0] ]),weight);
	  }
	  histos["MTsum" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[0] ])+mT(MET,lep_p4()[iSSv2[1] ]),weight);
	  histos["MTlsys"+ bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ]),weight);
	  histos["ptlsys"+ bb +skimFilePrefix]->Fill((lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ]).Pt(),weight);
	  histos["ptlsum"+ bb +skimFilePrefix]->Fill(lep_p4()[iSSv2[0] ].Pt()+lep_p4()[iSSv2[1] ].Pt(),weight);
	  histos["DPhillmin" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[iSSv2[0] ],lep_p4()[iSSv2[1] ]),weight);
	  int l1c(-1), l2c(-1.); LorentzVector jsum(0,0,0,0);
	  float HT = 0; float MT2ljmin(99999);
	  for(unsigned int j = 0; j<i2p5.size();++j){
	    if(j<=1&&jets_p4()[i2p5[j] ].Pt()<30) continue;
	    if(jets_p4()[i2p5[j] ].Pt()<20) continue;
	    if(l1c<0) l1c = i2p5[j];
	    else if(dR(lep_p4()[iSSv2[0] ],jets_p4()[i2p5[j] ])<dR(lep_p4()[iSSv2[0] ],jets_p4()[l1c ])) l1c = i2p5[j];
	    if(l2c<0) l2c = i2p5[j];
	    else if(dR(lep_p4()[iSSv2[1] ],jets_p4()[i2p5[j] ])<dR(lep_p4()[iSSv2[1] ],jets_p4()[l2c ])) l2c = i2p5[j];
	    jsum = jsum + jets_p4()[i2p5[j] ];
	    HT += jets_p4()[i2p5[j] ].Pt();
	    for(unsigned int k = j+1; k<i2p5.size();++k){
	      if(k<=1&&jets_p4()[i2p5[k] ].Pt()<30) continue;
	      if(jets_p4()[i2p5[k] ].Pt()<20) continue;
	      float mymt2 = MT2(lep_p4()[iSSv2[0] ]+jets_p4()[i2p5[j] ],lep_p4()[iSSv2[1] ]+jets_p4()[i2p5[k] ],MET,true);    if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	      mymt2 = MT2(lep_p4()[iSSv2[0] ]+jets_p4()[i2p5[k] ],lep_p4()[iSSv2[1] ]+jets_p4()[i2p5[j] ],MET,true);          if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	    }
	  }
	  if(l1c>=0&&l2c>=0){
	    if((lep_p4()[iSSv2[0] ]+jets_p4()[l1c ]).M()<(lep_p4()[iSSv2[1] ]+jets_p4()[l2c ]).M()) {
	      histos["Mljclosestmin" + bb +skimFilePrefix]->Fill((lep_p4()[iSSv2[0] ]+jets_p4()[l1c ]).M(),weight);
	      histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[iSSv2[1] ]+jets_p4()[l2c ]).M(),weight);
	    } else {
	      histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[iSSv2[0] ]+jets_p4()[l1c ]).M(),weight);
	      histos["Mljclosestmin" + bb +skimFilePrefix]->Fill((lep_p4()[iSSv2[1] ]+jets_p4()[l2c ]).M(),weight);
	    }
	  }
	  histos["DPhilsysMET" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[iSSv2[0] ]+lep_p4()[iSSv2[1] ],MET),weight);
	  histos["HT" + bb +skimFilePrefix]->Fill(HT,weight);
	  histos["ST" + bb +skimFilePrefix]->Fill(HT+lep_p4()[iSSv2[0] ].Pt()+lep_p4()[iSSv2[1] ].Pt()+met_pt(),weight);
	  histos["ptjsys"+bb+skimFilePrefix]->Fill(jsum.Pt(),weight);
	  histos["MT2ljmin"+bb+skimFilePrefix]->Fill(MT2ljmin,weight);
	}
      }//SS done

      if(SFOS[2]>=0){
	string bb = "";
	if(SFOS[2]==0) bb = "_0SFOS_";
	if(SFOS[2]==1) bb = "_1SFOS_";
	if(SFOS[2]==2) bb = "_2SFOS_";
	bool OS01 = (lep_pdgId()[i3lv2[0] ]*lep_pdgId()[i3lv2[1] ]<0);
	bool OS02 = (lep_pdgId()[i3lv2[0] ]*lep_pdgId()[i3lv2[2] ]<0);
	bool OS12 = (lep_pdgId()[i3lv2[1] ]*lep_pdgId()[i3lv2[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3lv2[0] ])==abs(lep_pdgId()[i3lv2[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3lv2[0] ])==abs(lep_pdgId()[i3lv2[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3lv2[1] ])==abs(lep_pdgId()[i3lv2[2] ]));
	if(passSelectionMll&&passSelectionDPhi) histos["MET" + bb +skimFilePrefix]->Fill(met_pt(),weight);
	if(passSelectionMET&&passSelectionDPhi){
	  if(OS01&&SF01) histos["MllSFOS" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]).M(),weight);
	  if(OS02&&SF02) histos["MllSFOS" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[2] ]).M(),weight);
	  if(OS12&&SF12) histos["MllSFOS" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M(),weight);
	}
	if(passSelectionMll&&passSelectionMET) histos["DPhilsysMET" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ],MET),weight);

	if(passSelectionMET&&passSelectionMll&&passSelectionDPhi){
	  histos["Mlsys" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).M(),weight);
	  int imin(-1);
	  if((dPhi(MET,(lep_p4()[i3lv2[0] ]))<dPhi(MET,(lep_p4()[i3lv2[1] ]))) && (dPhi(MET,(lep_p4()[i3lv2[0] ]))<dPhi(MET,(lep_p4()[i3lv2[2] ]))) ) histos["MTclosest" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[0] ]),weight);
	  else if(dPhi(MET,(lep_p4()[i3lv2[1] ]))<dPhi(MET,(lep_p4()[i3lv2[2] ]))) histos["MTclosest" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[1] ]),weight);
	  else                                                                     histos["MTclosest" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[2] ]),weight);
	  if((mT(MET,(lep_p4()[i3lv2[0] ]))<mT(MET,(lep_p4()[i3lv2[1] ]))) && (mT(MET,(lep_p4()[i3lv2[0] ]))<mT(MET,(lep_p4()[i3lv2[2] ]))) ) {
	    histos["MTmin" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[0] ]),weight);
	    if(mT(MET,(lep_p4()[i3lv2[2] ]))<mT(MET,(lep_p4()[i3lv2[1] ]))) histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[1] ]),weight);
	    else                                                            histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[2] ]),weight);
	  } else if(mT(MET,(lep_p4()[i3lv2[2] ]))<mT(MET,(lep_p4()[i3lv2[1] ]))) {
	    histos["MTmin" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[2] ]),weight);
	    if(mT(MET,(lep_p4()[i3lv2[0] ]))<mT(MET,(lep_p4()[i3lv2[1] ]))) histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[1] ]),weight);
	    else                                                            histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[0] ]),weight);
	  } else {
	    histos["MTmin" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[1] ]),weight);
	    if(mT(MET,(lep_p4()[i3lv2[0] ]))<mT(MET,(lep_p4()[i3lv2[2] ]))) histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[2] ]),weight);
	    else                                                            histos["MTmax" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[0] ]),weight);
	  }
	  histos["MTsum" + bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[0] ])+mT(MET,lep_p4()[i3lv2[1] ])+mT(MET,lep_p4()[i3lv2[2] ]),weight);
	  histos["MTlsys"+ bb +skimFilePrefix]->Fill(mT(MET,lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]),weight);
	  histos["ptlsys"+ bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[0] ]+lep_p4()[i3lv2[1] ]+lep_p4()[i3lv2[2] ]).Pt(),weight);
	  histos["ptlsum"+ bb +skimFilePrefix]->Fill(lep_p4()[i3lv2[0] ].Pt()+lep_p4()[i3lv2[1] ].Pt()+lep_p4()[i3lv2[2] ].Pt(),weight);
	  if( (dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[0] ])<dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[1] ])) && (dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[0] ])<dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[1] ])) ){
	    histos["DPhillmin" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[2] ]),weight);
	    if (dPhi(lep_p4()[i3lv2[1] ],lep_p4()[i3lv2[0] ])<dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[1] ]))  histos["DPhillmax" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[1] ],lep_p4()[i3lv2[2] ]),weight);
	    else                                                                                              histos["DPhillmax" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[1] ]),weight);
	  } else if(dPhi(lep_p4()[i3lv2[1] ],lep_p4()[i3lv2[0] ])<dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[1] ])){
	    histos["DPhillmin" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[1] ]),weight);
	    if (dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[0] ])<dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[1] ]))  histos["DPhillmax" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[1] ],lep_p4()[i3lv2[2] ]),weight);
	    else                                                                                              histos["DPhillmax" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[2] ]),weight);
	  } else {
	    histos["DPhillmin" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[1] ],lep_p4()[i3lv2[2] ]),weight);
	    if (dPhi(lep_p4()[i3lv2[2] ],lep_p4()[i3lv2[0] ])<dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[1] ]))  histos["DPhillmax" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[1] ]),weight);
	    else                                                                                              histos["DPhillmax" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[2] ]),weight);
	  }
	  if(OS01&&SF01) {
	    histos["DPhillSFOS" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[1] ]),weight);
	    histos["DRllSFOS"   + bb +skimFilePrefix]->Fill(dR(  lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[1] ]),weight);
	  }
	  if(OS02&&SF02) {
	    histos["DPhillSFOS" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[2] ]),weight);
	    histos["DRllSFOS"   + bb +skimFilePrefix]->Fill(dR(  lep_p4()[i3lv2[0] ],lep_p4()[i3lv2[2] ]),weight);
	  }
	  if(OS12&&SF12) {
	    histos["DPhillSFOS" + bb +skimFilePrefix]->Fill(dPhi(lep_p4()[i3lv2[1] ],lep_p4()[i3lv2[2] ]),weight);
	    histos["DRllSFOS"   + bb +skimFilePrefix]->Fill(dR(  lep_p4()[i3lv2[1] ],lep_p4()[i3lv2[2] ]),weight);
	  }
	  int l1c(-1), l2c(-1.), l3c(-1); LorentzVector jsum(0,0,0,0);
	  float HT = 0; float MT2ljmin(99999);
	  for(unsigned int j = 0; j<i2p5.size();++j){
	    if(j<=1&&jets_p4()[i2p5[j] ].Pt()<30) continue;
	    if(jets_p4()[i2p5[j] ].Pt()<20) continue;
	    if(l1c<0) l1c = i2p5[j];
	    else if(dR(lep_p4()[i3lv2[0] ],jets_p4()[i2p5[j] ])<dR(lep_p4()[i3lv2[0] ],jets_p4()[l1c ])) l1c = i2p5[j];
	    if(l2c<0) l2c = i2p5[j];
	    else if(dR(lep_p4()[i3lv2[1] ],jets_p4()[i2p5[j] ])<dR(lep_p4()[i3lv2[1] ],jets_p4()[l2c ])) l2c = i2p5[j];
	    if(l3c<0) l3c = i2p5[j];
	    else if(dR(lep_p4()[i3lv2[1] ],jets_p4()[i2p5[j] ])<dR(lep_p4()[i3lv2[1] ],jets_p4()[l3c ])) l3c = i2p5[j];
	    jsum = jsum + jets_p4()[i2p5[j] ];
	    HT += jets_p4()[i2p5[j] ].Pt();
	    for(unsigned int k = j+1; k<i2p5.size();++k){
	      if(k<=1&&jets_p4()[i2p5[k] ].Pt()<30) continue;
	      if(jets_p4()[i2p5[k] ].Pt()<20) continue;
	      float mymt2 = MT2(lep_p4()[i3lv2[0] ]+jets_p4()[i2p5[j] ],lep_p4()[i3lv2[1] ]+jets_p4()[i2p5[k] ],MET,true);    if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	      mymt2 = MT2(lep_p4()[i3lv2[0] ]+jets_p4()[i2p5[k] ],lep_p4()[i3lv2[1] ]+jets_p4()[i2p5[j] ],MET,true);          if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	      mymt2 = MT2(lep_p4()[i3lv2[0] ]+jets_p4()[i2p5[k] ],lep_p4()[i3lv2[2] ]+jets_p4()[i2p5[j] ],MET,true);          if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	      mymt2 = MT2(lep_p4()[i3lv2[0] ]+jets_p4()[i2p5[j] ],lep_p4()[i3lv2[2] ]+jets_p4()[i2p5[k] ],MET,true);          if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	      mymt2 = MT2(lep_p4()[i3lv2[1] ]+jets_p4()[i2p5[k] ],lep_p4()[i3lv2[2] ]+jets_p4()[i2p5[j] ],MET,true);          if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	      mymt2 = MT2(lep_p4()[i3lv2[1] ]+jets_p4()[i2p5[j] ],lep_p4()[i3lv2[2] ]+jets_p4()[i2p5[k] ],MET,true);          if(mymt2<MT2ljmin) MT2ljmin = mymt2;
	    }
	  }
	  if(MT2ljmin>99998) MT2ljmin = 0.;
	  if(l1c>=0&&l2c>=0&&l3c>=0){
	    if( ((lep_p4()[i3lv2[0] ]+jets_p4()[l1c ]).M()<(lep_p4()[i3lv2[1] ]+jets_p4()[l2c ]).M()) && ((lep_p4()[i3lv2[0] ]+jets_p4()[l1c ]).M()<(lep_p4()[i3lv2[2] ]+jets_p4()[l3c ]).M()) ) {
	      histos["Mljclosestmin" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[0] ]+jets_p4()[l1c ]).M(),weight);
	      if((lep_p4()[i3lv2[2] ]+jets_p4()[l3c ]).M()<(lep_p4()[i3lv2[1] ]+jets_p4()[l2c ]).M()) histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[1] ]+jets_p4()[l2c ]).M(),weight);
	      else                                                                                    histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[2] ]+jets_p4()[l3c ]).M(),weight);
	    } else if((lep_p4()[i3lv2[2] ]+jets_p4()[l3c ]).M()<(lep_p4()[i3lv2[1] ]+jets_p4()[l2c ]).M()){
	      histos["Mljclosestmin" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[2] ]+jets_p4()[l3c ]).M(),weight);
	      if((lep_p4()[i3lv2[0] ]+jets_p4()[l1c ]).M()<(lep_p4()[i3lv2[1] ]+jets_p4()[l2c ]).M()) histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[1] ]+jets_p4()[l2c ]).M(),weight);
	      else                                                                                    histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[0] ]+jets_p4()[l1c ]).M(),weight);
	    } else {
	      histos["Mljclosestmin" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[1] ]+jets_p4()[l2c ]).M(),weight);
	      if((lep_p4()[i3lv2[0] ]+jets_p4()[l1c ]).M()<(lep_p4()[i3lv2[2] ]+jets_p4()[l3c ]).M()) histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[2] ]+jets_p4()[l3c ]).M(),weight);
	      else                                                                                    histos["Mljclosestmax" + bb +skimFilePrefix]->Fill((lep_p4()[i3lv2[0] ]+jets_p4()[l1c ]).M(),weight);
	    }
	  }
	  histos["HT" + bb +skimFilePrefix]->Fill(HT,weight);
	  histos["ST" + bb +skimFilePrefix]->Fill(HT+lep_p4()[i3lv2[0] ].Pt()+lep_p4()[i3lv2[1] ].Pt()+lep_p4()[i3lv2[2] ].Pt()+met_pt(),weight);
	  histos["ptjsys"+bb+skimFilePrefix]->Fill(jsum.Pt(),weight);
	  histos["MT2ljmin"+bb+skimFilePrefix]->Fill(MT2ljmin,weight);
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
  
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    //add overflow
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
    //add underflow
    h->second->SetBinContent(1, h->second->GetBinContent(1)+ h->second->GetBinContent(0) );
    h->second->SetBinError(1, sqrt(pow(h->second->GetBinError(1),2)+pow(h->second->GetBinError(0),2) ) );
  }
  
  string filename = "rootfiles/StudyVariousObjects.root";
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
