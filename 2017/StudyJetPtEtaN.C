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

  histonames.push_back("Yields_mNB20");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_mNB25");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_mNB30");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_lNB20");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_lNB25");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_lNB30");                  hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt4040_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt4040_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt4040_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt4030_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt4030_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt4030_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3030_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3030_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3030_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3020_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3020_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3020_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2020_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2020_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2020_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3025_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3025_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt3025_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2525_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2525_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2525_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2520_2p5");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2520_3p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_SSjpt2520_5p0");          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_2p5_0j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_3p0_0j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_5p0_0j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_2p5_0j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_3p0_0j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_5p0_0j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_2p5_1j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_3p0_1j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_5p0_1j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_2p5_1j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_3p0_1j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_5p0_1j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_2p5_2j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_3p0_2j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt30_5p0_2j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_2p5_2j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_3p0_2j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_vetojpt20_5p0_2j");       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Yields_novetojptcut");           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  //ee,em,mm,0SFOS,1SFOS,2SFOS
  //e.g. select all standard cuts, vary NB20/25/50 (loose and medium) --> 6 regions
  //select standard cuts: for SS: 2j selection: pt 40/40, 40/30, 30/30, 30/20, 20/20 --> 5 region
  //select standard cuts: for SS: 2j: eta 2.5/3.0/5.0 (affects veto jets) --> 3 regions
  //select standard cuts: test veto jets:, 0j, <=1j, <=2j, pt>30/20, eta 2.5/3.0/5.0 --> 12 regions | 0j,<=1j,<=2j for 3l region (on top of sel.jets)
  //--> total study: 6 + 5*3 + 12 = 33 possibilities
  //after that: study dijet properties for SS.

  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
	string mapname = histonames[i] + "_"+skimFilePrefix;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
	//cout << mapname << endl;
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


      //cout << __LINE__ << endl;
      int nj(0),nb(0);
      int nj20(0),nj30(0);
      vector<int> i2p5;
      vector<int> i40_2p5, i30_2p5, i20_2p5, i40_3p0, i30_3p0, i20_3p0, i40_5p0, i30_5p0, i20_5p0, i25_2p5, i25_3p0, i25_5p0;
      vector<int> bl20, bl25, bl30, bm20, bm25, bm30;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	if(fabs(jets_p4()[n].Eta())<2.5) i2p5.push_back(n);
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<5.0) i20_5p0.push_back(n);
	if(jets_p4()[n].Pt()>25&&fabs(jets_p4()[n].Eta())<5.0) i25_5p0.push_back(n);
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<5.0) i30_5p0.push_back(n);
	if(jets_p4()[n].Pt()>40&&fabs(jets_p4()[n].Eta())<5.0) i40_5p0.push_back(n);
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<3.0) i20_3p0.push_back(n);
	if(jets_p4()[n].Pt()>25&&fabs(jets_p4()[n].Eta())<3.0) i25_3p0.push_back(n);
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<3.0) i30_3p0.push_back(n);
	if(jets_p4()[n].Pt()>40&&fabs(jets_p4()[n].Eta())<3.0) i40_3p0.push_back(n);
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.5) i20_2p5.push_back(n);
	if(jets_p4()[n].Pt()>25&&fabs(jets_p4()[n].Eta())<2.5) i25_2p5.push_back(n);
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.5) i30_2p5.push_back(n);
	if(jets_p4()[n].Pt()>40&&fabs(jets_p4()[n].Eta())<2.5) i40_2p5.push_back(n);
	if(jets_csv()[n]>0.5426&&jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.4) bl20.push_back(n);
	if(jets_csv()[n]>0.5426&&jets_p4()[n].Pt()>25&&fabs(jets_p4()[n].Eta())<2.4) bl25.push_back(n);
	if(jets_csv()[n]>0.5426&&jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.4) bl30.push_back(n);
	if(jets_csv()[n]>0.8484&&jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.4) bm20.push_back(n);
	if(jets_csv()[n]>0.8484&&jets_p4()[n].Pt()>25&&fabs(jets_p4()[n].Eta())<2.4) bm25.push_back(n);
	if(jets_csv()[n]>0.8484&&jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.4) bm30.push_back(n);
      }
      //cout << __LINE__ << endl;
      float Mjj_2p5 = -1; float Detajj_2p5 = -1;
      float Mjj_3p0 = -1; float Detajj_3p0 = -1;
      float Mjj_5p0 = -1; float Detajj_5p0 = -1;
      if(i20_2p5.size()>1){
	Mjj_2p5 = (jets_p4()[i20_2p5[0] ]+jets_p4()[i20_2p5[1] ]).M();
	Detajj_2p5 = dEta(jets_p4()[i20_2p5[0] ],jets_p4()[i20_2p5[1] ]);
      }
      if(i20_3p0.size()>1){
	Mjj_3p0 = (jets_p4()[i20_3p0[0] ]+jets_p4()[i20_3p0[1] ]).M();
	Detajj_3p0 = dEta(jets_p4()[i20_3p0[0] ],jets_p4()[i20_3p0[1] ]);
      }
      if(i20_5p0.size()>1){
	Mjj_5p0 = (jets_p4()[i20_5p0[0] ]+jets_p4()[i20_5p0[1] ]).M();
	Detajj_5p0 = dEta(jets_p4()[i20_5p0[0] ],jets_p4()[i20_5p0[1] ]);
      }
      bool passMDetajj_2p5 = false;
      bool passMDetajj_3p0 = false;
      bool passMDetajj_5p0 = false;
      if(i20_2p5.size()>=2&&fabs(Detajj_2p5)<1.5&&fabs(Mjj_2p5-85.)<20.) passMDetajj_2p5 = true;
      if(i20_3p0.size()>=2&&fabs(Detajj_3p0)<1.5&&fabs(Mjj_3p0-85.)<20.) passMDetajj_3p0 = true;
      if(i20_5p0.size()>=2&&fabs(Detajj_5p0)<1.5&&fabs(Mjj_5p0-85.)<20.) passMDetajj_5p0 = true;

      //cout << __LINE__ << endl;
     vector<int> i0, i1, i2, i3, v, iSS, i3l;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	if(abs(lep_pdgId()[i])==11&&lep_tightCharge()[i]==2&&lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSS.push_back(i);
	if(abs(lep_pdgId()[i])==13&&lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>30) iSS.push_back(i);
	if(abs(lep_pdgId()[i])==11&&/*lep_tightCharge()[i]==2&&*/lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3l.push_back(i);
	if(abs(lep_pdgId()[i])==13&&lep_relIso03EA()[i]<0.1 && fabs(lep_ip3d()[i])<0.015&&lep_p4()[i ].Pt()>20) i3l.push_back(i);	if(lep_relIso03EA()[i]<=0.1 && fabs(lep_ip3d()[i])<=0.015&&lep_p4()[i].Pt()>30) i0.push_back(i);
	if(lep_relIso03EA()[i]<=0.1 && fabs(lep_ip3d()[i])<=0.015) i1.push_back(i);
	if(lep_relIso03EA()[i]<=0.1 && fabs(lep_ip3d()[i])<=0.015 && lep_tightCharge()[i]==2){
	  i2.push_back(i);
	  if(abs(lep_pdgId()[i])==11&&lep_ptRatio()[i]>0.65) i3.push_back(i);
	  if(abs(lep_pdgId()[i])==13&&lep_ptRatio()[i]>0.5) i3.push_back(i);
	}
	if(lep_relIso03EA()[i]<=0.3) v.push_back(i);
      }
      //cout << __LINE__ << endl;

      bool preselectSS = false;//all cuts that are not ee/mumu/emu specific and unrelated to jets
      bool preselect3l = false;
      if(nisoTrack_mt2()==0/*&&v.size()==2&&*/&&nlep()==2&&iSS.size()==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&lep_p4()[iSS[0] ].Pt()>30&&lep_p4()[iSS[1] ].Pt()>30&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) preselectSS = true;
      if(/*nisoTrack_mt2()==0&&v.size()==3&&*/nlep()==3&&i3l.size()==3&&lep_p4()[i3l[0] ].Pt()>20&&lep_p4()[i3l[1] ].Pt()>20&&lep_p4()[i3l[2] ].Pt()>20&&dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET)>2.5) preselect3l = true;
      bool preselectee(false), preselectemu(false), preselectmumu(false);
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[0] ]==2&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>55.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) preselectee = true;
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&lep_tightCharge()[iSS[1] ]==2&&met_pt()>55.) preselectemu = true;
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&lep_tightCharge()[iSS[0] ]==2&&met_pt()>55.) preselectemu = true;
      if(preselectSS&&abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13) preselectmumu = true;
      bool preselect0SFOS(false), preselect1SFOS(false), preselect2SFOS(false);
      //	cout << __LINE__ << endl;
      if(preselect3l){
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
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-91.)<15.) SFOScounter = -1;
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-91.)<15.) SFOScounter = -1;
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-91.)<15.) SFOScounter = -1;
	}
	if(SFOScounter==1){
	  if(met_pt()<45.) SFOScounter = -1;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>56.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<111.) SFOScounter = -1;
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()>56.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<111.) SFOScounter = -1;
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()>56.&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<111.) SFOScounter = -1;
	}
	if(SFOScounter==2){
	  if(met_pt()<55.) SFOScounter = -1;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-91.)<20.) SFOScounter = -1;
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-91.)<20.) SFOScounter = -1;
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-91.)<20.) SFOScounter = -1;
	}
	if(SFOScounter==0) preselect0SFOS = true;
	if(SFOScounter==1) preselect1SFOS = true;
	if(SFOScounter==2) preselect2SFOS = true;
      }
      //	cout << __LINE__ << endl;

      //bool preselectee(false), preselectemu(false), preselectmumu(false);
      if(preselectSS){//less computational?
	//cout << __LINE__ << endl;
	if(i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_lNB20_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_lNB20_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_lNB20_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=2&&bl25.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_lNB25_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_lNB25_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_lNB25_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=2&&bl30.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_lNB30_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_lNB30_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_lNB30_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=2&&bm20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_mNB20_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_mNB20_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_mNB20_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=2&&bm25.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_mNB25_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_mNB25_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_mNB25_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=2&&bm30.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_mNB30_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_mNB30_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_mNB30_"+sn]->Fill(2.,weight);
	}
	if(i40_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt4040_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt4040_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt4040_2p5_"+sn]->Fill(2.,weight);
	}
	if(i40_2p5.size()>=1&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt4030_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt4030_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt4030_2p5_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt3030_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3030_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3030_2p5_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=1&&i20_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt3020_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3020_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3020_2p5_"+sn]->Fill(2.,weight);
	}
	if(i20_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt2020_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2020_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2020_2p5_"+sn]->Fill(2.,weight);
	}
	if(i40_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt4040_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt4040_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt4040_3p0_"+sn]->Fill(2.,weight);
	}
	if(i40_3p0.size()>=1&&i30_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt4030_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt4030_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt4030_3p0_"+sn]->Fill(2.,weight);
	}
	if(i30_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt3030_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3030_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3030_3p0_"+sn]->Fill(2.,weight);
	}
	if(i30_3p0.size()>=1&&i20_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt3020_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3020_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3020_3p0_"+sn]->Fill(2.,weight);
	}
	if(i20_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt2020_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2020_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2020_3p0_"+sn]->Fill(2.,weight);
	}
	if(i40_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt4040_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt4040_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt4040_5p0_"+sn]->Fill(2.,weight);
	}
	if(i40_5p0.size()>=1&&i30_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt4030_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt4030_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt4030_5p0_"+sn]->Fill(2.,weight);
	}
	if(i30_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt3030_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3030_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3030_5p0_"+sn]->Fill(2.,weight);
	}
	if(i30_5p0.size()>=1&&i20_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt3020_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3020_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3020_5p0_"+sn]->Fill(2.,weight);
	}
	if(i20_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt2020_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2020_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2020_5p0_"+sn]->Fill(2.,weight);
	}
	if(i30_5p0.size()>=1&&i25_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt3025_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3025_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3025_5p0_"+sn]->Fill(2.,weight);
	}
	if(i25_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt2525_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2525_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2525_5p0_"+sn]->Fill(2.,weight);
	}
	if(i25_5p0.size()>=1&&i20_5p0.size()>=2&&bl20.size()==0&&passMDetajj_5p0){
	  if(preselectee)   histos["Yields_SSjpt2520_5p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2520_5p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2520_5p0_"+sn]->Fill(2.,weight);
	}
	if(i30_3p0.size()>=1&&i25_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt3025_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3025_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3025_3p0_"+sn]->Fill(2.,weight);
	}
	if(i25_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt2525_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2525_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2525_3p0_"+sn]->Fill(2.,weight);
	}
	if(i25_3p0.size()>=1&&i20_3p0.size()>=2&&bl20.size()==0&&passMDetajj_3p0){
	  if(preselectee)   histos["Yields_SSjpt2520_3p0_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2520_3p0_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2520_3p0_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=1&&i25_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt3025_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt3025_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt3025_2p5_"+sn]->Fill(2.,weight);
	}
	if(i25_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt2525_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2525_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2525_2p5_"+sn]->Fill(2.,weight);
	}
	if(i25_2p5.size()>=1&&i20_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_SSjpt2520_2p5_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_SSjpt2520_2p5_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_SSjpt2520_2p5_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()<=2&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_2p5_0j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_2p5_0j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_2p5_0j_"+sn]->Fill(2.,weight);
	}
	if(i30_3p0.size()<=2&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_3p0_0j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_3p0_0j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_3p0_0j_"+sn]->Fill(2.,weight);
	}
	if(i30_5p0.size()<=2&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_5p0_0j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_5p0_0j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_5p0_0j_"+sn]->Fill(2.,weight);
	}
	if(i20_2p5.size()<=2&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_2p5_0j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_2p5_0j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_2p5_0j_"+sn]->Fill(2.,weight);
	}
	if(i20_3p0.size()<=2&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_3p0_0j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_3p0_0j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_3p0_0j_"+sn]->Fill(2.,weight);
	}
	if(i20_5p0.size()<=2&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_5p0_0j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_5p0_0j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_5p0_0j_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()<=3&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_2p5_1j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_2p5_1j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_2p5_1j_"+sn]->Fill(2.,weight);
	}
	if(i30_3p0.size()<=3&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_3p0_1j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_3p0_1j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_3p0_1j_"+sn]->Fill(2.,weight);
	}
	if(i30_5p0.size()<=3&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_5p0_1j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_5p0_1j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_5p0_1j_"+sn]->Fill(2.,weight);
	}
	if(i20_2p5.size()<=3&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_2p5_1j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_2p5_1j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_2p5_1j_"+sn]->Fill(2.,weight);
	}
	if(i20_3p0.size()<=3&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_3p0_1j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_3p0_1j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_3p0_1j_"+sn]->Fill(2.,weight);
	}
	if(i20_5p0.size()<=3&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_5p0_1j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_5p0_1j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_5p0_1j_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()<=4&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_2p5_2j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_2p5_2j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_2p5_2j_"+sn]->Fill(2.,weight);
	}
	if(i30_3p0.size()<=4&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_3p0_2j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_3p0_2j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_3p0_2j_"+sn]->Fill(2.,weight);
	}
	if(i30_5p0.size()<=4&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt30_5p0_2j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt30_5p0_2j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt30_5p0_2j_"+sn]->Fill(2.,weight);
	}
	if(i20_2p5.size()<=4&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_2p5_2j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_2p5_2j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_2p5_2j_"+sn]->Fill(2.,weight);
	}
	if(i20_3p0.size()<=4&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_3p0_2j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_3p0_2j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_3p0_2j_"+sn]->Fill(2.,weight);
	}
	if(i20_5p0.size()<=4&&i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_vetojpt20_5p0_2j_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_vetojpt20_5p0_2j_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_vetojpt20_5p0_2j_"+sn]->Fill(2.,weight);
	}
	if(i30_2p5.size()>=2&&bl20.size()==0&&passMDetajj_2p5){
	  if(preselectee)   histos["Yields_novetojptcut_"+sn]->Fill(0.,weight);
	  if(preselectemu)  histos["Yields_novetojptcut_"+sn]->Fill(1.,weight);
	  if(preselectmumu) histos["Yields_novetojptcut_"+sn]->Fill(2.,weight);
	}
	//cout << __LINE__ << endl;
      }//SS done
      //bool preselect0SFOS(false), preselect1SFOS(false), preselect2SFOS(false);
      if(preselect3l){
	//cout << __LINE__ << endl;
	if(i30_2p5.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_lNB20_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_lNB20_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_lNB20_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=1&&bl25.size()==0){
	  if(preselect0SFOS) histos["Yields_lNB25_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_lNB25_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_lNB25_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=1&&bl30.size()==0){
	  if(preselect0SFOS) histos["Yields_lNB30_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_lNB30_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_lNB30_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=1&&bm20.size()==0){
	  if(preselect0SFOS) histos["Yields_mNB20_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_mNB20_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_mNB20_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=1&&bm25.size()==0){
	  if(preselect0SFOS) histos["Yields_mNB25_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_mNB25_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_mNB25_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=1&&bm30.size()==0){
	  if(preselect0SFOS) histos["Yields_mNB30_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_mNB30_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_mNB30_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_SSjpt4040_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt4040_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt4040_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt4040_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt4040_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt4040_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt4040_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt4040_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt4040_5p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt4030_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt4030_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt4030_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt4030_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt4030_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt4030_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt4030_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt4030_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt4030_5p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3030_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3030_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3030_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3030_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3030_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3030_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3030_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3030_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3030_5p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3020_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3020_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3020_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3020_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3020_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3020_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3020_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3020_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3020_5p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2020_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2020_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2020_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2020_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2020_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2020_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2020_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2020_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2020_5p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3025_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3025_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3025_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3025_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3025_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3025_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt3025_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt3025_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt3025_5p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2525_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2525_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2525_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2525_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2525_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2525_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2525_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2525_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2525_5p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2520_2p5_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2520_2p5_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2520_2p5_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2520_3p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2520_3p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2520_3p0_"+sn]->Fill(5.,weight);
	  if(preselect0SFOS) histos["Yields_SSjpt2520_5p0_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_SSjpt2520_5p0_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_SSjpt2520_5p0_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=0&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_2p5_0j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_2p5_0j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_2p5_0j_"+sn]->Fill(5.,weight);
	}
	if(i30_3p0.size()<=0&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_3p0_0j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_3p0_0j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_3p0_0j_"+sn]->Fill(5.,weight);
	}
	if(i30_5p0.size()<=0&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_5p0_0j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_5p0_0j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_5p0_0j_"+sn]->Fill(5.,weight);
	}
	if(i20_2p5.size()<=0&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_2p5_0j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_2p5_0j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_2p5_0j_"+sn]->Fill(5.,weight);
	}
	if(i20_3p0.size()<=0&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_3p0_0j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_3p0_0j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_3p0_0j_"+sn]->Fill(5.,weight);
	}
	if(i20_5p0.size()<=0&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_5p0_0j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_5p0_0j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_5p0_0j_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_2p5_1j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_2p5_1j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_2p5_1j_"+sn]->Fill(5.,weight);
	}
	if(i30_3p0.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_3p0_1j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_3p0_1j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_3p0_1j_"+sn]->Fill(5.,weight);
	}
	if(i30_5p0.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_5p0_1j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_5p0_1j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_5p0_1j_"+sn]->Fill(5.,weight);
	}
	if(i20_2p5.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_2p5_1j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_2p5_1j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_2p5_1j_"+sn]->Fill(5.,weight);
	}
	if(i20_3p0.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_3p0_1j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_3p0_1j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_3p0_1j_"+sn]->Fill(5.,weight);
	}
	if(i20_5p0.size()<=1&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_5p0_1j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_5p0_1j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_5p0_1j_"+sn]->Fill(5.,weight);
	}
	if(i30_2p5.size()<=2&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_2p5_2j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_2p5_2j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_2p5_2j_"+sn]->Fill(5.,weight);
	}
	if(i30_3p0.size()<=2&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_3p0_2j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_3p0_2j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_3p0_2j_"+sn]->Fill(5.,weight);
	}
	if(i30_5p0.size()<=2&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt30_5p0_2j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt30_5p0_2j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt30_5p0_2j_"+sn]->Fill(5.,weight);
	}
	if(i20_2p5.size()<=2&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_2p5_2j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_2p5_2j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_2p5_2j_"+sn]->Fill(5.,weight);
	}
	if(i20_3p0.size()<=2&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_3p0_2j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_3p0_2j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_3p0_2j_"+sn]->Fill(5.,weight);
	}
	if(i20_5p0.size()<=2&&bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_vetojpt20_5p0_2j_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_vetojpt20_5p0_2j_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_vetojpt20_5p0_2j_"+sn]->Fill(5.,weight);
	}
	if(bl20.size()==0){
	  if(preselect0SFOS) histos["Yields_novetojptcut_"+sn]->Fill(3.,weight);
	  if(preselect1SFOS) histos["Yields_novetojptcut_"+sn]->Fill(4.,weight);
	  if(preselect2SFOS) histos["Yields_novetojptcut_"+sn]->Fill(5.,weight);
	}
	//cout << __LINE__ << endl;
      }//3l done

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
  string filename = "rootfiles/StupidStuff_jetptetaN.root";
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
