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
  
  if(skimFilePrefix!="WW"&&skimFilePrefix!="ttV"&&skimFilePrefix!="ZG"&&skimFilePrefix!="WG"&&skimFilePrefix!="WGstar") return 0;
  //if(skimFilePrefix=="Data"||skimFilePrefix=="Background") return 0;
  float wgtsum1(0), wgtsum2(0), wgtsum3(0);
  int counter = 0;
  bool applylepSF = false;

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

  histonames.push_back("SRyield"); hbins.push_back(6); hlow.push_back(    0); hup.push_back(6);
  
  
  cout << "booking histograms" << endl;
  for(unsigned int i = 0; i<histonames.size(); ++i){
    string mapname = histonames[i];
    if(skimFilePrefix=="WW"||skimFilePrefix=="ttV"){
    mapname = histonames[i] + "_trueSS"+skimFilePrefix;
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    mapname = histonames[i] + "_true3l"+skimFilePrefix;
    if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
    if(skimFilePrefix=="WG"||skimFilePrefix=="ZG"||skimFilePrefix=="WGstar"){
      cout << __LINE__ << endl;
      mapname = histonames[i] + "_photonfake"+skimFilePrefix;
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
    if(skimFilePrefix=="WW"){
	mapname = histonames[i] + "_WWRest";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_otherWWSS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_WWVBS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_WWVBSfake";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_WWDPS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_WWDPSFake";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
    if(skimFilePrefix.find("ttV")!=string::npos){
      	mapname = histonames[i] + "_OtherttVWW";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_ttWSS";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_OtherttVWWW";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_OtherttV";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_ttW3l";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_ttWfake";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
	mapname = histonames[i] + "_photonfakettG";
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1D(mapname.c_str(), "", hbins[i], hlow[i], hup[i]);
    }
  }

  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->Sumw2(); h->second->SetDirectory(rootdir);
    cout << h->first << endl;
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

      int nj(0),nb(0),nbmed(0);
      int nj20(0),nj30(0);
      vector<int> i2p5;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	if(fabs(jets_p4()[n].Eta())<2.5) i2p5.push_back(n);
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<5) ++nj;
	if(jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.5) ++nj30;
	if(jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.5) ++nj20;
	if(jets_csv()[n]>0.5426&&jets_p4()[n].Pt()>20&&fabs(jets_p4()[n].Eta())<2.4) ++nb;
	if(jets_csv()[n]>0.8484&&jets_p4()[n].Pt()>30&&fabs(jets_p4()[n].Eta())<2.5) ++nbmed;
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
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool isSS = false; bool is3l = false;
	bool isaSS = false; bool isa3l = false;
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
	if(fabs(lep_p4()[i].Eta())<2.4){
	  if(lep_pass_VVV_cutbased_fo_noiso()[i]&&abs(lep_pdgId()[i])==11&&lep_lostHits()[i]==0/*&&fabs(lep_ip3d()[i])<0.015*/){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.2){
	      if(lep_p4()[i ].Pt()>20&&!is3l)                          { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.2){
	      if(lep_p4()[i ].Pt()>20&&!is3l)                          { ia3l.push_back(i); isa3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2&&!isSS) { iaSS.push_back(i); isaSS = true; }
	    }
	  }
	  if(lep_pass_VVV_cutbased_fo_noiso()[i]&&abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.4&&fabs(lep_ip3d()[i])<0.015){
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

      if(n3l<2) continue;
      //if(nj30<2) continue;
      //if(nb==0) continue;
      if(!passofflineforTrigger) continue;
      if((nSS+naSS)==0&&lep_p4()[i3l[0] ].Pt()<25) continue;
      
      if(isData()){
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) { /*cout << "Event " << tas::run() << ":" << tas::lumi() << ":" << tas::run() << " is duplicated." << endl;*/ continue; }
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = false;
	if(nmu>=2&&(HLT_DoubleMu()) )                             passonlineTrigger = true;
	if(nmu25>=1&&nel>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nmu>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) passonlineTrigger = true;
	//if(nmu>=2&&(HLT_DoubleMu_noiso()) )                             passonlineTrigger = true;
	//if(nmu25>=1&&nel>=1&&(HLT_MuEG_noiso()||HLT_MuEG_noiso_2()))    passonlineTrigger = true;
	//if(nel25>=1&&nel>=2&&(HLT_DoubleEl_noiso()))                    passonlineTrigger = true;
	if(!passonlineTrigger) continue;
      }

      bool trueSS = false;
      bool true3l = false;
      bool photonfakeSS = false;
      bool photonfake3l = false;
      if((iSS.size()+iaSS.size())>=2){
	int l1(-1),l2(-1);
	if(iSS.size()>=2) { l1 = iSS[0]; l2 = iSS[1]; }
	else if(iSS.size()==1&&iaSS.size()>=1) { l1 = iSS[0]; l2 = iaSS[0]; }
	else if(iaSS.size()>=2) { l1 = iaSS[0]; l2 = iaSS[1]; }
	int gentype = gentype_v2(l1,l2,-1);
	if(     gentype==0) trueSS = true;
	else if(gentype==5) photonfakeSS = true;
      }
      if(i3l.size()>=3){
	int l1(-1), l2(-1), l3(-1);
	if(i3l.size()>=3) { l1 = i3l[0]; l2 = i3l[1]; l3 = i3l[2]; }
	int gentype = gentype_v2(l1,l2,l3);
	if(     gentype==0) true3l = true;
	else if(gentype==5) photonfake3l = true;
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
      
      bool ee[50],mm[50],em[50];
      for(int i = 0; i<50; ++i) { ee[i] = false; em[i] = false; mm[i] = false; }
      if(nj30>=2&&nb==0&&passMDetajj&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13&&met_pt()>40.&&MTmax>90.)                                               em[0] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                                        mm[0] = true;
	if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<40){ ee[0] = false; mm[0] = false; }
      }
      if(nj30>=2&&nb==0&&nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&nvetoSS==0&&nSS==2&&((lep_pdgId()[iSS[0] ]*lep_pdgId()[iSS[1] ])>0)){
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==11&&fabs((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()-90.)>10.) ee[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==11)                                                          em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==11&&abs(lep_pdgId()[iSS[1] ])==13)                                                          em[1] = true;
	if(abs(lep_pdgId()[iSS[0] ])==13&&abs(lep_pdgId()[iSS[1] ])==13)                                                          mm[1] = true;
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
      }

      int SFOS[50];
      double pTlll(-1), DPhilllMET(-1);
      bool passMET = false;
      double Mlll(-1), MSFOS(-1), MSFOS2(-1);
      if(i3l.size()>=3){
	DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	Mlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M();
      }
      for(int i = 0; i<50; ++i) { SFOS[i] = -1; }
      //0 : SR, 1: inverted Mll-Z, 2: SR but inverted DPhi,Pt, 3: inverted Mll-Z and inverted DPhi,Pt
      //cout << __LINE__ << " " << n3l << " nb " << nb << " nj " << nj << endl;
      if(nj<2&&nb==0&&nveto3l==0&&n3l==3){
	bool OS01 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[1] ]<0);
	bool OS02 = (lep_pdgId()[i3l[0] ]*lep_pdgId()[i3l[2] ]<0);
	bool OS12 = (lep_pdgId()[i3l[1] ]*lep_pdgId()[i3l[2] ]<0);
	bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	int SFOScounter = 0;
	if(OS01&&SF01) { ++SFOScounter; MSFOS = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M(); }
	if(OS02&&SF02) { ++SFOScounter; if(MSFOS<0) MSFOS = (lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(); else MSFOS2 = (lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M(); }
	if(OS12&&SF12) { ++SFOScounter; if(MSFOS<0) MSFOS = (lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(); else MSFOS2 = (lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M(); }
	bool pass0(false), pass1(false), pass2(false);
	bool pass0X(false), pass1X(false), pass2X(false);
	if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) SFOScounter = -1;
	if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) SFOScounter = -1;
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) SFOScounter = -1;
	if(SFOScounter==0){
	  pass0 = true;
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) pass0 = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) pass0 = false;
	  bool atleastoneSFOSZ = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<15.) { pass0 = false; if(OS01) atleastoneSFOSZ = true; }
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS02) atleastoneSFOSZ = true; }
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<15.) { pass0 = false; if(OS12) atleastoneSFOSZ = true; }
	  //else cout<<OS01<<" "<<SF01<<" ("<<(OS01&&SF01)<<") "<<OS02<<" "<<SF02<<" ("<<(OS02&&SF02)<<") "<<OS12<<" "<<SF12<<" ("<<(OS12&&SF12)<<") "<<endl;	  
	}
	if(SFOScounter==1){
	  pass1 = true;
	  pass1X = true;
	  if(met_pt()>45) passMET=true;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()>55.&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<110.) { pass1 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass1X = false;
	}
	if(SFOScounter==2){
	  pass2 = true;
	  pass2X = true;
	  if(met_pt()>55) passMET=true;
	  bool atleastoneSFOSZ = false;
	  if(OS01&&SF01&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS02&&SF02&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(OS12&&SF12&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-90.)<20.) { pass2 = false; atleastoneSFOSZ = true; }
	  if(!atleastoneSFOSZ) pass2X = false;
	}
	//SR
	if(DPhilllMET>2.7&&pTlll>60&&              pass0) SFOS[0] = 0;
	if(DPhilllMET>2.5&&pTlll>60&&met_pt()>45.&&pass1) SFOS[0] = 1;
	if(DPhilllMET>2.5&&pTlll>60&&met_pt()>55.&&pass2) SFOS[0] = 2;
	//SRpresel
	if(pass0) SFOS[1] = 0;
	if(pass1) SFOS[1] = 1;
	if(pass2) SFOS[1] = 2;
	//CR
	if(DPhilllMET>2.7&&pTlll>60&&              pass0X) SFOS[2] = 0;
	if(DPhilllMET>2.5&&pTlll>60&&met_pt()>45.&&pass1X) SFOS[2] = 1;
	if(DPhilllMET>2.5&&pTlll>60&&met_pt()>55.&&pass2X) SFOS[2] = 2;
	//CRpresel
	if(pass0X) SFOS[3] = 0;
	if(pass1X) SFOS[3] = 1;
	if(pass2X) SFOS[3] = 2;
      }
      for(int i = 0; i<50; ++i){
	if(SFOS[i]<0) continue;
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i] = -1; }
	if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[0] ])==11&&lep_motherIdSS()[i3l[0] ]==(-3)) { SFOS[i] = -1; }
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i] = -1; }
	if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[1] ])==11&&lep_motherIdSS()[i3l[1] ]==(-3)) { SFOS[i] = -1; }
	if(fname.find("wjets")!=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i] = -1; }
	if(fname.find("dy_")  !=string::npos&&abs(lep_pdgId()[i3l[2] ])==11&&lep_motherIdSS()[i3l[2] ]==(-3)) { SFOS[i] = -1; }
	if(abs(lep_charge()[i3l[0] ]+lep_charge()[i3l[1] ]+lep_charge()[i3l[2] ])==3) { SFOS[i] = -1; }
      }

      if(skimFilePrefix=="WW"||skimFilePrefix=="ttV"){
	if(ee[0]&&trueSS     ) histos["SRyield_trueSS"+skimFilePrefix]->Fill(0., weight);
	if(em[0]&&trueSS     ) histos["SRyield_trueSS"+skimFilePrefix]->Fill(1., weight);
	if(mm[0]&&trueSS     ) histos["SRyield_trueSS"+skimFilePrefix]->Fill(2., weight);
	if(SFOS[0]==0&&true3l) histos["SRyield_true3l"+skimFilePrefix]->Fill(3., weight);
	if(SFOS[0]==1&&true3l) histos["SRyield_true3l"+skimFilePrefix]->Fill(4., weight);
	if(SFOS[0]==2&&true3l) histos["SRyield_true3l"+skimFilePrefix]->Fill(5., weight);
      }
      if(ee[0]||em[0]||mm[0]){
	if((skimFilePrefix=="WW")&&(fname.find("wpwpjj")!=string::npos)&&trueSS){
	  if(ee[0]     ) histos["SRyield_WWVBS"]->Fill(0., weight);
	  if(em[0]     ) histos["SRyield_WWVBS"]->Fill(1., weight);
	  if(mm[0]     ) histos["SRyield_WWVBS"]->Fill(2., weight);
	} else if((skimFilePrefix=="WW")&&(fname.find("dbl")!=string::npos)&&trueSS){
	  if(ee[0]     ) histos["SRyield_WWDPS"]->Fill(0., weight);
	  if(em[0]     ) histos["SRyield_WWDPS"]->Fill(1., weight);
	  if(mm[0]     ) histos["SRyield_WWDPS"]->Fill(2., weight);
	} else if((skimFilePrefix=="WW")&&trueSS){
	  if(ee[0]     ) histos["SRyield_otherWWSS"]->Fill(0., weight);
	  if(em[0]     ) histos["SRyield_otherWWSS"]->Fill(1., weight);
	  if(mm[0]     ) histos["SRyield_otherWWSS"]->Fill(2., weight);
	} else if(skimFilePrefix=="WW"){
	  if(ee[0]     ) histos["SRyield_WWRest"]->Fill(0., weight);
	  if(em[0]     ) histos["SRyield_WWRest"]->Fill(1., weight);
	  if(mm[0]     ) histos["SRyield_WWRest"]->Fill(2., weight);
	}
	if((skimFilePrefix=="ttV")&&(fname.find("ttw")!=string::npos)&&trueSS){
	  if(ee[0]     ) histos["SRyield_ttWSS"]->Fill(0., weight);
	  if(em[0]     ) histos["SRyield_ttWSS"]->Fill(1., weight);
	  if(mm[0]     ) histos["SRyield_ttWSS"]->Fill(2., weight);
	} else if((skimFilePrefix=="ttV")&&trueSS){
	  if(ee[0]     ) histos["SRyield_OtherttVWW"]->Fill(0., weight);
	  if(em[0]     ) histos["SRyield_OtherttVWW"]->Fill(1., weight);
	  if(mm[0]     ) histos["SRyield_OtherttVWW"]->Fill(2., weight);
	} else if(skimFilePrefix=="ttV"){
	  if(ee[0]     ) histos["SRyield_OtherttV"]->Fill(0., weight);
	  if(em[0]     ) histos["SRyield_OtherttV"]->Fill(1., weight);
	  if(mm[0]     ) histos["SRyield_OtherttV"]->Fill(2., weight);
	}
      }
      if(SFOS[0]>=0){
	if((skimFilePrefix=="ttV")&&(fname.find("ttw")!=string::npos)&&true3l){
	  if(SFOS[0]==0) histos["SRyield_ttW3l"]->Fill(3., weight);
	  if(SFOS[0]==1) histos["SRyield_ttW3l"]->Fill(4., weight);
	  if(SFOS[0]==2) histos["SRyield_ttW3l"]->Fill(5., weight);
	} else if((skimFilePrefix=="ttV")&&true3l){
	  if(SFOS[0]==0) histos["SRyield_OtherttVWWW"]->Fill(3., weight);
	  if(SFOS[0]==1) histos["SRyield_OtherttVWWW"]->Fill(4., weight);
	  if(SFOS[0]==2) histos["SRyield_OtherttVWWW"]->Fill(5., weight);
	} else if(skimFilePrefix=="ttV"){
	  if(SFOS[0]==0) histos["SRyield_OtherttV"]->Fill(3., weight);
	  if(SFOS[0]==1) histos["SRyield_OtherttV"]->Fill(4., weight);
	  if(SFOS[0]==2) histos["SRyield_OtherttV"]->Fill(5., weight);
	}
      }
      if(skimFilePrefix=="WG"||skimFilePrefix=="ZG"||skimFilePrefix=="WGstar"||((skimFilePrefix=="ttV")&&(fname.find("ttg")!=string::npos))){
	string s = skimFilePrefix;
	if((skimFilePrefix=="ttV")&&(fname.find("ttg")!=string::npos)) s = "ttG";
	//if(skimFilePrefix=="WGstar") cout << "SRyield_photonfake"+s << endl;
	if((ee[0]||em[0]||mm[0])&&photonfakeSS){
	  if(ee[0]     ) histos["SRyield_photonfake"+s]->Fill(0.,weight);
	  if(em[0]     ) histos["SRyield_photonfake"+s]->Fill(1.,weight);
	  if(mm[0]     ) histos["SRyield_photonfake"+s]->Fill(2.,weight);
	}
	if(SFOS[0]>=0&&photonfake3l){
	  if(SFOS[0]==0) histos["SRyield_photonfake"+s]->Fill(3.,weight);
	  if(SFOS[0]==1) histos["SRyield_photonfake"+s]->Fill(4.,weight);
	  if(SFOS[0]==2) histos["SRyield_photonfake"+s]->Fill(5.,weight);	  
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
  
  string filename = "rootfiles/WWTTWSplitted.root";
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
  return 0;
  }
