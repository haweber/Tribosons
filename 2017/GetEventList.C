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

  vector<unsigned int> r,l,e;
  //r.push_back(1); l.push_back(); e.push_back();
  //r.push_back(1); l.push_back(1041); e.push_back(208123);
  //r.push_back(1); l.push_back(2895); e.push_back(578874);

r.push_back(1); l.push_back( 1152); e.push_back( 230244);
r.push_back(1); l.push_back( 6771); e.push_back( 1353835);
r.push_back(1); l.push_back( 2223); e.push_back( 444505);
r.push_back(1); l.push_back( 2269); e.push_back( 453566);
r.push_back(1); l.push_back( 2716); e.push_back( 542959);
r.push_back(1); l.push_back( 2943); e.push_back( 588446);
r.push_back(1); l.push_back( 4525); e.push_back( 904613);
r.push_back(1); l.push_back( 5955); e.push_back( 1190638);
r.push_back(1); l.push_back( 7434); e.push_back( 1486415);
r.push_back(1); l.push_back( 7574); e.push_back( 1514407);
r.push_back(1); l.push_back( 7862); e.push_back( 1571931);
r.push_back(1); l.push_back( 9423); e.push_back( 1884172);
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
  histonames.push_back("YieldsSR");                          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_raw");                      hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
  histonames.push_back("YieldsSR_rawweight");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);//blind the data
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

      bool checkevent = false;
      for(unsigned int n = 0; n<r.size();++n){
	if(tas::evt()==1048333) checkevent = true;
	if(r[n]==tas::run()&&l[n]==tas::lumi()&&e[n]==tas::evt()){
	  checkevent = true;
	  //break;
	}
      }
      if(checkevent) cout << "Check event " << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << endl;
    
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

      int nj(0),nb(0);
      int nj20(0),nj30(0);
      vector<int> i2p5;
      for(unsigned int n = 0; n<jets_csv().size();++n){
	if(checkevent) cout << "jet" << n << " pt " << jets_p4()[n].Pt() << " eta " << (jets_p4()[n].Eta()) << " csv " << jets_csv()[n] << endl;
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
      if(fabs(MjjL)>400.)   passDetajj = false;
      bool passDetajj2 = true;
      if(nj30<2)            passDetajj2 = false;
      if(fabs(Detajj)>1.5)  passDetajj2 = false;
     vector<int> vSS, v3l, v, iSS, i3l;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){

	bool isSS = false; bool is3l = false;
	bool isaSS = false; bool isa3l = false;
	bool islSS = false; bool isl3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015) {
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      //if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>20&&lep_tightCharge()[i]==2) { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      //if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>20&&lep_tightCharge()[i]==2) { i3l.push_back(i); is3l = true; }
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
	if(checkevent) cout << "lep" << i << " pdgID " << (lep_pdgId()[i]) << " pT " << lep_p4()[i ].Pt() << " eta " << (lep_p4()[i].Eta()) << " iso " << lep_relIso03EAv2()[i] << " IP " << (lep_ip3d()[i]) << " ID " << lep_pass_VVV_cutbased_tight_noiso()[i] << "/" << lep_pass_VVV_cutbased_veto_noiso()[i] << " 3q " << lep_tightCharge()[i] << " losthits " << lep_lostHits()[i] <<  " is3l " << is3l << endl;

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

      if(checkevent) cout << "n3l " << n3l << " nveto " << nveto3l << " nj " << nj << " nb " << nb << endl;

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
      if(checkevent) cout << "passofflinetrigger" << endl;

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
	if(checkevent) cout << "passonlinetrigger" << endl;

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
	  
	  // if(sn2=="others")
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

      if(checkevent&&n3l>=3) cout << "pTlll " << pTlll << " DPhilllMET " << DPhilllMET << " met " << met_pt() << " Mlll " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M() << " M01 " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M() << " M02 " << (lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M() << " M12 " << (lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M() << endl;
      if(checkevent) cout << "ee " << ee[0] << " emu " << em[0] << " mm " << mm[0] << " SFOS= " << SFOS[0] << endl;

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
  */
  
  ofstream eventstocheckEElog    ((string("logs/TESTeventListEE_"   +skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheckEMlog    ((string("logs/TESTeventListEM_"   +skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheckMMlog    ((string("logs/TESTeventListMM_"   +skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheck0SFOSlog ((string("logs/TESTeventList0SFOS_"+skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheck1SFOSlog ((string("logs/TESTeventList1SFOS_"+skimFilePrefix+".log")).c_str(), ios::trunc);
  ofstream eventstocheck2SFOSlog ((string("logs/TESTeventList2SFOS_"+skimFilePrefix+".log")).c_str(), ios::trunc);
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
