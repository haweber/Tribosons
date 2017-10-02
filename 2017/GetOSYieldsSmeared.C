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
#include "CMS3_WWW0116.cc"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.h"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.cc"
#include "/home/users/haweber/CORE/Tools/goodrun.h"
#include "/home/users/haweber/CORE/Tools/goodrun.cc"

#include "EnergyScaleCorrection_class.cc"
#include "EnergyScaleCorrection_class.h"

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

  EnergyScaleCorrection_class eScaler("data/Moriond17_23Jan_ele");
  eScaler.doScale     = true;
  eScaler.doSmearings = true;
  TRandom3 *rgen_ = new TRandom3();
  rgen_->SetSeed(123456);
  double ftarget = 1.294;


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

  histonames.push_back("Mee_raw");                hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedonly");        hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_scaledonly");         hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedandscaled");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_raw");              hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedonly");      hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_scaledonly");       hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedandscaled"); hbins.push_back(20); hlow.push_back(80); hup.push_back(100);

  histonames.push_back("Mmumu_raw");                hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_scaledonlyMM");         hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedandscaledMM");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MmumuSS_raw");              hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_scaledonlyMM");       hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedandscaledMM"); hbins.push_back(20); hlow.push_back(80); hup.push_back(100);

  histonames.push_back("Mee_scaledup_smearedonly");      hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_scaleddown_smearedonly");    hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedup_smearedonly");     hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smeareddown_smearedonly");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);

  histonames.push_back("MeeSS_scaledup_smearedonly");    hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_scaleddown_smearedonly");  hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedup_smearedonly");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smeareddown_smearedonly"); hbins.push_back(20); hlow.push_back(80); hup.push_back(100);

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

    string fname = currentFile->GetTitle();
    //if(fname.find("unmerged_mm")!=string::npos) continue;
    //if(fname.find("unmerged_em")!=string::npos) continue;

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

      double ff = (double)nEventsTotal/(double)nEventsChain;
      if(ff>=ftarget) cout << __LINE__ << endl;
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
      vector<int> vSSsc, v3lsc, vsc, iSSsc, i3lsc;
      vector<int> vSSsm, v3lsm, vsm, iSSsm, i3lsm;
      vector<int> vSSss, v3lss, vss, iSSss, i3lss;
      vector<int> vSSscMM, v3lscMM, vscMM, iSSscMM, i3lscMM;
      vector<int> vSSssMM, v3lssMM, vssMM, iSSssMM, i3lssMM;
      vector<LorentzVector> lep_rescaled, lep_rescaledMM, lep_smeared, lep_scaleup, lep_scaledown, lep_smearup, lep_smeardown, lep_newcentral, lep_smsc, lep_newcentralMM;
      if(ff>ftarget) cout << __LINE__ << endl;
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
       if(abs(lep_pdgId()[i])==13){
	 lep_rescaled   .push_back(lep_p4()[i]);
	 lep_rescaledMM .push_back(lep_p4()[i]);
	 lep_smsc       .push_back(lep_p4()[i]);
	 lep_smeared    .push_back(lep_p4()[i]);
	 lep_scaleup    .push_back(lep_p4()[i]);
	 lep_scaledown  .push_back(lep_p4()[i]);
	 lep_smearup    .push_back(lep_p4()[i]);
	 lep_smeardown  .push_back(lep_p4()[i]);
	 lep_newcentral .push_back(lep_p4()[i]);
	 lep_newcentralMM.push_back(lep_p4()[i]);
       } else if(abs(lep_pdgId()[i])==11){
	 bool isEBEle = ( (fabs(lep_etaSC()[i])<=1.479) ? true : false);
	 float scale_corr = 1.;
	 float sigma = 0.;
	 float sigmaData = 0.;
	 float error_scale = 0.;
	 float sigmaup = 0.;
	 float sigmadown = 0.;
	 //Moriond recipe
	 unsigned int gainswitch = 12;
	 if(lep_p4()[i].Et()>150.) gainswitch = 6;
	 if(lep_p4()[i].Et()>300.) gainswitch = 1;
	 /*if(isData())*/    scale_corr  = eScaler.ScaleCorrection(           tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch       );
	 if(!isData())   sigma       = eScaler.getSmearingSigma(          tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch, 0, 0);
	 if(!isData())   sigmaup     = eScaler.getSmearingSigma(          tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  1, 0);//smear up
	 if(!isData())   sigmadown   = eScaler.getSmearingSigma(          tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch, -1, 0);//smear down
	 if(!isData()){
	   float error_stat = eScaler.ScaleCorrectionUncertainty(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  1);//stat
	   float error_syst = eScaler.ScaleCorrectionUncertainty(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  2);//syst
	   float error_gain = eScaler.ScaleCorrectionUncertainty(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  4);//gain
	   error_scale = sqrt(error_stat*error_stat+error_syst*error_syst+error_gain*error_gain);
	 }
	 //ICHEP
	 //if(isData())    scale_corr=eScaler.ScaleCorrection(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et());
	 //if(!isData())   sigma= eScaler.getSmearingSigma(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), 0, 0);
	 //if(!isData())   error_scale=eScaler.ScaleCorrectionUncertainty(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et());
	 //if(!isData())   sigmaup= eScaler.getSmearingSigma(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), 1, 0);//smear up
	 //if(!isData())   sigmadown= eScaler.getSmearingSigma(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), -1, 0);//smear down
	 
	 float smearingfactor = rgen_->Gaus(0.,sigma);
	 if(isData()) smearingfactor = 0.;
	 if(isData()) lep_rescaled  .push_back(lep_p4()[i]*(scale_corr) );//applied to data only
	 else         lep_rescaled  .push_back(lep_p4()[i]*(1.) );//applied to data only
	 lep_rescaledMM.push_back(lep_p4()[i]*(scale_corr) );//applied to data only
	 lep_scaleup   .push_back(lep_p4()[i]*((1.+error_scale)+smearingfactor) );
	 lep_scaledown .push_back(lep_p4()[i]*((1.-error_scale)+smearingfactor) );
	 lep_smeardown .push_back(lep_p4()[i]*(1+smearingfactor*sigmadown/sigma) );
	 lep_smeared   .push_back(lep_p4()[i]*(1+smearingfactor) );//applied to mc only
	 lep_smearup   .push_back(lep_p4()[i]*(1+smearingfactor*sigmaup/sigma  ) );
	 lep_smsc      .push_back(lep_p4()[i]*(scale_corr+smearingfactor) );
	 if(isData()){
	   lep_newcentral .push_back(lep_rescaled[i]);
	   lep_newcentralMM.push_back(lep_rescaled[i]);
	 } else {
	   lep_newcentral .push_back(lep_smeared[i]);
	   lep_newcentralMM.push_back(lep_smsc[i]);
	 }

       }
     }
     if(ff>ftarget) cout << __LINE__ << endl;
   
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
	bool isSS = false; bool is3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_p4()[i ].Pt()>20)                          { i3l.push_back(i); is3l = true; }
	      if(lep_p4()[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSS.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06&&lep_isTriggerSafe_v1()[i]){
	    if(lep_p4()[i ].Pt()>20) { i3l.push_back(i); is3l = true; }
	    if(lep_p4()[i ].Pt()>30) { iSS.push_back(i); isSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_p4()[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4&&lep_isTriggerSafe_v1()[i]) {
	  v.push_back(i);
	  if(!isSS) vSS.push_back(i);
	  if(!is3l) v3l.push_back(i);
	}

	//scaled
	isSS = false; is3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_rescaled[i ].Pt()>20)                          { i3lsc.push_back(i); is3l = true; }
	      if(lep_rescaled[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSsc.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_rescaled[i ].Pt()>20)                          { i3lsc.push_back(i); is3l = true; }
	      if(lep_rescaled[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSsc.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06&&lep_isTriggerSafe_v1()[i]){
	    if(lep_rescaled[i ].Pt()>20) { i3lsc.push_back(i); is3l = true; }
	    if(lep_rescaled[i ].Pt()>30) { iSSsc.push_back(i); isSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_rescaled[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4&&lep_isTriggerSafe_v1()[i]) {
	  vsc.push_back(i);
	  if(!isSS) vSSsc.push_back(i);
	  if(!is3l) v3lsc.push_back(i);
	}

	//smeared
	isSS = false; is3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_smeared[i ].Pt()>20)                          { i3lsm.push_back(i); is3l = true; }
	      if(lep_smeared[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSsm.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_smeared[i ].Pt()>20)                          { i3lsm.push_back(i); is3l = true; }
	      if(lep_smeared[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSsm.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06&&lep_isTriggerSafe_v1()[i]){
	    if(lep_smeared[i ].Pt()>20) { i3lsm.push_back(i); is3l = true; }
	    if(lep_smeared[i ].Pt()>30) { iSSsm.push_back(i); isSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_smeared[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4&&lep_isTriggerSafe_v1()[i]) {
	  vsm.push_back(i);
	  if(!isSS) vSSsm.push_back(i);
	  if(!is3l) v3lsm.push_back(i);
	}
      
	//smeared+scaled
	isSS = false; is3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_newcentral[i ].Pt()>20)                          { i3lss.push_back(i); is3l = true; }
	      if(lep_newcentral[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSss.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_newcentral[i ].Pt()>20)                          { i3lss.push_back(i); is3l = true; }
	      if(lep_newcentral[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSss.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06&&lep_isTriggerSafe_v1()[i]){
	    if(lep_newcentral[i ].Pt()>20) { i3lss.push_back(i); is3l = true; }
	    if(lep_newcentral[i ].Pt()>30) { iSSss.push_back(i); isSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_newcentral[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4&&lep_isTriggerSafe_v1()[i]) {
	  vss.push_back(i);
	  if(!isSS) vSSss.push_back(i);
	  if(!is3l) v3lss.push_back(i);
	}

	if(ff>ftarget) cout << __LINE__ << endl;
	//scaled(incl.MC)
	isSS = false; is3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_rescaledMM[i ].Pt()>20)                          { i3lscMM.push_back(i); is3l = true; }
	      if(lep_rescaledMM[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSscMM.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_rescaledMM[i ].Pt()>20)                          { i3lscMM.push_back(i); is3l = true; }
	      if(lep_rescaledMM[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSscMM.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06&&lep_isTriggerSafe_v1()[i]){
	    if(lep_rescaledMM[i ].Pt()>20) { i3lscMM.push_back(i); is3l = true; }
	    if(lep_rescaledMM[i ].Pt()>30) { iSSscMM.push_back(i); isSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_rescaledMM[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4&&lep_isTriggerSafe_v1()[i]) {
	  vscMM.push_back(i);
	  if(!isSS) vSSscMM.push_back(i);
	  if(!is3l) v3lscMM.push_back(i);
	}
	if(ff>ftarget) cout << __LINE__ << endl;
	//smeared+scaled(incl.MC)
	isSS = false; is3l = false;
	if(lep_pass_VVV_cutbased_tight_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&fabs(lep_ip3d()[i])<0.015&&lep_isTriggerSafe_v1()[i]){
	  if(abs(lep_pdgId()[i])==11){
	    if(fabs(lep_etaSC()[i])<=1.479&&lep_relIso03EAv2()[i]<0.0588){
	      if(lep_newcentralMM[i ].Pt()>20)                          { i3lssMM.push_back(i); is3l = true; }
	      if(lep_newcentralMM[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSssMM.push_back(i); isSS = true; }
	    } else if(fabs(lep_etaSC()[i]) >1.479&&lep_relIso03EAv2()[i]<0.0571){
	      if(lep_newcentralMM[i ].Pt()>20)                          { i3lssMM.push_back(i); is3l = true; }
	      if(lep_newcentralMM[i ].Pt()>30&&lep_tightCharge()[i]==2) { iSSssMM.push_back(i); isSS = true; }
	    }
	  } else if(abs(lep_pdgId()[i])==13&&lep_relIso03EAv2()[i]<0.06&&lep_isTriggerSafe_v1()[i]){
	    if(lep_newcentralMM[i ].Pt()>20) { i3lssMM.push_back(i); is3l = true; }
	    if(lep_newcentralMM[i ].Pt()>30) { iSSssMM.push_back(i); isSS = true; }
	  }
	}
	if(lep_pass_VVV_cutbased_veto_noiso()[i]&&fabs(lep_p4()[i].Eta())<2.4&&lep_newcentralMM[i ].Pt()>10&&lep_relIso03EAv2()[i]<=0.4&&lep_isTriggerSafe_v1()[i]) {
	  vssMM.push_back(i);
	  if(!isSS) vSSssMM.push_back(i);
	  if(!is3l) v3lssMM.push_back(i);
	}
	if(ff>ftarget) cout << __LINE__ << endl;
      }
      
      if(ff>ftarget) cout << __LINE__ << endl;

      int nvetoSS = vSS.size();
      int nveto3l = v3l.size();
      int nSS = iSS.size();
      int n3l = i3l.size();

      int nvetoSSss = vSSss.size();
      int nveto3lss = v3lss.size();
      int nSSss = iSSss.size();
      int nvetoSSsm = vSSsm.size();
      int nveto3lsm = v3lsm.size();
      int nSSsm = iSSsm.size();
      int nvetoSSsc = vSSsc.size();
      int nveto3lsc = v3lsc.size();
      int nSSsc = iSSsc.size();

      int nvetoSSssMM = vSSssMM.size();
      int nveto3lssMM = v3lssMM.size();
      int nSSssMM = iSSssMM.size();
      int nvetoSSscMM = vSSscMM.size();
      int nveto3lscMM = v3lscMM.size();
      int nSSscMM = iSSscMM.size();
      
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
      if((nSS)>=2)         passofflineforTrigger = true;

      if(n3l<2&&(nSS)<2) continue;
      //if(nj30<2) continue;
      //if(nb!=0) continue;
      if(!passofflineforTrigger) continue;
      if((nSS)==0&&lep_p4()[i3l[0] ].Pt()<25) continue;
      
      if(isData()){
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) continue;
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = false;
	if(nmu>=2&&(HLT_DoubleMu()) )                             passonlineTrigger = true;
	if(nmu25>=1&&nel>=1&&HLT_MuEG())                          passonlineTrigger = true;
	if(nel25>=1&&nel>=2&&(HLT_DoubleEl()||HLT_DoubleEl_DZ())) passonlineTrigger = true;
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
	if((iSS.size())>=2||iSSsm.size()>=2||iSSsc.size()>=2||iSSss.size()>=2||iSSscMM.size()>=2||iSSssMM.size()>=2){
	  int l1(-1),l2(-1);
	  if(iSS.size()>=2) { l1 = iSS[0]; l2 = iSS[1]; }
	  else if(iSSss.size()>=2) { l1 = iSSss[0]; l2 = iSSss[1]; }
	  else if(iSSsm.size()>=2) { l1 = iSSsm[0]; l2 = iSSsm[1]; }
	  else if(iSSsc.size()>=2) { l1 = iSSsc[0]; l2 = iSSsc[1]; }
	  else if(iSSscMM.size()>=2) { l1 = iSSscMM[0]; l2 = iSSscMM[1]; }
	  else if(iSSssMM.size()>=2) { l1 = iSSssMM[0]; l2 = iSSssMM[1]; }
	  //else if(iSS.size()==1&&iaSS.size()>=1) { l1 = iSS[0]; l2 = iaSS[0]; }
	  //else if(iaSS.size()>=2) { l1 = iaSS[0]; l2 = iaSS[1]; }
	  int gentype = gentype_v2(l1,l2,-1);
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
	}
      }//bg

      if(ff>ftarget) cout << __LINE__ << endl;
      /*
      if(nSS>=1&&fname.find("wjets") !=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { continue; }
      if(nSS>=2&&fname.find("wjets") !=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { continue; }
      if(nSS>=1&&fname.find("dy_")   !=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { continue; }
      if(nSS>=2&&fname.find("dy_")   !=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { continue; }
      if(nSS>=1&&fname.find("ttbar_")!=string::npos&&abs(lep_pdgId()[iSS[0] ])==11&&lep_motherIdSS()[iSS[0] ]==(-3)) { continue; }
      if(nSS>=2&&fname.find("ttbar_")!=string::npos&&abs(lep_pdgId()[iSS[1] ])==11&&lep_motherIdSS()[iSS[1] ]==(-3)) { continue; }
      */
      if(nSS>=2&&fname.find("wjets") !=string::npos&&gentype_v2(iSS[0],iSS[1],-1)==5) continue;
      if(nSS>=2&&fname.find("dy_"  ) !=string::npos&&gentype_v2(iSS[0],iSS[1],-1)==5) continue;
      if(nSS>=2&&fname.find("ttbar_")!=string::npos&&gentype_v2(iSS[0],iSS[1],-1)==5) continue;
      if(ff>ftarget) cout << __LINE__ << endl;
	    
      float Mee = -1; float MeeSS = -1;
      if(nSS==2&&nvetoSS==0&&abs(lep_pdgId()[iSS[0] ])==11&&(lep_pdgId()[iSS[0] ]==(-1*lep_pdgId()[iSS[1] ]))) Mee   = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();
      if(nSS==2&&nvetoSS==0&&abs(lep_pdgId()[iSS[0] ])==11&&(lep_pdgId()[iSS[0] ]==(   lep_pdgId()[iSS[1] ]))) MeeSS = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();

      if(ff>ftarget) cout << __LINE__ << endl;

      float Meesc = -1; float MeescSS = -1;
      if(nSSsc==2&&nvetoSSsc==0&&abs(lep_pdgId()[iSSsc[0] ])==11&&(lep_pdgId()[iSSsc[0] ]==(-1*lep_pdgId()[iSSsc[1] ]))) Meesc   = (lep_rescaled[iSSsc[0] ]+lep_rescaled[iSSsc[1] ]).M();
      if(nSSsc==2&&nvetoSSsc==0&&abs(lep_pdgId()[iSSsc[0] ])==11&&(lep_pdgId()[iSSsc[0] ]==(   lep_pdgId()[iSSsc[1] ]))) MeescSS = (lep_rescaled[iSSsc[0] ]+lep_rescaled[iSSsc[1] ]).M();

      if(ff>ftarget) cout << __LINE__ << endl;

      float Meesm = -1; float MeesmSS = -1;
      float MeesmUpSm = -1; float MeesmSSUpSm = -1;
      float MeesmDnSm = -1; float MeesmSSDnSm = -1;
      float MeesmUpSc = -1; float MeesmSSUpSc = -1;
      float MeesmDnSc = -1; float MeesmSSDnSc = -1;
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(-1*lep_pdgId()[iSSsm[1] ]))) Meesm   = (lep_smeared[iSSsm[0] ]+lep_smeared[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(   lep_pdgId()[iSSsm[1] ]))) MeesmSS = (lep_smeared[iSSsm[0] ]+lep_smeared[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(-1*lep_pdgId()[iSSsm[1] ]))) MeesmUpSm   = (lep_smearup[iSSsm[0] ]+lep_smearup[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(   lep_pdgId()[iSSsm[1] ]))) MeesmSSUpSm = (lep_smearup[iSSsm[0] ]+lep_smearup[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(-1*lep_pdgId()[iSSsm[1] ]))) MeesmDnSm   = (lep_smeardown[iSSsm[0] ]+lep_smeardown[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(   lep_pdgId()[iSSsm[1] ]))) MeesmSSDnSm = (lep_smeardown[iSSsm[0] ]+lep_smeardown[iSSsm[1] ]).M();
      
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(-1*lep_pdgId()[iSSsm[1] ]))) MeesmUpSc   = (lep_scaleup[iSSsm[0] ]+lep_scaleup[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(   lep_pdgId()[iSSsm[1] ]))) MeesmSSUpSc = (lep_scaleup[iSSsm[0] ]+lep_scaleup[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(-1*lep_pdgId()[iSSsm[1] ]))) MeesmDnSc   = (lep_scaledown[iSSsm[0] ]+lep_scaledown[iSSsm[1] ]).M();
      if(nSSsm==2&&nvetoSSsm==0&&abs(lep_pdgId()[iSSsm[0] ])==11&&(lep_pdgId()[iSSsm[0] ]==(   lep_pdgId()[iSSsm[1] ]))) MeesmSSDnSc = (lep_scaledown[iSSsm[0] ]+lep_scaledown[iSSsm[1] ]).M();

      if(ff>ftarget) cout << __LINE__ << endl;

      float Meess = -1; float MeessSS = -1;
      if(nSSss==2&&nvetoSSss==0&&abs(lep_pdgId()[iSSss[0] ])==11&&(lep_pdgId()[iSSss[0] ]==(-1*lep_pdgId()[iSSss[1] ]))) Meess   = (lep_newcentral[iSSss[0] ]+lep_newcentral[iSSss[1] ]).M();
      if(nSSss==2&&nvetoSSss==0&&abs(lep_pdgId()[iSSss[0] ])==11&&(lep_pdgId()[iSSss[0] ]==(   lep_pdgId()[iSSss[1] ]))) MeessSS = (lep_newcentral[iSSss[0] ]+lep_newcentral[iSSss[1] ]).M();

      if(ff>ftarget) cout << __LINE__ << endl;

      
      float MeescMM = -1; float MeescSSMM = -1;
      if(nSSscMM==2&&nvetoSSscMM==0&&abs(lep_pdgId()[iSSscMM[0] ])==11&&(lep_pdgId()[iSSscMM[0] ]==(-1*lep_pdgId()[iSSscMM[1] ]))) MeescMM   = (lep_rescaled[iSSscMM[0] ]+lep_rescaled[iSSscMM[1] ]).M();
      if(nSSscMM==2&&nvetoSSscMM==0&&abs(lep_pdgId()[iSSscMM[0] ])==11&&(lep_pdgId()[iSSscMM[0] ]==(   lep_pdgId()[iSSscMM[1] ]))) MeescSSMM = (lep_rescaled[iSSscMM[0] ]+lep_rescaled[iSSscMM[1] ]).M();

      float MeessMM = -1; float MeessSSMM = -1;
      if(nSSssMM==2&&nvetoSSssMM==0&&abs(lep_pdgId()[iSSssMM[0] ])==11&&(lep_pdgId()[iSSssMM[0] ]==(-1*lep_pdgId()[iSSssMM[1] ]))) MeessMM   = (lep_newcentral[iSSssMM[0] ]+lep_newcentral[iSSssMM[1] ]).M();
      if(nSSssMM==2&&nvetoSSssMM==0&&abs(lep_pdgId()[iSSssMM[0] ])==11&&(lep_pdgId()[iSSssMM[0] ]==(   lep_pdgId()[iSSssMM[1] ]))) MeessSSMM = (lep_newcentral[iSSssMM[0] ]+lep_newcentral[iSSssMM[1] ]).M();
	    
      float Mmm = -1; float MmmSS = -1;
      if(nSS==2&&nvetoSS==0&&abs(lep_pdgId()[iSS[0] ])==13&&(lep_pdgId()[iSS[0] ]==(-1*lep_pdgId()[iSS[1] ]))) Mmm   = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();
      if(nSS==2&&nvetoSS==0&&abs(lep_pdgId()[iSS[0] ])==13&&(lep_pdgId()[iSS[0] ]==(   lep_pdgId()[iSS[1] ]))) MmmSS = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();

      if(ff>ftarget) cout << __LINE__ << " " << sn << endl;
      if(Mmm    >=80&&Mmm    <=100) histos["Mmumu_raw_"               +sn]->Fill(Mmm      ,weight);
      if(MmmSS  >=80&&MmmSS  <=100) histos["MmumuSS_raw_"             +sn]->Fill(MmmSS    ,weight);
      if(Mee    >=80&&Mee    <=100) histos["Mee_raw_"                 +sn]->Fill(Mee      ,weight);
      if(MeeSS  >=80&&MeeSS  <=100) histos["MeeSS_raw_"               +sn]->Fill(MeeSS    ,weight);
      if(Meesm  >=80&&Meesm  <=100) histos["Mee_smearedonly_"         +sn]->Fill(Meesm    ,weight);
      if(MeesmSS>=80&&MeesmSS<=100) histos["MeeSS_smearedonly_"       +sn]->Fill(MeesmSS  ,weight);
      if(Meesc  >=80&&Meesc  <=100) histos["Mee_scaledonly_"          +sn]->Fill(Meesc    ,weight);
      if(MeescSS>=80&&MeescSS<=100) histos["MeeSS_scaledonly_"        +sn]->Fill(MeescSS  ,weight);
      if(Meess  >=80&&Meess  <=100) histos["Mee_smearedandscaled_"    +sn]->Fill(Meess    ,weight);
      if(MeessSS>=80&&MeessSS<=100) histos["MeeSS_smearedandscaled_"  +sn]->Fill(MeessSS  ,weight);
      if(MeescMM  >=80&&MeescMM  <=100) histos["Mee_scaledonlyMM_"          +sn]->Fill(Meesc    ,weight);
      if(MeescSSMM>=80&&MeescSSMM<=100) histos["MeeSS_scaledonlyMM_"        +sn]->Fill(MeescSS  ,weight);
      if(MeessMM  >=80&&MeessMM  <=100) histos["Mee_smearedandscaledMM_"    +sn]->Fill(Meess    ,weight);
      if(MeessSSMM>=80&&MeessSSMM<=100) histos["MeeSS_smearedandscaledMM_"  +sn]->Fill(MeessSS  ,weight);

      if(MeesmUpSc  >=80&&MeesmUpSc  <=100) histos["Mee_scaledup_smearedonly_"         +sn]->Fill(MeesmUpSc    ,weight);
      if(MeesmSSUpSc>=80&&MeesmSSUpSc<=100) histos["MeeSS_scaledup_smearedonly_"       +sn]->Fill(MeesmSSUpSc  ,weight);
      if(MeesmDnSc  >=80&&MeesmDnSc  <=100) histos["Mee_scaleddown_smearedonly_"       +sn]->Fill(MeesmDnSc    ,weight);
      if(MeesmSSDnSc>=80&&MeesmSSDnSc<=100) histos["MeeSS_scaleddown_smearedonly_"     +sn]->Fill(MeesmSSDnSc  ,weight);
      if(MeesmUpSm  >=80&&MeesmUpSm  <=100) histos["Mee_smearedup_smearedonly_"        +sn]->Fill(MeesmUpSm    ,weight);
      if(MeesmSSUpSm>=80&&MeesmSSUpSm<=100) histos["MeeSS_smearedup_smearedonly_"      +sn]->Fill(MeesmSSUpSm  ,weight);
      if(MeesmDnSm  >=80&&MeesmDnSm  <=100) histos["Mee_smeareddown_smearedonly_"      +sn]->Fill(MeesmDnSm    ,weight);
      if(MeesmSSDnSm>=80&&MeesmSSDnSm<=100) histos["MeeSS_smeareddown_smearedonly_"    +sn]->Fill(MeesmSSDnSm  ,weight);
      if(ff>ftarget) cout << __LINE__ << " " << sn << endl;

 
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
  
  string filename = "rootfiles/OSChecks_smearedMoriond_gain12.root";
  TFile *f = new TFile(filename.c_str(),"UPDATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=histos.begin(); h!=histos.end();++h){
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
