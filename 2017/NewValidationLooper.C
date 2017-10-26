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
#include "MyFunctions.h"
#include "CMS3_WWW0116.cc" //--> is this loaded by MyFunctions.C ?
#include "/home/users/haweber/CORE/Tools/dorky/dorky.h"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.cc"
#include "/home/users/haweber/CORE/Tools/goodrun.h"
#include "/home/users/haweber/CORE/Tools/goodrun.cc"
//MT2 variants

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  bool btagreweighting = true;
  bool applylepSF = true;
  float muptmin = 20.1;                           float muptmax = 199.9; float muetamin =  0.01; float muetamax = 2.49;
  float elptmin = 10.1; float elptminReco = 25.1; float elptmax = 499.9; float eletamin = -2.49; float eletamax = 2.49;
  TFile *fSF;
  TH2F *hMu, *hX, *hElReco, *hElID;
  int SFloaded = -1;
  if(applylepSF) SFloaded = loadlepSFfile(fSF,hMu,hX,hElReco,hElID,"rootfiles/SF_TnP.root","muSF","","elSF_reco","elSF_ID");

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();

  //split ttX/ttV/ttZ
  histonames.push_back("NB_VR3l_Z_ge2j_ge1b");             hbins.push_back(3);  hlow.push_back( 1); hup.push_back(4);
  histonames.push_back("NB_VR3l_noZ_ge2j_ge1b");           hbins.push_back(3);  hlow.push_back( 1); hup.push_back(4);
  histonames.push_back("NB_VR3l_Z_ge2j_ge1bmed");          hbins.push_back(3);  hlow.push_back( 1); hup.push_back(4);
  histonames.push_back("NB_VR3l_noZ_ge2j_ge1bmed");        hbins.push_back(3);  hlow.push_back( 1); hup.push_back(4);
  histonames.push_back("NB_VRSS_MjjW_ge4j_ge1b");          hbins.push_back(3);  hlow.push_back( 1); hup.push_back(4);
  histonames.push_back("NB_VRSS_MjjW_ge4j_ge1bmed");       hbins.push_back(3);  hlow.push_back( 1); hup.push_back(4);

  histonames.push_back("ZPt_VR3l_Z_ge2j_ge1bmed");         hbins.push_back(15); hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("Ptlll_VR3l_Z_ge2j_ge1bmed");       hbins.push_back(15); hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("Ptlll_VR3l_noZ_ge2j_ge1bmed");     hbins.push_back(15); hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("MET_VR3l_Z_ge2j_ge1bmed");         hbins.push_back(15); hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("MET_VR3l_noZ_ge2j_ge1bmed");       hbins.push_back(15); hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("Ptll_VRSS_MjjW_ge4j_ge1bmed");     hbins.push_back(15); hlow.push_back( 0); hup.push_back(300);
  histonames.push_back("MET_VRSS_MjjW_ge4j_ge1bmed");      hbins.push_back(15); hlow.push_back( 0); hup.push_back(300);

  //WW VBS
  histonames.push_back("MjjL_VRSS_Detajj1p5");             hbins.push_back(10); hlow.push_back( 0); hup.push_back(500);
  histonames.push_back("Detajj_VRSS_MjjL400");             hbins.push_back(10); hlow.push_back( 0); hup.push_back(3.0);

  //Zg
  histonames.push_back("pTlll_VR3l_noZ_METleq50");         hbins.push_back(10); hlow.push_back( 0); hup.push_back(100);
  histonames.push_back("Mlll_VR3l_noZ_METleq50");          hbins.push_back(13); hlow.push_back(20); hup.push_back(150);

  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir,2);
  cout << "Loaded histograms" << endl;

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

      if(firstgoodvertex()!=0)   continue;
      if(nVert()<0)              continue;
      //if(nlep()<2)               continue;

      //weight = 1;
      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("www_2l_mia")!=string::npos)      weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")!=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(weight>100) cout << weight << " " << currentFile->GetTitle() << endl;
      if(isData()) weight = 1.;
      double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();

      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);

      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l);
      float lepSFSS(1.), lepSFerrSS(0.), lepSF3l(1.), lepSFerr3l(0.);
      if(applylepSF){
	lepSFSS = getlepSFWeightandError(lepSFerrSS,iSS,iaSS,hMu,hX,hElReco,hElID, muptmin,muptmax,muetamin,muetamax, muptmin,muptmax,muetamin,muetamax, elptminReco,elptmax,eletamin,eletamax, elptmin,elptmax,eletamin,eletamax);
	lepSF3l = getlepSFWeightandError(lepSFerr3l,i3l,ia3l,hMu,hX,hElReco,hElID, muptmin,muptmax,muetamin,muetamax, muptmin,muptmax,muetamin,muetamax, elptminReco,elptmax,eletamin,eletamax, elptmin,elptmax,eletamin,eletamax);
      }
      
      int nvetoSS = vSS.size();
      int nveto3l = v3l.size();
      int nSS = iSS.size();
      int n3l = i3l.size();
      int nvetoaSS = vaSS.size();
      int nvetoa3l = va3l.size();
      int naSS = iaSS.size();
      int na3l = ia3l.size();
      
      if((n3l+na3l)<2) continue;
      bool passofflineforTrigger = passofflineTriggers(i3l, ia3l);
      if(!passofflineforTrigger) continue;
      
      if(isData()){
	if(!passFilters()) continue;
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) { continue; }
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
	bool passonlineTrigger = passonlineTriggers(i3l, ia3l);//currently applied only to data
	if(!passonlineTrigger) continue;
      }

      string sn = skimFilePrefix;
      string sn2 = skimFilePrefix;
      string temp = process(fname,true ,iSS,iaSS);
      string temp2 = process(fname,false,i3l,ia3l);
      bool isphotonSS = (temp =="photonfakes");
      bool isphoton3l = (temp2=="photonfakes");
      if(sn.find("Other")!=string::npos&&splitVH(fname)){ sn = "WHtoWWW"; sn2 = sn; }
      else if(sn.find("Background")!=string::npos&&splitVH(fname)) continue;
      else if(sn.find("Background")!=string::npos&&!splitVH(fname)){
	sn  = temp;
	sn2 = temp2;
      }
      string extrasn1 = skimFilePrefix;
      string extrasn2 = skimFilePrefix;
      bool passextra1 = false;
      bool passextra2 = false;      
      if(sn=="WW"){
	extrasn1 = "WWVBS";
	sn  = "WWRest";
	sn2 = "WWRest";
	if(fname.find("wpwpjj_ewk")!=string::npos) passextra1 = true;
      }
      if(sn=="ttV"){
	extrasn1 = "ttW";
	extrasn2 = "ttZ";
	if(fname.find("ttw_")!=string::npos) passextra1 = true;
	if(fname.find("ttz_")!=string::npos) passextra2 = true;
      }


      float MTmax = -1;
      if(iSS.size()==2) MTmax = calcMTmax(iSS,MET);
      else if(iSS.size()==1&&iaSS.size()>=1){
	vector<int> temp; temp.push_back(iSS[0]); temp.push_back(iaSS[0]);
	MTmax = calcMTmax(temp,MET);
      }
      float MTmax3l = calcMTmax(i3l,MET,true);

      int VR[10]; bool selects3l[10];
      for(int i = 0; i<10; ++i) { VR[i] = -1; selects3l[i] = false; }
      int nbmed = numJ(20.,2.4,0.8484);
      //fake nj for 3lSR function
      int fakenj = 0;
      if(nj30<2) fakenj = 2;

      //0: VR3l,Z,ge2j,ge1b
      //1: VR3l,noZ,ge2j,ge1b
      checkbothSRCR3l(VR[1],VR[0],i3l,true,fakenj,nb,MET,0,true);
      VR[1] = std::max(VR[1],VR[0]);
      if(VR[1]>0){
	if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[0] ]+lep_p4()[i3l[0] ]).M()-MZ)<10.) { VR[1] = -1; }
	vector<float> MSFOSvec = allMSFOS(i3l);
	int countMZ = 0;
	for(unsigned int i = 0; i<MSFOSvec.size(); ++i) {
	  if(MSFOSvec[i]<20.) { VR[1] = -1; }
	  if(fabs(MSFOSvec[i]-MZ)<15.) ++countMZ;
	}
	if(VR[1]>0){
	  VR[1] = MSFOSvec.size()-countMZ;
	  VR[0] = countMZ;
	}
      }
      selects3l[0] = true; selects3l[1] = true;
      //2: VR3l,Z,ge2j,ge1bmed
      //3: VR3l,noZ,ge2j,ge1b
      checkbothSRCR3l(VR[3],VR[2],i3l,true,fakenj,nbmed,MET,0,true);
      VR[3] = std::max(VR[3],VR[2]);
      if(VR[3]>0){
	if(fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-MZ)<10.) { VR[3] = -1; }
	vector<float> MSFOSvec = allMSFOS(i3l);
	int countMZ = 0;
	for(unsigned int i = 0; i<MSFOSvec.size(); ++i) {
	  if(MSFOSvec[i]<20.) { VR[3] = -1; }
	  if(fabs(MSFOSvec[i]-MZ)<15.) ++countMZ;
	}
	if(VR[3]>0){
	  VR[3] = MSFOSvec.size()-countMZ;
	  VR[2] = countMZ;
	}
      }
      selects3l[2] = true; selects3l[3] = true;
      //4: VRSS,MjjW,ge4j,ge1b
      VR[4] = isSRSS(iSS,vSS,true,MTmax,nj30,nb,Mjj,MjjL,Detajj,MET,0,true);
      if(VR[4]>=0){
	if(nj30<4)            VR[4] = -1;
	if(fabs(Mjj-80.)>20.) VR[4] = -1;
      }
      //5: VRSS,MjjW,ge4j,ge1bmed
      VR[5] = isSRSS(iSS,vSS,true,MTmax,nj30,nbmed,Mjj,MjjL,Detajj,MET,0,true);
      if(VR[5]>=0){
	if(nj30<4)            VR[5] = -1;
	if(fabs(Mjj-80.)>20.) VR[5] = -1;
      }
      //6,7: VRSS,ge2j,0b
      VR[6] = isSRSS(iSS,vSS,true,MTmax,nj30,nb,Mjj,MjjL,Detajj,MET,0,false);
      VR[7] = VR[6];
      if(MjjL  <400.) VR[6] = -1;
      if(Detajj<1.5)  VR[7] = -1;
      //8: VR3l, le1j,0b
      VR[8] = isSR3l(i3l,true,nj,nb);
      if(VR[8]>=0){
	if(met_pt()>50.) VR[8] = -1;
      }
      selects3l[8] = true; 
      for(int i = 0; i<10; ++i) {
	if(!selects3l[i]){
	  if(vetophotonprocess(fname,isphotonSS))    { VR[i] = -1; }
	}
	else if(vetophotonprocess(fname,isphoton3l)){ VR[i] = -1; }
      }
      //cout << __LINE__ << " " << sn << " " << sn2 << " " << passextra1 << " " << passextra2 << " " << extrasn1 << " " << extrasn2 << endl;
      if(VR[0]>=1){
	if(     sn=="ttV"&&passextra1)      histos["NB_VR3l_Z_ge2j_ge1b_"+extrasn1]->Fill(nb,weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["NB_VR3l_Z_ge2j_ge1b_"+extrasn2]->Fill(nb,weight*lepSF3l);
	else                                histos["NB_VR3l_Z_ge2j_ge1b_"+     sn2]->Fill(nb,weight*lepSF3l);
      }
      if(VR[1]>=0){
	if(     sn=="ttV"&&passextra1)      histos["NB_VR3l_noZ_ge2j_ge1b_"+extrasn1]->Fill(nb,weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["NB_VR3l_noZ_ge2j_ge1b_"+extrasn2]->Fill(nb,weight*lepSF3l);
	else                                histos["NB_VR3l_noZ_ge2j_ge1b_"+     sn2]->Fill(nb,weight*lepSF3l);
      }
      if(VR[2]>=1){
	if(     sn=="ttV"&&passextra1)      histos["NB_VR3l_Z_ge2j_ge1bmed_"+extrasn1]->Fill(nb,weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["NB_VR3l_Z_ge2j_ge1bmed_"+extrasn2]->Fill(nb,weight*lepSF3l);
	else                                histos["NB_VR3l_Z_ge2j_ge1bmed_"+     sn2]->Fill(nb,weight*lepSF3l);
	float ptlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	if(     sn=="ttV"&&passextra1)      histos["Ptlll_VR3l_Z_ge2j_ge1bmed_"+extrasn1]->Fill(ptlll,weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["Ptlll_VR3l_Z_ge2j_ge1bmed_"+extrasn2]->Fill(ptlll,weight*lepSF3l);
	else                                histos["Ptlll_VR3l_Z_ge2j_ge1bmed_"+     sn2]->Fill(ptlll,weight*lepSF3l);
	if(     sn=="ttV"&&passextra1)      histos["MET_VR3l_Z_ge2j_ge1bmed_"  +extrasn1]->Fill(met_pt(),weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["MET_VR3l_Z_ge2j_ge1bmed_"  +extrasn2]->Fill(met_pt(),weight*lepSF3l);
	else                                histos["MET_VR3l_Z_ge2j_ge1bmed_"  +     sn2]->Fill(met_pt(),weight*lepSF3l);
	for(unsigned int i1 = 0; i1<3;++i1){
	  for(unsigned int i2 = i1+1; i2<3;++i2){
	    if((lep_pdgId()[i3l[i1] ])!=(-lep_pdgId()[i3l[i2] ])) continue;
	    if(fabs((lep_p4()[i3l[i1] ]+lep_p4()[i3l[i2] ]).M()-MZ)>15.) continue;
	    if(     sn=="ttV"&&passextra1)      histos["ZPt_VR3l_Z_ge2j_ge1bmed_"+extrasn1]->Fill((lep_p4()[i3l[i1] ]+lep_p4()[i3l[i2] ]).Pt(),weight*lepSF3l);
	    else if(sn=="ttV"&&passextra2)      histos["ZPt_VR3l_Z_ge2j_ge1bmed_"+extrasn2]->Fill((lep_p4()[i3l[i1] ]+lep_p4()[i3l[i2] ]).Pt(),weight*lepSF3l);
	    else                                histos["ZPt_VR3l_Z_ge2j_ge1bmed_"+     sn2]->Fill((lep_p4()[i3l[i1] ]+lep_p4()[i3l[i2] ]).Pt(),weight*lepSF3l);
	  }
	}
      }
      if(VR[3]>=0){
	if(     sn=="ttV"&&passextra1)      histos["NB_VR3l_noZ_ge2j_ge1bmed_"+extrasn1]->Fill(nb,weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["NB_VR3l_noZ_ge2j_ge1bmed_"+extrasn2]->Fill(nb,weight*lepSF3l);
	else                                histos["NB_VR3l_noZ_ge2j_ge1bmed_"+     sn2]->Fill(nb,weight*lepSF3l);
	float ptlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	if(     sn=="ttV"&&passextra1)      histos["Ptlll_VR3l_noZ_ge2j_ge1bmed_"+extrasn1]->Fill(ptlll,weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["Ptlll_VR3l_noZ_ge2j_ge1bmed_"+extrasn2]->Fill(ptlll,weight*lepSF3l);
	else                                histos["Ptlll_VR3l_noZ_ge2j_ge1bmed_"+     sn2]->Fill(ptlll,weight*lepSF3l);
	if(     sn=="ttV"&&passextra1)      histos["MET_VR3l_noZ_ge2j_ge1bmed_"  +extrasn1]->Fill(met_pt(),weight*lepSF3l);
	else if(sn=="ttV"&&passextra2)      histos["MET_VR3l_noZ_ge2j_ge1bmed_"  +extrasn2]->Fill(met_pt(),weight*lepSF3l);
	else                                histos["MET_VR3l_noZ_ge2j_ge1bmed_"  +     sn2]->Fill(met_pt(),weight*lepSF3l);
      }
      if(VR[4]>=0){
	if(     sn=="ttV"&&passextra1)      histos["NB_VRSS_MjjW_ge4j_ge1b_"+extrasn1]->Fill(nb,weight*lepSFSS);
	else if(sn=="ttV"&&passextra2)      histos["NB_VRSS_MjjW_ge4j_ge1b_"+extrasn2]->Fill(nb,weight*lepSFSS);
	else                                histos["NB_VRSS_MjjW_ge4j_ge1b_"+      sn]->Fill(nb,weight*lepSFSS);
      }
      if(VR[5]>=0){
	if(     sn=="ttV"&&passextra1)      histos["NB_VRSS_MjjW_ge4j_ge1bmed_"+extrasn1]->Fill(nb,weight*lepSFSS);
	else if(sn=="ttV"&&passextra2)      histos["NB_VRSS_MjjW_ge4j_ge1bmed_"+extrasn2]->Fill(nb,weight*lepSFSS);
	else                                histos["NB_VRSS_MjjW_ge4j_ge1bmed_"+      sn]->Fill(nb,weight*lepSFSS);
	float ptll = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt();
	if(     sn=="ttV"&&passextra1)      histos["Ptll_VRSS_MjjW_ge4j_ge1bmed_"+extrasn1]->Fill(ptll,weight*lepSFSS);
	else if(sn=="ttV"&&passextra2)      histos["Ptll_VRSS_MjjW_ge4j_ge1bmed_"+extrasn2]->Fill(ptll,weight*lepSFSS);
	else                                histos["Ptll_VRSS_MjjW_ge4j_ge1bmed_"+      sn]->Fill(ptll,weight*lepSFSS);
	if(     sn=="ttV"&&passextra1)      histos["MET_VRSS_MjjW_ge4j_ge1bmed_" +extrasn1]->Fill(met_pt(),weight*lepSFSS);
	else if(sn=="ttV"&&passextra2)      histos["MET_VRSS_MjjW_ge4j_ge1bmed_" +extrasn2]->Fill(met_pt(),weight*lepSFSS);
	else                                histos["MET_VRSS_MjjW_ge4j_ge1bmed_" +      sn]->Fill(met_pt(),weight*lepSFSS);
      }
      if(VR[6]>=0){
	if(sn=="WW"&&passextra1) histos["Detajj_VRSS_MjjL400_"+extrasn1]->Fill(Detajj,weight*lepSFSS);
	else                     histos["Detajj_VRSS_MjjL400_"+      sn]->Fill(Detajj,weight*lepSFSS);
      }
      if(VR[7]>=0){
	if(sn=="WW"&&passextra1) histos["MjjL_VRSS_Detajj1p5_"+extrasn1]->Fill(MjjL,weight*lepSFSS);
	else                     histos["MjjL_VRSS_Detajj1p5_"+      sn]->Fill(MjjL,weight*lepSFSS);
      }
      if(VR[8]>=0){
	histos["pTlll_VR3l_noZ_METleq50_"+sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt(),weight*lepSF3l);
	histos["Mlll_VR3l_noZ_METleq50_" +sn2]->Fill((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M() ,weight*lepSF3l);
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
  
  SaveHistosToFile("rootfiles/newValidation.root",histos,true,true);

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
