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
#include "CMS3_WWW0116.cc"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.h"
#include "/home/users/haweber/CORE/Tools/dorky/dorky.cc"
#include "/home/users/haweber/CORE/Tools/goodrun.h"
#include "/home/users/haweber/CORE/Tools/goodrun.cc"
//MT2 variants

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  bool blindSR = true;
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

  //YieldsCR_SSany_dropMjj.pdf in note (check test version)
  //Mjj_CRlike_SSany_allSS.pdf
  //MSFOS_CRlike_allSFOS.pdf
  //Mll_invertMETdPhiPt_allSFOS.pdf
  //Mll_inverteitherMETdPhiPt_allSFOS.pdf

  //SR
  histonames.push_back("YieldsSR");                                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_raw");                               hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_rawweight");                         hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR");                                   hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  //SS JEC uncertainty
  histonames.push_back("YieldsSR_jesup");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_jesdn");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj_jesup");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_dropMjj_jesdn");                     hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_jesup");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_jesdn");                             hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_cutonMjj");                          hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_cutonMjj_jesup");                    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_cutonMjj_jesdn");                    hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("Mjj_CRlike_allSS_jesup");                     hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_allSS_jesdn");                     hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  histonames.push_back("Mjj_CRlike_allSS");                           hbins.push_back(12);hlow.push_back(20);hup.push_back(260);
  //MSFOS/lep SF uncertainty
  histonames.push_back("YieldsCR_lepSFup");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsCR_lepSFdn");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lepSFup");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_lepSFdn");                           hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("MSFOS_CRlike_allSS");                         hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_CRlike_allSS_lepSFup");                 hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_CRlike_allSS_lepSFdn");                 hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_CRlike_allSS");                      hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_CRlike_allSS_lepSFup");              hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_CRlike_allSS_lepSFdn");              hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  //validation SS
  histonames.push_back("YieldsSR_Mjjsideband");                       hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_Mjjsideband_lowMET");                hbins.push_back(6); hlow.push_back(0); hup.push_back(6);
  histonames.push_back("YieldsSR_Mjjsideband_lowMTmax");              hbins.push_back(6); hlow.push_back(0); hup.push_back(6);

  //validation 3l
  histonames.push_back("MSFOS_all3l");                                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_lepSFup");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_lepSFdn");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_inverteitherMETdPhiPt");          hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_all3l_invertMETdPhiPt");                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l");                             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_lepSFup");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_lepSFdn");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_inverteitherMETdPhiPt");       hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_all3l_invertMETdPhiPt");             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_1SFOS");                                hbins.push_back(12);hlow.push_back(33);hup.push_back(165);
  histonames.push_back("MSFOS_1SFOS_lepSFup");                        hbins.push_back(12);hlow.push_back(33);hup.push_back(165);
  histonames.push_back("MSFOS_1SFOS_lepSFdn");                        hbins.push_back(12);hlow.push_back(33);hup.push_back(165);
  histonames.push_back("MSFOS_1SFOS_inverteitherMETdPhiPt");          hbins.push_back(12);hlow.push_back(33);hup.push_back(165);
  histonames.push_back("MSFOS_1SFOS_invertMETdPhiPt");                hbins.push_back(12);hlow.push_back(33);hup.push_back(165);
  histonames.push_back("MSFOS_2SFOS");                                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_lepSFup");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_lepSFdn");                        hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_inverteitherMETdPhiPt");          hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOS_2SFOS_invertMETdPhiPt");                hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS");                             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_lepSFup");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_lepSFdn");                     hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_inverteitherMETdPhiPt");       hbins.push_back(12);hlow.push_back(30);hup.push_back(150);
  histonames.push_back("MSFOSall_2SFOS_invertMETdPhiPt");             hbins.push_back(12);hlow.push_back(30);hup.push_back(150);

  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir);
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
      double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      
      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      //LorentzVector METx;   METx  .SetPxPyPzE(met_T1CHS_miniAOD_CORE_pt()   *TMath::Cos(met_T1CHS_miniAOD_CORE_phi()   ),met_T1CHS_miniAOD_CORE_pt()   *TMath::Sin(met_T1CHS_miniAOD_CORE_phi()   ),0,met_T1CHS_miniAOD_CORE_pt()   );
      LorentzVector MET_up; MET_up.SetPxPyPzE(met_T1CHS_miniAOD_CORE_up_pt()*TMath::Cos(met_T1CHS_miniAOD_CORE_up_phi()),met_T1CHS_miniAOD_CORE_up_pt()*TMath::Sin(met_T1CHS_miniAOD_CORE_up_phi()),0,met_T1CHS_miniAOD_CORE_up_pt());
      LorentzVector MET_dn; MET_dn.SetPxPyPzE(met_T1CHS_miniAOD_CORE_dn_pt()*TMath::Cos(met_T1CHS_miniAOD_CORE_dn_phi()),met_T1CHS_miniAOD_CORE_dn_pt()*TMath::Sin(met_T1CHS_miniAOD_CORE_dn_phi()),0,met_T1CHS_miniAOD_CORE_dn_pt());
      //if(fabs(MET.Pt()-METx.Pt())>0.01*MET.Pt()) cout << "Damn " << MET.Pt() << " " << METx.Pt() << endl;

     int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);

      int nj_up(0), nb_up(0), nj30_up(0);
      getalljetnumbers(nj_up,nj30_up,nb_up,1);
      float Mjj_up = -1;
      float MjjL_up = -1; float Detajj_up = -1;
      getMjjAndDeta(Mjj_up,MjjL_up,Detajj_up,1);

      int nj_dn(0), nb_dn(0), nj30_dn(0);
      getalljetnumbers(nj_dn,nj30_dn,nb_dn,-1);
      float Mjj_dn = -1;
      float MjjL_dn = -1; float Detajj_dn = -1;
      getMjjAndDeta(Mjj_dn,MjjL_dn,Detajj_dn,-1);

      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l);
      float lepSFSS(1.), lepSFerrSS(1.), lepSF3l(0.), lepSFerr3l(0.);
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

      if((n3l)<2) continue;
      bool passofflineforTrigger = passofflineTriggers(i3l, ia3l);
      if(!passofflineforTrigger) continue;
      
      if(isData()){
	if(!passFilters()) continue;
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        ) { /*cout << "Event " << tas::run() << ":" << tas::lumi() << ":" << tas::run() << " is duplicated." << endl;*/ continue; }
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
      
      float MTmax      = calcMTmax(iSS,MET);
      float MTmax3l    = calcMTmax(i3l,MET);
      float MTmax_up   = calcMTmax(iSS,MET_up);
      float MTmax3l_up = calcMTmax(i3l,MET_up);
      float MTmax_dn   = calcMTmax(iSS,MET_dn);
      float MTmax3l_dn = calcMTmax(i3l,MET_dn);

      int SRSS[15]; bool selects3l[15];
      int SR3l[15];
      for(int i = 0; i<15; ++i) { SRSS[i] = -1; SR3l[i] = -1; selects3l[i] = false; }

      //0: SR
      SRSS[ 0] = isSRSS(iSS,      vSS,false,MTmax,  nj30,nb,Mjj,MjjL,Detajj,MET,0);
      //1: SR preselect
      SRSS[ 1] = isSRSS(iSS,      vSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj,MET,0);
      //2: CR
      SRSS[ 2] = isCRSS(iSS,i3l,  v3l,false,MTmax3l,nj30,nb,Mjj,MjjL,Detajj,MET,0,false);
      //3: CR preselect
      SRSS[ 3] = isCRSS(iSS,i3l,  v3l,true ,MTmax3l,nj30,nb,Mjj,MjjL,Detajj,MET,0,false);
      //4: CR - noZ
      SRSS[ 4] = isCRSS(iSS,i3l,  v3l,false,MTmax3l,nj30,nb,Mjj,MjjL,Detajj,MET,0,true);
      //5: CR - noZ preselect
      SRSS[ 5] = isCRSS(iSS,i3l,  v3l,true ,MTmax3l,nj30,nb,Mjj,MjjL,Detajj,MET,0,true);
      //JEC UP
      //6: SR
      SRSS[ 6] = isSRSS(iSS,      vSS,false,MTmax_up,  nj30_up,nb_up,Mjj_up,MjjL_up,Detajj_up,MET_up,1);
      //7: SR preselect
      SRSS[ 7] = isSRSS(iSS,      vSS,true ,MTmax_up,  nj30_up,nb_up,Mjj_up,MjjL_up,Detajj_up,MET_up,1);
      //8: CR
      SRSS[ 8] = isCRSS(iSS,i3l,  v3l,false,MTmax3l_up,nj30_up,nb_up,Mjj_up,MjjL_up,Detajj_up,MET_up,1,false);
      //9: CR preselect
      SRSS[ 9] = isCRSS(iSS,i3l,  v3l,true ,MTmax3l_up,nj30_up,nb_up,Mjj_up,MjjL_up,Detajj_up,MET_up,1,false);
      //JEC DOWN
      //10: SR
      SRSS[10]= isSRSS(iSS,      vSS,false,MTmax_dn,  nj30_dn,nb_dn,Mjj_dn,MjjL_dn,Detajj_dn,MET_dn,-1);
      //11: SR preselect
      SRSS[11]= isSRSS(iSS,      vSS,true ,MTmax_dn,  nj30_dn,nb_dn,Mjj_dn,MjjL_dn,Detajj_dn,MET_dn,-1);
      //12: CR
      SRSS[12]= isCRSS(iSS,i3l,  v3l,false,MTmax3l_dn,nj30_dn,nb_dn,Mjj_dn,MjjL_dn,Detajj_dn,MET_dn,-1,false);
      //13: CR preselect
      SRSS[13]= isCRSS(iSS,i3l,  v3l,true ,MTmax3l_dn,nj30_dn,nb_dn,Mjj_dn,MjjL_dn,Detajj_dn,MET_dn,-1,false);
      selects3l[2] = true; selects3l[3] = true;
      selects3l[4] = true; selects3l[5] = true;
      selects3l[8] = true; selects3l[9] = true;
      selects3l[12]= true; selects3l[13]= true;
      //3l
      //0: SR, 2: CR
      checkbothSRCR3l(SR3l[ 0],SR3l[ 2],i3l,false,nj,nb,MET);
      //1: SR preselect, 3: CR preselect
      checkbothSRCR3l(SR3l[ 1],SR3l[ 3],i3l,true ,nj,nb,MET);
      //6: SR, 8: CR - JESUP
      checkbothSRCR3l(SR3l[ 6],SR3l[ 8],i3l,false,nj_up,nb_up,MET_up,1);
      //7: SR preselect, 9: CR preselect - JESUP
      checkbothSRCR3l(SR3l[ 7],SR3l[ 9],i3l,true ,nj_up,nb_up,MET_up,1);
      //10: SR, 12: CR - JESDOWN
      checkbothSRCR3l(SR3l[10],SR3l[12],i3l,false,nj_dn,nb_dn,MET_dn,-1);
      //11: SR preselect, 13: CR preselect - JESDOWN
      checkbothSRCR3l(SR3l[11],SR3l[13],i3l,true ,nj_dn,nb_dn,MET_dn,-1);

      for(int i = 0; i<6; ++i) {
	if(!selects3l[i]){
	  if(vetophotonprocess(fname,isphotonSS))    { SRSS[i] = -1; }
	}
	else if(vetophotonprocess(fname,isphoton3l)){ SRSS[i] = -1; }
	if(vetophotonprocess(fname,isphoton3l))     { SR3l[i] = -1; }
      }
	     
      float Mlll = -1;
      float pTlll = -1;
      float DPhilllMET = -1;
      if(SR3l[1]>=0||SR3l[5]>=5){
	Mlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M();
	pTlll = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	DPhilllMET = dPhi((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]),MET);
      }
      vector <float> MSFOSvec;
      if(SRSS[5]>=0||SR3l[3]>=0||SR3l[1]>=0) MSFOSvec = allMSFOS(i3l);
      
      //cout << sn << " " << sn2 << endl;
      if(!(blindSR&&isData())){
	fillSRhisto(histos, "YieldsSR",           sn, sn2, SRSS[ 0], SR3l[ 0], weight*lepSFSS, weight*lepSF3l);
	fillSRhisto(histos, "YieldsSR_raw",       sn, sn2, SRSS[ 0], SR3l[ 0], 1., 1.);
	fillSRhisto(histos, "YieldsSR_rawweight", sn, sn2, SRSS[ 0], SR3l[ 0], rawweight, rawweight);
	fillSRhisto(histos, "YieldsSR_jesup",     sn, sn2, SRSS[ 6], SR3l[ 6], weight*lepSFSS, weight*lepSF3l);
	fillSRhisto(histos, "YieldsSR_jesdn",     sn, sn2, SRSS[10], SR3l[10], weight*lepSFSS, weight*lepSF3l);
	if(SR3l[ 0]>=0||SRSS[ 1]>=0){
	  int t = SRSS[ 1];
	  if(MjjL>400.||Detajj>1.5) t = -1;
	  if(SRSS[ 1]==0&&(MET.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)            ) t = -1;
	  if(SRSS[ 1]==1&&(MET.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=30.||MTmax<=90.)) t = -1;
	  if(SRSS[ 1]==2&&(               (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)            ) t = -1;
	  fillSRhisto(histos, "YieldsSR_dropMjj",   sn, sn2, t, SR3l[ 0], weight*lepSFSS, weight*lepSF3l);
	}
	if(SR3l[ 6]>=0||SRSS[ 7]>=0){
	  int t = SRSS[ 7];
	  if(MjjL_up>400.||Detajj_up>1.5) t = -1;
	  if(SRSS[ 7]==0&&(MET_up.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)               ) t = -1;
	  if(SRSS[ 7]==1&&(MET_up.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=30.||MTmax_up<=90.)) t = -1;
	  if(SRSS[ 7]==2&&(                  (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)               ) t = -1;
	  fillSRhisto(histos, "YieldsSR_dropMjj_jesup",   sn, sn2, t, SR3l[ 6], weight*lepSFSS, weight*lepSF3l);
	}
	if(SR3l[10]>=0||SRSS[11]>=0){
	  int t = SRSS[11];
	  if(MjjL_dn>400.||Detajj_dn>1.5) t = -1;
	  if(SRSS[11]==0&&(MET_dn.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)               ) t = -1;
	  if(SRSS[11]==1&&(MET_dn.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=30.||MTmax_dn<=90.)) t = -1;
	  if(SRSS[11]==2&&(                  (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)               ) t = -1;
	  fillSRhisto(histos, "YieldsSR_dropMjj_jesdn",   sn, sn2, t, SR3l[10], weight*lepSFSS, weight*lepSF3l);
	}	
	fillSRhisto(histos, "YieldsSR_lepSFup",     sn, sn2, SRSS[ 0], SR3l[ 0], weight*(lepSFSS+lepSFerrSS), weight*(lepSF3l+lepSFerr3l));
	fillSRhisto(histos, "YieldsSR_lepSFdn",     sn, sn2, SRSS[ 0], SR3l[ 0], weight*(lepSFSS-lepSFerrSS), weight*(lepSF3l-lepSFerr3l));
	//MSFOS - SR blinded
	if(SR3l[ 0]>=0&&MSFOSvec.size()>0){
	  histos["MSFOS_all3l_"        +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  histos["MSFOS_all3l_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l+lepSFerr3l));
	  histos["MSFOS_all3l_lepSFdn_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l-lepSFerr3l));
	  if(SR3l[ 0]==1){
	    histos["MSFOS_1SFOS_"        +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	    histos["MSFOS_1SFOS_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l+lepSFerr3l));
	    histos["MSFOS_1SFOS_lepSFdn_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l-lepSFerr3l));
	  }
	  if(SR3l[ 0]==2){
	    histos["MSFOS_2SFOS_"        +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	    histos["MSFOS_2SFOS_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l+lepSFerr3l));
	    histos["MSFOS_2SFOS_lepSFdn_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l-lepSFerr3l));
	  }
	  for(unsigned int i = 0; i<MSFOSvec.size();++i){
	    histos["MSFOSall_all3l_"        +sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	    histos["MSFOSall_all3l_lepSFup_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l+lepSFerr3l));
	    histos["MSFOSall_all3l_lepSFdn_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l-lepSFerr3l));
	    if(SR3l[ 0]==2){
	      histos["MSFOSall_2SFOS_"        +sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	      histos["MSFOSall_2SFOS_lepSFup_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l+lepSFerr3l));
	      histos["MSFOSall_2SFOS_lepSFdn_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l-lepSFerr3l));
	    }
	  }
	}
      }
      //MSFOS - CR not blinded
      if(SR3l[ 2]>=0&&MSFOSvec.size()>0){
	histos["MSFOS_all3l_"        +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	histos["MSFOS_all3l_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l+lepSFerr3l));
	histos["MSFOS_all3l_lepSFdn_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l-lepSFerr3l));
	if(SR3l[ 2]==1){
	  histos["MSFOS_1SFOS_"        +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  histos["MSFOS_1SFOS_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l+lepSFerr3l));
	  histos["MSFOS_1SFOS_lepSFdn_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l-lepSFerr3l));
	}
	if(SR3l[ 2]==2){
	  histos["MSFOS_2SFOS_"        +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  histos["MSFOS_2SFOS_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l+lepSFerr3l));
	  histos["MSFOS_2SFOS_lepSFdn_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSF3l-lepSFerr3l));
	}
	for(unsigned int i = 0; i<MSFOSvec.size();++i){
	  histos["MSFOSall_all3l_"        +sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	  histos["MSFOSall_all3l_lepSFup_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l+lepSFerr3l));
	  histos["MSFOSall_all3l_lepSFdn_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l-lepSFerr3l));
	  if(SR3l[ 2]==2){
	    histos["MSFOSall_2SFOS_"        +sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	    histos["MSFOSall_2SFOS_lepSFup_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l+lepSFerr3l));
	    histos["MSFOSall_2SFOS_lepSFdn_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSF3l-lepSFerr3l));
	  }
	}
      }
      //Mjj sideband
      if(SRSS[ 1]>=0){
	int t = SRSS[ 1];
	if(MjjL>400.||Detajj>1.5) t = -1;
	if(fabs(Mjj   -80.)<20.)  t = -1;
	if(SRSS[ 1]==0&&(MET.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)            ) t = -1;
	if(SRSS[ 1]==1&&(MET.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=30.||MTmax<=90.)) t = -1;
	if(SRSS[ 1]==2&&(               (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)            ) t = -1;
	fillSRhisto(histos, "YieldsSR_Mjjsideband",   sn, sn2, t, -1, weight*lepSFSS, weight*lepSF3l);
	t = SRSS[ 1];
	if(MjjL>400.||Detajj>1.5) t = -1;
	if(fabs(Mjj   -80.)<20.)  t = -1;
	if(SRSS[ 1]==0&&(MET.Pt()>40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)            ) t = -1;
	if(SRSS[ 1]==1&&(MET.Pt()>40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=30.||MTmax<=90.)) t = -1;
	if(SRSS[ 1]==2&&(MET.Pt()>40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=40.)            ) t = -1;
	fillSRhisto(histos, "YieldsSR_Mjjsideband_lowMET",   sn, sn2, t, -1, weight*lepSFSS, weight*lepSF3l);
	t = SRSS[ 1];
	if(SRSS[ 1]!=1) t = -1;
	if(MjjL>400.||Detajj>1.5) t = -1;
	if(fabs(Mjj   -80.)<20.)  t = -1;
	if(SRSS[ 1]==1&&(MET.Pt()<=40.||(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()<=30.||MTmax>90.)) t = -1;
	fillSRhisto(histos, "YieldsSR_Mjjsideband_lowMTmax",   sn, sn2, t, -1, weight*lepSFSS, weight*lepSF3l);
      }
      
      fillSRhisto(histos, "YieldsCR",           sn2,sn2, SRSS[ 2], SR3l[ 2], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "YieldsCR_jesup",     sn2,sn2, SRSS[ 8], SR3l[ 8], weight*lepSFSS, weight*lepSF3l);
      fillSRhisto(histos, "YieldsCR_jesdn",     sn2,sn2, SRSS[12], SR3l[12], weight*lepSFSS, weight*lepSF3l);
      if(fabs(Mjj   -80.)<20.||SR3l[ 2]>=0) fillSRhisto(histos, "YieldsCR_cutonMjj",       sn2,sn2, SRSS[ 2], SR3l[ 2], weight*lepSFSS, weight*lepSF3l);
      if(fabs(Mjj_up-80.)<20.||SR3l[ 8]>=0) fillSRhisto(histos, "YieldsCR_cutonMjj_jesup", sn2,sn2, SRSS[ 8], SR3l[ 8], weight*lepSFSS, weight*lepSF3l);
      if(fabs(Mjj_dn-80.)<20.||SR3l[12]>=0) fillSRhisto(histos, "YieldsCR_cutonMjj_jesdn", sn2,sn2, SRSS[12], SR3l[12], weight*lepSFSS, weight*lepSF3l);
      if(SRSS[ 2]>=0) histos["Mjj_CRlike_allSS_"      +sn2]->Fill(Mjj,    weight*lepSFSS);
      if(SRSS[ 8]>=0) histos["Mjj_CRlike_allSS_jesup_"+sn2]->Fill(Mjj_up, weight*lepSFSS);
      if(SRSS[12]>=0) histos["Mjj_CRlike_allSS_jesdn_"+sn2]->Fill(Mjj_dn, weight*lepSFSS);
      fillSRhisto(histos, "YieldsCR_lepSFup",     sn2,sn2, SRSS[ 2], SR3l[ 2], weight*(lepSFSS+lepSFerrSS), weight*(lepSF3l+lepSFerr3l));
      fillSRhisto(histos, "YieldsCR_lepSFdn",     sn2,sn2, SRSS[ 2], SR3l[ 2], weight*(lepSFSS-lepSFerrSS), weight*(lepSF3l-lepSFerr3l));
      if(SRSS[ 4]>=0&&MSFOSvec.size()>0){
	histos["MSFOS_CRlike_allSS_"        +sn2]->Fill(MSFOSvec[0],    weight*lepSFSS);
	histos["MSFOS_CRlike_allSS_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSFSS+lepSFerrSS));
	histos["MSFOS_CRlike_allSS_lepSFup_"+sn2]->Fill(MSFOSvec[0],    weight*(lepSFSS-lepSFerrSS));
	for(unsigned int i = 0; i<MSFOSvec.size();++i){
	  histos["MSFOSall_CRlike_allSS_"        +sn2]->Fill(MSFOSvec[i],    weight*lepSFSS);
	  histos["MSFOSall_CRlike_allSS_lepSFup_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSFSS+lepSFerrSS));
	  histos["MSFOSall_CRlike_allSS_lepSFup_"+sn2]->Fill(MSFOSvec[i],    weight*(lepSFSS-lepSFerrSS));
	}
      }
      if((SR3l[ 1]>=0||SR3l[ 3]>=0)&&MSFOSvec.size()>0&&fabs(Mlll-90.)>10.){
	bool passMET = true;
	if(SR3l[ 1]==1&&MET.Pt()<45.) passMET = false;
	if(SR3l[ 1]==2&&MET.Pt()<55.) passMET = false;
	if(SR3l[ 3]==1&&MET.Pt()<45.) passMET = false;
	if(SR3l[ 3]==2&&MET.Pt()<55.) passMET = false;
	bool passPTlll = (pTlll>=60.);
	bool passDPhilllMET = (DPhilllMET>=2.7);
	if((SR3l[ 1]==1||SR3l[ 3]==1)&&(DPhilllMET>=2.5)) passDPhilllMET = true;
	bool passneither  = !passMET&&!passPTlll&&!passDPhilllMET;
	bool passnotall   = !(passMET&&passPTlll&&passDPhilllMET);
	
	if(passneither){
	  histos["MSFOS_all3l_invertMETdPhiPt_"       +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  if((SR3l[ 1]==1)||(SR3l[ 3]==1)){
	    histos["MSFOS_1SFOS_invertMETdPhiPt_"     +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  }
	  if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	    histos["MSFOS_2SFOS_invertMETdPhiPt_"     +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  }
	  for(unsigned int i = 0; i<MSFOSvec.size();++i){
	    histos["MSFOSall_all3l_invertMETdPhiPt_"  +sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	    if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	      histos["MSFOSall_2SFOS_invertMETdPhiPt_"+sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	    }
	  }
	}
	if(passnotall){
	  histos["MSFOS_all3l_inverteitherMETdPhiPt_"       +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  if((SR3l[ 1]==1)||(SR3l[ 3]==1)){
	    histos["MSFOS_1SFOS_inverteitherMETdPhiPt_"     +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  }
	  if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	    histos["MSFOS_2SFOS_inverteitherMETdPhiPt_"     +sn2]->Fill(MSFOSvec[0],    weight*lepSF3l);
	  }
	  for(unsigned int i = 0; i<MSFOSvec.size();++i){
	    histos["MSFOSall_all3l_inverteitherMETdPhiPt_"  +sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	    if((SR3l[ 1]==2)||(SR3l[ 3]==2)){
	      histos["MSFOSall_2SFOS_inverteitherMETdPhiPt_"+sn2]->Fill(MSFOSvec[i],    weight*lepSF3l);
	    }
	  }
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

  SaveHistosToFile("rootfiles/Check3lCRnew.root",histos,true,true);

  /*
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
  */
  // return
  deleteFiles((SFloaded>=0),fSF);
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
