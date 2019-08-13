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
#include "Functions.h"
//#include "CMS3_WWW0117.cc"
#include "CMS3_WWW100.cc"
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1) {

  bool blindSR = false;
  bool btagreweighting = true;
  bool applylepSF      = false;
  bool applytrigSF     = false;
  bool applyPUrewgt    = true;

  const char* json_file = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();

  histonames.push_back("SR_addone_SSee");                 hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMET_SSee");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMjj_SSee");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMll_SSee");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropDeta_SSee");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_SSem");                 hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMET_SSem");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMjj_SSem");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMll_SSem");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropDeta_SSem");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMTmax_SSem");       hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_SSmm");                 hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMET_SSmm");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMjj_SSmm");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMll_SSmm");         hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropDeta_SSmm");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);

  histonames.push_back("SR_additionaltests_SSee");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_additionaltests_SSem");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_additionaltests_SSmm");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);

  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir);
  cout << "Loaded histograms" << endl;

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
      float PUweight(1.), PUweightup(1.), PUweightdn(1.);
      if(applyPUrewgt&&!isData()){
	PUweight = getPUWeightAndError(PUweightdn,PUweightup);
	weight *= PUweight;
      }

      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      float Mjj = -1;
      float MjjL = -1; float Detajj = -1;
      getMjjAndDeta(Mjj,MjjL,Detajj);

      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      //getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l);//oldID
      //getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l,1,25,20);//newID
      getleptonindices_v2(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l,25,20);
      float lepSF(1.), lepSFerr(0.);//i3l and iSS have same ID
      if(applylepSF&&!isData()){
	lepSF = getlepSFWeightandError(lepSFerr,i3l,ia3l);
	weight *= lepSF;
      }
      float trigSF(1.), trigSFerr(1.);
      if(applytrigSF&&!isData()){
	trigSF    = getTriggerWeightandError(trigSFerr, i3l,ia3l);
	weight *= trigSF;
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

      string sample   = skimFilePrefix;
      string sn       = ((iSS.size()+iaSS.size())>=2) ? process(fname,true ,iSS,iaSS) : string("not2l");
      string sn2      = ((i3l.size()+ia3l.size())>=3) ? process(fname,false,i3l,ia3l) : string("not3l");
      bool isphotonSS = (sn =="photonfakes");
      bool isphoton3l = (sn2=="photonfakes");
      if(splitVH(fname)){ sample = "WHtoWWW"; }

      float MTmax = -1;
      if(iSS.size()==2) MTmax = calcMTmax(iSS,MET);
      else if(iSS.size()==1&&iaSS.size()>=1){
	vector<int> temp; temp.push_back(iSS[0]); temp.push_back(iaSS[0]);
	MTmax = calcMTmax(temp,MET);
      }
      vector<float> MTs;
      for(unsigned int i = 0; i<i3l.size(); ++i){
	MTs.push_back(mT(lep_p4()[i3l[i] ],MET));
      }
     sort(MTs.begin(),MTs.end(),sortDecreasing);

      int SRSSpresel = -1;
      int SR3lpresel = -1;
      int CR3lpresel = -1;
      //SS
      SRSSpresel = isSRSS(iSS,      vSS,true ,MTmax,  nj30,nb,Mjj,MjjL,Detajj,MET,0,false,false,1);//enter variables for quicker calculation
      //3l
      SR3lpresel = isSR3l(i3l,true ,nj,nb,MET,0,false,1);
      CR3lpresel = isCR3l(i3l,true ,nj,nb,MET,0,false,1);

      if(vetophotonprocess(fname,isphotonSS))    { SRSSpresel = -1; }
      if(vetophotonprocess(fname,isphoton3l))    { SR3lpresel = -1; }

      if(SRSSpresel>=0){
	string mysample = "";
	string mysn = "";
	if(SRSSpresel==0) { mysn = "SSee_"+sn; mysample = "SSee_"+sample; }
	if(SRSSpresel==1) { mysn = "SSem_"+sn; mysample = "SSem_"+sample; }
	if(SRSSpresel==2) { mysn = "SSmm_"+sn; mysample = "SSmm_"+sample; }
	float Mll  = (lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M();
	if(!isData()||!blindSR){
	  bool passMjj  = fabs(Mjj-80.)>=15.;//INVERTED MJJ !!!!
	  bool passMjjL = MjjL<400.;
	  bool passDeta = Detajj<1.5;
	  bool passMll  = ((SRSSpresel==1) ? ((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>30.) : ((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.));
	  bool passMET  = true;//mumu
	  if(SRSSpresel==0) passMET = (met_pt()>40.);//ee
	  if(SRSSpresel==1) passMET = (met_pt()>60.);//emu
	  bool passMT   = ((SRSSpresel==1) ? (MTmax>90.) : (true));
	  float minDR1(999.), minDRlj(999.), minDPhilj(999.);
	  float MTjjllMET(-1);
	  for(unsigned int n = 0; n<jets_p4().size();++n){
	    if(!(   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)) continue;
	    if(dR(  jets_p4()[n],lep_p4()[iSS[0] ])<minDRlj) minDRlj = dR(  jets_p4()[n],lep_p4()[iSS[0] ]);
	    if(dR(  jets_p4()[n],lep_p4()[iSS[1] ])<minDRlj) minDRlj = dR(  jets_p4()[n],lep_p4()[iSS[1] ]);
	    if(dPhi(jets_p4()[n],lep_p4()[iSS[0] ])<minDRlj) minDRlj = dPhi(jets_p4()[n],lep_p4()[iSS[0] ]);
	    if(dPhi(jets_p4()[n],lep_p4()[iSS[1] ])<minDRlj) minDRlj = dPhi(jets_p4()[n],lep_p4()[iSS[1] ]);
	    for(unsigned int m = n+1; m<jets_p4().size();++m){
	      if(!(   jets_p4()[m].Pt()>=30&&fabs(   jets_p4()[m].Eta())<=2.5)) continue;
	      if(dR(jets_p4()[n],jets_p4()[m])<minDR1){
		minDR1 = dR(jets_p4()[n],jets_p4()[m]);
		MTjjllMET = mT(jets_p4()[n]+jets_p4()[m]+lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ],MET);
	      }
	    }
	  }
	    
	  if(passMjj&&passMjjL&&passDeta&&passMll&&passMET&&passMT) {
	                                                      fillhisto(histos, "SR_addone",     mysample, mysn, 0.,      weight);
	    if((MTs[0]+MTs[1])>90.)                           fillhisto(histos, "SR_addone",     mysample, mysn, 1.,      weight);//MTsum
	    if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone",     mysample, mysn, 2.,      weight);//MTmax
	    if(MTs[1]>50.)                                    fillhisto(histos, "SR_addone",     mysample, mysn, 3.,      weight);//MTmin
	    if(dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5) fillhisto(histos, "SR_addone",     mysample, mysn, 4.,      weight);//DPhi(l,l)
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>20.) fillhisto(histos, "SR_addone",     mysample, mysn, 5.,      weight);//Mll
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) fillhisto(histos, "SR_addone",     mysample, mysn, 6.,      weight);//Mll
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.) fillhisto(histos, "SR_addone",     mysample, mysn, 7.,      weight);//Mll
	    if(Mjj<120.)                                      fillhisto(histos, "SR_addone",     mysample, mysn, 8.,      weight);//Mjj
	    if(fabs(Mjj-80.)<15.)                             fillhisto(histos, "SR_addone",     mysample, mysn, 9.,      weight);//Mjj
	    if(met_pt()>40.)                                  fillhisto(histos, "SR_addone",     mysample, mysn,10.,      weight);//MET
	    if(met_pt()>60.)                                  fillhisto(histos, "SR_addone",     mysample, mysn,11.,      weight);//MET
	    if(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5) fillhisto(histos, "SR_addone",     mysample, mysn,12.,      weight);//MET/pTll
	    if(MTjjllMET> 80.)                                          fillhisto(histos, "SR_addone",     mysample, mysn,13.,      weight);//MTjjllMET
	    if(MTjjllMET>100.)                                          fillhisto(histos, "SR_addone",     mysample, mysn,14.,      weight);//MTjjllMET
	    if(minDRlj<1.5)                                             fillhisto(histos, "SR_addone",     mysample, mysn,15.,      weight);//minDRlj
	    if(minDPhilj>0.5)                                           fillhisto(histos, "SR_addone",     mysample, mysn,16.,      weight);//minDPhilj
	    if(minDPhilj<1.5)                                           fillhisto(histos, "SR_addone",     mysample, mysn,17.,      weight);//minDPhilj
	    if(minDPhilj<1.5&&minDPhilj>0.5)                            fillhisto(histos, "SR_addone",     mysample, mysn,18.,      weight);//minDPhilj
	  }
	  if(passMjj&&passMjjL&&passDeta&&passMll         &&passMT){
	                                                      fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 0.,      weight);
	    if((MTs[0]+MTs[1])>90.)                           fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 1.,      weight);
	    if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 2.,      weight);
	    if(MTs[1]>50.)                                    fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 3.,      weight);
	    if(dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5) fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 4.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>20.) fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 5.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 6.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.) fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 7.,      weight);
	    if(Mjj<120.)                                      fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 8.,      weight);
	    if(fabs(Mjj-80.)<15.)                             fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 9.,      weight);
	    if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,10.,      weight);
	    if(met_pt()>60.)                                  fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,11.,      weight);
	    if(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5) fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,12.,      weight);//MET/pTll
	    if(MTjjllMET> 80.)                                          fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,13.,      weight);//MTjjllMET
	    if(MTjjllMET>100.)                                          fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,14.,      weight);//MTjjllMET
	    if(minDRlj<1.5)                                             fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,15.,      weight);//minDRlj
	    if(minDPhilj>0.5)                                           fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,16.,      weight);//minDPhilj
	    if(minDPhilj<1.5)                                           fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,17.,      weight);//minDPhilj
	    if(minDPhilj<1.5&&minDPhilj>0.5)                            fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,18.,      weight);//minDPhilj
	  }
	  if(         passMjjL&&passDeta&&passMll&&passMET&&passMT){
	                                                      fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 0.,      weight);
	    if((MTs[0]+MTs[1])>90.)                           fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 1.,      weight);
	    if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 2.,      weight);
	    if(MTs[1]>50.)                                    fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 3.,      weight);
	    if(dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5) fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 4.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>20.) fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 5.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 6.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.) fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 7.,      weight);
	    if(Mjj<120.)                                      fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 8.,      weight);
	    if(fabs(Mjj-80.)<15.)                             fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn, 9.,      weight);
	    if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,10.,      weight);
	    if(met_pt()>60.)                                  fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,11.,      weight);
	    if(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5) fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,12.,      weight);//MET/pTll
	    if(MTjjllMET> 80.)                                          fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,13.,      weight);//MTjjllMET
	    if(MTjjllMET>100.)                                          fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,14.,      weight);//MTjjllMET
	    if(minDRlj<1.5)                                             fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,15.,      weight);//minDRlj
	    if(minDPhilj>0.5)                                           fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,16.,      weight);//minDPhilj
	    if(minDPhilj<1.5)                                           fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,17.,      weight);//minDPhilj
	    if(minDPhilj<1.5&&minDPhilj>0.5)                            fillhisto(histos, "SR_addone_dropMjj",     mysample, mysn,18.,      weight);//minDPhilj
	  }
	  if(passMjj&&passMjjL&&passDeta         &&passMET&&passMT){
	                                                      fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 0.,      weight);
	    if((MTs[0]+MTs[1])>90.)                           fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 1.,      weight);
	    if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 2.,      weight);
	    if(MTs[1]>50.)                                    fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 3.,      weight);
	    if(dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5) fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 4.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>20.) fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 5.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 6.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.) fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 7.,      weight);
	    if(Mjj<120.)                                      fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 8.,      weight);
	    if(fabs(Mjj-80.)<15.)                             fillhisto(histos, "SR_addone_dropMll",     mysample, mysn, 9.,      weight);
	    if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,10.,      weight);
	    if(met_pt()>60.)                                  fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,11.,      weight);
	    if(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5) fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,12.,      weight);//MET/pTll
	    if(MTjjllMET> 80.)                                          fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,13.,      weight);//MTjjllMET
	    if(MTjjllMET>100.)                                          fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,14.,      weight);//MTjjllMET
	    if(minDRlj<1.5)                                             fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,15.,      weight);//minDRlj
	    if(minDPhilj>0.5)                                           fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,16.,      weight);//minDPhilj
	    if(minDPhilj<1.5)                                           fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,17.,      weight);//minDPhilj
	    if(minDPhilj<1.5&&minDPhilj>0.5)                            fillhisto(histos, "SR_addone_dropMll",     mysample, mysn,18.,      weight);//minDPhilj
	  }
	  if(passMjj&&passMjjL          &&passMll&&passMET&&passMT){
	                                                      fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 0.,      weight);
	    if((MTs[0]+MTs[1])>90.)                           fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 1.,      weight);
	    if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 2.,      weight);
	    if(MTs[1]>50.)                                    fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 3.,      weight);
	    if(dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5) fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 4.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>20.) fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 5.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 6.,      weight);
	    if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.) fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 7.,      weight);
	    if(Mjj<120.)                                      fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 8.,      weight);
	    if(fabs(Mjj-80.)<15.)                             fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn, 9.,      weight);
	    if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,10.,      weight);
	    if(met_pt()>60.)                                  fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,11.,      weight);
	    if(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5) fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,12.,      weight);//MET/pTll
	    if(MTjjllMET> 80.)                                          fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,13.,      weight);//MTjjllMET
	    if(MTjjllMET>100.)                                          fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,14.,      weight);//MTjjllMET
	    if(minDRlj<1.5)                                             fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,15.,      weight);//minDRlj
	    if(minDPhilj>0.5)                                           fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,16.,      weight);//minDPhilj
	    if(minDPhilj<1.5)                                           fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,17.,      weight);//minDPhilj
	    if(minDPhilj<1.5&&minDPhilj>0.5)                            fillhisto(histos, "SR_addone_dropDeta",     mysample, mysn,18.,      weight);//minDPhilj

	  }
	  if(SRSSpresel==1){
	    if(passMjj&&passMjjL&&passDeta&&passMll&&passMET      ){
	                                                        fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 0.,      weight);
	      if((MTs[0]+MTs[1])>90.)                           fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 1.,      weight);
	      if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 2.,      weight);
	      if(MTs[1]>50.)                                    fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 3.,      weight);
	      if(dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5) fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 4.,      weight);
	      if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>20.) fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 5.,      weight);
	      if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>40.) fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 6.,      weight);
	      if((lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.) fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 7.,      weight);
	      if(Mjj<120.)                                      fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 8.,      weight);
	      if(fabs(Mjj-80.)<15.)                             fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn, 9.,      weight);
	      if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,10.,      weight);
	      if(met_pt()>60.)                                  fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,11.,      weight);
	    if(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5) fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,12.,      weight);//MET/pTll
	    if(MTjjllMET> 80.)                                          fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,13.,      weight);//MTjjllMET
	    if(MTjjllMET>100.)                                          fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,14.,      weight);//MTjjllMET
	    if(minDRlj<1.5)                                             fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,15.,      weight);//minDRlj
	    if(minDPhilj>0.5)                                           fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,16.,      weight);//minDPhilj
	    if(minDPhilj<1.5)                                           fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,17.,      weight);//minDPhilj
	    if(minDPhilj<1.5&&minDPhilj>0.5)                            fillhisto(histos, "SR_addone_dropMTmax",     mysample, mysn,18.,      weight);//minDPhilj
	    }
	  }
	  //additional tests SSee
	  if(SRSSpresel==0){
	    if(passMjj&&passMjjL&&passDeta&&passMll&&passMT)                                                            fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//drop MET, Mjj
	    if(passMjj&&passMjjL&&passDeta&&passMll&&passMT&&MTs[1]>50.)                                                fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//drop MET, Mjj, add MTmin>50
	    if(passMjj&&passMjjL&&passDeta&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5)                      fillhisto(histos, "SR_additionaltests",     mysample, mysn, 2.,      weight);//drop MET, Mjj, Mll, add DPhi(l,l)>0.5
	    if(passMjj&&passMjjL&&passDeta&&passMET&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 3.,      weight);//drop Mjj, Mll, add DPhi(l,l)>0.5
	    if(passMjj&&passMjjL&&passDeta&&passMjj&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 4.,      weight);//drop MET, Mll, add DPhi(l,l)>0.5
	    if(passMjj&&passMjjL&&passDeta&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5&&MTs[1]>50.)          fillhisto(histos, "SR_additionaltests",     mysample, mysn, 5.,      weight);//drop MET, Mjj, Mll, add DPhi(l,l)>0.5, add MTmin>50.
	    if(passMjj&&passMjjL&&passDeta&&passMT&&passMET&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5&&MTs[1]>50.) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 6.,      weight);//drop Mjj, Mll, add DPhi(l,l)>0.5, add MTmin>50.
	    if(passMjj&&passMjjL&&passDeta&&passMT&&passMjj&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5&&MTs[1]>50.) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 7.,      weight);//drop MET, Mll, add DPhi(l,l)>0.5, add MTmin>50.
	  }
	  //additional tests SSem
	  if(SRSSpresel==1){
	    bool passMjj2 = fabs(Mjj-80.)>=15.;//INVERT MJJ
	    //if(passMjj&&passMjjL&&passDeta&&passMll&&passMET&&passMT)
	    if(passMjj2&&passMjjL&&passDeta&&passMll&&passMT)               fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//drop MET, tighten Mjj
	    if(passMjj2&&passMjjL&&passDeta&&passMll&&passMT&&met_pt()>60.) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//MET>60, tighten Mjj
	    if(passMjj2&&passMjjL&&passDeta&&passMll&&met_pt()>60.)         fillhisto(histos, "SR_additionaltests",     mysample, mysn, 2.,      weight);//MET>60, tighten Mjj, drop MTmax
	    if(passMjj2&&passMjjL&&passDeta&&passMT)                        fillhisto(histos, "SR_additionaltests",     mysample, mysn, 3.,      weight);//drop MET, tighten Mjj, drop Mll
	    if(passMjj2&&passMjjL&&passDeta&&passMT&&met_pt()>60.)          fillhisto(histos, "SR_additionaltests",     mysample, mysn, 4.,      weight);//MET>60, tighten Mjj, drop Mll
	    if(passMjj2&&passMjjL&&passDeta&&met_pt()>60.)                  fillhisto(histos, "SR_additionaltests",     mysample, mysn, 5.,      weight);//MET>60, tighten Mjj, drop Mll, drop MTmax
	    if(passMjj&&passMjjL&&passDeta&&met_pt()>60.)                   fillhisto(histos, "SR_additionaltests",     mysample, mysn, 6.,      weight);//MET>60, drop Mll, drop MTmax
	    if(passMjj2&&passMjjL&&passDeta&&passMll&&passMT&&(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5)) fillhisto(histos, "SR_additionaltests",     mysample, mysn,7.,      weight);//MET/pTll>0.5, tighten Mjj

	  }
	  //additional tests SSmm
	  if(SRSSpresel==2){
	    bool passMjj2 = fabs(Mjj-80.)>=15.;//INVERT MJJ
	    if(passMjj2&&passMjjL&&passDeta&&passMET&&passMT&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//tighten Mjj, Mll > 60
	    if(Mjj<120.&&passMjj2&&passMjjL&&passDeta&&passMET&&passMT&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//relax Mjj, Mll > 60
	    if(passMjj2&&passMjjL&&passDeta&&passMll&&passMT&&(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5)) fillhisto(histos, "SR_additionaltests",     mysample, mysn,3.,      weight);//MET/pTll>0.5, tighten Mjj
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
  
  SaveHistosToFile("rootfiles/ImproveStuffv2NewBabe_Mjjsideband.root",histos,true,true,(chainnumber==0));
  deleteHistograms(histos);

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
