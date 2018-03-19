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

  histonames.push_back("SR_addone_0SFOS");                hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMET_0SFOS");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_droppTlll_0SFOS");      hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropDPhilllMET_0SFOS"); hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_1SFOS");                hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMET_1SFOS");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_changeMSFOS_1SFOS");    hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_droppTlll_1SFOS");      hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropDPhilllMET_1SFOS"); hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMT3rd_1SFOS");      hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_2SFOS");                hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropMET_2SFOS");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_changeMSFOS_2SFOS");    hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_droppTlll_2SFOS");      hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_addone_dropDPhilllMET_2SFOS"); hbins.push_back(50); hlow.push_back(0); hup.push_back(50);

  histonames.push_back("SR_additionaltests_SSee");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_additionaltests_SSem");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_additionaltests_SSmm");        hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_additionaltests_0SFOS");       hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_additionaltests_1SFOS");       hbins.push_back(50); hlow.push_back(0); hup.push_back(50);
  histonames.push_back("SR_additionaltests_2SFOS");       hbins.push_back(50); hlow.push_back(0); hup.push_back(50);

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
	  bool passMjj  = fabs(Mjj-80.)<15.;
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
	    if(passMjjL&&passDeta&&passMll&&passMT)                                                            fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//drop MET, Mjj
	    if(passMjjL&&passDeta&&passMll&&passMT&&MTs[1]>50.)                                                fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//drop MET, Mjj, add MTmin>50
	    if(passMjjL&&passDeta&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5)                      fillhisto(histos, "SR_additionaltests",     mysample, mysn, 2.,      weight);//drop MET, Mjj, Mll, add DPhi(l,l)>0.5
	    if(passMjjL&&passDeta&&passMET&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 3.,      weight);//drop Mjj, Mll, add DPhi(l,l)>0.5
	    if(passMjjL&&passDeta&&passMjj&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 4.,      weight);//drop MET, Mll, add DPhi(l,l)>0.5
	    if(passMjjL&&passDeta&&passMT&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5&&MTs[1]>50.)          fillhisto(histos, "SR_additionaltests",     mysample, mysn, 5.,      weight);//drop MET, Mjj, Mll, add DPhi(l,l)>0.5, add MTmin>50.
	    if(passMjjL&&passDeta&&passMT&&passMET&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5&&MTs[1]>50.) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 6.,      weight);//drop Mjj, Mll, add DPhi(l,l)>0.5, add MTmin>50.
	    if(passMjjL&&passDeta&&passMT&&passMjj&&dPhi(lep_p4()[iSS[0] ],lep_p4()[iSS[1] ])>0.5&&MTs[1]>50.) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 7.,      weight);//drop MET, Mll, add DPhi(l,l)>0.5, add MTmin>50.
	  }
	  //additional tests SSem
	  if(SRSSpresel==1){
	    bool passMjj2 = fabs(Mjj-80.)<15.;
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
	    bool passMjj2 = fabs(Mjj-80.)<15.;
	    if(passMjj2&&passMjjL&&passDeta&&passMET&&passMT&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//tighten Mjj, Mll > 60
	    if(Mjj<120.&&passMjjL&&passDeta&&passMET&&passMT&&(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).M()>60.)             fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//relax Mjj, Mll > 60
	    if(passMjj2&&passMjjL&&passDeta&&passMll&&passMT&&(met_pt()/(lep_p4()[iSS[0] ]+lep_p4()[iSS[1] ]).Pt()>0.5)) fillhisto(histos, "SR_additionaltests",     mysample, mysn,3.,      weight);//MET/pTll>0.5, tighten Mjj
	  }
	}
      }
      bool SFOS12(false), SFOS13(false), SFOS23(false); int nSFOS = 0;
      vector<LorentzVector> l;
      int mZ1(-1), mZ2(-1), mZ3(-1);
      int lZ1(-1), lZ2(-1), lZ3(-1);
      float pTSFOSmax = -1;
      if(i3l.size()==3){
	l.push_back(lep_p4()[i3l[0] ]); l.push_back(lep_p4()[i3l[1] ]); l.push_back(lep_p4()[i3l[2] ]);
	if(lep_pdgId()[i3l[0] ]==(-lep_pdgId()[i3l[1] ])) { SFOS12 = true; ++nSFOS; }
	if(lep_pdgId()[i3l[0] ]==(-lep_pdgId()[i3l[2] ])) { SFOS13 = true; ++nSFOS; }
	if(lep_pdgId()[i3l[1] ]==(-lep_pdgId()[i3l[2] ])) { SFOS23 = true; ++nSFOS; }
	if(nSFOS==0) { mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	if(nSFOS==1) {
	  if(SFOS12) { mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	  if(SFOS13) { mZ1 = 0; mZ2 = 2; mZ3 = 1; lZ1 = 0; lZ2 = 2; lZ3 = 1; }
	  if(SFOS23) { mZ1 = 1; mZ2 = 2; mZ3 = 0; lZ1 = 1; lZ2 = 2; lZ3 = 0; }
	}
	if(nSFOS==2) {
	  if(SFOS12&&SFOS13){
	    if(fabs((l[0]+l[1]).M()-MZ)<fabs((l[0]+l[2]).M()-MZ)){ mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 0; lZ2 = 2; lZ3 = 1; }
	    else                                                 { mZ1 = 0; mZ2 = 2; mZ3 = 1; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	  }
	  if(SFOS12&&SFOS23){
	    if(fabs((l[0]+l[1]).M()-MZ)<fabs((l[1]+l[2]).M()-MZ)){ mZ1 = 0; mZ2 = 1; mZ3 = 2; lZ1 = 1; lZ2 = 2; lZ3 = 0; }
	    else                                                 { mZ1 = 1; mZ2 = 2; mZ3 = 0; lZ1 = 0; lZ2 = 1; lZ3 = 2; }
	  }
	  if(SFOS13&&SFOS23){
	    if(fabs((l[0]+l[2]).M()-MZ)<fabs((l[1]+l[2]).M()-MZ)){ mZ1 = 0; mZ2 = 2; mZ3 = 1; lZ1 = 1; lZ2 = 2; lZ3 = 0; }
	    else                                                 { mZ1 = 1; mZ2 = 2; mZ3 = 0; lZ1 = 0; lZ2 = 2; lZ3 = 1; }
	  }
	}
       pTSFOSmax = ((l[mZ1]+l[mZ2]).Pt() > (l[lZ1]+l[lZ2]).Pt() ? (l[mZ1]+l[mZ2]).Pt() : (l[lZ1]+l[lZ2]).Pt() );
      }
      if(SR3lpresel>=0){
	string mysample = "";
	string mysn = "";
	if(SR3lpresel==0) { mysn = "0SFOS_"+sn2; mysample = "0SFOS_"+sample; }
	if(SR3lpresel==1) { mysn = "1SFOS_"+sn2; mysample = "1SFOS_"+sample; }
	if(SR3lpresel==2) { mysn = "2SFOS_"+sn2; mysample = "2SFOS_"+sample; }
	bool passZero  = true;
	if(SR3lpresel==0){
	  bool SF01 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[1] ]));
	  bool SF02 = (abs(lep_pdgId()[i3l[0] ])==abs(lep_pdgId()[i3l[2] ]));
	  bool SF12 = (abs(lep_pdgId()[i3l[1] ])==abs(lep_pdgId()[i3l[2] ]));
	  if(SF01&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()<20.) passZero = false;
	  if(SF02&&(lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()<20.) passZero = false;
	  if(SF12&&(lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()<20.) passZero = false;
	  if(SF01&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]).M()-MZ)<15.) passZero = false;
	  if(SF02&&abs(lep_pdgId()[i3l[0] ])==11&&fabs((lep_p4()[i3l[0] ]+lep_p4()[i3l[2] ]).M()-MZ)<15.) passZero = false;
	  if(SF12&&abs(lep_pdgId()[i3l[1] ])==11&&fabs((lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M()-MZ)<15.) passZero = false;
	}
	if(passZero){
	  float DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	  float pTlll      = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	  float Mlll       = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M();
	  if(!isData()||!blindSR){
	    bool passMET   = (met_pt()>30.); if(SR3lpresel==1) passMET = (met_pt()>40.); if(SR3lpresel==2) passMET = (met_pt()>55.);
	    bool passDPhi  = ((SR3lpresel==0) ? (DPhilllMET>2.7) : (DPhilllMET>2.5));
	    bool passpTlll = ((SR3lpresel==0) ? (true) : (pTlll>60.));
	    bool passMlll  = ((SR3lpresel>=1) ? (fabs(Mlll-MZ)>10.) : (true));
	    bool passMT3rd = ((SR3lpresel==1) ? (mT(l[mZ3],MET)>90.) : (true));
	    if(passMET&&passDPhi&&passpTlll&&passMlll&&passMT3rd){
	                                                        fillhisto(histos, "SR_addone",     mysample, mysn, 0.,      weight);
	      if((l[mZ3]+MET).Pt()>40.&&nSFOS>=1)               fillhisto(histos, "SR_addone",     mysample, mysn, 1.,      weight);//pTW
	      if((l[mZ3]+MET).Pt()>60.&&nSFOS>=1)               fillhisto(histos, "SR_addone",     mysample, mysn, 2.,      weight);//pTW
	      if((l[mZ1]+l[mZ2]).Pt()>40.&&nSFOS>=1)            fillhisto(histos, "SR_addone",     mysample, mysn, 3.,      weight);//pTSFOS
	      if((l[mZ1]+l[mZ2]).Pt()>60.&&nSFOS>=1)            fillhisto(histos, "SR_addone",     mysample, mysn, 4.,      weight);//pTSFOS
	      if(pTSFOSmax>40.&&nSFOS>=2)                       fillhisto(histos, "SR_addone",     mysample, mysn, 5.,      weight);//pTSFOSmax
	      if(pTSFOSmax>60.&&nSFOS>=2)                       fillhisto(histos, "SR_addone",     mysample, mysn, 6.,      weight);//pTSFOSmax
	      if((MTs[0]+MTs[1]+MTs[2])>100.)                   fillhisto(histos, "SR_addone",     mysample, mysn, 7.,      weight);//MTsum
	      if((MTs[0]+MTs[1]+MTs[2])>120.)                   fillhisto(histos, "SR_addone",     mysample, mysn, 8.,      weight);//MTsum
	      if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone",     mysample, mysn,10.,      weight);//MTmax
	      if(MTs[2]>45.)                                    fillhisto(histos, "SR_addone",     mysample, mysn,11.,      weight);//MTmin
	      if(MTs[2]>50.)                                    fillhisto(histos, "SR_addone",     mysample, mysn,12.,      weight);//MTmin
	      if(mT(l[mZ3],MET)>50.&&nSFOS>=1)                  fillhisto(histos, "SR_addone",     mysample, mysn,13.,      weight);//MT3rd
	      if(mT(l[mZ3],MET)>90.&&nSFOS>=1)                  fillhisto(histos, "SR_addone",     mysample, mysn,14.,      weight);//MT3rd
	      if(mT(l[0]+l[1]+l[2],MET)>90.)                    fillhisto(histos, "SR_addone",     mysample, mysn,15.,      weight);//MTlll
	      if(dPhi(l[mZ3],l[mZ1]+l[mZ2])>90.&&nSFOS>=1)      fillhisto(histos, "SR_addone",     mysample, mysn,16.,      weight);//DPhi(ll,l)
	      if(pTlll>30.)                                     fillhisto(histos, "SR_addone",     mysample, mysn,17.,      weight);//pTlll
	      if(pTlll>50.)                                     fillhisto(histos, "SR_addone",     mysample, mysn,18.,      weight);//pTlll
	      if(met_pt()>30.)                                  fillhisto(histos, "SR_addone",     mysample, mysn,19.,      weight);//MET
	      if(met_pt()>40.)                                  fillhisto(histos, "SR_addone",     mysample, mysn,20.,      weight);//MET
	      if(met_pt()>50.)                                  fillhisto(histos, "SR_addone",     mysample, mysn,21.,      weight);//MET
	      if(DPhilllMET>2.1)                                fillhisto(histos, "SR_addone",     mysample, mysn,22.,      weight);//DPhilllMET
	      if(DPhilllMET>2.5)                                fillhisto(histos, "SR_addone",     mysample, mysn,23.,      weight);//DPhilllMET
	      if(DPhilllMET>2.7)                                fillhisto(histos, "SR_addone",     mysample, mysn,24.,      weight);//DPhilllMET
	      if(met_pt()/pTlll>0.5)                            fillhisto(histos, "SR_addone",     mysample, mysn,25.,      weight);//MET/pTlll
	      if((l[mZ3]).Pt()>40.&&nSFOS>=1)                   fillhisto(histos, "SR_addone",     mysample, mysn,26.,      weight);//pTl(3rd)
	    }
	    if(         passDPhi&&passpTlll&&passMlll&&passMT3rd){
	                                                        fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 0.,      weight);
	      if((l[mZ3]+MET).Pt()>40.&&nSFOS>=1)               fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 1.,      weight);
	      if((l[mZ3]+MET).Pt()>60.&&nSFOS>=1)               fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 2.,      weight);
	      if((l[mZ1]+l[mZ2]).Pt()>40.&&nSFOS>=1)            fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 3.,      weight);
	      if((l[mZ1]+l[mZ2]).Pt()>60.&&nSFOS>=1)            fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 4.,      weight);
	      if(pTSFOSmax>40.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 5.,      weight);
	      if(pTSFOSmax>60.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 6.,      weight);
	      if((MTs[0]+MTs[1]+MTs[2])>100.)                   fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 7.,      weight);
	      if((MTs[0]+MTs[1]+MTs[2])>120.)                   fillhisto(histos, "SR_addone_dropMET",     mysample, mysn, 8.,      weight);
	      if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,10.,      weight);
	      if(MTs[2]>45.)                                    fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,11.,      weight);
	      if(MTs[2]>50.)                                    fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,12.,      weight);
	      if(mT(l[mZ3],MET)>50.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,13.,      weight);
	      if(mT(l[mZ3],MET)>90.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,14.,      weight);
	      if(mT(l[0]+l[1]+l[2],MET)>90.)                    fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,15.,      weight);
	      if(dPhi(l[mZ3],l[mZ1]+l[mZ2])>90.&&nSFOS>=1)      fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,16.,      weight);
	      if(pTlll>30.)                                     fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,17.,      weight);
	      if(pTlll>50.)                                     fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,18.,      weight);
	      if(met_pt()>30.)                                  fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,19.,      weight);
	      if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,20.,      weight);
	      if(met_pt()>50.)                                  fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,21.,      weight);
	      if(DPhilllMET>2.1)                                fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,22.,      weight);
	      if(DPhilllMET>2.5)                                fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,23.,      weight);
	      if(DPhilllMET>2.7)                                fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,24.,      weight);
	      if(met_pt()/pTlll>0.5)                            fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,25.,      weight);//MET/pTlll
	      if((l[mZ3]).Pt()>40.&&nSFOS>=1)                   fillhisto(histos, "SR_addone_dropMET",     mysample, mysn,26.,      weight);//pTl(3rd)
	    }
	    if(passMET&&passDPhi           &&passMlll&&passMT3rd){
	                                                        fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 0.,      weight);
	      if((l[mZ3]+MET).Pt()>40.&&nSFOS>=1)               fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 1.,      weight);
	      if((l[mZ3]+MET).Pt()>60.&&nSFOS>=1)               fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 2.,      weight);
	      if((l[mZ1]+l[mZ2]).Pt()>40.&&nSFOS>=1)            fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 3.,      weight);
	      if((l[mZ1]+l[mZ2]).Pt()>60.&&nSFOS>=1)            fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 4.,      weight);
	      if(pTSFOSmax>40.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 5.,      weight);
	      if(pTSFOSmax>60.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 6.,      weight);
	      if((MTs[0]+MTs[1]+MTs[2])>100.)                   fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 7.,      weight);
	      if((MTs[0]+MTs[1]+MTs[2])>120.)                   fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn, 8.,      weight);
	      if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,10.,      weight);
	      if(MTs[2]>45.)                                    fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,11.,      weight);
	      if(MTs[2]>50.)                                    fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,12.,      weight);
	      if(mT(l[mZ3],MET)>50.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,13.,      weight);
	      if(mT(l[mZ3],MET)>90.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,14.,      weight);
	      if(mT(l[0]+l[1]+l[2],MET)>90.)                    fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,15.,      weight);
	      if(dPhi(l[mZ3],l[mZ1]+l[mZ2])>90.&&nSFOS>=1)      fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,16.,      weight);
	      if(pTlll>30.)                                     fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,17.,      weight);
	      if(pTlll>50.)                                     fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,18.,      weight);
	      if(met_pt()>30.)                                  fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,19.,      weight);
	      if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,20.,      weight);
	      if(met_pt()>50.)                                  fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,21.,      weight);
	      if(DPhilllMET>2.1)                                fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,22.,      weight);
	      if(DPhilllMET>2.5)                                fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,23.,      weight);
	      if(DPhilllMET>2.7)                                fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,24.,      weight);
	      if(met_pt()/pTlll>0.5)                            fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,25.,      weight);//MET/pTlll
	      if((l[mZ3]).Pt()>40.&&nSFOS>=1)                   fillhisto(histos, "SR_addone_droppTlll",     mysample, mysn,26.,      weight);//pTl(3rd)
	    }
	    if(passMET          &&passpTlll&&passMlll&&passMT3rd){
	                                                        fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 0.,      weight);
	      if((l[mZ3]+MET).Pt()>40.&&nSFOS>=1)               fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 1.,      weight);
	      if((l[mZ3]+MET).Pt()>60.&&nSFOS>=1)               fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 2.,      weight);
	      if((l[mZ1]+l[mZ2]).Pt()>40.&&nSFOS>=1)            fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 3.,      weight);
	      if((l[mZ1]+l[mZ2]).Pt()>60.&&nSFOS>=1)            fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 4.,      weight);
	      if(pTSFOSmax>40.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 5.,      weight);
	      if(pTSFOSmax>60.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 6.,      weight);
	      if((MTs[0]+MTs[1]+MTs[2])>100.)                   fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 7.,      weight);
	      if((MTs[0]+MTs[1]+MTs[2])>120.)                   fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn, 8.,      weight);
	      if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,10.,      weight);
	      if(MTs[2]>45.)                                    fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,11.,      weight);
	      if(MTs[2]>50.)                                    fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,12.,      weight);
	      if(mT(l[mZ3],MET)>50.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,13.,      weight);
	      if(mT(l[mZ3],MET)>90.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,14.,      weight);
	      if(mT(l[0]+l[1]+l[2],MET)>90.)                    fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,15.,      weight);
	      if(dPhi(l[mZ3],l[mZ1]+l[mZ2])>90.&&nSFOS>=1)      fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,16.,      weight);
	      if(pTlll>30.)                                     fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,17.,      weight);
	      if(pTlll>50.)                                     fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,18.,      weight);
	      if(met_pt()>30.)                                  fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,19.,      weight);
	      if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,20.,      weight);
	      if(met_pt()>50.)                                  fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,21.,      weight);
	      if(DPhilllMET>2.1)                                fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,22.,      weight);
	      if(DPhilllMET>2.5)                                fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,23.,      weight);
	      if(DPhilllMET>2.7)                                fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,24.,      weight);
	      if(met_pt()/pTlll>0.5)                            fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,25.,      weight);//MET/pTlll
	      if((l[mZ3]).Pt()>40.&&nSFOS>=1)                   fillhisto(histos, "SR_addone_dropDPhilllMET",     mysample, mysn,26.,      weight);//pTl(3rd)
	    }
	    if((SR3lpresel==1||CR3lpresel==1)&&(passMET&&passDPhi&&passpTlll&&passMlll)){
	                                                        fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 0.,      weight);
	      if((l[mZ3]+MET).Pt()>40.&&nSFOS>=1)               fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 1.,      weight);//pTW
	      if((l[mZ3]+MET).Pt()>60.&&nSFOS>=1)               fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 2.,      weight);//pTW
	      if((l[mZ1]+l[mZ2]).Pt()>40.&&nSFOS>=1)            fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 3.,      weight);//pTSFOS
	      if((l[mZ1]+l[mZ2]).Pt()>60.&&nSFOS>=1)            fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 4.,      weight);//pTSFOS
	      if(pTSFOSmax>40.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 5.,      weight);//pTSFOSmax
	      if(pTSFOSmax>60.&&nSFOS>=2)                       fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 6.,      weight);//pTSFOSmax
	      if((MTs[0]+MTs[1]+MTs[2])>100.)                   fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 7.,      weight);//MTsum
	      if((MTs[0]+MTs[1]+MTs[2])>120.)                   fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn, 8.,      weight);//MTsum
	      if(MTs[0]>90.)                                    fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,10.,      weight);//MTmax
	      if(MTs[2]>45.)                                    fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,11.,      weight);//MTmin
	      if(MTs[2]>50.)                                    fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,12.,      weight);//MTmin
	      if(mT(l[mZ3],MET)>50.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,13.,      weight);//MT3rd
	      if(mT(l[mZ3],MET)>90.&&nSFOS>=1)                  fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,14.,      weight);//MT3rd
	      if(mT(l[0]+l[1]+l[2],MET)>90.)                    fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,15.,      weight);//MTlll
	      if(dPhi(l[mZ3],l[mZ1]+l[mZ2])>90.&&nSFOS>=1)      fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,16.,      weight);//DPhi(ll,l)
	      if(pTlll>30.)                                     fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,17.,      weight);//pTlll
	      if(pTlll>50.)                                     fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,18.,      weight);//pTlll
	      if(met_pt()>30.)                                  fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,19.,      weight);//MET
	      if(met_pt()>40.)                                  fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,20.,      weight);//MET
	      if(met_pt()>50.)                                  fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,21.,      weight);//MET
	      if(DPhilllMET>2.1)                                fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,22.,      weight);//DPhilllMET
	      if(DPhilllMET>2.5)                                fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,23.,      weight);//DPhilllMET
	      if(DPhilllMET>2.7)                                fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,24.,      weight);//DPhilllMET
	      if(met_pt()/pTlll>0.5)                            fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,25.,      weight);//MET/pTlll
	      if((l[mZ3]).Pt()>40.&&nSFOS>=1)                   fillhisto(histos, "SR_addone_dropMT3rd",     mysample, mysn,26.,      weight);//pTl(3rd)
	    }
	    if(SR3lpresel==0){
	      if(met_pt()>30.&&passDPhi&&passMlll&&MTs[0]>90. ) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//drop pTlll, MET>30, MTmax>90
	      if(met_pt()>30.&&passDPhi&&passMlll&&MTs[2]>50. ) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//drop pTlll, MET>30, MTmin>50
	      if(met_pt()>50.&&passDPhi&&passMlll&&MTs[0]>90. ) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 2.,      weight);//drop pTlll, MET>50, MTmax>90
	      if(met_pt()>50.&&passDPhi&&passMlll&&MTs[2]>50. ) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 3.,      weight);//drop pTlll, MET>50, MTmin>50
	    }
	  }
	}
      }
      if(SR3lpresel>=1||CR3lpresel>=1){
	string mysample = "";
	string mysn = "";
	if(CR3lpresel==0) { mysn = "0SFOS_"+sn2; mysample = "0SFOS_"+sample; }
	if(CR3lpresel==1) { mysn = "1SFOS_"+sn2; mysample = "1SFOS_"+sample; }
	if(CR3lpresel==2) { mysn = "2SFOS_"+sn2; mysample = "2SFOS_"+sample; }
	if(SR3lpresel==0) { mysn = "0SFOS_"+sn2; mysample = "0SFOS_"+sample; }
	if(SR3lpresel==1) { mysn = "1SFOS_"+sn2; mysample = "1SFOS_"+sample; }
	if(SR3lpresel==2) { mysn = "2SFOS_"+sn2; mysample = "2SFOS_"+sample; }
	vector<float> MSFOSvec = allMSFOS(i3l);
	float MSFOSZ = MSFOSvec[0];
	float MSFOSZ2;
	if(MSFOSvec.size()>1) MSFOSZ2 = MSFOSvec[1];
	float DPhilllMET = dPhi(lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ],MET);
	float pTlll      = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).Pt();
	float Mlll       = (lep_p4()[i3l[0] ]+lep_p4()[i3l[1] ]+lep_p4()[i3l[2] ]).M();
	bool passMET   = true; if(SR3lpresel==1) passMET = (met_pt()>45.); if(SR3lpresel==2) passMET = (met_pt()>55.);
	bool passDPhi  = ((SR3lpresel==0) ? (DPhilllMET>2.7) : (DPhilllMET>2.5));
	bool passpTlll = (pTlll>60.);
	bool passMlll  = ((SR3lpresel>=1) ? (fabs(Mlll-MZ)>10.) : (true));
	bool passMT3rd = ((SR3lpresel==1) ? (mT(l[mZ3],MET)>90.) : (true));
	if(!isData()||!blindSR) {
	  if(passMET&&passDPhi&&passpTlll&&passMlll&&passMT3rd) {
	    if(SR3lpresel>=1)                                       fillhisto(histos, "SR_addone_changeMSFOS", mysample, mysn, 0., weight);
	    if(nSFOS==1&&fabs(MSFOSZ-MZ)>20.)                       fillhisto(histos, "SR_addone_changeMSFOS", mysample, mysn, 1., weight);
	    if(nSFOS==1&&fabs(MSFOSZ-MZ)>15.)                       fillhisto(histos, "SR_addone_changeMSFOS", mysample, mysn, 2., weight);
	    if(nSFOS==2&&fabs(MSFOSZ-MZ)>20.&&fabs(MSFOSZ2-MZ)>20.) fillhisto(histos, "SR_addone_changeMSFOS", mysample, mysn, 1., weight);
	    if(nSFOS==2&&fabs(MSFOSZ-MZ)>15.&&fabs(MSFOSZ2-MZ)>15.) fillhisto(histos, "SR_addone_changeMSFOS", mysample, mysn, 2., weight);
	  }
	}
	bool passMSFOSorg = true;
	bool passMSFOS20  = true;
	bool passMSFOS15  = true;
	if(SR3lpresel==1||CR3lpresel==1){
	  if(MSFOSZ >=55.&&MSFOSZ <=110.) passMSFOSorg = false;
	  if(fabs(MSFOSZ -MZ)<20.)        passMSFOS20  = false;
	  if(fabs(MSFOSZ -MZ)<15.)        passMSFOS15  = false;
	}
	if(SR3lpresel==2||CR3lpresel==2){
	  if(fabs(MSFOSZ2-MZ)<20.)        passMSFOSorg = false;
	  if(fabs(MSFOSZ2-MZ)<20.)        passMSFOS20  = false;
	  if(fabs(MSFOSZ2-MZ)<15.)        passMSFOS15  = false;
	  if(fabs(MSFOSZ -MZ)<20.)        passMSFOSorg = false;
	  if(fabs(MSFOSZ -MZ)<20.)        passMSFOS20  = false;
	  if(fabs(MSFOSZ -MZ)<15.)        passMSFOS15  = false;
	}
	if(SR3lpresel==1||CR3lpresel==1){
	  if(passMlll&&passMSFOSorg){
	    if(met_pt()>40.&&passDPhi&&passpTlll&&mT(l[mZ3],MET)>90.)       fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//met>40, MT3rd>90
	    if(met_pt()>40.&&passDPhi&&mT(l[mZ3],MET)>90.)                  fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//met>40, MT3rd>90, drop pTlll
	    if(met_pt()>40.&&passpTlll&&mT(l[mZ3],MET)>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn, 2.,      weight);//met>40, MT3rd>90, drop Dphi
	    if(met_pt()>40.&&DPhilllMET>2.1&&passpTlll&&mT(l[mZ3],MET)>90.) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 3.,      weight);//met>40, MT3rd>90, Dphi>2.1
	    if(met_pt()>40.&&DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)            fillhisto(histos, "SR_additionaltests",     mysample, mysn, 4.,      weight);//met>40, MT3rd>90, Dphi>2.1, drop pTlll
	    if(met_pt()>40.&&mT(l[mZ3],MET)>90.)                            fillhisto(histos, "SR_additionaltests",     mysample, mysn, 5.,      weight);//met>40, MT3rd>90, drop Dphi, drop pTlll
	    if(passMET&&mT(l[mZ3],MET)>90.)                                 fillhisto(histos, "SR_additionaltests",     mysample, mysn, 6.,      weight);//MT3rd>90, drop Dphi, drop pTlll
	    if(passMET&&DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn, 7.,      weight);//MT3rd>90, Dphi>2.1, drop pTlll
	  }
	  if(passMlll&&passMSFOS20){
	    if(met_pt()>40.&&passDPhi&&passpTlll&&mT(l[mZ3],MET)>90.)       fillhisto(histos, "SR_additionaltests",     mysample, mysn, 8.,      weight);//met>40, MT3rd>90
	    if(met_pt()>40.&&passDPhi&&mT(l[mZ3],MET)>90.)                  fillhisto(histos, "SR_additionaltests",     mysample, mysn, 9.,      weight);//met>40, MT3rd>90, drop pTlll
	    if(met_pt()>40.&&passpTlll&&mT(l[mZ3],MET)>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn,10.,      weight);//met>40, MT3rd>90, drop Dphi
	    if(met_pt()>40.&&DPhilllMET>2.1&&passpTlll&&mT(l[mZ3],MET)>90.) fillhisto(histos, "SR_additionaltests",     mysample, mysn,11.,      weight);//met>40, MT3rd>90, Dphi>2.1
	    if(met_pt()>40.&&DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)            fillhisto(histos, "SR_additionaltests",     mysample, mysn,12.,      weight);//met>40, MT3rd>90, Dphi>2.1, drop pTlll
	    if(met_pt()>40.&&mT(l[mZ3],MET)>90.)                            fillhisto(histos, "SR_additionaltests",     mysample, mysn,13.,      weight);//met>40, MT3rd>90, drop Dphi, drop pTlll
	    if(passMET&&mT(l[mZ3],MET)>90.)                                 fillhisto(histos, "SR_additionaltests",     mysample, mysn,14.,      weight);//MT3rd>90, drop Dphi, drop pTlll
	    if(passMET&&DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn,15.,      weight);// MT3rd>90, Dphi>2.1, drop pTlll
	  }
	  if(passMlll&&passMSFOS15){
	    if(met_pt()>40.&&passDPhi&&passpTlll&&mT(l[mZ3],MET)>90.)       fillhisto(histos, "SR_additionaltests",     mysample, mysn,16.,      weight);//met>40, MT3rd>90
	    if(met_pt()>40.&&passDPhi&&mT(l[mZ3],MET)>90.)                  fillhisto(histos, "SR_additionaltests",     mysample, mysn,17.,      weight);//met>40, MT3rd>90, drop pTlll
	    if(met_pt()>40.&&passpTlll&&mT(l[mZ3],MET)>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn,18.,      weight);//met>40, MT3rd>90, drop Dphi
	    if(met_pt()>40.&&DPhilllMET>2.1&&passpTlll&&mT(l[mZ3],MET)>90.) fillhisto(histos, "SR_additionaltests",     mysample, mysn,19.,      weight);//met>40, MT3rd>90, Dphi>2.1
	    if(met_pt()>40.&&DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)            fillhisto(histos, "SR_additionaltests",     mysample, mysn,20.,      weight);//met>40, MT3rd>90, Dphi>2.1, drop pTlll
	    if(met_pt()>40.&&mT(l[mZ3],MET)>90.)                            fillhisto(histos, "SR_additionaltests",     mysample, mysn,21.,      weight);//met>40, MT3rd>90, drop Dphi, drop pTlll
	    if(passMET&&mT(l[mZ3],MET)>90.)                                 fillhisto(histos, "SR_additionaltests",     mysample, mysn,22.,      weight);//MT3rd>90, drop Dphi, drop pTlll
	    if(passMET&&DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn,23.,      weight);// MT3rd>90, Dphi>2.1, drop pTlll
	    //if(passMET&&passDPhi&&passpTlll&&passMlll){
	  }
	}
	if(SR3lpresel==2||CR3lpresel==2){
	  if(passMET&&passMlll&&passMSFOSorg){
	    if(passDPhi&&passpTlll&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)       fillhisto(histos, "SR_additionaltests",     mysample, mysn, 0.,      weight);//add pTW>40, MTmax>90
	    if(passDPhi&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)                  fillhisto(histos, "SR_additionaltests",     mysample, mysn, 1.,      weight);//add pTW>40, MTmax>90, drop pTlll
	    if(passpTlll&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn, 2.,      weight);//add pTW>40, MTmax>90, drop DPhi
	    if(DPhilllMET>2.1&&passpTlll&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.) fillhisto(histos, "SR_additionaltests",     mysample, mysn, 3.,      weight);//add pTW>40, MTmax>90, DPhi>2.1
	    if(DPhilllMET>2.1&&passpTlll&&MTs[0]>90.)                        fillhisto(histos, "SR_additionaltests",     mysample, mysn, 4.,      weight);//add MTmax>90, DPhi>2.1
	    if((l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)                            fillhisto(histos, "SR_additionaltests",     mysample, mysn, 6.,      weight);//add pTW>40, MTmax>90, drop pTlll, drop DPhi
	    if(DPhilllMET>2.1&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)            fillhisto(histos, "SR_additionaltests",     mysample, mysn, 7.,      weight);//add pTW>40, MTmax>90, drop pTlll, DPhi>2.1
	    if(DPhilllMET>2.1&&MTs[0]>90.)                                   fillhisto(histos, "SR_additionaltests",     mysample, mysn, 8.,      weight);//add MTmax>90, drop pTlll, DPhi>2.1
	    if(DPhilllMET>2.1&&(l[mZ3]+MET).Pt()>40.)                        fillhisto(histos, "SR_additionaltests",     mysample, mysn, 9.,      weight);//add pTW>40, drop pTlll, DPhi>2.1
	    if(MTs[0]>90.)                                                   fillhisto(histos, "SR_additionaltests",     mysample, mysn,20.,      weight);//add MTmax>90, drop pTlll, drop DPhi
	    if((l[mZ3]+MET).Pt()>40.)                                        fillhisto(histos, "SR_additionaltests",     mysample, mysn,21.,      weight);//add pTW>40, drop pTlll, drop DPhi
	    if(true)                                                         fillhisto(histos, "SR_additionaltests",     mysample, mysn,22.,      weight);//drop pTlll, drop DPhi
	    if(DPhilllMET>2.1)                                               fillhisto(histos, "SR_additionaltests",     mysample, mysn,23.,      weight);//drop pTlll, DPhi>2.1
	    if(mT(l[mZ3],MET)>90.)                                           fillhisto(histos, "SR_additionaltests",     mysample, mysn, 5.,      weight);//drop pTlll, drop DPhi, add MT3rd>90
	    if(DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)                           fillhisto(histos, "SR_additionaltests",     mysample, mysn,28.,      weight);//drop pTlll, DPhi>2.1, add MT3rd>90
	  }
	  if(passMET&&passMlll&&passMSFOS15){
	    if(passDPhi&&passpTlll&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)       fillhisto(histos, "SR_additionaltests",     mysample, mysn,10.,      weight);//add pTW>40, MTmax>90
	    if(passDPhi&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)                  fillhisto(histos, "SR_additionaltests",     mysample, mysn,11.,      weight);//add pTW>40, MTmax>90, drop pTlll
	    if(passpTlll&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)                 fillhisto(histos, "SR_additionaltests",     mysample, mysn,12.,      weight);//add pTW>40, MTmax>90, drop DPhi
	    if(DPhilllMET>2.1&&passpTlll&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.) fillhisto(histos, "SR_additionaltests",     mysample, mysn,13.,      weight);//add pTW>40, MTmax>90, DPhi>2.1
	    if(DPhilllMET>2.1&&passpTlll&&MTs[0]>90.)                        fillhisto(histos, "SR_additionaltests",     mysample, mysn,14.,      weight);//add MTmax>90, DPhi>2.1
	    if((l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)                            fillhisto(histos, "SR_additionaltests",     mysample, mysn,16.,      weight);//add pTW>40, MTmax>90, drop pTlll, drop DPhi
	    if(DPhilllMET>2.1&&(l[mZ3]+MET).Pt()>40.&&MTs[0]>90.)            fillhisto(histos, "SR_additionaltests",     mysample, mysn,17.,      weight);//add pTW>40, MTmax>90, drop pTlll, DPhi>2.1
	    if(DPhilllMET>2.1&&MTs[0]>90.)                                   fillhisto(histos, "SR_additionaltests",     mysample, mysn,18.,      weight);//add MTmax>90, drop pTlll, DPhi>2.1
	    if(DPhilllMET>2.1&&(l[mZ3]+MET).Pt()>40.)                        fillhisto(histos, "SR_additionaltests",     mysample, mysn,19.,      weight);//add pTW>40, drop pTlll, DPhi>2.1
	    if(MTs[0]>90.)                                                   fillhisto(histos, "SR_additionaltests",     mysample, mysn,24.,      weight);//add MTmax>90, drop pTlll, drop DPhi
	    if((l[mZ3]+MET).Pt()>40.)                                        fillhisto(histos, "SR_additionaltests",     mysample, mysn,25.,      weight);//add pTW>40, drop pTlll, drop DPhi
	    if(true)                                                         fillhisto(histos, "SR_additionaltests",     mysample, mysn,26.,      weight);//drop pTlll, drop DPhi
	    if(DPhilllMET>2.1)                                               fillhisto(histos, "SR_additionaltests",     mysample, mysn,27.,      weight);//drop pTlll, DPhi>2.1
	    if(mT(l[mZ3],MET)>90.)                                           fillhisto(histos, "SR_additionaltests",     mysample, mysn,15.,      weight);//drop pTlll, drop DPhi, add MT3rd>90
	    if(DPhilllMET>2.1&&mT(l[mZ3],MET)>90.)                           fillhisto(histos, "SR_additionaltests",     mysample, mysn,29.,      weight);//drop pTlll, DPhi>2.1, add MT3rd>90
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
  
  SaveHistosToFile("rootfiles/ImproveStuffv2_NewBabe.root",histos,true,true,(chainnumber==0));
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
