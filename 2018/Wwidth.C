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
//#define USE_CMS3_WWW100 

//#include "Functions112.h"
//#include "CMS3_WWW112.cc"
#include "Functions.h"
#include "CMS3_WWW121.cc"
/*
#include "Functions112.h" 
#ifdef USE_CMS3_WWW100
#include "CMS3_WWW112.cc"
#else
#include "CMS3_WWW0118.cc"
#endif
*/
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1) {

  bool blindSR         = false;
  bool btagreweighting = true;
  bool applylepSF      = true;
  bool applytrigSF     = true;
  bool applyPUrewgt    = true;
  bool getJECunc       = true;
  
  const char* json_file = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int>    hbins;      hbins.clear();
  vector<float>  hlow;       hlow.clear();
  vector<float>  hup;        hup.clear();

  histonames.push_back("Wwidthfromdecay_SRSSpresel");  hbins.push_back(30); hlow.push_back(50); hup.push_back(110);
  histonames.push_back("WgenwidthLC_SRSSpresel");  hbins.push_back(30); hlow.push_back(50); hup.push_back(110);
  histonames.push_back("WgenwidthLC_SR3lpresel");  hbins.push_back(30); hlow.push_back(50); hup.push_back(110);
  histonames.push_back("WgenwidthLC_all");         hbins.push_back(30); hlow.push_back(50); hup.push_back(110);
  histonames.push_back("NWLC_all");                hbins.push_back(5); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("WgenwidthFC_SRSSpresel");  hbins.push_back(30); hlow.push_back(50); hup.push_back(110);
  histonames.push_back("WgenwidthFC_SR3lpresel");  hbins.push_back(30); hlow.push_back(50); hup.push_back(110);
  histonames.push_back("WgenwidthFC_all");         hbins.push_back(30); hlow.push_back(50); hup.push_back(110);
  histonames.push_back("NWFC_all");                hbins.push_back(5); hlow.push_back(0); hup.push_back(5);
  histonames.push_back("Mjj_SRSSpresel");          hbins.push_back(40); hlow.push_back(0); hup.push_back(400);
  histonames.push_back("Mjj_SR3lpresel");          hbins.push_back(40); hlow.push_back(0); hup.push_back(400);
  histonames.push_back("Mjj_all");                 hbins.push_back(40); hlow.push_back(0); hup.push_back(400);
  histonames.push_back("MjjL_SRSSpresel");         hbins.push_back(25); hlow.push_back(0); hup.push_back(1000);
  histonames.push_back("MjjL_SR3lpresel");         hbins.push_back(25); hlow.push_back(0); hup.push_back(1000);
  histonames.push_back("MjjL_all");                hbins.push_back(25); hlow.push_back(0); hup.push_back(1000);
  histonames.push_back("ST_SRSSpresel");           hbins.push_back(25); hlow.push_back(0); hup.push_back(2000);
  histonames.push_back("ST_SR3lpresel");           hbins.push_back(25); hlow.push_back(0); hup.push_back(2000);
  histonames.push_back("ST_all");                  hbins.push_back(25); hlow.push_back(0); hup.push_back(2000);
  
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
  TH1D* h_c;


  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    string fname = currentFile->GetTitle();
    h_c = (TH1D*)file->Get("h_neventsinfile");
    if(h_c!=NULL) h_c->SetDirectory(0);

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
      //if(nLlep()<2)              continue;
      //if(nTlepSS()<1&&nTlep()<2) continue;//preselection can be done already here
      //if(nb()!=0)                continue;//preselection can be done already here

      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_") !=string::npos && gen_ht()>100.) continue;
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos && gen_ht()>100.) continue;
      //if(string(currentFile->GetTitle()).find("www_2l_mia")        !=string::npos) weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")   !=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(isData()) weight = 1.;
      //double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      if(!isData()&&applyPUrewgt)    weight *= purewgt();
      if(!isData()&&applylepSF)      weight *= lepsf();
      if(!isData()&&applytrigSF)     weight *= trigsf();
      bool checkevent = false;
            
      if(isData()){
        if(!passFilters())                      continue;
        if(checkevent) cout << "pass filter"    << endl;
        duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
        if( is_duplicate(id) )                  continue; 
        if(checkevent) cout << "pass duplicate" << endl;
        if( !goodrun(tas::run(), tas::lumi()) ) continue;
        if(checkevent) cout << "pass goodrun"   << endl;
        weight = 1.;
      } 
      //if(!passTriggers(true,true)) continue;//pass trigger for data, and offline lepton kinematic cuts for data/simulation
      if(checkevent) cout << "pass online/offline triggers" << endl;
      
      string sample   = skimFilePrefix;
      if(splitVH(fname)){ sample = "WHtoWWW"; }
      string sn = string(bkgtype().Data());
      sn = sample;
      if(vetophoton()) continue;
      
      int SRSS[30]; 
      int SR3l[30];
      for(int i = 0; i<30; ++i) { SRSS[i] = -1; SR3l[i] = -1;  }

      //SS
      passAnySS(SRSS[1],SRSS[3],SRSS[5],true); //preselection: 1: SR, 3: AR, 5: CR
      //3l
      passAny3l(SR3l[1],SR3l[3],SR3l[5],true); //preselection: 1: SR, 3: AR, 5: CR

      //LorentzVector METlv; METlv.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());

      float ST   = 0;
      float STall = 0;
      ST    += met_pt();
      STall += met_pt();
      if(SRSS[1]>=0) {
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(         lep_pass_VVV_cutbased_tight   ()[i]&&lep_p4()[i].Pt()>25.)    ST += lep_p4()[i].Pt();
        }
      }
      if(SR3l[1]>=0) {
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(         lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>20.)    ST += lep_p4()[i].Pt();
        }
      }
      for(unsigned int i = 0; i<lep_pdgId().size();++i){
        if(           lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>20.) STall += lep_p4()[i].Pt();
      }
      for(unsigned int n = 0; n<jets_csv().size();++n){
        if(                    jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5)    STall += jets_p4()[n].Pt();
        if(  SR3l[1]>=0 &&     jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5)    ST    += jets_p4()[n].Pt();
        if(  SRSS[1]>=0 &&     jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)  ST    += jets_p4()[n].Pt();
        if(   jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=2.4&&   jets_csv()[n]>=0.5426) {
          if(                !(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5))   STall += jets_p4()[n].Pt();
          if(SR3l[1]>=0 &&   !(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5))   ST    += jets_p4()[n].Pt();
          if(SRSS[1]>=0 &&   !(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)) ST    += jets_p4()[n].Pt();
        }
      }
      vector<float> mWLC; mWLC.clear();
      vector<float> mWFC; mWFC.clear();
      for(unsigned i = 0; i<genPart_pdgId().size();++i){
        if(abs(genPart_pdgId()[i])!=24) continue;
        if(genPart_status()[i]==22) mWLC.push_back(genPart_p4()[i].M());
        if(genPart_status()[i]==62) mWFC.push_back(genPart_p4()[i].M());
      }

      vector<float> mWfromDecayEMU; mWLC.clear();
      if(SRSS[1]==1){
        LorentzVector e(0,0,0,0),mu(0,0,0,0),enu(0,0,0,0),munu(0,0,0,0),q(0,0,0,0),qp(0,0,0,0);
        int idq(0),idqp(0);
        for(unsigned i = 0; i<genPart_pdgId().size();++i){
          if(abs(genPart_motherId()[i])!=24) continue;
          if(abs(genPart_pdgId()[i])==11) e    = genPart_p4()[i];
          if(abs(genPart_pdgId()[i])==12) enu  = genPart_p4()[i];
          if(abs(genPart_pdgId()[i])==13) mu   = genPart_p4()[i];
          if(abs(genPart_pdgId()[i])==14) munu = genPart_p4()[i];
          if(abs(genPart_pdgId()[i])<=4&&genPart_pdgId()[i]!=0) {
            if(idq==0) { q = genPart_p4()[i]; idq = genPart_pdgId()[i]; }
            else if(idqp==0) { qp = genPart_p4()[i]; idqp = genPart_pdgId()[i]; }
            else cout << __LINE__ << " " << idq << " " << idqp <<  " " << genPart_pdgId()[i] << endl;
          }
        }
        if( e.Pt()>0&& enu.Pt()>0) mWfromDecayEMU.push_back((e+enu  ).M());
        if(mu.Pt()>0&&munu.Pt()>0) mWfromDecayEMU.push_back((mu+munu).M());
        if( q.Pt()>0&&  qp.Pt()>0) mWfromDecayEMU.push_back((q+qp   ).M());
      }

      //cout << sample << " " << sn << endl;
      if(!isData()||!blindSR){//SR is blinded
        if(passTriggers(true,true)){
          if(SRSS[1]>=0){
            for(unsigned int j = 0; j<mWLC.size();++j) fillhisto(histos, "WgenwidthLC_SRSSpresel", sample, sn, mWLC[j],    weight);
            for(unsigned int j = 0; j<mWFC.size();++j) fillhisto(histos, "WgenwidthFC_SRSSpresel", sample, sn, mWFC[j],    weight);
            for(unsigned int j = 0; j<mWfromDecayEMU.size();++j) fillhisto(histos, "Wwidthfromdecay_SRSSpresel", sample, sn, mWfromDecayEMU[j],    weight);
            fillhisto(histos,  "Mjj_SRSSpresel",  sample, sn,  Mjj(), weight);
            fillhisto(histos, "MjjL_SRSSpresel",  sample, sn, MjjL(), weight);
            fillhisto(histos,   "ST_SRSSpresel",  sample, sn, ST,     weight);
          }
          if(SR3l[1]>=0){
            for(unsigned int j = 0; j<mWLC.size();++j) fillhisto(histos, "WgenwidthLC_SR3lpresel", sample, sn, mWLC[j],    weight);
            for(unsigned int j = 0; j<mWFC.size();++j) fillhisto(histos, "WgenwidthFC_SR3lpresel", sample, sn, mWFC[j],    weight);
            fillhisto(histos,  "Mjj_SR3lpresel",  sample, sn,  Mjj(), weight);
            fillhisto(histos, "MjjL_SR3lpresel",  sample, sn, MjjL(), weight);
            fillhisto(histos,   "ST_SR3lpresel",  sample, sn, ST,     weight);
          }
        }
      }
      for(unsigned int j = 0; j<mWLC.size();++j) fillhisto(histos, "WgenwidthLC_all", sample, sn, mWLC[j],    weight);
      for(unsigned int j = 0; j<mWFC.size();++j) fillhisto(histos, "WgenwidthFC_all", sample, sn, mWFC[j],    weight);
      if(Mjj()>0)  fillhisto(histos,  "Mjj_all",  sample, sn,  Mjj(), weight);
      if(MjjL()>0) fillhisto(histos, "MjjL_all",  sample, sn, MjjL(), weight);
      fillhisto(             histos,   "ST_all",  sample, sn, STall,  weight);
      fillhisto(             histos, "NWLC_all", sample, sn, mWLC.size(),    weight);
      fillhisto(             histos, "NWFC_all", sample, sn, mWFC.size(),    weight);
      if(checkevent) cout << endl;
      
    }//event loop
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }


  SaveHistosToFile("rootfiles/Wwidth.root",histos,true,true,(chainnumber==0));
  deleteHistograms(histos);
  
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.02f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.02f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
