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
//#include "Functions112.h"
//#include "CMS3_WWW112.cc"
#include "Functions.h"
#include "CMS3_WWW121.cc"
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1, int year=2016) {
  
  int    counterSS(0), counter3l(0);
  int    cnt1(0), cnt2(0);
  double  weightSS(0),  weight3l(0);

  bool blindSR         = false;
  bool btagreweighting = true;
  bool applylepSF      = true;
  bool applytrigSF     = true;
  bool applyPUrewgt    = true;
  
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  const char* json_file = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int> hbins; hbins.clear();
  vector<float> hlow; hlow.clear();
  vector<float> hup; hup.clear();

  //split ttX/ttV/ttZ
  histonames.push_back("NJ_3lCRpreselect");          hbins.push_back(10);  hlow.push_back( 0); hup.push_back(10);
  histonames.push_back("ST_3lCRSSpreselect");        hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("STmod_3lCRSSpreselect");     hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_3lCR3lpreselect");        hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("STmod_3lCR3lpreselect");     hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL");            hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_preselect");         hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL");            hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_preselect");         hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("SRfull");           hbins.push_back(9);  hlow.push_back( 0); hup.push_back(9);
  histonames.push_back("SRfullsb");           hbins.push_back(9);  hlow.push_back( 0); hup.push_back(9);
  histonames.push_back("SRpresel");         hbins.push_back(9);  hlow.push_back( 0); hup.push_back(9);

  histonames.push_back("ST_SR3l_noMjjL_JESup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_JESdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_JESup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_JESdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_PDFup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_PDFdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_PDFup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_PDFdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_Qsqup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_Qsqdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_Qsqup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_Qsqdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);

  
  histonames.push_back("ST_SR3l_noMjjL_PUup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_PUdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_PUup");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_PUdn");      hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_lepSFup");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_lepSFdn");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_lepSFup");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_lepSFdn");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_trgSFup");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_trgSFdn");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_trgSFup");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_trgSFdn");   hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_bLSFup");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_bLSFdn");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_bLSFup");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_bLSFdn");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_bHSFup");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SR3l_noMjjL_bHSFdn");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_bHSFup");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);
  histonames.push_back("ST_SRSS_noMjjL_bHSFdn");    hbins.push_back(20);  hlow.push_back( 0); hup.push_back(5000);

  map<string, TH1D*> histos =  bookhistograms(skimFilePrefix, histonames,hbins, hlow, hup, rootdir,2);
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
    h_c->SetDirectory(0);

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
      if(nLlep()<2)              continue;
      if(nTlepSS()<1&&nTlep()<2) continue;//preselection can be done already here

      //weight = 1;
      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos){
	if(gen_ht()>100) continue;
      }
      //if(string(currentFile->GetTitle()).find("www_2l_mia")!=string::npos)      weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")!=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      if(weight>100) cout << weight << " " << currentFile->GetTitle() << endl;
      if(isData()) weight = 1.;
      double rawweight = weight;
      if(!isData()&&btagreweighting) weight *= weight_btagsf();
      if(!isData()&&applyPUrewgt)    weight *= purewgt();
      if(!isData()&&applylepSF)      weight *= lepsf();
      if(!isData()&&applytrigSF)     weight *= trigsf();

      if(isData()){
        if(!passFilters())                      continue;
        duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
        if( is_duplicate(id) )                  continue; 
        if( !goodrun(tas::run(), tas::lumi()) ) continue;
      } 
      if(!passTriggers()) continue;//pass trigger for data, and offline lepton kinematic cuts for data/simulation
      
      string sample   = skimFilePrefix;
      if(splitVH(fname)){ sample = "WHtoWWW"; }
      string sn = string(bkgtype().Data());
      if(vetophoton()) continue;
      
      string extrasn = "";
      bool passextra1 = false;
      bool passextra2 = false;
      if(fname.find("wpwpjj_ewk")!=string::npos) { passextra1 = true; extrasn  = "WWVBS"; }
      //if(fname.find("ttw_"      )!=string::npos) { passextra2 = true; extrasn  = "ttW";   }
      //if(fname.find("ttz_"      )!=string::npos) { passextra2 = true; extrasn  = "ttZ";   }
      //if(fname.find("vbswz")!=string::npos)      { passextra1 = true; extrasn  = "ttZ"; }//this is stupid but should work

      int SRSS = isSRSS(true);
      int SRSSju = isSRSS(true, 1.);
      int SRSSjd = isSRSS(true,-1.);
      int SR3l = isSR3l(true);
      int SRSSf = isSRSS(false);
      int SRSSfsb = isSRSS(false,0,false,true);
      int SR3lf = isSR3l(false);
      int SR3lfju = isSR3l(false, 1.);
      int SR3lfjd = isSR3l(false,-1.);
      int SRSSfull = SRSS;
      int SRSSfullju = SRSSju;
      int SRSSfulljd = SRSSjd;
      if(fabs(DetajjL())>1.5/*&&fabs(Mjj()-80.)>=15.*/) SRSSfull = -1;
      //if(SRSSfull>=0 && fabs(Mjj()-80.)>=15.) ++cnt1;
      //if(SRSSfull>=0 && fabs(Mjj()-80.)< 15.) ++cnt2;
      if(SRSSfull==0 && !(fabs(MeeSS()-MZ)>10.&&met_pt()>60.&&MllSS()>40.)) SRSSfull = -1;
      if(SRSSfull==1 && !(met_pt()>60.&&MllSS()>40.&&MTmax()>90.))          SRSSfull = -1;
      if(SRSSfull==2 && !(met_pt()>60.&&MllSS()>40.))                       SRSSfull = -1;
      if(fabs(DetajjL_up())>1.5/*&&fabs(Mjj_up()-80.)>=15.*/) SRSSfullju = -1;
      if(SRSSfullju==0 && !(fabs(MeeSS()-MZ)>10.&&met_up_pt()>60.&&MllSS()>40.)) SRSSfullju = -1;
      if(SRSSfullju==1 && !(met_up_pt()>60.&&MllSS()>40.&&MTmax_up()>90.))       SRSSfullju = -1;
      if(SRSSfullju==2 && !(met_up_pt()>60.&&MllSS()>40.))                       SRSSfullju = -1;
      if(fabs(DetajjL_dn())>1.5/*&&fabs(Mjj_dn()-80.)>=15.*/) SRSSfulljd = -1;
      if(SRSSfulljd==0 && !(fabs(MeeSS()-MZ)>10.&&met_dn_pt()>60.&&MllSS()>40.)) SRSSfulljd = -1;
      if(SRSSfulljd==1 && !(met_dn_pt()>60.&&MllSS()>40.&&MTmax_dn()>90.))       SRSSfulljd = -1;
      if(SRSSfulljd==2 && !(met_dn_pt()>60.&&MllSS()>40.))                       SRSSfulljd = -1;
      int CRSS = isCRSS(true);
      int CR3l = isCR3l(true);

      float ST = 0;
      float STju = 0;
      float STjd = 0;
      float STmod = 0;
      LorentzVector METlv; METlv.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      ST += met_pt();
      STju += met_up_pt();
      STjd += met_dn_pt();
      STmod = met_pt();

      if(SRSSfullju>=0) {
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_p4()[i].Pt()>25.) STju  += lep_p4()[i].Pt();
        }
      }
      if(SRSSfulljd>=0) {
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_p4()[i].Pt()>25.) STjd  += lep_p4()[i].Pt();
        }
      }
      if(SR3lfju>=0) {
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(         lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>20.)  STju += lep_p4()[i].Pt();
        }
      }
      if(SR3lfjd>=0) {
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(         lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>20.)  STjd += lep_p4()[i].Pt();
        }
      }
          
      if(SRSS>=0) {
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_p4()[i].Pt()>25.) ST    += lep_p4()[i].Pt();
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_p4()[i].Pt()>25.) STmod += lep_p4()[i].Pt();
        }
      }
      else if(SR3l>=0||CRSS>=0||CR3l>=0) {
        int nlepcount=0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(         lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>20.)    ST += lep_p4()[i].Pt();
          if(CRSS<=0&&lep_pass_VVV_cutbased_3l_tight()[i]&&lep_p4()[i].Pt()>20.) STmod += lep_p4()[i].Pt();
        }
        if(CRSS>=0){
          int lep3 = -1; int SS1 = lep_idx0_SS(); int SS2 = lep_idx1_SS();
          if(SS1<0||SS2<0) cout << "WTF " << __LINE__ << " " << SS1 << " " << SS2 << endl;
          if(SS1==0&&SS2==1) lep3 = 2;
          else if(SS1==0&&SS2==2) lep3 = 1;
          else if(SS1==1&&SS2==2) lep3 = 0;
          //if(SS1 >=0) STmod += lep_p4()[SS1].Pt();
          //if(SS2 >=0) STmod += lep_p4()[SS2].Pt();
          if(lep3>=0) STmod += (METlv+lep_p4()[lep3]).Pt();
          if(lep3>=0) STmod -= met_pt();
          if(lep3>=0) STmod -= (lep_p4()[lep3]).Pt();
        }
      }
      for(unsigned int n = 0; n<jets_up_csv().size();++n){
        if((SR3lfju>=0)    &&   jets_up_p4()[n].Pt()>=30&&fabs(   jets_up_p4()[n].Eta())<=5)        STju += jets_up_p4()[n].Pt();
        if((SRSSfullju>=0) &&   jets_up_p4()[n].Pt()>=30&&fabs(   jets_up_p4()[n].Eta())<=2.5)      STju += jets_up_p4()[n].Pt();
        if(   jets_up_p4()[n].Pt()>=20&&fabs(   jets_up_p4()[n].Eta())<=2.4&&   jets_up_csv()[n]>=0.5426) {
          if((SR3lfju>=0)    &&   !(jets_up_p4()[n].Pt()>=30&&fabs(   jets_up_p4()[n].Eta())<=5))   STju += jets_up_p4()[n].Pt();
          if((SRSSfullju>=0) &&   !(jets_up_p4()[n].Pt()>=30&&fabs(   jets_up_p4()[n].Eta())<=2.5)) STju += jets_up_p4()[n].Pt();
        }
      }
      for(unsigned int n = 0; n<jets_dn_csv().size();++n){
        if((SR3lfjd>=0)    &&   jets_dn_p4()[n].Pt()>=30&&fabs(   jets_dn_p4()[n].Eta())<=5)        STjd += jets_dn_p4()[n].Pt();
        if((SRSSfulljd>=0) &&   jets_dn_p4()[n].Pt()>=30&&fabs(   jets_dn_p4()[n].Eta())<=2.5)      STjd += jets_dn_p4()[n].Pt();
        if(   jets_dn_p4()[n].Pt()>=20&&fabs(   jets_dn_p4()[n].Eta())<=2.4&&   jets_dn_csv()[n]>=0.5426) {
          if((SR3lfjd>=0)    &&   !(jets_dn_p4()[n].Pt()>=30&&fabs(   jets_dn_p4()[n].Eta())<=5))   STjd += jets_dn_p4()[n].Pt();
          if((SRSSfulljd>=0) &&   !(jets_dn_p4()[n].Pt()>=30&&fabs(   jets_dn_p4()[n].Eta())<=2.5)) STjd += jets_dn_p4()[n].Pt();
        }
      }
      
      for(unsigned int n = 0; n<jets_csv().size();++n){
        if((SR3l>=0||CR3l>=0) &&   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5)        ST += jets_p4()[n].Pt();
        if((SRSS>=0||CRSS>=0) &&   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)      ST += jets_p4()[n].Pt();
        if(   jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=2.4&&   jets_csv()[n]>=0.5426) {
          if((SR3l>=0||CR3l>=0) &&   !(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5))   ST += jets_p4()[n].Pt();
          if((SRSS>=0||CRSS>=0) &&   !(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)) ST += jets_p4()[n].Pt();
        }
        if((SR3l>=0||CR3l>=0) &&   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5)        STmod += jets_p4()[n].Pt();
        if((SRSS>=0||CRSS>=0) &&   jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)      STmod += jets_p4()[n].Pt();
        if(   jets_p4()[n].Pt()>=20&&fabs(   jets_p4()[n].Eta())<=2.4&&   jets_csv()[n]>=0.5426) {
          if((SR3l>=0||CR3l>=0) &&   !(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=5))   STmod += jets_p4()[n].Pt();
          if((SRSS>=0||CRSS>=0) &&   !(jets_p4()[n].Pt()>=30&&fabs(   jets_p4()[n].Eta())<=2.5)) STmod += jets_p4()[n].Pt();
        }
      }
      bool fillboolean = true; string spln = sn;
      //if(passextra1){ fillboolean = false; spln = extrasn; } // uncomment this for first study

      if(CRSS>=0||CR3l>=0) fillhisto(histos, "NJ_3lCRpreselect",         sample, spln, nj30(),  weight, fillboolean);
      if(CRSS>=0) fillhisto(histos, "ST_3lCRSSpreselect",         sample, spln, ST,    weight, fillboolean);
      if(CRSS>=0) fillhisto(histos, "STmod_3lCRSSpreselect",      sample, spln, STmod, weight, fillboolean);
      if(CR3l>=0) fillhisto(histos, "ST_3lCR3lpreselect",         sample, spln, ST,    weight, fillboolean);
      if(CRSS>=0 && ST>=2000.) { ++counterSS; weightSS += weight; }
      if(CR3l>=0 && ST>=1500.) { ++counter3l; weight3l += weight; }
      if(!(blindSR&&isData())){
        if(SRSS>=0)     fillhisto(histos, "ST_SRSS_preselect",      sample, spln, ST,    weight, fillboolean);
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL",         sample, spln, ST,    weight, fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL",         sample, spln, ST,    weight, fillboolean);
        if(SR3l>=0)     fillhisto(histos, "ST_SR3l_preselect",      sample, spln, ST,    weight, fillboolean);
        fillSRhisto(histos, "SRpresel",         sample, sn, SRSS,  SR3l,  weight);
        fillSRhisto(histos, "SRfull",           sample, sn, SRSSf, SR3lf, weight);
        fillSRhisto(histos, "SRfullsb",           sample, sn, SRSSfsb, SR3lf, weight);

        
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_PDFup",         sample, spln, ST,    weight*weight_pdf_up()/weight_fr_r1_f1()      , fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_PDFup",         sample, spln, ST,    weight*weight_pdf_up()/weight_fr_r1_f1()      , fillboolean);
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_PDFdn",         sample, spln, ST,     weight*weight_pdf_down()/weight_fr_r1_f1()   , fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_PDFdn",         sample, spln, ST,     weight*weight_pdf_down()/weight_fr_r1_f1()   , fillboolean);
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_Qsqup",         sample, spln, ST,    weight*weight_fr_r2_f2()/weight_fr_r1_f1()    , fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_Qsqup",         sample, spln, ST,    weight*weight_fr_r2_f2()/weight_fr_r1_f1()    , fillboolean);
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_Qsqdn",         sample, spln, ST,    weight*weight_fr_r0p5_f0p5()/weight_fr_r1_f1(), fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_Qsqdn",         sample, spln, ST,    weight*weight_fr_r0p5_f0p5()/weight_fr_r1_f1(), fillboolean);
        /*
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_PDFup",         sample, spln, ST,    weight*weight_pdf_up()/weight_fr_r1_f1()      *h_c->GetBinContent(2)/h_c->GetBinContent(11), fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_PDFup",         sample, spln, ST,    weight*weight_pdf_up()/weight_fr_r1_f1()      *h_c->GetBinContent(2)/h_c->GetBinContent(11), fillboolean);
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_PDFdn",         sample, spln, ST,     weight*weight_pdf_down()/weight_fr_r1_f1()   *h_c->GetBinContent(2)/h_c->GetBinContent(12), fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_PDFdn",         sample, spln, ST,     weight*weight_pdf_down()/weight_fr_r1_f1()   *h_c->GetBinContent(2)/h_c->GetBinContent(12), fillboolean);
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_Qsqup",         sample, spln, ST,    weight*weight_fr_r2_f2()/weight_fr_r1_f1()    *h_c->GetBinContent(2)/h_c->GetBinContent(6), fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_Qsqup",         sample, spln, ST,    weight*weight_fr_r2_f2()/weight_fr_r1_f1()    *h_c->GetBinContent(2)/h_c->GetBinContent(6), fillboolean);
        if(SRSSfull>=0) fillhisto(histos, "ST_SRSS_noMjjL_Qsqdn",         sample, spln, ST,    weight*weight_fr_r0p5_f0p5()/weight_fr_r1_f1()*h_c->GetBinContent(2)/h_c->GetBinContent(10), fillboolean);
        if(SR3lf>=0)    fillhisto(histos, "ST_SR3l_noMjjL_Qsqdn",         sample, spln, ST,    weight*weight_fr_r0p5_f0p5()/weight_fr_r1_f1()*h_c->GetBinContent(2)/h_c->GetBinContent(10), fillboolean);  
        */  
        if(SRSSfullju>=0) fillhisto(histos, "ST_SRSS_noMjjL_JESup",         sample, spln, STju,    weight, fillboolean);
        if(SR3lfju>=0)    fillhisto(histos, "ST_SR3l_noMjjL_JESup",         sample, spln, STju,    weight, fillboolean);
        if(SRSSfulljd>=0) fillhisto(histos, "ST_SRSS_noMjjL_JESdn",         sample, spln, STjd,    weight, fillboolean);
        if(SR3lfjd>=0)    fillhisto(histos, "ST_SR3l_noMjjL_JESdn",         sample, spln, STjd,    weight, fillboolean);
        if(SRSSfull>=0&&trigsf()       !=0) fillhisto(histos, "ST_SRSS_noMjjL_trgSFup",         sample, spln, ST,    weight*trigsf_up()/trigsf(), fillboolean);
        if(SR3lf>=0   &&trigsf()       !=0) fillhisto(histos, "ST_SR3l_noMjjL_trgSFup",         sample, spln, ST,    weight*trigsf_up()/trigsf(), fillboolean);
        if(SRSSfull>=0&&trigsf()       !=0) fillhisto(histos, "ST_SRSS_noMjjL_trgSFdn",         sample, spln, ST,    weight*trigsf_dn()/trigsf(), fillboolean);
        if(SR3lf>=0   &&trigsf()       !=0) fillhisto(histos, "ST_SR3l_noMjjL_trgSFdn",         sample, spln, ST,    weight*trigsf_dn()/trigsf(), fillboolean);
        if(SRSSfull>=0&&purewgt()      !=0) fillhisto(histos, "ST_SRSS_noMjjL_PUup",            sample, spln, ST,    weight*purewgt_up()/purewgt(), fillboolean);
        if(SR3lf>=0   &&purewgt()      !=0) fillhisto(histos, "ST_SR3l_noMjjL_PUup",            sample, spln, ST,    weight*purewgt_up()/purewgt(), fillboolean);
        if(SRSSfull>=0&&purewgt()      !=0) fillhisto(histos, "ST_SRSS_noMjjL_PUdn",            sample, spln, ST,    weight*purewgt_up()/purewgt(), fillboolean);
        if(SR3lf>=0   &&purewgt()      !=0) fillhisto(histos, "ST_SR3l_noMjjL_PUdn",            sample, spln, ST,    weight*purewgt_up()/purewgt(), fillboolean);
        if(SRSSfull>=0&&lepsf()        !=0) fillhisto(histos, "ST_SRSS_noMjjL_lepSFup",         sample, spln, ST,    weight*lepsf_up()/lepsf(), fillboolean);
        if(SR3lf>=0   &&lepsf()        !=0) fillhisto(histos, "ST_SR3l_noMjjL_lepSFup",         sample, spln, ST,    weight*lepsf_up()/lepsf(), fillboolean);
        if(SRSSfull>=0&&lepsf()        !=0) fillhisto(histos, "ST_SRSS_noMjjL_lepSFdn",         sample, spln, ST,    weight*lepsf_dn()/lepsf(), fillboolean);
        if(SR3lf>=0   &&lepsf()        !=0) fillhisto(histos, "ST_SR3l_noMjjL_lepSFdn",         sample, spln, ST,    weight*lepsf_dn()/lepsf(), fillboolean);
        if(SRSSfull>=0&&weight_btagsf()!=0) fillhisto(histos, "ST_SRSS_noMjjL_bLSFup",          sample, spln, ST,    weight*weight_btagsf_light_UP()/weight_btagsf(), fillboolean);
        if(SR3lf>=0   &&weight_btagsf()!=0) fillhisto(histos, "ST_SR3l_noMjjL_bLSFup",          sample, spln, ST,    weight*weight_btagsf_light_UP()/weight_btagsf(), fillboolean);
        if(SRSSfull>=0&&weight_btagsf()!=0) fillhisto(histos, "ST_SRSS_noMjjL_bLSFdn",          sample, spln, ST,    weight*weight_btagsf_light_DN()/weight_btagsf(), fillboolean);
        if(SR3lf>=0   &&weight_btagsf()!=0) fillhisto(histos, "ST_SR3l_noMjjL_bLSFdn",          sample, spln, ST,    weight*weight_btagsf_light_DN()/weight_btagsf(), fillboolean);
        if(SRSSfull>=0&&weight_btagsf()!=0) fillhisto(histos, "ST_SRSS_noMjjL_bHSFup",          sample, spln, ST,    weight*weight_btagsf_heavy_UP()/weight_btagsf(), fillboolean);
        if(SR3lf>=0   &&weight_btagsf()!=0) fillhisto(histos, "ST_SR3l_noMjjL_bHSFup",          sample, spln, ST,    weight*weight_btagsf_heavy_UP()/weight_btagsf(), fillboolean);
        if(SRSSfull>=0&&weight_btagsf()!=0) fillhisto(histos, "ST_SRSS_noMjjL_bHSFdn",          sample, spln, ST,    weight*weight_btagsf_heavy_DN()/weight_btagsf(), fillboolean);
        if(SR3lf>=0   &&weight_btagsf()!=0) fillhisto(histos, "ST_SR3l_noMjjL_bHSFdn",          sample, spln, ST,    weight*weight_btagsf_heavy_DN()/weight_btagsf(), fillboolean);
      }
      if(isData()&&(CRSS>=0 || CR3l>=0)&&(ST>1500||STmod>1500)) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " " << ST << " " << STmod << endl;
      if(isData()&&(SRSS>=0 || SR3l>=0)&&(ST>1500)) cout << tas::run() << ":" << tas::lumi() << ":" << tas::evt() << " - " << ST << " " << STmod << endl;
    }//event loop  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  cout << "For " << skimFilePrefix << " have " << counterSS << " (" << counter3l << ") events in SS (3l) CR with total weight sum " << weightSS << " (" << weight3l << ")" << endl;
  cout << "cnt1 " << cnt1 << " cnt2 " << cnt2 << endl;
  
  SaveHistosToFile("rootfiles/CheckBGaQGC.root",histos,true,true,(chainnumber==0));
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
