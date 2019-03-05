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

//check code https://github.com/haweber/Tribosons/blob/master/2017/GetOSYieldsSmeared.C
#include "EnergyScaleCorrection_class.cc"
#include "EnergyScaleCorrection_class.h"

#include "RoccoR.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int chainnumber=1) {

  bool blindSR         = true;
  bool btagreweighting = true;
  bool applylepSF      = true;
  bool applytrigSF     = true;
  bool applyPUrewgt    = true;
  
  const char* json_file = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
  set_goodrun_file_json(json_file);

  
  EnergyScaleCorrection_class eScaler("data/Moriond17_23Jan_ele");
  eScaler.doScale     = true;
  eScaler.doSmearings = true;
  TRandom3 *rgen_ = new TRandom3();
  rgen_->SetSeed(123456);

  RoccoR  rc("rcdata.2016.v3"); //directory path as input for now; initialize only once, contains all variations


  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  vector<string> histonames; histonames.clear();
  vector<int>   hbins;       hbins.clear();
  vector<float> hlow;        hlow.clear();
  vector<float> hup;         hup.clear();

  histonames.push_back("AllSR_xcheck");           hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_raw");              hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_allsmearedscaled"); hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_muosmearedscaled"); hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_elesmearedscaled"); hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_allsmearedscaled_eleup"); hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_allsmearedscaled_eledn"); hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_allsmearedscaled_muoup"); hbins.push_back(9); hlow.push_back(0); hup.push_back(9);
  histonames.push_back("AllSR_allsmearedscaled_muodn"); hbins.push_back(9); hlow.push_back(0); hup.push_back(9);

  histonames.push_back("Mee_raw");                hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedandscaled");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_raw");              hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedandscaled"); hbins.push_back(20); hlow.push_back(80); hup.push_back(100);

  histonames.push_back("Mmumu_raw");                hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mmumu_smearedandscaled");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MmumuSS_raw");              hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MmumuSS_smearedandscaled"); hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  
  histonames.push_back("Mee_smearedandscaled_scaledup");    hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedandscaled_scaleddn");    hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedandscaled_smearup");     hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("Mee_smearedandscaled_smeardn");     hbins.push_back(20); hlow.push_back(80); hup.push_back(100);

  histonames.push_back("MeeSS_smearedandscaled_scaledup");  hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedandscaled_scaleddn");  hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedandscaled_smearup");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);
  histonames.push_back("MeeSS_smearedandscaled_smeardn");   hbins.push_back(20); hlow.push_back(80); hup.push_back(100);

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

      if(nVert()<0)              continue;
      if(firstgoodvertex()!=0)   continue;
      if(nLlep()<2)              continue;
      if(nTlepSS()<2)            continue;//preselection can be done already here
      if(nb()>1 )                continue;//preselection can be done already here
      //if(nj30()<1)               continue;//preselection can be done already here

      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_") !=string::npos && gen_ht()>100.) continue;
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos && gen_ht()>100.) continue;
      //if(string(currentFile->GetTitle()).find("www_2l_mia")        !=string::npos) weight *= 0.066805* 91900./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
      //if(string(currentFile->GetTitle()).find("www_2l_ext1_mia")   !=string::npos) weight *= 0.066805*164800./(91900.+164800.);//(208fb/1pb)*BR(WWW—> >=2l)*combineweight
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
      if(splitVH(fname)) sample = "WHtoWWW"; 
      string sn = string(bkgtype().Data());
      if(vetophoton()) continue;

      LorentzVector MET; MET.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());
      vector<int> lep_id;
      vector<int> lep_q;
      vector<LorentzVector> lep_raw;
      vector<LorentzVector> lep_mod;//muon + electron
      vector<LorentzVector> lep_modmu;//muon only, documentation on variations is poor
      vector<LorentzVector> lep_modmuup;//muon only, documentation on variations is poor
      vector<LorentzVector> lep_modmudn;//muon only, documentation on variations is poor
      //below: all electrons
      vector<LorentzVector> lep_model;//electron only, documentation on variations is quite ok
      vector<LorentzVector> lep_modelup;//electron only, documentation on variations is quite ok
      vector<LorentzVector> lep_modeldn;//electron only, documentation on variations is quite ok
      vector<LorentzVector> lep_scupel;
      vector<LorentzVector> lep_scdnel;
      vector<LorentzVector> lep_smupel;
      vector<LorentzVector> lep_smdnel;
      vector<LorentzVector> lep_scsmel;
      bool isee(false), ismm(false);
      bool isss(false);

      for(unsigned int i = 0; i<lep_pdgId().size(); ++i){
        //if(lep_p4()[i].Pt()<20) continue;
        //if(fabs(lep_p4()[i].Eta())>2.5) continue;
        //if(!lep_pass_VVV_cutbased_tight()[i]) continue;
        lep_id.push_back(lep_pdgId()[i]);
        lep_q.push_back(lep_charge()[i]);
        lep_raw.push_back(lep_p4()[i]);
        if(abs(lep_pdgId()[i])==13){
          lep_model.push_back(lep_p4()[i]);
          lep_scupel.push_back(lep_p4()[i]);
          lep_scdnel.push_back(lep_p4()[i]);
          lep_smupel.push_back(lep_p4()[i]);
          lep_smdnel.push_back(lep_p4()[i]);
          lep_modelup.push_back(lep_p4()[i]);
          lep_modeldn.push_back(lep_p4()[i]);
          if(isData()){
            double dataSF      = rc.kScaleDT(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), 0, 0);
            double dataSFerrSC = rc.kScaleDT(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), 0, 1);
            //cout << __LINE__ << " " << dataSF << " " << dataSFerrSC << endl;
            LorentzVector lv = lep_p4()[i];
            lv.SetPxPyPzE(dataSF*lep_p4()[i].Px(),dataSF*lep_p4()[i].Py(),dataSF*lep_p4()[i].Pz(), sqrt(pow(0.114,2)+pow(dataSF*lep_p4()[i].P(),2)) );
            lep_mod.push_back(lv);
            lep_modmu.push_back(lv);
            lep_modmuup.push_back(lv);
            lep_modmudn.push_back(lv);
          } else {
            double mcSF = 0;
            double mcSFerrSC = 0;
            int    Nerrors = 0;
            int match = -1;
            for(unsigned int j = 0; j<genPart_p4().size(); ++j){
              if(genPart_status()[j]!=23 && genPart_status()[j]!=1) continue;
              if(abs( genPart_pdgId()[j])!=13) continue;
              if(dR(genPart_p4()[j],lep_p4()[i])>0.05) continue;
              match = j;
              break;
            }
            if(match>=0){
              double u1 = gRandom->Rndm();
              mcSF      = rc.kScaleFromGenMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], genPart_p4()[match].Pt(), u1, 0, 0);
              double sum = 0;
              for(unsigned int k = 0; k<100; ++k)
                sum += rc.kScaleFromGenMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], genPart_p4()[match].Pt(), u1, 1, k)/100;
              mcSFerrSC += pow(sum-mcSF,2);
              mcSFerrSC += pow(rc.kScaleFromGenMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], genPart_p4()[match].Pt(), u1, 2, 0)-mcSF,2);
              sum = 0;
              for(unsigned int k = 0; k<5; ++k)
                sum += rc.kScaleFromGenMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], genPart_p4()[match].Pt(), u1, 4, k)/5;
              mcSFerrSC += pow(sum-mcSF,2);
              sum = 0;
              for(unsigned int k = 0; k<5; ++k)
                sum += rc.kScaleFromGenMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], genPart_p4()[match].Pt(), u1, 5, k)/5;
              mcSFerrSC += pow(sum-mcSF,2);
            } else {
              double u1 = gRandom->Rndm();
              double u2 = gRandom->Rndm();
              mcSF      = rc.kScaleAndSmearMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], u1, u2, 0, 0);
              double sum = 0;
              for(unsigned int k = 0; k<100; ++k)
                sum += rc.kScaleAndSmearMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], u1, u2, 1, k)/100;
              mcSFerrSC += pow(sum-mcSF,2);
              mcSFerrSC += pow(rc.kScaleAndSmearMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], u1, u2, 2, 0)-mcSF,2);
              sum = 0;
              for(unsigned int k = 0; k<5; ++k)
                sum += rc.kScaleAndSmearMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], u1, u2, 4, k)/5;
              mcSFerrSC += pow(sum-mcSF,2);
              sum = 0;
              for(unsigned int k = 0; k<5; ++k)
                sum += rc.kScaleAndSmearMC(lep_charge()[i], lep_p4()[i].Pt(), lep_p4()[i].Eta(), lep_p4()[i].Phi(), lep_nlayers()[i], u1, u2, 5, k)/5;
              mcSFerrSC += pow(sum-mcSF,2);
            }
            //cout << __LINE__ << " " << mcSF << " " << mcSFerrSC << endl;
            mcSFerrSC = sqrt(mcSFerrSC);
            LorentzVector lv = lep_p4()[i];
            lv.SetPxPyPzE(mcSF*lep_p4()[i].Px(),mcSF*lep_p4()[i].Py(),mcSF*lep_p4()[i].Pz(), sqrt(pow(0.114,2)+pow(mcSF*lep_p4()[i].P(),2)) );
            lep_mod.push_back(lv);
            lep_modmu.push_back(lv);
            lv.SetPxPyPzE((mcSF+mcSFerrSC)*lep_p4()[i].Px(),(mcSF+mcSFerrSC)*lep_p4()[i].Py(),(mcSF+mcSFerrSC)*lep_p4()[i].Pz(), sqrt(pow(0.114,2)+pow((mcSF+mcSFerrSC)*lep_p4()[i].P(),2)) );
            lep_modmuup.push_back(lv);
            lv.SetPxPyPzE((mcSF-mcSFerrSC)*lep_p4()[i].Px(),(mcSF-mcSFerrSC)*lep_p4()[i].Py(),(mcSF-mcSFerrSC)*lep_p4()[i].Pz(), sqrt(pow(0.114,2)+pow((mcSF-mcSFerrSC)*lep_p4()[i].P(),2)) );
            lep_modmudn.push_back(lv);
          }
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
          scale_corr = eScaler.ScaleCorrection( tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch );//r9 is missing
          if(!isData())   sigma       = eScaler.getSmearingSigma(          tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch, 0, 0);
          if(!isData())   sigmaup     = eScaler.getSmearingSigma(          tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  1, 0);//smear up
          if(!isData())   sigmadown   = eScaler.getSmearingSigma(          tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch, -1, 0);//smear down
          if(!isData()){
            float error_stat = eScaler.ScaleCorrectionUncertainty(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  1);//stat
            float error_syst = eScaler.ScaleCorrectionUncertainty(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  2);//syst
            float error_gain = eScaler.ScaleCorrectionUncertainty(tas::run(), isEBEle, lep_r9()[i], fabs(lep_etaSC()[i]), lep_p4()[i].Et(), gainswitch,  4);//gain
            error_scale = sqrt(error_stat*error_stat+error_syst*error_syst+error_gain*error_gain);
          }
          float smearingfactor = rgen_->Gaus(0.,sigma);
          if(isData()) smearingfactor = 0.;
          lep_modmu.push_back(lep_p4()[i]);
          if(isData()){
            lep_mod.push_back(lep_p4()[i]*(scale_corr));
            lep_model.push_back(lep_p4()[i]*(scale_corr));
            lep_modmuup.push_back(lep_p4()[i]*(scale_corr));
            lep_modmudn.push_back(lep_p4()[i]*(scale_corr));
            lep_scupel.push_back(lep_p4()[i]*(1.+error_scale));
            lep_scdnel.push_back(lep_p4()[i]*(1.-error_scale));
            lep_modelup.push_back(lep_p4()[i]*(1.+error_scale));
            lep_modeldn.push_back(lep_p4()[i]*(1.-error_scale));
            lep_smupel.push_back(lep_p4()[i]);
            lep_smdnel.push_back(lep_p4()[i]);
          }
          else {
            lep_mod.push_back(lep_p4()[i]*(1+smearingfactor));
            lep_model.push_back(lep_p4()[i]*(1+smearingfactor));
            lep_modmuup.push_back(lep_p4()[i]*(1+smearingfactor));
            lep_modmudn.push_back(lep_p4()[i]*(1+smearingfactor));
            lep_scupel.push_back(lep_p4()[i]);
            lep_scdnel.push_back(lep_p4()[i]);
            lep_smupel.push_back(lep_p4()[i]*(1+smearingfactor*sigmaup/sigma));
            lep_smdnel.push_back(lep_p4()[i]*(1-smearingfactor*sigmadown/sigma));
            lep_modelup.push_back(lep_p4()[i]*(1+smearingfactor*sigmaup/sigma));
            lep_modeldn.push_back(lep_p4()[i]*(1-smearingfactor*sigmadown/sigma));
          }
        }
      }

      if(nVlep()<2)  continue;
      if(lep_raw.size()<2) continue;

      LorentzVector METlv; METlv.SetPxPyPzE(met_pt()*TMath::Cos(met_phi()),met_pt()*TMath::Sin(met_phi()),0,met_pt());


      int SRSS[40], SR3l[40];
      int dummy = -1;
      for(int i = 0; i<40; ++i) { SRSS[i] = -1; SR3l[i] = -1; }
      //SS
      //raw - SR 0 /2
      //all   SR 10 / 12
      int nleps  = 0;
      int nleps2 = 0;
      int SRid = -1;
      passAnySS(SRSS[ 0],dummy,SRSS[ 1]);      //full:         0: SR, 1: CR
      passAny3l(SR3l[0], dummy, SR3l[1]);
      if(lep_raw.size()==2&&nVlep()==2&& nisoTrack_mt2_cleaned_VVV_cutbased_veto()==0&&lep_q[0]*lep_q[1]>0){
        passAnySS(SRSS[ 0],dummy,SRSS[ 1]);      //full:         0: SR, 1: CR
        passAnySS(SRSS[ 2],dummy,SRSS[ 3],false, 0,false,false,true,1);//mjj side: 2
        passAnySS(SRSS[ 5],dummy,SRSS[ 1]);      //full:         0: SR, 1: CR
        passAnySS(SRSS[ 6],dummy,SRSS[ 3],false, 0,false,false,true,1);//mjj side: 2
        if(lep_id[0]*lep_id[1]==121)    SRid = 0;
        if(lep_id[0]*lep_id[1]==143)    SRid = 1;
        if(lep_id[0]*lep_id[1]==169)    SRid = 2;
        if(SRid==0 && met_pt()<60.) SRid = -1;
        if(SRid==1 && met_pt()<60.) SRid = -1;
        //raw   SR 0 / 2
        SRSS[0] = -1; SRSS[2] = -1;
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_raw[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_raw[0]+lep_raw[1]).M()-MZ)>10. && (lep_raw[0]+lep_raw[1]).M()>40.) SRSS[0] = 0;
            if(SRid==1 && (lep_raw[0]+lep_raw[1]).M()>30. && TMath::Max(mT(lep_raw[0],METlv),mT(lep_raw[1],METlv))>90.) SRSS[0] = 1;
            if(SRid==2 && (lep_raw[0]+lep_raw[1]).M()>40.) SRSS[0] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_raw[0]+lep_raw[1]).M()-MZ)>10. && (lep_raw[0]+lep_raw[1]).M()>40.) SRSS[2] = 0;
            if(SRid==1 && (lep_raw[0]+lep_raw[1]).M()>30. && TMath::Max(mT(lep_raw[0],METlv),mT(lep_raw[1],METlv))>90.) SRSS[2] = 1;
            if(SRid==2 && (lep_raw[0]+lep_raw[1]).M()>40.&&met_pt()>60.) SRSS[2] = 2;
          }
        }
        //ele+muo   SR 10 / 12
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_mod[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_mod[0]+lep_mod[1]).M()-MZ)>10. && (lep_mod[0]+lep_mod[1]).M()>40.) SRSS[10] = 0;
            if(SRid==1 && (lep_mod[0]+lep_mod[1]).M()>30. && TMath::Max(mT(lep_mod[0],METlv),mT(lep_mod[1],METlv))>90.) SRSS[10] = 1;
            if(SRid==2 && (lep_mod[0]+lep_mod[1]).M()>40.) SRSS[10] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_mod[0]+lep_mod[1]).M()-MZ)>10. && (lep_mod[0]+lep_mod[1]).M()>40.) SRSS[12] = 0;
            if(SRid==1 && (lep_mod[0]+lep_mod[1]).M()>30. && TMath::Max(mT(lep_mod[0],METlv),mT(lep_mod[1],METlv))>90.) SRSS[12] = 1;
            if(SRid==2 && (lep_mod[0]+lep_mod[1]).M()>40.&&met_pt()>60.) SRSS[12] = 2;
          }
        }
        //ele+muo ele up  SR 13 / 14
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_modelup[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_modelup[0]+lep_modelup[1]).M()-MZ)>10. && (lep_modelup[0]+lep_modelup[1]).M()>40.) SRSS[13] = 0;
            if(SRid==1 && (lep_modelup[0]+lep_modelup[1]).M()>30. && TMath::Max(mT(lep_modelup[0],METlv),mT(lep_modelup[1],METlv))>90.) SRSS[13] = 1;
            if(SRid==2 && (lep_modelup[0]+lep_modelup[1]).M()>40.) SRSS[13] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_modelup[0]+lep_modelup[1]).M()-MZ)>10. && (lep_modelup[0]+lep_modelup[1]).M()>40.) SRSS[14] = 0;
            if(SRid==1 && (lep_modelup[0]+lep_modelup[1]).M()>30. && TMath::Max(mT(lep_modelup[0],METlv),mT(lep_modelup[1],METlv))>90.) SRSS[14] = 1;
            if(SRid==2 && (lep_modelup[0]+lep_modelup[1]).M()>40.&&met_pt()>60.) SRSS[14] = 2;
          }
        }
        //ele+muo ele dn  SR 15 / 16
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_modeldn[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_modeldn[0]+lep_modeldn[1]).M()-MZ)>10. && (lep_modeldn[0]+lep_modeldn[1]).M()>40.) SRSS[15] = 0;
            if(SRid==1 && (lep_modeldn[0]+lep_modeldn[1]).M()>30. && TMath::Max(mT(lep_modeldn[0],METlv),mT(lep_modeldn[1],METlv))>90.) SRSS[15] = 1;
            if(SRid==2 && (lep_modeldn[0]+lep_modeldn[1]).M()>40.) SRSS[15] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_modeldn[0]+lep_modeldn[1]).M()-MZ)>10. && (lep_modeldn[0]+lep_modeldn[1]).M()>40.) SRSS[16] = 0;
            if(SRid==1 && (lep_modeldn[0]+lep_modeldn[1]).M()>30. && TMath::Max(mT(lep_modeldn[0],METlv),mT(lep_modeldn[1],METlv))>90.) SRSS[16] = 1;
            if(SRid==2 && (lep_modeldn[0]+lep_modeldn[1]).M()>40.&&met_pt()>60.) SRSS[16] = 2;
          }
        }
        //ele+muo muo up  SR 23 / 24
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_modmuup[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_modmuup[0]+lep_modmuup[1]).M()-MZ)>10. && (lep_modmuup[0]+lep_modmuup[1]).M()>40.) SRSS[23] = 0;
            if(SRid==1 && (lep_modmuup[0]+lep_modmuup[1]).M()>30. && TMath::Max(mT(lep_modmuup[0],METlv),mT(lep_modmuup[1],METlv))>90.) SRSS[23] = 1;
            if(SRid==2 && (lep_modmuup[0]+lep_modmuup[1]).M()>40.) SRSS[23] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_modmuup[0]+lep_modmuup[1]).M()-MZ)>10. && (lep_modmuup[0]+lep_modmuup[1]).M()>40.) SRSS[24] = 0;
            if(SRid==1 && (lep_modmuup[0]+lep_modmuup[1]).M()>30. && TMath::Max(mT(lep_modmuup[0],METlv),mT(lep_modmuup[1],METlv))>90.) SRSS[24] = 1;
            if(SRid==2 && (lep_modmuup[0]+lep_modmuup[1]).M()>40.&&met_pt()>60.) SRSS[24] = 2;
          }
        }
        //ele+muo muo dn  SR 25 / 26
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_modmudn[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_modmudn[0]+lep_modmudn[1]).M()-MZ)>10. && (lep_modmudn[0]+lep_modmudn[1]).M()>40.) SRSS[25] = 0;
            if(SRid==1 && (lep_modmudn[0]+lep_modmudn[1]).M()>30. && TMath::Max(mT(lep_modmudn[0],METlv),mT(lep_modmudn[1],METlv))>90.) SRSS[25] = 1;
            if(SRid==2 && (lep_modmudn[0]+lep_modmudn[1]).M()>40.) SRSS[25] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_modmudn[0]+lep_modmudn[1]).M()-MZ)>10. && (lep_modmudn[0]+lep_modmudn[1]).M()>40.) SRSS[26] = 0;
            if(SRid==1 && (lep_modmudn[0]+lep_modmudn[1]).M()>30. && TMath::Max(mT(lep_modmudn[0],METlv),mT(lep_modmudn[1],METlv))>90.) SRSS[26] = 1;
            if(SRid==2 && (lep_modmudn[0]+lep_modmudn[1]).M()>40.&&met_pt()>60.) SRSS[26] = 2;
          }
        }
        //ele   SR 20 / 22
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_model[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_model[0]+lep_model[1]).M()-MZ)>10. && (lep_model[0]+lep_model[1]).M()>40.) SRSS[20] = 0;
            if(SRid==1 && (lep_model[0]+lep_model[1]).M()>30. && TMath::Max(mT(lep_model[0],METlv),mT(lep_model[1],METlv))>90.) SRSS[20] = 1;
            if(SRid==2 && (lep_model[0]+lep_model[1]).M()>40.) SRSS[20] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_model[0]+lep_model[1]).M()-MZ)>10. && (lep_model[0]+lep_model[1]).M()>40.) SRSS[22] = 0;
            if(SRid==1 && (lep_model[0]+lep_model[1]).M()>30. && TMath::Max(mT(lep_model[0],METlv),mT(lep_model[1],METlv))>90.) SRSS[22] = 1;
            if(SRid==2 && (lep_model[0]+lep_model[1]).M()>40.&&met_pt()>60.) SRSS[22] = 2;
          }
        }
        //muo   SR 30 / 32
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_tight   ()[i]&&lep_modmu[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(passJetSSstate(false, false, 0, false, false,1)){
            if(SRid==0 && fabs((lep_modmu[0]+lep_modmu[1]).M()-MZ)>10. && (lep_modmu[0]+lep_modmu[1]).M()>40.) SRSS[30] = 0;
            if(SRid==1 && (lep_modmu[0]+lep_modmu[1]).M()>30. && TMath::Max(mT(lep_modmu[0],METlv),mT(lep_modmu[1],METlv))>90.) SRSS[30] = 1;
            if(SRid==2 && (lep_modmu[0]+lep_modmu[1]).M()>40.) SRSS[30] = 2;
          } else if(passJetSSstate(false, false, 0, false, true,1)){
            if(SRid==0 && fabs((lep_modmu[0]+lep_modmu[1]).M()-MZ)>10. && (lep_modmu[0]+lep_modmu[1]).M()>40.) SRSS[32] = 0;
            if(SRid==1 && (lep_modmu[0]+lep_modmu[1]).M()>30. && TMath::Max(mT(lep_modmu[0],METlv),mT(lep_modmu[1],METlv))>90.) SRSS[32] = 1;
            if(SRid==2 && (lep_modmu[0]+lep_modmu[1]).M()>40.&&met_pt()>60.) SRSS[32] = 2;
          }
        }
      }
      if(lep_raw.size()==3&&passJet3lstate()&&nVlep()==3&&abs(lep_q[0]+lep_q[1]+lep_q[2])!=3){
        passAny3l(SR3l[0], dummy, SR3l[1]);
        passAny3l(SR3l[5], dummy, SR3l[1]);

        bool OS01 = (lep_id[0]*lep_id[1]<0); bool SF01 = (abs(lep_id[0])==abs(lep_id[1]));
        bool OS02 = (lep_id[0]*lep_id[2]<0); bool SF02 = (abs(lep_id[0])==abs(lep_id[2]));
        bool OS12 = (lep_id[1]*lep_id[2]<0); bool SF12 = (abs(lep_id[1])==abs(lep_id[2]));
        int numSFOS = 0;
        if(OS01&&SF01) ++numSFOS;
        if(OS02&&SF02) ++numSFOS;
        if(OS12&&SF12) ++numSFOS;
        SRid = numSFOS;
        if(SRid==0){
          if(met_pt()<30.) SRid = -1;
        }
        if(SRid==1){
          if(met_pt()<40.) SRid = -1;
        }
        if(SRid==2){
          if(met_pt()<55.) SRid = -1;
        }
        SR3l[0] = -1;
        //raw ele+muo
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_raw[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_raw[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_raw[0]+lep_raw[1]+lep_raw[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_raw[0]+lep_raw[1]).M()<20.) pass = false;
            if(SF02 && (lep_raw[0]+lep_raw[2]).M()<20.) pass = false;
            if(SF12 && (lep_raw[1]+lep_raw[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_raw[0]+lep_raw[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_raw[0]+lep_raw[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_raw[1]+lep_raw[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_raw[0]+lep_raw[1]+lep_raw[2],METlv)>2.5 && TMath::Max(mT(lep_raw[0],METlv),TMath::Max(mT(lep_raw[1],METlv),mT(lep_raw[2],METlv)))) SR3l[0] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1;
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_raw[0]+lep_raw[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_raw[0]+lep_raw[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_raw[1]+lep_raw[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_raw[0]+lep_raw[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_raw[0]+lep_raw[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_raw[1]+lep_raw[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_raw[0]+lep_raw[1]+lep_raw[2],METlv)>2.5&&mT(lep_raw[third],METlv)>90.&&(lep_raw[0]+lep_raw[1]+lep_raw[2]).Pt()>60.) SR3l[0] = 1;
            if(SRid==2 && pass && dPhi(lep_raw[0]+lep_raw[1]+lep_raw[2],METlv)>2.5&&(lep_raw[0]+lep_raw[1]+lep_raw[2]).Pt()>60.) SR3l[0] = 2;
          }
        }
        //mod ele+muo
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_mod[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_mod[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_mod[0]+lep_mod[1]+lep_mod[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_mod[0]+lep_mod[1]).M()<20.) pass = false;
            if(SF02 && (lep_mod[0]+lep_mod[2]).M()<20.) pass = false;
            if(SF12 && (lep_mod[1]+lep_mod[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_mod[0]+lep_mod[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_mod[0]+lep_mod[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_mod[1]+lep_mod[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_mod[0]+lep_mod[1]+lep_mod[2],METlv)>2.5 && TMath::Max(mT(lep_mod[0],METlv),TMath::Max(mT(lep_mod[1],METlv),mT(lep_mod[2],METlv)))) SR3l[10] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1;
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_mod[0]+lep_mod[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_mod[0]+lep_mod[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_mod[1]+lep_mod[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_mod[0]+lep_mod[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_mod[0]+lep_mod[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_mod[1]+lep_mod[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_mod[0]+lep_mod[1]+lep_mod[2],METlv)>2.5&&mT(lep_mod[third],METlv)>90.&&(lep_mod[0]+lep_mod[1]+lep_mod[2]).Pt()>60.) SR3l[10] = 1;
            if(SRid==2 && pass && dPhi(lep_mod[0]+lep_mod[1]+lep_mod[2],METlv)>2.5&&(lep_mod[0]+lep_mod[1]+lep_mod[2]).Pt()>60.) SR3l[10] = 2;
          }
        }
        //mod ele up ele+muo
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modelup[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modelup[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_modelup[0]+lep_modelup[1]+lep_modelup[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_modelup[0]+lep_modelup[1]).M()<20.) pass = false;
            if(SF02 && (lep_modelup[0]+lep_modelup[2]).M()<20.) pass = false;
            if(SF12 && (lep_modelup[1]+lep_modelup[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_modelup[0]+lep_modelup[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_modelup[0]+lep_modelup[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_modelup[1]+lep_modelup[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_modelup[0]+lep_modelup[1]+lep_modelup[2],METlv)>2.5 && TMath::Max(mT(lep_modelup[0],METlv),TMath::Max(mT(lep_modelup[1],METlv),mT(lep_modelup[2],METlv)))) SR3l[13] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1;
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_modelup[0]+lep_modelup[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_modelup[0]+lep_modelup[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_modelup[1]+lep_modelup[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_modelup[0]+lep_modelup[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_modelup[0]+lep_modelup[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_modelup[1]+lep_modelup[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_modelup[0]+lep_modelup[1]+lep_modelup[2],METlv)>2.5&&mT(lep_modelup[third],METlv)>90.&&(lep_modelup[0]+lep_modelup[1]+lep_modelup[2]).Pt()>60.) SR3l[13] = 1;
            if(SRid==2 && pass && dPhi(lep_modelup[0]+lep_modelup[1]+lep_modelup[2],METlv)>2.5&&(lep_modelup[0]+lep_modelup[1]+lep_modelup[2]).Pt()>60.) SR3l[13] = 2;
          }
        }
        //mod ele dn ele+muo
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modeldn[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modeldn[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_modeldn[0]+lep_modeldn[1]+lep_modeldn[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_modeldn[0]+lep_modeldn[1]).M()<20.) pass = false;
            if(SF02 && (lep_modeldn[0]+lep_modeldn[2]).M()<20.) pass = false;
            if(SF12 && (lep_modeldn[1]+lep_modeldn[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_modeldn[0]+lep_modeldn[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_modeldn[0]+lep_modeldn[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_modeldn[1]+lep_modeldn[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_modeldn[0]+lep_modeldn[1]+lep_modeldn[2],METlv)>2.5 && TMath::Max(mT(lep_modeldn[0],METlv),TMath::Max(mT(lep_modeldn[1],METlv),mT(lep_modeldn[2],METlv)))) SR3l[15] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1;
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_modeldn[0]+lep_modeldn[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_modeldn[0]+lep_modeldn[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_modeldn[1]+lep_modeldn[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_modeldn[0]+lep_modeldn[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_modeldn[0]+lep_modeldn[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_modeldn[1]+lep_modeldn[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_modeldn[0]+lep_modeldn[1]+lep_modeldn[2],METlv)>2.5&&mT(lep_modeldn[third],METlv)>90.&&(lep_modeldn[0]+lep_modeldn[1]+lep_modeldn[2]).Pt()>60.) SR3l[15] = 1;
            if(SRid==2 && pass && dPhi(lep_modeldn[0]+lep_modeldn[1]+lep_modeldn[2],METlv)>2.5&&(lep_modeldn[0]+lep_modeldn[1]+lep_modeldn[2]).Pt()>60.) SR3l[15] = 2;
          }
        }
        //mod muo up ele+muo
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modmuup[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modmuup[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_modmuup[0]+lep_modmuup[1]+lep_modmuup[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_modmuup[0]+lep_modmuup[1]).M()<20.) pass = false;
            if(SF02 && (lep_modmuup[0]+lep_modmuup[2]).M()<20.) pass = false;
            if(SF12 && (lep_modmuup[1]+lep_modmuup[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_modmuup[0]+lep_modmuup[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_modmuup[0]+lep_modmuup[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_modmuup[1]+lep_modmuup[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_modmuup[0]+lep_modmuup[1]+lep_modmuup[2],METlv)>2.5 && TMath::Max(mT(lep_modmuup[0],METlv),TMath::Max(mT(lep_modmuup[1],METlv),mT(lep_modmuup[2],METlv)))) SR3l[23] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1;
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_modmuup[0]+lep_modmuup[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_modmuup[0]+lep_modmuup[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_modmuup[1]+lep_modmuup[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_modmuup[0]+lep_modmuup[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_modmuup[0]+lep_modmuup[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_modmuup[1]+lep_modmuup[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_modmuup[0]+lep_modmuup[1]+lep_modmuup[2],METlv)>2.5&&mT(lep_modmuup[third],METlv)>90.&&(lep_modmuup[0]+lep_modmuup[1]+lep_modmuup[2]).Pt()>60.) SR3l[23] = 1;
            if(SRid==2 && pass && dPhi(lep_modmuup[0]+lep_modmuup[1]+lep_modmuup[2],METlv)>2.5&&(lep_modmuup[0]+lep_modmuup[1]+lep_modmuup[2]).Pt()>60.) SR3l[23] = 2;
          }
        }
        //mod muo dn ele+muo
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modmudn[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modmudn[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_modmudn[0]+lep_modmudn[1]+lep_modmudn[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_modmudn[0]+lep_modmudn[1]).M()<20.) pass = false;
            if(SF02 && (lep_modmudn[0]+lep_modmudn[2]).M()<20.) pass = false;
            if(SF12 && (lep_modmudn[1]+lep_modmudn[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_modmudn[0]+lep_modmudn[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_modmudn[0]+lep_modmudn[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_modmudn[1]+lep_modmudn[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_modmudn[0]+lep_modmudn[1]+lep_modmudn[2],METlv)>2.5 && TMath::Max(mT(lep_modmudn[0],METlv),TMath::Max(mT(lep_modmudn[1],METlv),mT(lep_modmudn[2],METlv)))) SR3l[25] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1;
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_modmudn[0]+lep_modmudn[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_modmudn[0]+lep_modmudn[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_modmudn[1]+lep_modmudn[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_modmudn[0]+lep_modmudn[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_modmudn[0]+lep_modmudn[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_modmudn[1]+lep_modmudn[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_modmudn[0]+lep_modmudn[1]+lep_modmudn[2],METlv)>2.5&&mT(lep_modmudn[third],METlv)>90.&&(lep_modmudn[0]+lep_modmudn[1]+lep_modmudn[2]).Pt()>60.) SR3l[25] = 1;
            if(SRid==2 && pass && dPhi(lep_modmudn[0]+lep_modmudn[1]+lep_modmudn[2],METlv)>2.5&&(lep_modmudn[0]+lep_modmudn[1]+lep_modmudn[2]).Pt()>60.) SR3l[25] = 2;
          }
        }
        //mod ele
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_model[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_model[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_model[0]+lep_model[1]+lep_model[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_model[0]+lep_model[1]).M()<20.) pass = false;
            if(SF02 && (lep_model[0]+lep_model[2]).M()<20.) pass = false;
            if(SF12 && (lep_model[1]+lep_model[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_model[0]+lep_model[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_model[0]+lep_model[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_model[1]+lep_model[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_model[0]+lep_model[1]+lep_model[2],METlv)>2.5 && TMath::Max(mT(lep_model[0],METlv),TMath::Max(mT(lep_model[1],METlv),mT(lep_model[2],METlv)))) SR3l[20] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1;
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_model[0]+lep_model[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_model[0]+lep_model[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_model[1]+lep_model[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_model[0]+lep_model[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_model[0]+lep_model[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_model[1]+lep_model[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_model[0]+lep_model[1]+lep_model[2],METlv)>2.5&&mT(lep_model[third],METlv)>90.&&(lep_model[0]+lep_model[1]+lep_model[2]).Pt()>60.) SR3l[20] = 1;
            if(SRid==2 && pass && dPhi(lep_model[0]+lep_model[1]+lep_model[2],METlv)>2.5&&(lep_model[0]+lep_model[1]+lep_model[2]).Pt()>60.) SR3l[20] = 2;
          }
        }
        //mod muo
        nleps  = 0;
        nleps2 = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modmu[i].Pt()>25.) ++nleps2;
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_modmu[i].Pt()>20.) ++nleps;
        }
        if(nleps2>=1&&nleps==3 && fabs((lep_modmu[0]+lep_modmu[1]+lep_modmu[2]).M()-MZ)>10.){
          if(SRid==0 ){
            bool pass = true;
            if(SF01 && (lep_modmu[0]+lep_modmu[1]).M()<20.) pass = false;
            if(SF02 && (lep_modmu[0]+lep_modmu[2]).M()<20.) pass = false;
            if(SF12 && (lep_modmu[1]+lep_modmu[2]).M()<20.) pass = false;
            if(SF01 && abs(lep_id[0])==11 && fabs((lep_modmu[0]+lep_modmu[1]).M()-MZ)<15.) pass = false;
            if(SF02 && abs(lep_id[0])==11 && fabs((lep_modmu[0]+lep_modmu[2]).M()-MZ)<15.) pass = false;
            if(SF12 && abs(lep_id[1])==11 && fabs((lep_modmu[1]+lep_modmu[2]).M()-MZ)<15.) pass = false;
            if(pass && dPhi(lep_modmu[0]+lep_modmu[1]+lep_modmu[2],METlv)>2.5 && TMath::Max(mT(lep_modmu[0],METlv),TMath::Max(mT(lep_modmu[1],METlv),mT(lep_modmu[2],METlv)))) SR3l[30] = 0;
          }
          if(SRid>=1){
            bool pass = true;
            int third = -1; 
            if(OS01&&SF01) third = 2;
            if(OS02&&SF02) third = 1;
            if(OS12&&SF12) third = 0;
            if(OS01&&SF01 && (lep_modmu[0]+lep_modmu[1]).M()<20.) pass = false;
            if(OS02&&SF02 && (lep_modmu[0]+lep_modmu[2]).M()<20.) pass = false;
            if(OS12&&SF12 && (lep_modmu[1]+lep_modmu[2]).M()<20.) pass = false;
            if(OS01&&SF01 && fabs((lep_modmu[0]+lep_modmu[1]).M()-MZ)<20.) pass = false;
            if(OS02&&SF02 && fabs((lep_modmu[0]+lep_modmu[2]).M()-MZ)<20.) pass = false;
            if(OS12&&SF12 && fabs((lep_modmu[1]+lep_modmu[2]).M()-MZ)<20.) pass = false;
            if(SRid==1 && pass && dPhi(lep_modmu[0]+lep_modmu[1]+lep_modmu[2],METlv)>2.5&&mT(lep_modmu[third],METlv)>90.&&(lep_modmu[0]+lep_modmu[1]+lep_modmu[2]).Pt()>60.) SR3l[30] = 1;
            if(SRid==2 && pass && dPhi(lep_modmu[0]+lep_modmu[1]+lep_modmu[2],METlv)>2.5&&(lep_modmu[0]+lep_modmu[1]+lep_modmu[2]).Pt()>60.) SR3l[30] = 2;
          }
        }
      }
      //raw - SR 0 /2 (2 only for SS: Mjjside)
      //all   SR 10 / 12
      //ele   SR 20 / 22
      //muo   SR 30 / 32
      //fill SR histos here
      if(SRSS[ 0]>=0) fillhisto(  histos, "AllSR_raw",                     sample, sn, SRSS[ 0],   weight);
      if(SR3l[ 0]>=0) fillhisto(  histos, "AllSR_raw",                     sample, sn, SR3l[ 0]+6, weight);
      if(SRSS[ 2]>=0) fillhisto(  histos, "AllSR_raw",                     sample, sn, SRSS[ 2]+3, weight);
      if(SRSS[ 5]>=0) fillhisto(  histos, "AllSR_xcheck",                  sample, sn, SRSS[ 5],   weight);
      if(SR3l[ 5]>=0) fillhisto(  histos, "AllSR_xcheck",                  sample, sn, SR3l[ 5]+6, weight); 
      if(SRSS[ 6]>=0) fillhisto(  histos, "AllSR_xcheck",                  sample, sn, SRSS[ 6]+3, weight);
      if(SRSS[10]>=0) fillhisto(  histos, "AllSR_allsmearedscaled",        sample, sn, SRSS[10],   weight);
      if(SR3l[10]>=0) fillhisto(  histos, "AllSR_allsmearedscaled",        sample, sn, SR3l[10]+6, weight);
      if(SRSS[12]>=0) fillhisto(  histos, "AllSR_allsmearedscaled",        sample, sn, SRSS[12]+3, weight);
      if(SRSS[20]>=0) fillhisto(  histos, "AllSR_elesmearedscaled",        sample, sn, SRSS[20],   weight);
      if(SR3l[20]>=0) fillhisto(  histos, "AllSR_elesmearedscaled",        sample, sn, SR3l[20]+6, weight);
      if(SRSS[22]>=0) fillhisto(  histos, "AllSR_elesmearedscaled",        sample, sn, SRSS[22]+3, weight);
      if(SRSS[30]>=0) fillhisto(  histos, "AllSR_muosmearedscaled",        sample, sn, SRSS[30],   weight);
      if(SR3l[30]>=0) fillhisto(  histos, "AllSR_muosmearedscaled",        sample, sn, SR3l[30]+6, weight);
      if(SRSS[32]>=0) fillhisto(  histos, "AllSR_muosmearedscaled",        sample, sn, SRSS[32]+3, weight);
      if(SRSS[13]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_eleup",  sample, sn, SRSS[13],   weight);
      if(SR3l[13]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_eleup",  sample, sn, SR3l[13]+6, weight);
      if(SRSS[14]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_eleup",  sample, sn, SRSS[14]+3, weight);
      if(SRSS[15]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_eledn",  sample, sn, SRSS[15],   weight);
      if(SR3l[15]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_eledn",  sample, sn, SR3l[15]+6, weight);
      if(SRSS[16]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_eledn",  sample, sn, SRSS[16]+3, weight);
      if(SRSS[23]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_muoup",  sample, sn, SRSS[23],   weight);
      if(SR3l[23]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_muoup",  sample, sn, SR3l[23]+6, weight);
      if(SRSS[24]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_muoup",  sample, sn, SRSS[24]+3, weight);
      if(SRSS[25]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_muodn",  sample, sn, SRSS[25],   weight);
      if(SR3l[25]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_muodn",  sample, sn, SR3l[25]+6, weight);
      if(SRSS[26]>=0) fillhisto(  histos, "AllSR_allsmearedscaled_muodn",  sample, sn, SRSS[26]+3, weight);

      //fill Mll here for both OS and SS and some variations.
      if(lep_raw.size()==2){
        int elemuo = 0;
        if((lep_id[0]*lep_id[1])==(-121)) elemuo =  1;//ee OS
        if((lep_id[0]*lep_id[1])==( 121)) elemuo = -1;//ee SS
        if((lep_id[0]*lep_id[1])==(-169)) elemuo =  2;//mm OS
        if((lep_id[0]*lep_id[1])==( 169)) elemuo = -2;//mm SS
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_raw[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(elemuo== 1) fillhisto(  histos, "Mee_raw",          sample, sn, (lep_raw[0]+lep_raw[1]).M(), weight);
          if(elemuo==-1) fillhisto(  histos, "MeeSS_raw",        sample, sn, (lep_raw[0]+lep_raw[1]).M(), weight);
          if(elemuo== 2) fillhisto(  histos, "Mmumu_raw",        sample, sn, (lep_raw[0]+lep_raw[1]).M(), weight);
          if(elemuo==-2) fillhisto(  histos, "MmumuSS_raw",      sample, sn, (lep_raw[0]+lep_raw[1]).M(), weight);
        }
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_mod[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(elemuo== 1) fillhisto(  histos, "Mee_smearedandscaled",          sample, sn, (lep_mod[0]+lep_mod[1]).M(), weight);
          if(elemuo==-1) fillhisto(  histos, "MeeSS_smearedandscaled",        sample, sn, (lep_mod[0]+lep_mod[1]).M(), weight);
          if(elemuo== 2) fillhisto(  histos, "Mmumu_smearedandscaled",        sample, sn, (lep_mod[0]+lep_mod[1]).M(), weight);
          if(elemuo==-2) fillhisto(  histos, "MmumuSS_smearedandscaled",      sample, sn, (lep_mod[0]+lep_mod[1]).M(), weight);
        }
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_scupel[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(elemuo== 1) fillhisto(  histos, "Mee_smearedandscaled_scaledup",   sample, sn, (lep_scupel[0]+lep_scupel[1]).M(), weight);
          if(elemuo==-1) fillhisto(  histos, "MeeSS_smearedandscaled_scaledup", sample, sn, (lep_scupel[0]+lep_scupel[1]).M(), weight);
        }
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_scdnel[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(elemuo== 1) fillhisto(  histos, "Mee_smearedandscaled_scaleddn",   sample, sn, (lep_scdnel[0]+lep_scdnel[1]).M(), weight);
          if(elemuo==-1) fillhisto(  histos, "MeeSS_smearedandscaled_scaleddn", sample, sn, (lep_scdnel[0]+lep_scdnel[1]).M(), weight);
        }
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_smupel[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(elemuo== 1) fillhisto(  histos, "Mee_smearedandscaled_smearup",    sample, sn, (lep_smupel[0]+lep_smupel[1]).M(), weight);
          if(elemuo==-1) fillhisto(  histos, "MeeSS_smearedandscaled_smearup",  sample, sn, (lep_smupel[0]+lep_smupel[1]).M(), weight);
        }
        nleps = 0;
        for(unsigned int i = 0; i<lep_pdgId().size();++i){
          if(lep_pass_VVV_cutbased_3l_tight   ()[i]&&lep_smdnel[i].Pt()>25.) ++nleps;
        }
        if(nleps==2){
          if(elemuo== 1) fillhisto(  histos, "Mee_smearedandscaled_smeardn",    sample, sn, (lep_smdnel[0]+lep_smdnel[1]).M(), weight);
          if(elemuo==-1) fillhisto(  histos, "MeeSS_smearedandscaled_smeardn",  sample, sn, (lep_smdnel[0]+lep_smdnel[1]).M(), weight);
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

  SaveHistosToFile("rootfiles/LeptonSmearingResults.root",histos,true,true,(chainnumber==0));
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
