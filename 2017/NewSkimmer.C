// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TSystem.h"
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

vector<string> ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test"/*, string outdir="/nfs-7/userdata/haweber/WWWskims/WWW_v0.1.16/"*/) {

  bool mergesample = false;//probably false is good
  string outdir = "/nfs-7/userdata/haweber/WWWskims/WWW_v0.1.16/";
  string newname = outdir+skimFilePrefix+".root";
  vector<string> filestomerge;

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

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
  bool firsttree = true;
  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    string fname = currentFile->GetTitle();

    TTree *tree = (TTree*)file->Get("t");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    cms3.Init(tree);
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("nVert",1);
    tree->SetBranchStatus("nisoTrack_mt2_cleaned_VVV_cutbased_veto",1);
    tree->SetBranchStatus("ngenLepFromTau",1);
    tree->SetBranchStatus("ngenLep",1);
    tree->SetBranchStatus("met_pt",1);
    tree->SetBranchStatus("met_phi",1);
    tree->SetBranchStatus("lumi",1);
    tree->SetBranchStatus("lep_relIso03EAv2",1);
    tree->SetBranchStatus("lep_relIso03EA",1);
    tree->SetBranchStatus("lep_pdgId",1);
    tree->SetBranchStatus("nlep_VVV_cutbased_veto",1);
    tree->SetBranchStatus("lep_genPart_index",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_veto_noiso",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_veto",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_tight_noiso",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_tight",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_fo_noiso",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_fo",1);
    tree->SetBranchStatus("lep_p4*",1);//
    tree->SetBranchStatus("lep_motherIdSS",1);
    tree->SetBranchStatus("lep_mc_Id",1);
    tree->SetBranchStatus("lep_lostHits",1);
    tree->SetBranchStatus("lep_isTriggerSafe_v2",1);
    tree->SetBranchStatus("lep_isTriggerSafe_v1",1);
    tree->SetBranchStatus("lep_isFromZ",1);
    tree->SetBranchStatus("lep_isFromW",1);
    tree->SetBranchStatus("lep_isFromLF",1);
    tree->SetBranchStatus("lep_isFromL",1);
    tree->SetBranchStatus("lep_isFromC",1);
    tree->SetBranchStatus("lep_isFromB",1);
    tree->SetBranchStatus("lep_ip3derr",1);
    tree->SetBranchStatus("lep_ip3d",1);
    tree->SetBranchStatus("lep_etaSC",1);
    tree->SetBranchStatus("lep_dz",1);
    tree->SetBranchStatus("lep_dxy",1);
    tree->SetBranchStatus("lep_charge",1);
    tree->SetBranchStatus("lep_3ch_agree",1);
    tree->SetBranchStatus("jets_p4*",1);
    tree->SetBranchStatus("jets_csv",1);
    tree->SetBranchStatus("isData",1);
    tree->SetBranchStatus("genPart_charge",1);
    tree->SetBranchStatus("genPart_status",1);
    tree->SetBranchStatus("genPart_pdgId",1);
    tree->SetBranchStatus("genPart_p4*",1);
    tree->SetBranchStatus("genPart_motherId",1);
    tree->SetBranchStatus("firstgoodvertex",1);
    tree->SetBranchStatus("evt_scale1fb",1);
    tree->SetBranchStatus("evt_passgoodrunlist",1);
    tree->SetBranchStatus("evt_dataset",1);
    tree->SetBranchStatus("evt",1);
    tree->SetBranchStatus("HLT_MuEG_noiso",1);
    tree->SetBranchStatus("HLT_MuEG",1);
    tree->SetBranchStatus("HLT_DoubleMu_noiso",1);
    tree->SetBranchStatus("HLT_DoubleMu",1);
    tree->SetBranchStatus("HLT_DoubleEl_noiso",1);
    tree->SetBranchStatus("HLT_DoubleEl_DZ_2",1);
    tree->SetBranchStatus("HLT_DoubleEl_DZ",1);
    tree->SetBranchStatus("HLT_DoubleEl",1);
    tree->SetBranchStatus("Flag_goodVertices",1);
    tree->SetBranchStatus("Flag_globalTightHalo2016",1);
    tree->SetBranchStatus("Flag_eeBadScFilter",1);
    tree->SetBranchStatus("Flag_badMuonFilter",1);
    tree->SetBranchStatus("Flag_badChargedCandidateFilter",1);
    tree->SetBranchStatus("Flag_HBHENoiseFilter",1);
    tree->SetBranchStatus("Flag_HBHEIsoNoiseFilter",1);
    tree->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
    tree->SetBranchStatus("nlep",1);
    tree->SetBranchStatus("lep_tightCharge",1);
    tree->SetBranchStatus("jets_up_p4*",1);
    tree->SetBranchStatus("jets_up_csv",1);
    tree->SetBranchStatus("jets_dn_p4*",1);
    tree->SetBranchStatus("jets_dn_csv",1);
    tree->SetBranchStatus("met_T1CHS_miniAOD_CORE_up_pt",1);
    tree->SetBranchStatus("met_T1CHS_miniAOD_CORE_up_phi",1);
    tree->SetBranchStatus("met_T1CHS_miniAOD_CORE_dn_pt",1);
    tree->SetBranchStatus("met_T1CHS_miniAOD_CORE_dn_phi",1);
    tree->SetBranchStatus("gen_ht",1);
    tree->SetBranchStatus("weight_btagsf",1);
    tree->SetBranchStatus("weight_btagsf_heavy_UP",1);
    tree->SetBranchStatus("weight_btagsf_heavy_DN",1);
    tree->SetBranchStatus("weight_btagsf_light_UP",1);
    tree->SetBranchStatus("weight_btagsf_light_DN",1);
    vector<string> tempvec = split(fname);
    if(tempvec.size()==0) continue;
    
    string nametemp = outdir+tempvec[tempvec.size()-1];
    TFile *ftemp = new TFile(nametemp.c_str(),"recreate");
    TTree *newtree = (TTree*)tree->CloneTree(0);

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

      int nj(0), nb(0), nj30(0);
      getalljetnumbers(nj,nj30,nb);
      vector<int> vSS,   v3l,   iSS,   i3l; //lepton indices for both the SS and 3l signal regions
      vector<int> vaSS,  va3l,  iaSS,  ia3l;//loose, but not tight leptons.
      getleptonindices(iSS, i3l, iaSS, ia3l, vSS, v3l, vaSS, va3l);
      
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
      if(((nSS+naSS)==2&&nSS>=1&&nvetoaSS==0&&nj30>=2)||(((n3l+na3l)==3||n3l>=3)&&n3l>=2&&lep_p4()[i3l[0] ].Pt()>25.)) {
	if((n3l+na3l)<3){
	  if(nSS==2&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[ iSS[1] ]))<0) continue;
	  if(nSS==1&&((lep_pdgId()[iSS[0] ])*(lep_pdgId()[iaSS[0] ]))<0) continue;
	}
	//need to activate all branches
	tas::run();  tas::nVert();  tas::nisoTrack_mt2_cleaned_VVV_cutbased_veto();  tas::ngenLepFromTau();  tas::ngenLep();  tas::met_pt();  tas::met_phi();  tas::lumi();  tas::lep_relIso03EAv2();  tas::lep_relIso03EA();  tas::lep_pdgId();  tas::nlep_VVV_cutbased_veto();  tas::lep_genPart_index();  tas::lep_pass_VVV_cutbased_veto_noiso();  tas::lep_pass_VVV_cutbased_veto();  tas::lep_pass_VVV_cutbased_tight_noiso();  tas::lep_pass_VVV_cutbased_tight();  tas::lep_pass_VVV_cutbased_fo_noiso();  tas::lep_pass_VVV_cutbased_fo();  tas::lep_p4();  tas::lep_motherIdSS();  tas::lep_mc_Id();  tas::lep_lostHits();  tas::lep_isTriggerSafe_v2();  tas::lep_isTriggerSafe_v1();  tas::lep_isFromZ();  tas::lep_isFromW();  tas::lep_isFromLF();  tas::lep_isFromL();  tas::lep_isFromC();  tas::lep_isFromB();  tas::lep_ip3derr();  tas::lep_ip3d();  tas::lep_etaSC();  tas::lep_dz();  tas::lep_dxy();  tas::lep_charge();  tas::lep_3ch_agree();  tas::jets_p4();  tas::jets_csv();  tas::isData();  tas::genPart_charge();  tas::genPart_status();  tas::genPart_pdgId();  tas::genPart_p4();  tas::genPart_motherId();  tas::firstgoodvertex();  tas::evt_scale1fb();  tas::evt_passgoodrunlist();  tas::evt_dataset();  tas::evt();  tas::HLT_MuEG_noiso();  tas::HLT_MuEG();  tas::HLT_DoubleMu_noiso();  tas::HLT_DoubleMu();  tas::HLT_DoubleEl_noiso();  tas::HLT_DoubleEl_DZ_2();  tas::HLT_DoubleEl_DZ();  tas::HLT_DoubleEl();  tas::Flag_goodVertices();  tas::Flag_globalTightHalo2016();  tas::Flag_eeBadScFilter();  tas::Flag_badMuonFilter();  tas::Flag_badChargedCandidateFilter();  tas::Flag_HBHENoiseFilter();  tas::Flag_HBHEIsoNoiseFilter();  tas::Flag_EcalDeadCellTriggerPrimitiveFilter();  tas::nlep();  tas::lep_tightCharge();  tas::jets_up_p4();  tas::jets_up_csv();  tas::jets_dn_p4();  tas::jets_dn_csv();  tas::met_T1CHS_miniAOD_CORE_up_pt();  tas::met_T1CHS_miniAOD_CORE_up_phi();  tas::met_T1CHS_miniAOD_CORE_dn_pt();  tas::met_T1CHS_miniAOD_CORE_dn_phi(); tas::gen_ht(); tas::weight_btagsf(); tas::weight_btagsf_heavy_UP(); tas::weight_btagsf_heavy_DN(); tas::weight_btagsf_light_UP(); tas::weight_btagsf_light_DN();

	newtree->Fill();
      }
    }//event loop

    filestomerge.push_back(nametemp);
    newtree->AutoSave();
    ftemp->Close();
    delete ftemp;
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }//file loop
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  if(mergesample){
    TChain *newchain = new TChain("t");
    newchain->SetMaxTreeSize(5000000000LL);
    //for(unsigned int i = 0; i<filestomerge.size();++i) cout << filestomerge[i].c_str() << endl;
    for(unsigned int i = 0; i<filestomerge.size();++i) newchain->Add(filestomerge[i].c_str());
    int nFiles_out = newchain->Merge(newname.c_str(), "fast");
    cout << "merged skim for " << newname << endl;
    for(unsigned int i = 0; i<filestomerge.size();++i){
      string command = "rm " + filestomerge[i];
      gSystem->Exec(command.c_str());
    }
  }

  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return filestomerge;
}
