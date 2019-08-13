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
#include "Functions_test.h"//update this line
#include "CMS3.cc"//update this line
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/dorky/dorky.cc"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/goodrun.cc"

using namespace std;
using namespace tas;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test", int dummy=-1) {
  //int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test"/*, string outdir="/nfs-7/userdata/haweber/WWWskims/WWW_v0.1.16/"*/) {
//vector<string> ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test"/*, string outdir="/nfs-7/userdata/haweber/WWWskims/WWW_v0.1.16/"*/) {

  //if(skimFilePrefix=="Background") return 0;
  bool mergesample = false;//probably false is good
  string outdir = "/nfs-7/userdata/haweber/WWWskims/WWW_v1.0.11/";
  string newname = outdir+skimFilePrefix+".root";
  vector<string> filestomerge;

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  const char* json_file = "data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
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
    tree->SetBranchStatus("lumi",1);
    tree->SetBranchStatus("evt",1);
    tree->SetBranchStatus("isData",1);
    tree->SetBranchStatus("evt_scale1fb",1);
    tree->SetBranchStatus("evt_passgoodrunlist",1);
    tree->SetBranchStatus("HLT_DoubleMu",1);
    tree->SetBranchStatus("HLT_DoubleEl",1);
    tree->SetBranchStatus("HLT_DoubleEl_DZ",1);
    tree->SetBranchStatus("HLT_MuEG",1);
    tree->SetBranchStatus("lep_p4*",1);
    tree->SetBranchStatus("lep_coneCorrPt",1);
    tree->SetBranchStatus("lep_ip3d",1);
    tree->SetBranchStatus("lep_isTriggerSafe_v1",1);
    tree->SetBranchStatus("lep_motherIdSS",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_3l_fo",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_3l_tight",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_fo",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_tight",1);
    tree->SetBranchStatus("lep_pass_VVV_cutbased_veto",1);
    tree->SetBranchStatus("lep_pdgId",1);
    tree->SetBranchStatus("lep_ptRatio",1);
    tree->SetBranchStatus("lep_tightCharge",1);
    tree->SetBranchStatus("lep_etaSC",1);
    tree->SetBranchStatus("lep_MVA",1);
    tree->SetBranchStatus("lep_isFromW",1);
    tree->SetBranchStatus("lep_isFromZ",1);
    tree->SetBranchStatus("lep_genPart_index",1);
    tree->SetBranchStatus("jets_p4*",1);
    tree->SetBranchStatus("jets_up_p4*",1);
    tree->SetBranchStatus("jets_dn_p4*",1);
    tree->SetBranchStatus("jets_csv",1);
    tree->SetBranchStatus("jets_up_csv",1);
    tree->SetBranchStatus("jets_dn_csv",1);
    tree->SetBranchStatus("met_pt",1);
    tree->SetBranchStatus("met_phi",1);
    tree->SetBranchStatus("met_up_pt",1);
    tree->SetBranchStatus("met_up_phi",1);
    tree->SetBranchStatus("met_dn_pt",1);
    tree->SetBranchStatus("met_dn_phi",1);
    tree->SetBranchStatus("firstgoodvertex",1);
    tree->SetBranchStatus("nTrueInt",1);
    tree->SetBranchStatus("nVert",1);
    tree->SetBranchStatus("nisoTrack_mt2_cleaned_VVV_cutbased_veto",1);
    tree->SetBranchStatus("weight_btagsf",1);
    tree->SetBranchStatus("weight_btagsf_heavy_DN",1);
    tree->SetBranchStatus("weight_btagsf_heavy_UP",1);
    tree->SetBranchStatus("weight_btagsf_light_DN",1);
    tree->SetBranchStatus("weight_btagsf_light_UP",1);
    tree->SetBranchStatus("gen_ht",1);
    tree->SetBranchStatus("genPart_p4*",1);
    tree->SetBranchStatus("genPart_motherId",1);
    tree->SetBranchStatus("genPart_pdgId",1);
    tree->SetBranchStatus("genPart_status",1);
    tree->SetBranchStatus("ngenLep",1);
    tree->SetBranchStatus("ngenLepFromTau",1);
    tree->SetBranchStatus("Flag_AllEventFilters",1);
    tree->SetBranchStatus("nVlep",1);
    tree->SetBranchStatus("nTlep",1);
    tree->SetBranchStatus("nTlepSS",1);
    tree->SetBranchStatus("nLlep",1);
    tree->SetBranchStatus("nSFOS",1);
    tree->SetBranchStatus("nSFOSinZ",1);
    tree->SetBranchStatus("nj",1);
    tree->SetBranchStatus("nj_up",1);
    tree->SetBranchStatus("nj_dn",1);
    tree->SetBranchStatus("nj30",1);
    tree->SetBranchStatus("nj30_up",1);
    tree->SetBranchStatus("nj30_dn",1);
    tree->SetBranchStatus("nb",1);
    tree->SetBranchStatus("nb_up",1);
    tree->SetBranchStatus("nb_dn",1);
    tree->SetBranchStatus("Mjj",1);
    tree->SetBranchStatus("Mjj_up",1);
    tree->SetBranchStatus("Mjj_dn",1);
    tree->SetBranchStatus("MjjL",1);
    tree->SetBranchStatus("MjjL_up",1);
    tree->SetBranchStatus("MjjL_dn",1);
    tree->SetBranchStatus("DetajjL",1);
    tree->SetBranchStatus("DetajjL_up",1);
    tree->SetBranchStatus("DetajjL_dn",1);
    tree->SetBranchStatus("MllSS",1);
    tree->SetBranchStatus("MeeSS",1);
    tree->SetBranchStatus("Mll3L",1);
    tree->SetBranchStatus("Mee3L",1);
    tree->SetBranchStatus("Mll3L1",1);
    tree->SetBranchStatus("M3l",1);
    tree->SetBranchStatus("Pt3l",1);
    tree->SetBranchStatus("DPhi3lMET",1);
    tree->SetBranchStatus("DPhi3lMET_up",1);
    tree->SetBranchStatus("DPhi3lMET_dn",1);
    tree->SetBranchStatus("MTmax",1);
    tree->SetBranchStatus("MTmax_up",1);
    tree->SetBranchStatus("MTmax_dn",1);
    tree->SetBranchStatus("MT3rd",1);
    tree->SetBranchStatus("MT3rd_up",1);
    tree->SetBranchStatus("MT3rd_dn",1);
    tree->SetBranchStatus("passSSee",1);
    tree->SetBranchStatus("passSSem",1);
    tree->SetBranchStatus("passSSmm",1);
    tree->SetBranchStatus("lep_idx0_SS",1);
    tree->SetBranchStatus("lep_idx1_SS",1);
    tree->SetBranchStatus("bkgtype",1);
    tree->SetBranchStatus("vetophoton",1);
    tree->SetBranchStatus("purewgt",1);
    tree->SetBranchStatus("purewgt_up",1);
    tree->SetBranchStatus("purewgt_dn",1);
    tree->SetBranchStatus("ffwgt",1);
    tree->SetBranchStatus("ffwgt_up",1);
    tree->SetBranchStatus("ffwgt_dn",1);
    tree->SetBranchStatus("ffwgtqcd",1);
    tree->SetBranchStatus("ffwgtqcd_up",1);
    tree->SetBranchStatus("ffwgtqcd_dn",1);
    tree->SetBranchStatus("lepsf",1);
    tree->SetBranchStatus("lepsf_up",1);
    tree->SetBranchStatus("lepsf_dn",1);
    tree->SetBranchStatus("trigeff",1);
    tree->SetBranchStatus("trigeff_up",1);
    tree->SetBranchStatus("trigeff_dn",1);
    tree->SetBranchStatus("M01",1);
    tree->SetBranchStatus("M02",1);
    tree->SetBranchStatus("M12",1);
    tree->SetBranchStatus("isSFOS01",1);
    tree->SetBranchStatus("isSFOS02",1);
    tree->SetBranchStatus("isSFOS12",1);
    
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

      if(firstgoodvertex()!=0)                     continue;
      if(nVert()<0)                                continue;
      if(nVlep()>3)                                continue;
      if(NtightSS()<1&&Ntight3l()<2)               continue;
      if(NlooseSS()<2&&Nloose3l()<3)               continue;
      if(NlooseSS()==2&&Nloose3l()<3){
	if(nj30()<1)                               continue;//for SS, require at least one jet
	if(!passSSee()&&!passSSem()&&!passSSmm())  continue;//require leptons to be SS
      }
      if(!passTriggers()) continue;//pass trigger for data, and offline lepton kinematic cuts for data/simulation
      //weight = 1;
      if(string(currentFile->GetTitle()).find("wjets_incl_mgmlm_" )!=string::npos && gen_ht()>100) continue;
      if(string(currentFile->GetTitle()).find("dy_m50_mgmlm_ext1_")!=string::npos && gen_ht()>100) continue;

      if(isData()){
	if(!passFilters()) continue;
	duplicate_removal::DorkyEventIdentifier id(tas::run(), tas::evt(), tas::lumi());
	if( is_duplicate(id)        )           continue;
	if( !goodrun(tas::run(), tas::lumi()) ) continue;
      }

	//need to activate all branches
      tas::run(); tas::lumi(); tas::evt(); tas::isData(); tas::evt_scale1fb(); tas::evt_passgoodrunlist(); tas::HLT_DoubleMu(); tas::HLT_DoubleEl(); tas::HLT_DoubleEl_DZ(); tas::HLT_MuEG(); tas::lep_p4(); tas::lep_coneCorrPt(); tas::lep_ip3d(); tas::lep_isTriggerSafe_v1(); tas::lep_motherIdSS(); tas::lep_pass_VVV_cutbased_3l_fo(); tas::lep_pass_VVV_cutbased_3l_tight(); tas::lep_pass_VVV_cutbased_fo(); tas::lep_pass_VVV_cutbased_tight(); tas::lep_pass_VVV_cutbased_veto(); tas::lep_pdgId(); tas::lep_ptRatio(); tas::lep_tightCharge(); tas::lep_etaSC(); tas::lep_MVA(); tas::lep_isFromW(); tas::lep_isFromZ(); tas::lep_genPart_index(); tas::jets_p4(); tas::jets_up_p4(); tas::jets_dn_p4(); tas::jets_csv(); tas::jets_up_csv(); tas::jets_dn_csv(); tas::met_pt(); tas::met_phi(); tas::met_up_pt(); tas::met_up_phi(); tas::met_dn_pt(); tas::met_dn_phi(); tas::firstgoodvertex(); tas::nTrueInt(); tas::nVert(); tas::nisoTrack_mt2_cleaned_VVV_cutbased_veto(); tas::weight_btagsf(); tas::weight_btagsf_heavy_DN(); tas::weight_btagsf_heavy_UP(); tas::weight_btagsf_light_DN(); tas::weight_btagsf_light_UP(); tas::gen_ht(); tas::genPart_p4(); tas::genPart_motherId(); tas::genPart_pdgId(); tas::genPart_status(); tas::ngenLep(); tas::ngenLepFromTau(); tas::Flag_AllEventFilters(); tas::nVlep(); tas::nTlep(); tas::nTlepSS(); tas::nLlep(); tas::nSFOS(); tas::nSFOSinZ(); tas::nj(); tas::nj_up(); tas::nj_dn(); tas::nj30(); tas::nj30_up(); tas::nj30_dn(); tas::nb(); tas::nb_up(); tas::nb_dn(); tas::Mjj(); tas::Mjj_up(); tas::Mjj_dn(); tas::MjjL(); tas::MjjL_up(); tas::MjjL_dn(); tas::DetajjL(); tas::DetajjL_up(); tas::DetajjL_dn(); tas::MllSS(); tas::MeeSS(); tas::Mll3L(); tas::Mee3L(); tas::Mll3L1(); tas::M3l(); tas::Pt3l(); tas::DPhi3lMET(); tas::DPhi3lMET_up(); tas::DPhi3lMET_dn(); tas::MTmax(); tas::MTmax_up(); tas::MTmax_dn(); tas::MT3rd(); tas::MT3rd_up(); tas::MT3rd_dn(); tas::passSSee(); tas::passSSem(); tas::passSSmm(); tas::lep_idx0_SS(); tas::lep_idx1_SS(); tas::bkgtype(); tas::vetophoton(); tas::purewgt(); tas::purewgt_up(); tas::purewgt_dn(); tas::ffwgt(); tas::ffwgt_up(); tas::ffwgt_dn(); tas::ffwgtqcd(); tas::ffwgtqcd_up(); tas::ffwgtqcd_dn(); tas::lepsf(); tas::lepsf_up(); tas::lepsf_dn(); tas::trigeff(); tas::trigeff_up(); tas::trigeff_dn(); tas::M01(); tas::M02(); tas::M12(); tas::isSFOS01(); tas::isSFOS02(); tas::isSFOS12();

      newtree->Fill();
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
  //return filestomerge;
  return 1;
}
