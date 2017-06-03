{

  //gROOT->ProcessLine(".L InvestigateTribosons.C+");
  //gROOT->ProcessLine(".L DeriveTreeForTrainBDTv2.C+");
  //gROOT->ProcessLine(".L InvestigateTribosonsAddition.C+");
  //gROOT->ProcessLine(".L StupidSync.C+");
  //gROOT->ProcessLine(".L StupidSyncCheck.C+");
  //gROOT->ProcessLine(".L StudyJetPtEtaN.C+");
  //gROOT->ProcessLine(".L StudyMjj.C+");
  //gROOT->ProcessLine(".L StupidMll.C+");
  gROOT->ProcessLine(".L StudyVariousObjects.C+");
  const unsigned int chainsize = 11;
  TChain *ch[chainsize];
  string dataset[chainsize];


  //string babylocation = "/hadoop/cms/store/user/haweber/AutoTwopler_babies/WWW/merged/ZMET/WWW_vX/skim/";
  string babylocation = "/nfs-7/userdata/bhashemi/WWW_babies/WWW_v0.1.4/skim/";
  string myhelper;
  
  dataset[0] = "WWW";
  ch[0] = new TChain("t");
  //myhelper = babylocation + "www_incl_amcnlo_*.root"; ch[0]->Add(myhelper.c_str());//240000
  myhelper = babylocation + "www_2l_mia_*.root"; ch[0]->Add(myhelper.c_str());//91900
  myhelper = babylocation + "www_2l_ext1_mia_*.root"; ch[0]->Add(myhelper.c_str());//164800

  dataset[1] = "TT2l";
  ch[1] = new TChain("t");
  myhelper = babylocation + "ttbar_dilep_mgmlm_ext1_*.root"; ch[1]->Add(myhelper.c_str());

  dataset[2] = "TT1l";
  ch[2] = new TChain("t");
  myhelper = babylocation + "ttbar_1ltop_mgmlm_ext1_*.root"; ch[2]->Add(myhelper.c_str());
  myhelper = babylocation + "ttbar_1ltbr_mgmlm_ext1_*.root"; ch[2]->Add(myhelper.c_str());

  dataset[3] = "Wjets";
  ch[3] = new TChain("t");
  myhelper = babylocation + "wjets_*.root";             ch[3]->Add(myhelper.c_str());
  //myhelper = babylocation + "wgjets_incl_mgmlm_*.root"; ch[3]->Add(myhelper.c_str());
  
  dataset[4] = "Zjets";
  ch[4] = new TChain("t");
  myhelper = babylocation + "dy_m50_mgmlm_*.root";   ch[4]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m1050_mgmlm_*.root"; ch[4]->Add(myhelper.c_str());

  dataset[5] = "WZ";
  ch[5] = new TChain("t");
  myhelper = babylocation + "wz_1l3n_amcnlo_*.root"; ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "wz_3lnu_powheg_*.root"; ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "wz_lnqq_amcnlo_*.root"; ch[5]->Add(myhelper.c_str());

  dataset[6] = "WW";
  ch[6] = new TChain("t");
  myhelper = babylocation + "ww_2l2nu_powheg_*.root"; ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "ww_lnuqq_powheg_*.root"; ch[6]->Add(myhelper.c_str());
  //myhelper = babylocation + "wmwm_powheg_*.root"; ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wpwpjj_ewk-qcd_madgraph_*.root"; ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "ww_2l2nu_dbl_scat_*.root"; ch[6]->Add(myhelper.c_str());

  dataset[7] = "ZZ";
  ch[7] = new TChain("t");
  myhelper = babylocation + "zz_2l2n_powheg_*.root"; ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "zz_2l2q_powheg_*.root"; ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "zz_4l_powheg_*.root"; ch[7]->Add(myhelper.c_str());

  dataset[8] = "ttV";
  ch[8] = new TChain("t");
  myhelper = babylocation + "ttg_incl_amcnlo_*.root";  ch[8]->Add(myhelper.c_str());
  myhelper = babylocation + "tth_bb_powheg_*.root";    ch[8]->Add(myhelper.c_str());
  myhelper = babylocation + "tth_nonbb_powheg_*.root"; ch[8]->Add(myhelper.c_str());
  myhelper = babylocation + "ttw_incl_mgmlm_*.root";   ch[8]->Add(myhelper.c_str());
  myhelper = babylocation + "ttz_incl_mgmlm_*.root";   ch[8]->Add(myhelper.c_str());
  //myhelper = babylocation + "tzq_ll_amcnlo_*.root";    ch[8]->Add(myhelper.c_str());

  dataset[9] = "singleT";
  ch[9] = new TChain("t");
  myhelper = babylocation + "sttw_*.root"; ch[9]->Add(myhelper.c_str());
  myhelper = babylocation + "stt_*.root";  ch[9]->Add(myhelper.c_str());

  dataset[10] = "WWWv2";
  ch[10] = new TChain("t");
  myhelper = babylocation + "www_incl_amcnlo_*.root"; ch[10]->Add(myhelper.c_str());//240000
  //myhelper = babylocation + "www_2l_mia_*.root"; ch[10]->Add(myhelper.c_str());//91900
  //myhelper = babylocation + "www_2l_ext1_mia_*.root"; ch[10]->Add(myhelper.c_str());//164800

  
  for(int i = 0; i<chainsize; ++i){
    //if(i!=5) continue;
    //if(i==0) continue;
    TChain *mych = ch[i];
    string mydataset = dataset[i];
    cout << "Now entering " << mydataset << endl;
    ScanChain(mych,true,-1,mydataset); 

  }
}
