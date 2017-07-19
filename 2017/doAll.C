{

  
  //gROOT->ProcessLine(".L Check3lCR.C+");
  //gROOT->ProcessLine(".L Nminus1Histos.C+");
  //gROOT->ProcessLine(".L CheckWGrelatedStuff.C+");
  gROOT->ProcessLine(".L CheckSingleEvents.C+");
  const unsigned int chainsize = 17;
  TChain *ch[chainsize];
  string dataset[chainsize];


  //string babylocation = "/hadoop/cms/store/user/haweber/AutoTwopler_babies/WWW/merged/ZMET/WWW_vX/skim/";
  //string babylocation = "/nfs-7/userdata/bhashemi/WWW_babies/WWW_v0.1.4/skim/";
  string babylocation = "/nfs-7/userdata/bhashemi/WWW_babies/WWW_v0.1.9/skim/";
  string myhelper;

  dataset[0] = "WWWv2";
  ch[0] = new TChain("t");
  myhelper = babylocation + "www_incl_amcnlo_*.root"; ch[0]->Add(myhelper.c_str());//240000
  dataset[1] = "WWW";
  ch[1] = new TChain("t");
  myhelper = babylocation + "www_2l_mia_*.root";      ch[1]->Add(myhelper.c_str());//91900
  myhelper = babylocation + "www_2l_ext1_mia_*.root"; ch[1]->Add(myhelper.c_str());//164800

  dataset[2] = "VVV";
  ch[2] = new TChain("t");
  myhelper = babylocation + "wwz_incl_amcnlo*.root"; ch[2]->Add(myhelper.c_str());
  myhelper = babylocation + "wzz_incl_amcnlo*.root"; ch[2]->Add(myhelper.c_str());
  myhelper = babylocation + "zzz_incl_amcnlo*.root"; ch[2]->Add(myhelper.c_str());

  dataset[3] = "tt1l";
  ch[3] = new TChain("t");
  myhelper = babylocation + "ttbar_1ltbr_mgmlm_ext1*.root"; ch[3]->Add(myhelper.c_str());
  myhelper = babylocation + "ttbar_1ltop_mgmlm_ext1*.root"; ch[3]->Add(myhelper.c_str());

  dataset[4] = "tt2l";
  ch[4] = new TChain("t");
  myhelper = babylocation + "ttbar_dilep_mgmlm_ext1*.root"; ch[4]->Add(myhelper.c_str());

  dataset[5] = "singleTop";
  ch[5] = new TChain("t");
  myhelper = babylocation + "stt_antitop_incdec_powheg*.root";          ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "stt_top_incdec_powheg*.root";              ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "sttw_antitop_nofullhaddecay_powheg*.root"; ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "sttw_top_nofullhaddecay_powheg*.root";     ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "sttwll_madgraph*.root";                    ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "sts_4f_leptonic_amcnlo*.root";             ch[5]->Add(myhelper.c_str());
  myhelper = babylocation + "tzq_ll_amcnlo*.root";                      ch[5]->Add(myhelper.c_str());
  
  dataset[6] = "Wjets";
  ch[6] = new TChain("t");
  myhelper = babylocation + "wjets_incl_mgmlm*.root";                ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wjets_ht100_mgmlm_ext1*.root";          ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wjets_ht200_mgmlm_ext1*.root";          ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wjets_ht400_mgmlm_ext1*.root";          ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wjets_ht600_mgmlm_ext1*.root";          ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wjets_ht800_mgmlm_ext1*.root";          ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wjets_ht1200_mgmlm_nonext*.root";       ch[6]->Add(myhelper.c_str());
  myhelper = babylocation + "wjets_ht2500_mgmlm_ext1*.root";         ch[6]->Add(myhelper.c_str());

  dataset[7] = "Zjets";
  ch[7] = new TChain("t");
  myhelper = babylocation + "dy_m1050_mgmlm*.root";             ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ext1*.root";          ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ht100_ext1*.root";    ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ht200_ext1*.root";    ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ht400_ext1*.root";    ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ht600_nonext*.root";  ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ht800_nonext*.root";  ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ht1200_nonext*.root"; ch[7]->Add(myhelper.c_str());
  myhelper = babylocation + "dy_m50_mgmlm_ht2500_nonext*.root"; ch[7]->Add(myhelper.c_str());

  dataset[8] = "WW";
  ch[8] = new TChain("t");
  myhelper = babylocation + "wpwpjj_ewk-qcd_madgraph*.root";  ch[8]->Add(myhelper.c_str());
  myhelper = babylocation + "ww_2l2nu_dbl_scat*.root";        ch[8]->Add(myhelper.c_str());
  myhelper = babylocation + "ww_2l2nu_powheg*.root";          ch[8]->Add(myhelper.c_str());
  myhelper = babylocation + "ww_lnuqq_powheg*.root";          ch[8]->Add(myhelper.c_str());

  dataset[9] = "WZ";
  ch[9] = new TChain("t");
  myhelper = babylocation + "wz_1l3n_amcnlo*.root";          ch[9]->Add(myhelper.c_str());
  myhelper = babylocation + "wz_3lnu_powheg*.root";          ch[9]->Add(myhelper.c_str());
  myhelper = babylocation + "wz_lnqq_amcnlo*.root";          ch[9]->Add(myhelper.c_str());

  dataset[10] = "ZZ";
  ch[10] = new TChain("t");
  myhelper = babylocation + "zz_2l2n_powheg*.root";        ch[10]->Add(myhelper.c_str());
  myhelper = babylocation + "zz_2l2q_powheg*.root";        ch[10]->Add(myhelper.c_str());
  myhelper = babylocation + "zz_2q2n_amcnlo*.root";        ch[10]->Add(myhelper.c_str());
  myhelper = babylocation + "zz_4l_powheg*.root";          ch[10]->Add(myhelper.c_str());

  dataset[11] = "ttV";
  ch[11] = new TChain("t");
  myhelper = babylocation + "ttg_incl_amcnlo*.root";        ch[11]->Add(myhelper.c_str());
  myhelper = babylocation + "tth_bb_powheg*.root";          ch[11]->Add(myhelper.c_str());
  myhelper = babylocation + "tth_nonbb_powheg*.root";       ch[11]->Add(myhelper.c_str());
  myhelper = babylocation + "ttw_incl_mgmlm*.root";         ch[11]->Add(myhelper.c_str());
  myhelper = babylocation + "ttz_incl_mgmlm*.root";         ch[11]->Add(myhelper.c_str());

  dataset[12] = "Other";
  ch[12] = new TChain("t");
  myhelper = babylocation + "vh_nonbb_amcnlo*.root";          ch[12]->Add(myhelper.c_str());

  dataset[13] = "WG";
  ch[13] = new TChain("t");
  myhelper = babylocation + "wgjets_incl_mgmlm_skim_1.root";          ch[13]->Add(myhelper.c_str());

  dataset[14] = "ZG";
  ch[14] = new TChain("t");
  myhelper = babylocation + "zgamma_2lG_amc_skim_1.root";          ch[14]->Add(myhelper.c_str());
  
  dataset[15] = "Background";
  ch[15] = new TChain("t");
  ch[15]->Add(ch[2]);
  ch[15]->Add(ch[3]);
  ch[15]->Add(ch[4]);
  ch[15]->Add(ch[5]);
  ch[15]->Add(ch[6]);
  ch[15]->Add(ch[7]);
  ch[15]->Add(ch[8]);
  ch[15]->Add(ch[9]);
  ch[15]->Add(ch[10]);
  ch[15]->Add(ch[11]);
  ch[15]->Add(ch[12]);
  ch[15]->Add(ch[13]);
  ch[15]->Add(ch[14]);

  babylocation = "/nfs-7/userdata/bhashemi/WWW_babies/WWW_v0.1.11/skim/";
  dataset[16] = "Data";
  ch[16] = new TChain("t");
  myhelper = babylocation + "data_*ee*.root";          ch[16]->Add(myhelper.c_str());
  myhelper = babylocation + "data_*em*.root";          ch[16]->Add(myhelper.c_str());
  myhelper = babylocation + "data_*mm*.root";          ch[16]->Add(myhelper.c_str());
  
  for(int i = 0; i<chainsize; ++i){
    //if(i<=14) continue;
    //if(i!=13&&i!=14) continue;
    if(i!=16) continue;
    //if(i==0) continue;
    TChain *mych = ch[i];
    string mydataset = dataset[i];
    cout << "Now entering " << mydataset << endl;
    ScanChain(mych,true,-1,mydataset); 

  }
}
