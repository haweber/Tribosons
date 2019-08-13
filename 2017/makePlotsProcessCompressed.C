void makePlotsProcessCompressed(){
    
    bool data = false;
    bool twosig = true;
    gStyle->SetOptStat(0);
    bool logy = false;
    map<string, TH1D*> h;
    map<string, TH1D*> rat;
    map<string, THStack*> stack;
    vector<string> histonames;
    vector<string> axisnames;
    //vector<Color_t> cols;
    vector<int> cols;
    vector<int> plottogether;
    vector<string> samples;
    vector<string> samplesleg;
    vector<int> onlySSSFOS;
    
    
    TFile *f = TFile::Open("Nminus1Histos.root");
    histonames.push_back("Detajj_SRlike_allSS");         axisnames.push_back("#Delta#eta_{jj}");                 onlySSSFOS.push_back(1);
    histonames.push_back("Detajj_SRlike_SSee");          axisnames.push_back("#Delta#eta_{jj}");                 onlySSSFOS.push_back(1);
    histonames.push_back("Detajj_SRlike_SSem");          axisnames.push_back("#Delta#eta_{jj}");                 onlySSSFOS.push_back(1);
    histonames.push_back("Detajj_SRlike_SSmm");          axisnames.push_back("#Delta#eta_{jj}");                 onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRlike_allSS");            axisnames.push_back("M_{jj} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRlike_SSee");             axisnames.push_back("M_{jj} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRlike_SSem");             axisnames.push_back("M_{jj} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRlike_SSmm");             axisnames.push_back("M_{jj} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRlike_allSS");            axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRlike_SSee");             axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRlike_SSem");             axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRlike_SSmm");             axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRlike_allSS");            axisnames.push_back("M_{ll} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRlike_SSee");             axisnames.push_back("M_{ll} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRlike_SSem");             axisnames.push_back("M_{ll} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRlike_SSmm");             axisnames.push_back("M_{ll} [GeV]");                    onlySSSFOS.push_back(1);
    histonames.push_back("MTmax_SRlike_SSem");           axisnames.push_back("M_{T}^{max} [GeV]");               onlySSSFOS.push_back(1);
    
    histonames.push_back("MET_SRlike_allSFOS");          axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(2);
    histonames.push_back("MET_SRlike_0SFOS");            axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(2);
    histonames.push_back("MET_SRlike_1SFOS");            axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(2);
    histonames.push_back("MET_SRlike_2SFOS");            axisnames.push_back("E_{T}^{miss} [GeV]");              onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_SRlike_allSFOS");        axisnames.push_back("M_{SFOS} [GeV]");                  onlySSSFOS.push_back(2);
    histonames.push_back("MSF_SRlike_0SFOS");            axisnames.push_back("M_{SF} [GeV]");                    onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_SRlike_1SFOS");          axisnames.push_back("M_{SFOS} [GeV]");                  onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_SRlike_2SFOS");          axisnames.push_back("M_{SFOS} [GeV]");                  onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRlike_allSFOS");        axisnames.push_back("p_{T}(lll) [GeV]");                onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRlike_0SFOS");          axisnames.push_back("p_{T}(lll) [GeV]");                onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRlike_1SFOS");          axisnames.push_back("p_{T}(lll) [GeV]");                onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRlike_2SFOS");          axisnames.push_back("p_{T}(lll) [GeV]");                onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRlike_allSFOS");   axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");    onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRlike_0SFOS");     axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");    onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRlike_1SFOS");     axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");    onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRlike_2SFOS");     axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");    onlySSSFOS.push_back(2);
    
    histonames.push_back("Detajj_SRloose_allSS");         axisnames.push_back("#Delta#eta_{jj}");                onlySSSFOS.push_back(1);
    histonames.push_back("Detajj_SRloose_SSee");          axisnames.push_back("#Delta#eta_{jj}");                onlySSSFOS.push_back(1);
    histonames.push_back("Detajj_SRloose_SSem");          axisnames.push_back("#Delta#eta_{jj}");                onlySSSFOS.push_back(1);
    histonames.push_back("Detajj_SRloose_SSmm");          axisnames.push_back("#Delta#eta_{jj}");                onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRloose_allSS");            axisnames.push_back("M_{jj} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRloose_SSee");             axisnames.push_back("M_{jj} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRloose_SSem");             axisnames.push_back("M_{jj} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_SRloose_SSmm");             axisnames.push_back("M_{jj} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRloose_allSS");            axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRloose_SSee");             axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRloose_SSem");             axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(1);
    histonames.push_back("MET_SRloose_SSmm");             axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRloose_allSS");            axisnames.push_back("M_{ll} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRloose_SSee");             axisnames.push_back("M_{ll} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRloose_SSem");             axisnames.push_back("M_{ll} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("Mll_SRloose_SSmm");             axisnames.push_back("M_{ll} [GeV]");                   onlySSSFOS.push_back(1);
    histonames.push_back("MTmax_SRloose_SSem");           axisnames.push_back("M_{T}^{max} [GeV]");              onlySSSFOS.push_back(1);
    
    histonames.push_back("MET_SRloose_allSFOS");          axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(2);
    histonames.push_back("MET_SRloose_0SFOS");            axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(2);
    histonames.push_back("MET_SRloose_1SFOS");            axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(2);
    histonames.push_back("MET_SRloose_2SFOS");            axisnames.push_back("E_{T}^{miss} [GeV]");             onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_SRloose_allSFOS");        axisnames.push_back("M_{SFOS} [GeV]");                 onlySSSFOS.push_back(2);
    histonames.push_back("MSF_SRloose_0SFOS");            axisnames.push_back("M_{SF} [GeV]");                   onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_SRloose_1SFOS");          axisnames.push_back("M_{SFOS} [GeV]");                 onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_SRloose_2SFOS");          axisnames.push_back("M_{SFOS} [GeV]");                 onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRloose_allSFOS");        axisnames.push_back("p_{T}(lll) [GeV]");               onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRloose_0SFOS");          axisnames.push_back("p_{T}(lll) [GeV]");               onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRloose_1SFOS");          axisnames.push_back("p_{T}(lll) [GeV]");               onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_SRloose_2SFOS");          axisnames.push_back("p_{T}(lll) [GeV]");               onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRloose_allSFOS");   axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");   onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRloose_0SFOS");     axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");   onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRloose_1SFOS");     axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");   onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_SRloose_2SFOS");     axisnames.push_back("#Delta#phi(lll,E_{T}^{miss})");   onlySSSFOS.push_back(2);
/*
    TFile *f = TFile::Open("ExcessPlots.root");
    histonames.push_back("SR_lMET");                   axisnames.push_back("SR | E_{T}^{miss} < 40 GeV, M_{jj} sideband");                                      onlySSSFOS.push_back(1);
    histonames.push_back("SR_hMET");                   axisnames.push_back("SR | E_{T}^{miss} > 40 GeV, M_{jj} sideband");                                      onlySSSFOS.push_back(1);
    histonames.push_back("SR_iMET");                   axisnames.push_back("SR | M_{jj} sideband");                                                             onlySSSFOS.push_back(1);
    histonames.push_back("SRpresel_lMET");             axisnames.push_back("SR (preselection) | E_{T}^{miss} < 40 GeV, M_{jj} sideband");                       onlySSSFOS.push_back(1);
    histonames.push_back("SRpresel_hMET");             axisnames.push_back("SR (preselection) | E_{T}^{miss} > 40 GeV, M_{jj} sideband");                       onlySSSFOS.push_back(1);
    histonames.push_back("SRpresel_iMET");             axisnames.push_back("SR (preselection) | M_{jj} sideband");                                              onlySSSFOS.push_back(1);
    
    histonames.push_back("SR_ee");                     axisnames.push_back("SR - e^{#pm}e^{#pm} | splitted by run range, M_{jj} sideband");                     onlySSSFOS.push_back(1);
    histonames.push_back("SR_em");                     axisnames.push_back("SR - e^{#pm}#mu^{#pm} | splitted by run range, M_{jj} sideband");                   onlySSSFOS.push_back(1);
    histonames.push_back("SR_mm");                     axisnames.push_back("SR - #mu^{#pm}#mu^{#pm} | splitted by run range, M_{jj} sideband");                 onlySSSFOS.push_back(1);
    histonames.push_back("SRpresel_ee");               axisnames.push_back("SR (preselection) - e^{#pm}e^{#pm} | splitted by run range, M_{jj} sideband");      onlySSSFOS.push_back(1);
    histonames.push_back("SRpresel_em");               axisnames.push_back("SR (preselection) - e^{#pm}#mu^{#pm} | splitted by run range, M_{jj} sideband");    onlySSSFOS.push_back(1);
    histonames.push_back("SRpresel_mm");               axisnames.push_back("SR (preselection) - #mu^{#pm}#mu^{#pm} | splitted by run range, M_{jj} sideband");  onlySSSFOS.push_back(1);
    
    histonames.push_back("NB_ee_lMET");                axisnames.push_back("N_{B} - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, M_{jj} sideband");                  onlySSSFOS.push_back(1);
    histonames.push_back("NB_ee_hMET");                axisnames.push_back("N_{B} - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, M_{jj} sideband");                  onlySSSFOS.push_back(1);
    histonames.push_back("NB_ee_iMET");                axisnames.push_back("N_{B} - e^{#pm}e^{#pm} | M_{jj} sideband");                                         onlySSSFOS.push_back(1);
    histonames.push_back("NB_em_lMET");                axisnames.push_back("N_{B} - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, M_{jj} sideband");                onlySSSFOS.push_back(1);
    histonames.push_back("NB_em_hMET");                axisnames.push_back("N_{B} - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, M_{jj} sideband");                onlySSSFOS.push_back(1);
    histonames.push_back("NB_em_iMET");                axisnames.push_back("N_{B} - e^{#pm}#mu^{#pm} | M_{jj} sideband");                                       onlySSSFOS.push_back(1);
    histonames.push_back("NB_mm_lMET");                axisnames.push_back("N_{B} - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("NB_mm_hMET");                axisnames.push_back("N_{B} - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("NB_mm_iMET");                axisnames.push_back("N_{B} - #mu^{#pm}#mu^{#pm} | M_{jj} sideband");                                     onlySSSFOS.push_back(1);
    
    histonames.push_back("MET_ee_0b");                 axisnames.push_back("E_{T}^{miss} [GeV] - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                        onlySSSFOS.push_back(1);
    histonames.push_back("MET_em_0b");                 axisnames.push_back("E_{T}^{miss} [GeV] - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                      onlySSSFOS.push_back(1);
    histonames.push_back("MET_mm_0b");                 axisnames.push_back("E_{T}^{miss} [GeV] - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                    onlySSSFOS.push_back(1);
    histonames.push_back("MTmax_ee_0b");               axisnames.push_back("M_{T}^{max} [GeV] - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                         onlySSSFOS.push_back(1);
    histonames.push_back("MTmax_em_0b");               axisnames.push_back("M_{T}^{max} [GeV] - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                       onlySSSFOS.push_back(1);
    histonames.push_back("MTmax_mm_0b");               axisnames.push_back("M_{T}^{max} [GeV] - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                     onlySSSFOS.push_back(1);
    
    histonames.push_back("q_ee_0b_lMET");              axisnames.push_back("charge - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");             onlySSSFOS.push_back(1);
    histonames.push_back("q_ee_0b_hMET");              axisnames.push_back("charge - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");             onlySSSFOS.push_back(1);
    histonames.push_back("q_ee_0b_iMET");              axisnames.push_back("charge - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                                    onlySSSFOS.push_back(1);
    histonames.push_back("q_em_0b_lMET");              axisnames.push_back("charge - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("q_em_0b_hMET");              axisnames.push_back("charge - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("q_em_0b_iMET");              axisnames.push_back("charge - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                  onlySSSFOS.push_back(1);
    histonames.push_back("q_mm_0b_lMET");              axisnames.push_back("charge - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");         onlySSSFOS.push_back(1);
    histonames.push_back("q_mm_0b_hMET");              axisnames.push_back("charge - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");         onlySSSFOS.push_back(1);
    histonames.push_back("q_mm_0b_iMET");              axisnames.push_back("charge - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                onlySSSFOS.push_back(1);
    
    histonames.push_back("NJ_ee_0b_lMET");             axisnames.push_back("N_{J} - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("NJ_ee_0b_hMET");             axisnames.push_back("N_{J} - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("NJ_ee_0b_iMET");             axisnames.push_back("N_{J} - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                                     onlySSSFOS.push_back(1);
    histonames.push_back("NJ_em_0b_lMET");             axisnames.push_back("N_{J} - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("NJ_em_0b_hMET");             axisnames.push_back("N_{J} - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("NJ_em_0b_iMET");             axisnames.push_back("N_{J} - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                   onlySSSFOS.push_back(1);
    histonames.push_back("NJ_mm_0b_lMET");             axisnames.push_back("N_{J} - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("NJ_mm_0b_hMET");             axisnames.push_back("N_{J} - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("NJ_mm_0b_iMET");             axisnames.push_back("N_{J} - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                 onlySSSFOS.push_back(1);
    
    histonames.push_back("Mll_ee_0b_lMET");            axisnames.push_back("M_{ll} [GeV] - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");       onlySSSFOS.push_back(1);
    histonames.push_back("Mll_ee_0b_hMET");            axisnames.push_back("M_{ll} [GeV] - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");       onlySSSFOS.push_back(1);
    histonames.push_back("Mll_ee_0b_iMET");            axisnames.push_back("M_{ll} [GeV] - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                              onlySSSFOS.push_back(1);
    histonames.push_back("Mll_em_0b_lMET");            axisnames.push_back("M_{ll} [GeV] - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("Mll_em_0b_hMET");            axisnames.push_back("M_{ll} [GeV] - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("Mll_em_0b_iMET");            axisnames.push_back("M_{ll} [GeV] - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                            onlySSSFOS.push_back(1);
    histonames.push_back("Mll_mm_0b_lMET");            axisnames.push_back("M_{ll} [GeV] - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("Mll_mm_0b_hMET");            axisnames.push_back("M_{ll} [GeV] - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("Mll_mm_0b_iMET");            axisnames.push_back("M_{ll} [GeV] - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                          onlySSSFOS.push_back(1);
    
    histonames.push_back("Mjj_ee_0b_lMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_ee_0b_hMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_ee_0b_iMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                            onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_em_0b_lMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_em_0b_hMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_em_0b_iMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                          onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_mm_0b_lMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband"); onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_mm_0b_hMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband"); onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_mm_0b_iMET");            axisnames.push_back("M_{jj} [GeV] (closest #DeltaR) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                        onlySSSFOS.push_back(1);
    
    histonames.push_back("MjjL_ee_0b_lMET");           axisnames.push_back("M_{jj} [GeV] (leading) - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_ee_0b_hMET");           axisnames.push_back("M_{jj} [GeV] (leading) - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_ee_0b_iMET");           axisnames.push_back("M_{jj} [GeV] (leading) - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                            onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_em_0b_lMET");           axisnames.push_back("M_{jj} [GeV] (leading) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_em_0b_hMET");           axisnames.push_back("M_{jj} [GeV] (leading) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_em_0b_iMET");           axisnames.push_back("M_{jj} [GeV] (leading) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                          onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_mm_0b_lMET");           axisnames.push_back("M_{jj} [GeV] (leading) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband"); onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_mm_0b_hMET");           axisnames.push_back("M_{jj} [GeV] (leading) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband"); onlySSSFOS.push_back(1);
    histonames.push_back("MjjL_mm_0b_iMET");           axisnames.push_back("M_{jj} [GeV] (leading) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                        onlySSSFOS.push_back(1);
    
    histonames.push_back("pTll_ee_0b_lMET");           axisnames.push_back("p_{T}(ll) [GeV] - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("pTll_ee_0b_hMET");           axisnames.push_back("p_{T}(ll) [GeV] - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("pTll_ee_0b_iMET");           axisnames.push_back("p_{T}(ll) [GeV] - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                                   onlySSSFOS.push_back(1);
    histonames.push_back("pTll_em_0b_lMET");           axisnames.push_back("p_{T}(ll) [GeV] - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("pTll_em_0b_hMET");           axisnames.push_back("p_{T}(ll) [GeV] - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("pTll_em_0b_iMET");           axisnames.push_back("p_{T}(ll) [GeV] - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                 onlySSSFOS.push_back(1);
    histonames.push_back("pTll_mm_0b_lMET");           axisnames.push_back("p_{T}(ll) [GeV] - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("pTll_mm_0b_hMET");           axisnames.push_back("p_{T}(ll) [GeV] - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("pTll_mm_0b_iMET");           axisnames.push_back("p_{T}(ll) [GeV] - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                               onlySSSFOS.push_back(1);
    
    histonames.push_back("dPhillMET_ee_0b_lMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");    onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_ee_0b_hMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");    onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_ee_0b_iMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                           onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_em_0b_lMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");  onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_em_0b_hMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");  onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_em_0b_iMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                         onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_mm_0b_lMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_mm_0b_hMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");onlySSSFOS.push_back(1);
    histonames.push_back("dPhillMET_mm_0b_iMET");      axisnames.push_back("#Delta#phi(ll,E_{T}^{miss}) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                       onlySSSFOS.push_back(1);
    
    histonames.push_back("dPhill_ee_0b_lMET");         axisnames.push_back("#Delta#phi(l,l) - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_ee_0b_hMET");         axisnames.push_back("#Delta#phi(l,l) - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_ee_0b_iMET");         axisnames.push_back("#Delta#phi(l,l) - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                                   onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_em_0b_lMET");         axisnames.push_back("#Delta#phi(l,l) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_em_0b_hMET");         axisnames.push_back("#Delta#phi(l,l) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_em_0b_iMET");         axisnames.push_back("#Delta#phi(l,l) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                 onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_mm_0b_lMET");         axisnames.push_back("#Delta#phi(l,l) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_mm_0b_hMET");         axisnames.push_back("#Delta#phi(l,l) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("dPhill_mm_0b_iMET");         axisnames.push_back("#Delta#phi(l,l) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                               onlySSSFOS.push_back(1);
    
    histonames.push_back("dEtall_ee_0b_lMET");         axisnames.push_back("#Delta#eta(l,l) - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_ee_0b_hMET");         axisnames.push_back("#Delta#eta(l,l) - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_ee_0b_iMET");         axisnames.push_back("#Delta#eta(l,l) - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                                   onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_em_0b_lMET");         axisnames.push_back("#Delta#eta(l,l) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_em_0b_hMET");         axisnames.push_back("#Delta#eta(l,l) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_em_0b_iMET");         axisnames.push_back("#Delta#eta(l,l) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                 onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_mm_0b_lMET");         axisnames.push_back("#Delta#eta(l,l) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_mm_0b_hMET");         axisnames.push_back("#Delta#eta(l,l) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("dEtall_mm_0b_iMET");         axisnames.push_back("#Delta#eta(l,l) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                               onlySSSFOS.push_back(1);
    
    histonames.push_back("minDRlj_ee_0b_lMET");        axisnames.push_back("min#DeltaR(l,j) - e^{#pm}e^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_ee_0b_hMET");        axisnames.push_back("min#DeltaR(l,j) - e^{#pm}e^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");            onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_ee_0b_iMET");        axisnames.push_back("min#DeltaR(l,j) - e^{#pm}e^{#pm} | 0b, M_{jj} sideband");                                   onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_em_0b_lMET");        axisnames.push_back("min#DeltaR(l,j) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_em_0b_hMET");        axisnames.push_back("min#DeltaR(l,j) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");          onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_em_0b_iMET");        axisnames.push_back("min#DeltaR(l,j) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                 onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_mm_0b_lMET");        axisnames.push_back("min#DeltaR(l,j) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_mm_0b_hMET");        axisnames.push_back("min#DeltaR(l,j) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("minDRlj_mm_0b_iMET");        axisnames.push_back("min#DeltaR(l,j) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                               onlySSSFOS.push_back(1);
    
    histonames.push_back("MuPt_em_0b_lMET");           axisnames.push_back("p_{T}(#mu) [GeV] - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");         onlySSSFOS.push_back(1);
    histonames.push_back("MuPt_em_0b_hMET");           axisnames.push_back("p_{T}(#mu) [GeV] - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");         onlySSSFOS.push_back(1);
    histonames.push_back("MuPt_em_0b_iMET");           axisnames.push_back("p_{T}(#mu) [GeV] - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                onlySSSFOS.push_back(1);
    histonames.push_back("MuPt_mm_0b_lMET");           axisnames.push_back("p_{T}(#mu) [GeV] - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");       onlySSSFOS.push_back(1);
    histonames.push_back("MuPt_mm_0b_hMET");           axisnames.push_back("p_{T}(#mu) [GeV] - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");       onlySSSFOS.push_back(1);
    histonames.push_back("MuPt_mm_0b_iMET");           axisnames.push_back("p_{T}(#mu) [GeV] - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                              onlySSSFOS.push_back(1);
    
    histonames.push_back("MuEta_em_0b_lMET");          axisnames.push_back("#eta(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");                onlySSSFOS.push_back(1);
    histonames.push_back("MuEta_em_0b_hMET");          axisnames.push_back("#eta(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");                onlySSSFOS.push_back(1);
    histonames.push_back("MuEta_em_0b_iMET");          axisnames.push_back("#eta(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                       onlySSSFOS.push_back(1);
    histonames.push_back("MuEta_mm_0b_lMET");          axisnames.push_back("#eta(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("MuEta_mm_0b_hMET");          axisnames.push_back("#eta(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("MuEta_mm_0b_iMET");          axisnames.push_back("#eta(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                     onlySSSFOS.push_back(1);
    
    histonames.push_back("MuPhi_em_0b_lMET");          axisnames.push_back("#phi(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");                onlySSSFOS.push_back(1);
    histonames.push_back("MuPhi_em_0b_hMET");          axisnames.push_back("#phi(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");                onlySSSFOS.push_back(1);
    histonames.push_back("MuPhi_em_0b_iMET");          axisnames.push_back("#phi(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                       onlySSSFOS.push_back(1);
    histonames.push_back("MuPhi_mm_0b_lMET");          axisnames.push_back("#phi(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("MuPhi_mm_0b_hMET");          axisnames.push_back("#phi(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");              onlySSSFOS.push_back(1);
    histonames.push_back("MuPhi_mm_0b_iMET");          axisnames.push_back("#phi(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                     onlySSSFOS.push_back(1);
    
    histonames.push_back("MuRelIso_em_0b_lMET");       axisnames.push_back("I_{rel}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");             onlySSSFOS.push_back(1);
    histonames.push_back("MuRelIso_em_0b_hMET");       axisnames.push_back("I_{rel}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");             onlySSSFOS.push_back(1);
    histonames.push_back("MuRelIso_em_0b_iMET");       axisnames.push_back("I_{rel}(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                    onlySSSFOS.push_back(1);
    histonames.push_back("MuRelIso_mm_0b_lMET");       axisnames.push_back("I_{rel}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("MuRelIso_mm_0b_hMET");       axisnames.push_back("I_{rel}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("MuRelIso_mm_0b_iMET");       axisnames.push_back("I_{rel}(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                  onlySSSFOS.push_back(1);
    
    histonames.push_back("MuIP3D_em_0b_lMET");         axisnames.push_back("IP_{3D}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");             onlySSSFOS.push_back(1);
    histonames.push_back("MuIP3D_em_0b_hMET");         axisnames.push_back("IP_{3D}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");             onlySSSFOS.push_back(1);
    histonames.push_back("MuIP3D_em_0b_iMET");         axisnames.push_back("IP_{3D}(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                    onlySSSFOS.push_back(1);
    histonames.push_back("MuIP3D_mm_0b_lMET");         axisnames.push_back("IP_{3D}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("MuIP3D_mm_0b_hMET");         axisnames.push_back("IP_{3D}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("MuIP3D_mm_0b_iMET");         axisnames.push_back("IP_{3D}(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                  onlySSSFOS.push_back(1);
    
    histonames.push_back("MuPtErrOvPt_em_0b_lMET");    axisnames.push_back("#sigma(p_{T})/p_{T}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("MuPtErrOvPt_em_0b_hMET");    axisnames.push_back("#sigma(p_{T})/p_{T}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");     onlySSSFOS.push_back(1);
    histonames.push_back("MuPtErrOvPt_em_0b_iMET");    axisnames.push_back("#sigma(p_{T})/p_{T}(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                            onlySSSFOS.push_back(1);
    histonames.push_back("MuPtErrOvPt_mm_0b_lMET");    axisnames.push_back("#sigma(p_{T})/p_{T}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("MuPtErrOvPt_mm_0b_hMET");    axisnames.push_back("#sigma(p_{T})/p_{T}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");   onlySSSFOS.push_back(1);
    histonames.push_back("MuPtErrOvPt_mm_0b_iMET");    axisnames.push_back("#sigma(p_{T})/p_{T}(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                          onlySSSFOS.push_back(1);
    
    histonames.push_back("MuValidFrac_em_0b_lMET");    axisnames.push_back("valid fraction(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");      onlySSSFOS.push_back(1);
    histonames.push_back("MuValidFrac_em_0b_hMET");    axisnames.push_back("valid fraction(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");      onlySSSFOS.push_back(1);
    histonames.push_back("MuValidFrac_em_0b_iMET");    axisnames.push_back("valid fraction(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                             onlySSSFOS.push_back(1);
    histonames.push_back("MuValidFrac_mm_0b_lMET");    axisnames.push_back("valid fraction(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");    onlySSSFOS.push_back(1);
    histonames.push_back("MuValidFrac_mm_0b_hMET");    axisnames.push_back("valid fraction(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");    onlySSSFOS.push_back(1);
    histonames.push_back("MuValidFrac_mm_0b_iMET");    axisnames.push_back("valid fraction(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                           onlySSSFOS.push_back(1);
    
    histonames.push_back("MuChi2N_em_0b_lMET");        axisnames.push_back("#chi^{2}_{N}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("MuChi2N_em_0b_hMET");        axisnames.push_back("#chi^{2}_{N}(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");        onlySSSFOS.push_back(1);
    histonames.push_back("MuChi2N_em_0b_iMET");        axisnames.push_back("#chi^{2}_{N}(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                               onlySSSFOS.push_back(1);
    histonames.push_back("MuChi2N_mm_0b_lMET");        axisnames.push_back("#chi^{2}_{N}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");      onlySSSFOS.push_back(1);
    histonames.push_back("MuChi2N_mm_0b_hMET");        axisnames.push_back("#chi^{2}_{N}(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");      onlySSSFOS.push_back(1);
    histonames.push_back("MuChi2N_mm_0b_iMET");        axisnames.push_back("#chi^{2}_{N}(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                             onlySSSFOS.push_back(1);
    
    histonames.push_back("MuLostHits_em_0b_lMET");     axisnames.push_back("lost hits(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("MuLostHits_em_0b_hMET");     axisnames.push_back("lost hits(#mu) - e^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");           onlySSSFOS.push_back(1);
    histonames.push_back("MuLostHits_em_0b_iMET");     axisnames.push_back("lost hits(#mu) - e^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                  onlySSSFOS.push_back(1);
    histonames.push_back("MuLostHits_mm_0b_lMET");     axisnames.push_back("lost hits(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} < 40 GeV, 0b, M_{jj} sideband");         onlySSSFOS.push_back(1);
    histonames.push_back("MuLostHits_mm_0b_hMET");     axisnames.push_back("lost hits(#mu) - #mu^{#pm}#mu^{#pm} | E_{T}^{miss} > 40 GeV, 0b, M_{jj} sideband");         onlySSSFOS.push_back(1);
    histonames.push_back("MuLostHits_mm_0b_iMET");     axisnames.push_back("lost hits(#mu) - #mu^{#pm}#mu^{#pm} | 0b, M_{jj} sideband");                                onlySSSFOS.push_back(1);
     */
     /*
    TFile *f = TFile::Open("Check3lCRv2.root");
    
    histonames.push_back("YieldsSR");                          axisnames.push_back("signal regions");                                                           onlySSSFOS.push_back(0);
    histonames.push_back("YieldsCR_SSany_dropMjj");            axisnames.push_back("control regions");                                                          onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_lowMET_SRlike_allSS");           axisnames.push_back("M_{jj} [GeV] | SR but low E_{T}^{miss}");                                   onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_lowMET_CRlike_allSS");           axisnames.push_back("M_{jj} [GeV] | CR but low E_{T}^{miss}");                                   onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRlike_allSS");                  axisnames.push_back("M_{jj} [GeV] | CR");                                                        onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRlike_nolowMll_allSS");         axisnames.push_back("M_{jj} [GeV] | CR and no low M_{ll} cut");                                  onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhiPt_1SFOS");            axisnames.push_back("M_{ll} [GeV] | 1SFOS, invert #Delta$phi(lll,E_{T}^{miss}) and p_{T}(lll)"); onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhi_1SFOS");              axisnames.push_back("M_{ll} [GeV] | 1SFOS, invert #Delta$phi(lll,E_{T}^{miss})");                onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertPt_1SFOS");                axisnames.push_back("M_{ll} [GeV] | 1SFOS, invert p_{T}(lll)");                                  onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertMET_1SFOS");               axisnames.push_back("M_{ll} [GeV] | 1SFOS, invert E_{T}^{miss}");                                onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertMETdPhiPt_1SFOS");         axisnames.push_back("M_{ll} [GeV] | 1SFOS, invert E_{T}^{miss}, #Delta$phi(lll,E_{T}^{miss}) and p_{T}(lll)");   onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhiPt_2SFOS");            axisnames.push_back("M_{ll} [GeV] | 2SFOS, invert #Delta$phi(lll,E_{T}^{miss}) and p_{T}(lll)"); onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhi_2SFOS");              axisnames.push_back("M_{ll} [GeV] | 2SFOS, invert #Delta$phi(lll,E_{T}^{miss})");                onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertPt_2SFOS");                axisnames.push_back("M_{ll} [GeV] | 2SFOS, invert p_{T}(lll)");                                  onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertMET_2SFOS");               axisnames.push_back("M_{ll} [GeV] | 2SFOS, invert E_{T}^{miss}");                                onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertMETdPhiPt_2SFOS");         axisnames.push_back("M_{ll} [GeV] | 2SFOS, invert E_{T}^{miss}, #Delta$phi(lll,E_{T}^{miss}) and p_{T}(lll)");   onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhiOrPt_1SFOS");          axisnames.push_back("M_{ll} [GeV] | 1SFOS, invert #Delta$phi(lll,E_{T}^{miss}) or p_{T}(lll)");  onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhiOrPt_2SFOS");          axisnames.push_back("M_{ll} [GeV] | 2SFOS, invert #Delta$phi(lll,E_{T}^{miss}) or p_{T}(lll)");  onlySSSFOS.push_back(2);
    
    histonames.push_back("Mll_inverteitherMETdPhiPt_1SFOS");   axisnames.push_back("M_{ll} [GeV] | 1SFOS, invert E_{T}^{miss}, #Delta$phi(lll,E_{T}^{miss}) or p_{T}(lll)");    onlySSSFOS.push_back(2);
    histonames.push_back("Mll_inverteitherMETdPhiPt_2SFOS");   axisnames.push_back("M_{ll} [GeV] | 2SFOS, invert E_{T}^{miss}, #Delta$phi(lll,E_{T}^{miss}) or p_{T}(lll)");    onlySSSFOS.push_back(2);
    histonames.push_back("Mll_SRlike_1SFOS");                  axisnames.push_back("M_{ll} [GeV] | 1SFOS, SR");                                                 onlySSSFOS.push_back(2);
    histonames.push_back("Mll_SRlike_2SFOS");                  axisnames.push_back("M_{ll} [GeV] | 2SFOS, SR");                                                 onlySSSFOS.push_back(2);
    histonames.push_back("Mll_SRlike_MllclosestZ_2SFOS");      axisnames.push_back("M_{ll} [GeV] | 2SFOS, SR, min|M_{ll}-M_{Z}|");                              onlySSSFOS.push_back(2);
    histonames.push_back("Mll_SRlike_MllfurthestZ_2SFOS");     axisnames.push_back("M_{ll} [GeV] | 2SFOS, SR, max|M_{ll}-M_{Z}|");                              onlySSSFOS.push_back(2);
    histonames.push_back("Mll_ge2j_1SFOS");                    axisnames.push_back("M_{ll} [GeV] | 2SFOS, #geq2 jets");                                         onlySSSFOS.push_back(2);
    histonames.push_back("Mll_ge2j_2SFOS");                    axisnames.push_back("M_{ll} [GeV] | 2SFOS, #geq2 jets");                                         onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRlike_allSS_SSpairnotlead2");   axisnames.push_back("M_{jj} [GeV] | CR, SS not 2 lead leps");                                    onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRlike_allSS_SSpairnotlead2v2"); axisnames.push_back("M_{jj} [GeV] | CR, SS not 2 lead leps, M(ll) on lead leps");                onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_SRlike_allSS");                  axisnames.push_back("M_{jj} [GeV] | SR");                                                        onlySSSFOS.push_back(1);
    histonames.push_back("Mjj_CRlike_SSany_allSS");            axisnames.push_back("M_{jj} [GeV] | CR-loose");                                                  onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRlike_SSany_allSS_jesup");      axisnames.push_back("M_{jj}^{JECup} [GeV] | CR-loose");                                          onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRlike_SSany_allSS_jesdown");    axisnames.push_back("M_{jj}^{JECdown} [GeV] | CR-loose");                                        onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRloose_butnotlike_SSany_allSS_wMjjL");axisnames.push_back("M_{jj} [GeV] | CR-preselection but not loose");                       onlySSSFOS.push_back(2);
    histonames.push_back("Mjj_CRloose_butnotlike_SSany_allSS_wMjjLDeta");axisnames.push_back("M_{jj} [GeV] | CR-preselection but not loose");                   onlySSSFOS.push_back(2);
    
    //histonames.push_back("Mll_CR_SSany_dropMjj");               axisnames.push_back("M_{ll} [GeV] | SS-CR");                                                  onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRlike_allSS");                   axisnames.push_back("M_{ll} [GeV] | SS-CR");                                                    onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRlike_allSS_SSpairnotlead2");    axisnames.push_back("M_{ll} [GeV] | SS-CR, SS not 2 lead leps");                                onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRlike_allSS_SSpairnotlead2v2");  axisnames.push_back("M_{ll} [GeV] | SS-CR, SS not 2 lead leps, M(ll) on lead leps");            onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRnoMjj_allSS");                  axisnames.push_back("M_{ll} [GeV] | SS-CR, no Mjj cut");                                        onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRnoMjj_allSS_SSpairnotlead2");   axisnames.push_back("M_{ll} [GeV] | SS-CR, no Mjj cut and SS not 2 lead leps");                 onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRnoMjj_allSS_SSpairnotlead2v2"); axisnames.push_back("M_{ll} [GeV] | SS-CR, no Mjj cut and SS not 2 lead leps, M(ll) on lead leps"); onlySSSFOS.push_back(2);
    
    histonames.push_back("Mjj_CRloose_allSS");                  axisnames.push_back("M_{jj} [GeV] | CR loose");                                                 onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRlike_allSS");                   axisnames.push_back("M_{ll} [GeV] | CR");                                                       onlySSSFOS.push_back(2);
    histonames.push_back("Mll_CRloose_allSS");                  axisnames.push_back("M_{ll} [GeV] | CR loose");                                                 onlySSSFOS.push_back(2);
    histonames.push_back("MET_CRlike_allSS");                   axisnames.push_back("E_{T}^{miss} [GeV] | CR");                                                 onlySSSFOS.push_back(2);
    histonames.push_back("MET_CRloose_allSS");                  axisnames.push_back("E_{T}^{miss} [GeV] | CR loose");                                           onlySSSFOS.push_back(2);
    histonames.push_back("Detajj_CRlike_allSS");                axisnames.push_back("Delta#eta_{jj} | CR");                                                     onlySSSFOS.push_back(2);
    histonames.push_back("Detajj_CRloose_allSS");               axisnames.push_back("Delta#eta_{jj} | CR loose");                                               onlySSSFOS.push_back(2);
    histonames.push_back("MTmax_CRlike_allSS");                 axisnames.push_back("M_{T}^{max} [GeV] | CR");                                                  onlySSSFOS.push_back(2);
    histonames.push_back("MTmax_CRloose_allSS");                axisnames.push_back("M_{T}^{max} [GeV] | CR loose");                                            onlySSSFOS.push_back(2);
    
    histonames.push_back("MET_CRlike_allSFOS");                 axisnames.push_back("E_{T}^{miss} [GeV] | CR");                                                 onlySSSFOS.push_back(2);
    histonames.push_back("MET_CRloose_allSFOS");                axisnames.push_back("E_{T}^{miss} [GeV] | CR loose");                                           onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_CRlike_allSFOS");               axisnames.push_back("M_{SFOS} [GeV] | CR");                                                     onlySSSFOS.push_back(2);
    histonames.push_back("MSFOS_CRloose_allSFOS");              axisnames.push_back("M_{SFOS} [GeV] | CR loose");                                               onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_CRlike_allSFOS");               axisnames.push_back("p_{T}(lll) [GeV] | CR");                                                   onlySSSFOS.push_back(2);
    histonames.push_back("pTlll_CRloose_allSFOS");              axisnames.push_back("p_{T}(lll) [GeV] | CR loose");                                             onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_CRlike_allSFOS");          axisnames.push_back("#Delta#phi(lll,E_{T}^{miss}) | CR");                                       onlySSSFOS.push_back(2);
    histonames.push_back("dPhiMETlll_CRloose_allSFOS");         axisnames.push_back("#Delta#phi(lll,E_{T}^{miss}) | CR loose");                                 onlySSSFOS.push_back(2);
    
    histonames.push_back("Mll_invertdPhiPt_allSFOS");           axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, #Delta$phi(lll,E_{T}^{miss}) and p_{T}(lll)");     onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhiOrPt_allSFOS");         axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, #Delta$phi(lll,E_{T}^{miss}) or p_{T}(lll)");      onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertdPhi_allSFOS");             axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, invert #Delta$phi(lll,E_{T}^{miss})");             onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertPt_allSFOS");               axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, invert p_{T}(lll)");                               onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertMET_allSFOS");              axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, invert E_{T}^{miss}");                             onlySSSFOS.push_back(2);
    histonames.push_back("Mll_invertMETdPhiPt_allSFOS");        axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, invert E_{T}^{miss}, #Delta$phi(lll,E_{T}^{miss}) and p_{T}(lll)");    onlySSSFOS.push_back(2);
    histonames.push_back("Mll_inverteitherMETdPhiPt_allSFOS");  axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, invert E_{T}^{miss}, #Delta$phi(lll,E_{T}^{miss}) or p_{T}(lll)");     onlySSSFOS.push_back(2);
    histonames.push_back("Mll_SRlike_allSFOS");                 axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, SR");                                              onlySSSFOS.push_back(2);
    histonames.push_back("Mll_ge2j_allSFOS");                   axisnames.push_back("M_{ll} [GeV] | 1+2SFOS, #geq2 jets");                                      onlySSSFOS.push_back(2);
    
    //histonames.push_back("el_mID_SRpresel");                    axisnames.push_back("lep_motherIdSS for e");                                                  onlySSSFOS.push_back(0);
    //histonames.push_back("mu_mID_SRpresel");                    axisnames.push_back("lep_motherIdSS for #mu");                                                onlySSSFOS.push_back(0);
     */
    TColor *lightblue  = new TColor(2001,91/255.,187/255.,241/255.);
    TColor *blue       = new TColor(2002,60/255.,144/255.,196/255.);
    TColor *orange     = new TColor(2003,230/255.,159/255.,0/255.);
    TColor *brown      = new TColor(2004,180/255.,117/255.,0/255.);
    TColor *yellow     = new TColor(2005,245/255.,236/255.,69/255.);
    TColor *darkyellow = new TColor(2006,215/255.,200/255.,0/255.);
    TColor *blueviolet = new TColor(2007,70/255.,109/255.,171/255.);
    TColor *violet     = new TColor(2008,70/255.,90/255.,134/255.);
    TColor *darkviolet = new TColor(2009,55/255.,65/255.,100/255.);
    TColor *lightgreen = new TColor(2010,120/255.,160/255.,0/255.);
    TColor *green      = new TColor(2011,0/255.,158/255.,115/255.);
    TColor *pink       = new TColor(2012,204/255.,121/255.,167/255.);


    vector<int> is3lSS;
    samples.push_back("WWW");         samplesleg.push_back("WWW");                          is3lSS.push_back(0);
    if(twosig) { samples.push_back("WHtoWWW");     samplesleg.push_back("WH#rightarrowWWW");is3lSS.push_back(0);}
    samples.push_back("others");      samplesleg.push_back("Other");                        is3lSS.push_back(0);
    
    samples.push_back("photonfakes");      samplesleg.push_back("#gamma #rightarrow l");    is3lSS.push_back(0);
    
    samples.push_back("chargeflips"); samplesleg.push_back("charge flip");                  is3lSS.push_back(0);
    //samples.push_back("doublefakes"); samplesleg.push_back("double fakes");                 is3lSS.push_back(0);
    samples.push_back("fakes");       samplesleg.push_back("jet fakes");                    is3lSS.push_back(0);
    //samples.push_back("FakePred");    samplesleg.push_back("fakes from data");              is3lSS.push_back(0);
    samples.push_back("3lLL");        samplesleg.push_back("lost lepton");                  is3lSS.push_back(2);
    samples.push_back("SSLL");        samplesleg.push_back("lost lepton");                  is3lSS.push_back(1);
    samples.push_back("true3L");      samplesleg.push_back("3l");                           is3lSS.push_back(2);
    samples.push_back("trueWWW");     samplesleg.push_back("WWW bg (ttW)");                 is3lSS.push_back(2);
    samples.push_back("trueSS");      samplesleg.push_back("SS bg (W^{#pm}W^{#pm},ttW)");   is3lSS.push_back(1);
    samples.push_back("bg");          samplesleg.push_back("background");                   is3lSS.push_back(-1);
    
    
    cols.push_back(632);//kRed
    if(twosig) cols.push_back(kBlue);//kRed
    cols.push_back(2012);
    
    cols.push_back(kGray);
    
    cols.push_back(2007);
    //cols.push_back(2006);
    cols.push_back(2005);
    cols.push_back(2011);
    cols.push_back(2003);
    cols.push_back(2003);
    cols.push_back(2001);
    cols.push_back(2001);
    cols.push_back(kBlack);
    /*
    samples.push_back("WWW");         samplesleg.push_back("WWW");                          onlySSSFOS.push_back(0);
    if(twosig) { samples.push_back("WHtoWWW");     samplesleg.push_back("WH#rightarrowWWW");             onlySSSFOS.push_back(0);}
    samples.push_back("others");      samplesleg.push_back("Other");                        onlySSSFOS.push_back(0);

    samples.push_back("photonfakes");      samplesleg.push_back("#gamma #rightarrow l");        onlySSSFOS.push_back(0);
    
    samples.push_back("chargeflips"); samplesleg.push_back("charge flip");                  onlySSSFOS.push_back(0);
    //samples.push_back("doublefakes"); samplesleg.push_back("double fakes");                 onlySSSFOS.push_back(0);
    samples.push_back("fakes");       samplesleg.push_back("jet fakes");                 onlySSSFOS.push_back(0);
    samples.push_back("3lLL");        samplesleg.push_back("lost lepton");                  onlySSSFOS.push_back(2);
    samples.push_back("SSLL");        samplesleg.push_back("lost lepton");                  onlySSSFOS.push_back(1);
    samples.push_back("true3L");      samplesleg.push_back("3l");                           onlySSSFOS.push_back(2);
    samples.push_back("trueWWW");     samplesleg.push_back("WWW bg (ttW)");                 onlySSSFOS.push_back(2);
    samples.push_back("trueSS");      samplesleg.push_back("SS bg (W^{#pm}W^{#pm},ttW)");   onlySSSFOS.push_back(1);
    samples.push_back("bg");          samplesleg.push_back("background");                   onlySSSFOS.push_back(-1);

  
    cols.push_back(632);//kRed
    if(twosig) cols.push_back(kBlue);//kRed
    cols.push_back(2012);

    cols.push_back(kGray);

    cols.push_back(2007);
    //cols.push_back(2006);
    cols.push_back(2005);
    cols.push_back(2011);
    cols.push_back(2003);
    cols.push_back(2003);
    cols.push_back(2001);
    cols.push_back(2001);
    cols.push_back(kBlack);
    */
    int nsig = 1;
    if(twosig) ++nsig;
    for(unsigned int n = 0; n<histonames.size();++n){
        for(int b = 0; b<3;++b){
            string SS = "";
            if(b!=0) continue;
            //if(b==0) SS = "ee";
            //if(b==1) SS = "em";
            //if(b==2) SS = "mm";
            string bgname = SS+histonames[n] + "_bg";
            string stackname = SS+histonames[n];
            //cout << stackname << endl;
            stack[stackname] = new THStack();
            stack[stackname]->SetName(stackname.c_str());
            string mapname = "";
            if(data) {
                mapname = SS+histonames[n] + "_Data";
                h[mapname ]=(TH1D*)f->Get(mapname.c_str());
                h[mapname]->SetLineColor(kBlack);
                h[mapname]->SetLineWidth(2);
                h[mapname]->SetMarkerStyle(20);
            }
            for(unsigned int s = 0; s<samples.size()-1; ++s){
                mapname = SS+histonames[n] + "_" + samples[s];
                cout << mapname << endl;
                h[mapname ]=(TH1D*)f->Get(mapname.c_str());
                h[mapname]->GetXaxis()->SetTitle(axisnames[n].c_str());
                h[mapname]->GetYaxis()->SetTitle("Events");
                h[mapname]->GetXaxis()->SetTitleSize(0.);
                //if(mapname.find("Mjj_")!=string::npos){
                //    if(!(mapname.find("MjjW")!=string::npos)){
                //        h[mapname]->SetBinContent(h[mapname]->GetNbinsX(),0.2*h[mapname]->GetBinContent(h[mapname]->GetNbinsX()));
                //        h[mapname]->SetBinError  (h[mapname]->GetNbinsX(),0.2*h[mapname]->GetBinError  (h[mapname]->GetNbinsX()));
                //    }
                //}
                //if(s==0&&(mapname.find("SSee")!=string::npos||mapname.find("1SFOS")!=string::npos||mapname.find("2SFOS")!=string::npos)) h[mapname]->Scale(5.);
                for(int b = 0; b<h[mapname]->GetNbinsX();++b){ if(h[mapname]->GetBinContent(b)<0){h[mapname]->SetBinContent(b,0); h[mapname]->SetBinError(b,0); } }
                if(s==nsig) h[bgname] = (TH1D*)h[mapname]->Clone(bgname.c_str());
                else if(s>nsig) {
                    h[bgname]->Add(h[mapname],1.);
                }
                if(s==0) {
                    h[mapname]->SetLineColor(cols[s]); h[mapname]->SetLineWidth(3);
                }
                else if(twosig&&s==1) {
                    h[mapname]->SetLineColor(cols[s]); h[mapname]->SetLineWidth(3); h[mapname]->SetLineStyle(7);
                }
                else { h[mapname]->SetFillColor(cols[s]); h[mapname]->SetLineColor(kBlack); }
                if(s==nsig){
                    h[bgname]->SetFillColor(cols[samples.size()-1]); h[bgname]->SetLineColor(cols[samples.size()-1]); h[bgname]->SetFillStyle(3544);
                }
                if(s>=nsig){
                    stack[stackname]->Add(h[mapname],"");
                }
                //cout << h[mapname]->Integral() << endl;
            }
            if(data&&mapname.find("Mjj_SRloose")!=string::npos){
                h[stackname + "_Data"]->SetBinContent(3,-999); h[stackname + "_Data"]->SetBinContent(4,-999);
            }
            if(data){
                for(int i = 1; i<=h[stackname + "_Data"        ]->GetNbinsX();++i){
                    //if(((h[stackname + "_"+samples[0] ]->GetBinContent(i)/h[stackname + "_bg"]->GetBinContent(i)>0.1) || (h[stackname + "_"+samples[0] ]->GetBinContent(i)/h[stackname + "_bg"]->GetBinContent(i)>0.5*h[stackname + "_bg"]->GetBinError(i)))&&h[stackname + "_"+samples[0] ]->GetBinContent(i)>0.5 ) h[stackname + "_Data"]->SetBinContent(i,-999);
                }
            }
            if(data) rat[stackname] = (TH1D*)h[stackname + "_Data"        ]->Clone(stackname.c_str());
            else     rat[stackname] = (TH1D*)h[stackname + "_"+samples[0] ]->Clone(stackname.c_str());
            rat[stackname]->Divide(h[stackname + "_bg"]);
            if(twosig&&!data){
                rat[stackname+"v2"] = (TH1D*)h[stackname + "_"+samples[1] ]->Clone(stackname.c_str());
                rat[stackname+"v2"]->Divide(h[stackname + "_bg"]);
            }
        }
    }
    
    TLatex *tLumi = new TLatex(0.95,0.955,"35.9 fb^{-1} (13 TeV)");
    tLumi->SetNDC();
    tLumi->SetTextAlign(31);
    tLumi->SetTextFont(42);
    tLumi->SetTextSize(0.042);
    tLumi->SetLineWidth(2);
    TLatex *tECM = new TLatex(0.95,0.955,"(13 TeV)");
    tECM->SetNDC();
    tECM->SetTextAlign(31);
    tECM->SetTextFont(42);
    tECM->SetTextSize(0.042);
    tECM->SetLineWidth(2);
    //tLumi->Draw();
    TLatex *tCMS = new TLatex(0.125,0.955,"CMS");
    tCMS->SetNDC();
    tCMS->SetTextAlign(11);
    tCMS->SetTextFont(61);
    tCMS->SetTextSize(0.0525);
    tCMS->SetLineWidth(2);
    //tCMS->Draw();
    TLatex *tSim = new TLatex(0.225,0.955,"Simulation");
    tSim->SetNDC();
    tSim->SetTextAlign(11);
    tSim->SetTextFont(52);
    tSim->SetTextSize(0.042);
    tSim->SetLineWidth(2);
    TLatex *tPrel = new TLatex(0.295,0.955,"Preliminary");
    tPrel->SetNDC();
    tPrel->SetTextAlign(11);
    tPrel->SetTextFont(52);
    tPrel->SetTextSize(0.042);
    tPrel->SetLineWidth(2);
    TLatex *tComment = new TLatex(0.185,0.90,"Last bin scaled x0.2");
    tComment->SetNDC();
    tComment->SetTextAlign(11);
    tComment->SetTextFont(42);
    tComment->SetTextSize(0.035);
    tComment->SetLineWidth(2);
    TLatex *tlx = new TLatex();
    tlx->SetTextFont(42);
    tlx->SetNDC();
    tlx->SetTextAlign(11);
    tlx->SetTextSize(0.042);
    tlx->SetTextColor(1);
    


    TLegend *leg1 = new TLegend(0.1667,0.63,0.5,0.9075,NULL,"brNDC");
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.033);
    leg1->SetLineColor(1);
    leg1->SetLineStyle(1);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1001);
    TLegend *leg1SS = new TLegend(0.2,0.67,0.5,0.89,NULL,"brNDC");
    leg1SS->SetBorderSize(0);
    leg1SS->SetTextSize(0.035);
    leg1SS->SetLineColor(1);
    leg1SS->SetLineStyle(1);
    leg1SS->SetLineWidth(2);
    leg1SS->SetFillColor(0);
    leg1SS->SetFillStyle(1001);
    TLegend *leg13l = new TLegend(0.2,0.67,0.5,0.89,NULL,"brNDC");
    leg13l->SetBorderSize(0);
    leg13l->SetTextSize(0.035);
    leg13l->SetLineColor(1);
    leg13l->SetLineStyle(1);
    leg13l->SetLineWidth(2);
    leg13l->SetFillColor(0);
    leg13l->SetFillStyle(1001);
    int legcount(0);
    int legcountSS(0);
    int legcount3l(0);
    for(unsigned int s = nsig; s<(samples.size()-1);++s){
        if(is3lSS[s]==0||is3lSS[s]==1) {
            leg1->AddEntry(h[histonames[0]+"_"+samples[s] ],samplesleg[s].c_str(), "f");
            legcount = s;
        }
    }
    for(unsigned int s = nsig; s<(samples.size()/2)+2;++s){
        if(is3lSS[s]==0||is3lSS[s]==1) {
            leg1SS->AddEntry(h[histonames[0]+"_"+samples[s] ],samplesleg[s].c_str(), "f");
            legcountSS = s;
        }
        if(is3lSS[s]==0||is3lSS[s]==2) {
            leg13l->AddEntry(h[histonames[0]+"_"+samples[s] ],samplesleg[s].c_str(), "f");
            legcount3l = s;
        }
    }
    TLegend *leg2 = new TLegend(0.55,0.63,0.85,0.9075,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.033);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1001);
    TLegend *leg2SS = new TLegend(0.55,0.67,0.85,0.89,NULL,"brNDC");
    leg2SS->SetBorderSize(0);
    leg2SS->SetTextSize(0.035);
    leg2SS->SetLineColor(1);
    leg2SS->SetLineStyle(1);
    leg2SS->SetLineWidth(2);
    leg2SS->SetFillColor(0);
    leg2SS->SetFillStyle(1001);
    TLegend *leg23l = new TLegend(0.55,0.67,0.85,0.89,NULL,"brNDC");
    leg23l->SetBorderSize(0);
    leg23l->SetTextSize(0.035);
    leg23l->SetLineColor(1);
    leg23l->SetLineStyle(1);
    leg23l->SetLineWidth(2);
    leg23l->SetFillColor(0);
    leg23l->SetFillStyle(1001);
    for(unsigned int s = nsig; s<(samples.size()-1);++s){
        if(is3lSS[s]==0||is3lSS[s]==2) {
            leg2->AddEntry(h[histonames[0]+"_"+samples[s] ],samplesleg[s].c_str(), "f");
            legcount = s;
        }
    }
    for(unsigned int s = legcountSS+1; s<samples.size()-1;++s){
        if(is3lSS[s]==0||is3lSS[s]==1) leg2SS->AddEntry(h[histonames[0]+"_"+samples[s] ],samplesleg[s].c_str(), "f");
    }
    for(unsigned int s = legcount3l+1; s<samples.size()-1;++s){
        if(is3lSS[s]==0||is3lSS[s]==2) leg23l->AddEntry(h[histonames[0]+"_"+samples[s] ],samplesleg[s].c_str(), "f");
    }
    leg1->AddEntry(h[histonames[0]+"_"+samples[0] ],samplesleg[0].c_str(), "l");
    leg2->AddEntry(h[histonames[0]+"_"+samples[0] ],samplesleg[0].c_str(), "l");
    leg2SS->AddEntry(h[histonames[0]+"_"+samples[0] ],samplesleg[0].c_str(), "l");
    leg23l->AddEntry(h[histonames[0]+"_"+samples[0] ],samplesleg[0].c_str(), "l");
    if(twosig){
        leg1->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
        leg2->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
        leg2SS->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
        leg23l->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
    }
    if(data){
        leg1->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
        leg2->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
        leg2SS->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
        leg23l->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
    }

    for(unsigned int n = 0; n<histonames.size();++n){
        for(int b =0; b<3; ++b){
            string SS="";
            if(b!=0) continue;
            //if(b==0) SS = "ee";
            //if(b==1) SS = "em";
            //if(b==2) SS = "mm";
            //cout << b << " " << SS << endl;
            TCanvas *c1 = new TCanvas("c1", "",334,192,600,600);
            c1->SetFillColor(0);
            c1->SetBorderMode(0);
            c1->SetBorderSize(2);
            //if(logy) c1->SetLogy();    // Log y
            c1->SetTickx(1);
            c1->SetTicky(1);
            c1->SetLeftMargin(0.18);
            c1->SetRightMargin(0.05);
            c1->SetTopMargin(0.07);
            c1->SetBottomMargin(0.15);
            c1->SetFrameFillStyle(0);
            c1->SetFrameBorderMode(0);
            c1->SetFrameFillStyle(0);
            c1->SetFrameBorderMode(0);
            TPad *plotpad = new TPad("plotpad", "Pad containing the overlay plot",0,0.165,1,1);//0,0.18,1,1);
            plotpad->Draw();
            plotpad->cd();
            plotpad->Range(-85.71429,-3.864499,628.5714,6.791402);//(133.1169,-3.101927,782.4675,0.7583922);
            plotpad->SetFillColor(0);
            plotpad->SetBorderMode(0);
            plotpad->SetBorderSize(2);
            if(logy) plotpad->SetLogy();
            plotpad->SetTickx(1);
            plotpad->SetTicky(1);
            plotpad->SetLeftMargin(0.12);
            plotpad->SetRightMargin(0.04);
            plotpad->SetTopMargin(0.05);
            // plotpad->SetBottomMargin(0.15);
            plotpad->SetFrameFillStyle(0);
            plotpad->SetFrameBorderMode(0);
            plotpad->SetFrameFillStyle(0);
            plotpad->SetFrameBorderMode(0);
        
            plotpad->cd();
    
            string stackname = SS+histonames[n];
            string bgname = stackname + "_bg";
            //cout << stackname << endl;
            double max;
            if(data) max = TMath::Max((h[bgname]->GetBinContent(h[bgname]->GetMaximumBin() )+0.5*h[bgname]->GetBinError(h[bgname]->GetMaximumBin() ) ),(h[stackname+"_Data"]->GetBinContent(h[stackname+"_Data"]->GetMaximumBin() )+0.5*h[stackname+"_Data"]->GetBinError(h[stackname+"_Data"]->GetMaximumBin() ) ) )*1.45;
            else max = (h[bgname]->GetBinContent(h[bgname]->GetMaximumBin() )+0.5*h[bgname]->GetBinError(h[bgname]->GetMaximumBin()))*1.45;
            stack[stackname]->SetMaximum(max);
            stack[stackname]->Draw("hist");
            h[bgname]->Draw("sameE2");
            if(data) h[stackname+"_Data"]->Draw("sameE0X0");
            string signame = stackname + "_"+ samples[0];
            h[signame]->Draw("histsame");
            signame = stackname + "_"+ samples[1];
            if(twosig) h[signame]->Draw("histsame");
            string outname = stackname + ".pdf";
            if(logy) outname = "plots/processsplit/log/" + outname;
            else     outname = "plots/processsplit/" + outname;
            //if(logy) outname = "plots/processsplit/excessDataFakes/log/" + outname;
            //else     outname = "plots/processsplit/excessDataFakes/" + outname;

            if(onlySSSFOS[n]==0){ leg1  ->Draw(); leg2  ->Draw(); }
            if(onlySSSFOS[n]==1){ leg1SS->Draw(); leg2SS->Draw(); }
            if(onlySSSFOS[n]==2){ leg13l->Draw(); leg23l->Draw(); }
            tLumi->Draw();
            //if(stackname.find("Mjj_")!=string::npos){
            //    if(!(stackname.find("MjjW")!=string::npos)){
            //        tComment->Draw();
            //    }
            //}
            tCMS->Draw();
            tSim->Draw();
            //if(b==0) tlx->DrawLatex(0.6,0.95,"e^{#pm}e^{#pm}");
            //if(b==1) tlx->DrawLatex(0.6,0.95,"e^{#pm}#mu^{#pm}");
            //if(b==2) tlx->DrawLatex(0.6,0.95,"#mu^{#pm}#mu^{#pm}");
            if(stackname.find("SSee" )  !=string::npos) tlx->DrawLatex(0.45,0.955,"e^{#pm}e^{#pm}");
            if(stackname.find("SSem" )  !=string::npos) tlx->DrawLatex(0.45,0.955,"e^{#pm}#mu^{#pm}");
            if(stackname.find("SSmm" )  !=string::npos) tlx->DrawLatex(0.45,0.955,"#mu^{#pm}#mu^{#pm}");
            if(stackname.find("allSS")  !=string::npos) tlx->DrawLatex(0.45,0.955,"SS region");
            if(stackname.find("0SFOS")  !=string::npos) tlx->DrawLatex(0.45,0.955,"0 SFOS pairs");
            if(stackname.find("1SFOS")  !=string::npos) tlx->DrawLatex(0.45,0.955,"1 SFOS pair");
            if(stackname.find("2SFOS")  !=string::npos) tlx->DrawLatex(0.45,0.955,"2 SFOS pairs");
            if(stackname.find("allSFOS")!=string::npos) tlx->DrawLatex(0.45,0.955,"3l region");
            //if(stackname.find("SSee")!=string::npos||stackname.find("1SFOS")!=string::npos||stackname.find("2SFOS")!=string::npos) tlx->DrawLatex(0.65,0.9,"WWW scaled by x5");

            c1->cd();
            TPad *ratiopad = new TPad("ratiopad", "Pad containing the ratio",0,0,1,0.21); //0,0,1,0.26);
            ratiopad->Draw();
            ratiopad->cd();
            //ratiopad->Range(-85.71429,-0.4,628.5714,2.266667);  //(133.1169,0.06923079,782.4675,1.607692);
            ratiopad->SetFillColor(0);
            ratiopad->SetBorderMode(0);
            ratiopad->SetBorderSize(2);
            ratiopad->SetTickx(1);
            ratiopad->SetTicky(1);
            ratiopad->SetLeftMargin(0.12);
            ratiopad->SetRightMargin(0.04);
            //  ratiopad->SetTopMargin(0.04);
            ratiopad->SetBottomMargin(0.3);
            ratiopad->SetFrameFillStyle(0);
            ratiopad->SetFrameBorderMode(0);
            ratiopad->SetFrameFillStyle(0);
            ratiopad->SetFrameBorderMode(0);
            if(data){
                rat[stackname]->SetMinimum(0.);
                rat[stackname]->SetMaximum(2.);
            }
            else {
                rat[stackname]->SetMinimum(0.);
                rat[stackname]->SetMaximum(0.2);
            }
            rat[stackname]->GetXaxis()->SetTitle(axisnames[n].c_str());
            rat[stackname]->GetXaxis()->SetTitleSize(0.16);
            rat[stackname]->GetXaxis()->SetTitleOffset(0.76);
            rat[stackname]->GetXaxis()->SetLabelSize(0.0);
            rat[stackname]->GetYaxis()->SetNdivisions(504);
            rat[stackname]->GetYaxis()->SetTitle("signal / bgd");
            if(data) rat[stackname]->GetYaxis()->SetTitle("data / bgd");
            rat[stackname]->GetYaxis()->SetTitleSize(0.14);
            rat[stackname]->GetYaxis()->SetTitleOffset(0.28);
            rat[stackname]->GetYaxis()->SetLabelSize(0.14);
            rat[stackname]->Draw();
            if(twosig&&!data){
                rat[stackname+"v2"]->SetMinimum(0.);
                rat[stackname+"v2"]->SetMaximum(0.2);
                rat[stackname+"v2"]->GetXaxis()->SetTitle(axisnames[n].c_str());
                rat[stackname+"v2"]->GetXaxis()->SetTitleSize(0.16);
                rat[stackname+"v2"]->GetXaxis()->SetTitleOffset(0.76);
                rat[stackname+"v2"]->GetXaxis()->SetLabelSize(0.0);
                rat[stackname+"v2"]->GetYaxis()->SetNdivisions(504);
                rat[stackname+"v2"]->GetYaxis()->SetTitle("signal / bgd");
                if(data) rat[stackname+"v2"]->GetYaxis()->SetTitle("data / bgd");
                rat[stackname+"v2"]->GetYaxis()->SetTitleSize(0.14);
                rat[stackname+"v2"]->GetYaxis()->SetTitleOffset(0.28);
                rat[stackname+"v2"]->GetYaxis()->SetLabelSize(0.14);
                rat[stackname+"v2"]->Draw("same");
            }
            if(data){
                TLine *rline = new TLine(rat[stackname]->GetXaxis()->GetBinLowEdge(1),1.,rat[stackname]->GetXaxis()->GetBinLowEdge(rat[stackname]->GetNbinsX()+1),1.);
                rline->SetLineWidth(2);
                rline->SetLineStyle(7);
                rline->Draw();
            }
            //    for(int i =1; i<=rat[stackname]->GetNbinsX();++i) cout << stackname << " " << i << " " << rat[stackname]->GetBinContent(i) << " +/- " << rat[stackname]->GetBinError(i) << "   " << h[signame]->GetBinContent(i) << " +/- " << h[signame]->GetBinError(i) << endl;
            c1->cd();
            c1->SaveAs(outname.c_str());
            //c1->Clear();
            c1->cd();
        }
    }
}
