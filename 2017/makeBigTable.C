float Zbi_(float sig, float bg, float bgrelunc=0.3){
    
    double bgunc = bgrelunc*bg;
    double tau = bg/pow(bgunc,2);//bgunc is absolute
    double n_on = sig+bg;//total yield in SR = sig + bg
    double n_off = tau*bg;
    double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);
    double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
    return Z_Bi;
}

void makeBigTable(){
    
    bool raw = false;
    bool SoverB = false;
    bool twosig = true;
    map<string, TH1D*> h;
    vector<string> histonames;
    vector<string> caption;
    vector<string> samples;
    vector<string> samples2;
    TFile *f = TFile::Open("Check3lCRv2.root");
    histonames.push_back("YieldsSR");                                      caption.push_back("SR");
    histonames.push_back("YieldsCR_using3lorMll");                         caption.push_back("CR: in SS add 3rd $\\ell$. Require SFOS pair within Z mass window: $\\pm$10 GeV (SS); invert SFOS veto (3l)");
    histonames.push_back("YieldsCR_dropMjj");                              caption.push_back("CR: in SS drop Mjj-W mass cut");
    histonames.push_back("YieldsCR_SSany_dropMjj");                              caption.push_back("CR: in SS drop Mjj-W mass cut, any leps SS pair");
    histonames.push_back("YieldsCR_SSany_dropMjj_test");                              caption.push_back("CRtest: in SS drop Mjj-W mass cut, any leps SS pair");
    histonames.push_back("YieldsCR_SSany_dropMjj_jesdown");                              caption.push_back("CR: in SS drop Mjj-W mass cut, any leps SS pair, JESDOWN");
    histonames.push_back("YieldsCR_SSany_dropMjj_jesup");                              caption.push_back("CR: in SS drop Mjj-W mass cut, any leps SS pair, JESUP");
    //histonames.push_back("YieldsSR_invertdPhiPtfor3Mjjfor2l");             caption.push_back("SR: SS inverted Mjj-W mass cut; 3l inverted $\\Delta\\phi(\\ell\\ell\\ell,E_\\mathrm{T}^\\mathrm{miss})$, $p_\\mathrm{T}(\\ell\\ell\\ell)$ cut");
    //histonames.push_back("YieldsCR_invertdPhiPtfor3Mjjfor2l");             caption.push_back("CR: SS inverted Mjj-W mass cut; 3l inverted $\\Delta\\phi(\\ell\\ell\\ell,E_\\mathrm{T}^\\mathrm{miss})$, $p_\\mathrm{T}(\\ell\\ell\\ell)$ cut");
    //histonames.push_back("YieldsSR_lowMET");                               caption.push_back("SR but invert MET cut");
    //histonames.push_back("YieldsCR_lowMET");                               caption.push_back("CR but invert MET cut");
    //histonames.push_back("YieldsSR_lowMETinvertdPhiPtMjj");                caption.push_back("SR but invert MET cut: SS inverted Mjj-W mass cut; 3l inverted $\\Delta\\phi(\\ell\\ell\\ell,E_\\mathrm{T}^\\mathrm{miss})$, $p_\\mathrm{T}(\\ell\\ell\\ell)$ cut");
    //histonames.push_back("YieldsCR_lowMETinvertdPhiPtMjj");                caption.push_back("CR but invert MET cut: SS inverted Mjj-W mass cut; 3l inverted $\\Delta\\phi(\\ell\\ell\\ell,E_\\mathrm{T}^\\mathrm{miss})$, $p_\\mathrm{T}(\\ell\\ell\\ell)$ cut");

    //histonames.push_back("YieldsCR_dropcutsbutNJisotr");                   caption.push_back("CR: add 3rd $\\ell$ forming a Z (SS). Keep only NJ/NB/isotrack");
    //histonames.push_back("YieldsCR_dropcutsbutNJvetolep");                 caption.push_back("CR: add 3rd $\\ell$ forming a Z (SS); invert SFOS veto (3l). Keep only NJ/NB/vetolep");
    //histonames.push_back("YieldsCR_dropcutsbutNJ");                        caption.push_back("CR: add 3rd $\\ell$ forming a Z (SS); invert SFOS veto (3l). Keep only NJ/NB/$\\geq$3lep");
    //histonames.push_back("YieldsCRnoMll_dropcutsbutNJisotr");              caption.push_back("CR: add 3rd $\\ell$ (SS). Keep only NJ/NB/isotrack");
    //histonames.push_back("YieldsCRnoMll_dropcutsbutNJvetolep");            caption.push_back("CR: add 3rd $\\ell$ (SS). Keep only NJ/NB/vetolep");
    //histonames.push_back("YieldsCRnoMll_dropcutsbutNJ");                   caption.push_back("CR: add 3rd $\\ell$ (SS). Keep only NJ/NB/$\\geq$3lep");

    histonames.push_back("YieldsSR_jesup");                   caption.push_back("SR - jesUP");
    histonames.push_back("YieldsSR_jesdown");                 caption.push_back("SR - jesDOWN");
    histonames.push_back("YieldsCR_dropMjj_jesup");           caption.push_back("CR w/o Mjj-W mass cut - jesUP");
    histonames.push_back("YieldsCR_dropMjj_jesdown");         caption.push_back("CR w/o Mjj-W mass cut - jesDOWN");
    histonames.push_back("YieldsCR_using3lorMll_jesup");      caption.push_back("CR - jesUP");
    histonames.push_back("YieldsCR_using3lorMll_jesdown");    caption.push_back("CR - jesDOWN");
    histonames.push_back("YieldsCR_dropcutsbutNJ_jesup");     caption.push_back("CR loose - jesUP");
    histonames.push_back("YieldsCR_dropcutsbutNJ_jesdown");   caption.push_back("CR loose - jesDOWN");
    histonames.push_back("YieldsSR_preselection");            caption.push_back("SR preselection - veto SR data");
    //histonames.push_back("YieldsSR_rawweight");                                caption.push_back("SR");
    //histonames.push_back("YieldsSR_raw");                                      caption.push_back("SR");

    //   TFile *f = TFile::Open("CheckWGrelatedstuff.root");
    /*
    histonames.push_back("YieldsSRmumue_0mhits_invertMETDPhiOrPt");        caption.push_back("SRmumue: 0 missing hits");
    histonames.push_back("YieldsCRmumue_0mhits_invertMETDPhiOrPt");        caption.push_back("CRmumue: 0 missing hits");
    histonames.push_back("YieldsSRmumue_1mhits_invertMETDPhiOrPt");        caption.push_back("SRmumue: 1 missing hits");
    histonames.push_back("YieldsCRmumue_1mhits_invertMETDPhiOrPt");        caption.push_back("CRmumue: 1 missing hits");
    histonames.push_back("YieldsSRmumue_2mhits_invertMETDPhiOrPt");        caption.push_back("SRmumue: 2 missing hits");
    histonames.push_back("YieldsCRmumue_2mhits_invertMETDPhiOrPt");        caption.push_back("CRmumue: 2 missing hits");
     */
    /*
    histonames.push_back("YieldsSR_lowMSFOS_invertMETDPhiOrPt");                             caption.push_back("SRinverted, added low MSFOS$>$20 GeV cut");
    histonames.push_back("YieldsSR_lowMSFOS_Mlll_invertMETDPhiOrPt");                        caption.push_back("SRinverted, added low MSFOS$>$20 GeV, Mlll cuts");
    histonames.push_back("YieldsSR_dRllmin_invertMETDPhiOrPt");                              caption.push_back("SRinverted, add minDeltaR(l,l)>0.2 cut");
    histonames.push_back("YieldsSR_invertMETDPhiOrPt");                                      caption.push_back("SRinverted");
    histonames.push_back("YieldsSR_lowMSFOS");                             caption.push_back("SR, added low MSFOS$>$20 GeV cut");
    histonames.push_back("YieldsSR_lowMSFOS_Mlll");                        caption.push_back("SR, added low MSFOS$>$20 GeV, Mlll cuts");
    histonames.push_back("YieldsSR_dRllmin");                              caption.push_back("SR, add minDeltaR(l,l)>0.2 cut");
    histonames.push_back("YieldsSR");                                      caption.push_back("SR");
    //histonames.push_back("YieldsCR");                                      caption.push_back("CR");
     */
    /*
    histonames.push_back("YieldsSR_tightcharge");                          caption.push_back("SR");
    histonames.push_back("YieldsCR_tightcharge");                          caption.push_back("CR");
    histonames.push_back("YieldsSR_nisotrack");                            caption.push_back("SR");
    histonames.push_back("YieldsCR_nisotrack");                            caption.push_back("CR");
    histonames.push_back("YieldsSR_nisotrack_tightcharge");                caption.push_back("SR");
    histonames.push_back("YieldsCR_nisotrack_tightcharge");                caption.push_back("CR");
    */
    /*
    TFile *f = TFile::Open("QuickAndDirtyChecks.root");
    histonames.push_back("SR");                            caption.push_back("SR");
    histonames.push_back("SR_Mjjsideband");                caption.push_back("SR - M_{jj} sideband");
    histonames.push_back("SR_fakesource");                 caption.push_back("SR fake source");
    histonames.push_back("SR_Mjjsideband_fakesource");     caption.push_back("SR fake source - M_{jj} sideband");
    histonames.push_back("AR_fakesource");                 caption.push_back("AR fake source");
    histonames.push_back("AR_Mjjsideband_fakesource");     caption.push_back("AR fake source - M_{jj} sideband");
    histonames.push_back("SR_presel_fakesource");          caption.push_back("SR preselection - fake source");
    histonames.push_back("AR_presel_fakesource");          caption.push_back("AR preselection - fake source");
    */
    /*
    TFile *f = TFile::Open("TryAdditionalValidationRegion.root");
    histonames.push_back("SR_ZplusLep_METle10");                             caption.push_back("SR_ZplusLep_METle10");
    histonames.push_back("SR_ZplusLep_MTle25OrMTge150");                     caption.push_back("SR_ZplusLep_MTle25OrMTge150");
    histonames.push_back("SR_ZplusLep_MTle25OrMTge150_METle10");             caption.push_back("SR_ZplusLep_MTle25OrMTge150_METle10");
    histonames.push_back("AR_ZplusLep_METle10");                             caption.push_back("AR_ZplusLep_METle10");
    histonames.push_back("AR_ZplusLep_MTle25OrMTge150");                     caption.push_back("AR_ZplusLep_MTle25OrMTge150");
    histonames.push_back("AR_ZplusLep_MTle25OrMTge150_METle10");             caption.push_back("AR_ZplusLep_MTle25OrMTge150");
    histonames.push_back("SR_ZplusLep_Mjjsideband_METle10");                 caption.push_back("SR_ZplusLep_Mjjsideband_METle10");
    histonames.push_back("SR_ZplusLep_Mjjsideband_MTle25OrMTge150");         caption.push_back("SR_ZplusLep_Mjjsideband_MTle25OrMTge150");
    histonames.push_back("SR_ZplusLep_Mjjsideband_MTle25OrMTge150_METle10"); caption.push_back("SR_ZplusLep_Mjjsideband_MTle25OrMTge150_METle10");
    histonames.push_back("AR_ZplusLep_Mjjsideband_METle10");                 caption.push_back("AR_ZplusLep_Mjjsideband_METle10");
    histonames.push_back("AR_ZplusLep_Mjjsideband_MTle25OrMTge150");         caption.push_back("AR_ZplusLep_Mjjsideband_MTle25OrMTge150");
    histonames.push_back("AR_ZplusLep_Mjjsideband_MTle25OrMTge150_METle10"); caption.push_back("AR_ZplusLep_Mjjsideband_MTle25OrMTge150_METle10");
    histonames.push_back("SR_2b_presel");                                    caption.push_back("SR_2b_presel");
    histonames.push_back("SR_2b");                                           caption.push_back("SR_2b");
    histonames.push_back("AR_2b_presel");                                    caption.push_back("AR_2b_presel");
    histonames.push_back("AR_2b");                                           caption.push_back("AR_2b");
    histonames.push_back("SR_2b_Mjjsideband_presel");                        caption.push_back("SR_2b_Mjjsideband_presel");
    histonames.push_back("SR_2b_Mjjsideband");                               caption.push_back("SR_2b_Mjjsideband");
    histonames.push_back("AR_2b_Mjjsideband_presel");                        caption.push_back("AR_2b_Mjjsideband_presel");
    histonames.push_back("AR_2b_Mjjsideband");                               caption.push_back("AR_2b_Mjjsideband");
    */
    samples.push_back("WWW");
    if(twosig) samples.push_back("WHtoWWW");
    samples.push_back("VVV");
    samples.push_back("WZ");
    samples.push_back("WW");
    samples.push_back("ZZ");
    samples.push_back("Wjets");
    samples.push_back("Zjets");
    samples.push_back("ttV");
    samples.push_back("tt2l");
    samples.push_back("tt1l");
    samples.push_back("singleTop");
    samples.push_back("WG");
    samples.push_back("ZG");
    samples.push_back("Other");
    samples.push_back("bg");
    /*
    samples2.push_back("Data");
    samples2.push_back("trueSS");
    samples2.push_back("trueWWW");
    samples2.push_back("true3L");
    samples2.push_back("SSLL");
    samples2.push_back("3lLL");
    samples2.push_back("chargeflips");
    samples2.push_back("fakes");
    samples2.push_back("doublefakes");
    samples2.push_back("photonfakes");
    samples2.push_back("photondoublefakes");
    samples2.push_back("photontriplefakes");
    samples2.push_back("fakesphotonfakes");
    samples2.push_back("otherphotonfakes");
    samples2.push_back("others");
    */
    samples2.push_back("Data");
    samples2.push_back("trueSS");
    samples2.push_back("trueWWW");
    samples2.push_back("true3L");
    samples2.push_back("SSLL");
    samples2.push_back("3lLL");
    samples2.push_back("chargeflips");
    samples2.push_back("fakes");
    samples2.push_back("photonfakes");
    samples2.push_back("others");

    int nsig = 1;
    if(twosig) ++nsig;

    for(unsigned int n = 0; n<histonames.size();++n){
        for(unsigned int s = 0; s<samples.size()-1; ++s){
            string mapname = histonames[n] + "_" + samples[s];
            string bgname = histonames[n] + "_bg";
            if(raw){
                mapname = "Raw" + mapname;
                bgname = "Raw" + bgname;
            }
            cout << mapname << endl;
            h[mapname ]=(TH1D*)f->Get(mapname.c_str());
            //if(samples[s]=="WWW") h[mapname ]->Scale(0.2086);
            for(int b = 0; b<h[mapname]->GetNbinsX();++b){ if(h[mapname]->GetBinContent(b)<0){h[mapname]->SetBinContent(b,0); h[mapname]->SetBinError(b,0); } }
            if(s==2) h[bgname] = (TH1D*)h[mapname]->Clone(bgname.c_str());
            else if(s>2) h[bgname]->Add(h[mapname],1.);
        }
        for(unsigned int s = 0; s<samples2.size(); ++s){
            string mapname = histonames[n] + "_" + samples2[s];
            h[mapname ]=(TH1D*)f->Get(mapname.c_str());
            for(int b = 0; b<h[mapname]->GetNbinsX();++b){ if(h[mapname]->GetBinContent(b)<0){h[mapname]->SetBinContent(b,0); h[mapname]->SetBinError(b,0); } }
        }
    }
    string ss;
    cout << "\\begin{table}" << endl
    << "\\centering" << endl
    << "\\tiny" << endl
    << "\\begin{tabular}{|r|ccc|ccc|}" << endl;
    for(unsigned int n = 0; n<histonames.size(); ++n){
    cout << "\\hline" << endl
    << "\\multicolumn{7}{|c|}{"<<caption[n] <<"} \\\\" << endl
    << " & \\multicolumn{3}{|c|}{SS} & \\multicolumn{3}{|c|}{3$\\ell$} \\\\" << endl
    << " & $ee$ & $e\\mu$ & $\\mu\\mu$ & 0SFOS & 1SFOS & 2SFOS \\\\" << endl
    << "\\hline" << endl;
    for(unsigned int s = nsig; s<samples.size(); ++s){
        ss = histonames[n] + "_" + samples[s];
        if(raw) ss = "Raw"+ss;
        cout << samples[s];
        for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b) << "$\\pm$" << h[ss]->GetBinError(b);
        cout << " \\\\" << endl;
    }
    cout << "\\hline" << endl;
    ss = histonames[n] + "_WWW";
    if(raw) ss = "Raw"+ss;
    cout << "WWW";
    for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b) << "$\\pm$" << h[ss]->GetBinError(b);
        cout << " \\\\" << endl;
        if(twosig){
        ss = histonames[n] + "_WHtoWWW";
       if(raw) ss = "Raw"+ss;
        cout << "WH$\\to$WWW";
        for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b) << "$\\pm$" << h[ss]->GetBinError(b);
            cout << " \\\\" << endl;
        }
    cout << "\\hline" << endl;
    if(SoverB){
    string ssb =histonames[n] + "_bg";
    if(raw) ssb = "Raw"+ssb;
    cout << "S/B ";
    for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b)/h[ssb]->GetBinContent(b) << "$\\pm$" << sqrt(pow(h[ss]->GetBinError(b)/h[ssb]->GetBinContent(b),2)+pow(h[ss]->GetBinContent(b)*h[ssb]->GetBinError(b)/pow(h[ssb]->GetBinContent(b),2),2));
    cout << "\\\\" << endl;
    cout << "S/$\\sqrt{\\mathrm{B}}$ ";
    for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b)/sqrt(h[ssb]->GetBinContent(b)) << "$\\pm$" << sqrt(pow(h[ss]->GetBinError(b)/sqrt(h[ssb]->GetBinContent(b)),2)+pow(h[ss]->GetBinContent(b)*h[ssb]->GetBinError(b)/pow(h[ssb]->GetBinContent(b),1.5),2));
    cout << "\\\\" << endl;
    cout << "$Z_\\mathrm{bi}$ ";
    for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << Zbi_(h[ss]->GetBinContent(b),h[ssb]->GetBinContent(b),sqrt(pow(h[ssb]->GetBinError(b)/h[ssb]->GetBinContent(b),2)+pow(0.3,2)) );
    cout << "\\\\" << endl
    << "\\hline" << endl;
    }
    ss = histonames[n] + "_Data";
    if(raw) ss = "Raw"+ss;
    cout << "Data";
    for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b) << "$\\pm$" << h[ss]->GetBinError(b);
    cout << " \\\\" << endl
    << "\\hline" << endl;
    for(unsigned int s = 1; s<samples2.size(); ++s){
        ss = histonames[n] + "_" + samples2[s];
        if(raw) ss = "Raw"+ss;
        cout << samples2[s];
        for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b) << "$\\pm$" << h[ss]->GetBinError(b);
        cout << " \\\\" << endl;
    }
        cout << "\\hline" << endl;
    }
    cout << "\\end{tabular}" << endl
    << "\\end{table}" << endl;


}
