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
    bool twosig = false;
    boll addsig = true;
    map<string, TH1D*> h;
    vector<string> histonames;
    vector<string> caption;
    vector<string> samples;
    vector<string> samples2;
    TFile *f = TFile::Open("rootfiles/newLooper.root");
    histonames.push_back("SignalRegion");                  caption.push_back("signal region");
    histonames.push_back("ApplicationRegion");             caption.push_back("application region");
    histonames.push_back("WZControlRegion");               caption.push_back("lost lepton / 3l-with-Z control region");
    histonames.push_back("SignalRegionPresel");            caption.push_back("signal region - preselection");
    histonames.push_back("ApplicationRegionPresel");       caption.push_back("application region - preselection");
    histonames.push_back("WZControlRegionPresel");         caption.push_back("lost lepton / 3l-with-Z control region - preselection");


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
            if(s==0&&addsig){
                string addname = histonames[n] + "_WHtoWWW";
                if(raw) addname = "Raw" + mapname;
                h[addname ]=(TH1D*)f->Get(addname.c_str());
                h[mapname ]->Add(h[addname],1.);
            }
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
