void makeQuickTable(){
    
    bool raw = false;
    bool SoverB = true;
    
    TFile *f = TFile::Open("StupidSync.root");
    map<string, TH1D*> h;
    vector<string> histonames;
    vector<string> caption;
    vector<string> samples;
    histonames.push_back("YieldsSync");                     caption.push_back("Sync Region");
    histonames.push_back("YieldsSyncTightIso");             caption.push_back("Sync Region, tight iso (e-POG style)");
    histonames.push_back("YieldsSyncTighttest");            caption.push_back("Sync Region, SUSY-SS iso");
    histonames.push_back("YieldsSync3025_252020");          caption.push_back("Sync Region, SS-pT=30,25 GeV, 3l-pT=25,20,20 GeV");
    histonames.push_back("YieldsSync3020_252015");          caption.push_back("Sync Region, SS-pT=30,20 GeV, 3l-pT=25,20,15 GeV");
    histonames.push_back("YieldsSync3015_252010");          caption.push_back("Sync Region, SS-pT=30,15 GeV, 3l-pT=25,20,10 GeV");
    histonames.push_back("YieldsSync3030_251515");          caption.push_back("Sync Region, SS-pT=30,30 GeV, 3l-pT=25,15,15 GeV");
    histonames.push_back("YieldsSync3030_251510");          caption.push_back("Sync Region, SS-pT=30,30 GeV, 3l-pT=25,15,10 GeV");
    histonames.push_back("YieldsSync2525_202015");          caption.push_back("Sync Region, SS-pT=25,25 GeV, 3l-pT=20,20,15 GeV");
    histonames.push_back("YieldsSync2520_202010");          caption.push_back("Sync Region, SS-pT=25,20 GeV, 3l-pT=20,20,10 GeV");
    histonames.push_back("YieldsSync2515_201515");          caption.push_back("Sync Region, SS-pT=25,15 GeV, 3l-pT=20,15,15 GeV");
    histonames.push_back("YieldsSync2525_201510");          caption.push_back("Sync Region, SS-pT=25,25 GeV, 3l-pT=20,15,10 GeV");
    histonames.push_back("YieldsSync2525_201010");          caption.push_back("Sync Region, SS-pT=25,25 GeV, 3l-pT=20,10,10 GeV");

    samples.push_back("WWW");
    samples.push_back("WZ");
    samples.push_back("WW");
    samples.push_back("ZZ");
    samples.push_back("Wjets");
    samples.push_back("Zjets");
    samples.push_back("ttV");
    samples.push_back("TT2l");
    samples.push_back("TT1l");
    samples.push_back("singleT");
    samples.push_back("bg");
    
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
            if(s==1) h[bgname] = (TH1D*)h[mapname]->Clone(bgname.c_str());
            else if(s>1) h[bgname]->Add(h[mapname],1.);
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
    for(unsigned int s = 1; s<samples.size(); ++s){
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
    cout << " \\\\" << endl
    << "\\hline" << endl;
    if(SoverB){
    string ssb =histonames[n] + "_bg";
    if(raw) ssb = "Raw"+ssb;
    cout << "S/B ";
    for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b)/h[ssb]->GetBinContent(b) << "$\\pm$" << sqrt(pow(h[ss]->GetBinError(b)/h[ssb]->GetBinContent(b),2)+pow(h[ss]->GetBinContent(b)*h[ssb]->GetBinError(b)/pow(h[ssb]->GetBinContent(b),2),2));
    cout << "\\\\" << endl
     << "\\hline" << endl;
    }
        
    }
    cout << "\\end{tabular}" << endl
    << "\\end{table}" << endl;
    /*
    cout << endl << endl << "\\begin{table}" << endl
    << "\\centering" << endl
    << "\\begin{tabular}{|r|c|c|}" << endl
    << "\\hline" << endl
    << "\\multicolumn{3}{|c|}{Bobak's sync} \\\\" << endl
    << " & $\\mu\\mu$ & 0SFOS \\\\" << endl
    << "\\hline" << endl;
    for(unsigned int s = 1; s<samples.size(); ++s){
        ss = "Bobak_" + samples[s];
        if(raw) ss = "Raw"+ss;
       cout << samples[s];
        for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b) << "$\\pm$" << h[ss]->GetBinError(b);
        cout << " \\\\" << endl;
    }
    cout << "\\hline" << endl;
    ss = "Bobak_WWW";
    if(raw) ss = "Raw"+ss;
    cout << "WWW";
    for(int b =1; b<=h[ss]->GetNbinsX();++b) cout << " & " << fixed << setprecision(2) << h[ss]->GetBinContent(b) << "$\\pm$" << h[ss]->GetBinError(b);
    cout << " \\\\" << endl;
    cout << "\\hline" << endl
    << "\\end{tabular}" << endl
    << "\\end{table}" << endl;
    */

}
