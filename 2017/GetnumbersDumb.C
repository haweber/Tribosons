void GetnumbersDumb(){
    
    map<string, TH1D*> h;
    vector<string> histonames;
    vector<int> group;
    vector<string> samples;
    TFile *f = TFile::Open("CheckMll.root");
    
    histonames.push_back("Mll_SSee");
    histonames.push_back("Mll_SSem");
    histonames.push_back("Mll_SSmm");
    histonames.push_back("Mll_ee_0SFOS");
    histonames.push_back("Mll_mm_0SFOS");
    histonames.push_back("Mll_SF_0SFOS");
    histonames.push_back("Mll_ee_1SFOS");
    histonames.push_back("Mll_mm_1SFOS");
    histonames.push_back("Mll_SF_1SFOS");
    histonames.push_back("Mll_SFOS_1SFOS");
    histonames.push_back("Mll_ee_2SFOS");
    histonames.push_back("Mll_mm_2SFOS");
    histonames.push_back("Mll_SF_2SFOS");
    histonames.push_back("Mll_SFOS_2SFOS");
    
    samples.push_back("WWW");
    //samples.push_back("WWWv2");
    samples.push_back("WZ");
    samples.push_back("WW");
    samples.push_back("ZZ");
    samples.push_back("ttV");
    samples.push_back("TT2l");
    samples.push_back("TT1l");
    samples.push_back("singleT");
    samples.push_back("Zjets");
    samples.push_back("Wjets");
    samples.push_back("bg");
    
    for(unsigned int n = 0; n<histonames.size();++n){
        string bgname = histonames[n] + "_bg";
        string stackname = histonames[n];
        for(unsigned int s = 0; s<samples.size()-1; ++s){
            string mapname = histonames[n] + "_" + samples[s];
            //cout << mapname << endl;
            h[mapname ]=(TH1D*)f->Get(mapname.c_str());
            for(int b = 0; b<h[mapname]->GetNbinsX();++b){ if(h[mapname]->GetBinContent(b)<0){h[mapname]->SetBinContent(b,0); h[mapname]->SetBinError(b,0); } }
            if(s==1) h[bgname] = (TH1D*)h[mapname]->Clone(bgname.c_str());
            else if(s>1) h[bgname]->Add(h[mapname],1.);
            //cout << h[mapname]->Integral() << endl;
        }
        h[stackname] = (TH1D*)h[stackname + "_"+samples[0] ]->Clone(stackname.c_str());
        h[stackname]->Divide(h[stackname + "_bg"]);
    }

    vector<string> testing;
    vector<float> cut1, cut2;
    testing.push_back("Mll_SSee");        cut1.push_back(0.);   cut2.push_back(-1.);
    testing.push_back("Mll_SSee");        cut1.push_back(20.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSee");        cut1.push_back(30.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSee");        cut1.push_back(40.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSee");        cut1.push_back(50.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSem");        cut1.push_back(0.);   cut2.push_back(-1.);
    testing.push_back("Mll_SSem");        cut1.push_back(20.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSem");        cut1.push_back(30.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSem");        cut1.push_back(40.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSem");        cut1.push_back(50.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSmm");        cut1.push_back(0.);   cut2.push_back(-1.);
    testing.push_back("Mll_SSmm");        cut1.push_back(20.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSmm");        cut1.push_back(30.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSmm");        cut1.push_back(40.);  cut2.push_back(-1.);
    testing.push_back("Mll_SSmm");        cut1.push_back(50.);  cut2.push_back(-1.);
    testing.push_back("Mll_SFOS_1SFOS");  cut1.push_back(0.);   cut2.push_back(-1.);
    testing.push_back("Mll_SFOS_1SFOS");  cut1.push_back(20.);  cut2.push_back(-1.);
    testing.push_back("Mll_SFOS_1SFOS");  cut1.push_back(40.);  cut2.push_back(-1.);
    testing.push_back("Mll_SFOS_2SFOS");  cut1.push_back(0.);   cut2.push_back(-1.);
    testing.push_back("Mll_SFOS_2SFOS");  cut1.push_back(20.);  cut2.push_back(-1.);
    testing.push_back("Mll_SFOS_2SFOS");  cut1.push_back(40.);  cut2.push_back(-1.);

    for(unsigned int n = 0; n<testing.size();++n){
        double errS, errB;
        double sig, bg;
        int b1, b2;
        string  bgname = testing[n] + "_bg";
        string ratname = testing[n];
        string signame = testing[n] + "_"+samples[0];
        b1 = h[signame]->GetXaxis()->FindBin(cut1[n]+0.001);
        if(cut2[n]>0) b2 = h[signame]->GetXaxis()->FindBin(cut2[n]-0.001);
        else          b2 = h[signame]->GetNbinsX();
        sig = h[signame]->IntegralAndError(b1,b2,errS);
        bg = h[bgname]->IntegralAndError(b1,b2,errB);
        cout << ratname << " " << cut1[n] << "-" << cut2[n] << setprecision(3) << " WWW " << sig << " +/- " << errS << " bg " << bg << " +/- " << errB << "   S/B " << sig/bg << " +/- " << sqrt(pow(errS/bg,2)+pow(sig*errB/(bg*bg),2)) << endl;
    }

}
