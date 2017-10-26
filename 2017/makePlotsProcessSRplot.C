void makePlotsProcessSRplot(){
    
    bool data = true;
    bool twosig = true;
    bool addlostlepto3l = true;
    gStyle->SetOptStat(0);
    bool logy = false;
    map<string, TH1D*> h;
    map<string, TH1D*> rat;
    map<string, THStack*> stack;
    vector<string> histonames;
    vector<string> axisnames;
    vector<int> cols;
    vector<int> plottogether;
    vector<string> samples;
    vector<string> samplesleg;
    vector<int> onlySSSFOS;

    TFile *f = TFile::Open("rootfiles/newLooper.root");
    histonames.push_back("SignalRegion");                  caption.push_back("signal region");                                         onlySSSFOS.push_back(0);
    histonames.push_back("ApplicationRegion");             caption.push_back("application region");                                    onlySSSFOS.push_back(0);
    histonames.push_back("WZControlRegion");               caption.push_back("lost lepton / 3l-with-Z control region");                onlySSSFOS.push_back(0);
    histonames.push_back("SignalRegionPresel");            caption.push_back("signal region - preselection");                          onlySSSFOS.push_back(0);
    histonames.push_back("ApplicationRegionPresel");       caption.push_back("application region - preselection");                     onlySSSFOS.push_back(0);
    histonames.push_back("WZControlRegionPresel");         caption.push_back("lost lepton / 3l-with-Z control region - preselection"); onlySSSFOS.push_back(0);


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
    samples.push_back("WWW");         samplesleg.push_back("WWW signal");                      is3lSS.push_back(0);
////    if(twosig) { samples.push_back("WHtoWWW");     samplesleg.push_back("WH#rightarrowWWW");   is3lSS.push_back(0); }
    //samples.push_back("others");      samplesleg.push_back("Other");                           is3lSS.push_back(0);
    samples.push_back("photonfakes"); samplesleg.push_back("#gamma #rightarrow l");            is3lSS.push_back(0);
    samples.push_back("chargeflips"); samplesleg.push_back("charge misassignment");            is3lSS.push_back(0);
    //samples.push_back("doublefakes"); samplesleg.push_back("double fakes");                    is3lSS.push_back(0);
    samples.push_back("fakes");       samplesleg.push_back("non-prompt leptons");              is3lSS.push_back(0);
    //samples.push_back("FakePred");    samplesleg.push_back("fakes from data");                 is3lSS.push_back(0);
    if(!addlostlepto3l) { samples.push_back("3lLL");        samplesleg.push_back("lost leptons");                     is3lSS.push_back(2);}
    samples.push_back("SSLL");        samplesleg.push_back("lost leptons");                     is3lSS.push_back(1);
    samples.push_back("true3L");      samplesleg.push_back("three leptons");                   is3lSS.push_back(2);
    samples.push_back("trueWWW");     samplesleg.push_back("WWW bg (ttW)");                    is3lSS.push_back(2);
    samples.push_back("trueSS");      samplesleg.push_back("W^{#pm}W^{#pm} bg (DPS,VBS,ttW)"); is3lSS.push_back(1);
    samples.push_back("bg");          samplesleg.push_back("background");                      is3lSS.push_back(-1);
    if(addlostlepto3l){
        for(unsigned int i = 0; i<samples.size();++i){
            if(samples[i].find("true3L")!=string::npos) samplesleg[i] = "three leptons/lost leptons";
        }
    }
    
    cols.push_back(632);//kRed
////    if(twosig) cols.push_back(kBlue);//kRed
    //cols.push_back(2012);//other
    cols.push_back(kGray);
    cols.push_back(2007);
    //cols.push_back(2006);
    cols.push_back(2005);
if(!addlostlepto3l)     cols.push_back(2011);
    cols.push_back(2003);
    cols.push_back(2003);
    cols.push_back(2001);
    cols.push_back(2001);
    cols.push_back(kBlack);
    
    int nsig = 1;
////    if(twosig) ++nsig;
    for(unsigned int n = 0; n<histonames.size();++n){
        string bgname = histonames[n] + "_bg";
        string stackname = histonames[n];
        stack[stackname] = new THStack();
        stack[stackname]->SetName(stackname.c_str());
        string mapname = "";
        if(data) {
            mapname = histonames[n] + "_Data";
            string  mapnameX = histonames[n] + "_X";
            h[mapnameX ]=(TH1D*)f->Get(mapname.c_str());
            h[mapnameX]->SetName(mapnameX.c_str());
            h[mapname] = new TH1D(mapname.c_str(),"",h[mapnameX]->GetNbinsX(), h[mapnameX]->GetXaxis()->GetBinLowEdge(1), h[mapnameX]->GetXaxis()->GetBinLowEdge(h[mapnameX]->GetNbinsX()+1));
            h[mapname]->SetBinErrorOption(TH1::kPoisson);
            for(int i = 1; i<=h[mapname]->GetNbinsX(); ++i){
                //cout << h[mapnameX]->GetBinContent(i) << endl;
                for(int n = 1; n<=h[mapnameX]->GetBinContent(i); ++n){
                    h[mapname]->Fill(i-0.5,1);
                }
            }
            h[mapname]->SetLineColor(kBlack);
            h[mapname]->SetLineWidth(2);
            h[mapname]->SetMarkerStyle(20);
            h[mapname]->GetXaxis()->SetTitleSize(0.0);
            h[mapname]->GetXaxis()->SetLabelSize(0.05);
            h[mapname]->GetYaxis()->SetTitle("Events");
            h[mapname]->GetYaxis()->SetTitleSize(0.05);
            if(onlySSSFOS[0]==1){ h[mapname]->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); }
            if(onlySSSFOS[0]==2){ h[mapname]->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); }
            if(onlySSSFOS[0]==0){ h[mapname]->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(4,"0 SFOS"); h[mapname]->GetXaxis()->SetBinLabel(5,"1 SFOS"); h[mapname]->GetXaxis()->SetBinLabel(6,"2 SFOS"); }
        }
        for(unsigned int s = 0; s<samples.size()-1; ++s){
            mapname = histonames[n] + "_" + samples[s];
            cout << mapname << endl;
            h[mapname ]=(TH1D*)f->Get(mapname.c_str());
            h[mapname]->GetXaxis()->SetTitle(axisnames[n].c_str());
            h[mapname]->GetYaxis()->SetTitle("Events");
            h[mapname]->GetXaxis()->SetTitleSize(0.00);
            h[mapname]->GetYaxis()->SetTitleSize(0.05);
            h[mapname]->GetXaxis()->SetLabelSize(0.05);
            if(addlostlepto3l&&samples[s].find("true3L")!=string::npos){
                string mapname2 =histonames[n] + "_3lLL";
                h[mapname2 ]=(TH1D*)f->Get(mapname2.c_str());
                for(int i = 1; i<=h[mapname]->GetNbinsX();++i) h[mapname2]->GetXaxis()->SetBinLabel(i,h[mapname]->GetXaxis()->GetBinLabel(i));
                h[mapname]->Add(h[mapname2],1.);
            }
         if(onlySSSFOS[0]==1){ h[mapname]->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); }
            if(onlySSSFOS[0]==2){ h[mapname]->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); }
            if(onlySSSFOS[0]==0){ h[mapname]->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); h[mapname]->GetXaxis()->SetBinLabel(4,"0 SFOS"); h[mapname]->GetXaxis()->SetBinLabel(5,"1 SFOS"); h[mapname]->GetXaxis()->SetBinLabel(6,"2 SFOS"); }
            for(int b = 0; b<h[mapname]->GetNbinsX();++b){ if(h[mapname]->GetBinContent(b)<0){h[mapname]->SetBinContent(b,0); h[mapname]->SetBinError(b,0); } }
            if(s==nsig) h[bgname] = (TH1D*)h[mapname]->Clone(bgname.c_str());
            else if(s>nsig) h[bgname]->Add(h[mapname],1.);
            if(s==0) {
                h[mapname]->SetLineColor(cols[s]); h[mapname]->SetLineWidth(3);
                //new twosig
                if(twosig){
                    string mapname2 =histonames[n] + "_WHtoWWW";
                    h[mapname2 ]=(TH1D*)f->Get(mapname2.c_str());
                    for(int i = 1; i<=h[mapname]->GetNbinsX();++i) h[mapname2]->GetXaxis()->SetBinLabel(i,h[mapname]->GetXaxis()->GetBinLabel(i));
                    h[mapname]->Add(h[mapname2],1.);
                }
            }
////            else if(twosig&&s==1) {
////                h[mapname]->SetLineColor(cols[s]); h[mapname]->SetLineWidth(3); h[mapname]->SetLineStyle(7);
////            }
            else { h[mapname]->SetFillColor(cols[s]); h[mapname]->SetLineColor(kBlack); }
            if(s==nsig){
                h[bgname]->SetFillColor(cols[samples.size()-1]); h[bgname]->SetLineColor(cols[samples.size()-1]); h[bgname]->SetFillStyle(3544);
            }
            if(s>=nsig) stack[stackname]->Add(h[mapname],"");
        }
        if(data) rat[stackname] = (TH1D*)h[stackname + "_Data"        ]->Clone(stackname.c_str());
        else     rat[stackname] = (TH1D*)h[stackname + "_"+samples[0] ]->Clone(stackname.c_str());
        rat[stackname]->Divide(h[stackname + "_bg"]);
////        if(twosig&&!data){
////            rat[stackname+"v2"] = (TH1D*)h[stackname + "_"+samples[1] ]->Clone(stackname.c_str());
////            rat[stackname+"v2"]->Divide(h[stackname + "_bg"]);
////        }
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
    TLatex *tPrel = new TLatex(0.225,0.955,"Preliminary");
    tPrel->SetNDC();
    tPrel->SetTextAlign(11);
    tPrel->SetTextFont(52);
    tPrel->SetTextSize(0.042);
    tPrel->SetLineWidth(2);
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
////    if(twosig){
////        leg1->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
////        leg2->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
////        leg2SS->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
////        leg23l->AddEntry(h[histonames[0]+"_"+samples[1] ],samplesleg[1].c_str(), "l");
////    }
    if(data){
        leg1->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
        leg2->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
        leg2SS->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
        leg23l->AddEntry(h[histonames[0]+"_Data" ],"Data", "ep");
    }
    
    for(unsigned int n = 0; n<histonames.size();++n){
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
    
        string stackname = histonames[n];
        string bgname = stackname + "_bg";
        //cout << stackname << endl;
        double max;
        if(data) max = TMath::Max((h[bgname]->GetBinContent(h[bgname]->GetMaximumBin() )+0.5*h[bgname]->GetBinError(h[bgname]->GetMaximumBin() ) ),(h[stackname+"_Data"]->GetBinContent(h[stackname+"_Data"]->GetMaximumBin() )+0.5*h[stackname+"_Data"]->GetBinError(h[stackname+"_Data"]->GetMaximumBin() ) ) )*1.667;
        else max = (h[bgname]->GetBinContent(h[bgname]->GetMaximumBin() )+0.5*h[bgname]->GetBinError(h[bgname]->GetMaximumBin()))*1.667;
        stack[stackname]->SetMaximum(max);
        //h[bgname]->Draw("axis");
        stack[stackname]->Draw("hist");
        stack[stackname]->SetHistogram(h[bgname]);
        stack[stackname]->Draw("hist");
        h[bgname]->Draw("sameE2");
        if(data) h[stackname+"_Data"]->Draw("sameE0X0");
        string signame = stackname + "_"+ samples[0];
        h[signame]->Draw("histsame");
        signame = stackname + "_"+ samples[1];
////        if(twosig) h[signame]->Draw("histsame");
        string outname = stackname + ".pdf";
        if(logy) outname = "plots/processsplit/testV2/log/" + outname;
        else     outname = "plots/processsplit/testV2/" + outname;
        if(onlySSSFOS[n]==0){ leg1  ->Draw(); leg2  ->Draw(); }
        if(onlySSSFOS[n]==1){ leg1SS->Draw(); leg2SS->Draw(); }
        if(onlySSSFOS[n]==2){ leg13l->Draw(); leg23l->Draw(); }
        if(onlySSSFOS[n]==0){
            TLine *splitter = new TLine(3,0,3,max);
            splitter->SetLineColor(kMagenta);
            splitter->SetLineWidth(5);
            splitter->SetLineStyle(3);
            splitter->Draw();
        }
        tLumi->Draw();

        tCMS->Draw();
        if(data) tPrel->Draw();
        else     tSim ->Draw();

        c1->cd();
        TPad *ratiopad = new TPad("ratiopad", "Pad containing the ratio",0,0,1,0.21); //0,0,1,0.26);
        ratiopad->Draw();
        ratiopad->cd();
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
            rat[stackname]->SetMinimum(0.5);
            rat[stackname]->SetMaximum(2.-0.5);
        }
        else {
            rat[stackname]->SetMinimum(0.);
            rat[stackname]->SetMaximum(0.65);
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
        rat[stackname]->Draw("E0X0");
/*        if(twosig&&!data){
            rat[stackname+"v2"]->SetMinimum(0.);
            rat[stackname+"v2"]->SetMaximum(0.65);
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
        }*/
        if(data){
            TLine *rline = new TLine(rat[stackname]->GetXaxis()->GetBinLowEdge(1),1.,rat[stackname]->GetXaxis()->GetBinLowEdge(rat[stackname]->GetNbinsX()+1),1.);
            rline->SetLineWidth(2);
            rline->SetLineStyle(7);
            rline->Draw();
        }
        c1->cd();
        c1->SaveAs(outname.c_str());
        c1->cd();
    }
}
