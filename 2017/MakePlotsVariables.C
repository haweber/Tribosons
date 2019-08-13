void MakePlotsVariables(){
    
    gStyle->SetOptStat(0);
    bool logy = false;
    
    map<string, TH1D*> h;
    map<string, TH1D*> rat;
    map<string, THStack*> stack;
    vector<string> histonamestemp;
    vector<string> histonames;
    vector<Color_t> cols;
    vector<int> plottogether;
    vector<string> samples;
    /*
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
     */
    TFile *f = TFile::Open("StudyVariousObjects.root");
    histonamestemp.push_back("MET");
    histonamestemp.push_back("MTclosest");
    histonamestemp.push_back("MTsum");
    histonamestemp.push_back("MTmin");
    histonamestemp.push_back("MTmax");
    histonamestemp.push_back("MTlsys");
    histonamestemp.push_back("ptlsum");
    histonamestemp.push_back("ptlsys");
    histonamestemp.push_back("Mlsys");
    histonamestemp.push_back("MllSFOS");
    histonamestemp.push_back("DPhillSFOS");
    histonamestemp.push_back("DRllSFOS");
    histonamestemp.push_back("DPhillmin");
    histonamestemp.push_back("DPhillmax");
    histonamestemp.push_back("Mljclosestmin");
    histonamestemp.push_back("Mljclosestmax");
    histonamestemp.push_back("DPhilsysMET");
    histonamestemp.push_back("HT");
    histonamestemp.push_back("ST");
    histonamestemp.push_back("ptjsys");
    histonamestemp.push_back("MT2ljmin");
    for(unsigned int i = 0; i<histonamestemp.size();++i){
        histonames.push_back(histonamestemp[i]+"_SSee");
        histonames.push_back(histonamestemp[i]+"_SSem");
        histonames.push_back(histonamestemp[i]+"_SSmm");
        histonames.push_back(histonamestemp[i]+"_0SFOS");
        histonames.push_back(histonamestemp[i]+"_1SFOS");
        histonames.push_back(histonamestemp[i]+"_2SFOS");
    }

    
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
    
    
    cols.push_back(kRed);
    //cols.push_back(kRed);
    cols.push_back(kBlue);
    cols.push_back(kCyan+1);
    cols.push_back(kCyan+2);
    cols.push_back(kGreen+2);
    cols.push_back(kGreen+1);
    cols.push_back(kMagenta+1);
    cols.push_back(kMagenta+2);
    cols.push_back(kYellow+1);
    cols.push_back(kYellow+2);
    cols.push_back(kBlack);
    
    /*
     samples.push_back("WWW");
     //samples.push_back("WWWv2");
     samples.push_back("Wjets");
     samples.push_back("TT1l");
     samples.push_back("ZZ");
     samples.push_back("Zjets");
     samples.push_back("TT2l");
     samples.push_back("singleT");
     samples.push_back("ttV");
     samples.push_back("WZ");
     samples.push_back("WW");
     samples.push_back("bg");
     
    cols.push_back(kRed);
    //cols.push_back(kRed);
    cols.push_back(kOrange-2);  //Wjets   //fake LF
    cols.push_back(kOrange-1);  //TT1l    //fake HF
    //cols.push_back(kYellow+2);  //ZZ      //SM Z+ZZ
    //cols.push_back(kYellow+1);  //Zjets   //SM Z//drop this
    cols.push_back(kRed-6); //TT2l    //SM OS2l w HF <-- incorporate sT here if tW
    cols.push_back(kRed-9); //sT      //SM OS2l w/o HF (WW)
    cols.push_back(kAzure-7); //TT2l    //SM OS2l w HF <-- incorporate sT here if tW
    cols.push_back(kAzure-2); //sT      //SM OS2l w/o HF (WW)
    cols.push_back(kGreen-2); //ttV     //SM LL/3l w HF
    cols.push_back(kGreen-3); //WZ      //SM LL/3l w/o HF   kOrange-2
    cols.push_back(kAzure+6);  //WW      //SM ++/-- w/o HF // kAzure + 9 for SM ++/-- with HF
    cols.push_back(kBlack);
*/
    for(unsigned int n = 0; n<histonames.size();++n){
        string bgname = histonames[n] + "_bg";
        string stackname = histonames[n];
        //cout << stackname << endl;
        stack[stackname] = new THStack();
        stack[stackname]->SetName(stackname.c_str());
        for(unsigned int s = 0; s<samples.size()-1; ++s){
            string mapname = histonames[n] + "_" + samples[s];
            h[mapname ]=(TH1D*)f->Get(mapname.c_str());
            h[mapname]->GetXaxis()->SetTitle(stackname.c_str());
            //h[mapname]->Rebin(2);
            if(s==0&&(mapname.find("1SFOS")!=string::npos||mapname.find("2SFOS")!=string::npos||mapname.find("SSee")!=string::npos)){
                h[mapname]->Scale(5.);
            } else if(s==0) h[mapname]->Scale(1.);
            if(mapname.find("_Mjj_")!=string::npos){
                if(!(mapname.find("MjjW")!=string::npos)){
                    h[mapname]->SetBinContent(h[mapname]->GetNbinsX(),0.2*h[mapname]->GetBinContent(h[mapname]->GetNbinsX()));
                    h[mapname]->SetBinError  (h[mapname]->GetNbinsX(),0.2*h[mapname]->GetBinError  (h[mapname]->GetNbinsX()));
                }
            }
            for(int b = 0; b<h[mapname]->GetNbinsX();++b){ if(h[mapname]->GetBinContent(b)<0){h[mapname]->SetBinContent(b,0); h[mapname]->SetBinError(b,0); } }
            if(s==1) h[bgname] = (TH1D*)h[mapname]->Clone(bgname.c_str());
            else if(s>1) h[bgname]->Add(h[mapname],1.);
            if(s==0) {
                h[mapname]->SetLineColor(cols[s]); h[mapname]->SetLineWidth(3);
            }
            else h[mapname]->SetFillColor(cols[s]);
            if(s==1){
                h[bgname]->SetFillColor(cols[samples.size()-1]); h[bgname]->SetFillStyle(3244);
            }
            if(s>=1){
                stack[stackname]->Add(h[mapname],"");
            }
            //cout << h[mapname]->Integral() << endl;
        }
        rat[stackname] = (TH1D*)h[stackname + "_"+samples[0] ]->Clone(stackname.c_str());
        rat[stackname]->Divide(h[stackname + "_bg"]);
    }
    
    TLatex *tLumi = new TLatex(0.95,0.95,"35.9 fb^{-1} (13 TeV)");
    tLumi->SetNDC();
    tLumi->SetTextAlign(31);
    tLumi->SetTextFont(42);
    tLumi->SetTextSize(0.042);
    tLumi->SetLineWidth(2);
    TLatex *tECM = new TLatex(0.95,0.95,"(13 TeV)");
    tECM->SetNDC();
    tECM->SetTextAlign(31);
    tECM->SetTextFont(42);
    tECM->SetTextSize(0.042);
    tECM->SetLineWidth(2);
    //tLumi->Draw();
    TLatex *tCMS = new TLatex(0.185,0.95,"CMS");
    tCMS->SetNDC();
    tCMS->SetTextAlign(11);
    tCMS->SetTextFont(61);
    tCMS->SetTextSize(0.0525);
    tCMS->SetLineWidth(2);
    //tCMS->Draw();
    TLatex *tSim = new TLatex(0.295,0.95,"Supplementary");
    tSim->SetNDC();
    tSim->SetTextAlign(11);
    tSim->SetTextFont(52);
    tSim->SetTextSize(0.042);
    tSim->SetLineWidth(2);
    TLatex *tPrel = new TLatex(0.295,0.95,"Preliminary");
    tPrel->SetNDC();
    tPrel->SetTextAlign(11);
    tPrel->SetTextFont(52);
    tPrel->SetTextSize(0.042);
    tPrel->SetLineWidth(2);
    TLatex *tComment = new TLatex(0.185,0.95,"Last bin scaled x0.2");
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
    

    TLegend *leg1 = new TLegend(0.2,0.67,0.5,0.89,NULL,"brNDC");
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.035);
    leg1->SetLineColor(1);
    leg1->SetLineStyle(1);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1001);
    int legcount(0);
    for(unsigned int s = 1; s<(samples.size()/2);++s){
        leg1->AddEntry(h[histonames[0]+"_"+samples[s] ],samples[s].c_str(), "f");
        legcount = s;
    }
    TLegend *leg2 = new TLegend(0.55,0.67,0.85,0.89,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.035);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1001);
    for(unsigned int s = legcount+1; s<samples.size()-1;++s){
        leg2->AddEntry(h[histonames[0]+"_"+samples[s] ],samples[s].c_str(), "f");
        legcount = s;
    }
    leg2->AddEntry(h[histonames[0]+"_"+samples[0] ],samples[0].c_str(), "l");
    
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
        //cout << stackname << endl;
        stack[stackname]->SetMaximum(1.75*stack[stackname]->GetMaximum());
        stack[stackname]->Draw("hist");
        string bgname = stackname + "_bg";
        h[bgname]->Draw("sameE1");
        string signame = stackname + "_"+ samples[0];
        h[signame]->Draw("histsame");
        string outname = stackname + ".pdf";
        if(logy) outname = "plotsCut/variables/Log/" + outname;
        else     outname = "plotsCut/variables/Lin/" + outname;
        leg1->Draw();
        leg2->Draw();
        tLumi->Draw();
        if(stackname.find("_Mjj_")!=string::npos){
            if(!(stackname.find("MjjW")!=string::npos)){
                tComment->Draw();
            }
        }
        c1->cd();
        TPad *ratiopad = new TPad("ratiopad", "Pad containing the ratio",0,0,1,0.17); //0,0,1,0.26);
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
        rat[stackname]->SetMinimum(0.);
        rat[stackname]->SetMaximum(0.33);
        rat[stackname]->Draw();
        //    for(int i =1; i<=rat[stackname]->GetNbinsX();++i) cout << stackname << " " << i << " " << rat[stackname]->GetBinContent(i) << " +/- " << rat[stackname]->GetBinError(i) << "   " << h[signame]->GetBinContent(i) << " +/- " << h[signame]->GetBinError(i) << endl;
        c1->cd();
        c1->SaveAs(outname.c_str());
        //c1->Clear();
        c1->cd();
    }
}
