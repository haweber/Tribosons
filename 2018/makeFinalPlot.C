void makeFinalPlot(){
    
    bool stacksignal = true;
    gStyle->SetOptStat(0);

    TH1D *WWW     = new TH1D("WWW",    "",6,0,6);
    TH1D *WWorWWW = new TH1D("WWorWWW","",6,0,6);
    TH1D *LL3l    = new TH1D("LL3l",   "",6,0,6);
    TH1D *fakes   = new TH1D("fakes",  "",6,0,6);
    TH1D *gamma   = new TH1D("gamma",  "",6,0,6);
    TH1D *qflips  = new TH1D("qflips", "",6,0,6);
    TH1D *Data    = new TH1D("Data",   "",6,0,6);
    TH1D *bg      = new TH1D("bg",     "",6,0,6);
    //TGraphAsymmErrors *gbgrat   = new TGraphAsymmErrors();
    //TGraphAsymmErrors *gdatarat = new TGraphAsymmErrors();
    Data->SetBinErrorOption(TH1::kPoisson);
    WWW    ->SetBinContent(1,0.61);     WWW    ->SetBinError(1,0.13);
    WWW    ->SetBinContent(2,2.57);     WWW    ->SetBinError(2,0.41);
    WWW    ->SetBinContent(3,3.00);     WWW    ->SetBinError(3,0.29);
    WWW    ->SetBinContent(4,0.82);     WWW    ->SetBinError(4,0.26);
    WWW    ->SetBinContent(5,1.96);     WWW    ->SetBinError(5,0.34);
    WWW    ->SetBinContent(6,0.80);     WWW    ->SetBinError(6,0.27);
    
    WWorWWW->SetBinContent(1,0.61);     WWorWWW->SetBinError(1,0.16);
    WWorWWW->SetBinContent(2,2.19);     WWorWWW->SetBinError(2,0.36);
    WWorWWW->SetBinContent(3,2.38);     WWorWWW->SetBinError(3,0.36);
    WWorWWW->SetBinContent(4,0.11);     WWorWWW->SetBinError(4,0.03);
    WWorWWW->SetBinContent(5,0.18);     WWorWWW->SetBinError(5,0.07);
    WWorWWW->SetBinContent(6,0.10);     WWorWWW->SetBinError(6,0.09);
    
    LL3l    ->SetBinContent(1,1.62);     LL3l   ->SetBinError(1,0.67);
    LL3l    ->SetBinContent(2,2.39);     LL3l   ->SetBinError(2,0.73);
    LL3l    ->SetBinContent(3,5.02);     LL3l   ->SetBinError(3,1.11);
    LL3l    ->SetBinContent(4,0.37);     LL3l   ->SetBinError(4,0.13);
    LL3l    ->SetBinContent(5,11.05);    LL3l   ->SetBinError(5,1.62);
    LL3l    ->SetBinContent(6,9.29);     LL3l   ->SetBinError(6,1.56);

    qflips  ->SetBinContent(1,0.17);     qflips ->SetBinError(1,0.17);
    qflips  ->SetBinContent(2,0.15);     qflips ->SetBinError(2,0.15);
    qflips  ->SetBinContent(3,0.00);     qflips ->SetBinError(3,0.00);
    qflips  ->SetBinContent(4,0.11);     qflips ->SetBinError(4,0.11);
    qflips  ->SetBinContent(5,0.00);     qflips ->SetBinError(5,0.00);
    qflips  ->SetBinContent(6,0.00);     qflips ->SetBinError(6,0.00);
    
    fakes   ->SetBinContent(1,0.66);     fakes  ->SetBinError(1,0.54);
    fakes   ->SetBinContent(2,1.17);     fakes  ->SetBinError(2,0.82);
    fakes   ->SetBinContent(3,4.20);     fakes  ->SetBinError(3,1.77);
    fakes   ->SetBinContent(4,0.48);     fakes  ->SetBinError(4,0.26);
    fakes   ->SetBinContent(5,0.45);     fakes  ->SetBinError(5,0.36);
    fakes   ->SetBinContent(6,0.38);     fakes  ->SetBinError(6,0.26);

    gamma   ->SetBinContent(1,0.32);     gamma  ->SetBinError(1,0.32);
    gamma   ->SetBinContent(2,0.57);     gamma  ->SetBinError(2,0.57);
    gamma   ->SetBinContent(3,0.00);     gamma  ->SetBinError(3,0.00);
    gamma   ->SetBinContent(4,0.00);     gamma  ->SetBinError(4,0.00);
    gamma   ->SetBinContent(5,0.03);     gamma  ->SetBinError(5,0.03);
    gamma   ->SetBinContent(6,0.00);     gamma  ->SetBinError(6,0.00);
    
    for(int i = 0; i< 0; ++i) Data->Fill(0.5);
    for(int i = 0; i< 0; ++i) Data->Fill(1.5);
    for(int i = 0; i< 0; ++i) Data->Fill(2.5);
    for(int i = 0; i< 0; ++i) Data->Fill(3.5);
    for(int i = 0; i< 0; ++i) Data->Fill(4.5);
    for(int i = 0; i< 0; ++i) Data->Fill(5.5);
    
    bg->Add(WWorWWW,1.);
    bg->Add(LL3l   ,1.);
    bg->Add(qflips ,1.);
    bg->Add(fakes  ,1.);
    bg->Add(gamma  ,1.);
    if(stacksignal){
        for(int i = 1; i<=6; ++i) WWW->SetBinContent(i,WWW->GetBinContent(i)+bg->GetBinContent(i));
    }
    
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

    WWW    ->SetLineColor(kRed  ); WWW    ->SetLineWidth(3     ); WWW    ->SetMarkerSize(0);
    WWorWWW->SetLineColor(kBlack); WWorWWW->SetFillColor(2001  ); WWorWWW->SetMarkerSize(0);
    LL3l   ->SetLineColor(kBlack); LL3l   ->SetFillColor(2003  ); LL3l   ->SetMarkerSize(0);
    qflips ->SetLineColor(kBlack); qflips ->SetFillColor(2007  ); qflips ->SetMarkerSize(0);
    gamma  ->SetLineColor(kBlack); gamma  ->SetFillColor(kGray ); gamma  ->SetMarkerSize(0);
    fakes  ->SetLineColor(kBlack); fakes  ->SetFillColor(2005  ); fakes  ->SetMarkerSize(0);
    bg     ->SetLineColor(kBlack); bg     ->SetFillColor(kBlack); bg     ->SetMarkerSize(0); bg     ->SetFillStyle(3544);
    Data   ->SetLineColor(kBlack); Data   ->SetLineWidth(2     ); Data   ->SetMarkerStyle(20);
    
    
    Data   ->GetXaxis()->SetTitleSize(0.0);
    Data   ->GetXaxis()->SetLabelSize(0.06);
    //Data   ->GetXaxis()->SetLabelOffset(1.0);
    Data   ->GetYaxis()->SetTitle("Events");
    Data   ->GetYaxis()->SetTitleSize(0.04);
    Data   ->GetYaxis()->SetLabelSize(0.04);
    Data   ->GetYaxis()->SetTitleOffset(1.10);
    Data   ->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); Data   ->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); Data   ->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); Data   ->GetXaxis()->SetBinLabel(4,"0 SFOS"); Data   ->GetXaxis()->SetBinLabel(5,"1 SFOS"); Data   ->GetXaxis()->SetBinLabel(6,"2 SFOS");
    bg     ->GetXaxis()->SetTitleSize(0.0);
    bg     ->GetXaxis()->SetLabelSize(0.06);
    //bg     ->GetXaxis()->SetLabelOffset(1.00);
    bg     ->GetYaxis()->SetTitle("Events");
    bg     ->GetYaxis()->SetLabelSize(0.04);
    bg     ->GetYaxis()->SetTitleSize(0.04);
    bg     ->GetYaxis()->SetTitleOffset(1.10);
    bg     ->GetXaxis()->SetBinLabel(1,"e^{#pm}e^{#pm}"); bg     ->GetXaxis()->SetBinLabel(2,"e^{#pm}#mu^{#pm}"); bg     ->GetXaxis()->SetBinLabel(3,"#mu^{#pm}#mu^{#pm}"); bg     ->GetXaxis()->SetBinLabel(4,"0 SFOS"); bg     ->GetXaxis()->SetBinLabel(5,"1 SFOS"); bg     ->GetXaxis()->SetBinLabel(6,"2 SFOS");
    
    
    TH1D *bgnoerr = (TH1D*)bg->Clone("bgnoerr");
    for(int i = 1; i<=6; ++i) bgnoerr->SetBinError(i,0.);
    
    TH1D *gbgrat   = (TH1D*)bg  ->Clone("gbgrat");
    for(int i = 1; i<=6; ++i) {
        gbgrat->SetBinContent(i,1);
        gbgrat->SetBinError(  i, bg->GetBinError(i)/bg->GetBinContent(i));
    }
    //gbgrat  ->Divide(bgnoerr);
    //TH1D *gdatarat = (TH1D*)Data->Clone("gdatarat");
    //gdatarat->Divide(bgnoerr);
    //gbgrat  -> TGraphAsymmErrors::Divide(bg,  bgnoerr, "cp");
    //gdatarat-> TGraphAsymmErrors::Divide(Data,bgnoerr, "cp");
    TGraphAsymmErrors *gdatarat = new TGraphAsymmErrors();
    gdatarat->SetName("gdatarat");
    for(int i = 1; i<=6; ++i) {
        gdatarat->SetPoint(      i-1,i-0.5,Data->GetBinContent( i)/bg->GetBinContent(i));
        gdatarat->SetPointEYhigh(i-1,      Data->GetBinErrorUp( i)/bg->GetBinContent(i));
        gdatarat->SetPointEYlow( i-1,      Data->GetBinErrorLow(i)/bg->GetBinContent(i));
    }
    gdatarat   ->SetLineColor(kBlack); gdatarat   ->SetLineWidth(2     ); gdatarat   ->SetMarkerStyle(20);

    THStack *stack = new THStack();
    
    stack->Add(qflips ,"");
    stack->Add(gamma  ,"");
    stack->Add(fakes  ,"");
    stack->Add(LL3l   ,"");
    stack->Add(WWorWWW,"");
 
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
    //TLegend *leg1 = new TLegend(0.1667,0.63,0.85,0.9075,NULL,"brNDC");
    //leg1->SetNColumns(2);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.033);
    leg1->SetLineColor(1);
    leg1->SetLineStyle(1);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1001);
    leg1->AddEntry(WWorWWW, "irreducible",   "f");
    leg1->AddEntry(LL3l   , "lost-leptons/three-leptons",   "f");
    leg1->AddEntry(fakes  , "non-prompt leptons",           "f");
    leg1->AddEntry(gamma  , "#gamma #rightarrow l",         "f");

 
    TLegend *leg2 = new TLegend(0.55,0.63,0.85,0.9075,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.033);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1001);
    
    leg2->AddEntry(qflips , "charge misassignment",         "f");
    leg2->AddEntry(WWW    , "W^{#pm}W^{#pm}W^{#mp} signal (stacked)", "l");
    leg2->AddEntry(Data   , "Data",                        "ep");
    
        TCanvas *c1 = new TCanvas("c1", "",334,192,600,600);
        c1->SetFillColor(0);
        c1->SetBorderMode(0);
        c1->SetBorderSize(2);
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
    

        double max = 35.;
        stack->SetMaximum(max);
        //h[bgname]->Draw("axis");
        stack->Draw("hist");
        stack->SetHistogram(bg);
        stack->Draw("hist");
        bg   ->Draw("sameE2");
        Data ->Draw("sameE0X0");
        WWW  ->Draw("histsame");
    leg1 ->Draw();
    leg2 ->Draw();
    //TLine *splitter = new TLine(3,0,3,max);
    //splitter->SetLineColor(kMagenta);
    //splitter->SetLineWidth(5);
    //splitter->SetLineStyle(3);
    //splitter->Draw();

        tLumi->Draw();

        tCMS->Draw();
        tPrel->Draw();

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
    gbgrat  ->SetMinimum(0.5);
    gbgrat  ->SetMaximum(3.0);
    gdatarat->SetMinimum(0.5);
    gdatarat->SetMaximum(3.0);
        gbgrat  ->GetXaxis()->SetTitle("signal regions");
        gbgrat  ->GetXaxis()->SetTitleSize(0.17);
        gbgrat  ->GetXaxis()->SetTitleOffset(0.70);
        gbgrat  ->GetXaxis()->SetLabelSize(0.0);
        gbgrat  ->GetYaxis()->SetNdivisions(504);
        gbgrat  ->GetYaxis()->SetTitle("data / pred.");
        gbgrat  ->GetYaxis()->SetTitleSize(0.15);
        gbgrat  ->GetYaxis()->SetTitleOffset(0.28);
        gbgrat  ->GetYaxis()->SetLabelSize(0.14);
    gdatarat->GetXaxis()->SetTitle("signal regions");
    gdatarat->GetXaxis()->SetTitleSize(0.17);
    gdatarat->GetXaxis()->SetTitleOffset(0.70);
    gdatarat->GetXaxis()->SetLabelSize(0.0);
    gdatarat->GetYaxis()->SetNdivisions(504);
    gdatarat->GetYaxis()->SetTitle("data / pred.");
    gdatarat->GetYaxis()->SetTitleSize(0.15);
    gdatarat->GetYaxis()->SetTitleOffset(0.28);
    gdatarat->GetYaxis()->SetLabelSize(0.14);
    gbgrat  ->Draw("E2");
    //gdatarat->Draw("sameE0X0");
    gdatarat->Draw("p0Z");
            TLine *rline = new TLine(0.,1.,6.,1.);
            rline->SetLineWidth(2);
            rline->SetLineStyle(7);
            rline->Draw();
        c1->cd();
    c1->SaveAs("ResultPlot.pdf");
        c1->cd();
}
