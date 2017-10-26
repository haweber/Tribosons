#from ROOT import TCanvas  # global import, easy but NOT recommended
import ROOT
from inspect import currentframe, getframeinfo

frameinfo = getframeinfo(currentframe())

addsig = True
twosig = False
data = True

ROOT.gStyle.SetOptStat(0);

f = ROOT.TFile("Check3lCRv2.root","open")

histonames = ["YieldsSR","YieldsCR_SSany_dropMjj"]
axisnames  = ["signal regions","control regions"]
SSorSFOS   = [ 0,               0]

samples  = ["WWW","others","photonfakes"         ,"chargeflips","fakes"    ,"3lLL",       "SSLL"       ,"true3L","trueSS"                 ,"trueWWW"]
legnames = ["WWW","Other" ,"#gamma #rightarrow l","charge flip","jet fakes","lost lepton","lost lepton","3l"    ,"W^{#pm}W^{#pm}(W^{#mp})","W^{#pm}W^{#pm}(W^{#mp})"]
colors   = [ 632,  2012,    920,                   2007,         2005,       2011,         2003,         2003,    2001,                     2001]  
isSS3l   = [ -1,    -1,      0,                     0,            0,          2,            1,            2,       1,                        2]  
if twosig :
    samples.insert(1,"WHtoWWW")
    legnames.insert(1,"WH#rightarrowWWW")
    colors.insert(1,600)
    isSS3l.insert(1,-1)

xx1 = ROOT.TColor(2001, 91/255.,187/255.,241/255.)
xx2 = ROOT.TColor(2003,230/255.,159/255.,  0/255.)
xx3 = ROOT.TColor(2005,245/255.,236/255., 69/255.)
xx4 = ROOT.TColor(2007, 70/255.,109/255.,171/255.)
xx5 = ROOT.TColor(2011,  0/255.,158/255.,115/255.)
xx6 = ROOT.TColor(2012,204/255.,121/255.,167/255.)

stacks = dict()
histos = dict()

for h,a in zip(histonames,axisnames) :
    stacks[h] = ROOT.THStack()
    for s,l,c in zip(samples,legnames,colors) :
        #print h+"_"+s
        histos[h+"_"+s] = f.Get(h+"_"+s)
        if s == "WWW" and addsig :
            histos[h+"_WHtoWWW"] = f.Get(h+"_WHtoWWW")
            histos[h+"_"+s].Add(histos[h+"_WHtoWWW",1.])
        histos[h+"_"+s].GetXaxis().SetTitle(a)
        histos[h+"_"+s].GetYaxis().SetTitle("Events")
        histos[h+"_"+s].GetXaxis().SetTitleSize(0.)
        for b in range(histos[h+"_"+s].GetNbinsX()) :
            if histos[h+"_"+s].GetBinContent(b) < 0 :
                histos[h+"_"+s].SetBinContent(b,0)
                histos[h+"_"+s].SetBinError(b,0)
        if s != "WWW" and s != "WHtoWWW" :
            #print h+"_"+s,c
            histos[h+"_"+s].SetLineColor(c)
            histos[h+"_"+s].SetFillColor(1)
            stacks[h].Add(histos[h+"_"+s])
            if h+"_bg" not in histos :
                histos[h+"_bg"] = histos[h+"_"+s].Clone(h+"_bg")
                histos[h+"_bg"].SetLineColor(ROOT.kBlack)
                histos[h+"_bg"].SetFillColor(ROOT.kBlack)
                histos[h+"_bg"].SetFillStyle(3544)
            else :
                histos[h+"_bg"].Add(histos[h+"_"+s])
        elif s == "WHtoWWW" :
            #print h+"_"+s,c
            histos[h+"_"+s].SetLineColor(c)
            histos[h+"_"+s].SetLineWidth(3)
            histos[h+"_"+s].SetLineStyle(7)
        else : #this is WWW
            #print h+"_"+s,c
            histos[h+"_"+s].SetLineColor(c)
            histos[h+"_"+s].SetLineWidth(3)
            histos[h+"_"+s].SetLineStyle(1)
    if data :
        histos[h+"_Data"] = f.Get(h+"_Data")
        histos[h+"_Data"].SetLineColor(ROOT.kBlack)
        histos[h+"_Data"].SetLineWidth(2)
        histos[h+"_Data"].SetMarkerStyle(20)
        histos[h+"_Ratio"] = histos[h+"_Data"].Clone(h+"_Ratio")
        histos[h+"_Ratio"].Divide(histos[h+"_bg"])
    else :
        histos[h+"_Ratio"] = histos[h+"_WWW"].Clone(h+"_Ratio")
        histos[h+"_Ratio"].Divide(histos[h+"_bg"])
        if twosig :
            histos[h+"_Ratio2"] = histos[h+"_WHtoWWW"].Clone(h+"_Ratio2")
            histos[h+"_Ratio2"].Divide(histos[h+"_bg"])
            
tLumi = ROOT.TLatex(0.95,0.955,"35.9 fb^{-1} (13 TeV)")
tLumi.SetNDC()
tLumi.SetTextAlign(31)
tLumi.SetTextFont(42)
tLumi.SetTextSize(0.042)
tLumi.SetLineWidth(2)
tECM = ROOT.TLatex(0.95,0.955,"(13 TeV)")
tECM.SetNDC()
tECM.SetTextAlign(31)
tECM.SetTextFont(42)
tECM.SetTextSize(0.042)
tECM.SetLineWidth(2)
tCMS = ROOT.TLatex(0.125,0.955,"CMS")
tCMS.SetNDC()
tCMS.SetTextAlign(11)
tCMS.SetTextFont(61)
tCMS.SetTextSize(0.0525)
tCMS.SetLineWidth(2)
tSim = ROOT.TLatex(0.295,0.955,"Simulation")
tSim.SetNDC()
tSim.SetTextAlign(11)
tSim.SetTextFont(52)
tSim.SetTextSize(0.042)
tSim.SetLineWidth(2)
tPrel = ROOT.TLatex(0.295,0.955,"Preliminary")
tPrel.SetNDC()
tPrel.SetTextAlign(11)
tPrel.SetTextFont(52)
tPrel.SetTextSize(0.042)
tPrel.SetLineWidth(2)

#one column 0.2 0.67 0.50 0.89
#two column 0.2 0.67 0.85 0.89
legSS = ROOT.TLegend(0.2,0.65,0.5,0.9025,"","brNDC")
legSS.SetBorderSize(0)
legSS.SetTextSize(0.033)
legSS.SetLineColor(1)
legSS.SetLineStyle(1)
legSS.SetLineWidth(2)
legSS.SetFillColor(0)
legSS.SetFillStyle(1001)
#one column 0.55 0.67 0.85 0.89
#two column 0.20 0.67 0.85 0.89
leg3l = ROOT.TLegend(0.2,0.65,0.85,0.9025,"","brNDC")
leg3l.SetBorderSize(0)
leg3l.SetTextSize(0.033)
leg3l.SetLineColor(1)
leg3l.SetLineStyle(1)
leg3l.SetLineWidth(2)
leg3l.SetFillColor(0)
leg3l.SetFillStyle(1001)
for s,l,i in zip(reversed(samples),reversed(legnames),reversed(isSS3l)) :
    if i == 0 or i == 1 :
        #print histonames[0]+"_"+s,legSS.GetNRows()
        legSS.AddEntry(histos[histonames[0]+"_"+s],l,"f")
    if i == 0 or i == 2 :
        leg3l.AddEntry(histos[histonames[0]+"_"+s],l,"f")
legSS.AddEntry(histos[histonames[0]+"_Data"],"Data","ep")
leg3l.AddEntry(histos[histonames[0]+"_Data"],"Data","ep")
legSS.AddEntry(histos[histonames[0]+"_WWW"],"WWW","l")
leg3l.AddEntry(histos[histonames[0]+"_WWW"],"WWW","l")
if twosig :
    legSS.AddEntry(histos[histonames[0]+"_WHtoWWW"],"WH#rightarrowWWW","l")
    leg3l.AddEntry(histos[histonames[0]+"_WHtoWWW"],"WH#rightarrowWWW","l")
        
for h,a,b in zip(histonames,axisnames,SSorSFOS) :
    c1 = ROOT.TCanvas("c1", "",334,192,600,600)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetBorderSize(2)
    c1.SetTickx(1)
    c1.SetTicky(1)
    c1.SetLeftMargin(0.18)
    c1.SetRightMargin(0.05)
    c1.SetTopMargin(0.07)
    c1.SetBottomMargin(0.15)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    plotpad = ROOT.TPad("plotpad", "Pad containing the overlay plot",0,0.165,1,1)
    plotpad.Draw()
    plotpad.cd()
    plotpad.SetFillColor(0)
    plotpad.SetBorderMode(0)
    plotpad.SetBorderSize(2)
    plotpad.SetTickx(1)
    plotpad.SetTicky(1)
    plotpad.SetLeftMargin(0.12)
    plotpad.SetRightMargin(0.04)
    plotpad.SetTopMargin(0.05)
    #plotpad.SetBottomMargin(0.15)
    plotpad.SetFrameFillStyle(0)
    plotpad.SetFrameBorderMode(0)
    plotpad.SetFrameFillStyle(0)
    plotpad.SetFrameBorderMode(0)
        
    plotpad.cd()
    maximum = max((histos[h+"_bg"].GetBinContent(histos[h+"_bg"].GetMaximumBin() )+0.5*histos[h+"_bg"].GetBinError(histos[h+"_bg"].GetMaximumBin() ) ),(histos[h+"_Data"].GetBinContent(histos[h+"_Data"].GetMaximumBin() )+0.5*histos[h+"_Data"].GetBinError(histos[h+"_Data"].GetMaximumBin() ) ) )*1.667 if data else (histos[h+"_bg"].GetBinContent(h[bgname].GetMaximumBin() )+0.5*histos[h+"_bg"].GetBinError(histos[h+"_bg"].GetMaximumBin()))*1.667
    #print  getframeinfo(currentframe()).lineno, h
    stacks[h].SetMaximum(maximum)
    histos[h+"_bg"].SetMaximum(maximum)
    histos[h+"_Data"].SetMaximum(maximum)
    stacks[h].Draw("hist")
    stacks[h].SetHistogram(histos[h+"_bg"])
    stacks[h].Draw("hist")
    histos[h+"_bg"].Draw("sameE2")
    histos[h+"_Data"].Draw("sameE0X0") if data else True
    histos[h+"_WWW"].Draw("histsame")
    histos[h+"_WHtoWWW"].Draw("histsame") if twosig else True
    outname = "plots/processsplit/" + h + ".pdf"
    #print h,b
    if b == 0 :
        legSS.SetNColumns(1)
        legSS.SetX2(0.5)
        legSS.SetX2NDC(0.5)
        legSS.Draw()
        leg3l.SetNColumns(1)
        leg3l.SetX1(0.55)
        leg3l.SetX1NDC(0.55)
        leg3l.Draw()
        leg3l.SetNColumns(1)
        leg3l.SetX1(0.55)
        leg3l.SetX1NDC(0.55)
        leg3l.Draw()
        splitter = ROOT.TLine(3,0,3,maximum)
        splitter.SetLineColor(616)
        splitter.SetLineWidth(5)
        splitter.SetLineStyle(3)
        splitter.Draw()
    elif b == 1 :
        legSS.SetNColumns(2)
        legSS.SetX2(0.85)
        legSS.SetX2NDC(0.85)
        legSS.Draw()
    elif b == 2 :
        leg3l.SetNColumns(2)
        leg3l.SetX1(0.20)
        leg3l.SetX1NDC(0.20)
        leg3l.Draw()

    tLumi.Draw()
    tCMS.Draw()
    tPrel.Draw()

    c1.cd()
    ratiopad = ROOT.TPad("ratiopad", "Pad containing the ratio",0,0,1,0.21)
    ratiopad.Draw()
    ratiopad.cd()
    ratiopad.SetFillColor(0)
    ratiopad.SetBorderMode(0)
    ratiopad.SetBorderSize(2)
    ratiopad.SetTickx(1)
    ratiopad.SetTicky(1)
    ratiopad.SetLeftMargin(0.12)
    ratiopad.SetRightMargin(0.04)
    ratiopad.SetBottomMargin(0.3)
    ratiopad.SetFrameFillStyle(0)
    ratiopad.SetFrameBorderMode(0)
    ratiopad.SetFrameFillStyle(0)
    ratiopad.SetFrameBorderMode(0)
    if data :
        histos[h+"_Ratio"].SetMinimum(0.)
        histos[h+"_Ratio"].SetMaximum(2.)
    else : 
        histos[h+"_Ratio"].SetMinimum(0.)
        histos[h+"_Ratio"].SetMaximum(0.65)
    histos[h+"_Ratio"].GetXaxis().SetTitle(a)
    histos[h+"_Ratio"].GetXaxis().SetTitleSize(0.16)
    histos[h+"_Ratio"].GetXaxis().SetTitleOffset(0.76)
    histos[h+"_Ratio"].GetXaxis().SetLabelSize(0.0)
    histos[h+"_Ratio"].GetYaxis().SetNdivisions(504)
    histos[h+"_Ratio"].GetYaxis().SetTitle("data / bgd") if data else histos[h+"_Ratio"].GetYaxis().SetTitle("signal / bgd")
    histos[h+"_Ratio"].GetYaxis().SetTitleSize(0.14)
    histos[h+"_Ratio"].GetYaxis().SetTitleOffset(0.28)
    histos[h+"_Ratio"].GetYaxis().SetLabelSize(0.14)
    histos[h+"_Ratio"].Draw()        
    if twosig and not data :
        histos[h+"_Ratio2"].SetMinimum(0.)
        histos[h+"_Ratio2"].SetMaximum(0.65)
        histos[h+"_Ratio2"].GetXaxis().SetTitle(axisnames[n].c_str())
        histos[h+"_Ratio2"].GetXaxis().SetTitleSize(0.16)
        histos[h+"_Ratio2"].GetXaxis().SetTitleOffset(0.76)
        histos[h+"_Ratio2"].GetXaxis().SetLabelSize(0.0)
        histos[h+"_Ratio2"].GetYaxis().SetNdivisions(504)
        histos[h+"_Ratio2"].GetYaxis().SetTitle("signal / bgd")
        histos[h+"_Ratio2"].GetYaxis().SetTitleSize(0.14)
        histos[h+"_Ratio2"].GetYaxis().SetTitleOffset(0.28)
        histos[h+"_Ratio2"].GetYaxis().SetLabelSize(0.14)
        histos[h+"_Ratio2"].Draw("same")
    if data :
        rline = ROOT.TLine(histos[h+"_Ratio"].GetXaxis().GetBinLowEdge(1),1.,histos[h+"_Ratio"].GetXaxis().GetBinLowEdge(histos[h+"_Ratio"].GetNbinsX()+1),1.)
        rline.SetLineWidth(2)
        rline.SetLineStyle(7)
        rline.Draw()
    c1.cd()
    c1.SaveAs(outname)
    c1.cd()
    
for s in stacks :
    stacks[s].Delete()
