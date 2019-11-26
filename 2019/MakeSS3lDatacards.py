import ROOT
from inspect import currentframe, getframeinfo
import csv
import math
from optparse import OptionParser
import re
import numpy as np
import sys
import os
import cStringIO
import subprocess
from BackgroundOtherPredictionHistogramAndPrinter import *

def getyieldanderror(hS,mybin):
    if mybin < 0:
        return -4.,-4
    if mybin > hS.GetNbinsX():
        return -5.,-5
    binnum = hS.FindBin(mybin)
    sig =  hS.GetBinContent(binnum)
    err =  hS.GetBinError(binnum)
    return sig,err


def GetYieldsAndError(histos,bins,corrsys,uncorrsys):
    yall,yerrC,yerrU = [],[],[]
    for h,b in zip(histos,bins):
        y,e = getyieldanderror(h,b)
        if y>0:
            yall.append(y)
            yerrC.append(1.+corrsys)
            yerrU.append(1.+math.sqrt(e*e+pow(uncorrsys*y,2))/y)
        else:
            yall.append(0)
            yerrC.append(1.)
            yerrU.append(1.)
    return yall,yerrC,yerrU

def PrintAllDatacards():
    histos = dict()
    LoadAllPredictions(histos)
    for i in range(1,13):
        PrintDatacards(histos,i)
        
def PrintSingleDatacards(histos,mybin):
    histos = dict()
    LoadAllPredictions(histos)
    PrintDatacards(histos,mybin)
    
def PrintDatacards(h,i):#h = histos, i = mybin
    #systematicsSignal      = ["LepSF","TrigSF","Pileup","Qsq","PDF","JES","AlphaS","BTagLF","BTagHF"]
    systematicsSignal      = ["Pileup","Qsq","PDF","JES","AlphaS","BTagLF","BTagHF"]
    #systematicsWZupdn      = ["LepSF","TrigSF","Pileup","Qsq","PDF","JES","OtherErr","DataErr"]
    systematicsWZupdncor   = ["Pileup","Qsq","PDF","JES"]
    systematicsWZupdn      = ["OtherErr","DataErr"]
    systematicsWZupcor     = ["WZscaled","ttZscaled"]
    systematicsWZup        = ["TFstat"]
    systematicsFakeupdncor = ["FakeRateEl","FakeRateMu","FakeClosureEl","FakeClosureMu"]
    systematicsFakeupdn    = ["Data","Purity"]
    systematicsVBSupcor    = ["Syst"]
    systematicsVBSup       = ["SimStat"]
    systematicsttWupcor    = ["Syst"]
    systematicsttWup       = ["SimStat"]
    systematicsQflipupcor  = ["Syst"]
    systematicsQflipup     = ["SimStat"]
    systematicsPhFakeupcor = ["Syst2016","Syst2017","Syst2018"]
    systematicsPhFakeup    = ["SimStat"]


    channelname = "ch"+('%02d' % i)
    r = "%.4f"
    l = 16
    datacard = cStringIO.StringIO()
    datacard.write("imax   1  number of channels\n")
    datacard.write("jmax   5  number of backgrounds\n")
    datacard.write("kmax   *  number of nuisance parameters\n")
    datacard.write("------------\n")
    datacard.write("shapes *            "+channelname+"  FAKE\n")
    datacard.write("------------\n")
    datacard.write("# observations \n")
    #datacard.write("bin                 "+channelname+"\n")
    datacard.write("bin".ljust(2*l)+channelname+"\n")
    #datacard.write("observation         ")
    datacard.write("observation".ljust(2*l))
    datacard.write("0")#dummy
    datacard.write(" \n")
    datacard.write("------------\n")
    datacard.write("#now we list all expected number of events\n")
    #datacard.write("bin                 ")
    datacard.write("bin".ljust(2*l))
    datacard.write(channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname+"            "+channelname)
    datacard.write(" \n")
    #datacard.write("process             ")
    datacard.write("process".ljust(2*l))
    datacard.write("WWWon           WWWh            WWZon           WWZh            WZZon           WZZh            ZZZon           ZZZh            LostLepton      Prompt          Fakes           ChargeFlip      GammaFakes  ")
    datacard.write(" \n")
    #datacard.write("process             ")
    datacard.write("process".ljust(2*l))
    datacard.write("-7              -6              -5              -4              -3              -2              -1              0               1               2               3               4               5               ")
    datacard.write(" \n")
    #datacard.write("rate                ")
    datacard.write("rate".ljust(2*l))
    datacard.write(    st(h["WWWonshellYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["WWWhiggsYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["WWZonshellYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["WWZhiggsYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["WZZonshellYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["WZZhiggsYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["ZZZonshellYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["ZZZhiggsYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["WZYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["ttWYield"].GetBinContent(i)+h["VBSYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["FakeYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["QFlipYield"].GetBinContent(i),r).ljust(l)+ \
                       st(h["PhotonFakeYield"].GetBinContent(i),r).ljust(l) )
    datacard.write(" \n")
    datacard.write("------------\n")
    datacard.write(printSignalUncertainty(h,i,"SimStat",False,False,r,l))
    for sys in systematicsSignal:
        datacard.write(printSignalUncertainty(h,i,sys,True,True,r,l))
    #make background function
    #def printBkgUncertainty(h,i,sys,WZ,IrrVBS,IrrttW,Fake,Qflip,PhFake,updown=True,cor=True,p="%.5f",sp=12):
    if i == 10:
        datacard.write(printBkgUncertainty(h,i,"Total",True ,False,False,False,False,False,True ,True ,r,l))   #WZ systematics for 0 SFOS     
    for sys in systematicsWZupdncor:
        datacard.write(printBkgUncertainty(h,i,sys,True ,False,False,False,False,False,True ,True ,r,l))
    for sys in systematicsWZupdn:
        datacard.write(printBkgUncertainty(h,i,sys,True ,False,False,False,False,False,True ,False,r,l))
    for sys in systematicsWZupcor:
        datacard.write(printBkgUncertainty(h,i,sys,True ,False,False,False,False,False,False,True ,r,l))
    for sys in systematicsWZup:
        datacard.write(printBkgUncertainty(h,i,sys,True ,False,False,False,False,False,False,False,r,l))
    for sys in systematicsVBSupcor:
        datacard.write(printBkgUncertainty(h,i,sys,False,True ,False,False,False,False,False,True ,r,l))
    for sys in systematicsVBSup:
        datacard.write(printBkgUncertainty(h,i,sys,False,True ,False,False,False,False,False,False,r,l))
    for sys in systematicsttWupcor:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,True ,False,False,False,False,True ,r,l))
    for sys in systematicsttWup:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,True ,False,False,False,False,False,r,l))
    for sys in systematicsFakeupdncor:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,False,True ,False,False,True ,True ,r,l))
    for sys in systematicsFakeupdn:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,False,True ,False,False,True ,False,r,l))
    for sys in systematicsQflipupcor:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,False,False,True ,False,False,True ,r,l))
    for sys in systematicsQflipup:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,False,False,True ,False,False,False,r,l))
    for sys in systematicsPhFakeupcor:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,False,False,False,True ,False,True ,r,l))
    for sys in systematicsPhFakeup:
        datacard.write(printBkgUncertainty(h,i,sys,False,False,False,False,False,True ,False,False,r,l))

    datacardname = "datacards/datacardSS3l_bin"+channelname+".txt"
    datacardtxt = open(datacardname,"w")
    datacardtxt.write(datacard.getvalue())
    datacardtxt.close()
    datacard.close()
    return True

def printBkgUncertainty(h,i,sys,WZ,IrrVBS,IrrttW,Fake,Qflip,PhFake,updown=True,cor=True,p="%.5f",sp=12):
    channelname = ('%02d' % i)
    process = "Bkg"
    outstring = ""
    Irr = (IrrVBS or IrrttW)
    if ( WZ) and (not Irr) and (not Fake) and (not Qflip) and (not PhFake): process = "LL"
    if (not WZ) and ( Irr) and (not Fake) and (not Qflip) and (not PhFake): process = "Prmpt"
    if (not WZ) and (not Irr) and ( Fake) and (not Qflip) and (not PhFake): process = "Fks"
    if (not WZ) and (not Irr) and (not Fake) and ( Qflip) and (not PhFake): process = "Qflp"
    if (not WZ) and (not Irr) and (not Fake) and (not Qflip) and ( PhFake): process = "PhFks"
    if process == "Prmpt":
        if IrrVBS:
            process = "VBS"
        else:
            process = "ttW"
    outstring = ""
    doreturn = False
    if cor:
        #outstring = (process+sys + " lnN     ").ljust(sp+4)
        outstring = (process+sys + " lnN ").ljust(sp*2)
    else:
        #outstring = (process+sys+channelname+" lnN   ").ljust(sp+4)
        outstring = (process+sys+channelname+" lnN ").ljust(sp*2)
    #return outstring
    if updown:
        outstring = outstring + "-               -               -               -               -               -               -               -               "
        if WZ:
            if getUpDnVar(h,"WZ","Yield",sys+"Up",sys+"Down",i,p,sp) != (st(1.,p)+"/"+st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpDnVar(h,"WZ","Yield",sys+"Up",sys+"Down",i,p,sp)
        else:
            outstring = outstring + "-               "
        if Irr:
            if IrrVBS:
                outstring = outstring #+ getUpDnVarPrompt(h,"VBS","ttW","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp)
            else:
                outstring = outstring #+ getUpDnVarPrompt(h,"ttW","VBS","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp)
        else:
            outstring = outstring + "-               "
        if Fake:
            if getUpDnVar(h,"Fake","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) != (st(1.,p)+"/"+st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpDnVar(h,"Fake","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp)
        else:
            outstring = outstring + "-               "
        if Qflip:
            if getUpVar(h,"QFlip","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) != (st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpVar(h,"QFlip","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp)
        else:
            outstring = outstring + "-               "
        if PhFake:
            if getUpVar(h,"PhotonFake","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) != (st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpVar(h,"PhotonFake","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp)
        else:
            outstring = outstring + "-               "
    else:
        outstring = outstring + "-               -               -               -               -               -               -               -               "
        if WZ:
            if getUpVar(h,"WZ","Yield",sys,i,p,sp) != (st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpVar(h,"WZ","Yield",sys,i,p,sp)
        else:
            outstring = outstring + "-               "
        if Irr:
            if IrrVBS:
                if getUpVarPrompt(h,"VBS","ttW","Yield",sys+"Err",i,p,sp) != (st(1.,p)).ljust(sp):
                    doreturn = True
                outstring = outstring + getUpVarPrompt(h,"VBS","ttW","Yield",sys+"Err",i,p,sp)
            else:
                if getUpVarPrompt(h,"ttW","VBS","Yield",sys+"Err",i,p,sp) != (st(1.,p)).ljust(sp):
                    doreturn = True
                outstring = outstring + getUpVarPrompt(h,"ttW","VBS","Yield",sys+"Err",i,p,sp)
        else:
            outstring = outstring + "-               "
        if Fake:
            if getUpVar(h,"Fake","Yield",sys+"Err",i,p,sp) != (st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpVar(h,"Fake","Yield",sys+"Err",i,p,sp)
        else:
            outstring = outstring + "-               "
        if Qflip:
            if getUpVar(h,"QFlip","Yield",sys+"Err",i,p,sp) != (st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpVar(h,"QFlip","Yield",sys+"Err",i,p,sp)
        else:
            outstring = outstring + "-               "
        if PhFake:
            if getUpVar(h,"PhotonFake","Yield",sys+"Err",i,p,sp) != (st(1.,p)).ljust(sp):
                doreturn = True
            outstring = outstring + getUpVar(h,"PhotonFake","Yield",sys+"Err",i,p,sp)
        else:
            outstring = outstring + "-               "
    outstring = outstring + "\n"
    if doreturn:
        return outstring
    return ""


def printSignalUncertainty(h,i,sys,updown=True,cor=True,p="%.5f",sp=12):
    channelname = ('%02d' % i)
    outstring = ""
    if sys == "SimStat":
        #outstring = ("SigSimStat"+channelname+" lnN   ").ljust(sp+4) + \
        outstring = ("SigSimStat"+channelname+" lnN ").ljust(sp*2) + \
          getUpVar(h,"WWWonshell","Yield","SimStatErr",i,p,sp) + \
          getUpVar(h,"WWWhiggs"  ,"Yield","SimStatErr",i,p,sp) + \
          getUpVar(h,"WWZonshell","Yield","SimStatErr",i,p,sp) + \
          getUpVar(h,"WWZhiggs"  ,"Yield","SimStatErr",i,p,sp) + \
          getUpVar(h,"WZZonshell","Yield","SimStatErr",i,p,sp) + \
          getUpVar(h,"WZZhiggs"  ,"Yield","SimStatErr",i,p,sp) + \
          getUpVar(h,"ZZZonshell","Yield","SimStatErr",i,p,sp) + \
          getUpVar(h,"ZZZhiggs"  ,"Yield","SimStatErr",i,p,sp) + \
          "-               -               -               -               -               \n"
    else:
        if cor:
            #outstring = ("Sig"+sys + " lnN     ").ljust(sp+4)
            outstring = ("Sig"+sys + " lnN ").ljust(sp*2)
        else:
            #outstring = ("Sig"+sys+channelname+" lnN   ").ljust(sp+4)
            outstring = ("Sig"+sys+channelname+" lnN ").ljust(sp*2)
        #return outstring
        if updown:
            outstring = outstring + \
              getUpDnVar(h,"WWWonshell","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              getUpDnVar(h,"WWWhiggs"  ,"Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              getUpDnVar(h,"WWZonshell","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              getUpDnVar(h,"WWZhiggs"  ,"Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              getUpDnVar(h,"WZZonshell","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              getUpDnVar(h,"WZZhiggs"  ,"Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              getUpDnVar(h,"ZZZonshell","Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              getUpDnVar(h,"ZZZhiggs"  ,"Yield",sys+"ErrUp",sys+"ErrDown",i,p,sp) + \
              "-               -               -               -               -               \n"
        else:
            outstring = outstring + \
              getUpVar(h,"WWWonshell","Yield",sys+"Err",i,p,sp) + \
              getUpVar(h,"WWWhiggs"  ,"Yield",sys+"Err",i,p,sp) + \
              getUpVar(h,"WWZonshell","Yield",sys+"Err",i,p,sp) + \
              getUpVar(h,"WWZhiggs"  ,"Yield",sys+"Err",i,p,sp) + \
              getUpVar(h,"WZZonshell","Yield",sys+"Err",i,p,sp) + \
              getUpVar(h,"WZZhiggs"  ,"Yield",sys+"Err",i,p,sp) + \
              getUpVar(h,"ZZZonshell","Yield",sys+"Err",i,p,sp) + \
              getUpVar(h,"ZZZhiggs"  ,"Yield",sys+"Err",i,p,sp) + \
              "-               -               -               -               -               \n"
    return outstring

        
def getUpVar(h,sample,central,up,mybin,prec,sp):
    #up is total up error
    c = h[sample+central].GetBinContent(mybin)
    u = h[sample+up].GetBinContent(mybin)
    if c <= 0:
        return ""
    else:
        return st(1.+u/c,prec).ljust(sp)

def getUpVarPrompt(h,sample,sampletwo,central,up,mybin,prec,sp):#hack
    #up is total up error
    c = h[sample+central].GetBinContent(mybin)+h[sampletwo+central].GetBinContent(mybin)
    u = h[sample+up].GetBinContent(mybin)
    if c <= 0:
        return ""
    else:
        return st(1.+u/c,prec).ljust(sp)
    
def getUpDnVar(h,sample,central,up,dn,mybin,prec,sp):
    #up is total up error
    #dn is total dn error
    c = h[sample+central].GetBinContent(mybin)
    u = h[sample+up].GetBinContent(mybin)
    d = h[sample+dn].GetBinContent(mybin)
    if c <= 0:
        return ""
    else:
        return (st(1./(1.+d/c),prec)+"/"+st(1.+u/c,prec)).ljust(sp)

def getUppVar(h,sample,central,upp,mybin,prec,sp):
    #upp is up variation = central + uperr
    c = h[sample+central].GetBinContent(mybin)
    u = h[sample+upp].GetBinContent(mybin)-c
    if c <= 0:
        return ""
    else:
        return st(1.+u/c,prec).ljust(12)

def getUppDnnVar(h,sample,central,upp,dnn,mybin,prec,sp):
    #upp is up variation = central + up error
    #dnn is dn variation = central - dn error
    c = h[sample+central].GetBinContent(mybin)
    u = h[sample+upp].GetBinContent(mybin) - c
    d = c - h[sample+dnn].GetBinContent(mybin)
    if c <= 0:
        return ""
    else:
        return (st(1./(1.+d/c),prec)+"/"+st(1.+u/c,prec)).ljust(sp)

PrintAllDatacards()

"""
def makeDummyDataCard(path,bins):
    ###hS,hP,hL,hF,hQ,hG = LoadHistos(path)
    fS = ROOT.TFile.Open(path+"signal.root")
    #fS = ROOT.TFile.Open(path+"signal_private.root")
    fP = ROOT.TFile.Open(path+"prompt.root")
    fL = ROOT.TFile.Open(path+"lostlep.root")
    fF = ROOT.TFile.Open(path+"fakes.root")
    #fF = ROOT.TFile.Open(path+"ddfakes.root")
    fQ = ROOT.TFile.Open(path+"qflip.root")
    fG = ROOT.TFile.Open(path+"photon.root")
    #print(fS,fP,fL,fF,fQ,fG)
    hS = []
    hP = []
    hL = []
    hF = []
    hQ = []
    hG = []
    histnames = []
    if len(bins)==12:
        histnames = ["SRSSeeMjjIn__CutsSS","SRSSemMjjIn__CutsSS","SRSSmmMjjIn__CutsSS","SRSSeeMjjOut__CutsSS","SRSSemMjjOut__CutsSS","SRSSmmMjjOut__CutsSS","SRSS1Jee1JPre__CutsSS","SRSS1Jem1JPre__CutsSS","SRSS1Jmm1JPre__CutsSS","SR0SFOS__CutsSFOS","SR1SFOS__CutsSFOS","SR2SFOS__CutsSFOS"]
    else:
        histnames = ["SRSSeeMjjIn__CutsSS","SRSSemMjjIn__CutsSS","SRSSmmMjjIn__CutsSS","SRSSeeMjjOut__CutsSS","SRSSemMjjOut__CutsSS","SRSSmmMjjOut__CutsSS","SRSS1Jee1JPre__CutsSS","SRSS1Jem1JPre__CutsSS","SRSS1Jmm1JPre__CutsSS","SR0SFOSeem__CutsSFOS","SR0SFOSemm__CutsSFOS","SR1SFOS__CutsSFOS","SR2SFOS__CutsSFOS"]
    for hn in histnames:
        h = fS.Get(hn)
        #print (hn,h)
        hh = h.Clone(hn+"Sig")
        #print (hn+"Sig",hh)
        hS.append(hh)
        #print (hS)
        h = fP.Get(hn)
        hh = h.Clone(hn+"Prompt")
        hP.append(hh)
        h = fL.Get(hn)
        hh = h.Clone(hn+"LL")
        hL.append(hh)
        h = fF.Get(hn)
        hh = h.Clone(hn+"Fake")
        hF.append(hh)
        h = fQ.Get(hn)
        hh = h.Clone(hn+"Qflip")
        hQ.append(hh)
        h = fG.Get(hn)
        hh = h.Clone(hn+"Gamma")
        hG.append(hh)
    #print hS
    sig, sigerrC, sigerrU  = GetYieldsAndError(hS,bins,0.10,0.10)
    pr,  prerrC,  prerrU   = GetYieldsAndError(hP,bins,0.15,0.10)
    LL,  LLerrC,  LLerrU   = GetYieldsAndError(hL,bins,0.15,0.15)
    fake,fakeerrC,fakeerrU = GetYieldsAndError(hF,bins,0.30,0.50)
    Qfl, QflerrC, QflerrU  = GetYieldsAndError(hQ,bins,0.01,0.50)
    gam, gamerrC, gamerrU  = GetYieldsAndError(hG,bins,0.11,0.50)
    dummydata = []
    for i in range(0,len(sig)):
        dummydata.append(int(sig[i]+pr[i]+LL[i]+fake[i]+Qfl[i]+gam[i]))
    
    fsp = 50
    scg = 25
    prc = "%.5f"
    datacard = cStringIO.StringIO()
    datacard.write("# My significance experiment\n")
    if len(bins)==12:
        datacard.write("imax   12  number of channels\n")
    else:
        datacard.write("imax   13  number of channels\n")
    datacard.write("jmax   5  number of backgrounds\n")
    datacard.write("kmax   *  number of nuisance parameters\n")
    datacard.write("------------\n")
    datacard.write("shapes *         ch01  FAKE\n")
    datacard.write("shapes *         ch02  FAKE\n")
    datacard.write("shapes *         ch03  FAKE\n")
    datacard.write("shapes *         ch04  FAKE\n")
    datacard.write("shapes *         ch05  FAKE\n")
    datacard.write("shapes *         ch06  FAKE\n")
    datacard.write("shapes *         ch07  FAKE\n")
    datacard.write("shapes *         ch08  FAKE\n")
    datacard.write("shapes *         ch09  FAKE\n")
    datacard.write("shapes *         ch10  FAKE\n")
    datacard.write("shapes *         ch11  FAKE\n")
    datacard.write("shapes *         ch12  FAKE\n")
    if len(bins)==13:
        datacard.write("shapes *         ch13  FAKE\n") 
    datacard.write("------------\n")
    datacard.write("# observations \n")
    if len(bins)==12:
        datacard.write("bin             ch01".rjust(12)+"ch02".rjust(12)+"ch03".rjust(12)+"ch04".rjust(12)+"ch05".rjust(12)+"ch06".rjust(12)+"ch07".rjust(12)+"ch08".rjust(12)+"ch09".rjust(12)+"ch10".rjust(12)+"ch11".rjust(12)+"ch12".rjust(12)+"\n")
    else:
        datacard.write("bin             ch01".rjust(12)+"ch02".rjust(12)+"ch03".rjust(12)+"ch04".rjust(12)+"ch05".rjust(12)+"ch06".rjust(12)+"ch07".rjust(12)+"ch08".rjust(12)+"ch09".rjust(12)+"ch10".rjust(12)+"ch11".rjust(12)+"ch12".rjust(12)+"ch13".rjust(12)+"\n")
    datacard.write("observation     ")
    for i in range(0,len(sig)):
        datacard.write(str(int(dummydata[i])).ljust(12))
    datacard.write(" \n")
    datacard.write("------------\n")
    datacard.write("#now we list all expected number of events\n")
    datacard.write("bin             ")
    for i in range(0,len(sig)):
        datacard.write("ch"+str(i+1).zfill(2)+"        ch"+str(i+1).zfill(2)+"        ch"+str(i+1).zfill(2)+"        ch"+str(i+1).zfill(2)+"        ch"+str(i+1).zfill(2)+"        ch"+str(i+1).zfill(2)+"        ")
    datacard.write(" \n")
    datacard.write("process         ")
    for i in range(0,len(sig)):
        datacard.write("sig         Prompt      LostLepton  Fakes       ChargeFlip  GammaFakes  ")
    datacard.write(" \n")
    datacard.write("process         ")
    for i in range(0,len(sig)):
        datacard.write("0           1           2           3           4           5           ")
    datacard.write(" \n")
    datacard.write("rate            ")
    for i in range(0,len(sig)):
        datacard.write(str(prc % sig[i]).ljust(12)+str(prc % pr[i]).ljust(12)+str(prc % LL[i]).ljust(12)+str(prc % fake[i]).ljust(12)+str(prc % Qfl[i]).ljust(12)+str(prc % gam[i]).ljust(12))
    datacard.write(" \n")
    datacard.write("------------\n")
    prc = "%.4f"
    for i in range(0,len(sig)):
        datacard.write("SigUnCorr"+str(i+1).zfill(2)+" lnN ")
        for j in range(0,len(sig)):
            if i == j:
                datacard.write(str(prc % sigerrU[i]).ljust(12)+"-           -           -           -           -           ")
            else:
                datacard.write("-           -           -           -           -           -           ")
        datacard.write(" \n")
    datacard.write("SigCorr     lnN ")
    for i in range(0,len(sig)):
        datacard.write(str(prc % sigerrC[i]).ljust(12)+"-           -           -           -           -           ")
    datacard.write(" \n")
    for i in range(0,len(sig)):
        datacard.write("PrUnCorr"+str(i+1).zfill(2)+"  lnN ")
        for j in range(0,len(sig)):
            if i == j:
                datacard.write("-           "+str(prc % prerrU[i]).ljust(12)+"-           -           -           -           ")
            else:
                datacard.write("-           -           -           -           -           -           ")
        datacard.write(" \n")
    datacard.write("PrCorr      lnN ")
    for i in range(0,len(sig)):
        datacard.write("-           "+str(prc % prerrC[i]).ljust(12)+"-           -           -           -           ")
    datacard.write(" \n")
    for i in range(0,len(sig)):
        datacard.write("LLUnCorr"+str(i+1).zfill(2)+"  lnN ")
        for j in range(0,len(sig)):
            if i == j:
                datacard.write("-           -           "+str(prc % LLerrU[i]).ljust(12)+"-           -           -           ")
            else:
                datacard.write("-           -           -           -           -           -           ")
        datacard.write(" \n")
    datacard.write("LLCorr      lnN ")
    for i in range(0,len(sig)):
        datacard.write("-           -           "+str(prc % LLerrC[i]).ljust(12)+"-           -           -           ")
    datacard.write(" \n")
    for i in range(0,len(sig)):
        datacard.write("FakUnCorr"+str(i+1).zfill(2)+" lnN ")
        for j in range(0,len(sig)):
            if i == j:
                datacard.write("-           -           -           "+str(prc % fakeerrU[i]).ljust(12)+"-           -           ")
            else:
                datacard.write("-           -           -           -           -           -           ")
        datacard.write(" \n")
    datacard.write("FakeCorr    lnN ")
    for i in range(0,len(sig)):
        datacard.write("-           -           -           "+str(prc % fakeerrC[i]).ljust(12)+"-           -           ")
    datacard.write(" \n")
    for i in range(0,len(sig)):
        datacard.write("QflUnCorr"+str(i+1).zfill(2)+" lnN ")
        for j in range(0,len(sig)):
            if i == j:
                datacard.write("-           -           -           -           "+str(prc % QflerrU[i]).ljust(12)+"-           ")
            else:
                datacard.write("-           -           -           -           -           -           ")
        datacard.write(" \n")
    datacard.write("QflCorr     lnN ")
    for i in range(0,len(sig)):
        datacard.write("-           -           -           -           "+str(prc % QflerrC[i]).ljust(12)+"-           ")
    datacard.write(" \n")
    for i in range(0,len(sig)):
        datacard.write("GamUnCorr"+str(i+1).zfill(2)+" lnN ")
        for j in range(0,len(sig)):
            if i == j:
                datacard.write("-           -           -           -           -           "+str(prc % gamerrU[i]).ljust(12))
            else:
                datacard.write("-           -           -           -           -           -           ")
        datacard.write(" \n")
    datacard.write("GamCorr     lnN ")
    for i in range(0,len(sig)):
        datacard.write("-           -           -           -           -           "+str(prc % gamerrC[i]).ljust(12))
    datacard.write(" \n")
    
    datacardname = "dummycard.txt"
    datacardtxt = open(datacardname,"w")
    datacardtxt.write(datacard.getvalue())
    datacardtxt.close()
    datacard.close()
    return True
        

path = "hists/combineyearsLoose_v5.3.0/reOptimize/"
#path = "hists/Loose2016_v5.3.0/reOptimize/"
#bins = [0,0,0,0,0,0,0,0,0,0,0,0] #preselection, no b-tagging
#bins = [1,1,1,1,1,1,1,1,1,1,1,1] #preselection, NB = 0
#bins = [2,2,2,2,2,2,2,2,2,2,2,2] #preselection, NB = NSB = 0
#bins = [193,199,121,193,199,193,   1,   1,   1,26032,28462,28461] #2016 selection
#bins = [ 56, 56, 50, 56, 56, 50,6056,6056,6056,    5,29813,29813] #2018 previous optimization
#bins = [18073,66386,2405,337,45031,14198,66,4065,6963,14,30356,13886] #current perfect selection, 0 SFOS one category
#bins = [18073,66386,2405,337,45031,14198,66,4065,6963,2,14,30356,13886] #current perfect selection, 0 SFOS eem + emm splitted

#bins = [296,296,290,296,296,266,3360,3366,3360,5240,17126,17126] #sensible selection, 0 SFOS one category
#bins = [296,296,290,290,290,266,3360,3366,3360,5240,17126,17126] #sensible selection, 0 SFOS one category
#bins = [296,296,266,296,296,266,3360,3366,3360,5240,17126,17126] #sensible selection, 0 SFOS one category
#bins = [296,296,266,290,290,266,3360,3366,3360,5240,17126,17126] #sensible selection, 0 SFOS one category

#bins = [296,296,290,296,296,266,3360,3366,3360,4970,17126,17126] #sensible selection, 0 SFOS one category
#bins = [296,296,290,290,290,266,3360,3366,3360,4970,17126,17126] #sensible selection, 0 SFOS one category
#bins = [296,296,266,296,296,266,3360,3366,3360,4970,17126,17126] #sensible selection, 0 SFOS one category
#bins = [296,296,266,290,290,266,3360,3366,3360,4970,17126,17126] #sensible selection, 0 SFOS one category

#bins = [296,296,290,296,296,266,3360,3366,3360,1892,2000,17126,17126] #sensible selection, 0 SFOS eem + emm splitted
#bins = [296,296,290,290,290,266,3360,3366,3360,1892,2000,17126,17126] #sensible selection, 0 SFOS eem + emm splitted
#bins = [296,296,266,296,296,266,3360,3366,3360,1892,2000,17126,17126] #sensible selection, 0 SFOS eem + emm splitted
#bins = [296,296,266,290,290,266,3360,3366,3360,1892,2000,17126,17126] #sensible selection, 0 SFOS eem + emm splitted

bins = [296,296,266,290,290,266,3360,3366,3360,14,26855,26855] #sensible selection, 0 SFOS one category
#bins = [296,296,266,290,290,266,3360,3366,3360,2,14,26855,26855] #sensible selection, 0 SFOS eem + emm splitted
#bins = [296,296,266,290,290,266,3360,3366,3360,3,14,26855,26855] #sensible selection, 0 SFOS eem + emm splitted
#bins = [296,296,266,290,290,266,3360,3366,3360,14,14,26855,26855] #sensible selection, 0 SFOS eem + emm splitted

makeDummyDataCard(path,bins)

#subprocess.call(["combine","-M","Significance dummycard.txt "])
os.system("combine -M Significance dummycard.txt  -t -1 --expectSignal=1")
"""
