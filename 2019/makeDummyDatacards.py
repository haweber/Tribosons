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

