import ROOT
import math


relunc=[0.25,0.30,1.0,0.5,1.0]
relunc=[0.25,0.25,0.5,0.5,0.5]

def printSelection(mybin,mytype="SFOS"):
    if "SFOS" in mytype:
        b = mybin
        d3 = int(b/8100)
        b = b%8100
        met = int(b/1620)
        b = b%1620
        pt3 = int(b/270)
        b = b%270
        mt = int(b/54)
        b = b%54
        lpt1 = int(b/27)        
        b = b%27
        lpt = int(b/9)
        b = b%9
        nj = int(b/3)
        b = b%3
        nb = int(b)

        snb  = "nb == 0       " if nb==1 else "nb == nsb == 0" if nb==2 else "              "
        snj  = "nj30 <= 1" if nj==1 else "nj30 == 0" if nj==2 else "         "
        spt  = "l-pt > 25,20,20" if lpt==0 else "l-pt > 25,25,25" if lpt==1 else "l-pt > 30,30,30" if lpt==2 else "               "
        spt1 = "l1-pt > 30" if lpt1==1 else "          "
        smt  = "MTmax >  60" if mt==1 else "MTmax >  90" if mt==2 else "MTmax > 120" if mt==3 else "MTmax > 150" if mt==4 else "           "
        spt3 = "Pt3l > 30" if pt3==1 else "Pt3l > 50" if pt3==2 else "Pt3l > 60" if pt3==3 else "Pt3l > 75" if pt3==4 else "Pt3l > 90" if pt3==5 else "         "
        smet  = "MET > 30" if met==1 else "MET > 45" if met==2 else "MET > 60" if met==3 else "MET > 75" if met==4 else "        "
        sd3  = "DPhi3lMET > 1.5" if d3==1 else "DPhi3lMET > 2.1" if d3==2 else "DPhi3lMET > 2.5" if d3==3 else "DPhi3lMET > 2.7" if d3==4 else "DPhi3lMET > 2.9" if d3==5 else "               "
        return spt+", "+spt1+", "+snb+", "+snj+", "+smt+", "+smet+", "+spt3+", "+sd3
    elif "SRSS" in mytype or "SS1J" in mytype:
        b = mybin
        ptll = int(b/21000)
        b = b%21000
        drlj = int(b/3000)
        b = b%3000
        mlj = int(b/600)
        b = b%600
        mll = int(b/120)
        b = b%120
        met = int(b/24)
        b = b%24
        mt = int(b/6)
        b = b%6
        lpt = int(b/3)
        b = b%3
        nb = int(b)

        snb  = "nb == 0       " if nb==1 else "nb == nsb == 0" if nb==2 else "              "
        spt  = "l-pt > 25,25" if lpt==0 else "l-pt > 30,30" if lpt==1 else "            "
        smt  = "MTmax >  90" if mt==1 else "MTmax > 120" if mt==2 else "MTmax > 150" if mt==3 else "           "
        smet  = "MET > 30" if met==1 else "MET > 45" if met==2 else "MET > 60" if met==3 else "MET > 75" if met==4 else "        "
        smll  = "Mll > 40" if mll==1 else "Mll > 60" if mll==2 else "Mll > 75" if mll==3 else "Mll > 90" if mll==4 else "Mll > 20"
        smlj  = "Mljmin <  75" if mlj==1 else "Mljmin < 100" if mlj==2 else "Mljmin < 125" if mlj==3 else "Mljmin < 150" if mlj==4 else "            "
        sdrlj  = "DRljmin < 1.2" if drlj==1 else "DRljmin < 1.5" if drlj==2 else "DRljmin < 1.8" if drlj==3 else "DRljmin < 2.0" if drlj==4 else "DRljmin < 2.1" if drlj==5 else "DRljmin < 2.5" if drlj==6 else "             "
        sptll  = "Ptll <  75" if ptll==1 else "Ptll < 100" if ptll==2 else "Ptll < 140" if ptll==3 else "          "

        return spt+", "+snb+", "+smt+", "+smet+", "+smll+", "+sdrlj+", "+smlj+", "+sptll
    else:
        print("Type "+mytype+" not recognized.")
        return ""



def getbkgsysterror(hbkg,relunc,mybin):
    SB = 0.
    if len(hbkg) != len(relunc):
        print("hbkg has different length than relunc",len(hbkg),len(relunc))
        return -3.
    if len(hbkg)==0:
        print("you need to provide signal and background histogram lists that are not empty")
        return -2.
    if mybin < 0:
        return -4.
    if mybin > hbkg[0].GetNbinsX():
        return -5.

    binnum = hbkg[0].FindBin(mybin)
    P =  hbkg[0].GetBinContent(binnum)
    L =  hbkg[1].GetBinContent(binnum)
    F =  hbkg[2].GetBinContent(binnum)
    Q =  hbkg[3].GetBinContent(binnum)
    G =  hbkg[4].GetBinContent(binnum)
    Pe =  hbkg[0].GetBinError(binnum)
    Le =  hbkg[1].GetBinError(binnum)
    Fe =  hbkg[2].GetBinError(binnum)
    Qe =  hbkg[3].GetBinError(binnum)
    Ge =  hbkg[4].GetBinError(binnum)

    Pe = math.sqrt(Pe*Pe+pow((relunc[0])*P,2))
    Le = math.sqrt(Le*Le+pow((relunc[1])*L,2))
    Fe = math.sqrt(Fe*Fe+pow((relunc[2])*F,2))
    Qe = math.sqrt(Qe*Qe+pow((relunc[3])*Q,2))
    Ge = math.sqrt(Ge*Ge+pow((relunc[4])*G,2))
    #print(P,L,F,Q,G,lowcut,upcut)

    Besq = pow(Pe,2)+pow(Le,2)+pow(Fe,2)+pow(Qe,2)+pow(Ge,2)
    return math.sqrt(Besq)

def GetErrHist(hbkg,relunc):
    h = hbkg[0].Clone("hBErr")
    for i in range(0,h.GetNbinsX()+1):
        h.SetBinError(i,0)
        h.SetBinContent(i,0)

    if len(hbkg) != len(relunc):
        print("hbkg has different length than relunc",len(hbkg),len(relunc))
        return h
    if len(hbkg)==0:
        print("you need to provide signal and background histogram lists that are not empty")
        return h

    for i in range(0,h.GetNbinsX()+1):
        P =  hbkg[0].GetBinContent(i)
        L =  hbkg[1].GetBinContent(i)
        F =  hbkg[2].GetBinContent(i)
        Q =  hbkg[3].GetBinContent(i)
        G =  hbkg[4].GetBinContent(i)
        Pe =  hbkg[0].GetBinError(i)
        Le =  hbkg[1].GetBinError(i)
        Fe =  hbkg[2].GetBinError(i)
        Qe =  hbkg[3].GetBinError(i)
        Ge =  hbkg[4].GetBinError(i)

        Pe = math.sqrt(Pe*Pe+pow((relunc[0])*P,2))
        Le = math.sqrt(Le*Le+pow((relunc[1])*L,2))
        Fe = math.sqrt(Fe*Fe+pow((relunc[2])*F,2))
        Qe = math.sqrt(Qe*Qe+pow((relunc[3])*Q,2))
        Ge = math.sqrt(Ge*Ge+pow((relunc[4])*G,2))
        h.SetBinContent(i,math.sqrt(pow(Pe,2)+pow(Le,2)+pow(Fe,2)+pow(Qe,2)+pow(Ge,2)))
    return h

def getyieldanderror(hS,mybin):
    if mybin < 0:
        return -4.,-4
    if mybin > hS.GetNbinsX():
        return -5.,-5
    binnum = hS.FindBin(mybin)
    sig =  hS.GetBinContent(binnum)
    err =  hS.GetBinError(binnum)
    return sig,err

def signif(sig,bkg,bkgerr):
    if bkg == 0:
        return -1;
    #print(sig,bkg,bkgerr,bkgerr*bkgerr,math.sqrt(bkg+pow(bkgerr,2)))
    if bkg+pow(bkgerr,2) <= 0.01:
        return -1
    return sig/math.sqrt(bkg+pow(bkgerr,2))

def GetSBplus(hsig,hbkg,relunc,mybin):
    if len(hbkg) != len(relunc):
        print("hbkg has different length than relunc",len(hbkg),len(relunc))
        return -3.,0.,0.,0.,0.
    if len(hsig)==0 or len(hbkg)==0:
        print("you need to provide signal and background histogram lists that are not empty")
        return -2.,0.,0.,0.,0.
    hS = hsig[0].Clone("hS")
    for i in range(1,len(hsig)):
        hS.Add(hsig[i])
    hB = hbkg[0].Clone("hB")
    for i in range(1,len(hbkg)):
        hB.Add(hbkg[i])
    if mybin < 0:
        return -4.,0.,0.,0.,0.
    if mybin > hS.GetNbinsX():
        return -5.,0.,0.,0.,0.
    s,serr = getyieldanderror(hS,mybin)
    if s < 0:
        return -1.5,0.,0.,0.,0.
    b,berr = getyieldanderror(hB,mybin)
    if b < 0:
        return -1.,0.,0.,0.,0.
    bsyst =  getbkgsysterror(hbkg,relunc,mybin)
    sb = signif(s,b,bsyst)
    pSB = 0
    if b>0:
        pSB = math.sqrt(2*((s+b)*math.log(1+s/b)-s))
    return sb,s,serr,b,berr,pSB

def GetSB(hsig,hbkg,hbkgunc,mybin):
    if mybin < 0:
        return -4.,0.,0.,0.,0.
    if mybin > hsig.GetNbinsX():
        return -5.,0.,0.,0.,0.
    binnum = hsig.FindBin(mybin)
    s    = hsig.GetBinContent(binnum)
    serr = hsig.GetBinError(binnum)
    b    = hbkg.GetBinContent(binnum)
    berr = hbkg.GetBinError(binnum)
    bsys = hbkgunc.GetBinContent(binnum)
    sb = signif(s,b,bsys)
    pSB = 0
    if b>0:
        pSB = math.sqrt(2*((s+b)*math.log(1+s/b)-s))
    return sb,s,serr,b,berr,pSB


def gettype(sel):
    if "SRSSeeMjjIn" in sel:
        return "SRSSeeMjjIn "
    if "SRSSeeMjjOut" in sel:
        return "SRSSeeMjjOut"
    if "SRSSemMjjIn" in sel:
        return "SRSSemMjjIn "
    if "SRSSemMjjOut" in sel:
        return "SRSSemMjjOut"
    if "SRSSmmMjjIn" in sel:
        return "SRSSmmMjjIn "
    if "SRSSmmMjjOut" in sel:
        return "SRSSmmMjjOut"
    if "SRSS1Jee" in sel:
        return "SS1Jee      "
    if "SRSS1Jem" in sel:
        return "SS1Jem      "
    if "SRSS1Jmm" in sel:
        return "SS1Jmm      "
    if "SR0SFOSeem" in sel:
        return "0SFOSeem    "
    if "SR0SFOSemm" in sel:
        return "0SFOSemm    "
    if "SR0SFOS" in sel:
        return "0SFOS       "
    if "SR1SFOS" in sel:
        return "1SFOS       "
    if "SR2SFOS" in sel:
        return "2SFOS       "
    return ""
def PrintSB(sb,s,serr,b,berr,mybin,sel):
    if b<=0:
        return True
    SBs = "{0:.3g}".format(sb)
    ss    = "{0:.2g}".format(s)
    serrs = "{0:.2g}".format(serr)
    bs    = "{0:.2g}".format(b)
    berrs = "{0:.2g}".format(berr)
    sobs  = "{0:.3g}".format(s/b)
    sosbs = "{0:.3g}".format(s/math.sqrt(b))
    pSB = 0
    if b>0:
        pSB = math.sqrt(2*((s+b)*math.log(1+s/b)-s))
    spsbs  = "{0:.3g}".format(pSB)
    mytype = gettype(sel)
    if mytype=="":
        print("type not recognized from histogram",sel)
        return False
    selection = printSelection(mybin,mytype)
    if selection != "":
        print ("for "+mytype+" and "+selection+" have S/sqrt(B+dB^2) = "+SBs.ljust(5)+", Philip's Sig = "+spsbs.ljust(5)+", S/sqrt(B) = "+sosbs.ljust(5)+", S/B = "+sobs.ljust(6)+", S="+ss.rjust(4)+"+/-"+serrs.ljust(4)+", B="+bs.rjust(5)+"+/-"+berrs.ljust(4)+" ("+str(int(mybin))+")")
    return True

def Printyields(S,SErr,B,BErr,P,PErr,L,LErr,F,FErr,Q,QErr,G,GErr,mybin,sel):

    if B<=0:
        return True
    ss    = "{0:.3g}".format(S)
    serrs = "{0:.3g}".format(SErr)
    bs    = "{0:.3g}".format(B)
    berrs = "{0:.3g}".format(BErr)
    ps    = "{0:.3g}".format(P)
    perrs = "{0:.3g}".format(PErr)
    ls    = "{0:.3g}".format(L)
    lerrs = "{0:.3g}".format(LErr)
    fs    = "{0:.3g}".format(F)
    ferrs = "{0:.3g}".format(FErr)
    qs    = "{0:.3g}".format(Q)
    qerrs = "{0:.3g}".format(QErr)
    gs    = "{0:.3g}".format(G)
    gerrs = "{0:.3g}".format(GErr)
    mytype = gettype(sel)
    if mytype=="":
        print("type not recognized from histogram",sel)
        return False
    selection = printSelection(mybin,mytype)
    if selection != "":
        print ("for "+mytype+" and "+selection+": S="+ss.rjust(4)+"+/-"+serrs.ljust(4)+", B="+bs.rjust(5)+"+/-"+berrs.ljust(4)+" (LL="+ls.rjust(5)+"+/-"+lerrs.ljust(4)+", Pr="+ps.rjust(4)+"+/-"+perrs.ljust(4)+", F="+fs.rjust(4)+"+/-"+ferrs.ljust(4)+", Qf="+qs.rjust(4)+"+/-"+qerrs.ljust(4)+", gf="+gs.rjust(4)+"+/-"+gerrs.ljust(4)+")")
    return True

def GetAllyields(files,selection,mybin,Smin=-1,SBmin=-1,Serrmax=0.1):

    histname = selection
    hS = files[0].Get(histname)
    hP = files[1].Get(histname)
    hL = files[2].Get(histname)
    hF = files[3].Get(histname)
    hQ = files[4].Get(histname)
    hG = files[5].Get(histname)
    
    hsig = [hS]
    hbkg = [hP,hL,hF,hQ,hG]
    #relunc=[0.25,0.3,1.,0.5,1.]
    P, Perr = getyieldanderror(hP,mybin)
    L, Lerr = getyieldanderror(hL,mybin)
    F, Ferr = getyieldanderror(hF,mybin)
    Q, Qerr = getyieldanderror(hQ,mybin)
    G, Gerr = getyieldanderror(hG,mybin)
    
    SB,S,SErr,B,BErr,pSB = GetSBplus(hsig,hbkg,relunc,mybin)
    if S >= Smin and SB >= SBmin and SErr < Serrmax*S:
        Printyields(S,SErr,B,BErr,P,Perr,L,Lerr,F,Ferr,Q,Qerr,G,Gerr,mybin,selection)
    return S,B,SB,pSB


def GetSignifForSelection(files,selection,mybin,Smin=-1,SBmin=-1,Serrmax=0.1):

    histname = selection
    hS = files[0].Get(histname)
    hP = files[1].Get(histname)
    hL = files[2].Get(histname)
    hF = files[3].Get(histname)
    hQ = files[4].Get(histname)
    hG = files[5].Get(histname)
    
    hsig = [hS]
    hbkg = [hP,hL,hF,hQ,hG]
    #relunc=[0.10,0.10,0.2,0.3,0.5]
    ##relunc=[0.25,0.3,1.,0.5,1.]
    ##relunc=[0.,0.,0.,0.,0.]
    
    SB,S,SErr,B,BErr,pSB = GetSBplus(hsig,hbkg,relunc,mybin)
    if S >= Smin and SB >= SBmin and SErr < Serrmax*S:
        #print(SB,S/B,S,SErr,B,BErr,selection,mybin)
        PrintSB(SB,S,SErr,B,BErr,mybin,selection)
    return S,B,SB,pSB



def GetSignifForSelectionScan(files,selection,Smin=-1,SBmin=-1,Serrmax=0.1):

    print("xx")
    hS = files[0].Get(selection)
    print(hS)
    hP = files[1].Get(selection)
    hL = files[2].Get(selection)
    hF = files[3].Get(selection)
    hQ = files[4].Get(selection)
    hG = files[5].Get(selection)
    
    hsig = [hS]
    hbkg = [hP,hL,hF,hQ,hG]
    #relunc=[0.25,0.3,1.,0.5,1.]
    ##relunc=[0.,0.,0.,0.,0.]
    
    hB = hP.Clone("hB")
    hB.Add(hL)
    hB.Add(hF)
    hB.Add(hQ)
    hB.Add(hG)
    hBerr = GetErrHist(hbkg,relunc)

    for i in range(0,hS.GetNbinsX()):
        SB,S,SErr,B,BErr,pSB = GetSB(hS,hB,hBerr,i)
        if S >= Smin and SB >= SBmin and SErr/S < Serrmax:
            PrintSB(SB,S,SErr,B,BErr,i,selection)
    return True

def Sort_ListedTuple(tup,sortwhat="S/sqrt(B+dB^2)"):  
  
    # reverse = None (Sorts in Ascending order)  
    # key is set to sort using second element of  
    # sublist lambda has been used
    if sortwhat == "S":
        tup.sort(key = lambda x: x[3], reverse = True)
    elif sortwhat == "B":
        tup.sort(key = lambda x: x[4])
    elif sortwhat == "SoB" or sortwhat == "SoverB":
        tup.sort(key = lambda x: x[1], reverse = True)
    elif sortwhat == "pSB":
        tup.sort(key = lambda x: x[2], reverse = True)
    else:
        tup.sort(key = lambda x: x[0], reverse = True)

    return tup  

def GetSignifForSelectionScanSorted(files,selection,Smin=-1,SBmin=-1,Serrmax=0.1,sortwhat="S/sqrt(B+dB^2)"):
    print("********* Below is optimized yields for "+selection+" *********")
    hS = files[0].Get(selection)
    hP = files[1].Get(selection)
    hL = files[2].Get(selection)
    hF = files[3].Get(selection)
    hQ = files[4].Get(selection)
    hG = files[5].Get(selection)
    
    hsig = [hS]
    hbkg = [hP,hL,hF,hQ,hG]
    #relunc=[0.25,0.3,1.,0.5,1.]
    ##relunc=[0.,0.,0.,0.,0.]

    hB = hP.Clone("hB")
    hB.Add(hL)
    hB.Add(hF)
    hB.Add(hQ)
    hB.Add(hG)
    hBerr = GetErrHist(hbkg,relunc)

    listSB, listSoB, listS, listB, listSErr, listBErr,listi, listpSB = [],[],[],[],[],[],[],[]
    for i in range(0,hS.GetNbinsX()):
    #for i in range(0,100):
        SB,S,SErr,B,BErr,pSB = GetSB(hS,hB,hBerr,i)
        SoB = -1
        if B>0: SoB = S/B
        listSB.append(SB)
        listSoB.append(SoB)
        listS.append(S)
        listB.append(B)
        listSErr.append(SErr)
        listBErr.append(BErr)
        listi.append(i)
        listpSB.append(pSB)
        #PrintSB(SB,S,SErr,B,BErr,i,selection)

    allcombined = zip(listSB,listSoB,listpSB,listS,listB,listSErr,listBErr,listi)
    #print(allcombined)
    allcombined = Sort_ListedTuple(allcombined,sortwhat)
    #print(len(allcombined))
    for i in range(0,len(allcombined)):
        mySB = allcombined[i][0]
        if sortwhat == "S":
            mySB = allcombined[i][3]
        elif sortwhat == "B":
            mySB = allcombined[i][4]
        elif sortwhat == "SoB" or sortwhat == "SoverB":
            mySB = allcombined[i][1]
        elif sortwhat == "pSB":
            mySB = allcombined[i][2]
        #if i in range(0,20):
        #    PrintSB(allcombined[i][0],allcombined[i][3],allcombined[i][5],allcombined[i][4],allcombined[i][6],allcombined[i][7],selection)
        #if S >= Smin and SB >= SBmin and SErr/S < Serrmax:
        if allcombined[i][3] >= Smin and mySB >= SBmin and allcombined[i][5] < Serrmax*allcombined[i][3]:
            #PrintSB(SB,S,SErr,B,BErr,i,selection)
            #continue
            PrintSB(allcombined[i][0],allcombined[i][3],allcombined[i][5],allcombined[i][4],allcombined[i][6],allcombined[i][7],selection)
    #print(allcombined)

def GetNumBins(files,histname):
    h = files[0].Get(histname)
    return h.GetNbinsX()



path = "hists/combineyearsLoose_v5.3.0/reOptimize/"
path = "hists/Loose2016_v5.3.0/reOptimize/"
fS = ROOT.TFile.Open(path+"signal.root")
#fS = ROOT.TFile.Open(path+"signal_private.root")
fP = ROOT.TFile.Open(path+"prompt.root")
fL = ROOT.TFile.Open(path+"lostlep.root")
fF = ROOT.TFile.Open(path+"fakes.root")
#fF = ROOT.TFile.Open(path+"ddfakes.root")
fQ = ROOT.TFile.Open(path+"qflip.root")
fG = ROOT.TFile.Open(path+"photon.root")
files = [fS,fP,fL,fF,fQ,fG]



print("************* ee-Mjj in *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",               0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",               1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",               2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",120*1+6*0+24*3+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",120*1+6*0+24*3+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",      6*1+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjIn__CutsSS",               0,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjIn__CutsSS",               1,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjIn__CutsSS",               2,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjIn__CutsSS",120*1+6*0+24*3+1,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjIn__CutsSS",120*1+6*0+24*3+2,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjIn__CutsSS",      6*1+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSSeeMjjIn__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSSeeMjjIn__CutsSS",1.5,SB,0.3)
print("************* em-Mjj in *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",               0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",               1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",               2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",120*1+6*1+24*3+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",120*1+6*1+24*3+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",      6*1+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjIn__CutsSS",               0,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjIn__CutsSS",               1,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjIn__CutsSS",               2,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjIn__CutsSS",120*1+6*1+24*3+1,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjIn__CutsSS",120*1+6*1+24*3+2,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjIn__CutsSS",      6*1+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSSemMjjIn__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSSemMjjIn__CutsSS",3.7,SB,0.15)
print("************* mm-Mjj in *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",               0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",               1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",               2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",120*1+6*0+24*0+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",120*1+6*0+24*0+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*0+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*0+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjIn__CutsSS",               0,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjIn__CutsSS",               1,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjIn__CutsSS",               2,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjIn__CutsSS",120*1+6*0+24*0+1,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjIn__CutsSS",120*1+6*0+24*0+2,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjIn__CutsSS",      6*0+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSSmmMjjIn__CutsSS",      6*0+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSSmmMjjIn__CutsSS",4,SB,0.15)


print("************* ee-Mjj out *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",               0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",               1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",               2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",120*1+6*0+24*3+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",120*1+6*0+24*3+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",      6*1+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjOut__CutsSS",               0,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjOut__CutsSS",               1,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjOut__CutsSS",               2,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjOut__CutsSS",120*1+6*0+24*3+1,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjOut__CutsSS",120*1+6*0+24*3+2,-1,-1,0.99)
#GetAllyields(files,"SRSSeeMjjOut__CutsSS",      6*1+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSSeeMjjOut__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSSeeMjjOut__CutsSS",1.5,SB,0.3)
print("************* em-Mjj out *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",               0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",               1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",               2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",120*1+6*1+24*3+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",120*1+6*1+24*3+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",      6*1+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjOut__CutsSS",               0,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjOut__CutsSS",               1,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjOut__CutsSS",               2,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjOut__CutsSS",120*1+6*1+24*3+1,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjOut__CutsSS",120*1+6*1+24*3+2,-1,-1,0.99)
#GetAllyields(files,"SRSSemMjjOut__CutsSS",      6*1+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSSemMjjOut__CutsSS",      6*1+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSSemMjjOut__CutsSS",3.7,SB,0.15)
print("************* mm-Mjj out *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",               0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",               1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",               2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",120*1+6*0+24*3+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",120*1+6*0+24*3+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",      6*0+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",      6*0+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjOut__CutsSS",               0,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjOut__CutsSS",               1,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjOut__CutsSS",               2,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjOut__CutsSS",120*1+6*0+24*3+1,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjOut__CutsSS",120*1+6*0+24*3+2,-1,-1,0.99)
#GetAllyields(files,"SRSSmmMjjOut__CutsSS",      6*0+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSSmmMjjOut__CutsSS",      6*0+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSSmmMjjOut__CutsSS",4,SB,0.15)


print("************* ee-1j *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jee1JPre__CutsSS",                0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jee1JPre__CutsSS",                1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jee1JPre__CutsSS",                2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jee1JPre__CutsSS",3000*2+6*1+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jee1JPre__CutsSS",3000*2+6*1+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jee1JPre__CutsSS",                0,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jee1JPre__CutsSS",                1,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jee1JPre__CutsSS",                2,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jee1JPre__CutsSS",3000*2+6*1+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSS1Jee1JPre__CutsSS",3000*2+6*1+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSS1Jee1JPre__CutsSS",2.5,SB,0.3)
print("************* em-1j *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jem1JPre__CutsSS",                0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jem1JPre__CutsSS",                1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jem1JPre__CutsSS",                2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jem1JPre__CutsSS",3000*2+6*1+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jem1JPre__CutsSS",3000*2+6*1+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jem1JPre__CutsSS",                0,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jem1JPre__CutsSS",                1,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jem1JPre__CutsSS",                2,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jem1JPre__CutsSS",3000*2+6*1+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSS1Jem1JPre__CutsSS",3000*2+6*1+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSS1Jem1JPre__CutsSS",4.,SB,0.3)
print("************* mm-1j *************")
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jmm1JPre__CutsSS",                0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jmm1JPre__CutsSS",                1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jmm1JPre__CutsSS",                2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jmm1JPre__CutsSS",3000*2+6*1+24*2+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSS1Jmm1JPre__CutsSS",3000*2+6*1+24*2+2,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jmm1JPre__CutsSS",                0,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jmm1JPre__CutsSS",                1,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jmm1JPre__CutsSS",                2,-1,-1,0.99)
#GetAllyields(files,"SRSS1Jmm1JPre__CutsSS",3000*2+6*1+24*2+1,-1,-1,0.99)
GetAllyields(files,"SRSS1Jmm1JPre__CutsSS",3000*2+6*1+24*2+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SRSS1Jmm1JPre__CutsSS",4.,SB,0.3)


print("************* 0 SFOS *************")
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",                       0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",                       1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",                       2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",8100*3+1620*1+54*2+3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",8100*3+1620*1+54*2+3*1+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",                   3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",                   3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOS__CutsSFOS",                       0,-1,-1,0.99)
#GetAllyields(files,"SR0SFOS__CutsSFOS",                       1,-1,-1,0.99)
#GetAllyields(files,"SR0SFOS__CutsSFOS",                       2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOS__CutsSFOS",8100*3+1620*1+54*2+3*1+1,-1,-1,0.99)
#GetAllyields(files,"SR0SFOS__CutsSFOS",8100*3+1620*1+54*2+3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOS__CutsSFOS",                   3*1+1,-1,-1,0.99)
GetAllyields(files,"SR0SFOS__CutsSFOS",                   3*1+2,-1,-1,0.99)
GetAllyields(files,"SR0SFOS__CutsSFOS",                 9+3*1+2,-1,-1,0.99)
GetAllyields(files,"SR0SFOS__CutsSFOS",                 9+3*0+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SR0SFOS__CutsSFOS",4.,SB,0.15)
print("************* 0 SFOS eem *************")
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",                       0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",                       1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",                       2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",8100*3+1620*1+54*2+3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",8100*3+1620*1+54*2+3*1+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",                   3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",                   3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSeem__CutsSFOS",                       0,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSeem__CutsSFOS",                       1,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSeem__CutsSFOS",                       2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSeem__CutsSFOS",8100*3+1620*1+54*2+3*1+1,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSeem__CutsSFOS",8100*3+1620*1+54*2+3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSeem__CutsSFOS",                   3*1+1,-1,-1,0.99)
GetAllyields(files,"SR0SFOSeem__CutsSFOS",                   3*1+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SR0SFOSeem__CutsSFOS",4.,SB,0.15)
print("************* 0 SFOS emm *************")
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",                       0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",                       1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",                       2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",8100*3+1620*1+54*2+3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",8100*3+1620*1+54*2+3*1+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",                   3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",                   3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSemm__CutsSFOS",                       0,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSemm__CutsSFOS",                       1,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSemm__CutsSFOS",                       2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSemm__CutsSFOS",8100*3+1620*1+54*2+3*1+1,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSemm__CutsSFOS",8100*3+1620*1+54*2+3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR0SFOSemm__CutsSFOS",                   3*1+1,-1,-1,0.99)
GetAllyields(files,"SR0SFOSemm__CutsSFOS",                   3*1+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SR0SFOSemm__CutsSFOS",4.,SB,0.15)

print("************* 1 SFOS *************")
S,B,SB,pSB=GetSignifForSelection(files,"SR1SFOS__CutsSFOS",                             0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR1SFOS__CutsSFOS",                             1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR1SFOS__CutsSFOS",                             2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR1SFOS__CutsSFOS",8100*3+1620*2+270*3+54*2+3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR1SFOS__CutsSFOS",8100*3+1620*2+270*3+54*2+3*1+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR1SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR1SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR1SFOS__CutsSFOS",                             0,-1,-1,0.99)
#GetAllyields(files,"SR1SFOS__CutsSFOS",                             1,-1,-1,0.99)
#GetAllyields(files,"SR1SFOS__CutsSFOS",                             2,-1,-1,0.99)
#GetAllyields(files,"SR1SFOS__CutsSFOS",8100*3+1620*2+270*3+54*2+3*1+1,-1,-1,0.99)
#GetAllyields(files,"SR1SFOS__CutsSFOS",8100*3+1620*2+270*3+54*2+3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR1SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+1,-1,-1,0.99)
GetAllyields(files,"SR1SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+2,-1,-1,0.99)
GetAllyields(files,"SR1SFOS__CutsSFOS",29813,-1,-1,0.99)
GetAllyields(files,"SR1SFOS__CutsSFOS",29816,-1,-1,0.99)
GetAllyields(files,"SR1SFOS__CutsSFOS",25226,-1,-1,0.99)
GetAllyields(files,"SR1SFOS__CutsSFOS",21986,-1,-1,0.99)
GetAllyields(files,"SR1SFOS__CutsSFOS",21716,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SR1SFOS__CutsSFOS",4,SB,0.1)
print("************* 2 SFOS *************")
S,B,SB,pSB=GetSignifForSelection(files,"SR2SFOS__CutsSFOS",                             0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR2SFOS__CutsSFOS",                             1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR2SFOS__CutsSFOS",                             2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR2SFOS__CutsSFOS",8100*3+1620*2+270*3+54*2+3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR2SFOS__CutsSFOS",8100*3+1620*2+270*3+54*2+3*1+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR2SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+1,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR2SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR2SFOS__CutsSFOS",                             0,-1,-1,0.99)
#GetAllyields(files,"SR2SFOS__CutsSFOS",                             1,-1,-1,0.99)
#GetAllyields(files,"SR2SFOS__CutsSFOS",                             2,-1,-1,0.99)
#GetAllyields(files,"SR2SFOS__CutsSFOS",8100*3+1620*3+270*3+54*0+3*1+1,-1,-1,0.99)
#GetAllyields(files,"SR2SFOS__CutsSFOS",8100*3+1620*3+270*3+54*0+3*1+2,-1,-1,0.99)
#GetAllyields(files,"SR2SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+1,-1,-1,0.99)
GetAllyields(files,"SR2SFOS__CutsSFOS",8100*3+1620*3+270*2+54*2+3*1+2,-1,-1,0.99)
#GetSignifForSelectionScanSorted(files,"SR2SFOS__CutsSFOS",4,SB,0.1)
"""
print("************* extras 1 *************")
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOS__CutsSFOS",                   14,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",                   14,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",                   14,-1,-1,0.99)
GetAllyields(files,"SR0SFOS__CutsSFOS",                   14,-1,-1,0.99)
GetAllyields(files,"SR0SFOSeem__CutsSFOS",                   14,-1,-1,0.99)
GetAllyields(files,"SR0SFOSemm__CutsSFOS",                   14,-1,-1,0.99)
"""

"""
print("************* extras 2 *************")
GetSignifForSelection(files,"SRSSeeMjjIn__CutsSS",    bins[ 0],-1,-1,0.99)
GetSignifForSelection(files,"SRSSemMjjIn__CutsSS",    bins[ 1],-1,-1,0.99)
GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",    bins[ 2],-1,-1,0.99)
GetSignifForSelection(files,"SRSSeeMjjOut__CutsSS",   bins[ 3],-1,-1,0.99)
GetSignifForSelection(files,"SRSSemMjjOut__CutsSS",   bins[ 4],-1,-1,0.99)
GetSignifForSelection(files,"SRSSmmMjjOut__CutsSS",   bins[ 5],-1,-1,0.99)
GetSignifForSelection(files,"SRSS1Jee1JPre__CutsSS",  bins[ 6],-1,-1,0.99)
GetSignifForSelection(files,"SRSS1Jem1JPre__CutsSS",  bins[ 7],-1,-1,0.99)
GetSignifForSelection(files,"SRSS1Jmm1JPre__CutsSS",  bins[ 8],-1,-1,0.99)
GetSignifForSelection(files,"SR0SFOS__CutsSFOS",      bins[ 9],-1,-1,0.99)
GetSignifForSelection(files,"SR0SFOSeem__CutsSFOS",   bins[10],-1,-1,0.99)
GetSignifForSelection(files,"SR0SFOSemm__CutsSFOS",   bins[11],-1,-1,0.99)
GetSignifForSelection(files,"SR1SFOS__CutsSFOS",      bins[12],-1,-1,0.99)
GetSignifForSelection(files,"SR2SFOS__CutsSFOS",      bins[13],-1,-1,0.99)
"""
"""
GetAllyields(files,"SRSSeeMjjIn__CutsSS",    bins[ 0],-1,-1,0.99)
GetAllyields(files,"SRSSemMjjIn__CutsSS",    bins[ 1],-1,-1,0.99)
GetAllyields(files,"SRSSmmMjjIn__CutsSS",    bins[ 2],-1,-1,0.99)
GetAllyields(files,"SRSSeeMjjOut__CutsSS",   bins[ 3],-1,-1,0.99)
GetAllyields(files,"SRSSemMjjOut__CutsSS",   bins[ 4],-1,-1,0.99)
GetAllyields(files,"SRSSmmMjjOut__CutsSS",   bins[ 5],-1,-1,0.99)
GetAllyields(files,"SRSS1Jee1JPre__CutsSS",  bins[ 6],-1,-1,0.99)
GetAllyields(files,"SRSS1Jem1JPre__CutsSS",  bins[ 7],-1,-1,0.99)
GetAllyields(files,"SRSS1Jmm1JPre__CutsSS",  bins[ 8],-1,-1,0.99)
GetAllyields(files,"SR0SFOS__CutsSFOS",      bins[ 9],-1,-1,0.99)
GetAllyields(files,"SR0SFOSeem__CutsSFOS",   bins[10],-1,-1,0.99)
GetAllyields(files,"SR0SFOSemm__CutsSFOS",   bins[11],-1,-1,0.99)
GetAllyields(files,"SR1SFOS__CutsSFOS",      bins[12],-1,-1,0.99)
GetAllyields(files,"SR2SFOS__CutsSFOS",      bins[13],-1,-1,0.99)
"""
"""
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*0+24*2+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*0+24*1+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*0+24*0+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",    3+6*0+24*0+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",3*120+6*0+24*0+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*1+24*0+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*0+24*3+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*1+24*2+2,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"SRSSmmMjjIn__CutsSS",      6*2+24*2+2,-1,-1,0.99)
"""
print("Path",path)
