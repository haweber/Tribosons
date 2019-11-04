import ROOT
import math


def printSelection(mybin,mytype="0SFOS"):
    if mytype=="0SFOS":
        b = mybin
        nb0 = int(b/108000)
        b = b%108000
        ssid0 = int(b/10800)
        b = b%10800
        fid0 = int(b/3600)
        b = b%3600
        d30 = int(b/720)
        b = b%720
        met0 = int(b/144)
        b = b%144
        pt30 = int(b/36)
        b = b%36
        mt0 = int(b/9)
        b = b%9
        nj0 = int(b/3)
        b = b%3
        lpt0 = int(b)
        #if mt0!=2:
        #    return ""
        #if dl0!=0:
        #    return ""
        #if nj0<2:
        #    return ""
        #if(fid0!=0): return ""
        #if(ssid0!=0): return ""
        spt0  = "l-pt > 25,20,20" if lpt0==0 else "l-pt > 25,25,25" if lpt0==1 else "l-pt > 30,30,30" if lpt0==2 else "               "
        snb0  = "nb == 0" if nb0==1 else "       "
        snj0  = "nj30 <= 1" if nj0==1 else "nj30 == 0" if nj0==2 else "         "
        smt0  = "MTmax >  60" if mt0==1 else "MTmax >  90" if mt0==2 else "MTmax > 120" if mt0==3 else "           "
        spt30 = "Pt3l > 30" if pt30==1 else "Pt3l > 50" if pt30==2 else "Pt3l > 60" if pt30==3 else "         "
        smet0 = "MET > 30" if met0==1 else "MET > 45" if met0==2 else "MET > 60" if met0==3 else "MET > 75" if met0==4 else "        "
        sd30  = "DPhi3lMET > 1.5" if d30==1 else "DPhi3lMET > 2.1" if d30==2 else "DPhi3lMET > 2.5" if d30==3 else "DPhi3lMET > 2.7" if d30==4 else "               "
        sfid0 = "eemu " if fid0==1 else "emumu" if fid0==2 else "     "
        #sssid0 = "e SS-tight     " if ssid0==1 else "mu SS-tight    " if ssid0==2 else "e,mu SS-tight  " if ssid0==3 else "sub-e SS-tight " if ssid0==4 else "sub-mu SS-tight" if ssid0==5 else "sub-l SS-tight " if ssid0==6 else "               "
        sssid0 = "e SS-tight     " if ssid0==1 else "mu SS-tight    " if ssid0==2 else "e,mu SS-tight  " if ssid0==3 else "sub-e SS-tight " if ssid0==4 else "sub-mu SS-tight" if ssid0==5 else "sub-l SS-tight " if ssid0==6 else "e Ph-tight     " if ssid0==7 else "mu Ph-tight    " if ssid0==8 else "e,mu Ph-tight  " if ssid0==9 else "               "

        #return spt0+("" if spt0=="" else ", ")+snb0+("" if snb0=="" else ", ")+snj0+("" if snj0=="" else ", ")+smt0+("" if smt0=="" else ", ")+spt30+("" if spt30=="" else ", ")+smet0+("" if smet0=="" else ", ")+sd30+("" if sd30=="" else ", ")+sfid0+("" if sfid0=="" else ", ")+sssid0
        return spt0+", "+snb0+", "+snj0+", "+smt0+", "+spt30+", "+smet0+", "+sd30+", "+sfid0+", "+sssid0
    elif mytype=="1SFOS" or mytype=="2SFOS":
        b = mybin
        nb1 = int(b/777600)
        b = b%777600
        ssid1 = int(b/77760)
        b = b%77760
        fid1 = int(b/25920)
        b = b%25920
        d31 = int(b/4320)
        b = b%4320
        pt31 = int(b/720)
        b = b%720
        met1 = int(b/144)
        b = b%144
        mt1 = int(b/36)
        b = b%36
        nj1 = int(b/12)
        b = b%12
        pt11 = int(b/3)
        b = b%3
        pt1 = int(b)
        spt1  = "l-pt > 25,20,20" if pt1==0 else "l-pt > 25,25,25" if pt1==1 else "l-pt > 30,30,30" if pt1==2 else "                "
        spt11 = "l1-pt > 30" if pt11==1 else "l1-pt > 40" if pt11==2 else "l1-pt > 50" if pt11==3 else "          "
        snb1  = "nb == 0" if nb1==1 else "       "
        snj1  = "nj30 <= 1" if nj1==1 else "nj30 == 0" if nj1==2 else "         "
        smt1  = "MT3rd >  90" if mt1==1 else "MT3rd > 120" if mt1==2 else "MT3rd > 150" if mt1==3 else "           "
        if mytype=="2SFOS":
            smt1  = "MTmax >  90" if mt1==1 else "MTmax > 120" if mt1==2 else "MTmax > 150" if mt1==3 else "           "
        smet1  = "MET > 45" if met1==1 else "MET > 60" if met1==2 else "MET > 75" if met1==3 else "MET > 30"
        spt31 = "Pt3l >  30" if pt31==1 else "Pt3l >  50" if pt31==2 else "Pt3l >  60" if pt31==3 else "Pt3l >  75" if pt31==4 else "Pt3l > 100" if pt31==5 else "          "
        sd31  = "DPhi3lMET > 1.5" if d31==1 else "DPhi3lMET > 2.1" if d31==2 else "DPhi3lMET > 2.5" if d31==3 else "DPhi3lMET > 2.7" if d31==4 else "DPhi3lMET > 2.9" if d31==5 else "               "
        sfid1 = "eemu " if fid1==1 else "emumu" if fid1==2 else "     "
        sssid1 = "e SS-tight     " if ssid1==1 else "mu SS-tight    " if ssid1==2 else "e,mu SS-tight  " if ssid1==3 else "sub-e SS-tight " if ssid1==4 else "sub-mu SS-tight" if ssid1==5 else "sub-l SS-tight " if ssid1==6 else "e Ph-tight     " if ssid1==7 else "mu Ph-tight    " if ssid1==8 else "e,mu Ph-tight  " if ssid1==9 else "               "

        #return spt1+("" if spt1=="" else ", ")+spt11+("" if spt11=="" else ", ")+snb1+("" if snb1=="" else ", ")+snj1+("" if snj1=="" else ", ")+smt1+("" if smt1=="" else ", ")+smet1+("" if smet1=="" else ", ")+spt31+("" if spt31=="" else ", ")+sd31+("" if sd31=="" else ", ")+sfid1+("" if sfid1=="" else ", ")+sssid1;
        return spt1+", "+spt11+", "+snb1+", "+snj1+", "+smt1+", "+smet1+", "+spt31+", "+sd31+", "+sfid1+", "+sssid1;
    elif "SRSS" in mytype:
        b = mybin
        nbs = int(b/54000)
        b = b%54000
        ssids = int(b/13500)
        b = b%13500
        ptlls = int(b/4500)
        b = b%4500
        drljs = int(b/900)
        b = b%900
        mljs = int(b/180)
        b = b%180
        mlls = int(b/36)
        b = b%36
        mets = int(b/9)
        b = b%9
        mts = int(b/3)
        b = b%3
        pts = int(b)
        #if (mljs!=0 and drljs!=0) or (mljs!=0 and ptlls!=0) or (mljs!=0 and dls!=0) or (drljs!=0 and ptlls!=0) or (drljs!=0 and dls!=0) or (ptlls!=0 and dls!=0):
        #    return ""
        spts = "l1-pt > 30" if pts==1 else "l1-pt > 40" if pts==2 else "          "
        snbs  = "nb == 0" if nbs==1 else "       "
        smts  = "MTmax >  90" if mts==1 else "MTmax > 120" if mts==2 else "           "
        smets  = "MET > 45" if mets==1 else "MET > 60" if mets==2 else "MET > 75" if mets==3 else "        "
        smlls  = "Mll > 40" if mlls==1 else "Mll > 60" if mlls==2 else "Mll > 75" if mlls==3 else "Mll > 90" if mlls==4 else "Mll > 20"
        smlj1  = "Mljmin <  75" if mljs==1 else "Mljmin < 100" if mljs==2 else "Mljmin < 125" if mljs==3 else "Mljmin < 150" if mljs==4 else "            "
        sdrlj1  = "DRljmin < 1.5" if drljs==1 else "DRljmin < 1.8" if drljs==2 else "DRljmin < 2.1" if drljs==3 else "DRljmin < 2.5" if drljs==4 else "             "
        sptlls  = "Ptll <  75" if ptlls==1 else "Ptll < 100" if ptlls==2 else "          "
        sssids = "e Ph-tight     " if ssids==1 else "mu Ph-tight    " if ssids==2 else "e,mu Ph-tight  " if ssids==3 else "               "
        #return spts+("" if spts=="" else ", ")+snbs+("" if snbs=="" else ", ")+smts+("" if smts=="" else ", ")+smets+("" if smets=="" else ", ")+smlls+("" if smlls=="" else ", ")+smlj1+("" if smlj1=="" else ", ")+sdrlj1+("" if sdrlj1=="" else ", ")+sptlls+("" if sptlls=="" else ", ")+sssids;
        return spts+", "+snbs+", "+smts+", "+smets+", "+smlls+", "+smlj1+", "+sdrlj1+", "+sptlls+", "+sssids;
    elif "SS1J" in mytype:
        b = mybin
        nbs = int(b/54000)
        b = b%54000
        ssids = int(b/13500)
        b = b%13500
        ptlls = int(b/4500)
        b = b%4500
        mljs = int(b/900)
        b = b%900
        drljs = int(b/180)
        b = b%180
        mts = int(b/45)
        b = b%45
        mlls = int(b/15)
        b = b%15
        mets = int(b/3)
        b = b%3
        pts = int(b)
        spts  = "l-pt > 25,25" if pts==1 else "l-pt > 30,30" if pts==2 else "l-pt > 40,40" if pts==3 else "            "
        snbs  = "nb == 0" if nbs==1 else "       "
        smets  = "MET > 40" if mets==1 else "MET > 60" if mets==2 else "MET > 75" if mets==3 else "MET > 90" if mets==4 else "        "
        smlls  = "Mll > 40" if mlls==1 else "Mll > 60" if mlls==2 else "Mll > 20"
        smts  = "MTmax >  90" if mts==1 else "MTmax > 120" if mts==2 else "MTmax > 150" if mts==3 else "           "
        sdrljs  = "dRljmin < 1.2" if drljs==1 else "dRljmin < 1.5" if drljs==2 else "dRljmin < 1.8" if drljs==3 else "dRljmin < 2.1" if drljs==4 else "             "
        smljs  = "Mljmin <  75" if mljs==1 else "Mljmin < 100" if mljs==2 else "Mljmin < 125" if mljs==3 else "Mljmin < 150" if mljs==4  else "            "
        sptlls  = "Ptll < 100" if ptlls==1 else "Ptll < 150" if ptlls==2 else "          "
        sssids = "e Ph-tight     " if ssids==1 else "mu Ph-tight    " if ssids==2 else "e,mu Ph-tight  " if ssids==3 else "               "
        #return spts+("" if spts=="" else ", ")+snbs+("" if snbs=="" else ", ")+smts+("" if smts=="" else ", ")+smets+("" if smets=="" else ", ")+smlls+("" if smlls=="" else ", ")+smljs+("" if smljs=="" else ", ")+sdrljs+("" if sdrljs=="" else ", ")+sptlls+("" if sptlls=="" else ", ")+sssids;
        return spts+", "+snbs+", "+smts+", "+smets+", "+smlls+", "+smljs+", "+sdrljs+", "+sptlls+", "+sssids;

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
    if "SRSSee" in sel:
        return "SRSSee"
    if "SRSSem" in sel:
        return "SRSSem"
    if "SRSSmm" in sel:
        return "SRSSmm"
    if "SR0SFOS" in sel:
        return "0SFOS"
    if "SR1SFOS" in sel:
        return "1SFOS"
    if "SR2SFOS" in sel:
        return "2SFOS"
    if "SRSS1Jee" in sel:
        return "SS1Jee"
    if "SRSS1Jem" in sel:
        return "SS1Jem"
    if "SRSS1Jmm" in sel:
        return "SS1Jmm"
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
    relunc=[0.25,0.3,1.,0.5,1.]
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
    relunc=[0.25,0.3,1.,0.5,1.]
    #relunc=[0.,0.,0.,0.,0.]
    
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
    relunc=[0.25,0.3,1.,0.5,1.]
    #relunc=[0.,0.,0.,0.,0.]
    
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
    relunc=[0.25,0.3,1.,0.5,1.]
    #relunc=[0.,0.,0.,0.,0.]

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
        #if S >= Smin and SB >= SBmin and SErr/S < Serrmax:
        if allcombined[i][3] >= Smin and mySB >= SBmin and allcombined[i][5] < Serrmax*allcombined[i][3]:
            #PrintSB(SB,S,SErr,B,BErr,i,selection)
           PrintSB(allcombined[i][0],allcombined[i][3],allcombined[i][5],allcombined[i][4],allcombined[i][6],allcombined[i][7],selection)
    #print(allcombined)

def GetNumBins(files,histname):
    h = files[0].Get(histname)
    return h.GetNbinsX()


#path = "new_vEval/"
#path = "hists/combineyears/new_vEval/"
path = "hists/combineyears_v5.2.0/new_vEval3/"
#path = "hists/WWW2016_v5.2.0/new2016_vEval2/"
#path = "hists/WWW2017_v5.2.0/new2017_vEval2/"
#path = "hists/WWW2018_v5.2.0/new2018_vEval2/"
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
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjIn__CutsSS",              0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjIn__CutsSS",        1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjIn__CutsSS",    2*9+1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjIn__CutsSS",              0,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjIn__CutsSS",        1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjIn__CutsSS",    2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSSeeMjjIn__CutsSS",1.5,0.01,0.3)
print("************* em-Mjj in *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjIn__CutsSS",              0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjIn__CutsSS",        1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjIn__CutsSS",              0,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjIn__CutsSS",        1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSSemMjjIn__CutsSS",3.7,SB,0.15)
print("************* mm-Mjj in *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSmmMjjIn__CutsSS",              0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSmmMjjIn__CutsSS",        1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSmmMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSmmMjjIn__CutsSS",              0,-1,-1,0.99)
GetAllyields(files,"V2SRSSmmMjjIn__CutsSS",        1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSmmMjjIn__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSSmmMjjIn__CutsSS",4,SB,0.15)

print("************* ee-Mjj out *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjOut__CutsSS",              0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjOut__CutsSS",        1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjOut__CutsSS",1*3    +1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSeeMjjOut__CutsSS",2*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjOut__CutsSS",              0,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjOut__CutsSS",        1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjOut__CutsSS",1*3    +1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSeeMjjOut__CutsSS",2*3+2*9+1*54000,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSSeeMjjOut__CutsSS",2,SB,0.2)
print("************* em-Mjj out *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjOut__CutsSS",              0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjOut__CutsSS",        1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjOut__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSemMjjOut__CutsSS",2*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjOut__CutsSS",              0,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjOut__CutsSS",        1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjOut__CutsSS",1*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSemMjjOut__CutsSS",2*3+2*9+1*54000,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSSemMjjOut__CutsSS",4,SB)
print("************* mm-Mjj out *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSmmMjjOut__CutsSS",              0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSmmMjjOut__CutsSS",        1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSmmMjjOut__CutsSS",    2*9+1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSSmmMjjOut__CutsSS",2*3+2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSmmMjjOut__CutsSS",              0,-1,-1,0.99)
GetAllyields(files,"V2SRSSmmMjjOut__CutsSS",        1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSmmMjjOut__CutsSS",    2*9+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSSmmMjjOut__CutsSS",2*3+2*9+1*54000,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSSmmMjjOut__CutsSS",4,SB,0.15)

print("************* ee-1j *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jee1JPre__Cuts1JSS",                         0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jee1JPre__Cuts1JSS",                   1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jee1JPre__Cuts1JSS",1*2+3*3+1*45+2*180+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jee1JPre__Cuts1JSS",                         0,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jee1JPre__Cuts1JSS",                   1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jee1JPre__Cuts1JSS",1*2+3*3+1*45+2*180+1*54000,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSS1Jee1JPre__Cuts1JSS",2.5,SB,0.3)
print("************* em-1j *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jem1JPre__Cuts1JSS",                         0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jem1JPre__Cuts1JSS",                   1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jem1JPre__Cuts1JSS",1*2+3*3+1*45+2*180+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jem1JPre__Cuts1JSS",                         0,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jem1JPre__Cuts1JSS",                   1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jem1JPre__Cuts1JSS",1*2+3*3+1*45+2*180+1*54000,-1,-1,0.99)
##GetAllyields(files,"V2SRSS1Jem1JPre__Cuts1JSS",0+2+24+3*360+1+4,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSS1Jem1JPre__Cuts1JSS",4,SB,0.1)
print("************* mm-1j *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jmm1JPre__Cuts1JSS",                         0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jmm1JPre__Cuts1JSS",                   1*54000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SRSS1Jmm1JPre__Cuts1JSS",1*2+3*3+1*45+2*180+1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jmm1JPre__Cuts1JSS",                         0,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jmm1JPre__Cuts1JSS",                   1*54000,-1,-1,0.99)
GetAllyields(files,"V2SRSS1Jmm1JPre__Cuts1JSS",1*2+3*3+1*45+2*180+1*54000,-1,-1,0.99)
##GetAllyields(files,"V2SRSS1Jmm1JPre__Cuts1JSS",0+2*3+24+3*360+1,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SRSS1Jmm1JPre__Cuts1JSS",4,SB,0.12)


print("************* 0 SFOS *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",                           0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",                    1*108000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+    1*3+2*9+3*720+1*108000,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+2*1+1*3+3*9+2*720+1*108000,-1,-1,0.99)
GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",                           0,-1,-1,0.99)
GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",                    1*108000,-1,-1,0.99)
GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+    1*3+2*9+3*720+1*108000,-1,-1,0.99)
GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+2*1+1*3+3*9+2*720+1*108000,-1,-1,0.99)
#GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0,-1,-1,0.99)
#GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+1*3+2*9+3*576,-1,-1,0.99)
#GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+2*1+1*3+3*9+2*576,-1,-1,0.99)
#S,B,SB,pSB=GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SR0SFOSDYVeto__Cuts0SFOS",4.,SB,0.15)
###GetSignifForSelectionScanSorted(files,"V2SR0SFOSDYVeto__Cuts0SFOS",2.5,1.5,0.2)
#GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",8659,-1,-1,0.99)
#GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",2880,-1,-1,0.99)
#GetAllyields(files,"V2SR0SFOSDYVeto__Cuts0SFOS",5760,-1,-1,0.99)

#print("*************")
#GetSignifForSelectionScanSorted(files,"V2SR0SFOSDYVeto__Cuts0SFOS",5.5,SB,0.15)
##S,B,SB,pSB=GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+2*1+1*3+3*9+2*576,-1,-1,0.99)
##GetSignifForSelectionScanSorted(files,"V2SR0SFOSDYVeto__Cuts0SFOS",4.0,1.5,0.2,"pSB")
#GetSignifForSelectionScanSorted(files,"V2SR0SFOSDYVeto__Cuts0SFOS",1.0,0.1,0.9,"pSB")
print("************* 1 SFOS *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SR1SFOSDYVeto__Cuts1SFOS",                                        0,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR1SFOSDYVeto__Cuts1SFOS",                                 1*777600,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR1SFOSDYVeto__Cuts1SFOS",    1*12+1*36+1*144+3*720+3*4320+1*777600,-1,-1,0.99)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR1SFOSDYVeto__Cuts1SFOS",2*1+1*12+2*36+1*144+      2*4320+1*777600,-1,-1,0.99)
GetAllyields(files,"V2SR1SFOSDYVeto__Cuts1SFOS",                                        0,-1,-1,0.99)
GetAllyields(files,"V2SR1SFOSDYVeto__Cuts1SFOS",                                 1*777600,-1,-1,0.99)
GetAllyields(files,"V2SR1SFOSDYVeto__Cuts1SFOS",    1*12+1*36+1*144+3*720+3*4320+1*777600,-1,-1,0.99)
GetAllyields(files,"V2SR1SFOSDYVeto__Cuts1SFOS",2*1+1*12+2*36+1*144+      2*4320+1*777600,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SR1SFOSDYVeto__Cuts1SFOS",4,SB,0.1)
print("************* 2 SFOS *************")
S,B,SB,pSB=GetSignifForSelection(files,"V2SR2SFOSDYVeto__Cuts1SFOS",                                        0,-1,-1,0.25)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR2SFOSDYVeto__Cuts1SFOS",                                 1*777600,-1,-1,0.25)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR2SFOSDYVeto__Cuts1SFOS",    1*12+     2*144+3*720+3*4320+1*777600,-1,-1,0.25)
S,B,SB,pSB=GetSignifForSelection(files,"V2SR2SFOSDYVeto__Cuts1SFOS",2*1+1*12+2*36+2*144+      2*4320+1*777600,-1,-1,0.25)
GetAllyields(files,"V2SR2SFOSDYVeto__Cuts1SFOS",                                        0,-1,-1,0.99)
GetAllyields(files,"V2SR2SFOSDYVeto__Cuts1SFOS",                                 1*777600,-1,-1,0.99)
GetAllyields(files,"V2SR2SFOSDYVeto__Cuts1SFOS",    1*12+     2*144+3*720+3*4320+1*777600,-1,-1,0.99)
GetAllyields(files,"V2SR2SFOSDYVeto__Cuts1SFOS",2*1+1*12+2*36+1*144+      2*4320+1*777600,-1,-1,0.99)
GetSignifForSelectionScanSorted(files,"V2SR2SFOSDYVeto__Cuts1SFOS",3.5,SB,0.17)

print("************* extras *************")
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0          ,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+    10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  2*10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  3*10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  4*10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  5*10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  6*10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  7*10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  8*10800,-1,-1,0.99)
GetSignifForSelection(files,"V2SR0SFOSDYVeto__Cuts0SFOS",0+  9*10800,-1,-1,0.99)
