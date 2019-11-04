import ROOT
import math

relunc=[0.25,0.30,1.0,0.5,1.0]
relunc=[0.25,0.25,0.5,0.5,0.5]


def getyieldanderror1D(hS,mybin):
    if mybin < 0:
        return -4.,-4
    if mybin > hS.GetNbinsX():
        return -5.,-5
    binnum = hS.FindBin(mybin-0.0001) #round down
    err = ROOT.Double(-1.)
    sig =  hS.IntegralAndError(1,binnum,err)
    return sig,err

def getyieldanderror(hS,mybinX,mybinY):
    if mybinX < 0:
        return -4.,-4
    if mybinX > hS.GetNbinsX():
        return -5.,-5
    if mybinY < 0:
        return -6.,-6
    if mybinY > hS.GetNbinsY():
        return -7.,-7
    binnumX = hS.GetXaxis().FindBin(mybinX-0.0001) #round down
    binnumY = hS.GetYaxis().FindBin(mybinY-0.0001) #round down
    err = ROOT.Double(-1.)
    sig =  hS.IntegralAndError(1,binnumX,1,binnumY,err)
    return sig,err

def getBkgHist1D(hbkg,relunc):
    h = hbkg[0].Clone("hB")
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
        h.SetBinContent(i,P+L+F+Q+G)
        h.SetBinError(i,math.sqrt(pow(Pe,2)+pow(Le,2)+pow(Fe,2)+pow(Qe,2)+pow(Ge,2)))
    return h

#def getBkgHist2D(hbkg,relunc):
def getBkgHist(hbkg):
    h = hbkg[0].Clone("hB")
    for i in range(1,len(hbkg)):
        h.Add(hbkg[i])
    #if len(hbkg) != len(relunc):
    #    print("hbkg has different length than relunc",len(hbkg),len(relunc))
    #    return h
    if len(hbkg)==0:
        print("you need to provide signal and background histogram lists that are not empty")
        return h
    """
    for i in range(0,h.GetNbinsX()+1):
        for j in range(0,h.GetNbinsY()+1):
            P =  hbkg[0].GetBinContent(i,j)
            L =  hbkg[1].GetBinContent(i,j)
            F =  hbkg[2].GetBinContent(i,j)
            Q =  hbkg[3].GetBinContent(i,j)
            G =  hbkg[4].GetBinContent(i,j)
            Pe =  hbkg[0].GetBinError(i,j)
            Le =  hbkg[1].GetBinError(i,j)
            Fe =  hbkg[2].GetBinError(i,j)
            Qe =  hbkg[3].GetBinError(i,j)
            Ge =  hbkg[4].GetBinError(i,j)
            
            Pe = math.sqrt(Pe*Pe+pow((relunc[0])*P,2))
            Le = math.sqrt(Le*Le+pow((relunc[1])*L,2))
            Fe = math.sqrt(Fe*Fe+pow((relunc[2])*F,2))
            Qe = math.sqrt(Qe*Qe+pow((relunc[3])*Q,2))
            Ge = math.sqrt(Ge*Ge+pow((relunc[4])*G,2))
            h.SetBinContent(i,j,P+L+F+Q+G)
            h.SetBinError(i,j,math.sqrt(pow(Pe,2)+pow(Le,2)+pow(Fe,2)+pow(Qe,2)+pow(Ge,2)))
    """
    return h

def signif(sig,bkg,bkgerr):
    if bkg == 0:
        return -1;
    #print(sig,bkg,bkgerr,bkgerr*bkgerr,math.sqrt(bkg+pow(bkgerr,2)))
    if bkg+pow(bkgerr,2) <= 0.01:
        return -1
    return sig/math.sqrt(bkg+pow(bkgerr,2))

#def GetSB(hsig,hbkg,hB,relunc,mybin):
def GetSB1D(hsig,hbkg,hB,relunc,mybin):
    if len(hbkg) != len(relunc):
        print("hbkg has different length than relunc",len(hbkg),len(relunc))
        return -3.,0.,0.,0.,0.
    if len(hsig)==0 or len(hbkg)==0:
        print("you need to provide signal and background histogram lists that are not empty")
        return -2.,0.,0.,0.,0.
    hS = hsig[0].Clone("hS")
    for i in range(1,len(hsig)):
        hS.Add(hsig[i])
    if mybin < 0:
        return -4.,0.,0.,0.,0.
    if mybin > hS.GetNbinsX():
        return -5.,0.,0.,0.,0.
    s,serr = getyieldanderror1D(hS,mybin)
    if s < 0:
        return -1.5,0.,0.,0.,0.
    b,berr = getyieldanderror1D(hB,mybin)
    for i in range(1,len(hbkg)):
        temp,temperr =  getyieldanderror1D(hbkg[i],mybin)
        berr = math.sqrt(pow(berr,2)+pow(relunc[i]*temp,2))
    sb = signif(s,b,berr)
    pSB = 0
    if b>0:
        pSB = math.sqrt(2*((s+b)*math.log(1+s/b)-s))
    return sb,s,serr,b,berr,pSB

def GetSB(hsig,hbkg,hB,relunc,mybinX,mybinY):
    if len(hbkg) != len(relunc):
        print("hbkg has different length than relunc",len(hbkg),len(relunc))
        return -3.,0.,0.,0.,0.
    if len(hsig)==0 or len(hbkg)==0:
        print("you need to provide signal and background histogram lists that are not empty")
        return -2.,0.,0.,0.,0.
    hS = hsig[0].Clone("hS")
    for i in range(1,len(hsig)):
        hS.Add(hsig[i])
    if mybinX < 0 or mybinY < 0:
        return -4.,0.,0.,0.,0.
    if mybinX > hS.GetNbinsX() or mybinY > hS.GetNbinsY():
        return -5.,0.,0.,0.,0.
    s,serr = getyieldanderror(hS,mybinX,mybinY)
    if s < 0:
        return -1.5,0.,0.,0.,0.
    b,berr = getyieldanderror(hB,mybinX,mybinY)#XX
    for i in range(1,len(hbkg)):
        temp,temperr =  getyieldanderror(hbkg[i],mybinX,mybinY)
        berr = math.sqrt(pow(berr,2)+pow(relunc[i]*temp,2))
    #binX = hS.GetXaxis().FindBin(mybinX-0.0001)
    #binY = hS.GetYaxis().FindBin(mybinY-0.0001)
    #print(mybinX,mybinY,binX,binY)
    #error = ROOT.Double(-1)
    #print("getsb","S",hS.Integral(1,binX,1,binY),"B",hbkg[0].Integral(1,binX,1,binY)+hbkg[1].Integral(1,binX,1,binY)+hbkg[2].Integral(1,binX,1,binY)+hbkg[3].Integral(1,binX,1,binY)+hbkg[4].Integral(1,binX,1,binY),"Bsyst",math.sqrt(pow(0.25*hbkg[0].Integral(1,binX,1,binY),2)+pow(0.3*hbkg[1].Integral(1,binX,1,binY),2)+pow(1.0*hbkg[2].Integral(1,binX,1,binY),2)+pow(0.5*hbkg[3].Integral(1,binX,1,binY),2)+pow(1.*hbkg[4].Integral(1,binX,1,binY),2)),"hB",hB.IntegralAndError(1,binX,1,binY,error),"hBerr",error,"berr",berr)
    
    sb = signif(s,b,berr)
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
    if "SR0SFOSeem" in sel:
        return "0SFOS,eem"
    if "SR0SFOSemm" in sel:
        return "0SFOS,emm"
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

def getcut1D(sel,mybin):
    if "ele_relIso" in sel:
        return "Iso(el) < "+str('%.4f' % mybin)
    if "muo_relIso" in sel:
        return "Iso(mu) < "+str('%.4f' % mybin)
    return ""
def getcut(sel,mybinX,mybinY):
    mytype = ""
    if "ele_relIso" in sel:
        mytype = mytype + "Iso(el) < "+str('%.4f' % mybinX)+ "  "
    if "muo_relIso" in sel:
        mytype = mytype + "Iso(mu) < "+str('%.4f' % mybinY)+ "  "
    return mytype

def getselection(sel):
    if (("Presel" in sel) or ("PreSel" in sel) or ("DYVeto" in sel)):
        return "Preselection, 0 b-tag(loose)             "
    if "NsoftbVeto" in sel:
        return "Preselection, 0 b-tag(loose) and 0 soft-b"
    return "Preselection, no cuts on b-tagging       "


def PrintSB1D(sb,s,serr,b,berr,f,ferr,pSB,mybin,sel):
    if b<=0:
        return True
    SBs = "{0:.3g}".format(sb)
    ss    = "{0:.2g}".format(s)
    serrs = "{0:.2g}".format(serr)
    bs    = "{0:.2g}".format(b)
    berrs = "{0:.2g}".format(berr)
    fs    = "{0:.2g}".format(f)
    ferrs = "{0:.2g}".format(ferr)
    sobs  = "{0:.3g}".format(s/b)
    sosbs = "{0:.3g}".format(s/math.sqrt(b))
    spsbs  = "{0:.3g}".format(pSB)
    mytype = gettype(sel)
    if mytype=="":
        print("type not recognized from histogram",sel)
        return False
    cut = getcut1D(sel,mybin)
    selection = getselection(sel)
    if cut != "":
        print (mytype+", "+selection+" and "+cut+"- S/sqrt(B+dB^2) = "+SBs.ljust(5)+", Philip's Sig = "+spsbs.ljust(5)+", S/sqrt(B) = "+sosbs.ljust(5)+", S/B = "+sobs.ljust(6)+", S="+ss.rjust(4)+"+/-"+serrs.ljust(4)+", B="+bs.rjust(5)+"+/-"+berrs.ljust(4)+", F="+fs.rjust(5)+"+/-"+ferrs.ljust(4))
    return True

def PrintSB(sb,s,serr,b,berr,f,ferr,pSB,mybinX,mybinY,sel):
    if b<=0:
        return True
    SBs = "{0:.3g}".format(sb)
    ss    = "{0:.2g}".format(s)
    serrs = "{0:.2g}".format(serr)
    bs    = "{0:.2g}".format(b)
    berrs = "{0:.2g}".format(berr)
    fs    = "{0:.2g}".format(f)
    ferrs = "{0:.2g}".format(ferr)
    sobs  = "{0:.3g}".format(s/b)
    sosbs = "{0:.3g}".format(s/math.sqrt(b))
    spsbs  = "{0:.3g}".format(pSB)
    mytype = gettype(sel)
    if mytype=="":
        print("type not recognized from histogram",sel)
        return False
    cut = getcut(sel,mybinX,mybinY)
    selection = getselection(sel)
    if cut != "":
        print (mytype+", "+selection+" and "+cut+"- S/sqrt(B+dB^2) = "+SBs.ljust(5)+", Philip's Sig = "+spsbs.ljust(5)+", S/sqrt(B) = "+sosbs.ljust(5)+", S/B = "+sobs.ljust(6)+", S="+ss.rjust(4)+"+/-"+serrs.ljust(4)+", B="+bs.rjust(5)+"+/-"+berrs.ljust(4)+", F="+fs.rjust(5)+"+/-"+ferrs.ljust(4))
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
    cut = getcut(sel,mybin)
    selection = getselection(sel)
    if cut != "":
        print (mytype+", "+selection+" and "+cut+"- S="+ss.rjust(4)+"+/-"+serrs.ljust(4)+", B="+bs.rjust(5)+"+/-"+berrs.ljust(4)+" (LL="+ls.rjust(5)+"+/-"+lerrs.ljust(4)+", Pr="+ps.rjust(4)+"+/-"+perrs.ljust(4)+", F="+fs.rjust(4)+"+/-"+ferrs.ljust(4)+", Qf="+qs.rjust(4)+"+/-"+qerrs.ljust(4)+", gf="+gs.rjust(4)+"+/-"+gerrs.ljust(4)+")")
    return True

def Printyields(S,SErr,B,BErr,P,PErr,L,LErr,F,FErr,Q,QErr,G,GErr,mybinX,mybinY,sel):

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
    cut = getcut(sel,mybinX,mybinY)
    selection = getselection(sel)
    if cut != "":
        print (mytype+", "+selection+" and "+cut+"- S="+ss.rjust(4)+"+/-"+serrs.ljust(4)+", B="+bs.rjust(5)+"+/-"+berrs.ljust(4)+" (LL="+ls.rjust(5)+"+/-"+lerrs.ljust(4)+", Pr="+ps.rjust(4)+"+/-"+perrs.ljust(4)+", F="+fs.rjust(4)+"+/-"+ferrs.ljust(4)+", Qf="+qs.rjust(4)+"+/-"+qerrs.ljust(4)+", gf="+gs.rjust(4)+"+/-"+gerrs.ljust(4)+")")
    return True

def GetAllyields1D(files,selection,mybin,Smin=-1,SBmin=-1,Serrmax=0.1):

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
    hB = getBkgHist(hbkg)
    P, Perr = getyieldanderror1D(hP,mybin)
    L, Lerr = getyieldanderror1D(hL,mybin)
    F, Ferr = getyieldanderror1D(hF,mybin)
    Q, Qerr = getyieldanderror1D(hQ,mybin)
    G, Gerr = getyieldanderror1D(hG,mybin)
    
    
    SB,S,SErr,B,BErr,pSB = GetSB1D(hsig,hbkg,hB,relunc,mybin)
    if S >= Smin and SB >= SBmin and SErr < Serrmax*S:
        Printyields(S,SErr,B,BErr,P,Perr,L,Lerr,F,Ferr,Q,Qerr,G,Gerr,mybin,selection)
    return S,B,SB,pSB

def GetAllyields(files,selection,mybinX,mybinY,Smin=-1,SBmin=-1,Serrmax=0.1):

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
    hB = getBkgHist(hbkg)
    P, Perr = getyieldanderror(hP,mybinX,mybinY)
    L, Lerr = getyieldanderror(hL,mybinX,mybinY)
    F, Ferr = getyieldanderror(hF,mybinX,mybinY)
    Q, Qerr = getyieldanderror(hQ,mybinX,mybinY)
    G, Gerr = getyieldanderror(hG,mybinX,mybinY)

    SB,S,SErr,B,BErr,pSB = GetSB(hsig,hbkg,hB,relunc,mybinX,mybinY)
    if S >= Smin and SB >= SBmin and SErr < Serrmax*S:
        Printyields(S,SErr,B,BErr,P,Perr,L,Lerr,F,Ferr,Q,Qerr,G,Gerr,mybinX,mybinY,selection)
    return S,B,SB,pSB



def GetSignifForSelection1D(files,selection,mybin,Smin=-1,SBmin=-1,Serrmax=0.1):

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
    ##relunc=[0.25,0.3,1.,0.5,1.]
    ##relunc=[0.,0.,0.,0.,0.]
    hB = getBkgHist(hbkg)

    SB,S,SErr,B,BErr,pSB = GetSB1D(hsig,hbkg,hB,relunc,mybin)
    F,Fe = getyieldanderror1D(hF,mybin)
    if S >= Smin and SB >= SBmin and SErr < Serrmax*S:
        PrintSB1D(SB,S,SErr,B,BErr,F,Fe,pSB,mybin,selection)
    return S,B,SB,pSB

def GetSignifForSelection(files,selection,mybinX,mybinY,Smin=-1,SBmin=-1,Serrmax=0.1):

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
    ##relunc=[0.25,0.3,1.,0.5,1.]
    ##relunc=[0.,0.,0.,0.,0.]
    hB = getBkgHist(hbkg)
    #mybin = hS.FindBin(mybinX-0.0001,mybinY-0.0001)
    #binX = hS.GetXaxis().FindBin(mybinX-0.0001)
    #binY = hS.GetYaxis().FindBin(mybinY-0.0001)
    #print(mybinX,mybinY,binX,binY)
    #error = ROOT.Double(-1)
    #print(selection,"S",hS.Integral(1,binX,1,binY),"B",hP.Integral(1,binX,1,binY)+hL.Integral(1,binX,1,binY)+hF.Integral(1,binX,1,binY)+hQ.Integral(1,binX,1,binY)+hG.Integral(1,binX,1,binY),"Bsyst",math.sqrt(pow(0.25*hP.Integral(1,binX,1,binY),2)+pow(0.3*hL.Integral(1,binX,1,binY),2)+pow(1.0*hL.Integral(1,binX,1,binY),2)+pow(0.5*hQ.Integral(1,binX,1,binY),2)+pow(1.*hG.Integral(1,binX,1,binY),2)),"hB",hB.IntegralAndError(1,binX,1,binY,error),"hBerr",error)
    
    SB,S,SErr,B,BErr,pSB = GetSB(hsig,hbkg,hB,relunc,mybinX,mybinY)
    F,Fe = getyieldanderror(hF,mybinX,mybinY)
    if S >= Smin and SB >= SBmin and SErr < Serrmax*S:
        PrintSB(SB,S,SErr,B,BErr,F,Fe,pSB,mybinX,mybinY,selection)
    return S,B,SB,pSB


def GetSignifForSelectionScan1D(files,selection,Smin=-1,SBmin=-1,Serrmax=0.1):

    hS = files[0].Get(selection)
    #print(hS)
    hP = files[1].Get(selection)
    hL = files[2].Get(selection)
    hF = files[3].Get(selection)
    hQ = files[4].Get(selection)
    hG = files[5].Get(selection)
    
    hsig = [hS]
    hbkg = [hP,hL,hF,hQ,hG]
    #relunc=[0.25,0.3,1.,0.5,1.]
    ##relunc=[0.,0.,0.,0.,0.]
    hB = getBkgHist(hbkg)

    for i in range(0,hS.GetNbinsX()):
        mybin = hS.GetXaxis().GetBinLowEdge(i)
        #SB,S,SErr,B,BErr,pSB = GetSB1D(hsig,hbkg,hB,relunc,mybin)
        SB,S,SErr,B,BErr,pSB = GetSB1D(hsig,hbkg,hB,relunc,mybin)
        pSB = 0
        F,Fe = getyieldanderror1D(hF,mybin)
        if S >= Smin and SB >= SBmin and SErr/S < Serrmax:
            PrintSB1D(SB,S,SErr,B,BErr,F,Fe,pSB,mybin,selection)
    return True

def GetSignifForSelectionScan2D(files,selection,Smin=-1,SBmin=-1,Serrmax=0.1):

    hS = files[0].Get(selection)
    #print(hS)
    hP = files[1].Get(selection)
    hL = files[2].Get(selection)
    hF = files[3].Get(selection)
    hQ = files[4].Get(selection)
    hG = files[5].Get(selection)
    
    hsig = [hS]
    hbkg = [hP,hL,hF,hQ,hG]
    #relunc=[0.25,0.3,1.,0.5,1.]
    ##relunc=[0.,0.,0.,0.,0.]
    hB = getBkgHist(hbkg)

    for i in range(0,hS.GetNbinsX()):
        mybinX = hS.GetXaxis().GetBinLowEdge(i)+hS.GetXaxis().GetBinWidth(i)
        if mybinX>0.1025:
            continue
        for j in range(0,hS.GetNbinsY()):
            mybinY = hS.GetYaxis().GetBinLowEdge(j)+hS.GetYaxis().GetBinWidth(j)
            if mybinY>0.1525:
                continue
            SB,S,SErr,B,BErr,pSB = GetSB(hsig,hbkg,hB,relunc,mybinX,mybinY)
            F,Fe = getyieldanderror(hF,mybinX,mybinY)
            if S >= Smin and SB >= SBmin and SErr/S < Serrmax:
                PrintSB(SB,S,SErr,B,BErr,F,Fe,pSB,mybinX,mybinY,selection)
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

def GetSignifForSelectionScanSorted1D(files,selection,Smin=-1,SBmin=-1,Serrmax=0.1,sortwhat="S/sqrt(B+dB^2)"):
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
    hB = getBkgHist(hbkg)


    listSB, listSoB, listS, listB, listSErr, listBErr,listi, listpSB,listF,listFErr = [],[],[],[],[],[],[],[],[],[]
    for i in range(0,hS.GetNbinsX()):
    #for i in range(0,100):
        SB,S,SErr,B,BErr,pSB = GetSB1D(hsig,hbkg,hB,relunc,i)
        F,Fe = getyieldanderror1D(hF,i)
        SoB = -1
        if B>0: SoB = S/B
        listSB.append(SB)
        listSoB.append(SoB)
        listS.append(S)
        listB.append(B)
        listF.append(F)
        listSErr.append(SErr)
        listBErr.append(BErr)
        listFErr.append(Fe)
        listi.append(i)
        listpSB.append(pSB)
        #PrintSB(SB,S,SErr,B,BErr,F,Fe,pSB,i,selection)

    allcombined = zip(listSB,listSoB,listpSB,listS,listB,listSErr,listBErr,listi,listF,listFErr)
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
        #    PrintSB(allcombined[i][0],allcombined[i][3],allcombined[i][5],allcombined[i][4],allcombined[i][6],allcombined[i][8],allcombined[i][9],allcombined[i][2],allcombined[i][7],selection)
        #if S >= Smin and SB >= SBmin and SErr/S < Serrmax:
        if allcombined[i][3] >= Smin and mySB >= SBmin and allcombined[i][5] < Serrmax*allcombined[i][3]:
            #PrintSB(SB,S,SErr,B,BErr,F,Fe,pSB,i,selection)
            #continue
            PrintSB(allcombined[i][0],allcombined[i][3],allcombined[i][5],allcombined[i][4],allcombined[i][6],allcombined[i][8],allcombined[i][9],allcombined[i][2],allcombined[i][7],selection)
    #print(allcombined)

def GetSignifForSelectionScanSorted2D(files,selection,Smin=-1,SBmin=-1,Serrmax=0.1,sortwhat="S/sqrt(B+dB^2)"):
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
    hB = getBkgHist(hbkg)


    listSB, listSoB, listS, listB, listSErr, listBErr,listi,listj, listpSB,listF,listFErr = [],[],[],[],[],[],[],[],[],[],[]
    for i in range(0,hS.GetNbinsX()):
        for j in range(0,hS.GetNbinsY()):
            #for i in range(0,10):
            #for j in range(0,10):
            SB,S,SErr,B,BErr,pSB = GetSB(hsig,hbkg,hB,relunc,i,j)
            F,Fe = getyieldanderror(hF,i,j)
            SoB = -1
            if B>0: SoB = S/B
            listSB.append(SB)
            listSoB.append(SoB)
            listS.append(S)
            listB.append(B)
            listF.append(F)
            listSErr.append(SErr)
            listBErr.append(BErr)
            listFErr.append(FErr)
            listi.append(i)
            listj.append(j)
            listpSB.append(pSB)
            #PrintSB(SB,S,SErr,B,BErr,F,Fe,pSB,i,j,selection)

    allcombined = zip(listSB,listSoB,listpSB,listS,listB,listSErr,listBErr,listi,listj,listF,listFErr)
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
        #    PrintSB(allcombined[i][0],allcombined[i][3],allcombined[i][5],allcombined[i][4],allcombined[i][6],allcombined[i][9],allcombined[i][10],allcombined[i][2],allcombined[i][7],allcombined[i][8],selection)
        #if S >= Smin and SB >= SBmin and SErr/S < Serrmax:
        if allcombined[i][3] >= Smin and mySB >= SBmin and allcombined[i][5] < Serrmax*allcombined[i][3]:
            #PrintSB(SB,S,SErr,B,BErr,F,Fe,pSB,i,selection)
            #continue
            PrintSB(allcombined[i][0],allcombined[i][3],allcombined[i][5],allcombined[i][4],allcombined[i][6],allcombined[i][9],allcombined[i][10],allcombined[i][2],allcombined[i][7],allcombined[i][8],selection)
    #print(allcombined)


def GetNumBins(files,histname):
    h = files[0].Get(histname)
    return h.GetNbinsX()


path = "hists/combineyears_v5.3.0/newIso/"
#path = "hists/Loose2017_v5.3.0/newIso/"
fS = ROOT.TFile.Open(path+"signal.root")
#fS = ROOT.TFile.Open(path+"signal_private.root")
fP = ROOT.TFile.Open(path+"prompt.root")
fL = ROOT.TFile.Open(path+"lostlep.root")
fF = ROOT.TFile.Open(path+"fakes.root")
#fF = ROOT.TFile.Open(path+"ddfakes.root")
fQ = ROOT.TFile.Open(path+"qflip.root")
fG = ROOT.TFile.Open(path+"photon.root")
files = [fS,fP,fL,fF,fQ,fG]
#GetSignifForSelection(files,"SR0SFOS__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)


print("************* ee-2j *************")
#GetSignifForSelectionScan1D(files,"SRSSee__ele_relIso03EAXXMax",              -1,-1,0.99)
#GetSignifForSelection(files,"SRSSee__ele_relIso03EAXXMax",0.10,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan1D(files,"SRSSeePreSel__ele_relIso03EAXXMax",        -1,-1,0.99)
#GetSignifForSelection(files,"SRSSeePreSel__ele_relIso03EAXXMax",0.10,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan1D(files,"SRSSeeNsoftbVeto__ele_relIso03EAXXMax",    -1,-1,0.99)
#GetSignifForSelection(files,"SRSSeeNsoftbVeto__ele_relIso03EAXXMax",0.10,Smin=-1,SBmin=-1,Serrmax=0.99)

print("************* em-2j *************")
##GetSignifForSelectionScan1D(files,"SRSSem__ele_relIso03EAXXMax",              -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSSemPreSel__ele_relIso03EAXXMax",        -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSSemNsoftbVeto__ele_relIso03EAXXMax",    -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSSem__muo_relIso03EAXXMax",              -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSSemPreSel__muo_relIso03EAXXMax",        -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSSemNsoftbVeto__muo_relIso03EAXXMax",    -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SRSSem__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",              -1,-1,0.99)
#GetSignifForSelection(files,"SRSSem__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SRSSemPreSel__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",        -1,-1,0.99)
#GetSignifForSelection(files,"SRSSemPreSel__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SRSSemNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",    -1,-1,0.99)
#GetSignifForSelection(files,"SRSSemNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)

print("************* mm-2j *************")
#GetSignifForSelectionScan1D(files,"SRSSmm__muo_relIso03EAXXMax",              -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SRSSmm__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",              -1,-1,0.99)
#GetSignifForSelection1D(files,"SRSSmm__muo_relIso03EAXXMax",0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan1D(files,"SRSSmmPreSel__muo_relIso03EAXXMax",        -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SRSSmmPreSel__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",              -1,-1,0.99)
#GetSignifForSelection1D(files,"SRSSmmPreSel__muo_relIso03EAXXMax",0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan1D(files,"SRSSmmNsoftbVeto__muo_relIso03EAXXMax",    -1,-1,0.99)
GetSignifForSelectionScan2D(files,"SRSSmmNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",              -1,-1,0.99)
GetSignifForSelection1D(files,"SRSSmmNsoftbVeto__muo_relIso03EAXXMax",0.15,Smin=-1,SBmin=-1,Serrmax=0.99)



print("************* ee-1j *************")
#GetSignifForSelectionScan1D(files,"SRSS1Jee1JPre__ele_relIso03EAXXMax",        -1,-1,0.99)
#GetSignifForSelection(files,"SRSS1Jee1JPre__ele_relIso03EAXXMax",0.10,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan1D(files,"SRSS1JeeNsoftbVeto__ele_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelection(files,"SRSS1JeeNsoftbVeto__ele_relIso03EAXXMax",0.10,Smin=-1,SBmin=-1,Serrmax=0.99)

print("************* em-1j *************")
##GetSignifForSelectionScan1D(files,"SRSS1Jem1JPre__ele_relIso03EAXXMax",        -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSS1JemNsoftbVeto__ele_relIso03EAXXMax",   -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSS1Jem1JPre__muo_relIso03EAXXMax",        -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SRSS1JemNsoftbVeto__muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SRSS1Jem1JPre__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",        -1,-1,0.99)
#GetSignifForSelection(files,"SRSS1Jem1JPre__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SRSS1JemNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelection(files,"SRSS1JemNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)

print("************* mm-1j *************")
#GetSignifForSelectionScan1D(files,"SRSS1Jmm1JPre__muo_relIso03EAXXMax",        -1,-1,0.99)
#GetSignifForSelection(files,"SRSS1Jmm1JPre__muo_relIso03EAXXMax",0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan1D(files,"SRSS1JmmNsoftbVeto__muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelection(files,"SRSS1JmmNsoftbVeto__muo_relIso03EAXXMax",0.15,Smin=-1,SBmin=-1,Serrmax=0.99)



print("************* 0 SFOS *************")
##GetSignifForSelectionScan1D(files,"SR0SFOS__ele_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSDYVeto__ele_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSNsoftbVeto__ele_relIso03EAXXMax",   -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOS__muo_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSDYVeto__muo_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSNsoftbVeto__muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOS__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",             -1,0.1,0.99)
#GetSignifForSelection(files,"SR0SFOS__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",       -1,0.45,0.99)
#GetSignifForSelection(files,"SR0SFOSDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",   -1,1.03,0.99)
#GetSignifForSelection(files,"SR0SFOSNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)

print("************* 0 SFOS,eem *************")
##GetSignifForSelectionScan1D(files,"SR0SFOSeem__ele_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSeemDYVeto__ele_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSeemNsoftbVeto__ele_relIso03EAXXMax",   -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSeem__muo_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSeemDYVeto__muo_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSeemNsoftbVeto__muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSeem__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",             -1,0.225,0.99)
#GetSignifForSelection(files,"SR0SFOSeem__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSeemDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",       -1,0.3,0.99)
#GetSignifForSelection(files,"SR0SFOSeemDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSeemNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",   -1,0.67,0.99)
#GetSignifForSelection(files,"SR0SFOSeemNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
print("************* 0 SFOS,emm *************")
##GetSignifForSelectionScan1D(files,"SR0SFOSemm__ele_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSemmDYVeto__ele_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSemmNsoftbVeto__ele_relIso03EAXXMax",   -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSemm__muo_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSemmDYVeto__muo_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR0SFOSemmNsoftbVeto__muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSemm__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",             -1,0.25,0.99)
#GetSignifForSelection(files,"SR0SFOSemm__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSemmDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",       -1,0.667,0.99)
#GetSignifForSelection(files,"SR0SFOSemmDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR0SFOSemmNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",   -1,0.75,0.99)
#GetSignifForSelection(files,"SR0SFOSemmNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)



print("************* 1 SFOS *************")
##GetSignifForSelectionScan1D(files,"SR1SFOS__ele_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR1SFOSDYVeto__ele_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR1SFOSNsoftbVeto__ele_relIso03EAXXMax",   -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR1SFOS__muo_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR1SFOSDYVeto__muo_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR1SFOSNsoftbVeto__muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SR1SFOS__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",             -1,-1,0.99)
#GetSignifForSelection(files,"SR1SFOS__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR1SFOSDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",       -1,0.15,0.99)
#GetSignifForSelection(files,"SR1SFOSDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR1SFOSNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",   -1,0.15,0.99)
#GetSignifForSelection(files,"SR1SFOSNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)

print("************* 2 SFOS *************")
##GetSignifForSelectionScan1D(files,"SR2SFOS__ele_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR2SFOSDYVeto__ele_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR2SFOSNsoftbVeto__ele_relIso03EAXXMax",   -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR2SFOS__muo_relIso03EAXXMax",             -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR2SFOSDYVeto__muo_relIso03EAXXMax",       -1,-1,0.99)
##GetSignifForSelectionScan1D(files,"SR2SFOSNsoftbVeto__muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelectionScan2D(files,"SR2SFOS__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",             -1,-1,0.99)
#GetSignifForSelection(files,"SR2SFOS__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR2SFOSDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",       -1,-1,0.99)
#GetSignifForSelection(files,"SR2SFOSDYVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)
#GetSignifForSelectionScan2D(files,"SR2SFOSNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",   -1,-1,0.99)
#GetSignifForSelection(files,"SR2SFOSNsoftbVeto__ele_relIso03EAXXMax_v_muo_relIso03EAXXMax",0.10,0.15,Smin=-1,SBmin=-1,Serrmax=0.99)


print("Path",path)
