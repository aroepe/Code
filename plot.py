from ROOT import *
import os


def filePath(inDir,channel):
    f=inDir+channel[1]
    if not os.path.exists(f): f=inDir+channel[1]
    return [TFile(f),channel[0]]

def filePathAttr():
    return "inDir <str>, channel <[str,str]> (first entry is name second entry is file path). returns [TFile(f) <TFile>,channel[0] <str>]"

def oneDplot(tree, name, distribution, bins, xMin, xMax, fiducial, cuts):
    if (name.find('data')!=-1 or not fiducial):
        tree.Draw("%s>>%s(%s,%s,%s)"%(distribution,name,bins,xMin,xMax),"(%s)"%cuts)  
    else:
        tree.Draw("%s>>%s(%s,%s,%s)"%(distribution,name,bins,xMin,xMax),"mcEventWeight*((Sum$(LLP_R_coord>4)>0) && (Sum$(LLP_R_coord<300)>0) && (Sum$(abs(LLP_Eta)<2.1)>0)\
 && (Sum$(abs(LLP_Z_coord)<800)>0))*(%s)"%cuts)
    h=gDirectory.Get(name).Clone(name)
    return h

def oneDplotAttr():
    return " tree <TTree>, name <str>, distribution <str>, bins <int/array>, xMin <float/int>, xMax <float/int>, fiducial <bool>, cuts <str> (set to 1 if there are no cuts). returns h <TObject>"

def twoDplot(tree, name, distribution, xBins, xMax, xMin, yBins, yMin, yMax, fiducial, cuts):
    if (name.find('data')!=-1 or not fiducial):
        t.Draw("%s>>%s(%s,%s,%s,%s,%s,%s)"%(distribution,name,xBins,xMin,xMax,yBins,yMin,yMax),"(%s)"%cuts)
    elif fiducial:
        t.Draw("%s>>%s(%s,%s,%s,%s,%s,%s)"%(distribution,name,xBins,xMin,xMax,yBins,yMin,yMax),"mcEventWeight*((Sum$(LLP_R_coord>4)>0) && (Sum$(LLP_R_coord<300)>0) && (Sum$(abs(LLP_Eta)<2.1)>0) && (Sum$(abs(LLP_Z_coord)<800)>0))*(%s)"%cuts)
    h=gDirectory.Get(name).Clone(name)
    return h

def twoDplotAttr():
    return " tree <TTree>, name <str>, distribution <str>, xBins <int/array>, xMin <float/int>, xMax <float/int>, yBins <int/array>, yMin <float/int>, yMax <float/int>, fiducial <bool>, cuts <str> (set to 1 if there are no cuts), returns h <TObject>"

def integral(histList,dimensions):
    for hist in histList:
        if (dimensions==1):
            hist.Scale(1./hist.Integral(0,hist.GetNbinsX()+1))
        else:
            hist.Scale(1./hist.Integral(0,hist.GetNbinsX()+1,0,hist.GetNbinsY()+1))
        hist.Draw()

def integralAttr():
    return "histList list of <TObject>, dimensions <int>"

def maxMin(histList,staticMax,staticMin):
    maximum=0
    for hist in histList:
        if hist.GetMaximum()>maximum: maximum=hist.GetMaximum()

    for hist in histList:
        if staticMax==-999: hist.SetMaximum(maximum*1.1)
        else: hist.SetMaximum(staticMax)
        hist.SetMinimum(staticMin)

def maxMinAttr():
    return "histList list of <TObject>, staticMax <float/int> (set to -999 for auto calculation of maximum), staticMin <float/int>"

def profile(hist,name,fit,color):
    hist.Draw("colz")
    hist.ProfileX(name)
    name.SetLineColor(color)
    name.Draw("same")
    name.Fit(fit)

def profileAttr():
    return "hist <TObject>, name <str>, fit <str>, color <TColor name>"

def logBins(nBins,xMin,xMax,yMin,yMax):
    binEdgesX=[]
    binEdgesY=[]
    delta=(log10(xMax)-log10(xMin))/nBins
    for i in xrange(0,nBins+1):
            edge=10**(log10(xMin)+delta*i)
            binEdgesX.append(edge)

    delta=(log10(yMax)-log10(yMin))/nBins
    for i in xrange(0,nBins+1):
        edge=10**(yMin+delta*i)                                                                                                                          
        binEdgesY.append(edge)

    ary = array("d",binEdgesY)                                                                          
    arx = array("d",binEdgesX)
    return [arx,ary]

def logBinsAttr():
    return "nBins <int>, xMin <float/int>, xMax <float/int>, yMin <float/int>, yMax <float/int>. returns [arx <array><float>,ary <array><float>]"
