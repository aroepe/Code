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
    return " tree <TTree>, name <str>, distribution <str>, bins <int/array>, xMin <float/int>, xMax <float/int>, fiducial <bool>, cuts <str> (set to 1 if there are no cuts). returns h <TH1>"

def twoDplot(tree, name, distribution, xBins, xMax, xMin, yBins, yMin, yMax, fiducial, cuts):
    if (name.find('data')!=-1 or not fiducial):
        t.Draw("%s>>%s(%s,%s,%s,%s,%s,%s)"%(distribution,name,xBins,xMin,xMax,yBins,yMin,yMax),"(%s)"%cuts)
    elif fiducial:
        t.Draw("%s>>%s(%s,%s,%s,%s,%s,%s)"%(distribution,name,xBins,xMin,xMax,yBins,yMin,yMax),"mcEventWeight*((Sum$(LLP_R_coord>4)>0) && (Sum$(LLP_R_coord<300)>0) && (Sum$(abs(LLP_Eta)<2.1)>0) && (Sum$(abs(LLP_Z_coord)<800)>0))*(%s)"%cuts)
    h=gDirectory.Get(name).Clone(name)
    return h

def twoDplotAttr():
    return " tree <TTree>, name <str>, distribution <str>, xBins <int/array>, xMin <float/int>, xMax <float/int>, yBins <int/array>, yMin <float/int>, yMax <float/int>, fiducial <bool>, cuts <str> (set to 1 if there are no cuts), returns h <TH2>"

def integral(histList,dimensions):
    for hist in histList:
        if (dimensions==1):
            hist.Scale(1./hist.Integral(0,hist.GetNbinsX()+1))
        else:
            hist.Scale(1./hist.Integral(0,hist.GetNbinsX()+1,0,hist.GetNbinsY()+1))
        hist.Draw()

def integralAttr():
    return "histList list of <TH1/TH2>, dimensions <int>"

def maxMin(histList,staticMax,staticMin):
    maximum=0
    for hist in histList:
        if hist.GetMaximum()>maximum: maximum=hist.GetMaximum()

    for hist in histList:
        if staticMax==-999: hist.SetMaximum(maximum*1.1)
        else: hist.SetMaximum(staticMax)
        hist.SetMinimum(staticMin)

def maxMinAttr():
    return "histList list of <TH1/TH2>, staticMax <float/int> (set to -999 for auto calculation of maximum), staticMin <float/int>"

def profile(hist,name,fit,color):
    hist.Draw("colz")
    hist.ProfileX(name)
    name.SetLineColor(color)
    name.Draw("same")
    name.Fit(fit)

def profileAttr():
    return "hist <TH1>, name <str>, fit <str>, color <TColor name>"

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

def optimizeOneD(signal1,signal2,data,target):
    signal1.Scale(1./signal1.Integral(0,signal1.GetNbinsX()+1))
    signal2.Scale(1./signal2.Integral(0,signal2.GetNbinsX()+1))
    data.Scale(1./data.Integral(0,data.GetNbinsX()+1))

    nBins=signal1.GetNbinsX()+1

    maxX=0

    maximum=0
    for i in range(1,nBins):
        fom=(signal1.Integral(0,i)+signal2.Integral(0,i))/max(0.00000000001,data.Integral(0,i))
        if fom>maximum and (data.Integral(0,i))/data.Integral() > target:
            maximum=fom
            maxX=i

    print "signal1 "+str(signal1.Integral(0,maxX)),maxX
    print "signal2 "+str(signal2.Integral(0,maxX))
    print "data "+str(data.Integral(0,maxX))

def optimizeOneDAttr():
    return "signal1 <TH1>, signal2 <TH1>, data <TH1>, target <float>"

def optimizeTwoD(signal1,signal2,data,target):
    signal1.Scale(1./signal1.Integral(0,signal1.GetNbinsX()+1,0,signal1.GetNbinsY()+1))
    signal2.Scale(1./signal2.Integral(0,signal2.GetNbinsX()+1,0,signal2.GetNbinsY()+1))
    data.Scale(1./data.Integral(0,data.GetNbinsX()+1,0,data.GetNbinsY()+1))

    nBinsX=signal1.GetNbinsX()+1
    nBinsY=signal1.GetNbinsY()+1

    maxX=0
    maxY=0

    maximum=0
    for i in range(1,nBinsX):
        fom=(signal1.Integral(0,i,0,nBinsY)+signal1.Integral(0,nBinsX,0,nBinsY)+signal2.Integral(0,i,0,nBinsY)+signal2.Integral(0,nBinsX,0,nBinsY)-signal1.Integral(0,i,0,nBinsY)-signal2.Integral(0,i,0,nBinsY))/max(0.00000000001,data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,nBinsY)-data.Integral(0,i,0,nBinsY))
        if fom>maximum and (data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,nBinsY)-data.Integral(0,i,0,nBinsY))/data.Integral() > target:
            maximum=fom
            maxX=i
            maxY=nBinsY

    print "signal1 "+str(signal1.Integral(0,maxX,0,nBinsY)+signal1.Integral(0,nBinsX,0,maxY)-signal1.Integral(0,maxX,0,maxY)),signal1.Integral(0,maxX,0,nBinsY),signal1.Integral(0,nBinsX,0,maxY),maxX,maxY
    print "signal2 "+str(signal2.Integral(0,maxX,0,nBinsY)+signal2.Integral(0,nBinsX,0,maxY)-signal2.Integral(0,maxX,0,maxY)),signal2.Integral(0,maxX,0,nBinsY),signal2.Integral(0,nBinsX,0,maxY)
    print "data "+str(data.Integral(0,maxX,0,nBinsY)+data.Integral(0,nBinsX,0,maxY)-data.Integral(0,maxX,0,maxY)),data.Integral(0,maxX,0,nBinsY),data.Integral(0,nBinsX,0,maxY)

def optimizeTwoDAttr():
    return "signal1 <TH2>, signal2 <TH2>, data <TH2>, target <float>"
