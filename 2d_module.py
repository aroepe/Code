from ROOT import *
from plot import *

import os
import pdb

gROOT.SetBatch(1)

inputDir='${EOS2}/run349451/'

channels=[#['mH125_mass_15_lt_10','rel21_mc_mH125_a15_ct10.v10/data-tree/mc16_13TeV.309853.PowhegPythia8EvtGen_WpH_H125_a15a15_4b_ctau10.merge.AOD.e6753_e5984_a875_r10201_r10210_tid14407181_00.root'],
	  ['mH125_mass_15_lt_100','rel21_mc_mH125_a15_ct100.nominal/data-tree/mc16_13TeV.309855.PowhegPythia8EvtGen_WpH_H125_a15a15_4b_ctau100.merge.AOD.e6753_e5984_a875_r10201_r10210_tid14407201_00.root'],
	  ['mH125_mass_55_lt_10','rel21_mc_mH125_a55_ct10.nominal/data-tree/mc16_13TeV.309854.PowhegPythia8EvtGen_WpH_H125_a55a55_4b_ctau10.merge.AOD.e6753_e5984_a875_r10201_r10210_tid14407191_00.root'],
	  #['mH125_mass_55_lt_100','rel21_mc_mH125_a55_ct100.v10/data-tree/mc16_13TeV.309856.PowhegPythia8EvtGen_WpH_H125_a55a55_4b_ctau100.merge.AOD.e6753_e5984_a875_r10201_r10210_tid14407211_00.root'],
	  ['data','rel21_mc_data.nominal/data-tree/data18_13TeV.root']]

observables=[#'jet_EMFrac[0]:jet_FracSamplingMax[0]']
	     #'jet_alpha_max[0]:jet_chf[0]']
             #'min(jet_alpha_max[0],jet_alpha_max[1]):min(abs(jet_Timing[1]),abs(jet_Timing[0]))']
	     'min(jet_alpha_max[0],jet_alpha_max[1]):min(jet_chf[1],jet_chf[0])']#,'jet_ip_sig_med', 'jet_chf']
	

#lifetimes=['(LLP_R_coord[0]<10 && LLP_R_coord[1]<10)','(LLP_R_coord[0]>10 && LLP_R_coord[1]>10) && (LLP_R_coord[0]<30 && LLP_R_coord[1]<30)','(LLP_R_coord[0]>30 && LLP_R_coord[1]>30)']

gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)

##############################################################################

def weight(tau, new, lifetime):

	return 'exp(- %s / %s)/exp(- %s / %s)'%(lifetime,str(new),lifetime,str(tau))





##############################################################################



#######################################                                                                                                                                            

def getCuts3(s1,s2,data):
	target=1.0
	s1.Scale(1./s1.Integral(0,s1.GetNbinsX()+1,0,s1.GetNbinsY()+1))
	s2.Scale(1./s2.Integral(0,s2.GetNbinsX()+1,0,s2.GetNbinsY()+1))
	data.Scale(1./data.Integral(0,data.GetNbinsX()+1,0,data.GetNbinsY()+1))

	nBinsX=s1.GetNbinsX()+1
	nBinsY=s1.GetNbinsY()+1

	maxX=0
	maxY=0
	
	maximum=0
	for i in range(1,nBinsX):
		fom=(s1.Integral(0,i,0,nBinsY)+s1.Integral(0,nBinsX,0,nBinsY)+s2.Integral(0,i,0,nBinsY)+s2.Integral(0,nBinsX,0,nBinsY)-s1.Integral(0,i,0,nBinsY)-s2.Integral(0,i,0,nBinsY))/max(0.00000000001,data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,nBinsY)-data.Integral(0,i,0,nBinsY))
		if fom>maximum and (data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,nBinsY)-data.Integral(0,i,0,nBinsY))/data.Integral() > target:
			maximum=fom
			maxX=i
			maxY=nBinsY

	print "s1 "+str(s1.Integral(0,maxX,0,nBinsY)+s1.Integral(0,nBinsX,0,maxY)-s1.Integral(0,maxX,0,maxY)),s1.Integral(0,maxX,0,nBinsY),s1.Integral(0,nBinsX,0,maxY),maxX,maxY
	print "s2 "+str(s2.Integral(0,maxX,0,nBinsY)+s2.Integral(0,nBinsX,0,maxY)-s2.Integral(0,maxX,0,maxY)),s2.Integral(0,maxX,0,nBinsY),s2.Integral(0,nBinsX,0,maxY)
	print "data "+str(data.Integral(0,maxX,0,nBinsY)+data.Integral(0,nBinsX,0,maxY)-data.Integral(0,maxX,0,maxY)),data.Integral(0,maxX,0,nBinsY),data.Integral(0,nBinsX,0,maxY)

        print "s1 "+str(s1.Integral(0,30,0,nBinsY)+s1.Integral(0,nBinsX,0,10)-s1.Integral(0,30,0,10)),s1.Integral(0,30,0,nBinsY),s1.Integral(0,nBinsX,0,10),30,10
        print "s2 "+str(s2.Integral(0,30,0,nBinsY)+s2.Integral(0,nBinsX,0,10)-s2.Integral(0,30,0,10)),s2.Integral(0,30,0,nBinsY),s2.Integral(0,nBinsX,0,10)
        print "data "+str(data.Integral(0,30,0,nBinsY)+data.Integral(0,nBinsX,0,10)-data.Integral(0,30,0,10)),data.Integral(0,30,0,nBinsY),data.Integral(0,nBinsY,0,10)

	print "s1 "+str(s1.Integral(0,30,0,nBinsY)+s1.Integral(0,nBinsX,0,3)-s1.Integral(0,30,0,3)),s1.Integral(0,30,0,nBinsY),s1.Integral(0,nBinsX,0,3),30,3
	print "s2 "+str(s2.Integral(0,30,0,nBinsY)+s2.Integral(0,nBinsX,0,3)-s2.Integral(0,30,0,3)),s2.Integral(0,30,0,nBinsY),s2.Integral(0,nBinsX,0,3)
	print "data "+str(data.Integral(0,30,0,nBinsY)+data.Integral(0,nBinsX,0,3)-data.Integral(0,30,0,3)),data.Integral(0,30,0,nBinsY),data.Integral(0,nBinsY,0,3)

	target=0.5

	maximum=0
	for i in range(1,nBinsX):
		for j in range(1,nBinsY):
			fom=(s1.Integral(0,i,0,nBinsY)+s1.Integral(0,nBinsX,0,j)+s2.Integral(0,i,0,nBinsY)+s2.Integral(0,nBinsX,0,j)-s1.Integral(0,i,0,j)-s2.Integral(0,i,0,j))/max(0.00000000001,data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,j)-data.Integral(0,i,0,j))
			if fom>maximum and (data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,j)-data.Integral(0,i,0,j))/data.Integral() > target:
				maximum=fom
				maxX=i
				maxY=j

	print "s1 "+str(s1.Integral(0,maxX,0,nBinsY)+s1.Integral(0,nBinsX,0,maxY)-s1.Integral(0,maxX,0,maxY)),s1.Integral(0,maxX,0,nBinsY),s1.Integral(0,nBinsX,0,maxY),maxX,maxY
	print "s2 "+str(s2.Integral(0,maxX,0,nBinsY)+s2.Integral(0,nBinsX,0,maxY)-s2.Integral(0,maxX,0,maxY)),s2.Integral(0,maxX,0,nBinsY),s2.Integral(0,nBinsX,0,maxY)
	print "data "+str(data.Integral(0,maxX,0,nBinsY)+data.Integral(0,nBinsX,0,maxY)-data.Integral(0,maxX,0,maxY)),data.Integral(0,maxX,0,nBinsY),data.Integral(0,nBinsX,0,maxY)

	target=0.25

	maximum=0
	for i in range(1,nBinsX):
		for j in range(1,nBinsY):
			fom=(s1.Integral(0,i,0,nBinsY)+s1.Integral(0,nBinsX,0,j)+s2.Integral(0,i,0,nBinsY)+s2.Integral(0,nBinsX,0,j)-s1.Integral(0,i,0,j)-s2.Integral(0,i,0,j))/max(0.00000000001,data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,j)-data.Integral(0,i,0,j))
			if fom>maximum and (data.Integral(0,i,0,nBinsY)+data.Integral(0,nBinsX,0,j)-data.Integral(0,i,0,j))/data.Integral() > target:
				maximum=fom
				maxX=i
				maxY=j

	print "s1 "+str(s1.Integral(0,maxX,0,nBinsY)+s1.Integral(0,nBinsX,0,maxY)-s1.Integral(0,maxX,0,maxY)),s1.Integral(0,maxX,0,nBinsY),s1.Integral(0,nBinsX,0,maxY),maxX,maxY
	print "s2 "+str(s2.Integral(0,maxX,0,nBinsY)+s2.Integral(0,nBinsX,0,maxY)-s2.Integral(0,maxX,0,maxY)),s2.Integral(0,maxX,0,nBinsY),s2.Integral(0,nBinsX,0,maxY)
	print "data "+str(data.Integral(0,maxX,0,nBinsY)+data.Integral(0,nBinsX,0,maxY)-data.Integral(0,maxX,0,maxY)),data.Integral(0,maxX,0,nBinsY),data.Integral(0,nBinsX,0,maxY)

        print "s1 "+str(s1.Integral(0,10,0,nBinsY)+s1.Integral(0,nBinsX,0,1)-s1.Integral(0,10,0,1)),s1.Integral(0,10,0,nBinsY),s1.Integral(0,nBinsX,0,1),10,1
        print "s2 "+str(s2.Integral(0,10,0,nBinsY)+s2.Integral(0,nBinsX,0,1)-s2.Integral(0,10,0,1)),s2.Integral(0,10,0,nBinsY),s2.Integral(0,nBinsX,0,1)
        print "data "+str(data.Integral(0,10,0,nBinsY)+data.Integral(0,nBinsX,0,1)-data.Integral(0,10,0,1)),data.Integral(0,10,0,nBinsY),data.Integral(0,nBinsY,0,1)


	return s1.GetXaxis().GetBinLowEdge(maxX)+s1.GetXaxis().GetBinWidth(maxX),s1.GetYaxis().GetBinLowEdge(maxY)+s1.GetYaxis().GetBinWidth(maxY)


###################  

def plot(distribution, channel):

	files=[]
	for c in channel:
		files.append(filePath(inputDir,c))

        trees=[]
	histograms=[]
	filterhistograms=[]
        for f in files:
		t=f[0].Get('outTree')
		trees.append(t)
		hName='h%s'%(f[1])
		histograms.append(twoDplot(t,hName,distribution,100,0,1,100,0,1,True,"1"))
		filterhistograms.append(twoDplot(t,hName+"filter",distribution,100,0,1,100,0,1,True,"passesFilter && (Length$(passedTriggers)>0)"))

	#pdb.set_trace()

        c=TCanvas()

        for h,l,f in zip(histograms,channel,filterhistograms):
		h.SetTitle('Minimum Jet Alpha Max vs CHF for '+l[0].replace('mH125_','').replace('_',' ')+';CHF;Alpha Max')
                f.SetTitle('Minimum Jet Alpha Max vs CHF for '+l[0].replace('mH125_','').replace('_',' ')+';CHF;Alpha Max')
		#h.SetTitle('Minimum Alpha Max vs. CHF for '+l[0].replace('mH125_','').replace('_',' ')+';CHF;Alpha Max')

	#print getCuts3(histograms[0],histograms[1],histograms[2])

	print "No Filter Optimize"
	
	optimizeTwoD(histograms[0],histograms[1],histograms[2],0.33)

	print "With Filter Optimize"

        optimizeTwoD(filterhistograms[0],filterhistograms[1],filterhistograms[2],0.33)

        print "No Filter Efficiency"

        twoDefficiency(histograms[0],histograms[1],histograms[2],0.3,0.1)

        print "With Filter Efficiency"

        twoDefficiency(filterhistograms[0],filterhistograms[1],filterhistograms[2],0.3,0.1)


	i=0
        for h,l in reversed(zip(histograms,channel)):
		#corr=h.GetCorrelationFactor()
		#cut=h.Integral(0,101,97,101)/h.Integral(0,101,0,101)
		#h.Scale(1./h.Integral(0,101,0,101))
		h.RebinX(5)
		h.RebinY(5)
		h.SetMinimum(0.1)
		h.SetMaximum(300)
                h.Draw('colz pfc')
                #latex1 = TLatex(0.52, 0.89, "#font[72]{EMFrac > 0.96: }"+str(round(cut,3)*100)+"%")
                #latex1.SetNDC()
                #latex1.SetTextSize(0.03)
                #latex1.SetTextAlign(13)
                #latex1.Draw()
		c.SetLogz()
		#c.SetLogy()
		#c.SaveAs('AlphaMax_vs_Chf_timing_%s_v10_log.pdf'%(l[0]))


################################################################################

for o in observables:
	plot(o,channels)


