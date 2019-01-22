#!/usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import *
import tdrstyle
import numpy as np 
import re
tdrstyle.setTDRStyle()

fPlotDir="plotsWtag/"

def drawFrame(iFrame,iData,iFuncs,iLegend,iColor,isPoisson=False):
    if isPoisson:
        iData.plotOn(iFrame,r.RooFit.DataError(r.RooAbsData.Poisson),r.RooFit.XErrorSize(0))
    else:
        iData.plotOn(iFrame)
    for pFunc in iFuncs:
        pFunc.plotOn(iFrame,r.RooFit.LineColor(iColor),r.RooFit.LineStyle(iColor != 50+1))
        pFunc.plotOn(iFrame,r.RooFit.Name("ErfExp"),r.RooFit.Components("erfExp*"),r.RooFit.LineColor(4))
        pFunc.plotOn(iFrame,r.RooFit.Name("Gaus"),r.RooFit.Components("gaus*"),r.RooFit.LineColor(3))
        pFunc.plotOn(iFrame,r.RooFit.Name("Poly"),r.RooFit.Components("poly*"),r.RooFit.LineColor(40))
        pFunc.plotOn(iFrame,r.RooFit.Name("Poly"),r.RooFit.Components("exp*"),r.RooFit.LineColor(6))
        pFunc.plotOn(iFrame,r.RooFit.Name("DoubleCB"),r.RooFit.Components("*CB*"),r.RooFit.LineColor(8))
    print iFrame.Print("")

def draw(iVar,iData,iFuncs,iRooFitResult,iLabel="A",iBinWidth=5):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    lLegend = getLegend()
    drawFrame(lFrame,iData,iFuncs,lLegend,50)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.2)
    lFrame.GetYaxis().SetTitle(" Events / %.1f GeV"%iBinWidth);
    lFrame.GetXaxis().SetTitle("PUPPI Softdrop Jet Mass (GeV) ");
    lchisq,lndof = getChi2NDOF(iData,iFuncs,iRooFitResult,iVar)
    addInfo = getPavetext()
    addInfo.AddText("#chi^{2}/nDOF = %.3f/%i"%(lchisq,lndof))
    print "chi^{2}/nDOF = %.3f/%i"%(lchisq,lndof)
    lFrame.Draw()
    addInfo.Draw()
    lCan.Modified()
    lCan.Update()
    lCan.SaveAs(fPlotDir+iLabel+".pdf")
    lCan.SaveAs(fPlotDir+iLabel+".png")

def getLegend():
    legend = TLegend(0.8,0.75,0.95,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.040)
    legend.SetTextAlign(12)
    return legend

def getPavetext():
    addInfo = TPaveText(0.2,0.75,0.5,0.9,"NDC")
    addInfo.SetFillColor(1)
    addInfo.SetLineColor(1)
    addInfo.SetFillStyle(1)
    addInfo.SetBorderSize(0)
    addInfo.SetTextFont(42)
    addInfo.SetTextSize(0.040)
    addInfo.SetTextAlign(12)
    return addInfo

def drawDataMc(iVar,iData,iFuncsData,iMc,iFuncsMc,iRooFitResult_mc,iRooFitResult_data,params,iLabel='A',iBinWidth=5):
    lCan   = r.TCanvas(str(iLabel),str(iLabel),800,600)
    lFrame = iVar.frame()
    lLegend = getLegend()
    drawFrame(lFrame,iData,iFuncsData,lLegend,4,True)
    data_chisq,data_ndof = getChi2NDOF(iData,iFuncsData,iRooFitResult_data,iVar)
    mc_chisq,mc_ndof     = getChi2NDOF(iMc,iFuncsMc,iRooFitResult_mc,iVar)
    data_chi2Prob        = r.TMath.Prob(data_chisq,data_ndof)
    mc_chi2Prob          = r.TMath.Prob(mc_chisq,mc_ndof)
    drawFrame(lFrame,iMc,iFuncsMc,lLegend,2)
    lFrame.GetYaxis().SetRangeUser(0,lFrame.GetMaximum()*1.3)
    lFrame.GetYaxis().SetTitle(" Events / %.1f GeV"%iBinWidth);
    lFrame.GetXaxis().SetTitle( "PUPPI Softdrop Jet Mass (GeV)")
    lLegend.AddEntry(lFrame.getCurve("model_mc_Norm[Puppijet0_msd]"),"MC","l")
    lLegend.AddEntry(lFrame.getCurve("model_data_Norm[Puppijet0_msd]"),"Data","l")
    addInfo = getPavetext()
    addInfo.AddText("Data #chi^{2}/nDOF = %.1f/%i, p = %.2f"%(data_chisq,data_ndof,data_chi2Prob))
    addInfo.AddText("MC #chi^{2}/nDOF = %.1f/%i, p = %.2f"%(mc_chisq,mc_ndof,mc_chi2Prob))
    lFrame.Draw()
    addInfo.Draw()
    lLegend.Draw()
    lCan.Modified()
    lCan.Update()
    lCan.SaveAs(fPlotDir+iLabel+".pdf")
    lCan.SaveAs(fPlotDir+iLabel+".png")

def getChi2NDOF(iData,iModel,iRooFitResult,iVar):
    lDataHist = iData.binnedClone(iData.GetName()+"_binnedClone",iData.GetName()+"_binnedClone");
    pChi2     = iModel[0].createChi2(lDataHist,r.RooFit.Range(iVar.getMin(),iVar.getMax()),r.RooFit.Extended(r.kTRUE),r.RooFit.DataError(r.RooAbsData.Poisson))
    return pChi2.getVal(), int(iVar.getBins() - iRooFitResult.floatParsFinal().getSize())

def getPull(iVar,iPlot):
    lHPull = iPlot.pullHist();
    pPull  = iVar.frame(r.RooFit.Title("Pull"), r.RooFit.Bins(int(iVar.getBins())));
    lLine  = r.TLine(iVar.getMin(),0,iVar.getMax(),0);
    lLine.SetLineWidth(2); lLine.SetLineColor(r.kRed);
    pPull.addObject(lLine);
    pPull.addPlotable(lHPull,"P");
    pPull.SetTitle("");
    pPull.GetXaxis().SetTitle("");
    pPull.GetYaxis().SetTitle("pull");
    pPull.GetYaxis().SetRangeUser(-4,4);
    pPull.GetXaxis().SetTitleSize(0.10);
    pPull.GetXaxis().SetLabelSize(0.10);
    pPull.GetYaxis().SetTitleSize(0.10);
    pPull.GetYaxis().SetLabelSize(0.10);
    return pPull
