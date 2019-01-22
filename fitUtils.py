#!/usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import *
import tdrstyle
import numpy as np 
import re

r.gSystem.Load("./PDFs/HWWLVJRooPdfs_cxx.so")

# add gaussian constraint                                                                                                   
def addConstraint(iWorkspace,iVar,iMean,iSigma,iList):
    print '---- Adding gaussian constraint to %s'%iVar.GetName()
    lMean = r.RooRealVar("%s_mean"%iVar.GetName(),"%s_mean"%iVar.GetName(),iMean);
    lSigma = r.RooRealVar("%s_sigma"%iVar.GetName(),"%s_sigma"%iVar.GetName(),iSigma);
    lConstraint = r.RooGaussian("constraint_pdf_%s"%iVar.GetName(),"constraint_pdf_%s"%iVar.GetName(),iVar,lMean,lSigma)
    lConstraint.Print()
    iList.append(lConstraint.GetName())
    getattr(iWorkspace,"import")(lConstraint,r.RooFit.RecycleConflictNodes())

# make Pdf from model, probably should put parameters in dictionary
def makePdf(iWorkspace,iVar,iLabel,iModel,iMc=False):
    print '---- Making pdf for %s with model %s'%(iLabel,iModel) 
    lVar = iWorkspace.var(iVar);
    lModelPdf = None
    lTag = "%s_%s"%(iModel,iLabel)
    pVarHigh       = r.RooRealVar("%sHi_%s"%(lVar.GetName(),lTag),"%sHi_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    pVarHigh1      = r.RooRealVar("%sHi1_%s"%(lVar.GetName(),lTag),"%sHi1_%s"%(lVar.GetName(),lTag),0.5,0.,1.);
    lConstraints = []
    gaus_mean = 82.
    gaus_sigma = 8.

    if iModel == "Exp_mc":
        lC_Exp_mc        = r.RooRealVar("c_mc"+lTag,"c_mc"+lTag,-0.01,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_mc);
        
    if iModel == "Exp_data":
        lC_Exp_data      = r.RooRealVar("c_data"+lTag,"c_data"+lTag,-0.01,-2.,0.05)
        lModelPdf = r.RooExponential("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_Exp_data);
       
    if "DoubleCB" in iModel:
        lMean1_gaus      = r.RooRealVar("mean_"+lTag,"mean_"+lTag,gaus_mean,gaus_mean*0.8,gaus_mean*1.2)
        lSigma1_gaus     = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,gaus_sigma,gaus_sigma*0.5,gaus_sigma*2)
        lAlpha1          = r.RooRealVar("alpha1_"+lTag,"alpha1_"+lTag,1.1,-10.,10.)
        lSign1           = r.RooRealVar("sign1_"+lTag,"sign1_"+lTag,37.,0.,48.)
        lAlpha2          = r.RooRealVar("alpha2_"+lTag,"alpha2_"+lTag,1.,-10.,10.)
        lSign2           = r.RooRealVar("sign2_"+lTag,"sign2_"+lTag,14.,0.,46.)
        if "pass" in iLabel:
            lAlpha1          = r.RooRealVar("alpha1_"+lTag,"alpha1_"+lTag,1.1,-10.,10.)
            lSign1           = r.RooRealVar("sign1_"+lTag,"sign1_"+lTag,17.,0.,48.)
            lAlpha2          = r.RooRealVar("alpha2_"+lTag,"alpha2_"+lTag,2.,-10.,10.)
            lSign2           = r.RooRealVar("sign2_"+lTag,"sign2_"+lTag,1.6,0.,10.)
        lModelPdf        = r.RooDoubleCrystalBall("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lMean1_gaus,lSigma1_gaus,lAlpha1,lSign1,lAlpha2,lSign2)
 
    if "ErfExpGaus_st" in iModel:
        lC_ErfExp        = r.RooRealVar("c_"+lTag,"c"+lTag,-0.04,-0.2,0.0)
        lOffSet_ErfExp   = r.RooRealVar("offset_"+lTag,"offset_"+lTag,82.,75.,90.)
        lWidth_ErfExp    = r.RooRealVar("width_"+lTag,"width_"+lTag,30.,10.,300.)
        lMean_Gaus       = r.RooRealVar("mean_"+lTag,"mean_"+lTag,82.,75.,90.)
        lSigma_Gaus      = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        pGaus            = r.RooGaussian("gaus_"+lTag,"gaus_"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp          = r.RooErfExpPdf("erfExp_"+lTag,"erfExp_"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp,pGaus),r.RooArgList(pVarHigh));

    if iModel == "ErfExp_st_mc":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,82.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if iModel == "ErfExp_st_data":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-5,1.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,82.,50.,1000.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,50.,20.,1000.)
        lModelPdf = r.RooErfExpPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);

    if "GausExp_tt" in iModel:
        if iLabel == "tt_mc" or iLabel == "tt_data":
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,80.,75.,90.)
        elif "fakeW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,80.,75.,90.)
        elif "realW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,80.,75.,90.)
        lSigma_Gaus        = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,40.)
        lC_Exp             = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-1.,1.)
        pGaus              = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pExp               = r.RooExponential("exp_"+lTag,"exp_%s"+lTag,lVar,lC_Exp)
        lModelPdf          = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pPoly,pExp),r.RooArgList(pVarHigh,pVarHigh1),1)

    if  "GausErfExp_tt" in iModel:
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,20.,0.,300.)
        if "realW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,82.,75.,90.)
        elif "fakeW" in iLabel:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,82.,75.,90.)
        else:
            lMean_Gaus     = r.RooRealVar("mean_"+lTag,"mean_"+lTag,82.,75.,90.)
        a1             = r.RooRealVar("a1_"+lTag,"a1_"+lTag,0.0001,-0.01,0.01)
        a2             = r.RooRealVar("a2_"+lTag,"a2_"+lTag,0.00001,-0.0001,0.01)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,20.,0.,60.)
        lSigma_Gaus    = r.RooRealVar("sigma_"+lTag,"sigma_"+lTag,7.,0.,60.)
        pPoly          = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(a1,a2),1)
        pErfExp        = r.RooErfExpPdf("erfExp_"+lTag,"erfExp_%s"+lTag,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pPoly          = r.RooPolynomial("poly_"+lTag,"poly_%s"+lTag,lVar,r.RooArgList(a1),1)
        pGaus          = r.RooGaussian("gaus_"+lTag,"gaus_%s"+lTag,lVar,lMean_Gaus,lSigma_Gaus)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pGaus,pErfExp,pPoly),r.RooArgList(pVarHigh,pVarHigh),1);

    if iModel == "ErfExp_tt_mc" or iModel == "ErfExp_tt_data":
        lC_ErfExp      = r.RooRealVar("c_"+lTag,"c_"+lTag,-0.02,-10,10.)
        lOffSet_ErfExp = r.RooRealVar("offset_"+lTag,"offset_"+lTag,200.,0.,300.)
        lWidth_ErfExp  = r.RooRealVar("width_"+lTag,"width_"+lTag,60.,0.,100.)
        pErfExp        = r.RooErfExpPdf("erfExp_%s"%iLabel,"erfExp_%s"%iLabel,lVar,lC_ErfExp,lOffSet_ErfExp,lWidth_ErfExp);
        lModelPdf = r.RooAddPdf("model_pdf_%s"%iLabel,"model_pdf_%s"%iLabel,r.RooArgList(pErfExp),r.RooArgList(pVarHigh),1);

    getattr(iWorkspace,"import")(lModelPdf,r.RooFit.RecycleConflictNodes())
    return iWorkspace.pdf("model_pdf_%s"%iLabel),lConstraints

# return RooExtendPdf and Constraints
def makeModel(iWorkspace,iVar,iLabel,iModel):
    print '---- Making model'
    lNumber = r.RooRealVar("number_%s"%iLabel,"number_%s"%iLabel,500,0.,1e7); 
    lModelPdf,lConstraints = makePdf(iWorkspace,iVar,iLabel,iModel)
    lModel = r.RooExtendPdf("model_%s"%iLabel,"model_%s"%iLabel,lModelPdf,lNumber)
    getattr(iWorkspace,"import")(lNumber,r.RooFit.RecycleConflictNodes())
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    iWorkspace.pdf("model_%s"%iLabel).Print()
    return iWorkspace.pdf("model_%s"%iLabel),lConstraints

# fit single rooDataset to model
# return list of constraints
def fitSingleMC(iWorkspace,iVar,iLabel,iModel):
    print '---- Fitting for %s'%iLabel
    lVar   = iWorkspace.var(iVar);
    lData  = iWorkspace.data(iLabel+"_D");
    lModel,lConstraints = makeModel(iWorkspace,iVar,iLabel,iModel)

    pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE));
    getattr(iWorkspace,"import")(lModel,r.RooFit.RecycleConflictNodes())
    pRooFitResult.Print()
    getattr(iWorkspace,"import")(pRooFitResult)
    
    draw(lVar,lData,[lModel],pRooFitResult,iLabel+"_only")
    return lConstraints
