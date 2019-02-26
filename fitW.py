#!/usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import *
import numpy as np 
import re
from fitUtils import *
from plotUtils import *

fVar="Puppijet0_msd"
fPt="Puppijet0_pt"
fPt_bins = [200,100000]
fCats = ['pass','fail']
fWCut = {'pass': "Puppijet1_lsfC_3>0.7",
         'fail': "Puppijet1_lsfC_3<=0.7"}
fWeight="weight"
fLumi=35.9
fBinWidth =5
fDir="/afs/cern.ch/user/c/cmantill/public/forSangEon/sklimWtag"
fOutput="w.root"

# files and names
fFiles = {'data': 'Data',
          'wlnu_mc' : 'WJets',
          'wlnu_data': 'WJets',
          'st_mc':   'ST',
          'st_data':   'ST',
#          'tt_fakeW_mc': 'TT',
          'tt_realW_mc': 'TT',
#          'tt_fakeW_data': 'TT',
          'tt_realW_data': 'TT',
          "tt_signal_mc": 'TT',
          "tt_signal_data": 'TT',
          "tt_bkg_mc": 'TT',
          "tt_bkg_data": 'TT',
          'mc': 'PseudoData'
          }

# pdfs
fPDFs = {
         'wlnu_mc':       'Exp_mc',
         'wlnu_data':     'Exp_data',
         'st_mc':         'ErfExpGaus_st_mc',
         'st_data':       'ErfExpGaus_st_data',
#         'tt_fakeW_mc':   'ErfExp_tt_mc',
         'tt_realW_mc':   'DoubleCB',
#         'tt_fakeW_data': 'ErfExp_tt_data',
         'tt_realW_data': 'DoubleCB'}
fPDFs['tt_signal_mc'] = fPDFs['tt_realW_mc']
fPDFs['tt_signal_data'] = fPDFs['tt_realW_data']
#fPDFs['tt_bkg_mc'] = fPDFs['tt_fakeW_mc']
#fPDFs['tt_bkg_data'] = fPDFs['tt_fakeW_data']

# floating params
fFloating = {}
fFloating['mc'] = ['mean_DoubleCB_tt_realW_mc_',
                   'sigma_DoubleCB_tt_realW_mc_',
                   'alpha1_DoubleCB_tt_realW_mc_',
                   'alpha2_DoubleCB_tt_realW_mc_',
                   'sign1_DoubleCB_tt_realW_mc_',
                   'sign2_DoubleCB_tt_realW_mc_',
                   'number_st_mc_',
                   'number_wlnu_mc_',
                   'number_tt_realW_mc_',
                   'eff_tt_realW_mc_',
                   'numbertotal_tt_realW_mc_',
                   ]
fFloating['data'] = ['mean_DoubleCB_tt_realW_data_',
                     'sigma_DoubleCB_tt_realW_data_',
                     'alpha1_DoubleCB_tt_realW_data_',
                     'alpha2_DoubleCB_tt_realW_data_',
                     'sign1_DoubleCB_tt_realW_data_',
                     'sign2_DoubleCB_tt_realW_data_',
                     'number_st_data_',
                     'number_wlnu_data_',
                     'number_tt_realW_data_',
                     'eff_tt_realW_data_',
                     'numbertotal_tt_realW_data_',
                     ]

r.gSystem.Load("./PDFs/HWWLVJRooPdfs_cxx.so")

def parser():
    parser = OptionParser()
    parser.add_option('--tag',     action='store', type='string', dest='tag',          default='AK8v12017WP026',  help='samples tag')
    parser.add_option('--xMin',    action='store', type='float',  dest='xmin',         default=50,                help='x-min')
    parser.add_option('--xMax',    action='store', type='float',  dest='xmax',         default=130,               help='x-max')
    parser.add_option('-b',        action='store_true',           dest='noX',          default=True,              help='no X11 windows')
    parser.add_option('--ind',     action='store_true',           dest='individual',   default=False,             help='individual pass and fail fit')
    parser.add_option('--comb',    action='store_true',           dest='combined',     default=False,             help='combine data and MC fit')
    parser.add_option('--fits',    action='store_true',           dest='fits',         default=True,              help='fits to each sample')
    (options,args) = parser.parse_args()
    return options

# main class
class WPeak():
    def __init__(self,options):
        self._lW      = r.RooWorkspace("w","w")

        self._lMSD_lo = fXMin
        self._lMSD_hi = fXMax
        self._lMSD    = r.RooRealVar(fVar,fVar,self._lMSD_lo,self._lMSD_hi)
        self._lMSD.setBins(fNBins)
        getattr(self._lW,"import")(self._lMSD)
        self._lPt     = r.RooRealVar(fPt,fPt,0,2000)
        self._lWeight = r.RooRealVar(fWeight,fWeight,-1e+25,1e+25)
        self._bFit    = options.fits
        
        self._lPDatas = {}
        self._lHPdfs  = {}
        self._lNPdfs  = {}
        self._lModels = {}
        self._lConstraints = {}
        self._lScaleNumber = {}
        self._TTsf = {};

        # get dataset for individual processes except for mc and tt_bkg(fakew), tt_signal(realW)
        for iLabel,iFile in fFiles.iteritems():
            if iLabel == "mc" or  "bkg" in iLabel or "signal" in iLabel: continue
            for Pt_it in range(len(fPt_bins)-1):
                pFile = fDir+fTag+"/"+iFile+"_"+fTag+".root"
                # pt cut and label
                fCut = "(Puppijet0_pt>" + str(fPt_bins[Pt_it]) + " && Puppijet0_pt<" + str(fPt_bins[Pt_it+1]) + ")"
                Pt_label = "pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
                # build dataset for each pt category
                for Cat_it in fCats:
                    self._TTsf[Pt_label+ '_' + Cat_it] = 1
                    print 'preparing sample %s in pT bin[%s,%s] and %s category'%(iLabel,fPt_bins[Pt_it],fPt_bins[Pt_it+1],Cat_it)
                    pName = iLabel + '_' + Pt_label + '_' + Cat_it
                    pCut = "(%s&&%s)"%(fCut,fWCut[Cat_it])
                    if 'wlnu' in iLabel: pCut = "(%s&&weight<0.2)"%pCut
                    if iLabel != "data":
                        self._lPDatas[pName] = self.getDataset(iLabel,pName,Pt_label+ '_' + Cat_it,pFile,pCut,True)
                    else:
                        self._lPDatas[pName] = self.getDataset(iLabel,pName,Pt_label+ '_' + Cat_it,pFile,pCut,False)
                    # fit individual mc pdfs (tt,st,wjets?)
                    if iLabel != "data" and self._bFit:
                        self._lConstraints[pName] = self.fitSingleMC(pName,fPDFs[iLabel])

        # get normalization SF for tt for pass and fail cats
        # sf_(cat)_pTbin = (data-st-wjets)/tt with cat=[pass,fail]
        for Pt_it in range(len(fPt_bins)-1):
            Pt_label = "pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
            for Cat_it in fCats:
                self._TTsf[Pt_label+'_'+Cat_it] = self.getTTSF(Pt_label,Cat_it)
        
        # get dataset for mc pass and fail and scale it by tt_sf norm
        iLabel = "mc"
        for Pt_it in range(len(fPt_bins)-1):
            pFile = fDir+fTag+"/"+fFiles[iLabel]+"_"+fTag+".root"
            Pt_label = "pT_" + str(fPt_bins[Pt_it]) + '_' + str(fPt_bins[Pt_it+1])
            for Cat_it in fCats:
                print 'preparing sample %s in pT bin[%s,%s] and %s category'%(iLabel,fPt_bins[Pt_it],fPt_bins[Pt_it+1],Cat_it)
                pName = iLabel + '_' + Pt_label + '_' + Cat_it
                pCut = "(%s&&%s)"%(fCut,fWCut[Cat_it])
                self._lPDatas[pName] = self.getDataset(iLabel,pName,Pt_label+ '_' + Cat_it,pFile,pCut,True)
                
        # print workspace
        self._lW.Print()
        
    # get dataset for each process from tree
    def getDataset(self,iLabel,iName,Pt_label,iFile,iCut="(1==1)",iWeight=False,iTree="otree2"):
        print 'preparing dataset for ',iFile,' with cut ',iCut, ' and TTsf ',self._TTsf[Pt_label]
        lFile   = r.TFile(iFile)
        lTree   = lFile.Get(iTree)
        lData = r.RooDataSet(iName+"_D",iName+"_D",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight))
        lCut = r.TTreeFormula("lCut",iCut,lTree)  
        lTree.SetNotify(lCut)       
        for i0 in range(lTree.GetEntriesFast()):
            lTree.LoadTree(i0) 
            lSel = False 
            for i1 in range(lCut.GetNdata()):  
                if (lCut.EvalInstance(i1)):  
                    lSel=True; break;
            if not lSel: continue     
            lTree.GetEntry(i0) 
            if iWeight:
                lWeight = getattr(lTree,fWeight)*fLumi*self._TTsf[Pt_label]
                #if 'wlnu' in iLabel: 
                #    print getattr(lTree,fWeight)
            else:
                lWeight = 1
            lMass = getattr(lTree,fVar)
            lMatched = 0
            jmatched = getattr(lTree,"Puppijet0_vMatching");
            jhadronic = getattr(lTree,"Puppijet0_isHadronicV");
            if jhadronic == 1.0 and jmatched < 0.8 and jmatched > 0.:
                lMatched = 1;
#####################################################
### REMOVED until we have GEN lepton info in bits ###
#####################################################
#            if 'realW' in iLabel and lMatched == 0: continue
#            if 'fakeW' in iLabel and lMatched == 1: continue
            if lMass < self._lMSD_hi and lMass > self._lMSD_lo:
                self._lMSD.setVal(lMass)
                lData.add(r.RooArgSet(self._lMSD), lWeight)
        getattr(self._lW,'import')(lData,r.RooFit.RecycleConflictNodes())
        lFile.Close()
        return self._lW.data(iName+"_D")
                     
    # fit single mc
    def fitSingleMC(self,iLabel,iModel):
        print '---- Fitting for %s'%iLabel
        lVar   = self._lW.var(fVar);
        lData  = self._lW.data(iLabel+"_D");
        lModel,lConstraints = makeModel(self._lW,fVar,iLabel,iModel)
        pRooFitResult = lModel.fitTo(lData,r.RooFit.Save(1),r.RooFit.SumW2Error(r.kTRUE),r.RooFit.Extended(r.kTRUE),r.RooFit.Minimizer("Minuit2"),r.RooFit.Verbose(r.kFALSE));
        pRooFitResult.Print()
        getattr(self._lW,"import")(lModel,r.RooFit.RecycleConflictNodes())
        getattr(self._lW,"import")(pRooFitResult)
        draw(lVar,lData,[lModel],pRooFitResult,iLabel+"_only",fBinWidth)
        return lConstraints

    # getModel
    def getModel(self, iLabel):
        print '---- Getting model for ', iLabel
        lModel    = self._lW.pdf("model_"+iLabel)
        lModel.Print()
        return lModel

    # return roohistpdf
    def histPdf(self, iLabel):
        lDataSet = self._lW.data(iLabel+"_D")
        lVar = self._lW.var(fVar)
        lReduced  = lDataSet.reduce(r.RooArgSet(lVar))
        lDataHist = lReduced.binnedClone(lReduced.GetName()+"_binnedClone",lReduced.GetName()+"_binnedClone");
        lHPdf     = r.RooHistPdf(lReduced.GetName()+"P",lReduced.GetName()+"P",r.RooArgSet(lVar),lDataHist,0)
        print lReduced.sumEntries()
        getattr(self._lW,"import")(lHPdf,r.RooFit.RecycleConflictNodes())
        return lHPdf

    # get TT normalization SF for pt label
    def getTTSF(self,iPtLabel,Cat_it):
        print '---- Getting TT SF for ', iPtLabel, Cat_it
        nData = nTT = nMinorBKG = 0;
        for Data_it in self._lPDatas:
            if Cat_it not in Data_it: continue
            if iPtLabel not in Data_it: continue
            nEntries = self._lW.data(Data_it+"_D").sumEntries()
            print Data_it,nEntries
            if Data_it == "data_"+iPtLabel+"_"+Cat_it:
                nData += nEntries
            if "_mc" in Data_it:
                if "realW" in Data_it or "fakeW" in Data_it:
                    nTT += nEntries
                else:
                    nMinorBKG += nEntries
        ttScalefactor = (nData-nMinorBKG)/nTT
        #ttScalefactor = 1
        scalefactor = r.RooRealVar("tt_scalefactor_"+Cat_it,"tt_scalefactor_"+Cat_it,ttScalefactor)
        getattr(self._lW,'import')(scalefactor,r.RooFit.RecycleConflictNodes())
        return scalefactor.getVal()

    # get Wtag SFs
    def getWtagSFs(self,iPtLabel,iModelLabel,iTTLabel):
        # efficiency: realW tt evts:  pass / pass+fail
        # N2 sf = efficiency in data / efficiency in mc 
        # mass sf = mean in data / mean in mc
        # mass shift = mean in data - mean in mc
        # mass smear = sigma in data / sigma in mc
        # sigma sf = sigma in data / sigma in mc
        pEff = {}; pMean = {}; pSigma = {}
        pEffErr = {}; pMeanErr = {}; pSigmaErr = {}

        for i0 in ['mc','data']:
            lTTLabel = iModelLabel.replace('mc',i0)+"_"+iTTLabel.replace('mc',i0)
            if "signal" in lTTLabel:
                pEff[i0] = self._lW.var("eff_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel).getVal(); 
                pEffErr[i0] = self._lW.var("eff_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel).getError();
            if "realW" in lTTLabel:
                pEff[i0] = self._lW.var("eff_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel).getVal();
                pEffErr[i0] = self._lW.var("eff_"+iTTLabel.replace('mc',i0)+"_"+iPtLabel).getError();

            print "mean_"+lTTLabel+"_"+iPtLabel+"_pass"
            pMean[i0] = self._lW.var("mean_"+lTTLabel+"_"+iPtLabel+"_pass").getVal();
            pMeanErr[i0] = self._lW.var("mean_"+lTTLabel+"_"+iPtLabel+"_pass").getError();
            
            pSigma[i0] = self._lW.var("sigma_"+lTTLabel+"_"+iPtLabel+"_pass").getVal(); 
            pSigmaErr[i0] = self._lW.var("sigma_"+lTTLabel+"_"+iPtLabel+"_pass").getError(); 

        pEffSf = pEff['data']/pEff['mc'];
        pEffSfErr = pEffSf * math.sqrt( (pEffErr['data']/pEff['data'])**2 + (pEffErr['mc']/pEff['mc'])**2 )

        pMeanShift = pMean['data'] - pMean['mc']
        pMeanShiftErr = math.sqrt(pMeanErr['data']**2 +pMeanErr['mc']**2)

        pSigmaSmear = pSigma['data']/pSigma['mc']
        pSigmaSmearErr = pSigmaSmear * math.sqrt( (pSigmaErr['data']/pSigma['data'])**2 + (pSigmaErr['mc']/pSigma['mc'])**2 )

        pMeanSf = pMean['data']/pMean['mc']
        pMeanSfErr = pMeanSf * math.sqrt( (pMeanErr['data']/pMean['data'])**2 + (pMeanErr['mc']/pMean['mc'])**2 )

        pSigmaSf = pSigmaSmear
        pSigmaSfErr =  pSigmaSmearErr

        print 'W-tagging SF: %.3f +/- %.3f '%(pEffSf,pEffSfErr)
        print 'Mass shift [GeV]: %.3f +/- %.3f '%(pMeanShift,pMeanShiftErr)
        print 'Mass SF: %.3f +/- %.3f '%(pMeanSf,pMeanSfErr)
        print 'Resolution SF: %.3f +/- %.3f '%(pSigmaSf,pSigmaSfErr)
        print '\n'
        print 'eff data: %.3f +/- %.3f, mc: %.3f +/- %.3f'%(pEff['data'],pEffErr['data'],pEff['mc'],pEffErr['mc'])
        print '<m> data: %.3f +/- %.3f, mc: %.3f +/- %.3f'%(pMean['data'],pMeanErr['data'],pMean['mc'],pMeanErr['mc'])
        print 'res data: %.3f +/- %.3f, mc: %.3f +/- %.3f'%(pSigma['data'],pSigmaErr['data'],pSigma['mc'],pSigmaErr['mc'])

    # fix parameters
    def fixParams(self):
        lFloat = []; lConst = [];
        args = self._lW.allVars()
        args.Print()
        iter = args.createIterator()
        var = iter.Next()
        lFloat = []
        lConst = []
        while var:
            if any(f in var.GetName() for f in fFloating['mc']+fFloating['data']):
                lFloat.append(var.GetName())
                pass
            else:
                lConst.append(var.GetName())
                var.setConstant(r.kTRUE)
            var = iter.Next()
        return lFloat,lConst

    # fit data to mc model
    def fit(self):
        lVar   = self._lW.var(fVar);
        for Pt_it in range(len(fPt_bins)-1):
            pPtLabel = "pT_" + str(fPt_bins[Pt_it]) + "_" + str(fPt_bins[Pt_it+1])            

            # define categories
            pCats  = r.RooCategory("pCats"  ,"pCats")
            for Cat_it in fCats:
                pCats.defineType(Cat_it)

            # link pass and fail categories for realW
            # tt_realw_pass = eff*(tt_realW_pass+tt_realW_fail) 
            # tt_realw_fail = (1-eff)*(tt_realW_pass+tt_realW_fail) 
            nrealW_pass = self._lW.var('number_tt_realW_mc_'+pPtLabel+"_pass").getVal()
            nrealW_fail = self._lW.var('number_tt_realW_mc_'+pPtLabel+"_fail").getVal()
            effS      = nrealW_pass/(nrealW_pass+nrealW_fail)
            print "Real MC signal efficiency = " ,effS; print "";

            # perform fit to data and mc to that category
            for Cat_it in fCats:
                print 'processing simult models for pT bin[%s,%s] and %s category'%(fPt_bins[Pt_it],fPt_bins[Pt_it+1],Cat_it)
                pPtLabelCat = pPtLabel + "_"+Cat_it

                # define normalizations and print ratio of tt/minorBKG
                # so that tt_fakeW_mc = ratio * (st_mc+wlnu_mc) - tt_realW_mc
                norm_fakeW = self._lW.arg("number_tt_fakeW_mc_" + pPtLabelCat).getVal()
                norm_realW = self._lW.arg("number_tt_realW_mc_" + pPtLabelCat).getVal()
                norm_st    = self._lW.arg("number_st_mc_" + pPtLabelCat).getVal()
                norm_wlnu  = self._lW.arg("number_wlnu_mc_" + pPtLabelCat).getVal()
                norm_ratio = (norm_fakeW + norm_realW) / (norm_st + norm_wlnu)
                print 'norm ratio TT/MINORBKG: %s'%norm_ratio

                # get normalizations vars from model
                st_norm_mc = self._lW.arg("number_st_mc_"+pPtLabelCat)
                wlnu_norm_mc = self._lW.arg("number_wlnu_mc_"+pPtLabelCat)
                norm_ratio_mc = r.RooRealVar("norm_ratio_mc_"+pPtLabelCat,"norm_ratio_mc_"+pPtLabelCat,norm_ratio)
                st_norm_data = self._lW.arg("number_st_data_"+pPtLabelCat)
                wlnu_norm_data = self._lW.arg("number_wlnu_data_"+pPtLabelCat)
                norm_ratio_data = r.RooRealVar("norm_ratio_data_"+pPtLabelCat,"norm_ratio_data_"+pPtLabelCat,norm_ratio)

                # define normalization
                # for pass:
                #  tt_realW = eff*tt_realW_mc
                #  tt_fakeW = ratio * (st_mc+wlnu_mc) - tt_realW_mc =  ratio * (st_mc+wlnu_mc) - eff*tt_realW_mc
                # for fail:
                #  tt realW = (1-eff)*tt_realW_mc
                #  tt_fakeW = ratio * (st_mc+wlnu_mc) - tt_realW_mc =  ratio * (st_mc+wlnu_mc) - eff*tt_realW_mc 
                # need to define everything twice 1.for mc model, 2.for data model
                if Cat_it == "pass":
                    print "defining realW"
                    tt_realW_total_mc = r.RooRealVar("numbertotal_tt_realW_mc_%s"%pPtLabel,"numbertotal_tt_realW_mc_%s"%pPtLabel,500,0.,1e7);
                    tt_realW_eff_mc = r.RooRealVar("eff_tt_realW_mc_%s"%pPtLabel,"eff_tt_realW_mc_%s"%pPtLabel,effS,effS*0.8,effS*1.2);
                    tt_realW_norm_mc = r.RooFormulaVar("number_tt_realW_mc_%s"%(pPtLabelCat), "@0*@1", r.RooArgList(tt_realW_eff_mc,tt_realW_total_mc));
                    tt_fakeW_norm_mc = r.RooFormulaVar("tt_fakeW_norm_mc_"+pPtLabelCat,"(@0*(@1+@2)-(@3*@4))",r.RooArgList(norm_ratio_mc,wlnu_norm_mc,st_norm_mc,tt_realW_eff_mc,tt_realW_total_mc))

                    print "defining realW for data"
                    tt_realW_total_data = r.RooRealVar("numbertotal_tt_realW_data_%s"%pPtLabel,"numbertotal_tt_realW_data_%s"%pPtLabel,500,0.,1e7);
                    tt_realW_eff_data = r.RooRealVar("eff_tt_realW_data_%s"%pPtLabel,"eff_tt_realW_data_%s"%pPtLabel,effS,effS*0.8,effS*1.2);
                    tt_realW_norm_data = r.RooFormulaVar("number_tt_realW_data_%s"%(pPtLabelCat), "@0*@1", r.RooArgList(tt_realW_eff_data,tt_realW_total_data));
                    tt_fakeW_norm_data = r.RooFormulaVar("tt_fakeW_norm_data_"+pPtLabelCat,"(@0*(@1+@2)-(@3*@4))",r.RooArgList(norm_ratio_data,wlnu_norm_data,st_norm_data,tt_realW_eff_data,tt_realW_total_data))
                else:
                    print "defining realW"
                    tt_realW_total_mc = self._lW.var("numbertotal_tt_realW_mc_%s"%pPtLabel)
                    tt_realW_eff_mc = self._lW.var("eff_tt_realW_mc_%s"%pPtLabel)
                    tt_realW_norm_mc = r.RooFormulaVar("number_tt_realW_mc_%s"%(pPtLabelCat), "(1-@0)*@1", r.RooArgList(tt_realW_eff_mc,tt_realW_total_mc));
                    tt_fakeW_norm_mc = r.RooFormulaVar("tt_fakeW_norm_mc_"+pPtLabelCat,"(@0*(@1+@2)-((1-@3)*@4))",r.RooArgList(norm_ratio_mc,wlnu_norm_mc,st_norm_mc,tt_realW_eff_mc,tt_realW_total_mc))

                    print "defining realW for data"
                    tt_realW_total_data = self._lW.var("numbertotal_tt_realW_data_%s"%pPtLabel)
                    tt_realW_eff_data = self._lW.var("eff_tt_realW_data_%s"%pPtLabel)
                    tt_realW_norm_data = r.RooFormulaVar("number_tt_realW_data_%s"%(pPtLabelCat), "(1-@0)*@1", r.RooArgList(tt_realW_eff_data,tt_realW_total_data));
                    tt_fakeW_norm_data = r.RooFormulaVar("tt_fakeW_norm_data_"+pPtLabelCat,"(@0*(@1+@2)-((1-@3)*@4))",r.RooArgList(norm_ratio_data,wlnu_norm_data,st_norm_data,tt_realW_eff_data,tt_realW_total_data))

                # tt = ratio *bkg = tt_realW + fakeW                                                                                                                                                          
                # tt_realW_pass/fail = eff*(tt_realW)                                                                                                                                                          
                # tt_fakeW = ratio *bkg - tt_realW                                                                                                                                                             
                # to mc 
                tt_realW_model_mc_ext = r.RooExtendPdf("ext_model_mc_total_tt_realW_"+pPtLabelCat,"ext_model_mc_total_tt_realW_"+pPtLabelCat,self._lW.pdf("model_pdf_tt_realW_mc_"+pPtLabelCat),tt_realW_norm_mc)
                tt_fakeW_model_mc_ext = r.RooExtendPdf("ext_model_mc_total_tt_fakeW_"+pPtLabelCat,"ext_model_mc_total_tt_fakeW_"+pPtLabelCat,self._lW.pdf("model_pdf_tt_fakeW_mc_"+pPtLabelCat),tt_fakeW_norm_mc)
                st_model_mc_ext   = r.RooExtendPdf("ext_model_mc_total_st_"+pPtLabelCat,"ext_model_mc_total_st_"+pPtLabelCat,self._lW.pdf("model_pdf_st_mc_"+pPtLabelCat),st_norm_mc)
                wlnu_model_mc_ext = r.RooExtendPdf("ext_model_mc_total_wlnu_"+pPtLabelCat,"ext_model_mc_total_wlnu_"+pPtLabelCat,self._lW.pdf("model_pdf_wlnu_mc_"+pPtLabelCat),wlnu_norm_mc)

                tt_realW_model_data_ext = r.RooExtendPdf("ext_model_data_total_tt_realW_"+pPtLabelCat,"ext_model_data_total_tt_realW_"+pPtLabelCat,self._lW.pdf("model_pdf_tt_realW_data_"+pPtLabelCat),tt_realW_norm_data)
                tt_fakeW_model_data_ext = r.RooExtendPdf("ext_model_data_total_tt_fakeW_"+pPtLabelCat,"ext_model_data_total_tt_fakeW_"+pPtLabelCat,self._lW.pdf("model_pdf_tt_fakeW_data_"+pPtLabelCat),tt_fakeW_norm_data)
                st_model_data_ext   = r.RooExtendPdf("ext_model_data_total_st_"+pPtLabelCat,"ext_model_data_total_st_"+pPtLabelCat,self._lW.pdf("model_pdf_st_data_"+pPtLabelCat),st_norm_data)
                wlnu_model_data_ext = r.RooExtendPdf("ext_model_data_total_wlnu_"+pPtLabelCat,"ext_model_data_total_wlnu_"+pPtLabelCat,self._lW.pdf("model_pdf_wlnu_data_"+pPtLabelCat),wlnu_norm_data)

                lFloat,lConst =self.fixParams()

                self._lModels['TotalMc_'+pPtLabelCat] = r.RooAddPdf(("model_total_mc_"+pPtLabelCat),("model_total_mc_"+pPtLabelCat),
                                                                    r.RooArgList(tt_realW_model_mc_ext,tt_fakeW_model_mc_ext,st_model_mc_ext,wlnu_model_mc_ext))

                self._lModels['TotalData_'+pPtLabelCat] = r.RooAddPdf(("model_total_data_"+pPtLabelCat),("model_total_data_"+pPtLabelCat),
                                                                    r.RooArgList(tt_realW_model_data_ext,tt_fakeW_model_data_ext,st_model_data_ext,wlnu_model_data_ext))

                self._lModels['TotalData_'+pPtLabelCat].Print()
                getattr(self._lW,"import")(self._lModels["TotalMc_"+pPtLabelCat],r.RooFit.RecycleConflictNodes())
                getattr(self._lW,"import")(self._lModels["TotalData_"+pPtLabelCat],r.RooFit.RecycleConflictNodes())

                print 'Floating ',lFloat
                print 'Constant ',lConst

            # combined data (pass and fail)
            combData_data = r.RooDataSet("combData_data","combData_data",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight),
                                         RooFit.Index(pCats),
                                         RooFit.Import("data_"+pPtLabel+"_pass",self._lPDatas["data_"+pPtLabel+"_pass"]),
                                         RooFit.Import("data_"+pPtLabel+"_fail",self._lPDatas["data_"+pPtLabel+"_fail"]))
            # combined mc (pass and fail)
            combData_mc = r.RooDataSet("combData_mc","combData_mc",r.RooArgSet(self._lMSD,self._lWeight),r.RooFit.WeightVar(self._lWeight),
                                       RooFit.Index(pCats),
                                       RooFit.Import("mc_"+pPtLabel+"_pass",self._lPDatas["mc_"+pPtLabel+"_pass"]),
                                       RooFit.Import("mc_"+pPtLabel+"_fail",self._lPDatas["mc_"+pPtLabel+"_fail"]))

            # simultaneous fit with tt model
            simPdf_total_data = r.RooSimultaneous("simPdf_total_data","simPdf_total_data",pCats)
            simPdf_total_mc   = r.RooSimultaneous("simPdf_total_mc"  ,"simPdf_total_mc"  ,pCats)

            for Cat_it in fCats:
                simPdf_total_data.addPdf(self._lW.pdf("model_total_data_"+pPtLabel+"_"+Cat_it),"data_"+pPtLabel+"_"+Cat_it)
                simPdf_total_mc  .addPdf(self._lW.pdf("model_total_mc_"  +pPtLabel+"_"+Cat_it),"mc_"  +pPtLabel+"_"+Cat_it)

            # do simult fit with tt model
            print "simultaneous pass and fail fit with tt model"
            simFit_total_mc   = simPdf_total_mc  .fitTo(combData_mc,r.RooFit.Save(r.kTRUE),r.RooFit.Verbose(r.kFALSE),r.RooFit.Minimizer("Minuit2"),r.RooFit.SumW2Error(r.kTRUE))
            simFit_total_mc.Print()

            simFit_total_data = simPdf_total_data.fitTo(combData_data,r.RooFit.Save(r.kTRUE),r.RooFit.Verbose(r.kFALSE),r.RooFit.Minimizer("Minuit2"),r.RooFit.SumW2Error(r.kTRUE))
            simFit_total_data.Print()

            # get Wtag with tt model
            self.getWtagSFs(pPtLabel,fPDFs["tt_realW_mc"],"tt_realW_mc");
            
            # draw simult fit with tt model
            params = {}                                                         
            drawDataMc(lVar,self._lPDatas["data_"+pPtLabel+"_pass"],
                       [self._lW.pdf("model_total_data_"+pPtLabel+"_pass")],
                       self._lPDatas["mc_"+pPtLabel+"_pass"],
                       [self._lW.pdf("model_total_mc_"+pPtLabel+"_pass")],
                       simFit_total_mc,simFit_total_data,params,"data_pass_simult_tt",
                       fBinWidth)

            drawDataMc(lVar,self._lPDatas["data_"+pPtLabel+"_fail"],
                       [self._lW.pdf("model_total_data_"+pPtLabel+"_fail")],
                       self._lPDatas["mc_"+pPtLabel+"_fail"],
                       [self._lW.pdf("model_total_mc_"+pPtLabel+"_fail")],
                       simFit_total_mc,simFit_total_data,params,"data_fail_simult_tt",
                       fBinWidth)

            # write workspace
            self._lW.Print()
            self._lW.writeToFile(fOutput)

if __name__ == "__main__":
    options = parser()
    print options
    global fTag,fXMin,fXMax,fNBins
    fTag   = options.tag
    fXMin  = options.xmin
    fXMax  = options.xmax
    fNBins = int( (fXMax - fXMin) / fBinWidth)
    # get roodataset and make individual fits
    lW = WPeak(options);
    # combined fit
    if options.combined:
        lW.fit();
