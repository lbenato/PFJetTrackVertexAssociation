#! /usr/bin/env python                                                          
                                                                                
import os, multiprocessing                                                      
import copy                                                                     
import math                                                                     
import numpy as np
from array import array                                                         
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory                 
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGraphErrors, TMultiGraph, TProfile, TF1
from ROOT import TStyle, TCanvas, TPad                                          
from ROOT import TLegend, TLatex, TText, TLine


#### IMPORT SAMPLES AND VARIABLES DICTIONARIES ####
from Analyzer.PFJetTrackVertexAssociation.samples import sample, samples
from Analyzer.PFJetTrackVertexAssociation.variables import *
from Analyzer.PFJetTrackVertexAssociation.selections import *
from Analyzer.PFJetTrackVertexAssociation.drawUtils import *

#### PARSER ####
import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-b", "--bash", action="store_true", default=False, dest="bash")
parser.add_option("-v", "--variable", action="store", type="string", dest="variable", default="")
parser.add_option("-c", "--cut", action="store", type="string", dest="cut", default="")
(options, args) = parser.parse_args()
if options.bash: gROOT.SetBatch(True)
gStyle.SetOptStat(0)

#### NTUPLE, PLOT DIRECTORIES ####
#NTUPLEDIR   = "/pnfs/desy.de/cms/tier2/store/user/lbenato/PFCHS_v0_debug_proxy_naf/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/crab_QCD_Pt_80to120_TuneCP5_13TeV_pythia8-v1/190604_134218/0000/"
#NTUPLEDIR   = "/afs/desy.de/user/l/lbenato/PF_track_vertex_association/CMSSW_9_4_1/src/Analyzer/PFJetTrackVertexAssociation/local_samples/"
#NTUPLEDIR   = "/pnfs/desy.de/cms/tier2/store/user/lbenato/PFCHS_v1_debug_proxy_naf/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/crab_QCD_Pt_15to30_TuneCP5_13TeV_pythia8-v1/190708_113615/0000/"
#NTUPLEDIR   = "/pnfs/desy.de/cms/tier2/store/user/lbenato/PFCHS_v1_debug_proxy_naf/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/crab_QCD_Pt_30to50_TuneCP5_13TeV_pythia8-v1/190708_122342/0000/"
#NTUPLEDIR   = "/pnfs/desy.de/cms/tier2/store/user/lbenato/PFCHS_v1_debug_proxy_naf/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/crab_QCD_Pt_80to120_TuneCP5_13TeV_pythia8-v1/190708_152325/0000/"

#NTUPLEDIR = "/nfs/dust/cms/user/lbenato/PFCHS_v1/"
NTUPLEDIR   = "/pnfs/desy.de/cms/tier2/store/user/lbenato/PFCHS_v0_debug_proxy_naf/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/crab_QCD_Pt_80to120_TuneCP5_13TeV_pythia8-v1/190604_134218/0000/"


PLOTDIR     = "/afs/desy.de/user/l/lbenato/PF_track_vertex_association/CMSSW_9_4_1/src/Analyzer/PFJetTrackVertexAssociation/macro/plots_v0/"
LUMI        = 1
SIGNAL      = 1
POISSON     = False

#### SAMPLES ####
#sign = ["QCD_Pt_15to30"]
#sign = ["QCD_Pt_30to50"]
#sign = ["QCD_Pt_80to120"]
sign = ["QCD"]


print "************************************"
print "Bringing these ntuples: ", NTUPLEDIR
print "And these samples: ", sign
print "************************************"


colors = [4, 410, 2, 881, 801, 856, 798, 634, 1, 798, 602, 921, 3, 5, 6, ]
lines = [1, 2, 3, 4, 5, 6, 7, 8]
markers = [24,25,26,32,36,24,25]
nomefile = {}
from collections import defaultdict
hist = defaultdict(dict)
fit0 = defaultdict(dict)
fit1 = defaultdict(dict)
hist_profX = defaultdict(dict)
uncert = defaultdict(dict)
mean_resolution = defaultdict(dict)
mean_uncertainty = defaultdict(dict)
RMS_resolution = defaultdict(dict)
sigma_resolution = defaultdict(dict)
sigma_uncertainty = defaultdict(dict)
chain = {}

gStyle.SetOptStat(0)

def compare_variables(sign, var_list, cut_gen, name, cut="", treename="ntuple/tree"):
    cv = TCanvas("cv","cv",1000,800)
    cv.SetGrid()
    leg = TLegend(0.62, 0.7, 0.95, 0.98)
    leg.SetHeader("Cut on PF candidates")
    massimo = 0.
    minimo  = 999.
    for b, s in enumerate(sign):
        for a, var in enumerate(var_list):
            if cut_gen:
                first, second = var.split('.')
                cutstring = first+".hasGenJ" 
                if first=="Jets":
                    label = "slimmedJets"
                elif first=="JetsNew1":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
                    label = "abs(dz)<x"
                elif first=="JetsNew2":
                    label = "fromPV()>0 || charge() ==0"
                    label = "fromPV()>0"
                elif first=="JetsNew3":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
                    label = "abs(dz)<x or abs(dz)<pt*dzError"
                elif first=="JetsNew4":
                    label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
                    label = "dzAssociatedPV()"
                elif first=="JetsNew5":
                    label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
                    label = "PUPPI-like"
                else:
                    label = ""
            else:
                first, second = var.split('n')
                if second=="Jets":
                    label = "slimmedJets"
                elif second=="JetsNew1":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
                    label = "abs(dz)<x"
                elif second=="JetsNew2":
                    label = "fromPV()>0 || charge() ==0"
                    label = "fromPV()>0"
                elif second=="JetsNew3":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
                    label = "abs(dz)<x or abs(dz)<pt*dzError"
                elif second=="JetsNew4":
                    label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
                    label = "dzAssociatedPV()"
                elif second=="JetsNew5":
                    label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
                    label = "PUPPI-like"
                else:
                    label = ""

                first = ""
                first = second
                cutstring = first+".hasGenJ"

            if cutstring=="":
                cutstring = cut
            else:
                if cut!="":
                    cutstring += " && " + first + "."+cut

            print "Variable: ", var
            print "Cut: ", cutstring
            print "*****"
            chain[s] = TChain(treename)                
            for j, ss in enumerate(samples[s]['files']):
                chain[s].Add(NTUPLEDIR + ss + ".root")
            if variable[var]['nbins']>0:
                hist[s][var] = TH1F(s+var, ";"+variable[var]['title'], variable[var]['nbins'], variable[var]['min'], variable[var]['max'])
            hist[s][var].Sumw2()
            chain[s].Project(s+var, var, cutstring)


            hist[s][var].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
            hist[s][var].Scale(samples[s]['weight'] if hist[s][var].Integral() >= 0 else 0)
            hist[s][var].SetLineColor(colors[a])
            hist[s][var].SetFillColor(colors[a])
            hist[s][var].SetLineWidth(2)
            hist[s][var].SetMarkerStyle(markers[a])
            hist[s][var].SetMarkerColor(colors[a])
            hist[s][var].SetTitle("")
            #hist[s][var].SetFillStyle(0)
            hist[s][var].SetLineStyle(lines[b])#original: lines[b] to have different styles for different backgrounds 
            addOverflow(hist[s][var], False)
            massimo = max(massimo, hist[s][var].GetMaximum())
            minimo = min(minimo, hist[s][var].GetMinimum())
            cv.cd()
            hist[s][var].Draw("SAMES, PL")# if (a>0 and b>0) else "HISTO")
            leg.AddEntry(hist[s][var],label,'PL')
            uncert[s][var] = hist[s][var].Clone(s+"_err")
            uncert[s][var].SetMarkerStyle(0)
            uncert[s][var].SetFillColor(colors[a])
            uncert[s][var].SetFillStyle(3001)
            uncert[s][var].Draw("SAME,E2")

        if variable[var]['log']:

            cv.SetLogy()


    for b, s in enumerate(sign):
        for a, var in enumerate(var_list):
            hist[s][var].SetMaximum(massimo*1.5)
            if variable[var]['log']==False:
                hist[s][var].SetMinimum(minimo*0.5)
    #leg.SetHeader(samples[sample]['leg_name'])
    #massimo = 0
    leg.Draw()
    drawCMS(-1,"Preliminary")
    cv.Update()
    cv.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+name+".pdf")
    cv.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+name+".png")
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")



def compare_2D_vs_eta(sign, var_list, cut_gen, name, cut="", treename="ntuple/tree"):
    cv = TCanvas("cv","cv",1000,800)
    cv.SetGrid()
    leg = TLegend(0.62, 0.3, 0.95, 0.58)
    leg.SetHeader("Cut on PF candidates")
    massimo = 0.
    minimo  = 999.
    
    for b, s in enumerate(sign):
        for a, var in enumerate(var_list):
            if cut_gen:
                first, second = var.split('.')
                cutstring = first+".hasGenJ" 
                if first=="Jets":
                    label = "slimmedJets"
                elif first=="JetsNew1":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
                    label = "abs(dz)<x"
                elif first=="JetsNew2":
                    label = "fromPV()>0 || charge() ==0"
                    label = "fromPV()>0"
                elif first=="JetsNew3":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
                    label = "abs(dz)<x or abs(dz)<pt*dzError"
                elif first=="JetsNew4":
                    label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
                    label = "dzAssociatedPV()"
                elif first=="JetsNew5":
                    label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
                    label = "PUPPI-like"
                else:
                    label = ""
            else:
                first, second = var.split('n')
                if second=="Jets":
                    label = "slimmedJets"
                elif second=="JetsNew1":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
                    label = "abs(dz)<x"
                elif second=="JetsNew2":
                    label = "fromPV()>0 || charge() ==0"
                    label = "fromPV()>0"
                elif second=="JetsNew3":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
                    label = "abs(dz)<x or abs(dz)<pt*dzError"
                elif second=="JetsNew4":
                    label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
                    label = "dzAssociatedPV()"
                elif second=="JetsNew5":
                    label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
                    label = "PUPPI-like"
                else:
                    label = ""
                
                    cutstring = ""

            if cutstring=="":
                cutstring = cut
            else:
                if cut!="":
                    cutstring += " && " + first + "."+cut

            var_x = first + ".eta"
            print "Variable y: ", var
            print "Variable x: ", var_x
            print "Cut: ", cutstring
            print "*****"
            chain[s] = TChain(treename)                
            for j, ss in enumerate(samples[s]['files']):
                chain[s].Add(NTUPLEDIR + ss + ".root")


            if (variable[var_x]['nbins']>0 and variable[var]['nbins']>0):
                print var_x, var
                hist[s][var] = TH2F(s+var,variable[var]['title']+" vs "+variable[var_x]['title'],variable[var_x]['nbins'],0, variable[var_x]['max'],variable[var]['nbins'],variable[var]['min'], variable[var]['max'])
                hist[s][var].Sumw2()
                chain[s].Project(s+var, var+":fabs("+var_x+")", cutstring)
                hist[s][var].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
                hist[s][var].GetXaxis().SetTitle("jets |#eta|")
                hist[s][var].SetTitle("")
                hist[s][var].GetYaxis().SetTitle(variable[var]['title'])
                hist[s][var].GetZaxis().SetRangeUser(0., hist[s][var].GetMaximum())
                hist[s][var].SetMarkerColor(colors[a])
            #hist[s].Draw("")
            #if variable[var]['log']:
            #    cv.SetLogy()

            cv.cd()
            hist_profX[s][var] = TProfile(hist[s][var].ProfileX(s+var+"prof"))
            massimo = max(massimo,hist_profX[s][var].GetMaximum())
            hist_profX[s][var].SetMaximum(1.5)
            hist_profX[s][var].SetLineColor(colors[a])
            hist_profX[s][var].SetFillColor(1)
            hist_profX[s][var].SetLineWidth(2)
            hist_profX[s][var].SetMarkerStyle(markers[a])
            hist_profX[s][var].SetMarkerColor(colors[a])
            hist_profX[s][var].GetXaxis().SetTitle("jets |#eta|")
            hist_profX[s][var].GetYaxis().SetTitle("<"+variable[var_list[0]]['title']+">")
            hist_profX[s][var].SetTitle("")
            hist_profX[s][var].Draw("PL,sames")
            hist_profX[s][var].SetMaximum(massimo*(1+0.1))
            leg.AddEntry(hist_profX[s][var],label,'PL')


    #for b, s in enumerate(sign):
        #for a, var in enumerate(var_list):
        #hist[s][var].SetMaximum(massimo*1.5)
        #if variable[var]['log']==False:
        #hist[s][var].SetMinimum(minimo*0.5)
        leg.Draw()
        drawCMS(-1,"Preliminary")
        cv.Update()
        cv.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_vs_"+var_x.replace('.', '_')+name+".pdf")
        cv.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_vs_"+var_x.replace('.', '_')+name+".png")
        if not gROOT.IsBatch(): raw_input("Press Enter to continue...")


def fractions_2D_vs_eta(sign, var_list, name, cut="", cut_gen=True, treename="ntuple/tree"):
    cv = TCanvas("cv","cv",1000,800)
    cv.SetGrid()
    leg = TLegend(0.62, 0.3, 0.95, 0.58)
    massimo = 0.
    minimo  = 999.
    
    for b, s in enumerate(sign):
        for a, var in enumerate(var_list):
            first, second = var.split('.')
            cutstring = first+".hasGenJ" 
            if first=="Jets":
                label = "slimmedJets"
            elif first=="JetsNew1":
                label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
                label = "abs(dz)<x"
            elif first=="JetsNew2":
                label = "fromPV()>0 || charge() ==0"
                label = "fromPV()>0"
            elif first=="JetsNew3":
                label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
                label = "abs(dz)<x or abs(dz)<pt*dzError"
            elif first=="JetsNew4":
                label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
                label = "dzAssociatedPV()"
            elif first=="JetsNew5":
                label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
                label = "PUPPI-like"
            else:
                label = ""

            leg_label = second
            if cutstring=="":
                cutstring = cut
            else:
                if cut!="":
                    cutstring += " && " + first + "."+cut

            var_x = first + ".eta"
            print "Variable y: ", var
            print "Variable x: ", var_x
            print "Cut: ", cutstring
            print "*****"
            chain[s] = TChain(treename)                
            for j, ss in enumerate(samples[s]['files']):
                chain[s].Add(NTUPLEDIR + ss + ".root")


            if (variable[var_x]['nbins']>0 and variable[var]['nbins']>0):
                print var_x, var
                hist[s][var] = TH2F(s+var,variable[var]['title']+" vs "+variable[var_x]['title'],variable[var_x]['nbins'],0, variable[var_x]['max'],variable[var]['nbins'],variable[var]['min'], variable[var]['max'])
                hist[s][var].Sumw2()
                chain[s].Project(s+var, var+":fabs("+var_x+")", cutstring)
                hist[s][var].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
                hist[s][var].GetXaxis().SetTitle("jets |#eta|")
                hist[s][var].SetTitle("")
                hist[s][var].GetYaxis().SetTitle(variable[var]['title'])
                hist[s][var].GetZaxis().SetRangeUser(0., hist[s][var].GetMaximum())
                hist[s][var].SetMarkerColor(colors[a])
            #hist[s].Draw("")
            #if variable[var]['log']:
            #    cv.SetLogy()

            cv.cd()
            hist_profX[s][var] = TProfile(hist[s][var].ProfileX(s+var+"prof"))
            hist_profX[s][var].SetLineColor(colors[a])
            hist_profX[s][var].SetFillColor(1)
            hist_profX[s][var].SetLineWidth(2)
            hist_profX[s][var].SetMarkerStyle(markers[a])
            hist_profX[s][var].SetMarkerColor(colors[a])
            hist_profX[s][var].GetXaxis().SetTitle("jets |#eta|")
            hist_profX[s][var].GetYaxis().SetTitle("<energy fractions>")
            hist_profX[s][var].SetTitle("")
            hist_profX[s][var].Draw("PL,sames")
            leg.AddEntry(hist_profX[s][var],leg_label,'PL')


    #for b, s in enumerate(sign):
        #for a, var in enumerate(var_list):
        #hist[s][var].SetMaximum(massimo*1.5)
        #if variable[var]['log']==False:
        #hist[s][var].SetMinimum(minimo*0.5)
        leg.SetHeader("PF candidates: "+label)
        leg.Draw()
        drawCMS(-1,"Preliminary")
        cv.Update()
        cv.Print(PLOTDIR+sign[0]+"_"+first+"_vs_"+var_x.replace('.', '_')+name+".pdf")
        cv.Print(PLOTDIR+sign[0]+"_"+first+"_vs_"+var_x.replace('.', '_')+name+".png")
        if not gROOT.IsBatch(): raw_input("Press Enter to continue...")



def multiplicities_2D_vs_eta(sign, var_list, name, cut="", cut_gen=True, treename="ntuple/tree"):
    cv = TCanvas("cv","cv",1000,800)
    cv.SetGrid()
    leg = TLegend(0.62, 0.3+0.4, 0.95, 0.58+0.4)
    massimo = 0.
    minimo  = 999.
    
    for b, s in enumerate(sign):
        for a, var in enumerate(var_list):
            first, second = var.split('.')
            cutstring = first+".hasGenJ" 
            if first=="Jets":
                label = "slimmedJets"
            elif first=="JetsNew1":
                label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
                label = "abs(dz)<x"
            elif first=="JetsNew2":
                label = "fromPV()>0 || charge() ==0"
                label = "fromPV()>0"
            elif first=="JetsNew3":
                label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
                label = "abs(dz)<x or abs(dz)<pt*dzError"
            elif first=="JetsNew4":
                label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
                label = "dzAssociatedPV()"
            elif first=="JetsNew5":
                label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
                label = "PUPPI-like"
            else:
                label = ""

            leg_label = second
            if cutstring=="":
                cutstring = cut
            else:
                if cut!="":
                    cutstring += " && " + first + "."+cut

            var_x = first + ".eta"
            print "Variable y: ", var
            print "Variable x: ", var_x
            print "Cut: ", cutstring
            print "*****"
            chain[s] = TChain(treename)                
            for j, ss in enumerate(samples[s]['files']):
                chain[s].Add(NTUPLEDIR + ss + ".root")


            if (variable[var_x]['nbins']>0 and variable[var]['nbins']>0):
                print var_x, var
                hist[s][var] = TH2F(s+var,variable[var]['title']+" vs "+variable[var_x]['title'],variable[var_x]['nbins'],0, variable[var_x]['max'],variable[var]['nbins'],variable[var]['min'], variable[var]['max'])
                hist[s][var].Sumw2()
                chain[s].Project(s+var, var+":fabs("+var_x+")", cutstring)
                hist[s][var].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
                hist[s][var].GetXaxis().SetTitle("jets |#eta|")
                hist[s][var].SetTitle("")
                hist[s][var].GetYaxis().SetTitle(variable[var]['title'])
                hist[s][var].GetZaxis().SetRangeUser(0., hist[s][var].GetMaximum())
                hist[s][var].SetMarkerColor(colors[a])

            cv.cd()
            hist_profX[s][var] = TProfile(hist[s][var].ProfileX(s+var+"prof"))
            hist_profX[s][var].SetLineColor(colors[a])
            hist_profX[s][var].SetMinimum(0)
            hist_profX[s][var].SetFillColor(1)
            hist_profX[s][var].SetLineWidth(2)
            hist_profX[s][var].SetMarkerStyle(markers[a])
            hist_profX[s][var].SetMarkerColor(colors[a])
            hist_profX[s][var].GetXaxis().SetTitle("jets |#eta|")
            hist_profX[s][var].GetYaxis().SetTitle("particle multiplicities")
            hist_profX[s][var].SetTitle("")
            hist_profX[s][var].Draw("PL,sames")
            leg.AddEntry(hist_profX[s][var],leg_label,'PL')


        leg.SetHeader("PF candidates: "+label)
        leg.Draw()
        drawCMS(-1,"Preliminary")
        cv.Update()
        cv.Print(PLOTDIR+sign[0]+"_"+first+"_vs_"+var_x.replace('.', '_')+name+".pdf")
        cv.Print(PLOTDIR+sign[0]+"_"+first+"_vs_"+var_x.replace('.', '_')+name+".png")
        if not gROOT.IsBatch(): raw_input("Press Enter to continue...")



def resolution(sign, var_list, cut_gen, name, cut="", treename="ntuple/tree"):
    cv = TCanvas("cv","cv",1000,800)
    cv.SetGrid()
    leg = TLegend(0.62, 0.7, 0.95, 0.98)
    leg.SetHeader("Cut on PF candidates")
    massimo = 0.
    minimo  = 999.
    for b, s in enumerate(sign):
        for a, var in enumerate(var_list):
            if cut_gen:
                first, second = var.split('.')
                cutstring = first+".hasGenJ" 
                if first=="Jets":
                    label = "slimmedJets"
                elif first=="JetsNew1":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
                    label = "abs(dz)<x"
                elif first=="JetsNew2":
                    label = "fromPV()>0 || charge() ==0"
                    label = "fromPV()>0"
                elif first=="JetsNew3":
                    label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
                    label = "abs(dz)<x or abs(dz)<pt*dzError"
                elif first=="JetsNew4":
                    label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
                    label = "dzAssociatedPV()"
                elif first=="JetsNew5":
                    label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
                    label = "PUPPI-like"
                else:
                    label = ""
            else:
                first, second = var.split('n')
                label = "aaa"
                first = ""
                first = second
                cutstring = first+".hasGenJ"

            if cutstring=="":
                cutstring = cut
            else:
                if cut!="":
                    cutstring += " && " + first + "."+cut

            print "Variable: ", var
            print "Cut: ", cutstring
            print "*****"
            chain[s] = TChain(treename)                
            for j, ss in enumerate(samples[s]['files']):
                chain[s].Add(NTUPLEDIR + ss + ".root")
            if variable[var]['nbins']>0:
                hist[s][var] = TH1F(s+var, ";"+variable[var]['title'], variable[var]['nbins'], variable[var]['min'], variable[var]['max'])
            hist[s][var].Sumw2()
            chain[s].Project(s+var, var, cutstring)


            hist[s][var].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
            hist[s][var].Scale(samples[s]['weight'] if hist[s][var].Integral() >= 0 else 0)
            hist[s][var].SetLineColor(colors[a])
            hist[s][var].SetFillColor(colors[a])
            hist[s][var].SetLineWidth(2)
            hist[s][var].SetMarkerStyle(markers[a])
            hist[s][var].SetMarkerColor(colors[a])
            hist[s][var].SetTitle("")
            #hist[s][var].SetFillStyle(0)
            hist[s][var].SetLineStyle(lines[b])#original: lines[b] to have different styles for different backgrounds 
            addOverflow(hist[s][var], False)
            massimo = max(massimo, hist[s][var].GetMaximum())
            minimo = min(minimo, hist[s][var].GetMinimum())
            cv.cd()
            hist[s][var].Draw("SAMES, PL")# if (a>0 and b>0) else "HISTO")
            leg.AddEntry(hist[s][var],label,'PL')
            uncert[s][var] = hist[s][var].Clone(s+"_err")
            uncert[s][var].SetMarkerStyle(0)
            uncert[s][var].SetFillColor(colors[a])
            uncert[s][var].SetFillStyle(3001)
            uncert[s][var].Draw("SAME,E2")
            fit0[s][var] = TF1("f0"+s+var,"gaus",variable[var]['min'],variable[var]['max'])
            print hist[s][var].GetMean()
            print hist[s][var].GetRMS()
            fit0[s][var].SetParameter(1,hist[s][var].GetMean())
            fit0[s][var].SetParameter(2,hist[s][var].GetRMS())
            hist[s][var].Fit("f0"+s+var,"PWMS","",hist[s][var].GetMean()-3*hist[s][var].GetRMS(),hist[s][var].GetMean()+3*hist[s][var].GetRMS())
            RMS_resolution[s][var] = hist[s][var].GetRMS()
            fit1[s][var] = TF1("f1"+s+var,"gaus",variable[var]['min'],variable[var]['max'])
            fit1[s][var].SetParameter(1,fit0[s][var].GetParameter(1))
            fit1[s][var].SetParameter(2,fit0[s][var].GetParameter(2))
            hist[s][var].Fit("f1"+s+var,"PWMS","",fit0[s][var].GetParameter(1)-1.5*fit0[s][var].GetParameter(2),fit0[s][var].GetParameter(1)+1.5*fit0[s][var].GetParameter(2))
            fit0[s][var].SetLineColor(colors[a])
            fit0[s][var].SetLineWidth(2)
            fit0[s][var].SetLineStyle(2)
            fit0[s][var].Draw("L,sames")
            fit1[s][var].SetLineColor(colors[a])
            fit1[s][var].SetLineWidth(2)
            fit1[s][var].SetLineStyle(1)
            fit1[s][var].Draw("L,sames")
            mean_resolution[s][var] = fit1[s][var].GetParameter(1) 
            mean_uncertainty[s][var] = fit1[s][var].GetParError(1)
            sigma_resolution[s][var] = fit1[s][var].GetParameter(2)
            sigma_uncertainty[s][var] = fit1[s][var].GetParError(2)
            print s, var, "sigma unc: ", sigma_uncertainty[s][var]

        if variable[var]['log']:

            cv.SetLogy()


    for b, s in enumerate(sign):
        for a, var in enumerate(var_list):
            hist[s][var].SetMaximum(massimo*1.5)
            if variable[var]['log']==False:
                hist[s][var].SetMinimum(minimo*0.5)
    #leg.SetHeader(samples[sample]['leg_name'])
    #massimo = 0
    leg.Draw()
    drawCMS(-1,"Preliminary")
    cv.Update()
    cv.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+name+".pdf")
    cv.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+name+".png")
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")


def resolution_vs_eta(sign, var, cut_gen, name, min_eta, max_eta, cut="", treename="ntuple/tree", drawDebugFit=False, col_index=-1):
    massimo = 0.
    minimo  = 999.
    if col_index==-1:
        col = 0
    else:
        col = col_index
    if cut_gen:
        first, second = var.split('.')
        cutstring = first+".hasGenJ" 
        if first=="Jets":
            label = "slimmedJets"
        elif first=="JetsNew1":
            label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
            label = "abs(dz)<x"
        elif first=="JetsNew2":
            label = "fromPV()>0 || charge() ==0"
            label = "fromPV()>0"
        elif first=="JetsNew3":
            label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
            label = "abs(dz)<x or abs(dz)<pt*dzError"
        elif first=="JetsNew4":
            label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
            label = "dzAssociatedPV()"
        elif first=="JetsNew5":
            label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
            label = "PUPPI-like"
        else:
            label = ""
    else:
        first, second = var.split('n')
        label = "aaa"
        first = ""
        first = second
        cutstring = first+".hasGenJ"

    if cutstring=="":
        cutstring = cut + ""
    else:
        if cut!="":
            cutstring += " && " + first + "."+cut + " && abs(" + first + ".eta)<"+str(max_eta)+" && abs(" + first + ".eta)>=" + str(min_eta)

    print "Variable: ", var
    print "Cut: ", cutstring
    print "*****"
    chain[sign] = TChain(treename)                
    for j, ss in enumerate(samples[sign]['files']):
        chain[sign].Add(NTUPLEDIR + ss + ".root")
    if variable[var]['nbins']>0:
        hist[sign][var] = TH1F(sign+var, ";"+variable[var]['title'], variable[var]['nbins'], variable[var]['min'], variable[var]['max'])
    hist[sign][var].Sumw2()
    chain[sign].Project(sign+var, var, cutstring)


    hist[sign][var].SetOption("%s" % chain[sign].GetTree().GetEntriesFast())
    hist[sign][var].Scale(samples[sign]['weight'] if hist[sign][var].Integral() >= 0 else 0)
    hist[sign][var].SetLineColor(colors[col])
    hist[sign][var].SetFillColor(colors[col])
    hist[sign][var].SetLineWidth(2)
    hist[sign][var].SetMarkerStyle(markers[1])
    hist[sign][var].SetMarkerColor(colors[col])
    hist[sign][var].SetTitle("") 
    hist[sign][var].SetLineStyle(lines[0])#original: lines[0] to have different styles for different backgrounds 
    addOverflow(hist[sign][var], False)
    massimo = max(massimo, hist[sign][var].GetMaximum())
    minimo = min(minimo, hist[sign][var].GetMinimum())
    if drawDebugFit:
        cv = TCanvas("cv","cv",1000,800)
        cv.SetGrid()
        leg = TLegend(0.62, 0.7, 0.95, 0.98)
        leg.SetHeader("|#eta|:["+str(min_eta)+","+str(max_eta)+"]")
        cv.cd()
        hist[sign][var].Draw("SAMES, PL")# if (a>0 and b>0) else "HISTO")
        uncert[sign][var] = hist[sign][var].Clone(sign+"_err")
        uncert[sign][var].SetMarkerStyle(markers[col])
        uncert[sign][var].SetFillColor(colors[col])
        uncert[sign][var].SetFillStyle(3001)
        uncert[sign][var].Draw("SAME,E2")
        leg.AddEntry(hist[sign][var],label,'PL')
    fit0[sign][var] = TF1("f0"+sign+var,"gaus",variable[var]['min'],variable[var]['max'])
    fit0[sign][var].SetParameter(1,hist[sign][var].GetMean())
    fit0[sign][var].SetParameter(2,hist[sign][var].GetRMS())
    hist[sign][var].Fit("f0"+sign+var,"PWMS","",hist[sign][var].GetMean()-4*hist[sign][var].GetRMS(),hist[sign][var].GetMean()+4*hist[sign][var].GetRMS())
    RMS_resolution = hist[sign][var].GetRMS()
    fit1[sign][var] = TF1("f1"+sign+var,"gaus",variable[var]['min'],variable[var]['max'])
    fit1[sign][var].SetParameter(1,fit0[sign][var].GetParameter(1))
    fit1[sign][var].SetParameter(2,fit0[sign][var].GetParameter(2))
    hist[sign][var].Fit("f1"+sign+var,"PWMS","",fit0[sign][var].GetParameter(1)-1.5*fit0[sign][var].GetParameter(2),fit0[sign][var].GetParameter(1)+1.5*fit0[sign][var].GetParameter(2))

    if drawDebugFit:
        fit0[sign][var].SetLineColor(colors[col])
        fit0[sign][var].SetLineWidth(2)
        fit0[sign][var].SetLineStyle(2)
        fit0[sign][var].Draw("L,sames")
        fit1[sign][var].SetLineColor(colors[col])
        fit1[sign][var].SetLineWidth(2)
        fit1[sign][var].SetLineStyle(1)
        fit1[sign][var].Draw("L,sames")
    mean_resolution = fit1[sign][var].GetParameter(1) 
    mean_uncertainty = fit1[sign][var].GetParError(1)
    sigma_resolution = fit1[sign][var].GetParameter(2)
    sigma_uncertainty = fit1[sign][var].GetParError(2)

    if variable[var]['log'] and drawDebugFit:
        cv.SetLogy()


    hist[sign][var].SetMaximum(massimo*1.5)
    if variable[var]['log']==False:
        hist[sign][var].SetMinimum(minimo*0.5)

    if drawDebugFit:
        leg.Draw()
        drawCMS(-1,"Preliminary")
        cv.Update()
        cv.Print(PLOTDIR+sign+"_"+var.replace('.', '_')+name+"_eta_"+str(min_eta).replace('.', 'p')+"_"+str(max_eta).replace('.', 'p')+".pdf")
        cv.Print(PLOTDIR+sign+"_"+var.replace('.', '_')+name+"_eta_"+str(min_eta).replace('.', 'p')+"_"+str(max_eta).replace('.', 'p')+".png")
        if not gROOT.IsBatch(): raw_input("Press Enter to continue...")
    return mean_resolution,mean_uncertainty,sigma_resolution,sigma_uncertainty,RMS_resolution

################################################
############### Variables ######################
################################################

'''
var_list = ["nJets","nJetsNew1","nJetsNew2","nJetsNew3","nJetsNew4","nJetsNew5"]
compare_variables(sign, var_list, cut_gen=False, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")

var_list = ["Jets.pt","JetsNew1.pt","JetsNew2.pt","JetsNew3.pt","JetsNew4.pt","JetsNew5.pt"]
compare_variables(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")

var_list = ["Jets.ptGenJ","JetsNew1.ptGenJ","JetsNew2.ptGenJ","JetsNew3.ptGenJ","JetsNew4.ptGenJ","JetsNew5.ptGenJ"]
compare_variables(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")

var_list = ["Jets.eta","JetsNew1.eta","JetsNew2.eta","JetsNew3.eta","JetsNew4.eta","JetsNew5.eta"]
compare_variables(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")

var_list = ["Jets.response","JetsNew1.response","JetsNew2.response","JetsNew3.response","JetsNew4.response","JetsNew5.response"]
compare_variables(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")
compare_2D_vs_eta(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15")
#exit()

#multi ratio
var_list = ["Jets.GenRecoNeuMulti","JetsNew1.GenRecoNeuMulti","JetsNew2.GenRecoNeuMulti","JetsNew3.GenRecoNeuMulti","JetsNew4.GenRecoNeuMulti","JetsNew5.GenRecoNeuMulti"]
compare_variables(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")
compare_2D_vs_eta(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15")

var_list = ["Jets.GenRecoChMulti","JetsNew1.GenRecoChMulti","JetsNew2.GenRecoChMulti","JetsNew3.GenRecoChMulti","JetsNew4.GenRecoChMulti","JetsNew5.GenRecoChMulti"]
compare_variables(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")
compare_2D_vs_eta(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15")

#exit()

#fractions

var_list = ["Jets.nHadEFrac","Jets.cHadEFrac","Jets.nEmEFrac","Jets.cEmEFrac","Jets.muEFrac"]
fractions_2D_vs_eta(sign, var_list, name="_ptGenJ_15_fractions", cut="ptGenJ>15")

var_list = ["JetsNew1.nHadEFrac","JetsNew1.cHadEFrac","JetsNew1.nEmEFrac","JetsNew1.cEmEFrac","JetsNew1.muEFrac"]
fractions_2D_vs_eta(sign, var_list, name="_ptGenJ_15_fractions", cut="ptGenJ>15")

var_list = ["JetsNew2.nHadEFrac","JetsNew2.cHadEFrac","JetsNew2.nEmEFrac","JetsNew2.cEmEFrac","JetsNew2.muEFrac"]
fractions_2D_vs_eta(sign, var_list, name="_ptGenJ_15_fractions", cut="ptGenJ>15")

var_list = ["JetsNew3.nHadEFrac","JetsNew3.cHadEFrac","JetsNew3.nEmEFrac","JetsNew3.cEmEFrac","JetsNew3.muEFrac"]
fractions_2D_vs_eta(sign, var_list, name="_ptGenJ_15_fractions", cut="ptGenJ>15")

var_list = ["JetsNew4.nHadEFrac","JetsNew4.cHadEFrac","JetsNew4.nEmEFrac","JetsNew4.cEmEFrac","JetsNew4.muEFrac"]
fractions_2D_vs_eta(sign, var_list, name="_ptGenJ_15_fractions", cut="ptGenJ>15")

var_list = ["JetsNew5.nHadEFrac","JetsNew5.cHadEFrac","JetsNew5.nEmEFrac","JetsNew5.cEmEFrac","JetsNew5.muEFrac"]
fractions_2D_vs_eta(sign, var_list, name="_ptGenJ_15_fractions", cut="ptGenJ>15")

#exit()

#multiplicites

var_list = ["Jets.npr","Jets.cMulti","Jets.nMulti","Jets.cHadMulti","Jets.nHadMulti","Jets.eleMulti","Jets.photonMulti"]
multiplicities_2D_vs_eta(sign, var_list, name="_ptGenJ_15_multiplicities", cut="ptGenJ>15")

var_list = ["JetsNew1.npr","JetsNew1.cMulti","JetsNew1.nMulti","JetsNew1.cHadMulti","JetsNew1.nHadMulti","JetsNew1.eleMulti","JetsNew1.photonMulti"]
multiplicities_2D_vs_eta(sign, var_list, name="_ptGenJ_15_multiplicities", cut="ptGenJ>15")

var_list = ["JetsNew2.npr","JetsNew2.cMulti","JetsNew2.nMulti","JetsNew2.cHadMulti","JetsNew2.nHadMulti","JetsNew2.eleMulti","JetsNew2.photonMulti"]
multiplicities_2D_vs_eta(sign, var_list, name="_ptGenJ_15_multiplicities", cut="ptGenJ>15")

var_list = ["JetsNew3.npr","JetsNew3.cMulti","JetsNew3.nMulti","JetsNew3.cHadMulti","JetsNew3.nHadMulti","JetsNew3.eleMulti","JetsNew3.photonMulti"]
multiplicities_2D_vs_eta(sign, var_list, name="_ptGenJ_15_multiplicities", cut="ptGenJ>15")

var_list = ["JetsNew4.npr","JetsNew4.cMulti","JetsNew4.nMulti","JetsNew4.cHadMulti","JetsNew4.nHadMulti","JetsNew4.eleMulti","JetsNew4.photonMulti"]
multiplicities_2D_vs_eta(sign, var_list, name="_ptGenJ_15_multiplicities", cut="ptGenJ>15")

var_list = ["JetsNew5.npr","JetsNew5.cMulti","JetsNew5.nMulti","JetsNew5.cHadMulti","JetsNew5.nHadMulti","JetsNew5.eleMulti","JetsNew5.photonMulti"]
multiplicities_2D_vs_eta(sign, var_list, name="_ptGenJ_15_multiplicities", cut="ptGenJ>15")

exit()
'''


###
###
# Resolution fit:

####Fitting reponse
#eta bins: 0.,1.,1.5,2.5,2.8,3.
var_list = ["Jets.response","JetsNew1.response","JetsNew2.response","JetsNew3.response","JetsNew4.response"]#,"JetsNew5.response"]


#compare_variables(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")
#compare_2D_vs_eta(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15")
#resolution(sign, var_list, cut_gen=True, name="_ptGenJ_15", cut="ptGenJ>15", treename="ntuple/tree")

eta_bins = defaultdict(dict)

eta_bins = {
    '0_1' :
        {
        'min_eta' : 0.,
        'max_eta' : 1.,
        },
    '1_1p5' :
        {
        'min_eta' : 1.,
        'max_eta' : 1.5,
        },
    '1p5_2p5' :
        {
        'min_eta' : 1.5,
        'max_eta' : 2.5,
        },
    '2p5_2p8' :
        {
        'min_eta' : 2.5,
        'max_eta' : 2.8,
        },
    '2p8_3' :
        {
        'min_eta' : 2.8,
        'max_eta' : 3,
        },
    '3_5p2' :
        {
        'min_eta' : 3,
        'max_eta' : 5.2,
        },
}

'''
eta_bins = {
    '0_1' :
        {
        'min_eta' : 0.,
        'max_eta' : 1.,
        },
}
'''

mean = defaultdict(dict)
sigma = defaultdict(dict)
mean_unc = defaultdict(dict)
sigma_unc = defaultdict(dict)
RMS = defaultdict(dict)
graph_mean = defaultdict(dict)
graph_sigma = defaultdict(dict)
graph_RMS = defaultdict(dict)
for a in sign:
    index = 0
    mg_mean = TMultiGraph()
    mg_sigma = TMultiGraph()
    mg_RMS = TMultiGraph()
    leg = TLegend(0.62, 0.7-0.6, 0.95, 0.98-0.6)

    color_index = 0
    for v in var_list:
        for e in eta_bins.keys():
            mean[v][e], mean_unc[v][e], sigma[v][e], sigma_unc[v][e], RMS[v][e] = resolution_vs_eta(a, v, cut_gen=True, name="_ptGenJ_15", min_eta = eta_bins[e]['min_eta'], max_eta=eta_bins[e]['max_eta'], cut="ptGenJ>15", treename="ntuple/tree", drawDebugFit=True, col_index=color_index)
        color_index +=1
    ##one resolution per sample, not to be mixed
    #print "mean: ", mean
    #print "mean unc: ", mean_unc
    #print sigma
    #print sigma_unc

    for v in var_list:
        print mean[v]
        print sigma[v]
        bins =       [0.,  1.,   1.5, 2.5,  2.8, 3.]
        bins =       [0.5, 1.25, 2.,  2.65, 2.9, 4.]
        bins_error = [0.5, 0.25, 0.5, 0.15, 0.1, 1.]
        sorted_mean = []
        sorted_sigma = []
        sorted_mean_unc = []
        sorted_sigma_unc = []
        sorted_RMS = []
        for a in sorted( mean[v].keys() ):
            sorted_mean.append(mean[v][a])
            sorted_sigma.append(sigma[v][a])
            sorted_mean_unc.append(mean_unc[v][a])
            sorted_sigma_unc.append(sigma_unc[v][a])
            sorted_RMS.append(RMS[v][a])
        print sorted_mean
        print sorted_sigma
        first, second = v.split('.')
        if first=="Jets":
            label = "slimmedJets"
        elif first=="JetsNew1":
            label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5)"
            label = "abs(dz)<x"
        elif first=="JetsNew2":
            label = "fromPV()>0 || charge() ==0"
            label = "fromPV()>0"
        elif first=="JetsNew3":
            label = "(abs(dz) < cosh(eta)*(0.02+0.01/pt)*5) ||  (hasTrackDetails() && abs(dz)<pt*dzError )"
            label = "abs(dz)<x or abs(dz)<pt*dzError"
        elif first=="JetsNew4":
            label = "vertexRef().key()==0 ||( dzAssociatedPV() < cosh(eta)*(0.003+0.01/pt)*5) || charge==0"
            label = "dzAssociatedPV()"
        elif first=="JetsNew5":
            label = "fromPV()==PVUsedInFit || ( (fromPV()==PVTight || fromPV()==PVLoose && abs(dz()) < 0.3 )"
            label = "PUPPI-like"
        else:
            label = ""

        graph_mean[v] = TGraphErrors(len(bins),np.array(bins),np.array(sorted_mean), np.array(bins_error), np.array(sorted_mean_unc) )
        graph_sigma[v] = TGraphErrors(len(bins),np.array(bins),np.array(sorted_sigma), np.array(bins_error), np.array(sorted_sigma_unc) )
        graph_RMS[v] = TGraphErrors(len(bins),np.array(bins),np.array(sorted_RMS), np.array(bins_error), np.array(sorted_sigma_unc) )

        graph_mean[v].SetMarkerStyle(markers[index])
        graph_mean[v].SetMarkerSize(1.)
        graph_mean[v].SetMarkerColor(colors[index])
        graph_mean[v].SetLineColor(colors[index])
        graph_mean[v].SetLineWidth(2)
        graph_mean[v].GetXaxis().SetTitle("jets |#eta|")
        graph_mean[v].SetTitle("")
        graph_mean[v].GetYaxis().SetTitle("mean resolution (fit)")
        graph_mean[v].SetMinimum(0.9)

        mg_mean.Add(graph_mean[v])
        leg.AddEntry(graph_mean[v],label,'PL')

        graph_sigma[v].SetMarkerStyle(markers[index])
        graph_sigma[v].SetMarkerSize(1.)
        graph_sigma[v].SetMarkerColor(colors[index])
        graph_sigma[v].SetLineColor(colors[index])
        graph_sigma[v].SetLineWidth(2)
        graph_sigma[v].GetXaxis().SetTitle("jets |#eta|")
        graph_sigma[v].SetTitle("")
        graph_sigma[v].GetYaxis().SetTitle("sigma resolution (fit)")
        graph_sigma[v].SetMinimum(0.9)
        mg_sigma.Add(graph_sigma[v])

        graph_RMS[v].SetMarkerStyle(markers[index])
        graph_RMS[v].SetMarkerSize(1.)
        graph_RMS[v].SetMarkerColor(colors[index])
        graph_RMS[v].SetLineColor(colors[index])
        graph_RMS[v].SetLineWidth(2)
        graph_RMS[v].GetXaxis().SetTitle("jets |#eta|")
        graph_RMS[v].SetTitle("")
        graph_RMS[v].GetYaxis().SetTitle("RMS of jet response")
        graph_RMS[v].SetMinimum(0.9)
        mg_RMS.Add(graph_RMS[v])

        index += 1
        print index


    cmg_mean = TCanvas("multi_graph_mean","A Simple Graph Example",1000,800)
    cmg_mean.cd()
    cmg_mean.SetGrid()
    #mg_mean.SetMinimum(0.9)
    mg_mean.Draw("AP")
    leg.Draw()
    mg_mean.GetXaxis().SetTitle("jets |#eta|")
    mg_mean.SetTitle("")
    mg_mean.GetYaxis().SetTitle("mean resolution (fit)")


    cmg_mean.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_fit_mean_vs_eta.pdf")
    cmg_mean.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_fit_mean_vs_eta.png")


    cmg_sigma = TCanvas("multi_graph_sigma","A Simple Graph Example",1000,800)
    cmg_sigma.cd()
    cmg_sigma.SetGrid()
    #mg_sigma.SetMinimum(0.9)
    mg_sigma.Draw("AP")
    leg.Draw()
    mg_sigma.GetXaxis().SetTitle("jets |#eta|")
    mg_sigma.SetTitle("")
    mg_sigma.GetYaxis().SetTitle("sigma resolution (fit)")

    cmg_sigma.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_fit_sigma_vs_eta.pdf")
    cmg_sigma.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_fit_sigma_vs_eta.png")

    cmg_RMS = TCanvas("multi_graph_RMS","A Simple Graph Example",1000,800)
    cmg_RMS.cd()
    cmg_RMS.SetGrid()
    #mg_RMS.SetMinimum(0.9)
    mg_RMS.Draw("AP")
    leg.Draw()
    mg_RMS.GetXaxis().SetTitle("jets |#eta|")
    mg_RMS.SetTitle("")
    mg_RMS.GetYaxis().SetTitle("RMS of jet response")

    cmg_RMS.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_RMS_vs_eta.pdf")
    cmg_RMS.Print(PLOTDIR+sign[0]+"_"+var_list[0].replace('.', '_')+"_RMS_vs_eta.png")

    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")


    exit()


exit()

