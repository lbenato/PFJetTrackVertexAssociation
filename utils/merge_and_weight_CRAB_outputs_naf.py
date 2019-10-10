#!/usr/bin/env python
import os, re
import multiprocessing
import logging
import commands
import math, time
import sys
from ROOT import TObject, TFile, TH1, TH1F
from Analyzer.PFJetTrackVertexAssociation.samples import sample
from array import array

print "*****************************************************************************"
print "\n"
print "Please input the correct lumi!"
print "\n"
LUMI = 1
print LUMI, " fb -1"
print "*****************************************************************************"



# use the following lists to include/exclude samples to be merged

blacklist = []
whitelist = []

#TIP = "/pnfs/desy.de/cms/tier2/store/user/lbenato/"
#DEST = "/nfs/dust/cms/user/lbenato/v1/"


########## DO NOT TOUCH BELOW THIS POINT ##########

import optparse

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-l", "--lists", action="store", type="string", dest="lists", default="")
parser.add_option("-i", "--input_folder", action="store", type="string", dest="input_folder", default="", help='the input folder containing the CRAB outputs')
parser.add_option("-o", "--output_folder", action="store", type="string", dest="output_folder", default="", help='the output folder containing the hadd of CRAB outputs')
parser.add_option("-g", "--groupofsamples", action="store", type="string", dest="groupofsamples", default="")
(options, args) = parser.parse_args()


from Analyzer.PFJetTrackVertexAssociation.samples import sample,samples

list_of_samples = ["QCD","QCD_Pt_15to30","QCD_Pt_30to50","QCD_Pt_80to120"]
print "Possible subgroups of samples:"
for a in list_of_samples:
    print a
    print "---------------"

from Analyzer.PFJetTrackVertexAssociation.crab_requests_lists import *

selected_requests = {}
if options.groupofsamples not in list_of_samples:
    print "Invalid subgroup of samples, aborting!"
    exit()

for b, k in enumerate(requests.keys()):
    if k in samples[options.groupofsamples]["files"]:
        print k
        selected_requests[k] = requests[k]


DEST = options.output_folder+'/'

if not os.path.exists(os.path.expandvars(options.input_folder)):
    print '--- ERROR: INPUT FOLDER ---'
    print '  \''+options.input_folder+'\' path not found'
    print '  please point to the correct path to the folder containing the CRAB output' 
    exit()

if not os.path.exists(os.path.expandvars(DEST)):
    print '--- ERROR: OUTPUT FOLDER ---'
    print '  \''+DEST+'\' path not found'
    print '  please point to the correct output path' 
    exit()


#########

jobs = []
names = []

def hadd_outputs(fold,name):
    if "_PRIVATE-MC" in name:
        short_name = name[:-11]
    else:
        short_name = name

######################This blocks naf machines
    #print name
    #os.system('hadd -k -f '+DEST+name+'.root '+fold+'/*/*/*/output_7*.root')# + ' ' +name+'/*/*/*/*_1.root')
    os.system('hadd -k -f '+DEST+name+'.root ' + fold + "/*/*/*/*.root")
pass

def weight(name):
    weight = 1.
    filename = TFile(DEST+name+'.root', "UPDATE")
    if ('Run2016') in name: weight = 1.
    ###
    # If you want to weight only one sample, specify
    #elif ('TT_TuneCUETP8M2T4_13TeV') in name:
    ###
    else:
        nevents = filename.Get("counter/c_nEvents").GetBinContent(1)
#        nevents = filename.Get("counter/c_nEvents").GetEntries()#try?!
        if ('VBFH_HToSSTobbbb') in name:
            xs = 1.
        else:
            xs = sample[name]['xsec']
        weight = LUMI * xs / nevents
        #tree = filename.Get("ntuple/tree")
        tree = filename.Get("trigger/tree")
        print name
        print weight
        tree.SetWeight(weight)
        tree.AutoSave()


subdirs = [x for x in os.listdir(options.input_folder) if os.path.isdir(os.path.join(options.input_folder, x))]
#subdirs have the names of the samples without v1, etc

#for naming purposes, we have to include v1, etc. Additional loop###
#print os.listdir(args.folder)
os.chdir(options.input_folder)

crab_subdirs = []
for l in subdirs:
    crab_subdirs += [x[5:] for x in os.listdir(l) if os.path.isdir(os.path.join(l, x))]
#here they have the proper names, including v1

os.chdir(options.input_folder)

for l in subdirs:
    fold = ""
    name = ""
    for a in crab_subdirs:
        #if l in a:
        if l in a and a in selected_requests.keys():
            fold = l
            name = a
            #print fold
            print "Being added...."
            print name
            #print "Not being added...."
            hadd_outputs(fold,name)
            print "##################################"

######################
#
#    Multiprocessing stucked naf machines, avoid - also, not tested with optparse
#
#    p = multiprocessing.Process(target=hadd_outputs, args=(fold,name))
#    jobs.append(p)
#    p.start()
######################
            
    #hadd_outputs(fold,name)

print "Ntuples ready in ", DEST
os.system('cd '+DEST+".. ")


#onlyfiles = [f for f in os.listdir(DEST) if (os.path.isfile(os.path.join(DEST, f)))]
#os.chdir(DEST)

#for b in onlyfiles:
#    #print b
#    if b[:-5] in selected_requests.keys():
#        print "I am going to weight:"
#        #print b
#        weight(b[:-5])
#        ##q = multiprocessing.Process(target=weight, args=(b[:-5],))
#        ##jobs.append(q)
#        ##q.start()                                                                             #                                                                                           
#print "Ntuples weighted in ", DEST
