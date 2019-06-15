#! /usr/bin/env python3

from hfile import *
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()

parser.add_argument('-pt',    type=float, default=2000, help='pt cut')    
parser.add_argument('-shape', type=str, default='plain_isd_1.0', help='name of the shape to process')
parser.add_argument('-lumi',  type=float, default=150.0, help='integrated lumi in fb^{-1}')
parser.add_argument('-fstep', type=float, default=0.05, help='steps in gg fraction')
parser.add_argument('-log',   action='store_true', help='when present, assume that the shape is binned in log(v)')

args = parser.parse_args()

script_dir = os.path.dirname(__file__)

# read in the cross-sections split in i->f channels
incl_qq2qq=get_array(os.path.join(script_dir,'2to2/res/lhc13-ct14nlo-ymax4.5-qq2qq.res'))
incl_qq2gg=get_array(os.path.join(script_dir,'2to2/res/lhc13-ct14nlo-ymax4.5-qq2gg.res'))
incl_qg2qg=get_array(os.path.join(script_dir,'2to2/res/lhc13-ct14nlo-ymax4.5-qg2qg.res'))
incl_gg2qq=get_array(os.path.join(script_dir,'2to2/res/lhc13-ct14nlo-ymax4.5-gg2qq.res'))
incl_gg2gg=get_array(os.path.join(script_dir,'2to2/res/lhc13-ct14nlo-ymax4.5-gg2gg.res'))

# convert this infto final-states admixtures
incl_pts=incl_qq2qq[:,0]
incl_sum=incl_qq2qq[:,1]+incl_qq2qq[:,2]+incl_gg2qq[:,1]+incl_gg2qq[:,2]+incl_qg2qg[:,1]+incl_qg2qg[:,2]+incl_qq2gg[:,1]+incl_qq2gg[:,2]+incl_gg2gg[:,1]+incl_gg2gg[:,2]
incl_qq=(incl_qq2qq[:,1]+incl_qq2qq[:,2]+incl_gg2qq[:,1]+incl_gg2qq[:,2])/incl_sum
incl_qg=(incl_qg2qg[:,1]+incl_qg2qg[:,2])/incl_sum
incl_gg=(incl_qq2gg[:,1]+incl_qq2gg[:,2]+incl_gg2gg[:,1]+incl_gg2gg[:,2])/incl_sum

# find the fractions at the requested pt
frac_qq = np.interp(args.pt, incl_pts, incl_qq)
frac_qg = np.interp(args.pt, incl_pts, incl_qg)
frac_gg = np.interp(args.pt, incl_pts, incl_gg)
#print (frac_qq, frac_qg, frac_gg, frac_qq+frac_qg+frac_gg)

# read the shape distributions
qq = get_array(os.path.join(script_dir,'../taggers-code/res/lhc13-py8.230_M13-qq2qq-ptmin'+str(int(args.pt))+'-full.res'), args.shape)
gg = get_array(os.path.join(script_dir,'../taggers-code/res/lhc13-py8.230_M13-gg2gg-ptmin'+str(int(args.pt))+'-full.res'), args.shape)

# buld corresponding "q" efficiencies
eff= np.empty([qq.shape[0]+1,3])

eff[0,0]=qq[0,0]
eff[0,1]=0.0
eff[0,2]=0.0
for i in range(qq.shape[0]):
    eff[i+1,0]=qq[i,2]
    eff[i+1,1]=eff[i,1]+(qq[i,2]-qq[i,0])*qq[i,3]
    eff[i+1,2]=eff[i,2]+(qq[i,2]-qq[i,0])*gg[i,3]

# deduce the fraction of gluons in the sample after tagging
frac_gg_cut=frac_gg*(1-eff[:,2])*(1-eff[:,2])/(frac_gg*(1-eff[:,2])*(1-eff[:,2])+frac_qg*(1-eff[:,1])*(1-eff[:,2])+frac_qq*(1-eff[:,1])*(1-eff[:,1]))

# we keep only the monotonically increasing part of the gg fraction
for max_index in range(eff[:,0].shape[0]-1):
    if frac_gg_cut[max_index]>1.0001*frac_gg_cut[max_index+1]:
        break

# get the cuts corresponding to some gg fractions
nsteps = int(0.5/args.fstep+0.1)
fracs=[args.fstep*i for i in range(1,nsteps+1)]
cuts=np.interp(fracs, frac_gg_cut[:max_index], eff[:max_index,0])

# compute other stats
#  - fraction of the x-section left after the cuts
q_effs_at_cuts=np.interp(cuts,eff[:,0],eff[:,1])
g_effs_at_cuts=np.interp(cuts,eff[:,0],eff[:,2])
incl_gg_fracs=(frac_gg*(1-g_effs_at_cuts)*(1-g_effs_at_cuts)+frac_qg*(1-q_effs_at_cuts)*(1-g_effs_at_cuts)+frac_qq*(1-q_effs_at_cuts)*(1-q_effs_at_cuts))

# - number of events in 150 fb^{-1} above the requested pt
lumi=args.lumi*1e+6
incl_sum_above=incl_sum.copy()
incl_sum_above[len(incl_sum_above)-1]=0.5*incl_sum[len(incl_sum_above)-1]*(6500-incl_pts[len(incl_sum_above)-1])
for i in range(len(incl_sum_above)-2,-1,-1):
    incl_sum_above[i]=incl_sum_above[i+1]+0.5*(incl_sum[i]+incl_sum[i+1])*(incl_pts[i+1]-incl_pts[i])
xs_pt=np.interp(args.pt, incl_pts, incl_sum_above)
nev_init=xs_pt*lumi
nevs=xs_pt*lumi*incl_gg_fracs

if args.log:
    cuts=np.exp(cuts)

print ('# number of events above pt=',args.pt,': ', nev_init)

print ("# gg frac    shape     kept    expected")
print ("# post cut    cut     gg frac   njets")
for frac, cut, incl_gg_frac, nev in zip(fracs, cuts,incl_gg_fracs, nevs):
    print (f"    {frac:.2f}   {cut:8.5f}   {incl_gg_frac:.4f}   {nev:7.2f}")

