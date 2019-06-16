#! /usr/bin/env python3

from hfile import *
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()

parser.add_argument('-pt',    type=float, default=2000, help='pt cut')    
parser.add_argument('-shape', type=str, default='plain_isd_1.0', help='name of the shape to process')
parser.add_argument('-lumi',  type=float, default=150.0, help='integrated lumi in fb^{-1}')
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
frac_qq_born = np.interp(args.pt, incl_pts, incl_qq)
frac_qg_born = np.interp(args.pt, incl_pts, incl_qg)
frac_gg_born = np.interp(args.pt, incl_pts, incl_gg)
#print (frac_qq, frac_qg, frac_gg, frac_qq+frac_qg+frac_gg)


# - number of events in 150 fb^{-1} above the requested pt
lumi=args.lumi*1e+6
incl_sum_above=incl_sum.copy()
incl_sum_above[len(incl_sum_above)-1]=0.5*incl_sum[len(incl_sum_above)-1]*(6500-incl_pts[len(incl_sum_above)-1])
for i in range(len(incl_sum_above)-2,-1,-1):
    incl_sum_above[i]=incl_sum_above[i+1]+0.5*(incl_sum[i]+incl_sum[i+1])*(incl_pts[i+1]-incl_pts[i])
xs_pt=np.interp(args.pt, incl_pts, incl_sum_above)
nev_init=xs_pt*lumi
nev_qq_born=nev_init*frac_qq_born
nev_qg_born=nev_init*frac_qg_born
nev_gg_born=nev_init*frac_gg_born

# read the shape distributions
q_distrib = get_array(os.path.join(script_dir,'../taggers-code/res/lhc13-py8.230_M13-qq2qq-ptmin'+str(int(args.pt))+'-full.res'), args.shape)
g_distrib = get_array(os.path.join(script_dir,'../taggers-code/res/lhc13-py8.230_M13-gg2gg-ptmin'+str(int(args.pt))+'-full.res'), args.shape)

# buld corresponding "q" efficiencies
eff= np.empty([q_distrib.shape[0]+1,3])

eff[0,0]=q_distrib[0,0]
eff[0,1]=0.0
eff[0,2]=0.0
for i in range(q_distrib.shape[0]):
    eff[i+1,0]=q_distrib[i,2]
    width=(q_distrib[i,2]-q_distrib[i,0])

    # efficienies fotr a cut at the current value of the shape
    eff[i+1,1]=eff[i,1]+width*q_distrib[i,3]
    eff[i+1,2]=eff[i,2]+width*g_distrib[i,3]

# weighted variant
positive=(q_distrib[:,3]+g_distrib[:,3])>0
weight  = g_distrib[positive,3]/(q_distrib[positive,3]+g_distrib[positive,3])
wgt_p_q = (q_distrib[positive,2]-q_distrib[positive,0])*q_distrib[positive,3]
wgt_p_g = (q_distrib[positive,2]-q_distrib[positive,0])*g_distrib[positive,3]

wgt_eff_q = np.sum(wgt_p_q*weight)
wgt_eff_g = np.sum(wgt_p_g*weight)

wgt2_eff_q = np.sum(wgt_p_q*weight*weight)
wgt2_eff_g = np.sum(wgt_p_g*weight*weight)

# fraction of the inclusive x-sect tagged as gluons
ggtag_frac_born_cut=(frac_gg_born*(1-eff[:,2])*(1-eff[:,2])+frac_qg_born*(1-eff[:,1])*(1-eff[:,2])+frac_qq_born*(1-eff[:,1])*(1-eff[:,1]))
ggtag_frac_born_wgt=(frac_gg_born*   wgt_eff_g*   wgt_eff_g+frac_qg_born*   wgt_eff_g*   wgt_eff_q+frac_qq_born*   wgt_eff_q*   wgt_eff_q)

# fraction of gluons in the sample after tagging i.e. purity
frac_gg_cut=frac_gg_born*(1-eff[:,2])*(1-eff[:,2])/ggtag_frac_born_cut
frac_gg_wgt=frac_gg_born*   wgt_eff_g*   wgt_eff_g/ggtag_frac_born_wgt

frac_qg_cut=frac_qg_born*(1-eff[:,1])*(1-eff[:,2])/ggtag_frac_born_cut
frac_qg_wgt=frac_qg_born*   wgt_eff_q*   wgt_eff_g/ggtag_frac_born_wgt

frac_qq_cut=frac_qq_born*(1-eff[:,1])*(1-eff[:,1])/ggtag_frac_born_cut
frac_qq_wgt=frac_qq_born*   wgt_eff_q*   wgt_eff_q/ggtag_frac_born_wgt

# number of events after cut
nevs_cut=nev_init*ggtag_frac_born_cut
nevs_wgt=nev_init*ggtag_frac_born_wgt

nevs_gg_cut=nevs_cut*frac_gg_cut
nevs_qg_cut=nevs_cut*frac_qg_cut
nevs_qq_cut=nevs_cut*frac_qq_cut

nevs_gg_wgt=nevs_wgt*frac_gg_wgt
nevs_qg_wgt=nevs_wgt*frac_qg_wgt
nevs_qq_wgt=nevs_wgt*frac_qq_wgt

# significance
sigs_gg_cut = nevs_gg_cut/np.sqrt(nevs_gg_cut+nevs_qg_cut+nevs_qq_cut)
sig_gg_wgt  = nevs_gg_wgt/np.sqrt(nev_gg_born*wgt2_eff_g*wgt2_eff_g+nev_qg_born*wgt2_eff_q*wgt2_eff_g+nev_qq_born*wgt2_eff_q*wgt2_eff_q)

print ('# number of events above pt=',args.pt,': ', nev_init)

print ("#   shape  gg frac    signific   expected")
print ("#    cut   post cut   post cut    njets")
for cut, frac, sig, nev in zip(eff[:,0], frac_gg_cut, sigs_gg_cut, nevs_cut):
    print (f"  {cut:7.2f}   {frac:.5f}   {sig:8.4f}  {nev:8.2f}")
print ("\n")
print (f" weighted   {frac_gg_wgt:.5f}   {sig_gg_wgt:8.4f}  {nevs_wgt:8.2f}")

