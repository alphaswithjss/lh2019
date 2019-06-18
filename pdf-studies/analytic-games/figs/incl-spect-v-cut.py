#! /usr/bin/env python3

from hfile import *
import argparse
import numpy as np
import os
import subprocess

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()

parser.add_argument('-pt',    type=float, default=2000, help='pt cut')    
parser.add_argument('-shape', type=str, default='nlh', help='name of the shape to process')
parser.add_argument('-alpha', type=float, default= 0.5, help='alpha parameter (for ECF)')
parser.add_argument('-ktcut', type=float, default= 1.0, help='kt cut (for NSD)')
parser.add_argument('-beta',  type=float, default= 0.0, help='cut on the shape')
parser.add_argument('-zcut',  type=float, default=-1.0, help='cut on the shape')
parser.add_argument('-out',   type=str,   default='incl-spect-v-cut.pdf', help='cut on the shape')

parser.add_argument('-vmin',  type=float, default= 0.0, help='min cut on the shape')
parser.add_argument('-vmax',  type=float, default=12.0, help='max cut on the shape')
parser.add_argument('-nv',    type=int,   default=101,  help='number of points')

#parser.add_argument('-subset', type=int, default=0, help='PDF mmeer set')

def central(array):
    return array[0]



args = parser.parse_args()

script_dir = os.path.dirname(__file__)

# read in the cross-sections split in i->f channels
#incl_qq_array=get_array(os.path.join(script_dir,'../../2to2/res/lhc13-CT14nnlo_"+str(args.subset)+"-ymax4.5-xx2qq.res'))
#incl_qg_array=get_array(os.path.join(script_dir,'../../2to2/res/lhc13-CT14nnlo_"+str(args.subset)+"-ymax4.5-qg2qg.res'))
#incl_gg_array=get_array(os.path.join(script_dir,'../../2to2/res/lhc13-CT14nnlo_"+str(args.subset)+"-ymax4.5-xx2gg.res'))

incl_qq_array=get_array(os.path.join(script_dir,'../../2to2/res/lhc13-NNPDF31_nnlo_as_0118_all-ymax4.5-xx2qq.res'))
incl_qg_array=get_array(os.path.join(script_dir,'../../2to2/res/lhc13-NNPDF31_nnlo_as_0118_all-ymax4.5-qg2qg.res'))
incl_gg_array=get_array(os.path.join(script_dir,'../../2to2/res/lhc13-NNPDF31_nnlo_as_0118_all-ymax4.5-xx2gg.res'))

# get the admixtures at the requested pt
incl_qq_born = np.interp(args.pt, incl_qq_array[:,0], incl_qq_array[:,1])
incl_qg_born = np.interp(args.pt, incl_qg_array[:,0], incl_qg_array[:,1])
incl_gg_born = np.interp(args.pt, incl_gg_array[:,0], incl_gg_array[:,1])

# now get the profile
cmdbase=['../get_cut_efficiency','-'+args.shape,'-R','0.4','-pt',str(args.pt),'-alpha',str(args.alpha),'-ktcut',str(args.ktcut),'-zcut',str(args.zcut),'-beta',str(args.beta)]

# we need acceptances for q & g varying the th params
muRs=[1.0,0.5,2.0,1.0,1.0]
muQs=[1.0,1.0,1.0,0.5,2.0]

# try to build a label for the shape
shape_label=''
if args.shape=='ecf':
    shape_label=r'ECF$^{('+str(args.alpha)+')}$'
else:
    shape_label=r'n$_{\mathrm{LH}}^{('+str(args.ktcut)+'\, \mathrm{GeV})}$'
shape_label = shape_label + ', '
if args.zcut>0:
    if args.beta<1e-4 and args.beta>-1e-4:
        shape_label = shape_label + 'mMDT('+str(args.zcut)+')'
    else:
        shape_label = shape_label + 'SD('+str(args.beta)+','+str(args.zcut)+')'
else:
    shape_label = shape_label + 'plain'

# for each cut we'll study acceptances and tagged x-sections
# We want the central scale and the uncertainties taken as the envelope of the scale variations
cuts=np.linspace(args.vmin,args.vmax,args.nv)
eps_q             = {'central':np.empty([args.nv]), 'min':np.empty([args.nv]), 'max':np.empty([args.nv])}
eps_g             = {'central':np.empty([args.nv]), 'min':np.empty([args.nv]), 'max':np.empty([args.nv])}
gg_purity         = {'central':np.empty([args.nv]), 'min':np.empty([args.nv]), 'max':np.empty([args.nv])}
incl_gg_tagged    = {'central':np.empty([args.nv]), 'min':np.empty([args.nv]), 'max':np.empty([args.nv])}
incl_ratio_tagged = {'central':np.empty([args.nv]), 'min':np.empty([args.nv]), 'max':np.empty([args.nv])}

print ('calculating distriutions')
for i,cut in zip(range(0,args.nv),cuts):
    print ('  cut=',str(cut))
    eps_q_v_mu=[]
    eps_g_v_mu=[]
    gg_purity_v_mu=[]
    incl_gg_tagged_v_mu=[]
    incl_ratio_tagged_v_mu=[]
    for muR, muQ in zip(muRs,muQs):
        epstag_q_as_g=float(subprocess.run(cmdbase+['-vcut',str(cut),'-muR',str(muR),'-muQ',str(muQ)],          stdout=subprocess.PIPE).stdout)
        epstag_g_as_g=float(subprocess.run(cmdbase+['-vcut',str(cut),'-muR',str(muR),'-muQ',str(muQ),'-gluon'], stdout=subprocess.PIPE).stdout)
        epstag_q_as_q=1-epstag_q_as_g
        epstag_g_as_q=1-epstag_g_as_g
        eps_q_v_mu.append(epstag_q_as_g)
        eps_g_v_mu.append(epstag_g_as_g)
        qq_tagged = incl_qq_born*epstag_q_as_q*epstag_q_as_q + incl_qg_born*epstag_q_as_q*epstag_g_as_q + incl_gg_born*epstag_g_as_q*epstag_g_as_q
        gg_tagged = incl_qq_born*epstag_q_as_g*epstag_q_as_g + incl_qg_born*epstag_q_as_g*epstag_g_as_g + incl_gg_born*epstag_g_as_g*epstag_g_as_g
        gg_purity_v_mu.append(incl_gg_born*epstag_g_as_g*epstag_g_as_g/gg_tagged)
        incl_gg_tagged_v_mu.append(gg_tagged/(incl_qq_born+incl_qg_born+incl_gg_born))
        incl_ratio_tagged_v_mu.append(gg_tagged/qq_tagged)

    eps_q['central'][i] = eps_q_v_mu[0]
    eps_q['min'][i] = min(eps_q_v_mu)
    eps_q['max'][i] = max(eps_q_v_mu)

    eps_g['central'][i] = eps_g_v_mu[0]
    eps_g['min'][i] = min(eps_g_v_mu)
    eps_g['max'][i] = max(eps_g_v_mu)

    gg_purity['central'][i] = gg_purity_v_mu[0]
    gg_purity['min'][i] = min(gg_purity_v_mu)
    gg_purity['max'][i] = max(gg_purity_v_mu)

    incl_gg_tagged['central'][i] = incl_gg_tagged_v_mu[0]
    incl_gg_tagged['min'][i] = min(incl_gg_tagged_v_mu)
    incl_gg_tagged['max'][i] = max(incl_gg_tagged_v_mu)

    incl_ratio_tagged['central'][i] = incl_ratio_tagged_v_mu[0]
    incl_ratio_tagged['min'][i] = min(incl_ratio_tagged_v_mu)
    incl_ratio_tagged['max'][i] = max(incl_ratio_tagged_v_mu)


# now we can plot things
pdfname = args.out
with PdfPages(pdfname) as pdf:
    fig = plt.figure(figsize=(6,6))

    #----------------------------------------------------------------------
    plt.title('fraction of $gg$-tagged events')
    ax=plt.gca()
    
    #ax=fig.add_subplot(111)
    plt.tick_params(which='both', direction='in')

    plt.xlabel(r'cut')
    plt.xlim(args.vmin,args.vmax)
    plt.ylabel(r'fraction of jets tagged as $gg$')
    plt.yscale('log')
    plt.ylim(1e-4,1)
    
    plt.grid(True,ls=':',lw=0.5)
    
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.95,0.9, shape_label, fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.05,0.11, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.05,0.05, 'Pythia8.230(Monash13)',       fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    
    plt.fill_between(cuts, incl_gg_tagged['min'], incl_gg_tagged['max'], color='blue', alpha=0.4, lw=0)
    plt.plot(cuts, incl_gg_tagged['central'], color='blue', ls='-', lw=2)
    #ax.legend(loc='upper left')
    pdf.savefig(bbox_inches='tight')
    plt.close()

    #----------------------------------------------------------------------
    fig = plt.figure(figsize=(6,6))
    plt.title('purity after $gg$ tagging')
    ax=plt.gca()

    plt.tick_params(which='both', direction='in')

    plt.xlabel(r'cut')
    plt.xlim(args.vmin,args.vmax)
    plt.ylabel(r'purity')
    plt.yscale('linear')
    plt.ylim(0.0,1.0)
    plt.yticks(np.linspace(0.0,1.0,11))
    
    plt.grid(True,ls=':',lw=0.5)
    
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.05,0.95, shape_label, fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.11, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.05, 'Pythia8.230(Monash13)',       fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)

    plt.fill_between(cuts, gg_purity['min'], gg_purity['max'], color='blue', alpha=0.4, lw=0)
    plt.plot(cuts, gg_purity['central'], color='blue', ls='-', lw=2)
    #ax.legend(loc='upper left')
    pdf.savefig(bbox_inches='tight')
    plt.close()

    #----------------------------------------------------------------------
    fig = plt.figure(figsize=(6,6))
    plt.title('relative uncertainty after $gg$ tagging')
    ax=plt.gca()

    plt.tick_params(which='both', direction='in')

    plt.xlabel(r'cut')
    plt.xlim(args.vmin,args.vmax)
    plt.ylabel(r'theory uncertainty')
    plt.yscale('linear')
    plt.ylim(0.0,1.0)
    plt.yticks(np.linspace(0.0,1.0,11))
    
    plt.grid(True,ls=':',lw=0.5)
    
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.05,0.95, shape_label, fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.11, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.05, 'Pythia8.230(Monash13)',       fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)

    plt.plot(cuts, 0.5*(incl_gg_tagged['max']-incl_gg_tagged['min'])/incl_gg_tagged['central'], color='blue', ls='-', lw=2)
    #ax.legend(loc='upper left')
    pdf.savefig(bbox_inches='tight')
    plt.close()

    # #----------------------------------------------------------------------
    # plt.title('ratio of $gg$-tagged over $qq$-tagged rates')
    # ax=plt.gca()
    # 
    # #ax=fig.add_subplot(111)
    # plt.tick_params(which='both', direction='in')
    # 
    # plt.xlabel(r'cut')
    # plt.xlim(args.vmin,args.vmax)
    # plt.ylabel(r'fraction of jets tagged as $gg$')
    # plt.yscale('log')
    # plt.ylim(1e-4,1e4)
    # 
    # plt.grid(True,ls=':',lw=0.5)
    # 
    # plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    # plt.text(0.95,0.9, shape_label, fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    # plt.text(0.05,0.11, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    # plt.text(0.05,0.05, 'Pythia8.230(Monash13)',       fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    # 
    # plt.fill_between(cuts, incl_ratio_tagged['min'], incl_ratio_tagged['max'], color='blue', alpha=0.4, lw=0)
    # plt.plot(cuts, incl_ratio_tagged['central'], color='blue', ls='-', lw=2)
    # #ax.legend(loc='upper left')
    # pdf.savefig(bbox_inches='tight')
    # plt.close()
    # 
    # #----------------------------------------------------------------------
    # fig = plt.figure(figsize=(6,6))
    # plt.title('relative uncertainty for the $gg/qq$-tagged results')
    # ax=plt.gca()
    # 
    # plt.tick_params(which='both', direction='in')
    # 
    # plt.xlabel(r'cut')
    # plt.xlim(args.vmin,args.vmax)
    # plt.ylabel(r'theory uncertainty')
    # plt.yscale('linear')
    # plt.ylim(0.0,1.0)
    # plt.yticks(np.linspace(0.0,1.0,11))
    # 
    # plt.grid(True,ls=':',lw=0.5)
    # 
    # plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    # plt.text(0.05,0.95, shape_label, fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    # plt.text(0.95,0.11, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    # plt.text(0.95,0.05, 'Pythia8.230(Monash13)',       fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    # 
    # print (0.5*(incl_ratio_tagged['max']-incl_ratio_tagged['min'])/incl_ratio_tagged['central'])
    # plt.plot(cuts, 0.5*(incl_ratio_tagged['max']-incl_ratio_tagged['min'])/incl_ratio_tagged['central'], color='blue', ls='-', lw=2)
    # #ax.legend(loc='upper left')
    # pdf.savefig(bbox_inches='tight')
    # plt.close()
    
    
