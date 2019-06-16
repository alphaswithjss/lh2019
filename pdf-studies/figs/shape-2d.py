#! /usr/bin/env python3
#
# For a given shape (and pt), make a plot of the qq qg and gg
# distributyions in the lambda1 lambda2 plane

from hfile import *
import argparse
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm

def plot_density(density, name, xmin, xmax, pdf):
    fig = plt.figure(figsize=(6,8))
    plt.title('density of shape values for '+name)
    ax=fig.add_subplot(111)
    ax.set_xlabel('$v_1$')
    ax.set_ylabel('$v_2$')
    plt.imshow(density, origin='lower', aspect='auto', cmap=cm.seismic, extent=[xmin, xmax, xmin, xmax])
    cbar = plt.colorbar(orientation='horizontal', pad=0.12, label=name+' density')
    pdf.savefig(bbox_inches='tight')
    plt.close()


def plot_density_ratio(weight, xmin, xmax, pdf):
    fig = plt.figure(figsize=(6,8))
    plt.title('weight = $n_{gg}/(n_{qq}+n_{qg}+n_{gg})$')
    ax=fig.add_subplot(111)
    ax.set_xlabel('$v_1$')
    ax.set_ylabel('$v_2$')
    plt.imshow(weight, norm=LogNorm(), origin='lower', aspect='auto', cmap=cm.seismic, vmin=1e-4, vmax=1.0, extent=[xmin, xmax, xmin, xmax])
    cbar = plt.colorbar(orientation='horizontal', pad=0.12, label='gg/all density')
    pdf.savefig(bbox_inches='tight')
    plt.close()

def plot_nevs(cuts, nevs_gg, nevs_tot, pdf):
    fig = plt.figure(figsize=(6,8))
    plt.title('Number of events kept after cut on weight')
    ax=fig.add_subplot(111)
    ax.set_xlabel('cut on the weight')
    ax.set_xscale('log')
    ax.set_ylabel('$N_{\\rm kept}$')
    ax.set_yscale('log')
    plt.ylim(1,11000)
    ax.grid(linestyle=':', zorder=1)
    ax.tick_params(which='both', direction='in')
    plt.plot(cuts, nevs_gg,  color='red',   ls='-', label='gg')
    plt.plot(cuts, nevs_tot, color='black', ls='-', label='all')
    pdf.savefig(bbox_inches='tight')
    plt.close()

def plot_nevs_v_purity(fraction_ratio, nevs_tot, wgt_fraction_ratio, wgt_nevs, pdf):
    fig = plt.figure(figsize=(6,8))
    plt.title('Purity (fraction of gg) after cut on weight')
    ax=fig.add_subplot(111)
    ax.set_xlabel('purity')
    plt.ylim(1,11000)
    plt.xlim(0,0.5)
    ax.set_ylabel('$N_{\\rm kept}$')
    ax.set_yscale('log')
    ax.grid(linestyle=':', zorder=1)
    ax.tick_params(which='both', direction='in')
    plt.plot(fraction_ratio, nevs_tot,  color='red',   ls='-', label='gg')
    plt.scatter(wgt_fraction_ratio, wgt_nevs, color='blue', marker='.', label='weighted')
    pdf.savefig(bbox_inches='tight')
    plt.close()



parser = argparse.ArgumentParser()

parser.add_argument('-pt',    type=float, default=2000, help='pt cut')    
parser.add_argument('-shape', type=str, default='plain_isd_1.0', help='name of the shape to process')
parser.add_argument('-lumi',  type=float, default=150.0, help='integrated lumi in fb^{-1}')
parser.add_argument('-vmin',  type=float, default=-0.5, help='min v to plot')
parser.add_argument('-vmax',  type=float, default=15.5, help='max v to plot')


args = parser.parse_args()

script_dir = os.path.dirname(__file__)


# read in the cross-sections split in i->f channels
incl_qq2qq=get_array(os.path.join(script_dir,'../2to2/res/lhc13-ct14nlo-ymax4.5-qq2qq.res'))
incl_qq2gg=get_array(os.path.join(script_dir,'../2to2/res/lhc13-ct14nlo-ymax4.5-qq2gg.res'))
incl_qg2qg=get_array(os.path.join(script_dir,'../2to2/res/lhc13-ct14nlo-ymax4.5-qg2qg.res'))
incl_gg2qq=get_array(os.path.join(script_dir,'../2to2/res/lhc13-ct14nlo-ymax4.5-gg2qq.res'))
incl_gg2gg=get_array(os.path.join(script_dir,'../2to2/res/lhc13-ct14nlo-ymax4.5-gg2gg.res'))

# convert this infto final-states admixtures
incl_pts=incl_qq2qq[:,0]
incl_sum=incl_qq2qq[:,1]+incl_qq2qq[:,2]+incl_gg2qq[:,1]+incl_gg2qq[:,2]+incl_qg2qg[:,1]+incl_qg2qg[:,2]+incl_qq2gg[:,1]+incl_qq2gg[:,2]+incl_gg2gg[:,1]+incl_gg2gg[:,2]
incl_qq=(incl_qq2qq[:,1]+incl_qq2qq[:,2]+incl_gg2qq[:,1]+incl_gg2qq[:,2])/incl_sum
incl_qg=(incl_qg2qg[:,1]+incl_qg2qg[:,2])/incl_sum
incl_gg=(incl_qq2gg[:,1]+incl_qq2gg[:,2]+incl_gg2gg[:,1]+incl_gg2gg[:,2])/incl_sum

# - number of events in 150 fb^{-1} above the requested pt
lumi=args.lumi*1e+6
incl_sum_above=incl_sum.copy()
incl_sum_above[len(incl_sum_above)-1]=0.5*incl_sum[len(incl_sum_above)-1]*(6500-incl_pts[len(incl_sum_above)-1])
for i in range(len(incl_sum_above)-2,-1,-1):
    incl_sum_above[i]=incl_sum_above[i+1]+0.5*(incl_sum[i]+incl_sum[i+1])*(incl_pts[i+1]-incl_pts[i])
xs_pt=np.interp(args.pt, incl_pts, incl_sum_above)
nev_init=xs_pt*lumi

# find the fractions at the requested pt
frac_qq_born = np.interp(args.pt, incl_pts, incl_qq)
frac_qg_born = np.interp(args.pt, incl_pts, incl_qg)
frac_gg_born = np.interp(args.pt, incl_pts, incl_gg)

# read the shape distributions
q_distrib = get_array(os.path.join(script_dir,'../../taggers-code/res/lhc13-py8.230_M13-qq2qq-ptmin'+str(int(args.pt))+'-full.res'), args.shape)
g_distrib = get_array(os.path.join(script_dir,'../../taggers-code/res/lhc13-py8.230_M13-gg2gg-ptmin'+str(int(args.pt))+'-full.res'), args.shape)

inrange=(q_distrib[:,1]>args.vmin) * (q_distrib[:,1]<args.vmax)

p_q = q_distrib[inrange,3]*(q_distrib[inrange,2]-q_distrib[inrange,0])
p_g = g_distrib[inrange,3]*(g_distrib[inrange,2]-g_distrib[inrange,0])

# now create 2-d distributions
density_qq=np.outer(p_q, p_q)
density_qg=0.5*(np.outer(p_q, p_g)+np.outer(p_g, p_q))
density_gg=np.outer(p_g, p_g)

# fractions of Born  events in each bin
n_qq=frac_qq_born*density_qq
n_qg=frac_qg_born*density_qg
n_gg=frac_gg_born*density_gg
n_tot=n_qq+n_qg+n_gg

# the "optimal" weight
weight=n_gg/(n_tot+1e-50)

# number of gg and total events
#
# note that for total number n_tot cancels
wgt_frac_tot=frac_gg_born
wgt_frac_gg=np.sum(n_gg*weight)

print(wgt_frac_gg,wgt_frac_tot)

# try to see the number of events and purity we get when imposing a
# cut on the weight (i.e. on isolines in the 2-d plane)
fs=np.logspace(-4,0,50)
frac_gg =numpy.asarray([np.sum(n_gg [(weight>f)]) for f in fs])
frac_tot=numpy.asarray([np.sum(n_tot[(weight>f)]) for f in fs])
# for f in :
#     condition=(weight>f)
#     
#     print (f,,np.sum(n_tot[condition]))

pdfname = 'shape-2d.pdf'
print("Creating "+pdfname)
with PdfPages(pdfname) as pdf:
    plot_density(density_qq, 'qq', args.vmin, args.vmax, pdf)
    plot_density(density_qg, 'qg', args.vmin, args.vmax, pdf)
    plot_density(density_gg, 'gg', args.vmin, args.vmax, pdf)

    plot_density_ratio(weight, args.vmin, args.vmax, pdf)

    plot_nevs(fs, nev_init*frac_gg, nev_init*frac_tot, pdf)
    plot_nevs_v_purity(frac_gg/frac_tot, nev_init*frac_tot, wgt_frac_gg/wgt_frac_tot, nev_init*wgt_frac_tot, pdf)
    
