#! /usr/bin/env python3
#
# For a given shape (and pt), make a plot of the qq qg and gg
# distributyions in the lambda1 lambda2 plane

from hfile import *
import argparse
import numpy as np
import os
import math

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm


parser = argparse.ArgumentParser()

parser.add_argument('-pt',  type=float, default=2000, help='pt cut')    
parser.add_argument('-pdf', type=str, default='NNPDF31_nnlo_as_0118', help='name of the shape to process')
parser.add_argument('-lumi',  type=float, default=300.0, help='integrated lumi in fb^{-1}')
parser.add_argument('-epsg',  type=float, default=0.42, help='gluon tagging efficiency')


args = parser.parse_args()

script_dir = os.path.dirname(__file__)

# read in the cross-sections split in final-state channels
born_array={}
born_array['qq']=get_array(os.path.join(script_dir,'../res/lhc13-'+args.pdf+'_all-ymax4.5-xx2qq.res'))
born_array['qg']=get_array(os.path.join(script_dir,'../res/lhc13-'+args.pdf+'_all-ymax4.5-qg2qg.res'))
born_array['gg']=get_array(os.path.join(script_dir,'../res/lhc13-'+args.pdf+'_all-ymax4.5-xx2gg.res'))
channels={'qq', 'qg', 'gg'}

nsets=born_array['qq'].shape[1]
pts=born_array['qq'][:,0]

born_xs={}        # all BBorn-level x-sections at pt
born_xs_avg={}    # average
born_xs_delta={}  # RELATIVE dispersion
print ('raw cross-sections')
print ('channel   average    rel_unc')
for c in channels:
    born_xs[c] = numpy.asarray([ np.interp(args.pt, pts, born_array[c][:,i]) for i in range(1,nsets)])
    born_xs_avg[c] = born_xs[c][0]
    born_xs_delta[c] = math.sqrt(np.sum((born_xs[c][1:]-born_xs_avg[c])*(born_xs[c][1:]-born_xs_avg[c]))/nsets)/born_xs_avg[c]
    print (f'  {c:7}{born_xs_avg[c]:.5g} {born_xs_delta[c]:9.6f}')
print ('----------------------------------------')

cmapdict={'red' :  ((0.0,  1.0, 1.0),
                    (0.25, 0.0, 0.0),
                    (0.5,  0.0, 0.0),
                    (0.75, 1.0, 1.0),
                    (1.0 , 1.0, 1.0)),
          'green': ((0.0,  1.0, 1.0),
                    (0.25, 0.0, 0.0),
                    (0.5,  1.0, 1.0),
                    (0.75, 1.0, 1.0),
                    (1.0 , 0.0, 0.0)),
          'blue':  ((0.0,  1.0, 1.0),
                    (0.25, 1.0, 1.0),
                    (0.5,  0.0, 0.0),
                    (0.75, 0.0, 0.0),
                    (1.0 , 0.0, 0.0))}
my_cmap = LinearSegmentedColormap('user-defined', cmapdict)


pdfname = 'performance-plots.pdf'
with PdfPages(pdfname) as pdf:
    # plot the (relative) PDF uncertainty as a function of purities
    fig = plt.figure(figsize=(6,8))
    plt.title('PDF uncertainty v (ideal) tagging purity')
    ax=fig.add_subplot(111)
    ax.set_xlabel('$gg$ purity')
    ax.set_ylabel('$qq$ contamination')
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.95,0.92,r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.85,args.pdf, fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.imshow([[math.sqrt((pqq*born_xs_delta['qq'])**2+((1-pqq-pgg)*born_xs_delta['qg'])**2+(pgg*born_xs_delta['gg'])**2) if pqq+pgg<1 else 0 for pgg in np.linspace(0,1,100) ] for pqq in np.linspace(0,1,100)], origin='lower', aspect='auto', vmin=0.02, vmax=0.14, cmap=my_cmap, extent=[0, 1, 0, 1])
    cbar = plt.colorbar(orientation='horizontal', pad=0.12, label='PDF uncertainty')
    pdf.savefig(bbox_inches='tight')
    plt.close()

    #now assume factorisation
    pgg=np.linspace(0.0,1.0,100)
    epsq_over_epsg=0.5*(np.sqrt((born_xs_avg['qg']/born_xs_avg['qq'])**2+4*born_xs_avg['gg']/born_xs_avg['qq']*(1/pgg-1))-born_xs_avg['qg']/born_xs_avg['qq'])
    delta_pdf=pgg*np.sqrt(born_xs_delta['gg']**2 + (epsq_over_epsg*born_xs_avg['qg']/born_xs_avg['gg']*born_xs_delta['qg'])**2  + (epsq_over_epsg**2*born_xs_avg['qq']/born_xs_avg['gg']*born_xs_delta['qq'])**2)

    # number of events per lumi*espg^2
    #
    # we need the cross-section above pt
    born_above={}
    for c in channels:
        b=born_array[c][:,1]
        born_above_array=b.copy()
        born_above_array[len(b)-1]=0.5*b[len(b)-1]*(6500-pts[len(b)-1])
        for i in range(len(b)-2,-1,-1):
            born_above_array[i]=born_above_array[i+1]+0.5*(b[i]+b[i+1])*(pts[i+1]-pts[i])
        born_above[c]=np.interp(args.pt, pts, born_above_array)

    nev_per_lepsg2=born_above['gg']+epsq_over_epsg*born_above['qg']+epsq_over_epsg*epsq_over_epsg*born_above['qq']
    delta_stat=1/(args.epsg*np.sqrt(args.lumi*1e6*nev_per_lepsg2))
    
    fig = plt.figure(figsize=(6,6))
    plt.title('PDF uncertainty v (ideal) tagging purity')
    ax=fig.add_subplot(111)
    ax.tick_params(which='both', direction='in')
    ax.set_xlabel('$gg$ purity')
    ax.set_ylabel('PDF uncertainty')
    plt.xlim(0,1)
    plt.ylim(0,0.14)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,0.14,15))
    plt.grid(True,ls='--',lw=0.5)
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.95,0.16, r'${\cal L}$='+f'{args.lumi:g}'+r' fb$^{-1}$, $\varepsilon_g=$'+f'{args.epsg:g}', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.10, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.04, args.pdf,                      fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.plot(pgg, delta_pdf, color='blue', ls='-', label='PDF uncertainty', lw=2)
    plt.plot(pgg, 0.05+0.0*pgg, color='red', ls=':', label='syst. uncert', lw=2)
    plt.plot(pgg, delta_stat, color='green', ls='--', label='stat. uncert', lw=2)
    ax.legend(loc='upper left')
    pdf.savefig(bbox_inches='tight')
    plt.close()
    
    # try to see what quark mistag rate this all corresponds to
    fig = plt.figure(figsize=(6,6))
    plt.title('quark mistag rate')
    ax=fig.add_subplot(111)
    ax.tick_params(which='both', direction='in')
    ax.set_xlabel('$gg$ purity')
    ax.set_ylabel(r'$\varepsilon_q$')
    plt.xlim(0,1)
    plt.ylim(0,1.0)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,1.0,11))
    plt.grid(True,ls='--',lw=0.5)
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.95,0.95, r'${\cal L}$='+f'{args.lumi:g}'+r' fb$^{-1}$, $\varepsilon_g=$'+f'{args.epsg:g}', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.89, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.95,0.83, args.pdf,                      fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.plot(pgg, epsq_over_epsg*args.epsg, color='blue', ls='-', label='$q$ mistag rate', lw=2)
    pdf.savefig(bbox_inches='tight')
    plt.close()
    
    # now we can see how taggers perform
    # We study several cases:
    #  - Casimir scaling,
    #  - ECF
   
    taggers=[]
    efficiencies={}
    colours={}

    # Casimir scaling: eps_g=eps_g^(CA/CF)
    tagger_name=r'Casimir: $(1-\varepsilon_g)^{C_F}=(1-\varepsilon_q)^{C_A}$'
    taggers.append(tagger_name)
    colours[tagger_name] = 'red'
    efficiencies[tagger_name]={'g' : np.linspace(0.0, 1.0, 101)}
    efficiencies[tagger_name]['q'] =  1.0-np.power(1.0-efficiencies[tagger_name]['g'], 4.0/9.0)

    # ECF
    tagger_name=r'ECF$^{(0.5)}$ [LH15]'
    taggers.append(tagger_name)
    colours[tagger_name] = 'green'
    efficiencies[tagger_name]={}
    for f in ['q', 'g']:
        distrib=get_array('../../../taggers-code/res/lhc13-py8.230_M13-'+f+f+'2'+f+f+'-ptmin'+str(args.pt)+'-full.res', 'plain_ecf_0.5')
        #eps_v_cut[f]=np.max(1.0-numpy.cumsum(distrib[:,3]*(distrib[:,2]-distrib[:,0])), 0.0) # avoid rounding close to enf of spectrum
        efficiencies[tagger_name][f]=1.0-numpy.cumsum(distrib[:,3]*(distrib[:,2]-distrib[:,0]))
        efficiencies[tagger_name][f] = numpy.insert(efficiencies[tagger_name][f],0, 1.0)

    # ISD
    tagger_name=r'SD$_{(\beta=2,z_{\mathrm{cut}}=0.05)}$ + $n_{LH}^{(1\,\mathrm{ GeV})}$ [LH17,19]'
    taggers.append(tagger_name)
    colours[tagger_name] = 'blue'
    efficiencies[tagger_name]={}
    for f in ['q', 'g']:
        distrib=get_array('../../../taggers-code/res/lhc13-py8.230_M13-'+f+f+'2'+f+f+'-ptmin'+str(args.pt)+'-full.res', 'sd_isd_1.0')
        #eps_v_cut[f]=np.max(1.0-numpy.cumsum(distrib[:,3]*(distrib[:,2]-distrib[:,0])), 0.0) # avoid rounding close to enf of spectrum
        efficiencies[tagger_name][f]=1.0-numpy.cumsum(distrib[:,3]*(distrib[:,2]-distrib[:,0]))
        efficiencies[tagger_name][f] = numpy.insert(efficiencies[tagger_name][f],0, 1.0)

    # NN
    tagger_name=r'Neural Network-like'
    taggers.append(tagger_name)
    colours[tagger_name] = 'black'
    efficiencies[tagger_name]={'q' : np.linspace(0.0, 1.0, 101)}
    efficiencies[tagger_name]['g'] =  1.0-(np.exp(8.0*(1.0-efficiencies[tagger_name]['q']))-1)/(math.exp(8.0)-1)

    # once we have the efficiencies we can compute all the rest: purity, pdf uncert, stat uncert, ...
    stat_uncerts={}
    pdf_uncerts={}
    for tagger in taggers:
        nev_qq = args.lumi*1e6*born_above['qq']*efficiencies[tagger]['q']*efficiencies[tagger]['q']
        nev_qg = args.lumi*1e6*born_above['qg']*efficiencies[tagger]['q']*efficiencies[tagger]['g']
        nev_gg = args.lumi*1e6*born_above['gg']*efficiencies[tagger]['g']*efficiencies[tagger]['g']
        nev_tot = nev_qq + nev_qg + nev_gg
        frac_qq = nev_qq/nev_tot
        frac_qg = nev_qg/nev_tot
        frac_gg = nev_gg/nev_tot
        pdf_uncerts[tagger]=np.sqrt((frac_gg*born_xs_delta['gg'])**2 + (frac_qg*born_xs_delta['qg'])**2  + (frac_qq*born_xs_delta['qq'])**2)
        stat_uncerts[tagger]=1.0/numpy.sqrt(nev_tot)

        
    # plot ROC curves
    fig = plt.figure(figsize=(6,6))
    plt.title('ROC curves for a few selected taggers')
    ax=fig.add_subplot(111)
    ax.tick_params(which='both', direction='in')
    ax.set_xlabel(r'$\varepsilon_g$')
    ax.set_ylabel(r'$\varepsilon_q$')
    plt.xlim(0,1)
    plt.ylim(0,1.0)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,1.0,11))
    plt.grid(True,ls='--',lw=0.5)
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.05,0.59, r'$p_t$='+str(args.pt)+' GeV', fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.05,0.53, 'Pythia8.230(Monash13)',       fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    for tagger in taggers:
        plt.plot(efficiencies[tagger]['g'], efficiencies[tagger]['q'], color=colours[tagger], ls='-', label=tagger, lw=2)
    ax.legend(loc='upper left')
    pdf.savefig(bbox_inches='tight')
    plt.close()
        
    # plot uncertainties
    fig = plt.figure(figsize=(6,6))
    plt.title('ROC curves for a few selected taggers')
    ax=fig.add_subplot(111)
    ax.tick_params(which='both', direction='in')
    ax.set_xlabel(r'$\varepsilon_g$')
    ax.set_ylabel(r'$\varepsilon_q$')
    plt.xlim(0,1)
    plt.ylim(0,0.1)
    plt.xticks(np.linspace(0,1,11))
    plt.yticks(np.linspace(0,plt.ylim()[1],11))
    plt.grid(True,ls='--',lw=0.5)
    plt.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.text(0.04,0.17, r'$p_t$='+str(args.pt)+' GeV, ${\cal L}$='+f'{args.lumi:g}'+r' fb$^{-1}$', fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.04,0.115, args.pdf,                                                                  fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.04,0.06, 'Pythia8.230(Monash13)',                                                   fontsize=13, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    plt.text(0.865,0.677, r'solid: $\delta_{\mathrm{PDF}}$',   fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.865,0.622, r'dashed: $\delta_{\mathrm{stat}}$', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.865,0.567, r'dotted: $\delta_{\mathrm{syst}}$', fontsize=13, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)

    plt.plot(np.linspace(0,1,2), 0.05+0.0*np.linspace(0,1,2), color='black', ls=':', lw=2)
    for tagger in taggers:
        plt.plot(efficiencies[tagger]['g'], pdf_uncerts [tagger], color=colours[tagger], ls='-', label=tagger, lw=2)
        plt.plot(efficiencies[tagger]['g'], stat_uncerts[tagger], color=colours[tagger], ls='--', lw=2)
    ax.legend(loc='upper right')
    pdf.savefig(bbox_inches='tight')
    plt.close()
    
