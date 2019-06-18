#!/usr/bin/env python3
# 

"""Plot NLO cross-sections for a given config."""

import sys
import os

import matplotlib as mpl
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator,FixedLocator,NullFormatter,FormatStrFormatter
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
import numpy as np
from hfile import *
import argparse

shapes=['ecf', 'isd']
shapelabels = {'ecf': r'ECF$^{(0.5)}$',
              'isd': r'$N_{LH}^{(1\, \mathrm{ GeV})}$'}

groomers=['plain', 'loose', 'tight']
groomerlabels={'plain':'plain', 'loose':'SD(2,0.05)', 'tight':'mMDT(0.1)'}

cuts = {'very_loose': 'very loose',
        'loose' : 'loose',
        'tight': 'tight',
        'very_tight': 'very tight'}

cut_choices=['very_loose', 'loose', 'tight']
#cut_choices=['tight']


#----------------------------------------------------------------------
def search(filehandle, regexp, return_line=False):
    """ looks through the file described by handle until it finds the regexp
    that's mentioned"""
    regexp_compiled = re.compile(regexp)
    while True:
        line = decode(filehandle.readline())
        if (line == "") : break        # empty line = end-of-file
        #if (re.search(regexp,line)) : 
        if (regexp_compiled.search(line)) : 
            if (return_line): return line
            else            : return filehandle
    return None

def decode(line):
    if (isinstance(line,bytes)): return line.decode(default_encoding)
    else: return line
    
#--------------------------------------------------------------------------
def get_array_from_yoda_file(yoda, regexp):
    # if file is a string, assume it's a filename, otherwise a filehandle
    if isinstance(yoda,basestring): yoda = open_any(yoda, 'r')
    if regexp is not None         : search(yoda,regexp)

    lines = []    # temporary store of lines, before conversion to array
    started = False
    while True:
        line = decode(yoda.readline())
        if not line: break        # empty line = end-of-file
        line = line.rstrip().lstrip() # strips trailing blanks, \n; then leading blanks
        if not line: continue
        if re.match('^END YODA', line): break        # yoda signal for end-of-file
        if re.match('^---', line): continue          # some versions have lines like this!
        if not re.match('[-0-9]', line) : continue
        lines.append(line)                        # collect the line

    # do some basic error checking
    if (len(lines) < 1):
        raise Error("Block in get_array had 0 useful lines")

    # now we know the size, transfer the information to a numpy ndarray
    ncol = len(lines[0].split('\t'))                
    num_array = numpy.empty( (len(lines), ncol) )
    for i in range(len(lines)):
        num_array[i,:] = lines[i].split('\t')

    # # rebin to avoid current low stat
    # result = numpy.empty( (int(len(lines)/2), ncol) )
    # for i in range(int(len(lines)/2)):
    #     result[i,0] = num_array[2*i,0]
    #     result[i,1] = num_array[2*i,1]
    #     result[i,2] = num_array[2*i+1,1]
    #     result[i,3] = ((num_array[2*i,1]-num_array[2*i,0])*num_array[2*i,2]+(num_array[2*i+1,1]-num_array[2*i+1,0])*num_array[2*i+1,2])/(num_array[2*i+1,1]-num_array[2*i,0])
    # 
    # return result
    return num_array

#--------------------------------------------------------------------------
def plot_xs_variations(shape, groomer, files, tags, cut, pdf):
    # prepare the plot
    fig, ax = pyplot.subplots(figsize=(4.5,3.6))

    pyplot.title('gg tag acceptance: '+shapelabels[shape]+', '+groomer+' jet')
    ax.tick_params(which='both', direction='in')
    
    ax.set_xlabel('$p_t$ [GeV]')
    ax.set_xscale('log')
    ticks=[100,200,500,1000,2000,5000]
    pyplot.xticks(ticks, [str(tic) for tic in ticks])

    ax.set_ylabel('ratio to '+tags[0])
    ax.set_yscale('linear')
    pyplot.ylim(0.0,2)
        
    ax.grid(linestyle=':', zorder=1)

    pyplot.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90,
                horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)


    pyplot.text(0.05,0.07,cuts[cut]+' cut',  fontsize=10, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    
    first=True
    at_2tev=[]
    for f,tag,colour in zip(files, tags,['red','red','blue','green']):
        incl = get_array_from_yoda_file(f, 'inclusive')
        data = get_array_from_yoda_file(f, groomer+'_'+shape+'_'+cut+'_gg')
        if first:
            reffrac=(data[:,2]/incl[:,2])
        else:
            pyplot.plot(0.5*(data[:,0]+data[:,1]), (data[:,2]/incl[:,2])/reffrac, color=colour, ls='-', label=tag)
            # rough (intergation between 1.75 and 2.5 TeV)
            # Note that the average is taken the "wrong" way just to get rid of fluctuations
            at_2tev.append(np.sum(data[39:44,2]/incl[39:44,2]/reffrac[39:44])/5)
        first=False

    ax.legend(loc='upper left')
    pyplot.tight_layout()
    pdf.savefig()
    pyplot.close()
    return max(np.max(at_2tev), 1/np.min(at_2tev))-1

#--------------------------------------------------------------------------
def plot_uncertainties(uncertainties, pdf):
    # prepare the plot
    ax={}
    pyplot.title('MC uncertainty')
    fig, (ax[shapes[0]],ax[shapes[1]]) = pyplot.subplots(2,1,sharex=True,figsize=(4.5,3.6))

    ax[shapes[1]].set_xlim(-0.5,len(cut_choices)-0.5)
    ax[shapes[1]].set_xticks(range(len(cut_choices)))
    ax[shapes[1]].set_xticklabels([cuts[cut_choice] for cut_choice in cut_choices])
    ax[shapes[1]].set_xlabel('gluon tagging tightness')

    ax['ecf'].set_ylim(0.0,0.5)
    ax['ecf'].set_yticks(np.linspace(0.0,0.5,6))
    ax['isd'].set_ylim(0.0,1.0)
    ax['isd'].set_yticks(np.linspace(0.0,1.0,6))
    
    for shape in shapes:
        ax[shape].tick_params(which='both', direction='in')
        ax[shape].set_ylabel(r'$\delta_{\mathrm{MC}}$')
        
        ax[shape].grid(linestyle=':', zorder=1)

        ax[shape].text(0.95,0.9,shapelabels[shape], fontsize=13, 
                       horizontalalignment='right', verticalalignment='center', transform=ax[shape].transAxes)

        ax[shape].plot(np.linspace(-0.5,len(cut_choices)-0.5,2),[0.05,0.05], ls=':', color='black')
        for groomer,colour,pt in zip(groomers,['red','green','blue'],['.','^','d']):
            ax[shape].plot(range(0,len(cut_choices)), [uncertainties[shape][groomer][cut] for cut in cut_choices], color=colour, ls='-', marker=pt, label=groomerlabels[groomer])
            ax[shape].legend(loc='upper left')
    pyplot.tight_layout()
    pdf.savefig()
    pyplot.close()

#--------------------------------------------------------------------------
def main():
    #----------------------------------------------------------------------
    # parse the command line
    #----------------------------------------------------------------------
    parser = argparse.ArgumentParser()

    parser.add_argument('-level', type=str, default='full', help='full or parton')    
    parser.add_argument('-pdf', action='store_true', help='do PDF variations instead of MC variations')    
    parser.add_argument('-out', type=str, default='plot-variations.pdf', help='output PDF file name')

    args = parser.parse_args()
    
    print("# "+" ".join(sys.argv))
    print ("# args=",args)

    # parse the command line 
    #----------------------------------------------------------------------
    parser = argparse.ArgumentParser()

    # get the data to plot
    #----------------------------------------------------------------------
    if args.pdf:
        candidates=['pythia8230_M13','pythia8230_M13-CT14nnlo', 'pythia8230_M13-MMHT2014nnlo68cl', 'pythia8230_M13-NNPDF31_nnlo_as_0118']
        cand_labels=['Pythia8.230(def)','Pythia8.230(CT14nnlo)', 'Pythia8.230(MMHT2014nnlo)', 'Pythia8.230(NNPDF31nnlo)']
    else:
        candidates=['pythia8230_M13', 'pythia8230_M13', 'herwig7.1.1', 'sherpa2.2.7']
        cand_labels=['Pythia8.230(def)', 'Pythia8.230(def)', 'Herwig7.1.1(def)', 'Sherpa2.2.7(def)']
    files=[]
    tags=[]
    for candidate,tag in zip(candidates, cand_labels):
        f = '../res/lhc13-'+candidate+'-'+args.level+'.yoda'
        if os.path.isfile(f):
            files.append(f)
            tags.append(tag)
    print (files,tags)

    # plot
    #----------------------------------------------------------------------
    print ("plotting to "+args.out)
    with PdfPages(args.out) as pdf:
        all_unc={}
        for shape in shapes:
            all_unc_v={}
            all_unc[shape]=all_unc_v
            for grooming in groomers: 
                all_unc_vg={}
                all_unc_v[grooming]=all_unc_vg
                for cut in cut_choices:
                    all_unc_vg[cut] = plot_xs_variations(shape,grooming,files,tags, cut, pdf)
        print(all_unc)
        plot_uncertainties(all_unc,pdf)

#------------------------------------------------------------------------
if __name__ == '__main__':
    main()
