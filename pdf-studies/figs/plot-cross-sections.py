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
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from hfile import *
import argparse

shapelabels = {'ecf': r'ECF$^{(0.5)}$',
              'isd': r'$N_{LH}^{(1\, \mathrm{ GeV})}$'}

cuts = {'very_loose': 'very_loose',
        'loose' : 'loose',
        'tight': 'tight',
        'very_tight': 'very tight'}

cut_choices=['very_loose', 'loose', 'tight']

styles = {'very_loose': 'red', 'loose': 'green', 'tight': 'blue', 'very_tight': 'black', 'total': 'magenta'}

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
    #     result[i,1] = num_array[2*i,2]
    #     result[i,2] = num_array[2*i+1,2]
    #     result[i,3] = ((num_array[2*i,2]-num_array[2*i,0])*num_array[2*i,3]+(num_array[2*i+1,2]-num_array[2*i+1,0])*num_array[2*i+1,3])/(num_array[2*i+1,2]-num_array[2*i,0])
    # 
    # return result
    return num_array

#--------------------------------------------------------------------------
def plot_xs(shape, groomer, files, tags, pdf):
    # prepare the plot
    fig, ax = pyplot.subplots(figsize=(4.5,3.6))

    pyplot.title('cross-section: '+shapelabels[shape]+', '+groomer+' jet')
    ax.tick_params(which='both', direction='in')
    
    ax.set_xlabel('$p_t$ [GeV]')
    ax.set_xscale('log')
    ticks=[100,200,500,1000,2000,5000]
    pyplot.xticks(ticks, [str(tic) for tic in ticks])

    ax.set_ylabel('$d\sigma/dp_t$ [pb/GeV]')
    ax.set_yscale('log')
        
    ax.grid(linestyle=':', zorder=1)

    pyplot.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90,
                horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    ypos=0.93
    for style,tag in zip(['solid', 'dashed', 'dotted'],tags):
        pyplot.text(0.95,ypos,style+':'+tag,  fontsize=8, horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
        ypos=ypos-0.06
        
    first=True
    for f,s in zip(files, ['-','--',':']):
        data = get_array_from_yoda_file(f, 'inclusive')
        pyplot.plot(0.5*(data[:,0]+data[:,1]), data[:,2], color=styles['total'], ls=s, label=('incl.' if first else None))
        first=False
    for cut in cut_choices:
        first=True
        for f,s in zip(files, ['-','--',':']):
            data = get_array_from_yoda_file(f, groomer+'_'+shape+'_'+cut+'_gg')
            pyplot.plot(0.5*(data[:,0]+data[:,1]), data[:,2], color=styles[cut], ls=s, label=(cuts[cut] if first else None))
            first=False
        
    ax.legend(loc='lower left')
    pyplot.tight_layout()
    pdf.savefig()
    pyplot.close()

def plot_xs_fractions(shape, groomer, files, tags, pdf):
    # prepare the plot
    fig, ax = pyplot.subplots(figsize=(4.5,3.6))

    pyplot.title('gg tag acceptance: '+shapelabels[shape]+', '+groomer+' jet')
    ax.tick_params(which='both', direction='in')
    
    ax.set_xlabel('$p_t$ [GeV]')
    ax.set_xscale('log')
    ticks=[100,200,500,1000,2000,5000]
    pyplot.xticks(ticks, [str(tic) for tic in ticks])

    ax.set_ylabel('ratio to inclusive')
    ax.set_yscale('log')
    pyplot.ylim(0.001,1)
        
    ax.grid(linestyle=':', zorder=1)

    pyplot.text(1.015,0.01,r'$\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4), |y|<4.5', fontsize=8, rotation=90,
                horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)

    ypos=0.56
    for style,tag in zip(['solid', 'dashed', 'dotted'],tags):
        pyplot.text(0.05,ypos,style+':'+tag,  fontsize=8, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        ypos=ypos-0.06
        
    first=True
    incldata=[]
    for f,s in zip(files, ['-','--',':']):
        data = get_array_from_yoda_file(f, 'inclusive')
        incldata.append(data[:,2])
        first=False
    for cut in cut_choices:
        first=True
        for f,idata,s in zip(files, incldata, ['-','--',':']):
            data = get_array_from_yoda_file(f, groomer+'_'+shape+'_'+cut+'_gg')
            pyplot.plot(0.5*(data[:,0]+data[:,1]), data[:,2]/idata, color=styles[cut], ls=s, label=(cuts[cut] if first else None))
            first=False



    # f_py='../res/lhc13-pythia8230_M13-full.yoda'
    # f_hw='../res/lhc13-herwig7.1.1-full.yoda'
    # f_sp='../res/lhc13-sherpa2.2.7-parton.yoda'
    # incl_py = get_array_from_yoda_file(f_py, 'inclusive')
    # incl_hw = get_array_from_yoda_file(f_hw, 'inclusive')
    # incl_sp = get_array_from_yoda_file(f_sp, 'inclusive')
    # for cut in cut_choices:
    #     data = get_array_from_yoda_file(f_py, groomer+'_'+shape+'_'+cut+'_gg')
    #     pyplot.plot(0.5*(data[:,0]+data[:,1]), data[:,2]/incl_py[:,2], color=styles[cut], ls='-', label=cuts[cut])
    #     data = get_array_from_yoda_file(f_hw, groomer+'_'+shape+'_'+cut+'_gg')
    #     pyplot.plot(0.5*(data[:,0]+data[:,1]), data[:,2]/incl_hw[:,2], color=styles[cut], ls='--')
    #     data = get_array_from_yoda_file(f_sp, groomer+'_'+shape+'_'+cut+'_gg')
    #     pyplot.plot(0.5*(data[:,0]+data[:,1]), data[:,2]/incl_sp[:,2], color=styles[cut], ls=':')

    ax.legend(loc='lower left')
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
    parser.add_argument('-out', type=str, default='plot-cross-sections.pdf', help='output PDF file name')

    args = parser.parse_args()
    
    print("# "+" ".join(sys.argv))
    print ("# args=",args)

    # parse the command line 
    #----------------------------------------------------------------------
    parser = argparse.ArgumentParser()

    # get the data to plot
    #----------------------------------------------------------------------
    if args.pdf:
        candidates=['pythia8230_M13-CT14nnlo', 'pythia8230_M13-MMHT2014nnlo68cl', 'pythia8230_M13-NNPDF31_nnlo_as_0118']
        cand_labels=['CT14nnlo', 'MMHT2014nnlo', 'NNPDF31nnlo']
    else:
        candidates=['pythia8230_M13', 'herwig7.1.1', 'sherpa2.2.7']
        cand_labels=['Pythia8.230', 'Herwig7.1.1', 'Sherpa2.2.7']
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
        for shape in ['ecf', 'isd']:
            for grooming in ['plain', 'loose', 'tight']:
                plot_xs(shape,grooming,files,tags, pdf)
        for shape in ['ecf', 'isd']:
            for grooming in ['plain', 'loose', 'tight']:
                plot_xs_fractions(shape,grooming,files,tags, pdf)


#------------------------------------------------------------------------
if __name__ == '__main__':
    main()
