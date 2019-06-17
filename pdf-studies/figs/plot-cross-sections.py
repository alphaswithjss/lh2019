#!/usr/bin/env python3
# 

"""Plot NLO cross-sections for a given config."""

import sys

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
              'isd': r'$N_{LH}^{(1 GeV)}$'}

cuts = {'very_loose': 'very_loose',
        'loose' : 'loose',
        'tight': 'tight',
        'very_tight': 'very tight'}

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
    return num_array

#--------------------------------------------------------------------------
def plot_xs(shape, groomer, pdf):
    # prepare the plot
    fig, ax = pyplot.subplots(figsize=(4.5,3.6))

    pyplot.title('cross-section for gg-tagged jets')
    ax.tick_params(which='both', direction='in')
    
    ax.set_xlabel('$p_t$ [GeV]')
    ax.set_xscale('log')
    ticks=[100,200,500,1000,2000,5000]
    pyplot.xticks(ticks, [str(tic) for tic in ticks])

    ax.set_ylabel('$d\sigma/dp_t$ [pb/GeV]')
    ax.set_yscale('log')
        
    ax.grid(linestyle=':', zorder=1)

    pyplot.text(0.05,0.06,r'$pp$, $\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4)', fontsize=8,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    pyplot.text(0.05,0.20,r'solid: Pythia8.230', fontsize=8,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    pyplot.text(0.05,0.13,r'dashed: Herwig7.1.1', fontsize=8,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    pyplot.text(0.05,0.6, shapelabels[shape], fontsize=11,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    pyplot.text(0.05,0.51, groomer+' jet', fontsize=11,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    f_py='../res/lhc13-pythia8230_M13-full.yoda'
    f_hw='../res/lhc13-herwig7.1.1-full.yoda'
    data = get_array_from_yoda_file(f_py, 'inclusive')
    pyplot.plot(data[:,1], data[:,3], color=styles['total'], ls='-', label='incl.')
    data = get_array_from_yoda_file(f_hw, 'inclusive')
    pyplot.plot(data[:,1], data[:,3], color=styles['total'], ls=':')
    for cut in ['very_loose', 'loose', 'tight', 'very_tight']:
        data = get_array_from_yoda_file(f_py, groomer+'_'+shape+'_'+cut+'_gg')
        pyplot.plot(data[:,1], data[:,3], color=styles[cut], ls='-', label=cuts[cut])
        data = get_array_from_yoda_file(f_hw, groomer+'_'+shape+'_'+cut+'_gg')
        pyplot.plot(data[:,1], data[:,3], color=styles[cut], ls=':')
        
    ax.legend(loc='upper right')
    pyplot.tight_layout()
    pdf.savefig()
    pyplot.close()

def plot_xs_fractions(shape, groomer, pdf):
    # prepare the plot
    fig, ax = pyplot.subplots(figsize=(4.5,3.6))

    pyplot.title('cross-section for different flavour channels')
    ax.tick_params(which='both', direction='in')
    
    ax.set_xlabel('$p_t$ [GeV]')
    ax.set_xscale('log')
    ticks=[100,200,500,1000,2000,5000]
    pyplot.xticks(ticks, [str(tic) for tic in ticks])

    ax.set_ylabel('ratio to inclusive')
    ax.set_yscale('linear')
    pyplot.ylim(0,1)
        
    pyplot.text(0.05,0.06,r'$pp$, $\sqrt{s}$ = 13 TeV, anti-$k_t$(0.4)', fontsize=8,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    pyplot.text(0.05,0.20,r'solid: Pythia8.230', fontsize=8,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    pyplot.text(0.05,0.13,r'dashed: Herwig7.1.1', fontsize=8,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    pyplot.text(0.05,0.6, shapelabels[shape], fontsize=11,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    pyplot.text(0.05,0.51, groomer+' jet', fontsize=11,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    ax.grid(linestyle=':', zorder=1)

    f_py='../res/lhc13-pythia8230_M13-full.yoda'
    f_hw='../res/lhc13-herwig7.1.1-full.yoda'
    incl_py = get_array_from_yoda_file(f_py, 'inclusive')
    incl_hw = get_array_from_yoda_file(f_hw, 'inclusive')
    for cut in ['very_loose', 'loose', 'tight', 'very_tight']:
        data = get_array_from_yoda_file(f_py, groomer+'_'+shape+'_'+cut+'_gg')
        pyplot.plot(data[:,1], data[:,3]/incl_py[:,3], color=styles[cut], ls='-', label=cuts[cut])
        data = get_array_from_yoda_file(f_hw, groomer+'_'+shape+'_'+cut+'_gg')
        pyplot.plot(data[:,1], data[:,3]/incl_hw[:,3], color=styles[cut], ls=':')

    ax.legend(loc='upper right')
    pyplot.tight_layout()
    pdf.savefig()
    pyplot.close()


#--------------------------------------------------------------------------
def main():
    #----------------------------------------------------------------------
    # parse the command line
    #----------------------------------------------------------------------
    parser = argparse.ArgumentParser()

    #parser.add_argument('-shape', type=str, default='loose_isd', help='choice of shape {plain,loose,tight}_{ecf,isd}')    
    parser.add_argument('-out', type=str, default='plot-cross-sections.pdf', help='output PDF file name')

    args = parser.parse_args()
    
    print("# "+" ".join(sys.argv))
    print ("# args=",args)

    # parse the command line 
    #----------------------------------------------------------------------
    parser = argparse.ArgumentParser()

    # get the data to plot
    #----------------------------------------------------------------------

    # plot
    #----------------------------------------------------------------------
    print ("plotting to "+args.out)
    with PdfPages(args.out) as pdf:
        for shape in ['ecf', 'isd']:
            for grooming in ['plain', 'loose', 'tight']:
                plot_xs(shape,grooming, pdf)
        for shape in ['ecf', 'isd']:
            for grooming in ['plain', 'loose', 'tight']:
                plot_xs_fractions(shape,grooming, pdf)


#------------------------------------------------------------------------
if __name__ == '__main__':
    main()
