# -*- GNUPlot -*-
reset
set term pdfcairo enhanced size 9cm, 8.5cm
set out 'gg-fractions.pdf'
set colors classic

shape='ecf_0.5'
multi='isd_1.0'

shapelabel='ECF(0.5)'
multilabel='N_{LH}(1 GeV)'

groomers='plain sd mmdt'
groomerlabels="'plain' 'SD(2,0.05)' 'mMDT(0.1)'"

file(flav,pt)='-f ../../taggers-code/res/lhc13-py8.230_M13-'.flav.'2'.flav.'-ptmin'.sprintf("%d",pt).'-full.res '
r=' | rebin.pl -arc 1'

#--------------------------------------------------------------------------------
# for each shape, plot  q/g efficiencies as a function of the cut
#--------------------------------------------------------------------------------

get(flav,pt)='< mergeidx.pl '.file(flav,pt)

set xlabel 'v'
set ylabel '{/Symbol S}(v)'
set yrange [0:1]

set log x
set xrange [1e-6:1]
set format x '10^{%T}'
set key top left reverse Left
set label 11 'solid: qq2qq' at graph 0.05,0.7
set label 12 'dashed: gg2gg' at graph 0.05,0.63
do for [ig=1:words(groomers)]{
    set title '{} {/:Bold '.shapelabel.' - '.word(groomerlabels,ig).'} {}'
    plot get('qq', 500).word(groomers,ig).'_'.shape.r u (exp($2)):4 w l dt 1      lc 1 lw 3 t '500 GeV',\
         get('gg', 500).word(groomers,ig).'_'.shape.r u (exp($2)):4 w l dt (16,8) lc 1 lw 3 not,\
         get('qq',1000).word(groomers,ig).'_'.shape.r u (exp($2)):4 w l dt 1      lc 3 lw 3 t '1 TeV',\
         get('gg',1000).word(groomers,ig).'_'.shape.r u (exp($2)):4 w l dt (16,8) lc 3 lw 3 not,\
         get('qq',2000).word(groomers,ig).'_'.shape.r u (exp($2)):4 w l dt 1      lc 7 lw 3 t '2 TeV',\
         get('gg',2000).word(groomers,ig).'_'.shape.r u (exp($2)):4 w l dt (16,8) lc 7 lw 3 not
}

unset log x
set xrange [0:15]
set format x '%g'

do for [ig=1:words(groomers)]{
    set title '{} {/:Bold '.multilabel.' - '.word(groomerlabels,ig).'} {}'
    plot get('qq', 500).word(groomers,ig).'_'.multi.r u 2:4 w l dt 1      lc 1 lw 3 t '500 GeV',\
         get('gg', 500).word(groomers,ig).'_'.multi.r u 2:4 w l dt (16,8) lc 1 lw 3 not,\
         get('qq',1000).word(groomers,ig).'_'.multi.r u 2:4 w l dt 1      lc 3 lw 3 t '1 TeV',\
         get('gg',1000).word(groomers,ig).'_'.multi.r u 2:4 w l dt (16,8) lc 3 lw 3 not,\
         get('qq',2000).word(groomers,ig).'_'.multi.r u 2:4 w l dt 1      lc 7 lw 3 t '2 TeV',\
         get('gg',2000).word(groomers,ig).'_'.multi.r u 2:4 w l dt (16,8) lc 7 lw 3 not
}

#--------------------------------------------------------------------------------
# try to see what fraction of gg events remain in the sample after a given cut
#--------------------------------------------------------------------------------

get(pt)='< mergeidx.pl '.file('qq',pt).file('gg',pt)

# rough fractions obtained by looking at the 2to2 results
x500qq=0.2
x500qg=0.54
x500gg=0.26

x1000qq=0.33
x1000qg=0.53
x1000gg=0.14

x2000qq=0.55
x2000qg=0.40
x2000gg=0.05

set xlabel 'v_{cut}'
set ylabel 'f_{gluon}'

set log x
set xrange [1e-6:1]
set format x '10^{%T}'
set key top left reverse Left
unset label 11
unset label 12
do for [ig=1:words(groomers)]{
    set title '{} {/:Bold '.shapelabel.' - '.word(groomerlabels,ig).'} {}'
    plot get( 500).word(groomers,ig).'_'.shape.r u (exp($2)):(x500gg *(1-$8)*(1-$8)/(x500gg *(1-$8)*(1-$8)+x500qg *(1-$4)*(1-$8)+x500qq *(1-$4)*(1-$4))) w l dt 1 lc 1 lw 3 t '500 GeV',\
         get(1000).word(groomers,ig).'_'.shape.r u (exp($2)):(x1000gg*(1-$8)*(1-$8)/(x1000gg*(1-$8)*(1-$8)+x1000qg*(1-$4)*(1-$8)+x1000qq*(1-$4)*(1-$4))) w l dt 1 lc 3 lw 3 t '1 TeV',\
         get(2000).word(groomers,ig).'_'.shape.r u (exp($2)):(x2000gg*(1-$8)*(1-$8)/(x2000gg*(1-$8)*(1-$8)+x2000qg*(1-$4)*(1-$8)+x2000qq*(1-$4)*(1-$4))) w l dt 1 lc 7 lw 3 t '2 TeV'
}


unset log x
set xrange [0:15]
set format x '%g'

do for [ig=1:words(groomers)]{
    set title '{} {/:Bold '.multilabel.' - '.word(groomerlabels,ig).'} {}'
    plot get( 500).word(groomers,ig).'_'.multi.r u 2:(x500gg *(1-$8)*(1-$8)/(x500gg *(1-$8)*(1-$8)+x500qg *(1-$4)*(1-$8)+x500qq *(1-$4)*(1-$4))) w l dt 1 lc 1 lw 3 t '500 GeV',\
         get(1000).word(groomers,ig).'_'.multi.r u 2:(x1000gg*(1-$8)*(1-$8)/(x1000gg*(1-$8)*(1-$8)+x1000qg*(1-$4)*(1-$8)+x1000qq*(1-$4)*(1-$4))) w l dt 1 lc 3 lw 3 t '1 TeV',\
         get(2000).word(groomers,ig).'_'.multi.r u 2:(x2000gg*(1-$8)*(1-$8)/(x2000gg*(1-$8)*(1-$8)+x2000qg*(1-$4)*(1-$8)+x2000qq*(1-$4)*(1-$4))) w l dt 1 lc 7 lw 3 t '2 TeV'
}
set out
