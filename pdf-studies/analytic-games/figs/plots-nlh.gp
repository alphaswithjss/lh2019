# -*- GNUPlot -*-
reset
set term pdfcairo enhanced size 8.5cm, 9cm
set out 'plots-nlh.pdf'
set colors classic

pt='2000'
beta='0'
zcut='0.1'
ktcut='1.0'

#--------------------------------------------------------------------------------
# distribution
#--------------------------------------------------------------------------------
set xlabel 'n_{LH}'

set key bottom right

set label 1 'p_t=2 TeV, R=0.4'                       right at graph 0.95, 0.57
set label 2 'k_t>'.ktcut.' GeV, {/Symbol b}='.beta.', z_{cut}='.zcut right at graph 0.95, 0.50

set grid

base='< ../test-lh-multiplicity -pt '.pt.' -R 0.4 -zcut '.zcut.' -beta '.beta.' -ktcut '.ktcut.' -cumul '

set label 11 'quarks' at graph 0.05,0.92
set ylabel '{/Symbol e}_q'
plot base.' -muR 1.0 -muQ 1.0' u 1:2 w l dt 1      lc 1 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=1.0',\
     base.' -muR 0.5 -muQ 1.0' u 1:2 w l dt (16,8) lc 3 lw 3 t '{/Symbol m}_R=0.5, {/Symbol m}_Q=1.0',\
     base.' -muR 2.0 -muQ 1.0' u 1:2 w l dt (4,4)  lc 3 lw 3 t '{/Symbol m}_R=2.0, {/Symbol m}_Q=1.0',\
     base.' -muR 1.0 -muQ 0.5' u 1:2 w l dt (16,8) lc 7 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=0.5',\
     base.' -muR 1.0 -muQ 2.0' u 1:2 w l dt (4,4)  lc 7 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=2.0'

set ylabel '{/Symbol e}_g'
set label 11 'gluon' at graph 0.05,0.92
plot base.' -muR 1.0 -muQ 1.0' u 1:3 w l dt 1      lc 1 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=1.0',\
     base.' -muR 0.5 -muQ 1.0' u 1:3 w l dt (16,8) lc 3 lw 3 t '{/Symbol m}_R=0.5, {/Symbol m}_Q=1.0',\
     base.' -muR 2.0 -muQ 1.0' u 1:3 w l dt (4,4)  lc 3 lw 3 t '{/Symbol m}_R=2.0, {/Symbol m}_Q=1.0',\
     base.' -muR 1.0 -muQ 0.5' u 1:3 w l dt (16,8) lc 7 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=0.5',\
     base.' -muR 1.0 -muQ 2.0' u 1:3 w l dt (4,4)  lc 7 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=2.0'

#--------------------------------------------------------------------------------
# ROC curve
#--------------------------------------------------------------------------------
set xlabel '{/Symbol e}_q'
set xrange [0:1]
set xtics 0.1

set ylabel '{/Symbol e}_g'
set yrange [0:1]
set ytics 0.1

set key top left reverse Left
set size square

unset label 11

set label 1 'p_t=2 TeV, R=0.4'                       left at graph 0.05, 0.57
set label 2 left at graph 0.05, 0.50

plot base.' -muR 1.0 -muQ 1.0' u 2:3 w l dt 1      lc 1 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=1.0',\
     base.' -muR 0.5 -muQ 1.0' u 2:3 w l dt (16,8) lc 3 lw 3 t '{/Symbol m}_R=0.5, {/Symbol m}_Q=1.0',\
     base.' -muR 2.0 -muQ 1.0' u 2:3 w l dt (4,4)  lc 3 lw 3 t '{/Symbol m}_R=2.0, {/Symbol m}_Q=1.0',\
     base.' -muR 1.0 -muQ 0.5' u 2:3 w l dt (16,8) lc 7 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=0.5',\
     base.' -muR 1.0 -muQ 2.0' u 2:3 w l dt (4,4)  lc 7 lw 3 t '{/Symbol m}_R=1.0, {/Symbol m}_Q=2.0'
     

set out
