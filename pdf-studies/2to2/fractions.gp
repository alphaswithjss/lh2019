# -*- GNUPlot -*-
reset
set term pdfcairo enhanced size 9cm,8.5cm
set out 'fractions.pdf'
set colors classic

set xlabel 'p_t'
set xrange [100:6500]
set log x
set xtics add (200, 500, 2000, 5000)

set ylabel 'channel fractions'
set yrange [0:1]

set grid


set label 10 '{/*0.8 {/Symbol \326}s=13 TeV, CT14nlo, anti-k_t(0.4)}' at graph 1.03,0.03 rotate by 90

do for [ymax in "0.5 4.5"]{
    set label 11 '|y|<'.ymax at graph 0.03, 0.07
    f(c)=' -f res/lhc13-ct14nlo-ymax'.ymax.'-'.c.'.res'
    m='< mergeidx.pl "^# pt"'.f('qq2qq').f('qq2gg').f('qg2qg').f('gg2qq').f('gg2gg')

    set label 1 'quark jet' at graph 0.03,0.92
    plot m u 1:( $3/($3+$6+$9+$12+$15)) w l lt 1 lc 1             lw 3 t 'qq2qq',\
         m u 1:( $6/($3+$6+$9+$12+$15)) w l lt 1 lc rgb "#00dd00" lw 3 t 'qq2gg',\
         m u 1:( $9/($3+$6+$9+$12+$15)) w l lt 1 lc 3             lw 3 t 'qg2qg',\
         m u 1:($12/($3+$6+$9+$12+$15)) w l lt 1 lc 4             lw 3 t 'gg2qq',\
         m u 1:($15/($3+$6+$9+$12+$15)) w l lt 1 lc 7             lw 3 t 'gg2gg'

    set label 1 'gluon jet' at graph 0.03,0.92
    plot m u 1:( $2/($2+$5+$8+$11+$14)) w l lt 1 lc 1             lw 3 t 'qq2qq',\
         m u 1:( $5/($2+$5+$8+$11+$14)) w l lt 1 lc rgb "#00dd00" lw 3 t 'qq2gg',\
         m u 1:( $8/($2+$5+$8+$11+$14)) w l lt 1 lc 3             lw 3 t 'qg2qg',\
         m u 1:($11/($2+$5+$8+$11+$14)) w l lt 1 lc 4             lw 3 t 'gg2qq',\
         m u 1:($14/($2+$5+$8+$11+$14)) w l lt 1 lc 7             lw 3 t 'gg2gg'

    set label 1 'final state qq' at graph 0.03,0.92
    plot m u 1:(( $2+ $3)/($2+$3+$11+$12)) w l lt 1 lc 1 lw 3 t 'init qq',\
         m u 1:(($11+$12)/($2+$3+$11+$12)) w l lt 1 lc 3 lw 3 t 'init gg'

    set label 1 'final state gg' at graph 0.03,0.92
    plot m u 1:(( $5+ $6)/($5+$6+$14+$15)) w l lt 1 lc 1 lw 3 t 'init qq',\
         m u 1:(($14+$15)/($5+$6+$14+$15)) w l lt 1 lc 3 lw 3 t 'init gg'

    set label 1 'any initial state' at graph 0.03,0.92
    plot m u 1:(($2+$3+$11+$12)/($2+$3+$5+$6+$8+$9+$11+$12+$14+$15)) w l lt 1 lc 1 lw 3 t 'final qq',\
         m u 1:(($8+$9)        /($2+$3+$5+$6+$8+$9+$11+$12+$14+$15)) w l lt 1 lc 3 lw 3 t 'final qg',\
         m u 1:(($5+$6+$14+$15)/($2+$3+$5+$6+$8+$9+$11+$12+$14+$15)) w l lt 1 lc 7 lw 3 t 'final gg'
}
    
set out
