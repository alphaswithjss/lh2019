# -*- GNUPlot -*-
reset
set term pdfcairo enhanced
set out 'cut-effects.pdf'
set colors classic

set xlabel 'fraction of gg events in gg-tagged sample'
set xrange [0:0.5]

set grid

groomers='plain sd mmdt'
groomerlabels="'plain' 'SD(2,0.05)' 'mMDT(0.1)'"

set style line 1 dt 1 lc 1             lw 3
set style line 2 dt 1 lc rgb "#00dd00" lw 3
set style line 3 dt 1 lc 3             lw 3

set ylabel 'cut on ECF(0.5)'
set yrange [*:*]
set key top left reverse Left
plot for [ig=1:words(groomers)] '< ../get-cuts-and-stats.py -lumi 150 -log -fstep 0.01 -shape '.word(groomers,ig).'_ecf_0.5' u 1:2 w l ls ig t word(groomerlabels,ig)

set ylabel 'number of events left after cut on ECF(0.5)'
set key top right noreverse Right
set log y
set yrange [1:11000]
plot for [ig=1:words(groomers)] '< ../get-cuts-and-stats.py -lumi 150 -log -fstep 0.01 -shape '.word(groomers,ig).'_ecf_0.5' u 1:4 w l ls ig t word(groomerlabels,ig)
unset log y

set ylabel 'cut on N_{LH}(1 GeV)'
set yrange [*:*]
set key top left reverse Left
plot for [ig=1:words(groomers)] '< ../get-cuts-and-stats.py -lumi 150 -fstep 0.01 -shape '.word(groomers,ig).'_isd_1.0' u 1:2 w l ls ig t word(groomerlabels,ig)

set ylabel 'number of events left after cut on N_{LH}(1 GeV)'
set key top right noreverse Right
set log y
set yrange [1:11000]
plot for [ig=1:words(groomers)] '< ../get-cuts-and-stats.py -lumi 150 -fstep 0.01 -shape '.word(groomers,ig).'_isd_1.0' u 1:4 w l ls ig t word(groomerlabels,ig)
unset log y

set ylabel 'cut on N_{LH}(2 GeV)'
set yrange [*:*]
set key top left reverse Left
plot for [ig=1:words(groomers)] '< ../get-cuts-and-stats.py -lumi 150 -fstep 0.01 -shape '.word(groomers,ig).'_isd_2.0' u 1:2 w l ls ig t word(groomerlabels,ig)

set ylabel 'number of events left after cut on N_{LH}(2 GeV)'
set key top right noreverse Right
set log y
set yrange [1:11000]
plot for [ig=1:words(groomers)] '< ../get-cuts-and-stats.py -lumi 150 -fstep 0.01 -shape '.word(groomers,ig).'_isd_2.0' u 1:4 w l ls ig t word(groomerlabels,ig)
unset log y


set out
