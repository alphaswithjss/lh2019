# -*- GNUPlot -*-
reset
set term pdfcairo enhanced size 9cm, 8.5cm
set out 'plots-distribs.pdf'
set colors classic

#--------------------------------------------------------------------------------
# plot distributions
#--------------------------------------------------------------------------------

set xlabel 'v'
set ylabel '1/N dN/dv'

vs='plain_ang_0.5 plain_ang_1.0 plain_ang_2.0 sd_ang_0.5 sd_ang_1.0 sd_ang_2.0 mmdt_ang_0.5 mmdt_ang_1.0 mmdt_ang_2.0 plain_ecf_0.5 plain_ecf_1.0 plain_ecf_2.0 sd_ecf_0.5 sd_ecf_1.0 sd_ecf_2.0 mmdt_ecf_0.5 mmdt_ecf_1.0 mmdt_ecf_2.0'
labs="'Angularity(0.5), plain' 'Angularity(1.0), plain' 'Angularity(2.0), plain' 'Angularity(0.5), SD(2,0.05)' 'Angularity(1.0), SD(2,0.05)' 'Angularity(2.0), SD(2,0.05)' 'Angularity(0.5), mMDT(0.1)' 'Angularity(1.0), mMDT(0.1)' 'Angularity(2.0), mMDT(0.1)' 'ECF(0.5), plain' 'ECF(1.0), plain' 'ECF(2.0), plain' 'ECF(0.5), SD(2,0.05)' 'ECF(1.0), SD(2,0.05)' 'ECF(2.0), SD(2,0.05)' 'ECF(0.5), mMDT(0.1)' 'ECF(1.0), mMDT(0.1)' 'ECF(2.0), mMDT(0.1)'"

set log x
set xrange [1e-6:1]
set format x '10^{%T}'
set key top left reverse Left
do for [iv=1:words(vs)]{
    set title '{} {/:Bold '.word(labs,iv).'} {}'
    plot '< mergeidx.pl -f ../res/lhc13-py8.230_M13-qq2qq-ptmin500-full.res '.word(vs,iv) u (exp($2)):4 w l dt 1 lc 1 lw 3 t 'qq2qq',\
         '< mergeidx.pl -f ../res/lhc13-py8.230_M13-gg2gg-ptmin500-full.res '.word(vs,iv) u (exp($2)):4 w l dt 1 lc 3 lw 3 t 'gg2gg'
}

vs='plain_isd_1.0 plain_isd_2.0 plain_isd_5.0 sd_isd_1.0 sd_isd_2.0 sd_isd_5.0 mmdt_isd_1.0 mmdt_isd_2.0 mmdt_isd_5.0'
labs="'ISD(1.0), plain' 'ISD(2.0), plain' 'ISD(5.0), plain' 'ISD(1.0), SD(2,0.05)' 'ISD(2.0), SD(2,0.05)' 'ISD(5.0), SD(2,0.05)' 'ISD(1.0), mMDT(0.1)' 'ISD(2.0), mMDT(0.1)' 'ISD(5.0), mMDT(0.1)'"

unset log x
set xrange [0:15]
set format x '%g'
set key top right noreverse Right

do for [iv=1:words(vs)]{
    set title '{} {/:Bold '.word(labs,iv).'} {}'
    plot '< mergeidx.pl -f ../res/lhc13-py8.230_M13-qq2qq-ptmin500-full.res '.word(vs,iv) u 2:4 w l dt 1 lc 1 lw 3 t 'qq2qq',\
         '< mergeidx.pl -f ../res/lhc13-py8.230_M13-gg2gg-ptmin500-full.res '.word(vs,iv) u 2:4 w l dt 1 lc 3 lw 3 t 'gg2gg'
}

#--------------------------------------------------------------------------------
# ROC curves
#--------------------------------------------------------------------------------

set style line 1 dt 1 lc 1 lw 3
set style line 2 dt 1 lc 3 lw 3
set style line 3 dt 1 lc 7 lw 3
set style line 4 dt (16,8) lc 1 lw 3
set style line 5 dt (16,8) lc 3 lw 3
set style line 6 dt (16,8) lc 7 lw 3
set style line 7 dt (4,4) lc 1 lw 3
set style line 8 dt (4,4) lc 3 lw 3
set style line 9 dt (4,4) lc 7 lw 3

set grid
set key top left reverse Left

set xrange [0:1]
set xlabel '{/Symbol e}_q'
set yrange [0:1]
set ylabel '{/Symbol e}_g'

# fix observable, vary the rest
m='< mergeidx.pl -f ../res/lhc13-py8.230_M13-qq2qq-ptmin500-full.res -f ../res/lhc13-py8.230_M13-gg2gg-ptmin500-full.res '
r=' | rebin.pl -arc 1'

set title '{} {/:Bold Angularities} {}'
plot m.'plain_ang_0.5'.r u 4:8 w l ls 1 t 'plain, {/Symbol a}=0.5',\
     m.'plain_ang_1.0'.r u 4:8 w l ls 2 t 'plain, {/Symbol a}=1.0',\
     m.'plain_ang_2.0'.r u 4:8 w l ls 3 t 'plain, {/Symbol a}=2.0',\
     m.'sd_ang_0.5'   .r u 4:8 w l ls 4 t 'SD, {/Symbol a}=0.5',\
     m.'sd_ang_1.0'   .r u 4:8 w l ls 5 t 'SD, {/Symbol a}=1.0',\
     m.'sd_ang_2.0'   .r u 4:8 w l ls 6 t 'SD, {/Symbol a}=2.0',\
     m.'mmdt_ang_0.5' .r u 4:8 w l ls 7 t 'mMDT, {/Symbol a}=0.5',\
     m.'mmdt_ang_1.0' .r u 4:8 w l ls 8 t 'mMDT, {/Symbol a}=1.0',\
     m.'mmdt_ang_2.0' .r u 4:8 w l ls 9 t 'mMDT, {/Symbol a}=2.0'
     
set title '{} {/:Bold ECFs} {}'
plot m.'plain_ecf_0.5'.r u 4:8 w l ls 1 t 'plain, {/Symbol a}=0.5',\
     m.'plain_ecf_1.0'.r u 4:8 w l ls 2 t 'plain, {/Symbol a}=1.0',\
     m.'plain_ecf_2.0'.r u 4:8 w l ls 3 t 'plain, {/Symbol a}=2.0',\
     m.'sd_ecf_0.5'   .r u 4:8 w l ls 4 t 'SD, {/Symbol a}=0.5',\
     m.'sd_ecf_1.0'   .r u 4:8 w l ls 5 t 'SD, {/Symbol a}=1.0',\
     m.'sd_ecf_2.0'   .r u 4:8 w l ls 6 t 'SD, {/Symbol a}=2.0',\
     m.'mmdt_ecf_0.5' .r u 4:8 w l ls 7 t 'mMDT, {/Symbol a}=0.5',\
     m.'mmdt_ecf_1.0' .r u 4:8 w l ls 8 t 'mMDT, {/Symbol a}=1.0',\
     m.'mmdt_ecf_2.0' .r u 4:8 w l ls 9 t 'mMDT, {/Symbol a}=2.0'
     
set title '{} {/:Bold ISDs} {}'
plot m.'plain_isd_1.0'.r u 4:8 w l ls 1 t 'plain, k_t>=1.0',\
     m.'plain_isd_2.0'.r u 4:8 w l ls 2 t 'plain, k_t>=2.0',\
     m.'plain_isd_5.0'.r u 4:8 w l ls 3 t 'plain, k_t>=5.0',\
     m.'sd_isd_1.0'   .r u 4:8 w l ls 4 t 'SD, k_t>1.0',\
     m.'sd_isd_2.0'   .r u 4:8 w l ls 5 t 'SD, k_t>2.0',\
     m.'sd_isd_5.0'   .r u 4:8 w l ls 6 t 'SD, k_t>5.0',\
     m.'mmdt_isd_1.0' .r u 4:8 w l ls 7 t 'mMDT, k_t>1.0',\
     m.'mmdt_isd_2.0' .r u 4:8 w l ls 8 t 'mMDT, k_t>2.0',\
     m.'mmdt_isd_5.0' .r u 4:8 w l ls 9 t 'mMDT, k_t>5.0'

set out
