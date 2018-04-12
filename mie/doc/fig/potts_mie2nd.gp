#!/usr/bin/env gnuplot



# Potts model



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced color size 9, 7 font "Helvetica, 36"
set output "potts_mie2nd.eps"

reset

set style line 1 lt 1 lc rgb "#002080" lw 1.0 pt  7 ps 1.7
set style line 2 lt 2 lc rgb "#a00000" lw 1.5 pt  4 ps 1.4
set style line 3 lt 3 lc rgb "#a00040" lw 1.5 pt 12 ps 2.0
set style line 4 lt 4 lc rgb "#00a020" lw 1.5 pt  8 ps 2.0
set style line 5 lt 5 lc rgb "#a0a000" lw 1.5 pt 10 ps 2.0
set style line 7 lt 1 lc rgb "#404040"

set multiplot

set size 1, 1
set origin 0, 0

set xtics 1 offset 0, 0.2
set mxtics 2
set xlabel "Simulation time, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 4})" offset 0, 0.5

set ytics 5
set mytics 5
set ylabel 'Estimated entropy, ~{/Times-Italic S}{0.5\^}' offset 1.0, 0

#set title "(a) Second-order MIE"
set key left bottom Left reverse spacing 1.5 at 2.0, 13.5 maxrows 2

#set arrow from 6.0e3, alpha0*0.5 to 2e4, alpha0*0.5 lt 1 lw 1 filled size screen 0.02,10,35
#set label "{/Times {/Times-Italic a}_0&{/*0.5 i}/2}" at 2.0e3, alpha0*0.5

fn = "../../data/potts/potts_q6_n10_T0.5.log"
ref = `head -n 1 @fn | cut -f15`
#print ref

plot [0:10][-3:18.5] \
    fn u ($1/1e4): 9:10 every 1 w lp ls 1 ps 2 t "Uncorrected", \
    fn u ($1/1e4):11:12 every 1 w lp ls 2 ps 2 t "Linear", \
    fn u ($1/1e4):13:14 every 1 w lp ls 4 ps 2 t "Exponential", \
    ref w l lt 1 t "Reference"


set size 0.6, 0.52
set origin 0.35, 0.15

unset key

set logscale x
set logscale y

set format x "10^{/*0.7 %T}"
set xtics auto offset 0, 0.2
set mxtics default
set xlabel "Simulation time, {/Times-Italic  t}" offset 0, 0.3

set format y "10^{/*0.7 %T}"
set ytics auto
set mytics default
set ylabel '{/*1.8 |} ~{/Times-Italic S}{0.5\^} - {/Times-Italic S}^{ {/Times ref} }{/*1.8 |}' offset -1, 0

set style arrow 3 head filled size screen 0.01,15,45 ls 1
set arrow 1 from 6e3, 10 to 3e3, 7 arrowstyle 3
set label 1 at 6.5e3, 11 "2.14{/Symbol \264}10^{/*0.7 4}{/Times /{/Times-Italic t}}"

set arrow 2 from 1e4, 0.02 to 2.56e4, 0.02 arrowstyle 3
set label 2 at 2.1e3, 0.023 "1.33{/Symbol \264}10^{/*0.7 7}{/Times / {/Times-Italic t}}^{/*0.7 2}"

plot [1e3:1e5][0.01:30] \
    fn u ($1):(abs($9-$15))  every 1 w lp ls 1 ps 0.7 t "Uncorrected", \
    fn u ($1):(abs($11-$15)) every 1 w lp ls 2 ps 0.7 t "Linearly", \
    fn u ($1):(abs($13-$15)) every 1 w lp ls 4 ps 0.7 t "Exponential", \
    2.14e4/x w l lt 7, \
    1.33e7/(x*x) w l lt 6


unset multiplot
unset output
set terminal pop
reset
