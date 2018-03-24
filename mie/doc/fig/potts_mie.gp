#!/usr/bin/env gnuplot



# Potts model



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced size 10, 5 font "Helvetica, 36"
set output "potts_mie.eps"

reset

set multiplot

wl = 0.52
wr = 1 - wl

set size wl, 1
set origin 0, 0
#set lmargin 8

set rmargin 2

#set logscale x
set xtics 2 offset 0, 0.2
set mxtics 2
#set format x "10^{/*0.8 %T}"
#set mxtics 10
#set xrange [1e4:1e8]
set xlabel "Simulation time, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 4})" offset 0, 0.5

#set logscale y
#set format y "10^{/*0.8 %T}"
set ytics 2
set mytics 2
#set ytics add ("{/Times {/Symbol a}_0/2}" alpha0*0.5)
#set mytics 10
set yrange [-4:20]
set ylabel 'Estimated entropy, ~{/Times-Italic S}{0.5\^}' offset 1.2, 0

set title "(a) Second-order MIE"
set key left bottom Left reverse spacing 1.5 at 0.4, -3

#set arrow from 6.0e3, alpha0*0.5 to 2e4, alpha0*0.5 lt 1 lw 1 filled size screen 0.02,10,35
#set label "{/Times {/Times-Italic a}_0&{/*0.5 i}/2}" at 2.0e3, alpha0*0.5

fn = "../../data/potts/potts_q6_n10_T0.5.log"

plot [:4.2][:] \
    fn u ($1/1e4):9:10  every 1 w lp lt 1 pt 5 t "Uncorrected", \
    fn u ($1/1e4):11:12 every 1 w lp lt 2 pt 6 t "Linearly corrected", \
    fn u ($1/1e4):13:14 every 1 w lp lt 4 pt 8 t "Exponentially corrected", \
    fn u ($1/1e4):15 w l lt 1 t "Reference"


set size wr, 1
set origin wl, 0

set title "(b) Third-order MIE"
set key at 2, -3

set lmargin 0
unset ylabel
unset format y ""

plot [:10][:] \
    fn u ($1/1e4):16:17 every 1 w lp lt 1 pt 5 t "Uncorrected", \
    fn u ($1/1e4):18:19 every 1 w lp lt 2 pt 6 t "Linearly corrected", \
    fn u ($1/1e4):20:21 every 1 w lp lt 4 pt 8 t "Exponentially corrected", \
    fn u ($1/1e4):22 w l lt 1 t "Reference"


unset multiplot
unset output
set terminal pop
reset
