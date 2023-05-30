set terminal pdf
set output 'phonon-bandplot.pdf'
maxy=65.0
set xrange[0.0:4.7462]
set yrange[0.0:maxy]
set xtics (" {/Symbol G} "  0.00000000," X "  1.0000," W "  1.5000, " L "  2.2071, " K "  2.8195," {/Symbol G} "  3.8801," L"   4.7462)
set arrow from  1.0000,  0.0 to  1.0000,  maxy nohead
set arrow from  1.5000,  0.0 to  1.5000,  maxy nohead
set arrow from  2.2071,  0.0 to  2.2071,  maxy nohead
set arrow from  2.8195,  0.0 to  2.8195,  maxy nohead
set arrow from  3.8801,  0.0 to  3.8801,  maxy nohead
set key outside below
plot 'Si.freq.gp' using 1:($2*0.123984) w l lt -1 lc rgb "red" lw 3 noti , \
'Si.freq.gp' using  1:($3*0.123984) w l lt -1 lc rgb "red" lw 3 noti , \
'Si.freq.gp' using  1:($4*0.123984) w l lt -1 lc rgb "red" lw 3 noti , \
'Si.freq.gp' using  1:($5*0.123984) w l lt -1 lc rgb "red" lw 3 noti , \
'Si.freq.gp' using  1:($6*0.123984) w l lt -1 lc rgb "red" lw 3 noti , \
'Si.freq.gp' using  1:($7*0.123984) w l lt -1 lc rgb "red" lw 3 noti

