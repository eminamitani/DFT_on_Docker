set term png 
set output 'dos.png'
# highest occupied, lowest unoccupied level (ev):    10.9263   13.1150
set xrange [9:18]
plot "GaN.dos" using 1:2 w l