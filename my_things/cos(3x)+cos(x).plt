set xrange[0:360]
set xlabel "{/Symbol w}[°]" 
set xtics 0, 60, 360
unset ytics
set key right bottom
plot [0:360] 0.5*cos(3*x*pi/180) + 0.5*cos(x*pi/180) lt rgb "blue" title "n = 1, n = 3", 0.5*cos(3*x*pi/180) lt rgb "red" dt 5 title "n = 3", 0.5*cos(x*pi/180) lt rgb "green" dt 5 title "n = 1"