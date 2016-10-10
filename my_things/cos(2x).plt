set xrange[0:360]
set xlabel "{/Symbol w}[°]" 
set xtics 0, 60, 360
unset key
unset ytics
plot [0:360] -cos(2*x*pi/180) lt rgb "blue"