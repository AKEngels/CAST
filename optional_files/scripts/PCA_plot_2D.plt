# USER INPUT
number_of_ticks = 4
inputfile = "pca_histogrammed.dat"
outputfile = "pca_hist.png"   

# OVERALL STUFF ABOUT PICTURE
set terminal png
set output outputfile
set pm3d map
set size square 1,1
set palette model CMY rgbformulae 7,5,15

# GETTING STATISTICS
stats inputfile using 1 name 'x' nooutput
stats inputfile using 2 name 'y' nooutput
x_diff = x_max - x_min
y_diff = y_max - y_min

# SET TICS STUFF
set tics out
set tics format "%.1f%%"
set xtics format "%.2t*10^{%S}"
set ytics format "%.2t*10^{%S}"
set ztics format "%.1f%%"
set format z "%.0f%%"
set xrange [x_min : x_max]
set yrange [y_min : y_max]
set xtics x_diff/number_of_ticks
set ytics y_diff/number_of_ticks

# DO PLOTTING
splot inputfile using 1:2:($3*100) title ""
