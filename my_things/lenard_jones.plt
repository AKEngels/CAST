set xlabel "R_{AB}" 
set ylabel "Energie"
set xtics ('0' 0, "R_{AB,0}" 1)
set ytics ('0' 0, "-{/Symbol e}_{AB}" -1)
unset key
set border
set yrange[-1.5:3]
set xzeroaxis
plot [0:3] (1/x)**12 - 2*(1/x)**6 lt rgb "blue"