
set title "Temperature distribution" font "default, 16"

set size square

set xtics font "default, 12"
set ytics font "default, 12"
set cbtics font "default,12"

set xlabel "x, [m]" font "default, 14"
set ylabel "y, [m]" font "default, 14"
set cblabel "T, [C]" font "default, 14" offset 2,0

p "output.dat" u 1:2:3 w image noti