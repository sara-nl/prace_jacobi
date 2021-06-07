
set title "Temperature distribution" font "default, 16"

set size square

set xtics font "default, 12"
set ytics font "default, 12"
set cbtics font "default,12"

set xlabel "x, [m]" font "default, 14"
set ylabel "y, [m]" font "default, 14"
set cblabel "T, [C]" font "default, 14" offset 2,0

# Convert binary to ASCII
`hexdump -v -e '3/8 "%6f "' -e '"\n"' output.dat > output_ascii.dat`

# Sort by "y" coordinate first, then by "x" coordinate
p "< sort -nk2,2 -nk1,1 output_ascii.dat" u 1:2:3 w image noti