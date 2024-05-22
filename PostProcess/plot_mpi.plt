# Copyright (c) 2024 Maksim Masterov, SURF
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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