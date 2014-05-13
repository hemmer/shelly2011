filelist = "flow.dat"
set palette rgbformulae 30,31,32
set terminal png size 735,395
set terminal png size 835,395
set size ratio 1
set nokey
set nogrid
set xlabel "x"
set ylabel "y"
set xrange [0:1]
set yrange [0:1]
set parametric
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

# set cbrange [9.5 to 48]

set pointsize 5

# useful functions
S2(e, s) = 0.25 * (4.0 * e * e + s * s)
O2(as) = 0.25 * as * as
f(e, s, as) = (S2(e, s) - O2(as)) / (S2(e, s) + O2(as))

do for [i=1:words(filelist)] {

    data_file = "flow.dat"
    set title "vorticity"
    set output "vorticity.png"
    set pm3d map
    splot data_file using 1:2:3

    set title "vx"
    set output "vx.png"
    set pm3d map
    splot data_file using 1:2:5
}

