#!/bin/bash

if [ $# -eq 1 ]; then
    quantity=$1
else
    echo "Requires 1 arg (given $#):" >&2
    echo "* the quantity of interest (e.g. Sxx, v etc)" >&2
    exit -1
fi


if [ `find . -maxdepth 1 -name 'st-*' | wc -l` -ne 0 ]; then

    # -r to get most recent images first
    filelist=`ls -r st-* | tr '\n' ' ' | tr -d 'st-'`
    num_files=`ls st-* | wc -w`

    if [ $num_files -eq 0 ]; then
        echo "No files found!" >&2
        exit -1;
    fi
else
    echo "No files found!" >&2
    exit -1;
fi


# get system size from output
Lx=`awk '/#L#/ {print $2}' output-data | tr -d ','`
Ly=`awk '/#L#/ {print $3}' output-data | tr -d ','`
ratio=`echo "scale=4; $Ly / $Lx" | bc`

if [ "$quantity" == "Sxx" ]; then
    quantity_col=3
elif [ "$quantity" == "Sxy" ]; then
    quantity_col=4
elif [ "$quantity" == "Syy" ]; then
    quantity_col=5
elif [ "$quantity" == "trS" ]; then
    quantity_col=6
elif [ "$quantity" == "v" ]; then
    quantity_col=7; quantity_x_col=8; quantity_y_col=9
    vector_str="vectors head size 0.05,20,60 filled"
elif [ "$quantity" == "vort" ]; then
    quantity_col=10

elif [ "$quantity" == "f" ]; then
    quantity_col=11; quantity_x_col=12; quantity_y_col=13
    vector_str="vectors head size 0.05,20,60 filled"


# calculated from existing
elif [ "$quantity" == "Q" ]; then
    quantity_col="(f(\$9, \$10, \$11))"
elif [ "$quantity" == "D" ]; then
    quantity_col="(S2(\$9, \$10))"
elif [ "$quantity" == "QD" ]; then
    quantity_col="(S2(\$9, \$10) * f(\$9, \$10, \$11))"

else
    echo "quantity not recognised!"
    exit -1
fi


if [ ! -z "${quantity_col##*[!0-9]*}" ]; then
    # find the extrema (so as to stop the colorscale jumping)
    cb_string=`awk '!/^(#|$)/ {print $'$quantity_col'}' st-* | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1};} END {printf "set cbrange [%.3g to %.3g]\n", min, max}'`
    echo "col = $quantity_col"
    echo $cb_string
fi

scale=0.25

mkdir -p images
rm -f images/${quantity}*

if [ -n "$quantity_x_col" ] && [ -n $quantity_y_col ]; then
    vector_overlay="data_file every 8:8 using (\$1 - 0.5*\$$quantity_x_col):(\$2 - 0.5*\$$quantity_y_col):(0):($scale*\$$quantity_x_col):($scale*\$$quantity_y_col):(0) with $vector_str linecolor rgb \"red\""

    if [ "$quantity" == "q" ]; then
        # in a seperate file (starting dt-00...) are a list of defects:
        #   +1/2 strength have 3rd col 1 (green)
        #   -1/2 strength have 3rd col 0 (red)
        vector_overlay="$vector_overlay, defect_file using 1:(\$3==1?\$2:1/0):(0) with points pt 7 ps 1 linecolor rgb \"green\""
        vector_overlay="$vector_overlay, defect_file using 1:(\$3==0?\$2:1/0):(0) with points pt 7 ps 1 linecolor rgb \"red\""
    fi
fi



gnuplotscript=$'
set palette rgbformulae 30,31,32
set terminal png size 635,695
set size ratio '$ratio'
set nokey
set nogrid
set xlabel "x"
set ylabel "y"
set xrange [0:'$Lx']
set yrange [0:'$Ly']
set parametric
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

'$cb_string'

set pointsize 5

# useful functions
S2(e, s) = 0.25 * (4.0 * e * e + s * s)
O2(as) = 0.25 * as * as
f(e, s, as) = (S2(e, s) - O2(as)) / (S2(e, s) + O2(as))

do for [i=1:words(filelist)] {

    data_file = "st-".word(filelist, i)
    defect_file = "dt-".word(filelist, i)
    set title sprintf("'$quantity': t = %09.3f ", word(filelist, i) + 0)
    set output "images/'$quantity'_".word(filelist, i).".png"
    set pm3d map
    set contour
    splot data_file using 1:2:'$quantity_col', '$vector_overlay'
    print sprintf("%09.3f ", word(filelist, i) + 0)
}
'



echo "filelist = \"$filelist\" $gnuplotscript" > plotting_$quantity
gnuplot plotting_$quantity

if [ $num_files -ge 10 ]; then
    echo "Stitching $num_files files..."
    fps=15; mov_w=835; mov_h=395
    mencoder "mf://./images/${quantity}*.png" -mf w=$mov_w:h=$mov_h:fps=$fps:type=png -ovc copy -oac copy -o ./images/output_${quantity}.avi
fi

echo "Finished!"
