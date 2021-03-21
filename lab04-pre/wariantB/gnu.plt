set term png
set view map
set pm3d interpolate 4,4
set xlabel "x"
set ylabel "y"

do for [t=3:12] {
    set out 'vec_'.(t-2).'.png'
    splot 'dane.dat' u 1:2:t w pm3d
}