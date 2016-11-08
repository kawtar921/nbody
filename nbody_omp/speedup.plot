set terminal png
set style line 1 linetype 1 linewidth 2 pointsize 0
set style line 2 linetype 2 linewidth 2
set style line 3 linetype 3 linewidth 2
set style line 4 linetype 4 linewidth 2
set style line 5 linetype 5 linewidth 2
set style line 6 linetype 6 linewidth 2
set style line 7 linetype 7 linewidth 2
set style line 8 linetype 8 linewidth 2
set style line 9 linetype 9 linewidth 2
set ylabel "speedup"
set xlabel "number of cores"
set xtics 1
set title "Mandelbrot Speedup for OpenMP execution"
set output "SpeedUpMandelbrot-OpenMP.png"

plot x ls 1 t 'Speedup max' with lines, 'speedup_omp.data' using 1:($2/$3) t 'Mandelbrot' ls 3 with linespoints;
