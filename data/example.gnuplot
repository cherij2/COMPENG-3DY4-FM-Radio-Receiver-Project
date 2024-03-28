reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 3,1 scale 1.0,1.0  # set two plots for this figure

# time domain
set ylabel 'Sample value'               # set y-axis label
set xlabel 'dem_mixer output'                   # set x-axis label
set yrange [-1:1]                       # set y plot range
set xrange [0:30]                      # sset x plot range


plot '../data/bits.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

# freq domain (Fourier)
set ylabel 'Spectrum (Mag)'              # set y-axis label
set xlabel 'rrc outp'               # set x-axis label
set yrange [-1.2:1.2]                    # set y plot range
set xrange [0:1000] 
set y2range [-0.2:1.2]                      # set x plot range
set x2range[0:33]
plot '../data/outp_rrc.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle \
     '../data/bits.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#32a852' notitle axes x2y2

# freq domain (PSD)
set ylabel 'Spectrum (dB/Hz)'            # set y-axis label
set xlabel 'PLl Input'             # set x-axis label
set autoscale
set yrange [-1:1]                       # set y plot range
set y2range [-1:1]
set xrange [0:119]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/rds_NCO_outp.dat' using 1:2 with lines lt 1 lw 3 lc rgb '#880000' notitle, \
     '../data/CR_rds_filtered.dat' using 1:2 with lines lt 1 lw 3 lc rgb '#880000' notitle axes x1y2

unset multiplot
