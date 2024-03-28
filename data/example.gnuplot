reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 3,1 scale 1.0,1.0  # set two plots for this figure

# time domain
set ylabel 'Sample value'               # set y-axis label
set xlabel 'dem_mixer output'                   # set x-axis label
set yrange [-2:2]                       # set y plot range
set xrange [0:20]                      # set x plot range
plot '../data/dem_mixer.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

# freq domain (Fourier)
set ylabel 'Spectrum (Mag)'              # set y-axis label
set xlabel 'rrc outp'               # set x-axis label
set yrange [-1.2:1.2]                    # set y plot range
set xrange [0:1000]                       # set x plot range
plot '../data/outp_rrc.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

# freq domain (PSD)
set ylabel 'Spectrum (dB/Hz)'            # set y-axis label
set xlabel 'PLl Input'             # set x-axis label
set autoscale
#set yrange [:]                       # set y plot range
set xrange [0:119]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/rds_NCO_outp.dat' using 1:2 with lines lt 1 lw 3 lc rgb '#880000' notitle

unset multiplot
