#!/bin/bash
#
# Version: $Id: plot_all.sh 566 2011-02-27 22:25:19Z fefo $
#

expBaseName=$1
outdir=$2

# Set title and output type 
title=${expBaseName}
font=/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans.ttf
fontsize=9
width=850
height=850

outputFile=results/${expBaseName}_all_results.png
gnustr="set terminal png font '${font},${fontsize}' size ${width}, ${height}"
gnustr="${gnustr}; set output '${outputFile}'"
# gnustr="${gnustr}; set output '| display png:-'"
gnustr="${gnustr}; set multiplot layout 2, 2 title '${title} results'"

# Plot
title="Mean Phase Residual"
xlabel="noise"
ylabel="ini"
str="set title '${title}'"
str="${str}; set zrange [*:*]"
str="${str}; set xlabel '${xlabel}'"
# str="${str}; set xtics (\"0\" 0, \"15\" 1, \"30\" 2, \"45\" 3, \"60\" 4) "
str="${str}; set ylabel '${ylabel}'"
# str="${str}; set ytics (\"0\" 0, \"1\" 1, \"3\" 2, \"5\" 3) "
# str="${str}; set hidden3d"
str="${str}; set view 70,330"
str="${str}; splot"
str="${str} '${outdir}/${expBaseName}_SCILAB_CSP_END_PR' matrix with lines ti \"CSP END\", "
str="${str} '${outdir}/${expBaseName}_SCILAB_USP_END_PR' matrix with lines ti \"USP END\" "
gnustr="${gnustr}; ${str}"
 
# Plot
title="Frobenius of difference of Transformation Matrices"
xlabel="noise"
ylabel="ini"
str="set title '${title}'"
str="${str}; set zrange [*:*]"
str="${str}; set xlabel '${xlabel}'"
# str="${str}; set xtics (\"0\" 0, \"15\" 1, \"30\" 2, \"45\" 3, \"60\" 4) "
str="${str}; set ylabel '${ylabel}'"
# str="${str}; set ytics (\"0\" 0, \"1\" 1, \"3\" 2, \"5\" 3) "
# str="${str}; set hidden3d"
str="${str}; set view 70,330"
str="${str}; splot"
str="${str} '${outdir}/${expBaseName}_SCILAB_CSP_END_Frob' matrix with lines ti \"CSP END\", "
str="${str} '${outdir}/${expBaseName}_SCILAB_USP_END_Frob' matrix with lines ti \"USP END\" "
gnustr="${gnustr}; ${str}"

# Plot
title="FSC"
xlabel="noise"
ylabel="ini"
str="set title '${title}'"
str="${str}; set zrange [*:*]"
str="${str}; set xlabel '${xlabel}'"
# str="${str}; set xtics (\"0\" 0, \"15\" 1, \"30\" 2, \"45\" 3, \"60\" 4) "
str="${str}; set ylabel '${ylabel}'"
# str="${str}; set ytics (\"0\" 0, \"1\" 1, \"3\" 2, \"5\" 3) "
# str="${str}; set hidden3d"
str="${str}; set view 70,330"
str="${str}; splot"
str="${str} '${outdir}/${expBaseName}_SCILAB_CSP_FSC' matrix with lines ti \"CSP END\", "
str="${str} '${outdir}/${expBaseName}_SCILAB_USP_FSC' matrix with lines ti \"USP END\" "
gnustr="${gnustr}; ${str}"

# Plot
title="FSC (vs. GT)"
xlabel="noise"
ylabel="ini"
str="set title '${title}'"
str="${str}; set zrange [*:*]"
str="${str}; set xlabel '${xlabel}'"
# str="${str}; set xtics (\"0\" 0, \"15\" 1, \"30\" 2, \"45\" 3, \"60\" 4) "
str="${str}; set ylabel '${ylabel}'"
# str="${str}; set ytics (\"0\" 0, \"1\" 1, \"3\" 2, \"5\" 3) "
# str="${str}; set hidden3d"
str="${str}; set view 70,330"
str="${str}; splot"
str="${str} '${outdir}/${expBaseName}_SCILAB_CSP_GTFSC' matrix with lines ti \"CSP END\", "
str="${str} '${outdir}/${expBaseName}_SCILAB_USP_GTFSC' matrix with lines ti \"USP END\" "
gnustr="${gnustr}; ${str}"

# Run gnuplot
command=`echo "${gnustr}" | gnuplot`
# echo ${command}
$command 

display ${outputFile} &
    
