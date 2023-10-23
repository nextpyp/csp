#!/bin/bash
#
# This script accept eighteen input parameters corresponding with the
# fields of the following structure and format them with the
# corresponding format.
#
# This is the structure that should be readed form the .par file, one
# line per image:
# C*********************************************************************
#       type img_prop
# C*********************************************************************
#       integer   img_index       ! Image number: Image index in the data
#                                 ! stack.
#       real      psi             ! First Euler angle: psi.
#       real      theta           ! Second Euler anlge: theta.
#       real      phi             ! Third Euler angle: phi.
#       real      shift_x         ! Horizontal translation in the image.
#       real      shift_y         ! Vertical trasnlation in the image.
#       real      mag             ! Micrograph magnification.
#       integer   film            ! Film (micrograph) number: identifies
#                                 ! which images are from same film.
#       real      df1             ! Defocus 1 (in \AA)
#       real      df2             ! Defocus 2 (in \AA)
#       real      angast          ! Astigmatism
#       real      presa           ! Phase residual: see publication #2 in
#                                 ! the manual for definition.
#       real      dpres           ! Difference in phase residual with
#                                 ! previous iteration
# C     New parameters
#       integer   ptl_index       ! Particle number: spike index (this
#                                 ! particle is a view of this spike).
#       real      tilt_angle      ! Angle of the micrograph where belongs
#                                 ! this image.
#       real      dose            ! Dose used in scanning the micrograph.
#       integer   scan_order      ! Number of the micrograph in the
#                                 ! scanned sequence.
#       real      confidence      ! Confidence in this particle.
#       real      ptl_CCX         ! Cross correlation for this spike with
#                                 ! the model (all the set of view with
#                                 ! respect to the model).
#       end type img_prop
# C*********************************************************************
#
# This is the fortran formatting line expected by frealign in each
# line of the .par file:
# 7011  FORMAT(I7,5F8.2,F8.0,I6,2F9.1,F8.2,F7.2,F8.2,I9,2F9.2,I9,2F9.2)
#
# $Id: format_dotpar_line.sh 173 2010-06-16 14:40:01Z fefo $

#
# Usage: format_dotpar_line.sh img_index psi theta phi shift_x shift_y
#                              mag film df1 df2 angast presa dpres
#                              ptl_index tilt_angle dose scan_order
#                              confidence ptl_CCX tilt_axis_angle
#                              normal[3] matrix[16]
#

tarray=( $@ )
#echo ${tarray[*]} $#

if [ $# -ne 39 ]; then

    echo "Usage: format_dotpar_line.sh img_index psi theta phi shift_x shift_y mag film df1 df2 angast \\"
    echo "                             presa dpres ptl_index tilt_angle dose scan_order confidence \\"
    echo "                             ptl_CCX tilt_axis_angle normal[0] normal[1] normal[2] \\"
    echo "                             matrix[00] matrix[01] matrix[02] matrix[03] \\"
    echo "                             matrix[10] matrix[11] matrix[12] matrix[13] \\"
    echo "                             matrix[20] matrix[21] matrix[22] matrix[23] \\"
    echo "                             matrix[30] matrix[31] matrix[32] matrix[33]"

    echo "e$#"
    exit 1
fi 

img_index=${tarray[0]}
psi=${tarray[1]}
theta=${tarray[2]}
phi=${tarray[3]}
shift_x=${tarray[4]}
shift_y=${tarray[5]}
mag=${tarray[6]}
film=${tarray[7]}
df1=${tarray[8]}
df2=${tarray[9]}
angast=${tarray[10]}
presa=${tarray[11]}
dpres=${tarray[12]}
ptl_index=${tarray[13]}
tilt_angle=${tarray[14]}
dose=${tarray[15]}
scan_order=${tarray[16]}
confidence=${tarray[17]}
ptl_CCX=${tarray[18]}
#
tilt_axis_angle=${tarray[19]}
#
normal_0=${tarray[20]}
normal_1=${tarray[21]}
normal_2=${tarray[22]}
#
matrix_00=${tarray[23]}
matrix_01=${tarray[24]}
matrix_02=${tarray[25]}
matrix_03=${tarray[26]}
matrix_10=${tarray[27]}
matrix_11=${tarray[28]}
matrix_12=${tarray[29]}
matrix_13=${tarray[30]}
matrix_20=${tarray[31]}
matrix_21=${tarray[32]}
matrix_22=${tarray[33]}
matrix_23=${tarray[34]}
matrix_30=${tarray[35]}
matrix_31=${tarray[36]}
matrix_32=${tarray[37]}
matrix_33=${tarray[38]}

# Print one formatted row 
echo -n "`printf %7d $img_index`"
echo -n "`printf %8.2f $psi`"
echo -n "`printf %8.2f $theta`"
echo -n "`printf %8.2f $phi`"
echo -n "`printf %8.2f $shift_x`"
echo -n "`printf %8.2f $shift_y`"
echo -n "`printf %8.0f $mag`"
echo -n "`printf %6d $film`"
echo -n "`printf %9.1f $df1`"
echo -n "`printf %9.1f $df2`"
echo -n "`printf %8.2f $angast`"
echo -n "`printf %7.2f $presa`"
echo -n "`printf %8.2f $dpres`"
echo -n "`printf %9d $ptl_index`"
echo -n "`printf %9.2f $tilt_angle`"
echo -n "`printf %9.2f $dose`"
echo -n "`printf %9d $scan_order`"
echo -n "`printf %9.2f $confidence`"
#echo "`printf %9.2f $ptl_CCX`"
echo -n "`printf %9.2f $ptl_CCX`"
#
# The following fields are not yet considered in our modified version
# of Frealign. Should run without problem (reading the new .par file)
# but will not write the new fields in the output .par
#
echo -n "`printf %9.2f $tilt_axis_angle`"
#
echo -n "`printf %9.2f $normal_0`"
echo -n "`printf %9.2f $normal_1`"
echo -n "`printf %9.2f $normal_2`"
#
echo -n "`printf %9.2f $matrix_00`"
echo -n "`printf %9.2f $matrix_01`"
echo -n "`printf %9.2f $matrix_02`"
echo -n "`printf %9.2f $matrix_03`"
echo -n "`printf %9.2f $matrix_10`"
echo -n "`printf %9.2f $matrix_11`"
echo -n "`printf %9.2f $matrix_12`"
echo -n "`printf %9.2f $matrix_13`"
echo -n "`printf %9.2f $matrix_20`"
echo -n "`printf %9.2f $matrix_21`"
echo -n "`printf %9.2f $matrix_22`"
echo -n "`printf %9.2f $matrix_23`"
echo -n "`printf %9.2f $matrix_30`"
echo -n "`printf %9.2f $matrix_31`"
echo -n "`printf %9.2f $matrix_32`"
echo -n "`printf %9.2f $matrix_33`"
#
echo ""
