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

if [ $# -ne 13 ]; then

    echo "Usage: format_dotpar_line.sh img_index psi theta phi shift_x shift_y mag film df1 df2 angast \\"
    echo "                             presa dpres"

    echo "e$#"
    exit 1
fi 

# img_index=${tarray[0]}
# psi=${tarray[1]}
# theta=${tarray[2]}
# phi=${tarray[3]}
# shift_x=${tarray[4]}
# shift_y=${tarray[5]}
# mag=${tarray[6]}
# film=${tarray[7]}
# df1=${tarray[8]}
# df2=${tarray[9]}
# angast=${tarray[10]}
# presa=${tarray[11]}
# dpres=${tarray[12]}

# # Print one formatted row 
# echo -n "`printf %7d $img_index`"
# echo -n "`printf %8.2f $psi`"
# echo -n "`printf %8.2f $theta`"
# echo -n "`printf %8.2f $phi`"
# echo -n "`printf %8.2f $shift_x`"
# echo -n "`printf %8.2f $shift_y`"
# echo -n "`printf %8.0f $mag`"
# echo -n "`printf %6d $film`"
# echo -n "`printf %9.1f $df1`"
# echo -n "`printf %9.1f $df2`"
# echo -n "`printf %8.2f $angast`"
# echo -n "`printf %7.2f $presa`"
# echo -n "`printf %8.2f $dpres`"
# #
# echo ""

echo -n "`printf %7d ${tarray[0]}``printf %8.2f ${tarray[1]}``printf %8.2f ${tarray[2]}``printf %8.2f ${tarray[3]}``printf %8.2f ${tarray[4]}``printf %8.2f ${tarray[5]}``printf %8.0f ${tarray[6]}``printf %6d ${tarray[7]}``printf %9.1f ${tarray[8]}``printf %9.1f ${tarray[9]}``printf %8.2f ${tarray[10]}``printf %7.2f ${tarray[11]}``printf %8.2f ${tarray[12]}`"
echo ""
