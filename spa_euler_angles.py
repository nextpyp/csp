def spa_euler_angles( tilt_angle, tilt_axis_angle, tilt_correction, normal, m, cutOffset):

    # 3DAVG refinement
    # This transformation matrix is different from the matrix[0-15] in the 3DAVG refinement file.
    # The ordering of the rotations has to be changed to match the new geometry.
    # 	- The ordering in 3DAVG is: RotZ1 * RotX * RotZ2.
    #	- The ordering here is: RotZ2 * RotX * RotZ1
    # The new transformation matrix has the following expression:

    # 3DAVG refinement transformation
    refinement = numpy.matrix( [ [m[0], -m[4], m[8], m[3]],
                               [-m[1], m[5],-m[9], m[7]],
                               [m[2],-m[6], m[10], m[11]],
                               [m[12], m[13], m[14], m[15]] ])

    # 3DAVG refinement rotation
    refinement_rotation = numpy.matrix( [ [m[0], -m[4], m[8], 0],
                                         [-m[1], m[5],-m[9], 0],
                                         [m[2],-m[6], m[10], 0],
                                         [m[12], m[13], m[14], m[15]] ])

    ## 3DAVG refinement transformation
    #refinement = numpy.matrix( [ [m[0], m[4], m[8], m[3]],
    #                           [m[1], m[5], m[9], m[7]],
    #                           [m[2], m[6], m[10], m[11]],
    #                           [m[12], m[13], m[14], m[15]] ])
    #
    ## 3DAVG refinement rotation
    #refinement_rotation = numpy.matrix( [ [m[0], m[4], m[8], 0],
    #                                     [m[1], m[5],m[9], 0],
    #                                     [m[2], m[6], m[10], 0],
    #                                     [m[12], m[13], m[14], m[15]] ])

    # RotateZ( tilt_axis_angle );
    r2D = vtk.rotation_matrix( numpy.radians(tilt_axis_angle), [0,0,1] )

    # correction in the direction perpendicular to the tilt axis
    # correction = [ -.5, 0, 0, 1 ]
    # tcorrection = numpy.dot( r2D, correction )

    # print 'tcorrection', tcorrection

    # apply .box coordinate discretization error
    t = vtk.translation_matrix( [tilt_correction[0], tilt_correction[1], 0] )

    # t = vtk.translation_matrix( [tcorrection[0], tcorrection[1], 0] )

    # tilt axis angle rotation: t->RotateZ( tilt_axis_angle )
    t = numpy.dot( r2D, t )

    # Difference vector between 2D rotation origins (IMOD vs. FREALIGN)
    diff2D = [ 1, 1, 0, 1 ]

    # Compute: t = Rot(C1) wrt C2 - C1 = Rot * ( C1 - C2 ) - ( C1 - C2 )
    tdiff2D = numpy.dot( r2D, diff2D ) - diff2D

    # print 'tdiff2D', tdiff2D

    # apply rotation center and refinement translations: t->Translate( -tdiff2D[0], -tdiff2D[1], -tdiff2D[2] );
    t = numpy.dot( vtk.translation_matrix( -tdiff2D ), t )

    # correct for center of tilt axis: t->Translate( -1, 0, 0 );
    t = numpy.dot( vtk.translation_matrix( [-1,0,0] ), t )

    # apply tilt angle rotation: t->RotateY( -tilt_angle );
    t = numpy.dot( vtk.rotation_matrix( numpy.radians(-tilt_angle), [0,1,0] ), t )

    # Convert to IMOD's tilt axis location: t->Translate( .5, 0, .5 );
    t = numpy.dot( vtk.translation_matrix( [.5,0,.5] ), t )

    # The remaining transformations are the spike normal and the 3DAVG refinement transformation.
    # Spike normals are a pure rotation R1. 3DAVG refinement is a full rotation and translation matrix F = R2 * T2
    # If the two transformations are composed, the net rotation is then: R = R1 * R2 and the translation component is T2.
    # The origin of the net rotation R is different from the image rotation, so to account for this we express
    # as R * T3, where T3 corrects for the difference in the origin of the rotation.

    # apply spike euler angles (pure rotation R1)
    t = numpy.dot( vtk.rotation_matrix( numpy.radians(-normal[2]), [0,0,1] ), t )
    t = numpy.dot( vtk.rotation_matrix( numpy.radians(-normal[0]), [1,0,0] ), t )
    t = numpy.dot( vtk.rotation_matrix( numpy.radians(-normal[1]), [0,0,1] ), t )

    # apply 3DAVG refinement transformation (rotation only, R2)
    t = numpy.dot( refinement_rotation, t )

    #  compute translation due to change in rotation origin for R1 * R2
    r = vtk.rotation_matrix( numpy.radians(-normal[2]), [0,0,1] )
    r = numpy.dot( vtk.rotation_matrix( numpy.radians(-normal[0]), [1,0,0] ), r )
    r = numpy.dot( vtk.rotation_matrix( numpy.radians(-normal[1]), [0,0,1] ), r )
    r = numpy.dot( refinement_rotation, r )

    # Difference vector between rotation origins C1=[51,50,50] and C2=[51,51,50]
    diff = [ 0, 1, 0, 1 ]

    # Compute: t = Rot(C1) wrt C2 - C1 = Rot * ( C1 - C2 ) - ( C1 - C2 )
    tdiff = numpy.dot( r, diff ) - diff

    # print 'tdiff', tdiff

    # compute post-multiplying translation component from refinement matrix, T2
    f = numpy.dot( numpy.linalg.inv( refinement_rotation ), refinement )

    # apply rotation center and refinement translations
    t = numpy.dot( vtk.translation_matrix( [-tdiff[0,0]+f[0,3], -tdiff[0,1]-f[1,3], -tdiff[0,2]+f[2,3] ] ), t )
    # t = numpy.dot( vtk.translation_matrix( -tdiff[0,:3] ), t )

    # t = numpy.dot( vtk.translation_matrix([0,0,-cutOffset]), t )

    # extract euler angles from transformation matrix
    if abs(t[2,2]) < 1 - numpy.nextafter(0,1):
        theta = math.acos( t[2,2] )
        psi = math.atan2( t[2,1] / math.sin(theta), -t[2,0] / math.sin(theta) )
        phi = math.atan2( t[1,2] / math.sin(theta), t[0,2] / math.sin(theta) )
    else:
        theta = math.pi * ( 1 - numpy.sign( t[2,2] ) ) / 2
        phi = 0
        psi = math.acos( t[0,0] / math.cos(theta) )

    frealign = numpy.degrees( numpy.array( [ psi, theta, phi ] ) )

    # frealign does not use negative angles, so we add 360 to each negative angle
    frealign = numpy.where( frealign < 0, frealign + 360.0, frealign )

    # Now project xyz shifts onto view plane for use in FREALIGN
    # rt = vtk.rotation_matrix( phi,[0,0,1] )
    # rt = numpy.dot( rt, vtk.rotation_matrix( theta,[0,1,0]) )
    # rt = numpy.dot( rt, vtk.rotation_matrix( psi,[0,0,1]) )
    # nt = numpy.dot( numpy.linalg.inv(rt), t )

    nt = - numpy.linalg.inv(t)

    return [ frealign[2], frealign[1], frealign[0], nt[0,3], nt[1,3] ]

