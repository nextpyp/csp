#!/opt/apps/rhel7/anaconda2/bin/python -B

import argparse, os, sys, numpy, subprocess, mrc, csv, commands

# Script for running frealing refinement called from frealign_ali.sh

def param(parameter,iteration):
    listed = parameter.split(':')
    return listed[ min(iteration-2,len(listed)-1) ]

def load_parameters( parameter_file):
    if os.path.isfile(parameter_file):
        reader = csv.reader( open( parameter_file, 'rb'), delimiter='\t' )
        return dict(x for x in reader)
    else:
        print 'File not found', parameter_file
        return 0

if __name__ == '__main__':

    # check and parse arguments    

    if len(sys.argv) != 6:
        print 'Usage: {0} dataset image_stack iteration reference dot_par_in'.format( sys.argv[0] )
        sys.exit()

    dataset = sys.argv[1]
    image_stack = sys.argv[2]
    i = int(sys.argv[3])
    reference = sys.argv[4]
    dot_par_in = sys.argv[5]
    
    # frealign parameters    
    fp = load_parameters('/scratch/frealign.config')

    # microscope parameters    
    mp = load_parameters('/scratch/.pyp_config')

    pixel_size = float( mp['scope_pixel'] ) * float( mp['data_bin'] ) * float( mp['particle_bin'] )
    dstep = pixel_size * float(mp['scope_mag']) / 10000.0
    
    # retrieve number of current inner iteration
    inner_iter = int( open('/scratch/.csp_current_inner_iteration').read() )
    
    # retrieve resolution limits for refinement
    cp = load_parameters('/scratch/parameters.config')
    rlowref = param( cp['LowResolutionForRefinement'], inner_iter )
    rref = param( cp['HighResolutionForRefinement'], inner_iter )
        
    # Append resolution table to end of file 
    res_table = '/scratch/resolution_table'
    if os.path.exists( res_table ):
        com = 'cat %s >> %s' % ( res_table, dot_par_in ) 
        # commands.getoutput(com)

    particles = int( mrc.readHeaderFromFile( image_stack )['nz'] )
        
    if float(param(fp['rrec'],i)) > 0:
	res_rec = param(fp['rrec'],i)
    else:
	res_rec = 2 * pixel_size
        
    label = dataset + '_%02d' % i    
    logfile = label + '_msearch_n.log1'
    logfile = '/dev/null'
    
    # run FREALIGN

    if 'cc' in param( fp['metric'], inner_iter ).lower():

        command="""
%s/frealign_v9.10/bin/frealign_v9.exe << eot >>%s 2>&1
M, 1, %s, %s, %s, %s, %s, T, %s, %s, %s, 0, F, %s, %s				!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FDUMP,IMEM,INTERP
%s, 0., %s, %s, %s, %s, %s, %s, %s, %s, %s					!RO,RI,PSIZE,MW,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
0,0,0,1,1									!MASK
1, %i										!IFIRST,ILAST 
C1
1.0, %s, %s, %s, %s, %s, 0., 0.
%s, %s, %s, %s, %s, %s
%s                                                                              !IMAGE STACK
%s_match.mrc
%s
%s
/dev/null
-100., 0., 0., 0., 0., 0., 0., 0.						!terminator with RELMAG=-100.0 to skip 3D reconstruction
%s
%s_weights
%s_map1.mrc
%s_map2.mrc
%s_phasediffs
%s_pointspread
eot
"""  % ( os.environ['PYP_DIR'], logfile,
    param(fp['fmag'],i), param(fp['fdef'],i), param(fp['fastig'],i), param(fp['fpart'],i), param(fp['iewald'],i), param(fp['ffilt'],i), param(fp['fbfact'],i), param(fp['fmatch'],i), param(fp['imem'],i), param(fp['interp'],i),
    mp['particle_rad'], pixel_size, mp['particle_mw'], mp['scope_wgh'], param(fp['xstd'],i), param(fp['pbc'],i), param(fp['boff'],i), param(fp['dang'],i), param(fp['itmax'],i), param(fp['ipmax'],i),
    particles,
    dstep, param(fp['target'],i), param(fp['threc'],i), mp['scope_cs'], mp['scope_voltage'],
    res_rec, rlowref, rref, rref, param(fp['dfsig'],i), param(fp['rbfact'],i),
    image_stack,
    label,
    dot_par_in,
    dot_par_in,
    reference,
    label,
    label,
    label,
    label,
    label )

    else:

        ref = 1

        my_mask = param( fp['mask'], i ).split(',')
        if my_mask[0] == '0':
            psib = 'no'
        else:
            psib = 'yes'
        if my_mask[1] == '0':
            thetab = 'no'
        else:
            thetab = 'yes'
        if my_mask[1] == '0':
            phib = 'no'
        else:
            phib = 'yes'
        if my_mask[3] == '0':
            shxb = 'no'
        else:
            shxb = 'yes'
        if my_mask[4] == '0':
            shyb = 'no'
        else:
            shyb = 'yes'

        if 't' in param( fp['fssnr'], i ).lower() and os.path.exists( '/scratch/statistics_r%02d.txt' % ref ) and len( open( '/scratch/statistics_r%02d.txt' % ref ).read() ) > 0:
            stats = 'yes'
        else:
            stats = 'no'

        if param( fp['focusmask'], i ) == '0,0,0,0':
            masking = 'no'
        else:
            masking = 'yes'
        if 't' in param( fp['fdef'], i ).lower():
            defocus = 'yes'
        else:
            defocus = 'no'
        if 't' in param( fp['fmatch'], i ).lower():
            match = 'yes'
        else:
            match = 'no'
        if 't' in param( fp['invert'], i ).lower():
            invert = 'yes'
        else:
            invert = 'no'

        globally = 'no'
        locally = 'yes'

        # set radius for global search
        if param(fp['srad'],i) == '0':
            srad = str( 1.5 * float(mp['particle_rad']) )
        else:
            srad = param(fp['srad'],i)

        boost = '30.0'
        if 't' in param( fp['fboost'], i ).lower():
            boost = param( fp['fboostlim'], i )

        if 'new' in param( fp['metric'], inner_iter ).lower():

            command = \
'{0}/frealign_v9.11/bin/refine3d << eot >>{1} 2>&1\n'.format( os.environ['PYP_DIR'], logfile ) + \
image_stack + "\n" \
'{0}\n{1}\n'.format( dot_par_in, reference ) + \
'/scratch/statistics_r%02d.txt\n' % ref + \
stats + '\n' + \
'{0}_match.mrc\n{1}\n/dev/null\n'.format( label, dot_par_in ) + \
param(fp['symmetry'],i) + '\n' + \
'{0}\n{1}\n'.format( 1, particles ) + \
'{0}\n{1}\n{2}\n{3}\n{4}\n'.format( pixel_size, mp['scope_voltage'], mp['scope_cs'], mp['scope_wgh'], mp['particle_mw'] ) + \
mp['particle_rad'] + '\n' + \
rlowref + '\n' + \
rref + '\n' + \
boost + '\n' + \
param(fp['rhcls'],i) + '\n' + \
srad + '\n' + \
rref + '\n' + \
param(fp['dang'],i) + '\n' + \
'20\n' + \
param(fp['searchx'],i) + '\n' + \
param(fp['searchy'],i) + '\n' + \
'\n'.join( param(fp['focusmask'],i).split(',') ) + '\n' + \
'500.0\n' + \
'50.0\n' + \
param(fp['iblow'],i) + '\n' + \
globally + '\n' + locally + '\n' + \
psib + '\n' + thetab + '\n' + phib + '\n' + shxb + '\n' + shyb + '\n' + \
match + '\n' + masking + '\n' + defocus + '\n' + invert + '\n' + \
'eot\n'

        else:
            
            # frealignx

            if 't' in param( fp['norm'], i ).lower():
                normalize = 'yes'
            else:
                normalize = 'no'

            exclude = 'no'
            normalize_input = 'no'
            threshold_input = 'no'

            command = \
'{0}/frealignx/refine3d << eot >>{1} 2>&1\n'.format( os.environ['PYP_DIR'], logfile ) + \
image_stack + "\n" + \
'{0}\n{1}\n'.format( dot_par_in, reference ) + \
'/scratch/statistics_r%02d.txt\n' % ref + \
stats + '\n' + \
'{0}_match.mrc\n{1}\n/dev/null\n'.format( label, dot_par_in ) + \
param(fp['symmetry'],i) + '\n' + \
'{0}\n{1}\n1\n'.format( 1, particles ) + \
'{0}\n{1}\n{2}\n{3}\n{4}\n0\n'.format( pixel_size, mp['scope_voltage'], mp['scope_cs'], mp['scope_wgh'], mp['particle_mw'] ) + \
mp['particle_rad'] + '\n' + \
rlowref + '\n' + \
rref + '\n' + \
boost + '\n' + \
param(fp['rhcls'],i) + '\n' + \
srad + '\n' + \
rref + '\n' + \
param(fp['dang'],i) + '\n' + \
'20\n' + \
param(fp['searchx'],i) + '\n' + \
param(fp['searchy'],i) + '\n' + \
'\n'.join( param(fp['focusmask'],i).split(',') ) + '\n' + \
'500.0\n' + \
'50.0\n' + \
param(fp['iblow'],i) + '\n' + \
globally + '\n' + locally + '\n' + \
psib + '\n' + thetab + '\n' + phib + '\n' + shxb + '\n' + shyb + '\n' + \
match + '\n' + masking + '\n' + defocus + '\n' + normalize + '\n' + invert + '\n' + exclude + "\n" + normalize_input + '\n' + threshold_input + "\n" \
'eot\n'

    # with open(logfile,'a') as f:
    #    f.write(command)

    # print command
    subprocess.Popen( command, shell=True).wait()
    
