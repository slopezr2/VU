
"""
Build LEIP executables.

Usage:

  ! build exectubale for selecting grid cells,
  ! path to exectuble is returned:
  select_x = Build( rcf, 'select' )
  
  ! run executable:
  subprocess.check_call( [select_x,...] )
  
"""

########################################################################
###
### common
###
########################################################################

# common stuff:
from leip_common import *


#####################################################################
###
### main
###
#####################################################################

def Build( rcf, target ) :

    """
    Build LEIP executable with nickname 'target' .
    
    Description:
    
        The routine will read settings from the rcfile to build
        a requested executable. The settings are either general,
        or specific for the executable or for the host where the
        script is running.

        First a temporary build directory is created,
        for example on a scratch space.

        The source files are then copied into the build directory.

        The name of the executable to be build is read from the rc file,
        for example 'target.x'

        One of the source files is 'Makefile', which will be used
        to compile the executable:
            make target.x
        The 'Makefile' should include host-specific settings from
        an include file with a pre-defined name:
            include Makefile_settings
        See below for the content, which should be used in the Makefile.
        
        
    Makefile settings
    
        Sample content of "Makefile_settings" written by this routine.
        The content is written based on settings read from the rcfile.
            
            # compiler
            FC = gfortran -std=f2003
            LD = gfortran

            # compiler
            FFLAGS = -ffree-line-length-none -fimplicit-none -fdefault-real-8 -O0 -fcheck=all
            LDFLAGS = 

            NETCDF_INCS = -I/opt/netcdf4/4.1.3/include
            NETCDF_LIBS = -L/opt/netcdf4/4.1.3/LP64/lib -lnetcdff -lnetcdf

    
    Return value:
    
        '/path/to/target.x' : fulll path to executable that was build
      
    Rcfile settings:
    
      Main settings:
      
        ! source location:
        leip.prefix                   :  ${LE}/tools/leip/vtrunk

        ! build directory:
        leip.build.prefix             :  ${PWD}/build

        ! subdirs to be copied:
        leip.build.copy.dirs          :  src

        ! names of executables to be used ;
        ! will be searched in <prefix>/bin :
        leip.build.x.select           :  c3po_select.x
        leip.build.x.spg2grad         :  c3po_spectral_gradient.x
        leip.build.x.spg2rgg          :  c3po_spg2rgg.x

        ! TNO HPC settings:
        #include ${leip.prefix}/rc/host-tno-hpc.rc

        ! flag groups:
        !leip.build.compiler.flags           :  wp-r8 optim-none check
        leip.build.compiler.flags           :  wp-r8 optim-fast     


      Host specific settings in for example "host-tno-hpc.rc" :

        ! compiler settings:
        #include ${leip.prefix}/rc/compiler-gfortran.rc

        ! names for libraries to be used,
        ! for each define library specific fortran and linker flags:
        build.libs                     :  netcdf

        ! lib specific flags:
        build.lib.netcdf.fflags        :  -I${NETCDF4_FORTRAN_HOME}/include
        build.lib.netcdf.ldflags       :  -L${NETCDF4_FORTRAN_HOME}/lib -lnetcdff -L${NETCDF4_HOME}/lib -lnetcdf


      Compiler settings in for example "compiler-gfortran.rc" :
      
        compiler.fc                             :  gfortran
        compiler.ld                             :  gfortran

        compiler.flags.default.fflags           :  -ffree-line-length-none -fimplicit-none
        compiler.flags.default.ldflags          :  

        compiler.flags.wp-r8.fflags             :  -fdefault-real-8
        compiler.flags.wp-r8.ldflags            :  

        compiler.flags.check.fflags             :  -fbounds-check -ffpe-trap=invalid,zero -Wall
        compiler.flags.check.ldflags            :  

        compiler.flags.optim-none.fflags        :  -O0
        compiler.flags.optim-none.ldflags       :  

        compiler.flags.optim-fast.fflags        :  -O3
        compiler.flags.optim-fast.ldflags       :  

    """
    
    # modules:
    import os
    import shutil
    import subprocess
    
    # info ...
    logger.info( 'build "%s" ...' % target )
    
    # original:
    leip_prefix = rcf.get( 'leip.prefix' )
    
    # build directory:
    build_prefix = rcf.get( 'leip.build.prefix' )
    
    # create if necessary:
    if not os.path.isdir(build_prefix) :
        logger.info( '  create build directory %s ...' % build_prefix )
        os.makedirs( build_prefix )
    #endif
    
    # source dirs to be copied:
    subdirs = rcf.get( 'leip.build.copy.dirs' ).split()
    # loop:
    for subdir in subdirs :
        # full source path:
        sdir = os.path.join( leip_prefix, subdir )
        # check ...
        if not os.path.isdir( sdir ) :
            logger.error( 'source dir not found : %s' % sdir )
            raise Exception
        #endif
        # files:
        fnames = os.listdir( sdir )
        # loop:
        for fname in fnames :
            # full path to source:
            sfile = os.path.join( sdir, fname )
            # skip directories:
            if os.path.isdir(sfile) : continue
            # full path to target file:
            tfile = os.path.join( build_prefix, fname )
            # copy if source is newer than destination:
            if (not os.path.isfile(tfile)) or (os.path.getmtime(sfile) > os.path.getmtime(tfile)) :
                # info ...
                logger.info( '  copy %s ...' % os.path.join(subdir,fname) )
                # copy:
                shutil.copy( sfile, tfile )
            #endif
        #endfor  # files
    #endfor  # subdirs
    
    # save current dir:
    cwd = os.getcwd()
    # goto build dir:
    os.chdir( build_prefix)
    
    # makefile with machine specific settings:
    mfile = 'Makefile_settings'
    
    # init content:
    lines = []
    lines.append( '# Make settings' )
    lines.append( '' )

    # add lines:
    lines.append( '# compiler' )
    lines.append( 'FC = %s' % rcf.get('compiler.fc') )
    lines.append( 'LD = %s' % rcf.get('compiler.ld') )
    lines.append( '' )

    # flag descriptions:
    flagnames = ['default'] + rcf.get( 'leip.build.compiler.flags' ).split()
    # init flags:
    fflags = ''
    ldflags = ''
    # add extra:
    for flagname in flagnames :
        # extend:
        fflags  = fflags  + ' ' + rcf.get( 'compiler.flags.%s.fflags' % flagname )
        ldflags = ldflags + ' ' + rcf.get( 'compiler.flags.%s.ldflags' % flagname )
    #endfor
    # add lines:
    lines.append( '# compiler' )
    lines.append( 'FFLAGS = %s' % fflags )
    lines.append( 'LDFLAGS = %s' % ldflags )
    lines.append( '' )
    
    # libs to be defined:
    libs = rcf.get( 'build.libs' ).split()
    # loop:
    for lib in libs :
        # extract:
        fflags  = rcf.get( 'build.lib.%s.fflags' % lib )
        ldflags = rcf.get( 'build.lib.%s.ldflags' % lib )
        # fill:
        lines.append( '%s_INCS = %s' % (lib.upper(),fflags) )
        lines.append( '%s_LIBS = %s' % (lib.upper(),ldflags) )
        lines.append( '' )
    #endfor
    
    # write:
    f = open( mfile, 'w' )
    for line in lines : f.write( line+'\n' )
    f.close()
    
    # target executalble
    target_x = rcf.get( 'leip.build.x.%s' % target )
    # make:
    command = [ 'make', target_x ]
    # run:
    subprocess.check_call( command )
    
    # back:
    os.chdir( cwd )
    
    # ok
    return os.path.join(build_prefix,target_x)

#enddef Build
