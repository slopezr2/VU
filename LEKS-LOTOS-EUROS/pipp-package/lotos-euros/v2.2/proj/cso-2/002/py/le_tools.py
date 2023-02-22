#! /usr/bin/env python3

"""
Classes for various LOTOS-EUROS tasks.
"""

########################################################################
###
### output postprocessing
###
########################################################################

class PostProc( object ) :

    """
    Postprocess, in output directory:
    
    * create collection of source files and settings;
    * copy runtime settings;
    * pack netcdf files.
    """
    
    def __init__( self, rcfile ) :
    
        """
        Initialize class.
        """
        
        # modules:
        import os
        import shutil
        import subprocess
        
        # tools:
        import rc
        
        # info ...
        print( 'postprocess ...' )
        
        # read settings:
        rcf = rc.RcFile( rcfile )
        
        # output dir:
        outdir = rcf.get( 'le.output.outdir' )

        #
        # * collect sources and settings
        #
        
        # target file:
        tarfile = os.path.join( outdir, 'build.tar.gz' )
        # info ...
        print( 'collect source files ...' )
        print( '  %s ...' % tarfile )
        
        # build directory:
        bdir = rcf.get( 'le.copy.prefix' )
        # parent directory:
        pdir = os.path.join(bdir,os.pardir)
        # info ..
        print( '  change to %s ...' % pdir )
        # one level back:
        os.chdir( pdir )
        
        # collection command:
        command = [ 'tar', '--create', '--gzip' ]
        command.append( '--file=%s' % tarfile )
        command.append( '--exclude=__pycache__' )
        command.append( '--exclude="*.o"' )
        command.append( '--exclude="*.mod"' )
        command.append( '--exclude="*.x"' )
        command.append( 'build/bin' )
        command.append( 'build/data' )
        command.append( 'build/py' )
        command.append( 'build/rc' )
        command.append( 'build/src' )
        
        # run:
        subprocess.check_call( command )

        
        #
        # * copy runtime rcfile
        #
        
        # info ...
        print( 'copy runtime settings ..' )
        
        # info ...
        print( '  change to %s ...' % outdir )
        # goto output dir:
        os.chdir( outdir )

        # runtime settings:
        fname = rcf.get( 'le.copy.rcwrite' )
        # info ...
        print( '  copy %s ...' % fname )
        
        # copy:
        shutil.copy( fname, os.path.basename(fname) )
        
        
        #
        # * pack output
        #
        
        # info ..
        print( 'pack netcdf output ...' )

        # info ...
        print( '  change to %s ...' % outdir )
        # goto output dir:
        os.chdir( outdir )
        
        # info ...
        print( '  pack ...' )
        # list:
        fnames = os.listdir(os.curdir)
        # loop:
        for fname in fnames :
            # filter:
            if not fname.endswith('.nc') : continue
            # info ...
            print( '    %s ...' % fname )
            # temporary name:
            tmpfile = fname+'_'
            # packing including deflation, this requires netcdf4:
            command = [ 'ncpdq', '--netcdf4', '--deflate', '1', '--output', tmpfile, fname ]
            # ALTERNATIVE?
            # Use 'Precision-Preserving Compression', see:
            #    http://nco.sourceforge.net/nco.html#Precision_002dPreserving-Compression
            #command = [ 'ncks', '--ppc', 'default=3', '--output', tmpfile, fname ]
            # run:
            subprocess.check_call( command )
            # rename:
            os.rename( tmpfile, fname )
        #endfor
        
    #enddef __init__
    
#endclass

