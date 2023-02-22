

"""
PYCASSO User Scripts

Do not change the arguments of the defined routines, only their content !
"""

#-------------------------------------------------
# arguments and options 
#-------------------------------------------------

def DefineOptions( parser ) :
    
    """
    Usage : DefineOptions( parser )
    Define options supported by user scripts.
    Arugments:
       parser   : optparse object 
    """
    
    ## define flag:
    #parser.add_option( "-c", "--clean", 
    #                     help="remove object and module files before compilation",
    #                     dest="clean", action="store_true", default=False )

    # ok
    return
    
#enddef

# *

def StoreOptions( settings, values ) :
        
    """
    Add the parsed flag to a dictionairy.
    The values have data fields with names defined by 'dest'
    in the previous calls to 'parser.add_option' .
    """
    
    # translate options into a dictionairy if they should
    # replace rcfile values:
    #if values.clean     : settings['build.make.clean'    ] = True

    # ok
    return
    
#enddef


#-------------------------------------------------
# source 
#-------------------------------------------------


def Build_FlagGroups( rcf, basic=False, nocheck=False ) :

    """
    Return list of compiler flag groups to be used.
    """
    
    # default :
    flaggroups = ['default']
    
    # add openmp ?
    if rcf.get('par.openmp','bool') : flaggroups.append('openmp')
    
    # add MPI ?
    if rcf.get('par.mpi','bool') : flaggroups.append('mpi')

    # include not standard flags ?
    if not basic :
        # read other:
        line = rcf.get('build.configure.flags')
        # add if necessary:
        if len(line) > 0 :
            # split at white space:
            fields = line.split()
            # loop over new groups:
            for field in fields :
                # skip ?
                if nocheck and ('check' in field) : continue
                # add:
                flaggroups.append( field )
            #endfor
        #endif
    #endif

    # ok
    return flaggroups

#enddef


# ***


def Build_Define( rcf, macro_group, mdefs ) :

    """
    Edit a list with macro's to be defined given
    the name of group of macro's and the other setings
    in the rcfile.
    """
    
    # external:
    import logging
    
    # info ...
    logging.info( '  user script Build_Define for macro group %s ...' % macro_group )
        
    # ok
    return mdefs
    
#enddef


# ***


def Build_Configure( rcf ) :

    """
    Configure a source code.
    This script is called from the source directory.
    Arguments:
      rcf  : dictionairy with settings from rc file
    """
    
    # external
    import logging
    
    # info ...
    logging.info( 'User script Build_Configure ...' )
    
    # generate source files if necessary:
    Build_Configure_Generate( rcf )
    
    # check on depricated stuff ...
    Build_Configure_Check( rcf )

    # ok
    return
    
#enddef   # Build_Configure


# *


def Build_Configure_Generate( rcF ) :

    """
    Generate source files.
    """
    
    # modules:
    import logging

    # generate a source ?
    flag = rcF.get('build.configure.genes','bool',default=False)
    if flag :
        # info ...
        logging.info( '  generate source files ...' )
        # tools:
        import genes
        # generate source:
        genes.GeneS( rcF, 'genes' )
    else :
        # info ...
        logging.info( '  no source files to be generated ...' )
    #endif
    
    # ok
    return

#enddef    # Build_Configure_Generate

    
# *


def Build_Configure_Check( rcf ) :

    """
    Check source file for undesired features.
    """
    
    # external:
    import logging
    import os
    import fnmatch

    # info ...
    logging.info( '    user script Build_Configure_Check ...' )
    
    # keywords for checks to be performed:
    checknames = rcf.get( 'build.configure.checks', default='' )
    # some defined ?
    if len(checknames) == 0 : return
    
    # error or just warnings ?
    with_error = rcf.get( 'build.configure.checks.error', 'bool', default=False )
    
    # set flag:
    any_warning = False
    
    # list files:
    srcfiles = os.listdir( os.curdir )
    srcfiles.sort()
    
    # loop over checks:
    for checkname in checknames.split() :

        # some specials are done below:
        if checkname == 'unknown_macro' : continue
    
        # paterns:
        test_msg   = rcf.get( 'build.configure.check.%s.msg'   % checkname )
        test_files = rcf.get( 'build.configure.check.%s.files' % checkname ).split()
        test_skip  = rcf.get( 'build.configure.check.%s.skip'  % checkname ).split()
        test_line  = rcf.get( 'build.configure.check.%s.test'  % checkname, 'str' )
        test_help  = rcf.get( 'build.configure.check.%s.help'  % checkname )

        # set flags:
        matching_files = False

        # loop over files:
        for srcfile in srcfiles :

            # match with patern ?
            match = False
            for pat in test_files :
                match = match or fnmatch.fnmatch(srcfile,pat)
                if match : break
            #endfor
            if not match : continue 

            # ... except if the name matches other patterns:
            match = False
            for pat in test_skip :
                match = match or fnmatch.fnmatch(srcfile,pat)
                if match : break
            #endfor
            if match : continue
            
            # always true ? trap this immediately,
            # this is probably a check on certain file names:
            if test_line == 'True' :

                # match found:
                match = True

            else: 

                # read file:
                f = open( srcfile )
                lines = f.readlines()
                f.close()

                #~ complain if none of the lines matches:
                if test_line.startswith('not') :
                    # by default no match:
                    match = False
                    # loop over lines:
                    for line in lines :
                        # test on this line:
                        match = match or (not eval( test_line ))
                        # try next file after first match ...
                        if match : break
                    #endfor
                    # check next file if the requested code was found:
                    if match : continue
                    # revert:
                    match = not match
                #~ complain if one of the lines matches:
                else :
                    # by default no match:
                    match = False
                    # loop over lines:
                    for line in lines[0:20] :
                        # test on this line:
                        match = match or eval( test_line )
                        # leave after first match ...
                        if match : break
                    #endfor
                #endfor
                
            #endif  # filename test or content test

            # found something ?
            if match :
                # info ...
                if not matching_files : logging.warning( '      %s : [found]' % test_msg )
                logging.warning( '        %s' % srcfile )
                # reset flags:
                matching_files = True
                any_warning = True
            #endif

        #endfor   # source files

        # info ...
        if matching_files :
            # display error message ?
            if with_error :
                # display help text; split at '\n' for newlines:
                for helpline in test_help.split('\\n') : logging.warning(helpline)
            #endif
        else :
            # no warnings for this test ...
            logging.info( '      %s [none ]' % test_msg )
        #endif

    #endfor  # checks
    
    # check for unknown macro's ?
    checkname = 'unknown_macro'
    if checkname in checknames :

        # settings:    
        test_msg   = rcf.get( 'build.configure.check.%s.msg' % checkname )

        # names of macro groups:
        macgroups = rcf.get( 'build.configure.macro.groups' ).split()
        # collect all supported macro's:
        macall = []
        for macgroup in macgroups :
            macs = rcf.get( 'build.configure.macro.%s.all' % macgroup ).split()
            macall = macall + macs
        #endfor

        # flag ...
        logged_msg = False

        # loop over files:
        for srcfile in srcfiles :

            # skip non-source files:
            if os.path.isdir(srcfile) : continue
            if fnmatch.fnmatch(srcfile,"*.mod"   ) : continue
            if fnmatch.fnmatch(srcfile,"*.o"     ) : continue
            if fnmatch.fnmatch(srcfile,"*.x"     ) : continue
            if fnmatch.fnmatch(srcfile,"Makefile") : continue

            # read file:
            f = open( srcfile, 'r' )
            lines = f.readlines()
            f.close()

            # flags:
            logged_srcfile = False

            # loop over lines:
            for iline in range(len(lines)) :
                # current:
                line = lines[iline].strip()
                # macro test ?
                if line.startswith('#ifdef') or line.startswith('#ifndef') :
                    # second element of line is macro name:
                    mac = line.split()[1].strip()
                    # not supported ?
                    if mac not in macall :
                        # test description if not done yet:
                        if not logged_msg :
                            logging.info( '      %s' % test_msg )
                            logged_msg = True
                        #endif
                        # intro if necessary:
                        if not logged_srcfile :
                            logging.error( '        unsupported macro(s) in %s :' % srcfile )
                            logged_srcfile = True
                        #endif
                        # line number and content:
                        logging.error( '        %6i : %s' % (iline,line) )
                        # set flag:
                        any_warning = True
                    #endif   # unsuported macro
                #endif  # line with macro test  
            #endfor   # lines

        #endfor   # source files

        # jippy ...
        if not logged_msg :
            # no warnings for this test ...
            logging.info( '      %s [none ]' % test_msg )
        #endif

    #endif   # test on unsupported macro's

    # break ?
    if any_warning and with_error :
        logging.error( 'some source code checks failed; break' )
        logging.error( '(set "build.configure.checks.error : False" in the expert.rc to avoid this error)' )
        raise Exception
    #endif
    
    # ok
    return

#enddef

    
# ***


def Build_Compiler( rcf ) :

    """
    Set compiler and linker names.
    Usually it is enough to read the name from the rcfile,
    but some compiler families have aliases for compilation with
    MPI or OpenMP enabled.
    Arguments:
      rcf  : dictionairy with settings from rc file
    Return values:
      fc,linker
    """
    
    # external
    import logging
    
    # info ...
    logging.info( '  user script Build_Compiler ...' )
    
    # extract compiler name:
    fc = rcf.get('compiler.fc')
    # or supporting openmp ?
    if rcf.get('par.openmp','bool') : fc = rcf.get( 'compiler.fc.openmp' )

    # or with mpi support ?
    if rcf.get('par.mpi','bool') : 
        # extract compiler name:
        fc = rcf.get( 'mpi.compiler.fc' )
        # or supporting openmp ?
        if rcf.get('par.openmp','bool') : fc = rcf.get( 'mpi.compiler.fc.openmp' )
    #endif

    # f77 compiler, by default the same as fc:
    f77 = rcf.get( 'compiler.f77', default=fc )
    
    # assume linker is the same:
    linker = fc

    # info ...
    logging.debug( '    fortran compiler : %s' % fc )
    logging.debug( '    f77 compiler     : %s' % f77 )
    logging.debug( '    linker           : %s' % linker )
    
    # ok
    return fc,f77,linker
    
#enddef


# ***


def Build_Make( rcf ) :

    """
    Make and install an executable.
    This script is called from the source directory.
    Arguments:
      rcf  : dictionairy with settings from rc file
    """
    
    # external
    import sys
    import os
    import logging
    
    # tools:
    import go
    
    # module dir ?
    mdir = rcf.get('compiler.mdir',default='None')
    if mdir != 'None' :
        # not present yet ? then create:
        if not os.path.exists( mdir ) : os.makedirs( mdir )
    #endif

    # info ...
    logging.debug( '  make ...' )

    # number of jobs available for make:
    build_jobs = rcf.get( 'build.jobs', default='' )

    # get maker command; replace some keys:
    maker = rcf.get('maker').replace('%{build.jobs}',build_jobs)
    
    # get name of Makefile supporting pycasso stuff,
    # e.g. including Makefile_flags, Makefile_deps etc:
    pycasso_makefile = rcf.get( 'build.make.makefile' )

    # get list of executables(mainfile) items, e.g.:
    #       tm5.x tracer ; boe.x
    items = rcf.get('build.make.exec')
    # loop over items:
    for item in items.split(';') :
        # might be of form 'exec mainfile', executable is always first item:
        exe = item.split()[0]
        # full make command:
        command = maker.split()+['-f',pycasso_makefile,exe]
        # info ...
        logging.debug( '  run command: %s' % str(command) )
        # run:
        p = go.subprocess.watch_call( command )
    #endfor
    
    # ok
    return
    
#enddef

