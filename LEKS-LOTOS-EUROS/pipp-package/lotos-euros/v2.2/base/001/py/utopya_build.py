
#-------------------------------------------------
# help
#-------------------------------------------------

"""
***********************
``utopya_build`` module
***********************

The :py:mod:`utopya_build` module provides classes to copy a source
code to a work directory and build (compile) an executable.
All classes are accesible through the top level :py:mod:`utopya` module.

* The :py:class:`.UtopyaCopy` class is used
to collect source code files and copy them into a work directory.
* The :py:class:`.UtopyaConfigure` is used to configure the source
code and create make files.
This is done using instances from two other classes defined in this module:

* The :py:class:`.UtopyaCompilerSettings` class is used to provide information on the
  compiler to be used, for example the compiler executable name(s) and
  useful compiler flags.
* The :py:class:`.UtopyaLibSettings` class is used to provide library names and paths
  on the target machine.

* The :py:class:`UtopyaMake` class compiles the source code and creates an executable.

The configuration of the objects is done using rcfile settings.


.. Label between '.. _' and ':' ; use :ref:`text <label>` for reference
.. _example-build:

Example script and settings
===========================

A script is provided to show how to build an executable from a Fortran source.
The example builds the ``SLM`` (Simple Linear Model) provided with UTOPyA.
The script uses UTOPyA Python classes from the :py:mod:`utopya_build` module
to first copy the source code to a work directory,
and then to configure and compile the code.

The exampe script is available as:

  `bin/run-example-build <../../../bin/run-example-build>`_

The example settings to build SLM are provided in:

  `rc/example-build-slm.rc <../../../rc/example-build-slm.rc>`_

These settings are loaded by default if the script is run without arguments::

  ./bin/run-example-build

where an example with explicit arguments would be::

  ./bin/run-example-build rc/example-build-slm.rc --rcbase=slm

By changing the name of the rcfile and the base key for the settings
the script could be used to build another souce code too.

Class hierarchy
===============

The classes provided by this module have been derived with the following hierchy:

* :py:class:`.UtopyaRc`

  * :py:class:`.UtopyaCopy`
  * :py:class:`.UtopyaConfigure`
  * :py:class:`.UtopyaCompilerSettings`

    * :py:class:`.UtopyaCompilerSettingsFortran`
    * :py:class:`.UtopyaCompilerSettingsF2py`

  * :py:class:`.UtopyaLibSettings`
  * :py:class:`.UtopyaDependencies`
  * :py:class:`.UtopyaMake`

Classes
=======

"""


#-------------------------------------------------
# modules
#-------------------------------------------------

import utopya_rc


#-------------------------------------------------
# copy
#-------------------------------------------------

class UtopyaCopy( utopya_rc.UtopyaRc ) :

    """
    Base class for objects that collect and copy source codes into a work directory.

    For initialization, provide the name of the rcfile with the settings
    that define the tasks to be done.
    Optionally, provide an 'rcbase' as prefix for the settings to be used,
    and an environment dictionairy 'env' for variable substitutions;
    see the initialization of the :py:class:`.UtopyaRc` class for details
    and the necessary configuration of the logging system.

    If the 'new' flag is present it overwrites the '\*.copy.new' flag;
    set to 'True' to remove an existing build directory if present.

    The following example shows the required settings for an rcbase 'optm';
    the comment should describe the purpose of each setting::

        ! List of source directories to be copied;
        ! files in later directories overwrite previously copied versions:
        optm.copy.dirs                :  sources/model/base  sources/model/special

        ! Sub directories to be maintained ;
        ! leave empty if no subdirs are defined:
        optm.copy.subdirs             :  src py

        ! directories to be inlcuded in copy,
        ! otherwise only files are copied:
        optm.copy.incdirs             :  include

        ! Skip files matching these filename patterns
        ! (tempoary editor files, compiled modules, etc)
        optm.copy.skip                :  .#* *~ *.pyc

        ! An optional flag is present to remove particular extensions
        ! from the copied filename. These extensions start with '__',
        ! and are used in certain application to distinguish different
        ! source files from eachother. With this flage enabled,
        ! the source file "work__special.F90" will be copied to "work.F90":
        optm.copy.remove__ext         :  True

        ! Prefix for destination of source and script files
        ! (base path for subdirectories src, py, etc) ;
        ! this should be an absolute path,
        ! use ${PWD} for the present dir if necessary:
        optm.copy.prefix              :  ${SCRATCH}/model/build

        ! The above defined refix could be extended with extra keywords,
        ! fore example 'build_optim-all_mpi' instead of just 'build'.
        ! While develloping it is sometimes useful to have different
        ! builds next to eachother to avoid a complete rebuild when
        ! a simple switch from optimal to debug compiler flags is needed,
        ! or one would like to test compilation with and without MPI.
        ! Provide here a space separated list of extension keywords,
        ! usually these are copied from other settings sush as compile options:
        optm.copy.prefix.extensions   :  optim-all mpi

        ! remove existing build directory if present ?
        optm.copy.new                 :  False

        ! optionally define a filename to which the evaluated rcfile settings will be written:
        optm.copy.rcwrite             :  ${SCRATCH}/model/build/optm-runtime.rc

    """

    def __init__( self, rcfile, rcbase='', env={}, new=None ) :

        # modules:
        import os

        # tools:
        import utopya_tools

        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )

        # info ...
        self.logger.info( 'copy source codes ...' )
        # add indent:
        self.LoggerSet( indent='  ' )

        # target directory:
        prefix = self.GetSetting( 'copy.prefix' )
        # info ...
        self.logger.info( 'destination prefix : %s ' % prefix )

        # remove existing build ?
        remove_existing_build = self.GetSetting( 'copy.new', 'bool' )
        # over-rule?
        if new is not None : remove_existing_build = new

        # extensions defined ?
        exts = self.GetSetting('copy.prefix.extensions').split()
        # defined ?
        if len(exts) > 0 :
            # start with original prefix:
            prefix_ext = prefix
            # extend:
            for ext in exts :
                prefix_ext = prefix_ext+'_'+ext
            #endfor
            # info ...
            self.logger.info( 'destination prefix extended : %s ' % prefix_ext )
            # create extended prefix directory if necessary:
            utopya_tools.CreateDirs( prefix_ext, renew=remove_existing_build, logger=self.logger )
            # test if link already exists; use 'lexists' to trap broken link:
            if os.path.lexists(prefix) :
                # is a link ?
                if os.path.islink(prefix) :
                    # remove current link:
                    os.remove( prefix )
                else :
                    self.logger.error( 'could not replace "%s" by a symbolic link; remove first' % prefix )
                    raise Exception
                #endif
            #endif
            # create link relative to current directory:
            os.symlink( os.path.basename(prefix_ext), prefix )
        else :
            # create prefix directory if necessary:
            utopya_tools.CreateDirs( prefix, renew=remove_existing_build, logger=self.logger )
        #endif

        # get list of sub directories to be scanned:
        subdirs = self.GetSetting('copy.subdirs').split()
        # eventually loop over sub directories:
        for subdir in subdirs :
            # full path:
            sdir = os.path.join(prefix,subdir)
            # create if necessary:
            utopya_tools.CreateDirs( sdir, logger=self.logger )
        #endfor

        # directories to be included:
        incdirs = self.GetSetting('copy.incdirs').split()

        # files to be skipped:
        skip_patterns = self.GetSetting('copy.skip').split()

        # format for listing files:
        fmt = '        %-100s     %-30s [%-8s]'

        # loop over source directories:
        sourcedirs = self.GetSetting('copy.dirs')
        if len(sourcedirs) > 0 :

            # info ...
            self.logger.info( 'copy files from source directories...' )
            # indent:
            self.LoggerSet( indent='  ' )

            # loop over source directories:
            for sourcedir in sourcedirs.split() :

                # set flag:
                found_some_files = False
                # info ...
                self.logger.info( 'scanning %s ...' % sourcedir)
                # should be a directory ...
                if not os.path.isdir(sourcedir) :
                    self.logger.error( 'specified source dir is not an existing directory : %s' % sourcedir )
                    raise Exception
                #endif
                # empty ? then add something for current directory:
                if len(subdirs) == 0 : subdirs = [os.path.curdir]
                # loop over sub directories:
                for subdir in  subdirs :
                    # add sub directory:
                    sourcepath = os.path.join(sourcedir,subdir)
                    # might be missing ...
                    if os.path.isdir(sourcepath) :
                        self.logger.debug( 'scanning %s ...' % sourcepath )
                    else :
                        self.logger.debug( 'missing  %s ; continue ...' % sourcepath )
                        continue
                    #endif

                    # copy all files from source path:
                    found = self._CopyFiles( sourcepath, prefix, subdir, incdirs=incdirs )
                    # update flag:
                    found_some_files = found_some_files or found

                #endfor   # subdirectories

                # check ...
                if not found_some_files :
                    self.logger.error('No source files found in standard subdirs "%s" of "%s" .  Mistake in "source.dirs" setting?' % (subdirs,sourcedir) )
                    raise Exception
                #endif

            #endfor   # source directories

            # reset indent:
            self.LoggerReset()

        else :
            # info ...
            self.logger.info( 'no source directories specified ...' )
        #endif

        # optional destination of evaluated rcfile:
        rcwrite = self.GetSetting( 'copy.rcwrite', default='None' )
        # defined?
        if rcwrite != 'None' :
            # write:
            self.WriteSettings( rcwrite )
        #endif

        # remove indent:
        self.LoggerReset()

    #enddef __init__

    # *

    def _CopyFiles( self, sourcepath, prefix, subdir, incdirs=[] ) :

        """
        Copy all files from source path to destination.
        """

        # modules:
        import os
        import shutil
        import fnmatch
        import filecmp

        # files to be skipped:
        skip_patterns = self.GetSetting('copy.skip').split()

        # remove extensions?
        remove__ext = self.GetSetting('copy.remove__ext',totype='bool',default=False)

        # format for listing files:
        fmt = '        %-100s     %-30s [%-8s]'

        # init flag:
        found_some_files = False

        # list all files in this directory, might be empty:
        sourcefiles = os.listdir(sourcepath)
        # loop over source files:
        for sfile in sourcefiles :
            # full filename:
            sourcefile = os.path.join(sourcepath,sfile)

            # directory?
            if os.path.isdir(sourcefile) :
                # to be copied?
                if sfile in incdirs :
                    # recursive call:
                    self._CopyFiles( os.path.join(sourcepath,sfile),
                                     prefix, os.path.join(subdir,sfile), incdirs=[] )
                #endif
                # skip, remainder is for files only:
                continue
            #endif

            # skip if matches pattern; init flag:
            skipit = False
            # loop over patterns to be skipped:
            for pattern in skip_patterns :
                # match?
                if fnmatch.fnmatch( sfile, pattern ) :
                    # info:
                    self.logger.debug( fmt % (sourcefile,'','skip') )
                    # set flag:
                    skipit = True
                    # leave loop:
                    break
                #endif
            #endfor # patterns
            # file skipped ? next:
            if skipit : continue

            # full target dir:
            targetdir = os.path.join( prefix, subdir )
            # create?
            if not os.path.isdir(targetdir) : os.makedirs( targetdir )

            # full target file:
            targetfile = os.path.join( targetdir, sfile )

            # remove '__xxx' extensions?
            if remove__ext :
                # has extension?
                if '__' in sfile :
                    # split:
                    bname,ext = os.path.splitext(sfile)
                    # new base:
                    bname = bname.split('__',-1)[0]
                    # new target name:
                    targetfile = os.path.join( targetdir, bname+ext )
                #endif
            #endif

            # already present ?
            if os.path.exists(targetfile) :
                # same content ?
                if filecmp.cmp(sourcefile,targetfile) :
                    stat = '...'
                else:
                    stat = 'differs'
                #endif
            else:
                stat = 'new'
            #endif
            # info ...
            self.logger.debug( fmt % (sourcefile,os.path.join(subdir,sfile),stat) )
            # copy source to target, preserve times etc:
            shutil.copy2( sourcefile, targetfile )
            # reset flag:
            found_some_files = True
        #endfor  # source files

        # ok
        return found_some_files

    #enddef _CopyFiles

#endclass UtopyaCopy


#-------------------------------------------------
# configure
#-------------------------------------------------


class UtopyaConfigure( utopya_rc.UtopyaRc ) :

    """
    Base class to perform tasks that configure a source code,
    which is mainly the creation of a suitable Makefile.
    The tasks are partly performed by the following methods:

    * :py:meth:`ConfigureCompilerFlags`
    * :py:meth:`ConfigureMacros`
    * :py:meth:`ConfigureChecks`
    * :py:meth:`ConfigureLibs`
    * :py:meth:`ConfigureMakeDependencies`

    Derived classes might re-define some of these if the standard versions
    are not sufficient.
    The following overview shows the rcfile settings needed for the configuration,
    and the methods that perform parts of the taks.

    * First specify the location of the source directory;
      typically this is the place to which :py:class:`UtopyaCopy`
      has copied the source files::

        ! where to configure ?
        [<rcbase>.]configure.dir            :  /work/model/src

    * Specify the name of the class that should provide compiler settings such
      as the compiler executable name(s) and the compilation flags.
      This class should be derived from UTOPyA's :py:class:`UtopyaCompilerSettings` class,
      for example the :py:class:`UtopyaCompilerSettingsFortran` class.
      Also specify the arguments for initialization; for the base class this
      is a rcfile name and keyword prefix::

        ! class for compiler settings, and arguments for initialization:
        [<rcbase>.]configure.compiler.class     :  utopya.UtopyaCompilerSettingsFortran
        [<rcbase>.]configure.compiler.args      :  rcfile='configure.compiler.rc', rcbase=''

      An instance of the specified class is created, and the :py:meth:`UtopyaCompilerSettingsFortran.GetCompilers`
      method of the instance is called to obtain the compiler(s) and linker names.
      The instance is also passed as argument to the :py:meth:`ConfigureCompilerFlags` methode.

    * The following flags specify if a source should be configured
      for running in parallel in an MPI or OpenMP environment.
      Their value is for example passed as arguments to the
      :py:meth:`UtopyaCompilerSettingsFortran.GetCompilers` method of (class derived from)
      the :py:class:`UtopyaCompilerSettings` class
      and to the :py:meth:`ConfigureCompilerFlags` method::

        ! enable parallel computing?
        [<rcbase>.]configure.enable.mpi     :  False
        [<rcbase>.]configure.enable.openmp  :  False

    * See the description of the :py:meth:`ConfigureCompilerFlags` method for
      the setings needed to select the appropriate compiler flags.

    * See the :py:meth:`ConfigureMacros` method for the settings needed
      for definition of preprocessing macro's.

    * If the macro definitions should be included in the compiler flags, e.g.::

        f90 .. -Dwith_test -D__VALUE__=12 ..

      then enable the following flag in the rcfile::

        [<rcbase>.]configure.macro.define_flags  :  True

      A line with the compiler flags is obtained from a call to the
      :py:meth:`UtopyaCompilerSettings.GetMacroDefinitionFlags` method of the
      :py:class:`UtopyaCompilerSettings` class (or one of its derivatives),
      which takes the the list of macro's
      returned by the :py:meth:`ConfigureMacros` method as argument.
      
    * Optionally specify a list of depricated or un-used files that will
      be removed from the build directory::
      
        [<rcbase>.]configure.remove  :  file1.F90 file2.F90 ...
        
      If files should only be removed if a certain macro is defined, use::
      
        [<rcbase>.]configure.remove.ifdef.with_test  :  file3.F90 file4.F90 ...
        
      If files should be removed if a certain macro is not defined, use:

        [<rcbase>.]configure.remove.ifndef.with_test :  test1.F90 test2.F90 ...
        
    * See the :py:meth:`ConfigureChecks` method for the settings needed
      to perform some quick checks on the source code.

    * The :py:meth:`ConfigureLibs` method is called to obtain the
      compiler and linker flags related to use of external libraries.

    The methode will write two textfiles that could be included in the final Makefile.
    It is up to the user to write a suitable Makefile that actually does this.
    The name of the included files should be spefied in the settings:

    * The compiler flags will be written to a textfile with a name specified by::

        configure.flags.includefile   :  Makefile_compiler_and_flags

      The written file might look like::

        #
        # include file with compiler names and flags for Makefile.
        #

        # compiler and linker:
        FC = f90
        F77 = f77
        LINKER = f90

        # compile flags:
        FFLAGS = --implicit-none -O2 -g -Dwith_test -I/opt/netcdf4/include

        # linker flags:
        LDFLAGS = -L/opt/netcdf4/lib -lnetcdf -L/usr/lib -ljasper

    * Optionally, explicit compilation rules are written to a file specified by::

        configure.rules.includefile   :  Makefile_rules

      If this file is defined, it will be created and filled
      with compilation rule(s) that differ from the default rule(s).
      This could be used to specify that for some files implicit typing is allowed
      while the default is to not allow that, or to disable optimizations
      if that causes failues for some files.
      The created include file might look like::

        #
        # include file with explicit compile rules for Makefile.
        #

        oldstuff.o:
          $(FC) -c -o $@ $(FFLAGS) --implicit-yes $<

        olderstuff.o:
          $(FC) -c -o $@ $(FFLAGS) --implicit-yes -O0 $<

      In this example, the '--implicit-yes' and '-O0' flags are inserted after
      the expansion of the default '$(FFLAGS)', which usually means that similar
      setting in '$(FFLAGS)' are over-rulled.

      To specify which objects require an explicit rule, and which flags
      should be added, use the following configuration::

        ! specify a list with objects and additional flag groups:
        !
        !    objectfile1 : flaggroup1 ;  objectfile2 : flaggroup2a flaggroup2b ; ...
        !
        configure.rules.explicit : \\
             oldstuff.o    :  implicit \\
             olderstuff.o  :  implict optim-none

      See the description of the :py:meth:`ConfigureCompilerFlags` method for
      how the flag group keywords 'implicit' etc are expanded to actual flags
      '--implicit-yes' etc.

    * Makefile dependencies might be created automatically if the following
      flag is enabled::

        ! create depencencies (True|False) ?
        configure.makedep        :  True

      If enabled, the method :py:meth:`ConfigureMakeDependencies` will be called.
      Depending on the settings and the sources found, this might for example
      write a file 'Makefile_deps' that will be included in the main 'Makefile'.

    The user defined Makefile could now look like::

      #
      # Makefile
      #

      # compiler name and flags:
      include "Makefile_compiler_and_flags"

      # default rule for how to form object files from f90 source:
      %.o: %.f90
        $(FC) -c -o $@ $(FFLAGS) $<

      # explicit rules for some files:
      include "Makefile_rules"

      # dependencies:
      #   myprog.o: mytool.o
      # here include dependencies from a file (created by configure?)
      include "Makefile_deps"

      # objects:
      OBJS = myprog.o mytool.o

      # target executable:
      myprog.x: myprog.o
        $(LINKER)  -o $@ $(OBJS) $(LDFLAGS)

    """

    def __init__( self, rcfile, rcbase='', env={} ) :

        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )

        # modules:
        import os
        import glob

        # tools:
        import utopya_tools

        # info ...
        self.logger.info( 'configure source code ...' )

        # work directory:
        configure_dir = self.GetSetting( 'configure.dir' )
        # info ...
        self.logger.info( '  change to %s ...' % configure_dir )
        # store old:
        owd = os.getcwd()
        # change:
        os.chdir( configure_dir )

        # environment defined?
        if len(env) > 0 :
            # info ...
            self.logger.info( '  setup environment ..' )
            # loop:
            for key in env.keys() :
                # info ...
                self.logger.info( '    %s = %s' % (key,env[key]) )
                # export to system envionment:
                os.environ[key] = env[key]
            #endfor
        #endif

        #
        # * user defined methods
        #

        # list of user methods:
        tasks = self.GetSetting( 'configure.tasks', default='' ).split()
        # loop:
        for task in tasks :
            # info ...
            self.logger.info( '  task "%s" ...' % task )
            # import class for compiler settings:
            cls = self.ImportClass( 'configure.task.%s.class' % task )
            # arguments for initialization:
            args = self.GetSetting( 'configure.task.%s.args' % task )
            # replace templates:
            args = args.replace('%{rcfile}',self.rcfile)
            # create object:
            tsk = eval( 'cls( %s )' % args )
        #endfor # tasks


        #
        # * compiler name(s) and flags
        #

        # import class for compiler settings:
        cls = self.ImportClass( 'configure.compiler.class' )
        # arguments for initialization:
        args = self.GetSetting( 'configure.compiler.args' )
        # replace templates:
        args = args.replace('%{rcfile}',self.rcfile)
        # create object:
        compiler = eval( 'cls( %s )' % args )

        # configure for parallel computing?
        mpi    = self.GetSetting( 'configure.enable.mpi'   , 'bool' )
        openmp = self.GetSetting( 'configure.enable.openmp', 'bool' )

        # call user script to set compilers and linker;
        # this returns a dictionairy with compiler names that will
        # be used in the Makefile and their value:
        #   dict( FC='f90', F77='f77', LINKER='f90' )
        compilers = compiler.GetCompilers( mpi=mpi, openmp=openmp )

        # obtain flags:
        fflags,ldflags = self.ConfigureCompilerFlags( compiler, mpi=mpi, openmp=openmp )

        #
        # * macro's
        #

        # obtain list with macro's to be defined,
        # and list with all supported macros':
        macros_defined,macros_supported = self.ConfigureMacros()

        # add macro definition flags ?
        define_flags = self.GetSetting( 'configure.macro.define_flags', 'bool' )
        # apply?
        if define_flags :
            # create list of macro definitions as command line arguments,
            # e.g. -Dwith_this_flag etc:
            Dflags = compiler.GetMacroDefinitionFlags( macros_defined )
            # add definitions to flags:
            fflags = fflags.strip()+' '+Dflags
        #endif


        #
        # * checks
        #

        # apply soure code checks:
        self.ConfigureChecks()


        #
        # * libraries
        #

        # obtain f and ld flags related to external libraries:
        fflags_libs,ldflags_libs = self.ConfigureLibs( macros_defined )

        # add:
        fflags  =  fflags+' '+fflags_libs
        ldflags = ldflags+' '+ldflags_libs


        #
        # * write compiler flags to include file incorporated into Makefile
        #

        # name of include file with dependencies, to be written by this script:
        makefile_flags = self.GetSetting('configure.flags.includefile')
        # info ...
        self.logger.info( '  write %s (compiler and flags) ...' % makefile_flags )
        self.logger.info( '    compiler fflags : '+fflags  )
        self.logger.info( '            ldflags : '+ldflags )
        # fill content; each line should end with the newline expression '\n' :
        src = []
        src.append( '#\n' )
        src.append( '# include file with compiler flags for Makefile.\n' )
        src.append( '#\n' )
        src.append( '\n' )
        src.append( '# compiler and linker:\n' )
        for key in compilers.keys() :
            src.append( '%s = %s\n' % (key,compilers[key]) )
        #endfor
        src.append( '\n' )
        src.append( '# compile flags:\n' )
        src.append( 'FFLAGS = %s\n' % fflags )
        src.append( '\n' )
        src.append( '# linker flags:\n' )
        src.append( 'LDFLAGS = %s\n' % ldflags )
        src.append( '\n' )

        # write:
        utopya_tools.WriteTextFile( makefile_flags, src )


        #
        # * write compile rules to makefile include file:
        #

        # name of include file with compiler rules, to be written by this script:
        makefile_rules = self.GetSetting('configure.rules.includefile',default='None')
        # defined?
        if makefile_rules != 'None' :
            # info ...
            self.logger.info( '  write %s (compile rules) ...' % makefile_rules )
            # init content; each line should end with the newline expression '\n' :
            src = []
            # header:
            src.append( '#\n' )
            src.append( '# include file with explicit compile rules for Makefile.\n' )
            src.append( '#\n' )
            src.append( '\n' )
            # explicit rules for some files
            # in rcfile, a list is specified with objects and additional flag groups:
            #      objectfile1  :  flaggroup1 ; objectfile2 : flaggroup2 ; ...
            # read list:
            rlist = self.GetSetting( 'configure.rules.explicit' )
            # split in ojbect/flaggroup parts:
            objgrps = rlist.split(';')
            # loop:
            for objgrp in objgrps :
                # cleanup:
                objgrp = objgrp.strip()
                # empty ?
                if len(objgrp) == 0 : continue
                # split into object and groups:
                obj,grps = objgrp.split(':')
                # cleanup:
                obj = obj.strip()
                # expand groups into flags:
                fflags_expl = ''
                for grp in grps.split() :
                    ff,ldf = compiler.GetFlags( grp )
                    fflags_expl = fflags_expl+' '+ff.strip()
                #endfor
                fflags_expl = fflags_expl.strip()
                # add rule:
                src.append( '%s:\n' % obj )
                src.append( '\t$(FC) -c -o $@ $(FFLAGS) %s $<\n' % fflags_expl )
                src.append( '\t@echo " "\n' )
                src.append( '\n' )
            #endfor

            # write:
            utopya_tools.WriteTextFile( makefile_rules, src )

        #endif # include file for rules defined


        #
        # * cleanup
        #

        # remove files ?
        self.logger.info( '  remove files ...' )
        # obtain file list, might not be present:
        sfiles = self.GetSetting('configure.remove').split()
        # loop:
        for sfile in sfiles :
            # info ...
            self.logger.debug( +'    %s ...' % sfile )
            # remove file if present:
            if os.path.exists(sfile) : os.remove( sfile )
        #endfor

        # clean for defined macro's (read list of pattern)
        self.logger.info( '  remove files for defined macro\'s ...' )
        for macr in macros_defined :
            # try to read list with files to be removed if macro is defined:
            patterns = self.GetSetting( 'configure.remove.ifdef.%s' % macr, default='' )
            # line with list of patterns was found ?
            if len(patterns) > 0:
                # loop over patterns:
                for pattern in patterns.split():
                    # loop over files matching the pattern:
                    for sfile in glob.glob(pattern) :
                        # info ...
                        self.logger.debug( '    %s ...' % sfile )
                        # remove:
                        os.remove(sfile)
                    #endfor # matching files
                #endfor # patterns
            #endif  # found ?
        #endfor   # loop over defined macros

        # clean for undefined macro's (read list of pattern)
        self.logger.info( '  remove files for undefined macro\'s ...' )
        for macr in macros_supported :
            # defined ? then next:
            if macr in macros_defined : continue
            # try to read list with files to be removed if macro is NOT defined:
            patterns = self.GetSetting( 'configure.remove.ifndef.%s' % macr, default='' )
            # line with list of patterns was found ?
            if len(patterns) > 0:
                # testing ...
                self.logger.info( '    macro "%s" undefined ...' % macr )
                # loop over patterns:
                for pattern in patterns.split():
                    # loop over files matching the pattern:
                    for sfile in glob.glob(pattern) :
                        # info ...
                        self.logger.info( '      remove %s ...' % sfile )
                        # remove:
                        os.remove( sfile )
                    #endfor # matching files
                #endfor # patterns
            #endif  # found ?
        #endfor   # loop over supported macros


        #
        # * deps
        #

        # read main flag:
        enabled = self.GetSetting('configure.makedep','bool',default=False)
        # enabled ?
        if enabled :
            # create Makefile dependencies based on settings and sources:
            self.ConfigureMakeDependencies( define_flags=define_flags )
        #endif

        #
        # *
        #

        # info:
        self.logger.debug( '  change back to %s ...' % owd )
        # back:
        os.chdir( owd )

    #enddef __init__

    # *

    def ConfigureCompilerFlags( self, compiler, mpi=False, openmp=False ) :

        """
        Return compiler flags to be used. Arguments:

        * compiler : instance of the :py:class:`UtopyaCompilerSettings` class

        Optional arguments:

        * mpi : bool flags to enable compilation with MPI compilers;.
        * openmp : bool flags to enable compilation with OpenMP support.

        Return values:

        * fflags : str object with Fortran compiler flags, for example '-O2 -g' ;
        * ldflags : str object with linker flags

        In the settings, specify a space seperated list of what is called 'compiler flags groups'::

            [<rcbase>.]configure.flags.groups     :  optim-none check-all

        Each group name is passed as argument to the :py:meth:`UtopyaCompilerSettings.GetFlags` method,
        which should then return the the actual flags to be passed to the compiler
        as defined in the compiler specific rcfile settings.

        The following group names are added automatically:

        * 'default' : to insert default flags, for example to disable implicit typing;
        * 'mpi' if MPI is enabled;
        * 'openmp' if OpenMP is enabled.

        """

        # info...
        self.logger.info( '  collect compiler flags ...' )

        # start without any flags:
        fflags  = ''
        ldflags = ''

        # list with names of compiler flag groups to be used:
        flaggroups = ['default']
        if mpi    : flaggroups.append( 'mpi' )
        if openmp : flaggroups.append( 'openmp' )
        flaggroups.extend( self.GetSetting( 'configure.flags.groups' ).split() )
        # loop over groups:
        for flaggroup in flaggroups :
            # info ...
            self.logger.info( '    add flag group "%s" ...' % flaggroup )
            # get flags:
            ff,ldf = compiler.GetFlags( flaggroup )
            # add to lists:
            fflags  =  fflags.strip()+' '+ff
            ldflags = ldflags.strip()+' '+ldf
        #endfor
        # info ...
        self.logger.info( '    collected:' )
        self.logger.info( '      fflags       : %s' % fflags )
        self.logger.info( '      ldflags      : %s' % ldflags )

        # ok
        return fflags,ldflags

    #enddef ConfigureCompilerFlags

    # *

    def ConfigureMacros( self ) :

        """
        The source files might contain pre-processing statements such as::

          #ifdef with_test
          print *, 'this is test output'
          #endif

          #if __VERSION__ == 12
          call Work_v12()
          #else
          call Work()
          #endif

        In this example, the print statement is only included in the actual code
        if the pre-processing macro 'with_test' is defined, and which version of
        the work routine is called depends on the value of the '__VERSION__' macro
        (if defined).
        Pre-processing macros could be defined in the following ways:

        * The compiler might automatically define certain macro's.
          Typically some macro's are defined to hold the compiler version,
          and if OpenMP is enabled then the macro '_OPENMP' is defined.

        * A macro definition flag could be added to the compiler flags::

            f90 .. -Dwith_test -D__VALUE__=12 ..

        * A macro definition statement could be included in the source file::

            #define with_test
            #define __VALUE__ 12

          Since preprocessing macro's are often used in more than one source file,
          it is common practice to collect all macro definitions in a header file
          (in Fortran called an 'include' file), and incorporate this in each source file::

            #include macros.inc

        The method using header (include) files is prefered since it ensures that a source code
        could be compiled later on again using exactly the same macro definitions.
        UTOPyA therefore provides a mechanism to create macro definition include files.

        First specify a list of macro group names, for example one name
        for each sub-program in a source::

          [<rcbase>.]configure.macro.groups       :  model optimizer

        For each of the groups, specify the following settings:

        * A list with all macro's in the group.
        * Optionally the name of a header file; if provided,
          this file is created in the source directory and filled
          with appropriate '#define' commands. The corresponding
          '#include' statement should be added to each source file
          that uses the macro's in the group.

        An example for the settings for the above macro groups ::

          [<rcbase>.]configure.macro.model.all         :  with_test without_error
          [<rcbase>.]configure.macro.model.hfile       :  model.inc

          [<rcbase>.]configure.macro.opimizer.all      :  DEBUG __VERSION__
          [<rcbase>.]configure.macro.optimizer.hfile   :  optim.inc

        Finally specify a list of the macro's that actually need to be
        defined in a certain application. The macro's defined should appear
        in at least one of the lists with all macro's for a group::

          [<rcbase>.]configure.macro.define       :  with_test __VERSION__=12

        With the above settings, this method will:

        * check if the requested macro's are supported;
        * write header files if necessary.

        Return values:

        * list with character strings holding the defined macro's (and values)::

          [ 'with_test', '__VERSION__=12' ]
          
        * list with character strings holding the supported macro's:

          [ 'with_test', 'with_fix', '__VERSION__' ]

        """

        # tools:
        import utopya_tools

        # info...
        self.logger.info( '  macro definitions ...' )

        # defined macro's:
        macros_def   = self.GetSetting( 'configure.macro.define' ).split()

        # create a 'clean' list of defined macro's without double definitions:
        macros_defined   = []
        for m in macros_def :
            if m not in macros_defined :
                macros_defined = macros_defined + macros_def
            #endfor
        #endfor # defined macro's

        # info ..
        if len(macros_defined) > 0 :
            self.logger.info( '    defined macros:' )
            for m in macros_defined :
              self.logger.info( '      %s' % m )
            #endfor
        else :
            self.logger.info( '    no macros defined ...' )
        #endif

        # macro groups:
        groups =  self.GetSetting( 'configure.macro.groups' ).split()

        ## apply user changes to default list:
        #for group in groups :
        #    # add (or remove!) macro defintions by user script:
        #    macros_def = pycus.Build_Define( rcf, group, macros_def )
        ##endfor

        # initialize a list to store all supported macro's without double definitions:
        macros_supported = []

        # loop over groups:
        for group in groups :
            # info...
            self.logger.info( '    process macro group "%s" ...' % group )
            # start of the rc keys:
            keybase = 'configure.macro.'+group
            # values for this group:
            macros_all   = self.GetSetting( keybase+'.all'   ).split()
            macros_hfile = self.GetSetting( keybase+'.hfile' )
            # write header file ?
            if len(macros_hfile) > 0 :
                # info ...
                self.logger.info( '      write header file "%s" ...' % macros_hfile )
                # fill text for include file in list of strings,
                # each line should end with the newline expression '\n' :
                src = []
                src.append( '!\n' )
                src.append( '! Include file with macro definitions.\n' )
                src.append( '!\n' )
                # loop over macro's to be defined:
                for mdef in macros_def :
                    # vallue assigned ?
                    if '=' in mdef :
                        # split in name and value:
                        mname,mval = mdef.split('=')
                        # in this group ?
                        if mname in macros_all :
                            # add line to file:
                            src.append( '#define %s %s\n' % (mname,mval) )
                        #endif
                    else :
                        # in this group ?
                        if mdef in macros_all :
                            # add line to file:
                            src.append( '#define %s\n' % mdef )
                        #endif
                    #endif
                #endfor
                # write the source, replace existing file only if it was different:
                utopya_tools.UpdateTextFile( macros_hfile, src )
            #endif
            # extend list:
            for m in macros_all :
                if m not in macros_supported : macros_supported = macros_supported + macros_all
            #endfor
        #endif

        # check for macro's that are not supported yet:
        any_error = False
        for macr_def in macros_defined :
            # split if necessary:
            if '=' in macr_def :
                macr = macr_def.split('=',1)[0]
            else :
                macr = macr_def
            #endif
            # not supported ?
            if macr not in macros_supported :
                # any unsupported macro's found yet ?
                if not any_error :
                    # initial error message:
                    self.logger.error( "one or more macro's have been defined that are not listed" )
                    self.logger.error( "in any 'configure.macro.*.all' lists:" )
                #endif
                # display problematic macro:
                self.logger.error( "  %s" % macr )
                # reset flag:
                any_error = True
            #endif
        #endfor
        # any error ? then leave:
        if any_error : raise Exception

        # ok
        return macros_defined,macros_supported

    #enddef ConfigureMacros

    # *

    def ConfigureChecks( self ) :

        """
        Apply checks on source code files.

        Provide a space seperated list with keywords
        for the checks to be applied::

            ! list with keywords for checks to be applied:
            [<rcbase>.]configure.checks      :  go_inc unknown_macro strange_file

        For each keyword, provide settings with:

        * a short message to describe the test;
        * the filename pattern to select the files on
          which the test should be applied;
        * a filename pattern to exclude some files;
        * a python boolean expression applied to each line of the file;
          the line itself is stored in a string variable named 'line';
        * a help text to be displayed in case a warning is issued;
          include \n for newlines.

        Example for check if all ``go*.F90`` files do include ``go.inc``::

            [<rcbase>.]configure.check.go_inc.msg       :  Test on GO files that do no not include "go.inc" ...
            [<rcbase>.]configure.check.go_inc.files     :  go*.F90
            [<rcbase>.]configure.check.go_inc.skip      :
            [<rcbase>.]configure.check.go_inc.test      :  not line.startswith('#include "go.inc"')
            [<rcbase>.]configure.check.go_inc.help      :  \n\
                  All GO files should include "go.inc" in the header.\n\
                  This include file is filled by the scripting with '#define' pre-processing macros.

        If for one of the lines in a source file the test evaluates to 'True',
        a warning is issued.
        If the test expression starts with ``not``, a warning is issued if none
        of the lines evaluates to ``True`` for the test after the ``not``.

        The following check will warn for unexpected file names::
        
            ! source directories should not be polluted with strang files
            [<rcbase>.]configure.check.strange_file.msg       :  Test on strange file names ...
            [<rcbase>.]configure.check.strange_file.files     :  *
            [<rcbase>.]configure.check.strange_file.skip      :  *.inc *.f *.F *.f90 *.F90 Makefile* *.o *.mod *.x
            [<rcbase>.]configure.check.strange_file.test      :  True
            [<rcbase>.]configure.check.strange_file.help      :  \n\
                  Source directories should not be polutted with strange files\
                  such as backup files created by editors.
                  
        Some checks are too complicated to be defined with the above settings.
        The following special test is used to check the macro definitions::

            ! check on unknown macros used in '#if[n]def' lines,
            ! but not listed in the configuration settings:
            [<rcbase>.]configure.check.unknown_macro.msg :  Test for unsupported macros ...

        Use the following flag to raise an error instead of a warning::

            [<rcbase>.]configure.check.error   :  True

        """
        
        # modules:
        import os
        import fnmatch

        # info...
        self.logger.info( '  source code checks ...' )

        # keywords for checks to be performed:
        checknames = self.GetSetting( 'configure.checks', default='' )
        # none defined ?
        if len(checknames) == 0 : return

        # error or just warnings ?
        with_error = self.GetSetting( 'configure.checks.error', 'bool' )

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
            test_msg   = self.GetSetting( 'configure.check.%s.msg'   % checkname )
            test_files = self.GetSetting( 'configure.check.%s.files' % checkname ).split()
            test_skip  = self.GetSetting( 'configure.check.%s.skip'  % checkname ).split()
            test_line  = self.GetSetting( 'configure.check.%s.test'  % checkname, 'str' )
            test_help  = self.GetSetting( 'configure.check.%s.help'  % checkname )

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
                    if not matching_files : self.logger.warning( '    %s : [found]' % test_msg )
                    self.logger.warning( '      %s' % srcfile )
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
                    for helpline in test_help.split('\\n') : self.logger.warning(helpline)
                #endif
            else :
                # no warnings for this test ...
                self.logger.info( '    %s [none ]' % test_msg )
            #endif

        #endfor  # checks

        # check for unknown macro's ?
        checkname = 'unknown_macro'
        if checkname in checknames :

            # settings:
            test_msg   = self.GetSetting( 'configure.check.%s.msg' % checkname )

            # names of macro groups:
            macgroups = self.GetSetting( 'configure.macro.groups' ).split()
            # collect all supported macro's:
            macall = []
            for macgroup in macgroups :
                macs = self.GetSetting( 'configure.macro.%s.all' % macgroup ).split()
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
                                self.logger.info( '    %s' % test_msg )
                                logged_msg = True
                            #endif
                            # intro if necessary:
                            if not logged_srcfile :
                                self.logger.error( '      unsupported macro(s) in %s :' % srcfile )
                                logged_srcfile = True
                            #endif
                            # line number and content:
                            self.logger.error( '      %6i : %s' % (iline,line) )
                            # set flag:
                            any_warning = True
                        #endif   # unsuported macro
                    #endif  # line with macro test
                #endfor   # lines

            #endfor   # source files

            # jippy ...
            if not logged_msg :
                # no warnings for this test ...
                self.logger.info( '      %s [none ]' % test_msg )
            #endif

        #endif   # test on unsupported macro's

        # break ?
        if any_warning and with_error :
            self.logger.error( 'some source code checks failed; break' )
            self.logger.error( '(set "configure.checks.error : False" in the expert.rc to avoid this error)' )
            raise Exception
        #endif

    #enddef ConfigureChecks

    # *

    def ConfigureLibs( self, macros=[] ) :

        """
        Return compiler flags related to external libraries.

        Optional arguments:

        * macros : list object holding character strings with defined macros
          as returned by the :py:meth:`ConfigureMacros` methode.

        Return values:

        * fflags : str object with Fortran compiler flags, for example::

            '-I/opt/netcdf4/include'

        * ldflags : str object with linker flags, for example::

            '-L/opt/netcdf4/lib -lnetcdf -L/usr/lib -ljasper'

        In the settings, first specify the name of the class that should provide
        library settings such as the compile and link flags.
        This class should be derived from UTOPyA's :py:class:`UtopyaLibSettings` class.
        Also specify the arguments for initialization, for the base class this
        is an rcfile name and keyword prefix.
        See the documentation of the class for the settings that need to be present
        in this rcfile; these are usually machine dependend::

          ! class for library settings, and arguments for initialization:
          [<rcbase>.]configure.libs.class     :  utopya.UtopyaLibSettings
          [<rcbase>.]configure.libs.args      :  rcfile='machine.rc', rcbase=''

        An instance of the specified class is created, and the :py:meth:`UtopyaLibSettings.GetFlags`
        method of the instance is called to obtain the compile and link flags.

        Then specify a list with names (keywords) for all supported external libraries::

            configure.libs.all    :  netcdf4 jasper hdf5 sz z

        The order of the names should be from 'high' level to 'low' level (a high level library
        might depend on a low level library). This order is important
        for linking, and should therefore be specified.

        Finaly specify a list with the names of the libraries that should actually be linked;
        the order is not important here::

            [<rcbase>.]configure.libs      :  jasper netcdf4

        If the optional 'macros' list with macro definitions is passed as argument,
        the settings are scanned for a specification of libraries that
        should be linked if a certain macro is defined::

            [<rcbase>.]configure.libs.ifdef.with_lapack   :  lapack blas
            [<rcbase>.]configure.libs.ifdef.with_mkl      :  mkl

        This feature could be used for codes that depend on libraries that
        are not always available in the same way on different machines.
        For example, when using Intel compilers, optimizer versions of BLAS
        and LAPACK are available through the Math Kernel Library (MKL).
        """

        # ~ all libraries

        # get list of all supported library names, in the right order for linking:
        libnames_all = self.GetSetting( 'configure.libs.all' ).split()

        # ~ default libraries

        # get list of default library names to be linked:
        libnames = self.GetSetting( 'configure.libs' )

        # info ...
        self.logger.debug( '  libraries to be used:' )
        if len(libnames) > 0 : self.logger.debug( '    %s    (default)' % libnames )

        # convert to list:
        libnames = libnames.split()

        # ~ macro depended libs

        # loop over defined macro's:
        for macro in macros :
            # read list of libraries that should be added if this macro is defined:
            libnames_ifdef = self.GetSetting( 'configure.libs.ifdef.%s' % macro, default='' )
            # defined ?
            if len(libnames_ifdef) > 0 :
                # add:
                libnames = libnames+libnames_ifdef.split()
                # info ...
                self.logger.debug( '    %s    (macro `%s` defined)' % (libnames_ifdef,macro) )
            #endif
        #enddef

        # ~ check ...

        # check if all libraries specified to be used are actually supported ...
        for libname in libnames :
            if libname not in libnames_all :
                self.logger.error( 'library name `%s` not in `configure.libs.all` list ...' % libname )
                raise Exception
            #endif
        #endfor

        # ~ flags

        # import class for libs settings:
        cls = self.ImportClass( 'configure.libs.class' )
        # arguments for initialization:
        args = self.GetSetting( 'configure.libs.args' )
        # replace templates:
        args = args.replace('%{rcfile}',self.rcfile)
        # create object:
        libset = eval( 'cls( %s )' % args )

        # info ...
        self.logger.debug( '  libraries linked (in this order):' )
        # now add compiler and linker flags ; initialize flag lines:
        fflags = ''
        ldflags = ''
        # loop over all supported libraries (this specfifies the linking order!)
        for libname in libnames_all :
            # not in use ? then skip:
            if libname not in libnames : continue
            # info ...
            self.logger.debug( '    %s' % libname )
            # obtain flags:
            ff,ldf = libset.GetFlags( libname )
            # add include, module, and link flags:
            fflags  =  fflags.strip()+' '+ff
            ldflags = ldflags.strip()+' '+ldf
        #endfor

        # ~

        # ok
        return fflags,ldflags

    #enddef ConfigureLibs

    # *

    def ConfigureMakeDependencies( self, define_flags=False, macros_defined=[] ) :

        """
        Create makefile dependencies based on settings and sources that are found.
        The dependencies will be written to file that should be included in the
        main `Makefile`.

        The following (optional) arguments could be set by the calling method:

        * Set the 'define_flags' bool to 'True' if source code contains pre-processing
          macro's that should be defined by compiler flags, e.g.::

              f90 .. -Dwith_test -D__VALUE__=12

          The depency generator needs to be aware of this to ensure that dependencies
          are defined for the final pre-processed source.

        * The 'macros_defined' list should contain the macro definitions that are needed
          when the 'define_flags' argument is 'True'.

        The following rcfile setting defines the name of the file that will contain the dependencies::

          ! include file to be written for Makefile:
          configure.makedep.includefile     :  Makefile_deps

        An ojbect derived rom the :py:class:`UtopyaDependencies` class
        is used to create the dependency lines.
        The name of the class and the initilization arguments are read
        from the settings::

          ! class for library settings, and arguments for initialization:
          [<rcbase>.]configure.makedep.class     :  utopya.UtopyaDependencies
          [<rcbase>.]configure.makedep.args      :  rcfile='build.rc', rcbase=''

        An instance of the specified class is created, and the :py:meth:`UtopyaDependencies.GetLines`
        method of the instance is called to obtain the dependency lines.
        The arguments for this method are read from the following settings::

          ! dependencies are created for files matching the pattern:
          configure.makedep.files           :  *.[Ff] *.[Ff]90

          ! define explicit dependencies for some objects,
          ! in case these are not recoqnized by 'makedep',
          ! for example dependencies on object files that are no modules;
          ! specifiy ';' seperated lines of target and object list seperated by ':' :
          !    target-object1 : object1 object2 ... ;
          !    target-object2 : object3 object4 ...
          configure.makedep.explicit        :

        After the dependency lines are obtained, the lines that define the
        linking of the executable(s) can be formed.
        The name(s) of the executables are read from the following setting::

          ! Which executables to be build ?
          ! Provide ';' seperated list of executable names,
          ! eventually followed by the basename of the source file
          ! with the main program if this is different from the
          ! basename of the executable. For example, the following
          ! tells to make executables 'test.x' and 'boe.x',
          ! where the first is to be compiled from a 'main.F90' or so,
          ! and the seccond is assumed to be compiled from 'boe.F90':
          configure.makedep.exec            : test.x main ; boe.x

        With this information, the make-includefile is written, which could look like::

            boe.o: boe.F90 tools.o
            driver.o: driver.F90 tools.o
            main.o: main.F90 driver.o tools.o
            tools.o: tools.F90

            test.x: main.o
                    $(LINKER) -o $@ main.o driver.o tools.o $(LDFLAGS)

            boe.x: boe.o
                    $(LINKER) -o $@ boe.o tools.o $(LDFLAGS)

        """

        # modules:
        import os

        # name of include file with dependencies to be written:
        makefile_deps = self.GetSetting('configure.makedep.includefile')
        # info ...
        self.logger.info( '  create %s ...' % makefile_deps )

        # import class for libs settings:
        cls = self.ImportClass( 'configure.makedep.class' )
        # arguments for initialization:
        args = self.GetSetting( 'configure.makedep.args' )
        # replace templates:
        args = args.replace('%{rcfile}',self.rcfile)
        # create object:
        makedep = eval( 'cls( %s )' % args )

        # file filter for dependencies:
        files = self.GetSetting('configure.makedep.files')

        # explicit dependencies:
        #    target-object1 : object1 object2 ... ;
        #    target-object2 : object3 object4 ...
        explicit = self.GetSetting( 'configure.makedep.explicit' )
        # convert into dictionairy:
        expl = {}
        targlists = explicit.split(';')
        for targlist in targlists :
            if len(targlist) > 0 :
                target,objs = targlist.split(':')
                target = target.strip()
                expl[target] = []
                for obj in objs.split() : expl[target].append( obj.strip() )
            #endif
        #endfor

        # create list of macro definitions as command line arguments,
        # e.g. -Dwith_this_flag etc; init empty:
        mk_defs = ''
        # add macro definitions ?
        if define_flags :
            # compiler depended prefix for macro definition:
            mk_D = '-D'
            # loop over macro's to be defined:
            for mdef in macros_defined :
                # add definition to command line argument list:
                mk_defs = mk_defs.strip()+(' %s%s' % (mk_D,mdef) )
            #endfor
        #endif

        # get dependency lines:
        deplines = makedep.GetLines( flags=mk_defs, files=files, expl=expl )

        # no link lines yet:
        linklines = ['']
        # names of target(mainfile)s to be build:
        items = self.GetSetting( 'configure.makedep.exec' )
        # loop over ';' seperated items:
        for item in items.split(';') :

            # split into list, might be single element:
            values = item.split()
            # check length to determine format:
            # ~ form is : exec
            if len(values) == 1 :
                # extract target:
                target = values[0]
                # main source file has same basename:
                mainfile,ext = os.path.splitext(target)
            # ~ form is : exec mainfile
            elif len(values) == 2 :
                # extract values:
                target,mainfile = values
            # ~ strange ...
            else :
                self.logger.error( 'items in exec list should have form :  exec [mainfile]' )
                self.logger.error( 'found item : %s' % item )
                raise Exception
            #endif

            # get list of objects for this target:
            objs = makedep.GetObjects( mainfile )
            # link lines:
            linklines.append( '%s: %s.o' % (target,mainfile) )
            linklines.append( '\t$(LINKER) -o $@ %s $(LDFLAGS)' % objs )
            linklines.append( '' )
        #endif

        # write result:
        f = open( makefile_deps, 'w' )
        for line in deplines+linklines : f.write(line+'\n')
        f.close()
        # add to log file:
        self.logger.debug( '' )
        self.logger.debug( '---[%s]------------------------------------------------------------' % makefile_deps )
        for line in deplines+linklines : self.logger.debug(line)
        self.logger.debug( '-------------------------------------------------------------------' )
        self.logger.debug( '' )

    #enddef ConfigureMakeDependencies

#endclass UtopyaConfigure


#-------------------------------------------------
# compiler settings
#-------------------------------------------------

class UtopyaCompilerSettings( utopya_rc.UtopyaRc ) :

    """
    Base class for objects that should return compiler specific settings
    defined in a rcfile.
    """

    def GetCompilers( self, mpi=False, openmp=False ) :

        """
        Returns a dictionairy with compiler names that will
        be used in the Makefile and their value::

            { 'FC' : 'f90', 'F77' : 'f77', 'LINKER' : 'f90' }

        The names of the compilers/linkers might be defined based
        on the ``mpi`` and/or ``openmp`` flags.

        This template returns some default settings.
        """

        # ok
        return { 'FC' : 'f90', 'F77' : 'f77', 'LINKER' : 'f90' }

    #enddef GetCompilers

    # *

    def GetFlags( self, name ) :

        """
        Return compiler flags for specified name.
        Empty for this base class.

        Returns two str objects:

        * fflags   : flags passed to fortran compiler
        * ldflags  : flags passed to linker
        """

        # read from settings:
        fflags  = ''
        ldflags = ''

        # ok:
        return fflags,ldflags

    #enddef GetFlags

    # *

    def GetMacroDefinitionFlags( self, macros ) :

        """
        Return a line with command line arguments to define
        preprocessing macro's. Such a line probably looks like::

          -Dwith_this_flag -Dwith_that_flag ...

        Pass a list with macro's to be defined as argument.

        The following setting should contain the macro definition flag
        of the compiler::

          [<rcbase>.]configure.compiler.defineflag    :  -D

        """

        # compiler depended prefix for macro definition:
        Dflag = '-D'

        # init empty:
        flags = ''
        # loop over macro's to be defined:
        for macro in macros :
            # add definition to command line argument list:
            flags = flags.strip()+(' %s%s' % (Dflag,macro) )
        #endfor

        # ok:
        return flags

    #enddef GetMacroDefinitionFlags

#endclass UtopyaCompilerSettings


# **


class UtopyaCompilerSettingsFortran( UtopyaCompilerSettings ) :

    """
    Class with methods that return Fortran compiler settings.
    """

    def GetCompilers( self, mpi=False, openmp=False ) :

        """
        Return actual names of compiler and linker to be
        used in make file.
        These names might depend on the flags that enable
        MPI and/or OpenMP enabled.
        The dictionairy that is returned has elements for the
        Fortran compiler, the linker (probably the same value),
        and for backwards compatability a Fortran77 specific version::

            { 'FC' : 'f90', 'F77' : 'f77', 'LINKER' : 'f90' }

        The names for compiler and linker are assumed to be the same,
        but his might be changed by derived classes.

        The name of the compiler/linker with MPI enabled is often
        defined by the MPI wrapper (MPICH, OpenMPI),
        but some compiler manufactures provide their own wrappers.

        OpenMP is either enabled by passing a flag to the standard
        compiler command (`--openmp`), but some compiler families
        use an alias instead.

        To support all these different styles, this method requires
        that the rcfile contains explicit names for all four possible
        combinations with/without MPI/OpenMP::

          ! fortran compiler:
          [<rcbase>.]configure.compiler.fc                       :  gfortran
          [<rcbase>.]configure.compiler.fc.openmp                :  gfortran

          ! compilers for MPI programs:
          [<rcbase>.]configure.compiler.fc.mpi                   :  mpif90
          [<rcbase>.]configure.compiler.fc.mpi.openmp            :  mpif90

          ! optional: f77 compiler (default same as fc):
          [<rcbase>.]configure.compiler.f77                      :  g77

        """

        # info ...
        self.logger.info( 'get compiler names ...' )

        # with mpi support ?
        if mpi :
            # supporting openmp ?
            if openmp :
                fc = self.GetSetting( 'configure.compiler.fc.mpi.openmp' )
            else :
                fc = self.GetSetting( 'configure.compiler.fc.mpi' )
            #endif
        else :
            # openmp enabled ?
            if openmp :
                fc = self.GetSetting( 'configure.compiler.fc.openmp' )
            else :
                fc = self.GetSetting( 'configure.compiler.fc' )
            #endif
        #endif

        # f77 compiler, by default the same as fc:
        f77 = self.GetSetting( 'configure.compiler.f77', default=fc )

        # assume linker is the same:
        linker = fc

        # info ...
        self.logger.debug( '  fortran compiler : %s' % fc )
        self.logger.debug( '  f77 compiler     : %s' % f77 )
        self.logger.debug( '  linker           : %s' % linker )

        # ok
        return { 'FC' : fc, 'F77' : f77, 'LINKER' : linker }

    #enddef GetCompilers

    # *

    def GetFlags( self, name ) :

        """
        Return compiler flags for specified name.
        Example settings for name ``optim-fast``::

          [<rcbase>.]configure.compiler.optim-fast.fflags        :  -O3
          [<rcbase>.]configure.compiler.optim-fast.ldflags       :

        Returns two str objects:

        * fflags   : flags passed to fortran compiler
        * ldflags  : flags passed to linker
        """

        # read from settings:
        fflags  = self.GetSetting( 'configure.compiler.%s.fflags'  % name )
        ldflags = self.GetSetting( 'configure.compiler.%s.ldflags' % name )

        # ok:
        return fflags,ldflags

    #enddef GetFlags

    # *

    def GetMacroDefinitionFlags( self, macros ) :

        """
        Return a line with command line arguments to define
        preprocessing macro's. Such a line probably looks like::

          -Dwith_this_flag -Dwith_that_flag ...

        Pass a list with macro's to be defined as argument.

        The following setting should contain the macro definition flag
        of the compiler::

          [<rcbase>.]configure.compiler.defineflag    :  -D

        """

        # compiler depended prefix for macro definition:
        Dflag = self.GetSetting( 'configure.compiler.defineflag' )

        # init empty:
        flags = ''
        # loop over macro's to be defined:
        for macro in macros :
            # add definition to command line argument list:
            flags = flags.strip()+(' %s%s' % (Dflag,macro) )
        #endfor

        # ok:
        return flags

    #enddef GetMacroDefinitionFlags

#endclass UtopyaCompilerSettingsFortran


# **


class UtopyaCompilerSettingsF2Py( UtopyaCompilerSettings ) :

    """
    Class with methods that return
    `F2Py <https://docs.scipy.org/doc/numpy/f2py/>`_
    compiler settings.
    The purpose of the ``F2Py`` (`Fortran to Python interface generator`)
    is to compile Fortran programs and create a Python interface
    to them.
    ``F2Py`` is an element of the Python's ``Numpy`` module,
    and therefore by default available in most distributions.

    """

    def GetCompilers( self, mpi=False, openmp=False ) :

        """
        Return actual names of compiler and linker to be
        used in make file.
        The dictionairy that is returned has only one element::

            { 'F2PY' : 'f2py' }

        The name of f2py compiler is read from
        the following rcfile setting::

          ! f2py compiler:
          [<rcbase>.]configure.compiler.f2py   :  f2py

        Flags for compilation with MPI and/or OpenMP are ignored.
        """

        # info ...
        self.logger.info( 'get compiler name ...' )

        # f2py compiler:
        f2py = self.GetSetting( 'configure.compiler.f2py' )

        # info ...
        self.logger.debug( '  f2py compiler : %s' % f2py )

        # ok
        return { 'F2PY' : 'f2py' }

    #enddef GetCompilers

    # *

    def GetFlags( self, name ) :

        """
        Return compiler flags for specified name.
        Example settings for name ``optim-fast``::

          [<rcbase>.]configure.compiler.optim-fast.fflags        :  -O3
          [<rcbase>.]configure.compiler.optim-fast.ldflags       :

        Returns two str objects:

        * fflags   : flags passed to fortran compiler
        * ldflags  : flags passed to linker
        """

        # read from settings:
        fflags  = self.GetSetting( 'configure.compiler.%s.fflags'  % name )
        ldflags = self.GetSetting( 'configure.compiler.%s.ldflags' % name )

        # ok:
        return fflags,ldflags

    #enddef GetFlags

    # *

    def GetMacroDefinitionFlags( self, macros ) :

        """
        Return a line with command line arguments to define
        preprocessing macro's. Such a line probably looks like::

          -Dwith_this_flag -Dwith_that_flag ...

        Pass a list with macro's to be defined as argument.

        The following setting should contain the macro definition flag
        of the compiler::

          [<rcbase>.]configure.compiler.defineflag    :  -D

        """

        # compiler depended prefix for macro definition:
        Dflag = self.GetSetting( 'configure.compiler.defineflag' )

        # init empty:
        flags = ''
        # loop over macro's to be defined:
        for macro in macros :
            # add definition to command line argument list:
            flags = flags.strip()+(' %s%s' % (Dflag,macro) )
        #endfor

        # ok:
        return flags

    #enddef GetMacroDefinitionFlags

#endclass UtopyaCompilerSettingsF2Py


#-------------------------------------------------
# library settings
#-------------------------------------------------

class UtopyaLibSettings( utopya_rc.UtopyaRc ) :

    """
    Base class for objects that should return library settings defined in a rcfile.
    """

    def GetFlags( self, name ) :

        """
        Return library flags for specified name. Example settings for name 'netcdf'::

          [<rcbase>.]configure.lib.netcdf.fflags     : -I${NETCDF_PREFIX}/include
          [<rcbase>.]configure.lib.netcdf.ldflags    : -L${NETCDF_PREFIX}/lib -lnetcdff -lnetcdf

        Return values:

        * fflags   : flags passed to Fortran compiler
        * ldflags  : flags passed to linker
        """

        # read from settings:
        fflags  = self.GetSetting( 'configure.lib.%s.fflags'  % name )
        ldflags = self.GetSetting( 'configure.lib.%s.ldflags' % name )

        # ok:
        return fflags,ldflags

    #enddef GetFlags

#endclass UtopyaLibSettings

#-------------------------------------------------
# makedep
#-------------------------------------------------

class UtopyaDependencies( utopya_rc.UtopyaRc ) :

    """
    Base class for automatic generation of Fortran source file dependencies.

    """

    def GetLines( self, flags='', files='*.[Ff]*', expl={} ) :

        """
        Return list with dependency lines for Fortran source.

        This method is actually an interface to the 'makedepf90' tool:

        * `MakeDepF90 Homepage <http://personal.inet.fi/private/erikedelmann/makedepf90/>`_

        To check if the tool is available, and for help on it's usage, try::

          makedepf90 --help

        The location of the tool could be specified (optionally) in the rcfile::

          makedepf90     :   /opt/makedepf90/bin/makedepf90

        Optional arguments:

        * 'flags': String with flags passed to 'makedepf90' tool,
           for example macro definition::

             -Dwith_some_library -D__VALUE__=12 ...

        * 'files': String with (list of) filename filters for source files
          for which dependecies should be made.
          Default is to use all Fortran files in the current directory.

        * 'expl': Dictionairy with explicit dependencies that are otherwise
          not recoqnized by 'makedepf90', for example files that
          are no modules. Example::

            { 'main.o' : ['tools.o','other.o',..], ... }

        Returns a list with dependency lines for a Makefile.
        Each line has the form::

            'file.o : dependency.o ...'

        """

        # tools:
        import utopya_tools

        # external program:
        prog = 'makedepf90'
        # eventually reset by optional path:
        prog = self.rcf.get( 'makedepf90', default=prog )

        # Check for the existence of makedepf90 on this machine:
        try :
            p = utopya_tools.Call( [prog,'--version'], logger=self.logger, verbose=False )
        except utopya_tools.CallingError as err :
            self.logger.error( err )
            self.logger.error( '')
            self.logger.error( '*******************************************************')
            self.logger.error( 'Cannot access a version of %s' % prog )
            self.logger.error( '(or something wrong with call to subprocess).')
            self.logger.error( 'If it is installed, please check your PATH variable.' )
            self.logger.error( 'If it is not installed, please install %s .' % prog )
            self.logger.error( 'For the source code and instruction, see:'  )
            self.logger.error( '  http://personal.inet.fi/private/erikedelmann/makedepf90/' )
            self.logger.error( '*******************************************************')
            self.logger.error( '')
            raise IOError( '%s not found' % prog )
        #endtry

        # command to generate a makefile:
        #command = '%s -Wconfused %s %s' % (prog,flags,files)
        # do not warn about confusion, gives too much noise ..
        command = '%s %s %s' % (prog,flags,files)
        # info ...
        self.logger.debug( '    ' )
        self.logger.debug( '    WARNING - Running makedepf90 ...' )
        self.logger.debug( '    WARNING - If the next command takes too much time and' )
        self.logger.debug( '    WARNING - seems to be in an infinite loop, try to run:' )
        self.logger.debug( '    WARNING -   dos2unix' )
        self.logger.debug( '    WARNING - on all your source files.' )
        self.logger.debug( '    WARNING - You probably used a Windows editor,' )
        self.logger.debug( '    WARNING - and makedepf90 does not like that ...' )
        self.logger.debug( '    ' )
        self.logger.debug( '    run command: %s' % command )
        # run command:
        try :
            # run as a shell command since the file list is probably '*.F90' :
            p = utopya_tools.Call( command, shell=True, logger=self.logger, verbose=False )
        except utopya_tools.CallingError as err :
            self.logger.error( err )
            raise Exception
        except utopya_tools.StatusError as err :
            for line in err.stderr : self.logger.error(line)
            self.logger.error( err )
            raise Exception
        #endtry

        # extract output lines:
        lines = p.stdout
        # translate into dictionary with lists,
        # start with empty:
        self.deps = {}
        # loop over lines:
        for line in lines :
            # should have a ':' as seperation:
            if ':' not in line :
                self.logger.error( 'found dependency line without ":" :' )
                self.logger.error( line )
                raise Exception
            #endif
            # split:
            target,objs = line.split(':',2)
            # cleanup:
            target = target.strip()
            # store in dictionairy as list:
            self.deps[target] = objs.split()
            # also explicitly defined dependencies ?
            if target in expl.keys() :
                # extend list:
                self.deps[target].extend( expl[target] )
            #endif
        #endfor

        # rewrite lines, including the explicitly added dependencies:
        lines = []
        for target in self.deps.keys() :
            # start:
            line = '%s :' % target
            # add objects:
            for obj in self.deps[target] : line = line+' '+obj
            # store:
            lines.append( line )
        #endfor

        # ok
        return lines

    #enddef   # GetLines

    # *

    def GetObjects( self, main ) :

        """
        Return character string of objects on which a main program depends.

        Arguments:

        * 'main': Basename of the main program, e.g. 'myprog' .
          The corresponding objectfile 'myprog.o' is the first value
          on the returned list. Other object names are added given
          the dependency information.
        """

        # first object file:
        obj = '%s.o' % main
        # init list of object files:
        objs = [obj]
        # loop until length is not changing anymore:
        nobj = len(objs)
        iloop = 1
        while True :
            # loop over current objects:
            for target in objs :
                # has dependencies ?
                if target in self.deps.keys() :
                    # loop over dependencies:
                    for obj in self.deps[target] :
                        # only object files ...
                        if not obj.endswith('.o') : continue
                        # add if not in list yet:
                        if obj not in objs : objs.append( obj )
                    #endfor
                #endif  # target with dependencies
            #endfor  # loop over current objects
            # new length:
            n = len(objs)
            # same ? then finished:
            if n == nobj : break
            # reset length:
            nobj = n
        #endwhile

        # sort the values:
        objs.sort()

        # write to string:
        line = ''
        for obj in objs : line = line+(' %s' % obj)

        # return clean list:
        return line.strip()

    #enddef

#endclass  # UtopyaDependencies

#-------------------------------------------------
# make
#-------------------------------------------------

class UtopyaMake( utopya_rc.UtopyaRc ) :

    """
    Base class to perform tasks to make an executable:

    * Change to source directory.
    * Make target.

    Example of rcfile settings::

        ! where to make?
        [<rcbase>.]make.dir             :  appl/model/src/

        ! targets to be build:
        [<rcbase>.]make.targets         :  model.x install

        ! Make command to be used for each of the targets;
        ! template '%{target}' is filled with the target name,
        ! and '%{njob}' is filled from the initialization argument:
        [<rcbase>.]make.target.model.x  :  make --file=Makefile --jobs=%{njob} %{target}
        [<rcbase>.]make.target.install  :  make --file=Makefile %{target}

    """

    def __init__( self, rcfile, rcbase='', env={}, njob=1 ) :

        # external:
        import os

        # tools:
        import utopya_tools

        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )

        # info ...
        self.logger.info( 'make executable ...' )

        # change directory:
        wd = self.GetSetting( 'make.dir' )
        # check ...
        if not os.path.isdir( wd ) :
            self.logger.error( 'workdir for make not found: %s' % wd )
            raise Exception
        #endif
        # info ...
        self.logger.info( '  change to %s ...' % wd )
        # store old:
        owd = os.getcwd()
        # change:
        os.chdir( wd )

        # build targets:
        targets = self.GetSetting( 'make.targets' ).split()
        # loop:
        for target in targets :

            # info ...
            self.logger.info( '  make "%s" ...' % target )

            # get maker command:
            maker = self.GetSetting( 'make.target.%s' % target )
            # replace templates:
            maker = maker.replace('%{target}',target)
            maker = maker.replace('%{njob}',str(njob))

            # command as list:
            command = maker.split()
            # info ...
            self.logger.info( '    run command: %s' % str(command) )
            # run, send output to logger ;
            # send error messages to output since sometimes warnings are
            # issued as errors and these could be better be shown immediately:
            p = utopya_tools.Call( command, logger=self.logger, err2out=True )

        #endfor # targets

        # back:
        os.chdir( owd )

    #enddef __init__

#endclass UtopyaMake


#-------------------------------------------------
# end
#-------------------------------------------------

