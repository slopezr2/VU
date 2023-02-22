#! /usr/bin/env python

"""
*            PYCASSO - PYthon Compile And Setup Scripts Organizer            *

INSTALL FOR YOUR APPLICATION

    Assume you have an application named 'yourmodel'.
    To create a compile/setup script for this application,
    copy the template setup script to a name that you and other users
    will recoqnize as the main setup script:

       cp  pycasso/py/pycasso_setup_template  setup_yourmodel

    Edit the 'Settings' section in this new script to set the location 
    of the PYCASSO scripts. Also specify the name of the module file
    with user scripts. If these are in a subdirectory './base/py', use fore example:

       # path to the PYCASSO modules:
       base_py = os.path.join('.','pycasso','py')

       # name of file with PYCASSO user scripts for this setup:
       pycasso_user_scripts = 'pycasso_user_scripts_tm5'

    Copy the template settings file 'pycasso__template.rc'
    to a suitable name:

        cp  pycasso/rc/pycasso_template.rc  yourmodel.rc

    Edit the settings if necessary, and start with:

      ./setup_yourmodel yourmodel.rc


LOGGING

    Logging is implemented via the standard python 'logging' module.
    
    By default, one logger is defined that writes to the standard output;
    usually this is the terminal window, but in a batch environment this
    could be redirected to a file.
    In addition, a second logger could be setup if the 'logfile'
    option in the rc file is set.
    
    To add a new message to one of the logs, follow the examples in the code:

        logging.info( 'rc-file read successfully' )

    This would create a log-message of level INFO.
    Other options for log messages are:

        logging.debug    (msg)     # extra information to debug the scripts
        logging.warning  (msg)     # warnings about undesired behaviour
        logging.info     (msg)     # the standard messages
        logging.exception(msg)     # something wrong, not fatal
        logging.error    (msg)     # something wrong, nearly fatal
        logging.critical (msg)     # something wrong, probably fatal

    The threshold level for messages to be included is DEBUG for the file output.
    The threshold level for messages to be included is INFO for the screen, 
    unless 'options.verbose' is True (set by '-v' or '--verbose' on the 
    command line).

"""


# ------------------------------------------------
# errors
# ------------------------------------------------


class PycassoError( Exception ) :
    """
    Base class for exceptions raised by Pycasso scripts.
    Raised after displaying proper error messages.
    """
#endclass


# ------------------------------------------------
# begin
# ------------------------------------------------


def Main( args, pycasso_user_scripts ) :

    """
    Start the compilation and setup of the application.
        args
                List of unparsed arguments, probalby from 'sys.argv[1:]'.
        pycasso_user_scripts
                Name of file with user scripts for this particular setup.
    """
    
    #
    # * arguments
    #
    
    # modules:
    import sys
    import logging
    
    # first parse the arguments to avoid that other messages
    # end up in the help text ..
    
    # parse the arguments;
    # returns a dictionairy 'options' and a single name of the rcfile:
    # trap error status if help text is requested:
    try :
        options,rcfile = ParseArguments( args, pycasso_user_scripts )
    except SystemExit :
        # help text was printed; return without exception:
        return
    except :
        # unknown exception ...
        logging.error( sys.exc_info()[1] )
        logging.error( 'exception from pycasso.Main' )
        raise PycassoError
    #endif

    
    #
    # * modules
    #
    
    # external:
    import os
    import subprocess
    import shutil
    import traceback
    import datetime
    
    #
    # * user scripts
    #
    
    # try to load user scripts and display proper error messages:
    try :
        pycus = __import__( pycasso_user_scripts )
    except ImportError :
        logging.error( 'could not import module : %s' % pycasso_user_scripts )
        logging.error( 'check definition of "pycasso_user_scripts" in "setup" script' )
        raise PycassoError
    except :
        logging.error( traceback.format_exc() )
        logging.error( 'error in user scripts module ?' )
        logging.error( 'check content of module : %s' % pycasso_user_scripts )
        raise PycassoError
    #endtry
    
    #
    # * setup logging to standard output
    #
    
    # initialise logging system:
    logger = Start_Logger()
    
    # initialise messages to standard output:
    stdout_handler = Start_StdOut_Log(logger)
    
    # info ...
    logging.info( '' )
    tnow = datetime.datetime.now().isoformat(' ').split('.',1)[0]
    logging.info( 'Started script at %s ...' % tnow )
    logging.info( '' )
    
    # info ...
    logging.info( 'options and arguments ...' )
    logging.info( '  parsed options         : %s' % options )
    logging.info( '  parsed rcfile argument : %s' % rcfile )
    
    #
    # * settings from rcfile
    #

    # tools:    
    import rc

    # info ...
    logging.info( 'read settings from %s ...' % rcfile )
    
    # read settings:
    rcf = rc.RcFile(rcfile)
    
    # options provided at command line ?
    if len(options) == 0 :
        # info ...
        logging.info( 'no rcfile settings to be added or changed by options' )
    else :
        # info ...
        logging.info( 'change rcfile settings given options ...' )
        # add or replace settings from the arguments passed to the script:
        for key in options.keys() :
            # value:
            val = options[key]
            # add or replace ?
            if rcf.has_key(key) :
                # replace existing value:
                logging.info( '  reset "%s" to "%s" ...' % (key,str(val)) )
                rcf.replace( key, str(val) )
            else :
                # add new key:
                logging.info( '  set "%s" to "%s" ...' % (key,str(val)) )
                rcf.add( key, str(val), comment='new key added after options passed to setup script' )
            #endif
        #endfor
    #endif
    
    #
    # * apply settings for logging
    #
    
    logging.info( 'setup logging according to settings ...')
    
    # debug messages to standard output ?
    flag = rcf.get( 'verbose', 'bool', default=False )
    if flag : 
        logging.info( '  verbose mode for standard output; print all messages ...' )
        stdout_handler.setLevel(logging.DEBUG)
    else :
        logging.info( '  quiet mode for standard output; print info messages only ...' )
    #endif
    
    # setup logfile ?
    key = 'logfile'
    if rcf.has_key(key) :
        logfile = rcf.get(key)
        logging.info( '  open additional log file : %s' % logfile )
        logfile_handler = Start_File_Log( logger, logfile )
    else :
        logging.info( '  no log file; print to standard output only ...' )
        logfile_handler = None
    #endif

    # test logging:
    logging.info    ('  test messages ...')
    logging.debug   ('    testing debug    message ...')
    logging.info    ('    testing info     message ...')
    logging.warning ('    testing warning  message ...')
    logging.error   ('    testing error    message ...')
    logging.critical('    testing critical message ...')
    
    #
    # * environment module
    #
    
    # module name ? might be empty, this is also the default if setting is not present:
    modfile = rcf.get( 'env.module.file', default='' )
    # defined ?
    if len(modfile) > 0 :
        # load:
        Env_Module( rcf, modfile )
    #endif
    
    #
    # * display all settings
    #
    
    # display settings:
    logging.debug( '' )
    logging.debug( 'Content of rc dictionary:' )
    for key in rcf.keys() :
        # get value as a str to avoid warnings about bool values:
        val = rcf.get(key,'str')
        # display:
        logging.debug( '  %s : %s' % (key,val) )
    #endfor
    logging.debug( '[done]' )
    
    #
    # * tasks
    #
    
    # copy a source ?
    flag = rcf.get('build.copy','bool')
    if flag :
        # info ...
        logging.info( 'copy source to build directory ...' )
        # create source:
        Build_Copy( rcf, pycasso_user_scripts )
    else :
        # info ...
        logging.info( 'no source to be copied ...' )
    #endif

    # configure a source ?
    flag = rcf.get('build.configure','bool')
    if flag :
        # info ...
        logging.info( 'configure source ...' )
        # configure source:
        Build_Configure( rcf, pycasso_user_scripts )
    else :
        # info ...
        logging.info( 'no source to be configured ...' )
    #endif

    # make an executable ?
    flag = rcf.get('build.make','bool')
    if flag :
        # info ...
        logging.info( 'make source ...' )
        # compile source:
        Build_Make( rcf, pycasso_user_scripts )
    else :
        # info ...
        logging.info( 'no source to be made ...' )
    #endif

    # change to destination directory ?  
    rundir = rcf.get('rundir')
    if len(rundir) > 0 :
        # create target directory if necessary:
        CreateDirs( rundir )
        # info ...
        logging.info( 'change to run directory %s ...' % rundir )
        # goto this directory:
        os.chdir(rundir)
    else :
        # info ...
        logging.info( 'no run directory specified; stay here ...' )
    #endif

    # copy files ?
    ifiles = rcf.get('install.copy')
    if len(ifiles) > 0 :
        # info ...
        logging.info( 'install files ...' )
        # loop over files to be installed:
        for ifile in ifiles.split() :
            # if the file contains a colon ':' ...
            if ':' in ifile :
                # ... the name after the colon is the target name
                sour,targ = ifile.split(':')
            else :
                # ... otherwise, copy to current directory:
                sour,targ = ifile,os.path.curdir
            #endif
            # info ...
            logging.info( '  copy "%s" to "%s" ...' % (sour,targ) )
            # check ...
            if not os.path.exists(sour) :
                logging.error( 'source file not found : %s' % sour )
                raise IOError
            #endif
            # use shell utility:
            shutil.copy( sour, targ )
        #endif
    else :
        # info ...
        logging.info( 'no files to be installed ...' )
    #endif
    
    # write rcfile ...
    newrc = rcf.get('install.rc')
    if len(newrc) > 0 :
        # info ...
        logging.info( '  install processed rcfile ...' )
        # write pre-processed rcfile:
        rcf.WriteFile( newrc )
    #endif

    #
    # * submit script / info / submitting
    #
    
    # tools, imported from build/bin directory:
    import pycasso_tools
    
    # name of submit script:
    submit_script = rcf.get( 'submit.script' )

    # defined ?
    if len(submit_script) > 0 :

        # full command:
        submit_command = rcf.get( 'submit.command' )

        # where to search for scripts ?    
        build_prefix = rcf.get('build.prefix')
        # use paths relative to run directory ?
        if rcf.get('submit.relpaths','bool') :
            # only available in recent versions ...
            if sys.version_info[0]+0.1*sys.version_info[1] >= 2.6 :
                build_prefix = os.path.relpath(build_prefix,start=rundir)
            else :
                logging.warning( "Option `submit.relpaths` requires at least python version 2.6 ; use absolute path instead" )
            #endif
        #endif

        # full path to submit script:
        submit_script_path = os.path.join( rundir, submit_script )
        # should exist ...
        if not os.path.isfile(submit_script_path) :
            logging.error( 'submit script not found:' )
            logging.error( '  %s' % submit_script_path )
            logging.error( 'not included in "install.copy" list ?' )
            raise PycassoError
        #endif
        # insert path to submit modules:
        pycasso_tools.modify_text_file( submit_script_path, 
                         "pypath_default = os.path.join( os.pardir, 'build', 'bin' )", 
                         "pypath_default = os.path.join('%s','bin')" % build_prefix )

        # info ...
        logging.info( '' )
        logging.info( 'To submit a job:' )
        logging.info( '' )
        indent = '  '
        # need to change to run directory ?
        if len(rundir) > 0 :
            logging.info( indent+'# change to the run directory:' )
            logging.info( indent+'cd %s' % rundir )
            logging.info( '' )
        #endif
        # disply usage text:
        logging.info( indent+'# submit a job (use --help for more information):' )
        logging.info( indent+'%s' % submit_command )

        # extra message ?
        submit_msg = rcf.get( 'submit.msg', default='None' )
        if submit_msg != 'None' :
            for line in submit_msg.split('\\n') :
                logging.info( line.strip().replace('\\t',indent) )
            #endfor
        #endif

        # submit automatically ?
        flag = rcf.get( 'submit.auto', 'bool' )
        if flag :
            # change to run directory:
            os.chdir( rundir )
            # info ...
            logging.info( '' )
            logging.info( 'submit automatically ...' )
            logging.info( '  %s' % submit_command )
            logging.info( '' )
            # call script:
            retcode = subprocess.call( submit_command.split() )
            if retcode != 0 :
                logging.error( 'from submission of : %s' % submit_command )
                raise PycassoError
            #endif
        #endif

    #endif   # submit script defined
    
    
    #
    # * end
    #
    
    # info ...
    logging.info( '' )
    tnow = datetime.datetime.now().isoformat(' ').split('.',1)[0]
    logging.info( 'End of script at %s .' % tnow )
    logging.info( '' )
    
    # close logs:
    if logfile_handler != None : logfile_handler.close()
    stdout_handler.close()
    
#enddef   # Main


# ***


def ParseArguments( args, pycasso_user_scripts ) :

    """
    Define the arguments accepted by a pycasso run script.

    Usage:

        options,rcfile = ParseArguments( args )

    Arguments:

        args
                Passed from the calling script, probably equal to : sys.argv[1:]
        pycasso_user_scripts
                Name of file with user scripts for this particular setup.

    Return values:

        options       # object with following data fields:
          .verbose    # (boolean) gives extra logging messages to the screen

        rcfile    # name of settings file
      
    """

    # external:
    import logging
    import optparse
    
    # load user scripts:
    pycus = __import__( pycasso_user_scripts )
    
    # set text for 'usage' help line:
    usage = "\n    %prog [options] rcfile\n    %prog -h|--help"
    
    # descriptive text
    description = "Driver script to compile and setup a model application. "\
        "The 'rcfile' is a textfile with settings read by the scripts or "\
        "the application, a template should be available with this script."
    
    # initialise the option parser:
    parser = optparse.OptionParser( usage=usage, description=description )
    
    # define verbose option:
    parser.add_option( "-v", "--verbose", 
                         help="""Print extra logging messages to standard output.
                            This option will set rcfile key 'verbose' to 'True'.""",
                         dest="verbose", action="store_true", default=False )
    # new build ?
    parser.add_option( "-n", "--new", 
                         help="""Create new build; remove old build directory.
                            This option will set rcfile key 'build.new' to 'True'.""",
                         dest="build_new", action="store_true", default=False )
    # how many jobs used for make etc ?
    parser.add_option( "-j", "--jobs", 
                         help="""Number of jobs (commands) to run simultaneously.
                            Empty value '' indicates unlimitted number.
                            Now only used to (re)set the number of jobs used by the maker ('gmake -j JOBS') .
                            This flag will replace the value of 'build.jobs' in the rcfile.""",
                         dest="jobs", action="store", default=None )
    # submit job ?
    parser.add_option( "-s", "--submit", 
                         help="""Submit the job after setup. 
                            See also the section on 'Submit options' below for 
                            options passed directly to the submit script.
                            This option will set rcfile key 'submit.auto' to 'True'.""",
                         dest="submit_auto", action="store_true", default=False )

    # options for submitting the job:
    group = optparse.OptionGroup( parser, "Submit options",
                         description="Options passed directly to the submit script, see its help text for details." )
    # where to submit to ?
    group.add_option( "-f", "--foreground", 
                         help="Run job in foreground.",
                         dest="submit_to", action="store_const", const='foreground' )
    group.add_option( "-b", "--background", 
                         help="Run job in background.",
                         dest="submit_to", action="store_const", const='background' )
    group.add_option( "-q", "--queue", 
                         help="Submit job to a queue system.",
                         dest="submit_to", action="store_const", const='queue' )
    # when submitted, run in debugger ?
    group.add_option( "-d", "--debugger",
                         help="Run executable in a debugger.",
                         dest="submit_debugger", action="store_true" )
    # add group:
    parser.add_option_group( group )
   
    # add the user options:
    group = optparse.OptionGroup( parser, "User model options",
                         description="Extra options are defined and handled in the 'pycasso_user_scripts_*' module, might be none." )
    pycus.DefineOptions( group )
    parser.add_option_group( group )

    # now parse the actual arguments:
    values,args = parser.parse_args( args=args )
    
    # at least rcfile should be specified as argument,
    # and no other arguments:
    if len(args) != 1 :
        parser.error("incorrect number of arguments\n")
    #endif
    
    # translate options into a dictionairy 
    # if they should replace rcfile values:
    opts = {}
    if values.verbose                 : opts['verbose'        ] = True
    if values.build_new               : opts['build.new'      ] = True
    if values.jobs != None            : opts['build.jobs'     ] = values.jobs
    if values.submit_auto             : opts['submit.auto'    ] = True
    if values.submit_to != None       : opts['submit.to'      ] = values.submit_to
    if values.submit_debugger != None : opts['submit.debugger'] = values.submit_debugger

    # add the parsed user options:
    pycus.StoreOptions( opts, values )
    
    # copy name of rcfile:
    rcfile = args[0]

    # return values:
    return opts,rcfile
    
#enddef


# ***


def Start_Logger() :

    """
    logger = Start_Logger()
    """

    # external:
    import logging

    # create logger
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    
    # ok
    return logger

#enddef


# ***


def Start_StdOut_Log(logger) :

    """
    stdout_handler = Start_StdOut_Log(logger)
    """

    # external:
    import sys
    import logging

    # set a format for screen use:
    logformat = '[%(levelname)-8s] %(message)s'
    # create formatter:
    formatter = logging.Formatter(logformat)
    
    # handler for standard output:
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    ## first messages:
    #logging.debug('testing debug message after start stdout log ...')
    #logging.info ('testing info  message after start stdout log ...')
    
    # ok
    return stdout_handler
    
#endif


# ***


def Start_File_Log(logger,logfile) :

    """
    logfile_handler = Start_File_Log(logger,logfile)
    
    Create handler for log file.
    Initial level is 'DEBUG', thus all messages will be written to the file.
    """
    
    # external:
    import sys
    import logging
    
    # set format of lines written to logging:
    if sys.version_info < (2, 5):
        logformat   = "%(asctime)s %(name)-10s, line %(lineno)4i  [%(levelname)-10s] %(message)s"
    else:
        logformat   = '%(lineno)4i %(filename)-30s -> %(funcName)-30s [%(levelname)-8s] %(message)s'
    #endif
    # create formatter:
    formatter = logging.Formatter(logformat)
    
    # now create a handler for the log file;
    # mode 'w' will cause the log file to be re-written:
    logfile_handler = logging.FileHandler(filename=logfile,mode='w')
    logfile_handler.setLevel(logging.DEBUG)
    logfile_handler.setFormatter(formatter)
    logger.addHandler(logfile_handler)
    
    ## first messages:
    #logging.debug('testing debug message after start file log ...')
    #logging.info ('testing info  message after start file log ...')
    
    # ok
    return logfile_handler

#enddef


# ***


def Env_Module( rcf, modname ) :

    # modules:
    import logging
    import os
    import subprocess
    
    # tools:
    import EnvMod
    
    # info ...
    logging.info( 'load environment module "%s" ...' % modname )

    # flag:
    found = False
    # get list of sub directories to be scanned:
    sourcedirs = rcf.get('build.copy.dirs').split()
    # loop:
    for sourcedir in sourcedirs :
        # requested module file is in the 'env' subdir:
        envdir = os.path.join( sourcedir, 'env' )
        # present ?
        if os.path.isfile( os.path.join(envdir,modname) ) :
            # store the current directory,  might be the latest found:
            usedir = envdir
            # reset flag:
            found = True
        #endif
    #endfor
    # not found ? problem ...
    if not found :
        logging.error( 'could not find module file "%s" in any of the environment directories:' % modname )
        for sourcedir in sourcedirs :
            logging.error( '  - %s' % os.path.join(sourcedir,'env') )
        #endfor
        raise PycassoError
    #endif
    
    # info ...
    logging.info( '  use module file from subdir "%s" ...' % usedir )
    # add to search path, load, and remove from path:
    EnvMod.Module( 'use', usedir )
    EnvMod.Module( 'load', modname )
    EnvMod.Module( 'unuse', usedir )
    
#enddef  # Env_Module


# ***


def Build_Copy( rcf, pycasso_user_scripts ) :

    # modules:
    import os
    import sys
    import shutil
    import logging
    import re
    
    # load user scripts:
    pycus = __import__( pycasso_user_scripts )
    
    # remove existing build ?
    remove_existing_build = rcf.get('build.new','bool')

    # target directory:
    prefix = rcf.get('build.prefix')
    # info ...
    logging.info( '  build prefix : %s ' % prefix )
    
    # extend name ?
    if rcf.get('build.prefix.extend','bool') :
        # start with original prefix:
        prefix_ext = prefix
        # get list with names of compiler flag groups to be used:
        flags = pycus.Build_FlagGroups( rcf )
        # add to prefix:
        if len(flags) > 0 :
            for flag in flags : prefix_ext = prefix_ext+'_'+flag
        else :
            prefix_ext = prefix_ext+'_'
        #endif
        # info ...
        logging.info( '  build prefix extended : %s ' % prefix_ext )
        # create extended prefix directory if necessary:
        CreateDirs( prefix_ext, forceclean=remove_existing_build )
        # test if link already exists; use 'lexists' to trap broken link:
        if os.path.lexists(prefix) :
            # is a link ?
            if os.path.islink(prefix) :
                # remove current link:
                os.remove( prefix )
            else :
                logging.error( 'could not replace "'+prefix+'" by a symbolic link; remove first' )
                raise PycassoError
            #endif
        #endif
        # create link relative to current directory:
        os.symlink( os.path.basename(prefix_ext), prefix )
    else :
        # create prefix directory if necessary:
        CreateDirs( prefix, forceclean=remove_existing_build )
    #endif
    
    # get list of sub directories to be scanned:
    subdirs = rcf.get('build.copy.subdirs').split()
    # eventually loop over sub directories:
    for subdir in subdirs :
        # full path:
        sdir = os.path.join(prefix,subdir)
        # create if necessary:
        CreateDirs( sdir )
    #endfor
    
    # loop over source directories:
    sourcedirs = rcf.get('build.copy.dirs')
    if len(sourcedirs) > 0 :
        # info ...
        logging.info( '  copy files from source directories...' )
        # some flags ...
        flag_remove__part = rcf.get('build.copy.remove.__part','bool')
        if flag_remove__part :
            logging.info( '    remove "__<name>" parts from sources files ...' )
        #endif
        # loop over source directories:
        for sourcedir in sourcedirs.split() :
            found_some_files = False
            # info ...
            logging.info( '    scanning %s ...' % sourcedir)
            # should be a directory ...
            if not os.path.isdir(sourcedir) :
                logging.error( 'specified source dir is not an existing directory : %s' % sourcedir )
                raise IOError
            #endif
            # empty ? then add something for current directory:
            if len(subdirs) == 0 : subdirs = os.path.curdir
            # loop over sub directories:
            for subdir in  subdirs :
                # full path:
                if subdir == '.' :
                    # just the directory:
                    sourcepath = sourcedir
                else :
                    # add sub directory:
                    sourcepath = os.path.join(sourcedir,subdir)
                    # info ...
                    logging.debug( '      scanning %s ...' % sourcepath )
                #endif
                # list all files in this directroy:
                try :
                    # list all files, might be empty:
                    sourcefiles = os.listdir(sourcepath)
                except :
                    # empty ? then next
                    continue
                #endtry
                # loop over source files:
                for sfile in sourcefiles :
                    found_some_files = True
                    # full filename:
                    sourcefile = os.path.join(sourcepath,sfile)
                    # skip directories:
                    if os.path.isdir(sourcefile) : continue

                    # trap unsaved emacs files:
                    # if file 'xxx.F90' is changed in emacs but not saved,
                    # a temporary unexisting link '.#xxx.F90' is made;
                    # the copy command will fail for this file, so therefore
                    # issue an error on forehand:
                    if sfile.startswith('.#') :
                        logging.error( 'file "%s" is being editted by Emacs; please save!' % sfile.lstrip('.#') )
                        raise PycassoError
                    #endif

                    # skip files with certain extensions:
                    name,ext = os.path.splitext(sfile)
                    skipit = False
                    for pattern in rcf.get('build.copy.skip.ext').split() :
                        if re.search(pattern, ext) :
                            logging.debug( '        %-50s     %-30s [%-8s]' % (sourcefile,'','skip') )
                            skipit = True
                            break
                        #endif
                    #endfor
                    if skipit : continue

                    # default target file:
                    outfile  = sfile
                    # remove '__name' part ?
                    if flag_remove__part :
                        # sourcefile contains '__' ?
                        if '__' in sfile :
                            # strip '__xxx' part from file name (used to identify versions):
                            name,ext = os.path.splitext(sfile)
                            outfile  = name.split('__')[0]+ext
                        #endif
                    #endif   # remove __part
                    # skip file ?
                    if outfile in rcf.get('build.copy.skip.file').split() :
                        logging.debug( '        %-50s     %-30s [%-8s]' % (sourcefile,'','skip') )
                        continue
                    #endif
                    # full target file:
                    targetfile = os.path.join( prefix, subdir, outfile )
                    # already present ?
                    if os.path.exists(targetfile) :
                        # different ?
                        if diff_text_files(sourcefile,targetfile) :
                            stat = 'differs'
                        else: 
                            stat = '...'
                        #endif
                    else: 
                        stat = 'new'
                    #endif
                    # info ...
                    logging.debug( '        %-50s  -> %-30s [%-8s]' % (sourcefile,os.path.join(subdir,outfile),stat) )
                    # copy source to target, preserve times etc:
                    shutil.copy2( sourcefile, targetfile )
                #endfor  # source files
            #endfor   # subdirectories
            # check ...
            if not found_some_files :
                logging.error('  found no source files in standard subdirs %s of %s.  Mistake in source.dirs?' % (subdirs,sourcedir) )
                raise Exception
            #endif

        #endfor   # source directories
    else :
        # info ...
        logging.info( '  no source directories specified ...' )
    #endif

    # add a new formed subdirs to pythonpath if necessary:
    subdirs = rcf.get( 'build.copy.subdirs.pypath' ).split()
    # loop over all subdirs:
    for subdir in subdirs :
        # info ...
        logging.info( '  add subdir <prefix>/%s to python path ...' % subdir )
        # extra directory:
        newdir = os.path.join(prefix,subdir)
        # add as first directory to search path:
        sys.path.insert(0,newdir)
    #endfor

#enddef


# ***


def Build_Configure( rcf, pycasso_user_scripts ) :

    # external:
    import os
    import sys
    import logging
    import subprocess
    import glob
    
    # tools:
    import go
    import pycasso_tools
    import MakeDep
    
    # load user scripts; try to reload first,
    # the path has been extended with subdirs from the build directory
    # (build/bin for example), so this version should replace an existing one:
    try :
        reload( pycus )
    except :
        pycus = __import__( pycasso_user_scripts )
    #endtry
    
    # change directory:
    current_dir = os.getcwd()
    configure_dir = rcf.get('build.sourcedir')
    logging.debug( '  change to %s ...' % configure_dir )
    os.chdir( configure_dir )

    #
    # * compiler and flags
    #
    
    # call user script to set compilers and linker:
    fc,f77,linker = pycus.Build_Compiler( rcf )
    
    # info...
    logging.info( '  collect compiler flags ...' )
    
    # start without any flags:
    fflags  = ''
    ldflags = ''
    libs    = ''
    
    # get list with names of compiler flag groups to be used:
    flaggroups = pycus.Build_FlagGroups( rcf )
    # loop over groups:
    for flaggroup in flaggroups :
        # info ...
        logging.info( '    add flag group "%s" ...' % flaggroup )
        # add flags for this group:
        fflags  =  fflags.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.fflags' )
        ldflags = ldflags.strip()+' '+rcf.get('compiler.flags.'+flaggroup+'.ldflags')
    #endfor
    logging.info( '      fflags       : %s' % fflags )
    logging.info( '      ldflags      : %s' % ldflags )
    

    #
    # * macro's
    #

    # default defined macro's:    
    macros_def   = rcf.get('build.configure.macro.define' ).split()

    # macro groups:
    groups =  rcf.get('build.configure.macro.groups').split()

    # apply user changes to default list:
    for group in groups :
        # add (or remove!) macro defintions by user script:
        macros_def = pycus.Build_Define( rcf, group, macros_def )
    #endfor
    
    # create a 'clean' list of defined macro's without double definitions:
    macros_defined   = []
    for m in macros_def :
        if m not in macros_defined   : macros_defined   = macros_defined   + macros_def
    #endfor

    # initialize a list to store all supported macro's without double definitions:
    macros_supported = []

    # loop over groups:
    for group in groups :
        # start of the rc keys:
        keybase = 'build.configure.macro.'+group
        # values for this group:
        macros_all   = rcf.get(keybase+'.all'    ).split()
        macros_hfile = rcf.get(keybase+'.hfile'  )
        # write header file ?
        if len(macros_hfile) > 0 :
            # info ...
            logging.info( '  update %s (include file with macro definitions) ...' % macros_hfile )
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
            pycasso_tools.update_text_file( macros_hfile, src )
        #endif
        # extend list:
        for m in macros_all :
            if m not in macros_supported : macros_supported = macros_supported + macros_all
        #endfor
    #endif
    
    # check for macro's that are not supported yet:
    any_error = False
    for macr in macros_defined :
        # not supported ?
        if macr not in macros_supported :
            # any unsupported macro's found yet ?
            if not any_error :
                # initial error message:
                logging.error( "one or more macro's have been defined that are not listed" )
                logging.error( "in any 'build.configure.macro.*.all' lists:" )
            #endif
            # display problematic macro:
            logging.error( "  %s" % macr )
            # reset flag:
            any_error = True
        #endif
    #endfor
    # any error ? then leave:
    if any_error : raise PycassoError

    # create list of macro definitions as command line arguments, e.g. -Dwith_this_flag etc:
    fc_defs = ''   # for fortran compiler
    mk_defs = ''   # for makedepf90
    # add macro definitions ?
    define_D = rcf.get( 'build.configure.define.D', 'bool' )
    if define_D :
        # compiler depended prefix for macro definition:
        fc_D = rcf.get('compiler.defineflag',default='-D')
        mk_D = '-D'
        # loop over macro's to be defined:
        for mdef in macros_defined :
            # add definition to command line argument list:
            fc_defs = fc_defs.strip()+(' %s%s' % (fc_D,mdef) )
            mk_defs = mk_defs.strip()+(' %s%s' % (mk_D,mdef) )
        #endfor
    #endif
    # add definitions to flags:
    fflags = fflags.strip()+' '+fc_defs
    

    #
    # * libraries
    #
    
    # get list of default library names to be linked:
    libnames = rcf.get( 'build.configure.libs' )
    
    # info ...
    logging.debug( '  libraries to be used:' )
    if len(libnames) > 0 : logging.debug( '    %s    (default)' % libnames )
    
    # convert to list:
    libnames = libnames.split()

    # loop over defined macro's:
    for mdef in macros_defined :
        # read list of libraries that should be added if this macro is defined:
        libnames_ifdef = rcf.get( 'build.configure.libs.ifdef.%s' % mdef, default='' )
        # defined ?
        if len(libnames_ifdef) > 0 :
            # add:
            libnames = libnames+libnames_ifdef.split()
            # info ...
            logging.debug( '    %s    (macro `%s` defined)' % (libnames_ifdef,mdef) )
        #endif
    #enddef
    
    # get list of all supported library names:
    libnames_all = rcf.get( 'build.configure.libs.all' ).split()
    
    # check if all libraries specified to be used are actually supported ...
    for libname in libnames :
        if libname not in libnames_all :
            logging.error( 'library name `%s` not in `build.configure.libs.all` list ...' % libname )
            raise PycassoError
        #endif
    #endfor
    
    # info ...
    logging.debug( '  libraries linked (in this order):' )
    # now add compiler and linker flags ;
    # loop over all supported libraries (this specfifies the linking order!)
    for libname in libnames_all :
        # not in use ? then skip:
        if libname not in libnames : continue
        # info ...
        logging.debug( '    %s' % libname )
        # add include, module, and link flags:
        fflags = fflags.strip()+' '+rcf.get('compiler.lib.'+libname+'.fflags')
        libs   =   libs.strip()+' '+rcf.get('compiler.lib.'+libname+'.libs')
    #endfor
    

    #
    # * write compiler flags to makefile include file:
    #
    
    # name of include file with dependencies, to be written by this script:
    makefile_flags = rcf.get('build.configure.flags.includefile')
    # info ...
    logging.info( '  write %s (compiler and flags) ...' % makefile_flags )
    logging.info( '    compiler fflags : '+fflags  )
    logging.info( '            ldflags : '+ldflags )
    logging.info( '               libs : '+libs    )
    # fill content; each line should end with the newline expression '\n' :
    src = []
    src.append( '#\n' )
    src.append( '# include file with compiler flags for Makefile.\n' )
    src.append( '#\n' )
    src.append( '\n' )
    src.append( '# compiler and linker:\n' )
    src.append( 'FC = %s\n' % fc )
    src.append( 'F77 = %s\n' % f77 )
    src.append( 'LINKER = %s\n' % linker )
    src.append( '\n' )
    src.append( '# compile flags:\n' )
    src.append( 'FFLAGS = %s\n' % fflags )
    src.append( '\n' )
    src.append( '# linker flags:\n' )
    src.append( 'LDFLAGS = %s\n' % ldflags )
    src.append( '\n' )
    src.append( '# library flags:\n' )
    src.append( 'LIBS = %s\n' % libs )
    src.append( '\n' )

    # write:
    pycasso_tools.write_text_file( makefile_flags, src )
    

    #
    # * write compile rules to makefile include file:
    #

    # name of include file with compiler rules, to be written by this script:
    makefile_rules = rcf.get('build.configure.rules.includefile')
    # info ...
    logging.info( '  write %s (compile rules) ...' % makefile_rules )
    # init content; each line should end with the newline expression '\n' :
    src = []
    # header:
    src.append( '#\n' )
    src.append( '# include file with compile rules for Makefile.\n' )
    src.append( '#\n' )
    src.append( '\n' )
    # explicit rules for some files 
    # in rcfile, a list is specified with objects and additional flag groups:
    #      objectfile1  :  flaggroup1 ; objectfile2 : flaggroup2 ; ...
    # read list:
    rlist = rcf.get( 'build.configure.rules.explicit' )
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
            fflags_expl = fflags_expl+' '+rcf.get('compiler.flags.'+grp+'.fflags' )
        #endfor
        fflags_expl = fflags_expl.strip()
        # add rule:
        src.append( '%s:\n' % obj )
        src.append( '\t$(FC) -c -o $@ $(FFLAGS) %s $<\n' % fflags_expl )
        src.append( '\t@echo " "\n' )
        src.append( '\n' )
    #endfor 

    # write:
    pycasso_tools.write_text_file( makefile_rules, src )
    

    #
    # * cleanup
    #
    
    # remove files ?
    logging.info( '  remove files ...' )
    for sfile in rcf.get('build.configure.remove').split() :
        # info ...
        logging.debug( '    %s ...' % sfile )
        # remove file if present:
        if os.path.exists(sfile) : os.remove( sfile )
    #endfor

    # clean for defined macro's (read list of pattern)
    logging.info( '  remove files for defined macro\'s ...' )
    for macr in macros_defined :
        # try to read list with files to be removed if macro is defined:
        patterns = rcf.get( 'build.configure.remove.ifdef.%s' % macr, default='' )
        # line with list of patterns was found ?
        if len(patterns) > 0:
            # loop over patterns:
            for pattern in patterns.split():
                # loop over files to be removed:
                for sfile in glob.glob(pattern) :
                    if os.path.exists(sfile) :
                        logging.debug( '    %s ...' % sfile )
                        os.remove(sfile)
                    #endif
                #endfor
            #endfor
        #endif  # found ?        
    #endfor   # loop over defined macros

    # clean for undefined macro's (read list of pattern)
    logging.info( '  remove files for undefined macro\'s ...' )
    for macr in macros_supported :
        # defined ? then next:
        if macr in macros_defined : continue
        # try to read list with files to be removed if macro is NOT defined:
        patterns = rcf.get( 'build.configure.remove.ifndef.%s' % macr, default='' )
        # line with list of patterns was found ?
        if len(patterns) > 0:
            # loop over patterns:
            for pattern in patterns.split():
                # loop over files to be removed
                for sfile in glob.glob(pattern) :
                    if os.path.exists(sfile) :
                        logging.debug( '    %s ...' % sfile )
                        os.remove(sfile)
                    #endif
                #endfor
            #endfor
        #endif  # found ?
    #endfor   # loop over supported macros


    #
    # * user script
    #

    # call user script:
    pycus.Build_Configure( rcf )


    #
    # * deps
    #
    
    # create deps ?
    flag = rcf.get('build.configure.makedep','bool')
    if flag :
        
        # name of include file with dependencies to be written:
        makefile_deps = rcf.get('build.configure.makedep.includefile')
        # info ...
        logging.info( '  create %s ...' % makefile_deps )

        # file filter for dependencies:
        files = rcf.get('build.configure.makedep.files')
        
        # explicit dependencies:
        #    target-object1 : object1 object2 ... ;
        #    target-object2 : object3 object4 ...
        explicit = rcf.get( 'build.configure.makedep.explicit' )
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
        
        # get dependency object:
        deps = MakeDep.Dependencies( flags=mk_defs, files=files, expl=expl )
        
        # no link lines yet:
        linklines = ['']
        # names of target(mainfile)s to be build:
        items = rcf.get( 'build.make.exec' )
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
                logging.error( 'items in exec list should have form :  exec [mainfile]' )
                logging.error( 'found item : %s' % item )
                raise Exception
            #endif
            
            # get list of objects for this target:
            objs = deps.GetObjects( mainfile )
            # link lines:
            linklines.append( '%s: %s.o' % (target,mainfile) )
            linklines.append( '\t$(FC) -o $@ $(FFLAGS) $(LDFLAGS) %s $(LIBS)' % objs )
            linklines.append( '' )
        #endif

        # write result:
        f = open( makefile_deps, 'w' )
        for line in deps.lines+linklines : f.write(line+'\n')
        f.close()
        # add to log file:
        logging.debug( '' )
        logging.debug( '---[%s]------------------------------------------------------------' % makefile_deps )
        for line in deps.lines+linklines : logging.debug(line)
        logging.debug( '-------------------------------------------------------------------' )
        logging.debug( '' )
        
    #endif
    
    
    #
    # * 
    #
    
    # back:
    
    logging.debug( '  change back to %s ...' % current_dir )
    os.chdir( current_dir )
    
    # ok
    return
    
#endif


# ***


def Build_Make( rcf, pycasso_user_scripts ) :

    # external:
    import os
    import logging
    
    # load user scripts:
    pycus = __import__( pycasso_user_scripts )
    
    # info ...
    logging.info( 'make executable ...' )
    
    # change directory:
    Make_dir = rcf.get('build.make.dir')
    logging.debug( '  change to %s ...' % Make_dir )
    os.chdir( Make_dir )

    # call user script:
    pycus.Build_Make( rcf )
    
#endif


# ------------------------------------------------
# tools
# ------------------------------------------------

#
# Functions copied from 'pycasso_tools.py' to ensure
# that that module is not imported before the sources
# are copied to the build directory.
#

def diff_text_files( fname1, fname2 ) :

    """
    Return True if the files are different, False if they are the same.
    """
    
    # external:
    import os
    import logging
    
    # read files:
    f = open( fname1, 'r' )
    lines1 = f.readlines()
    f.close()
    f = open( fname2, 'r' )
    lines2 = f.readlines()
    f.close()
    
    # compare:
    return lines1 != lines2
    
#enddef


# ***


def CreateDirs( dirname, forceclean=False ) :

    """
    Create a directory and report success.
    """

    # external:
    import os
    import shutil
    import logging
    
    # already present ?
    if os.path.isdir( dirname ) :
        # remove existing directory ?
        if forceclean:
            # info ...
            logging.info( 'remove existing %s ...' % dirname )
            # use the chainsaw:
            shutil.rmtree( dirname )
        #endif
    #endif
    
    # not present (anymore) ?
    if not os.path.isdir( dirname ) :
        # info ...
        logging.info( 'create new directory %s ...' % dirname )
        # create directory including parents if necesary:
        os.makedirs(dirname)
    #endif
    
    # ok
    return None

#enddef

            
# ------------------------------------------------
# end
# ------------------------------------------------

