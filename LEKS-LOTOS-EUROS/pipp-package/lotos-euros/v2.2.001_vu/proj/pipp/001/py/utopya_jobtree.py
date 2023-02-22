#! /usr/bin/env python

"""
*************************
``utopya_jobtree`` module
*************************

This module defines classes to create a sequence of job scripts.

* The :py:class:`UtopyaJobStep` class could be used to create a single element
  of a sequence. 
  It can create and submit a stan-alone job file performing some task configured
  by rcfile settings.
  In addition, methods are included to create and submit a following job,
  which are used by the :py:class:`.UtopyaJobTree` class to generate a sequence
  of jobs. 
* The :py:class:`UtopyaJobTree` class creates a sequence of job, configured
  as a tree with a main job and sub-jobs, which could have sub-jobs too, etc.
* The :py:class:`UtopyaJobIteration` class creates a job-tree with sub-jobs
  defined by an iteration number. The iteration over sub-jobs is stopped if certain
  conditions are reached.
  
The following classes are defined to perform actual tasks:

* The dummy class :py:class:`UtopyaJobTask` is used in the examples 
  to show where a user defined task should be defined.
* The task :py:class:`UtopyaJobTaskSubmit` could be used to create
  and submit a (job)tree (as defined in this module); this is useful
  to let a job submit other jobs and continue with doing other things.
  
Class hierarchy
===============

The classes are defined according to following hierarchy:

* :py:class:`.UtopyaRc`

  * :py:class:`.UtopyaJobStep`

    * :py:class:`.UtopyaJobTree`

      * :py:class:`.UtopyaJobIteration`

        * :py:class:`.UtopyaJobIteration_CheckFile`

      * JobTimeSeries
      
* :py:class:`UtopyaJobTask`

  * :py:class:`UtopyaJobTaskSubmit`

Classes
=======

"""


# ======================================================================
# ===
# === modules
# ===
# ======================================================================

# tools:
import utopya_rc


# ======================================================================
# ===
# === UtopyaJobStep
# ===
# ======================================================================


class UtopyaJobStep( utopya_rc.UtopyaRc ) :

    """
    Base class for single job step.
    
    The 'UtopyaJobStep' class and its derivatives are used to form a chain of jobs,
    where each job contains lines that create and submit the next job.
    The settings in the rcfile define for each of the jobs the class,
    where to submit the job file to (foreground or batch system),
    and a task to be performed.
    
    Because this is the base class for the :py:class:`.UtopyaJobTree` class
    that defines an actual sequence, a number of methods are included to create and
    submit a next job. Without using these, the class can be used to create a
    single stand-alone job file performing some task configured by rcfile settings.
    
    **Simple usage**
    
    Example of simple usage::
          
        # init UtopyaJobStep with name 'appl', and read settings for this name:
        jbs = UtopyaJobStep( 'appl', 'settings.rc', rcbase='', env={} )
        # write first job and submit:
        jbs.Start()
          
    A shown in the example, the following arguments are used on initialization:
    
    * name : Job step name, used to read settings and form the job file name.
    * rcfile : Name of settings file.    
    * rcbase : Optional prefix for keywords in rcfile.
    * env    : Optional dictionairy with variables that will be exported to the environment.
    
    An example of the rcfile settings needed for a job named 'appl' that should
    be submitted to an LSF queue::
      
        ! setup logging:
        *.logging.verbose         :  False
        *.logging.format          :  %(asctime)s [%(levelname)-8s] %(message)s

        ! class to create and submit this job within the jobtree:
        appl.class                     :  utopya.UtopyaJobStep

        ! class with the job script creator:
        appl.script.class              :  utopya.UtopyaJobScriptLSF
        
        ! work directory, here formed from job name "appl":
        appl.workdir                   :  ${my.work}/__NAME2PATH__

        ! (optional) add line to change to work directory,
        ! for example because job scheduler does not do that:
        appl.cwd                       :  True
        
        ! (optional) python search path:
        appl.pypath                    :  ./py
        
        ! (optional) environment modules:
        appl.modules                   :  load netcdf/4.4.1 ; load python/2.7
        
        ! (optional) extra environment variables that will be used
        !   to setup the next jobscript:
        appl.env                       :  HOMEDIR="${PWD}", WORKDIR="/work/me/appl"
        
        ! (optional) extra lines:
        appl.lines                     :  os.environ['OMP_NUM_THREADS'] = os.environ['JOB_NTHREAD']
        
        ! task:
        appl.task.class                :  utopya.UtopyaJobTask
        appl.task.args                 :  msg='Perform application task.'
    
    The rcfile should contain at least the settings required by the parent 
    class :py:class:`.UtopyRC`, in this case to setup default logging parameters.
    
    The first specific setting is the class of the job step in ``<name>.class``.
    This is a necessary setting for automatic creation of jobs by
    derived classes such as :py:class:`UtopyaJobTree` and :py:class:`UtopyaJobIteration`,
    where the creating job needs to know how to form the next job.

    The next setting is ``<name>.script.class``, which should specify the name of 
    the class that should be used to create and submit the job script. 
    This class should be derived from the :py:class:`.UtopyaJobScript` class, 
    see the :py:mod:`utopya_jobscript` module for default available classes.
    The above example settings will run the script in foreground.
    If the job should be submitted to a queue,
    use a different script class and specify the job options;
    an example for LSF queue system is ::

        ! Specify the module and class with the job script creator:
        appl.script.class              :  utopya.UtopyaJobScriptLSF

        ! batch options for "UtopyaJobScriptLSF" class:
        appl.batch.lsf.format          :  batch.lsf.format
        appl.batch.lsf.options         :  name output error
        appl.batch.lsf.option.name     :  J %(env:name)
        appl.batch.lsf.option.output   :  oo %(name).out
        appl.batch.lsf.option.error    :  eo %(name).err
    
        ! Define format of lsf options, here:
        !   #BSUB -flag value
        batch.lsf.format.comment           :  #
        batch.lsf.format.prefix            :  BSUB
        batch.lsf.format.arg               :  -
        batch.lsf.format.assign            :  ' '
        batch.lsf.format.template          :  %(key)
        batch.lsf.format.envtemplate       :  %(env:key)
        
    The python search path (for the utopya modules) could be specified
    using the ``<name>.pypath`` setting, read on initialization.
    Multiple search directories could be defined using a ``:``-seperated list.
    Note that this setting is not required and might be ommitted.

    Another optional setting is the definition of extra environment variables
    using the ``<name>.env`` setting, read on initialization. 
    The specified content will be converted to a dictionairy type, 
    and the key/value pairs will be exported to the environment.
        
    The actual work should be done by one or more objects,
    for which the class and initialization arguments are defined
    by ``<name>.<task>.*`` settings.
    For the above example, the following lines will be inserted in the job script::
        
        tskclass = utopya.ImportClass( "utopya.UtopyaJobTask" )
        tsk = tskclass(msg='Perform application task.')
        
    The arguments could include the template '%{rcfile}' to insert
    the name of the settings file::

        # use "settings.rc" from the example:
        appl.task.args      :  rcfile='%{rcfile}', rcbase='appl'
        
    Similar the templates '%{workdir}' and '%{name}' could be used to insert the 
    work directory and job name respectively;
    note that the work directory will always contain the pre-processed settings that
    were used to create the job in the file '%{name}.rc'::

        # use "task.rc" from the example:
        appl.task.args      :  rcfile='%{workdir}/%{name}.rc', rcbase='appl'
        
    The ``task`` part in the rc keys is actually a value out of a list.
    Use the following settings to define a list of three tasks::
        
        ! tasks:
        appl.tasks                      :  wakeup work sleep
        
        ! task:
        appl.wakeup.class               :  utopya.UtopyaJobTask
        appl.wakeup.args                :  msg='Wake up!'
        
        ! task:
        appl.work.class                 :  utopya.UtopyaJobTask
        appl.work.args                  :  msg='Work ...'
        
        ! task:
        appl.sleep.class                :  utopya.UtopyaJobTask
        appl.sleep.args                 :  msg='... and go to sleep'
        
    If no task list is specified, default list has just a single element named ``task``.
   
    **Job chain settings**
         
    In case of a job chain, a finished job sets up and submits the next one.
    The following optional setting are read to specify the work directory 
    where to write the next job file (default is current directory)::
    
        ! (optional) working directory:
        appl.wakeup.workdir        :  /scratch/you/appl-dir
        
    A special template ``__NAME2PATH__`` could be included in the path to insert the 
    job name (``appl.wakeup``) with the dots replaced by path-seperation characters;
    for example::

        ! working directory incl. subdirs for name elements:
        appl.wakeup.workdir        :  /scratch/you/__NAME2PATH__
        
    will be expanded to::
    
        /scratch/you/appl/wakeup
        
    The rcfile that should be used to initialize the object of the next job
    (default is the rcfile used for the finished job)::
    
        ! (optional) settings to initialize the job step:
        appl.wakeup.rcfile         :  rc/my-appl-wakeup.rc
    
    **Overview of methods**
    
    With the above settings it should be possible to define
    all application depended features of the jobs.
    If necessary to extend the class, the following is an
    overview of the underlying methods.
    
    * Start the job step (stand-alone), or the first job in a sequence (tree),
      by calling the following method:
        
      * :py:meth:`Start`

    * The job scripts are created and submitted with method:

      * :py:meth:`Run`

    * The 'Run' method uses the following methodes to fill the content
      of the job fils:

      * :py:meth:`AddHeader`
      * :py:meth:`AddOptions`
      * :py:meth:`AddModules`
      * :py:meth:`AddEnvModules`
      * :py:meth:`AddVariables`
      * :py:meth:`AddTasks`
      * :py:meth:`AddNextJob`
      * :py:meth:`AddFooter`

    * A job step might define special variables that can be passed
      in the task arguments. The following method is used to define these:

      * :py:meth:`GetVariables`

    * The job file is written and submitted using:

      * :py:meth:`GetJobFileName`
      * :py:meth:`WriteAndSubmit`

    * The following methods will be used to obtain the next jobname in the chain:

      * :py:meth:`GetNextElement`
      * :py:meth:`GetNextName`

    * In case the job chain consists of an iteration sequence,
      a method should be provided to decide on continuation or termination of the sequence:

      * :py:meth:`CheckContinuation`
        
    """
    
    def __init__( self, name, rcfile, rcbase='', env={} ) :
    
        """
        Initialize job step script from rcfile settings.
        """
        
        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )        
        # store arguments:
        self.name   = name
        self.rcfile = rcfile
        self.rcbase = rcbase
        
        # get generic name:
        self.name_generic = self._GetGenericName( self.name )
        
        # get optional path:
        pypath = self.GetSetting( self.name_generic+'.pypath', default='' )
        if len(pypath) > 0 :
            self.pypath = pypath.split(':')
        else :
            self.pypath = []
        #endif
        
        # environment to be set:
        self.env = env
        # eventually from settings:
        line = self.GetSetting( self.name_generic+'.env', default='' )
        if len(line) > 0 :
            # convert to dictionairy and store:
            self.env = eval( 'dict( %s )' % line )
        #endif
        
        # main log file for entire job chain?
        # defined by base of name only:
        name_base = self.name_generic.split('.')[0]
        # filename:
        self.logfile = self.GetSetting( name_base+'.logfile', default='' )
        # translation dictionairy for messages:
        self.logfile_translate = self.GetSetting( name_base+'.logfile.translate', 'dict', default={} )
 
        # empty content:
        self.lines = []

    #enddef __init__
    
    # *
    
    def _LogMessage( self, message, create=False ) :
    
        """
        Add message to log.
        """
        
        # modules:
        import datetime
        
        # tools:
        import utopya_tools
        
        # defined ?
        if len(self.logfile) > 0 :
        
            # translate following dictionairy:
            for key in self.logfile_translate.keys() :
                # replace key by value if present:
                message = message.replace( key, self.logfile_translate[key] )
            #endfor

            # current time:
            t = datetime.datetime.now()
            # message line:
            line = '%s: %s\n' % (t.strftime('%Y-%m-%d %H:%M:%S'),message)

            # create directory if necessary:
            utopya_tools.CreateFilePath( self.logfile )

            # file creation mode:
            if create :
                # info ...
                self.logger.info( 'create log file: %s' % self.logfile )
                # write, new file:
                cmode = 'w'
            else :
                # append:
                cmode = 'a'
            #endif
            # open:
            f = open( self.logfile, cmode )
            # write:
            f.write( line )
            # close:
            f.close()

        #endif
        
    #enddef _LogMessage
    
    # *
    
    def _SplitName( self, name ) :
        
        """
        Split names in roots and elements.
        Return values:
            root,element
        """

        # remove leading dots if present:
        name = name.lstrip('.')

        # split in root and element if possible:
        if '.' in name :
            root,element = name.rsplit('.',1)
        else :
            root,element = '',name
        #endif
        
        # ok
        return root,element

    #enddef _SplitName
    
    # *
    
    def _GetGenericName( self, name, indent='' ) :
        
        """
        Return generic version of actual name.
        
        The name is splitted at the dots into elements,
        and a generic name is build element by element.
        For each new element, the settings are scanned
        for presence of a key '<root>.generic' with '<root>'
        the generic name build up to then. If the key is
        not defined, the element is assumed to be non-generic.
        """
        
        # init result:
        name_generic = ''
        # split name:
        elements = name.split('.')
        # loop over elements:
        for element in elements :
            # first cannot be generic ..
            if len(name_generic) == 0 :
                # fill:
                name_generic = element
            else :
                # try to read generic name of this element,
                # if not defined, just use the element itself:
                element_generic = self.GetSetting( name_generic+'.generic', default=element )
                # add:
                name_generic = name_generic+'.'+element_generic
            #endif
        #endfor # elements
        
        # ok
        return name_generic
        
    #enddef _GetGenericName
    
    # *
    
    def GetGenericName( self ) :
    
        """
        Return own generic name, usefull to obtain extra settings by derived classes.
        """
        
        return self._GetGenericName( self.name )
        
    #enddef GetGenericName
    
    # *
    
    def GetVariables( self, element ) :
        
        """
        Return dictionairy with job variables for this class.
        This is used by the :py:meth:`AddVariables` method to add definition lines 
        to the job script.
        """
        
        # empty:
        var = {}
        
        # ok
        return var
        
    #enddef GetVariables
    
    # *
    
    def Append( self, line ) :
    
        """
        Add line to script content.
        A newline is added automatically.
        """

        # add line:
        self.lines.append( line+'\n' )
        
    #enddef GetFileName
    
    # *
    
    def AddHeader( self ) :
    
        """
        Add interprator line to script content, here for a python script.
        """
    
        # modules:
        import os
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # default shell:
        default = '/usr/bin/env python'
        # replace with settings if present:
        shell = self.GetSetting( name_generic+'.shell', default=default )
        
        # init script lines:
        self.Append( '#! %s' % shell )
        self.Append( '' )

    #enddef AddHeader
    
    # *
    
    def AddOptions( self ) :
    
        """
        Add batch job options to script content.
        
        If and how job options are formed is controlled by the script class,
        defined in the rcfile for each job.
        Example for a job name ``appl``::
            
            ! Specify the module and class with the job script creator:
            appl.script.class              :  utopya.JobScriptLSF

        The script class is derived from the :py:class:`.JobScript` class,
        and has a method :py:meth:`.JobScript.GetOptionsRc` to form lines 
        that can be added to the script. 
        The method reads setting from the rcfile for the
        provided generic name::

            ! batch options:
            appl.batch.lsf.format          :  lsf_format
            appl.batch.lsf.options         :  name output error workdir
            appl.batch.lsf.option.name     :  J %(env:name)
            appl.batch.lsf.option.output   :  oo %(name).out
            appl.batch.lsf.option.error    :  eo %(name).err
            appl.batch.lsf.option.workdir  :  cwd %(env:cwd)

        An enviroment is passed with pre-defined values that can be
        subsituted in the option values, here for the job actual name and the
        current working dirctory::
            
            env = { 'name' : 'step0001',
                    'cwd'  : '/scratch/me/test' }
                  
        With the following 'lsf_format' option formatting::
            
            ! Define format of lsf options, here:
            !   #BSUB -flag value
            lsf_format.comment           :  #
            lsf_format.prefix            :  BSUB
            lsf_format.arg               :  -
            lsf_format.assign            :  ' '
            lsf_format.template          :  %(key)
            lsf_format.envtemplate       :  %(env:key)

        this will lead to the following option lines::
            
            #BSUB -J step0001
            #BSUB -oo step0001.out
            #BSUB -eo step0001.err
            #BSUB -cwd /scratch/me/test              
                    
        """
    
        # modules:
        import os
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # import a module class derived from 'JobScript' :
        jbs_cls = self.ImportClass( name_generic+'.script.class' )        
        # job script object:
        jbs = jbs_cls()

        # environment for job options:
        env = {}
        env['name'] = self.name
        env['cwd' ] = os.getcwd()
        # form line with job options from rcfile settings:
        options = jbs.GetOptionsRc( self.rcfile, rcbase=self.rcbase+'.'+name_generic, env=env )
        # add:
        self.Append( options )
        self.Append( '' )
        
    #enddef AddOptions

    # *
    
    def AddModules( self ) :
    
        """
        Add lines to import standard modules (os, sys)
        and tool modules (utopya).
        """
        
        # info:
        self.Append( '#' )
        self.Append( '# *** modules' )
        self.Append( '#' )
        self.Append( '' )

        # default modules:
        self.Append( '# modules:' )
        self.Append( 'import os' )
        
        # extend path if necessary:
        if len(self.pypath) > 0 :
            self.Append( 'import sys' )
            self.Append( '' )
            self.Append( '# extend path:' )
            for i in range(len(self.pypath)) :
                self.Append( 'sys.path.insert( %i, "%s" )' % (i,self.pypath[i]) )
            #endfor
        #endif

        # add lines to import tools:
        self.Append( '' )
        self.Append( '# tools:' )
        self.Append( 'import utopya' )
        self.Append( '' )
        self.Append( '' )
        
    #enddef AddModules
    
    # *
    
    def AddEnvModules( self ) :
    
        """
        Add script lines for GNU Environment modules.
        On many computing platforms, the environment for running applications
        is managed using 'module' commands, e.g.::
        
            module load netcdf
            module load python
            
        The module commands to be performed are optionally defined in the
        recfile for the current job as a semi-colon seperated list::
            
            <name>.modules        :  load netcdf ; load python
            
        Appropriate job lines for these settings will be inserted in the script.
        
        The location of the GNU module scripts should be available in the environment::
        
            MODULESHOME=/opt/modules/3.2.10.4
            
        """
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # module commands defined?
        cline = self.GetSetting( name_generic+'.modules', default='' )
        # any?
        if len(cline) > 0 :
        
            # info:
            self.Append( '#' )
            self.Append( '# *** environment modules' )
            self.Append( '#' )
            self.Append( '' )
            self.Append( '# check ...' )
            self.Append( 'if "MODULESHOME" not in os.environ.keys() :' )
            self.Append( '    print( "ERROR - could not load environment modules without MODULESHOME" )' )
            self.Append( '    raise Exception' )
            self.Append( '#endif' )
            self.Append( '# initialization scripts for environment modules:' )
            self.Append( 'initfiles = []' )
            self.Append( 'initfiles.append( os.path.join( os.environ["MODULESHOME"], "init", "python.py" ) )' )
            self.Append( 'initfiles.append( os.path.join( os.environ["MODULESHOME"], "init", "python" ) )' )
            self.Append( 'initfiles.append( os.path.join( os.environ["MODULESHOME"], "init", "env_modules_python.py" ) )' )
            self.Append( '# search first one that is present:' )
            self.Append( 'found = False' )
            self.Append( 'for initfile in initfiles :' )
            self.Append( '    found = os.path.isfile(initfile)' )
            self.Append( '    if found : break' )
            self.Append( '#endfor' )
            self.Append( 'if not found :' )
            self.Append( '    print( "ERROR - could not find any environment modules initialization file: %s" % initfile )' )
            self.Append( '    for initfile in initfiles :' )
            self.Append( '        print( "ERROR -   %s" % initfile )' )
            self.Append( '    #endfor' )
            self.Append( '    raise Exception' )
            self.Append( '#endif' )
            self.Append( '# load module functions:' )
            self.Append( 'exec( open(initfile).read() )' )
            self.Append( '' )
            
            # add info:
            self.Append( '# module commands:' )     
            # split into individual commands:
            commands = cline.split(';')
            # loop:
            for command in commands :
                # start new line:
                line = 'module( '
                # split into elements:
                elements = command.split()
                # loop over elements:
                for k in range(len(elements)) :
                    # current:
                    element = elements[k]
                    # add seperation?
                    if k > 0 : line = line+', '
                    # add quoted version:
                    line = line+("'%s'" % element)
                #endfor # elements
                # close:
                line = line+' )'
                # add command:
                self.Append( line )
            #endfor  # commands
            
            # add lines to list currently loaded modules:
            self.Append( '' )
            self.Append( '# show current modules:' )
            self.Append( "module( 'list' )" )

            # add lines to update pythonpath:
            self.Append( '' )
            self.Append( '# modules might have extended the PYTHONPATH environment variable,' )
            self.Append( '# but these are not automatically added to the python search path;' )
            self.Append( '# check enviornment:' )
            self.Append( 'if "PYTHONPATH" in os.environ.keys() :' )
            self.Append( '    # extract as list:' )
            self.Append( '    ppaths = os.environ["PYTHONPATH"].split(":")' )
            self.Append( '    # reverse order to have correct order for inserting:' )
            self.Append( '    ppaths.reverse()' )
            self.Append( '    # add to search path' )
            self.Append( '    for ppath in ppaths :' )
            self.Append( '        # insert at start if not present yet:' )
            self.Append( '        if ppath not in sys.path : sys.path.insert(0,ppath)' )
            self.Append( '    #endfor # paths' )
            self.Append( '#endif  # PYTHONPATH defined' )
            self.Append( '' )

        #endif # module commands needed

    #enddef AddEnvModules
    
    # *
    
    def AddLines( self ) :
    
        """
        Add user defined script lines.
        Could be used to setup special environment.

        For example, if the job scheduler defines an environment variable
        'JOB_NTHREAD' for the number of OpenMP threads,
        the jobscript could use this to define the 'OMP_NUM_THREADS'
        variable needed by OpenMP code.
        In a python script, this looks like::
        
            os.environ['OMP_NUM_THREADS'] = os.environ['JOB_NTHREAD']
            
        Specify this code in the settings::
            
            <name>.lines    :  os.environ['OMP_NUM_THREADS'] = os.environ['JOB_NTHREAD']
        
        For a complete block of code, use '\\\\n' marks for line breaks, '\\\\t' for
        leading indents, and a multi-line rc value for a better readible definition.
        For example, the following definition::
        
            <name>.lines    : \\n\\
                print( 'environment:' )\\n\\
                for key in os.environ.keys() :\\n\\
                \\tprint( '  %=%' % (key,os.environ[key]) )
            
        will be expanded to::
        
            print( 'environment:' )
            for key in os.environ.keys() :
                print( '  %=%' % (key,os.environ[key]) )

        """
        
        # modules:
        import os
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # read flag, default False:
        lines = self.GetSetting( name_generic+'.lines', default='' )
        # enabled?
        if len(lines) > 0 :
        
            # info:
            self.Append( '#' )
            self.Append( '# *** user defined lines' )
            self.Append( '#' )
            self.Append( '' )
            for line in lines.split('\\n') :
                self.Append( line.strip().replace('\\t','    ') )
            #endfor    
            self.Append( '' )

        #endif # add user lines

    #enddef AddLines
    
    # *
    
    def _GetWorkdir( self ) :
    
        """
        Return work directory as specified by '<name>.workdir' in rcfile,
        or otherwise current directory.
        """
        
        # modules:
        import os
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )

        # work directory:
        wdir = self.GetSetting( name_generic+'.workdir', default='' )
        # workdir specified ?
        if len(wdir) > 0 :
            # special: if the special template is in the work dir,
            # insert the next name with '.' replaced by '/':
            wdir = wdir.replace( '__NAME2PATH__', self.name.replace('.','/') )
        else :
            # current:
            wdir = os.getcwd()
        #endif
        
        # ok
        return wdir
        
    #enddef _GetWorkdir

    # *
    
    def AddCwd( self ) :
    
        """
        Add script lines to change to work directory.
        This is necessary in case the job schedular does not change to the direcotry
        of the job script, or if the job options do not have a flag to specify it.
        These lines are only included if the following flag is set::
            
            <name>.cwd        :  True
            
        """
        
        # modules:
        import os
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # read flag, default False:
        enabled = self.GetSetting( name_generic+'.cwd', totype='bool', default=False )
        # enabled?
        if enabled :

            # work directory:
            wdir = self._GetWorkdir()

            # info:
            self.Append( '#' )
            self.Append( '# *** change workdir' )
            self.Append( '#' )
            self.Append( '' )            
            self.Append( '# change to workdir:' )
            self.Append( 'os.chdir( "%s" )' % wdir )
            self.Append( '' )

        #endif # change workdir

    #enddef AddCwd

    # *
    
    def AddVariables( self ) :
    
        """
        Add commands to set job variables that might be used by the tasks.
        
        Always a variable 'name' is defined with the job name.
        
        Extra variables might be defined by specific classes,
        for example an iteration class might define the iteration step.
        These are collected in a dictionairy named "env" with keys
        formed from the previous jobs in the tree.
        An iteration job named 'appl.run' stores for example::
            
            env["appl.run.__step__"] = 4

        """
        
        # add current values:
        self.Append( '#' )
        self.Append( '# *** variables' )
        self.Append( '#' )
        self.Append( '' )
        self.Append( '# current:' )
        self.Append( 'name = "%s"' % self.name )
        self.Append( '' )
        self.Append( '# job variables:' )
        self.Append( 'env = {}' )
        
        # keep list of keys defined:
        envkeys = []
        
        # walk through job tree up to current name ;
        # for each branch, define variables if necessary:
        root = ''
        for element in self.name.split('.') :
            # root defined ?
            if len(root) > 0 :
                # generic version:
                root_generic = self._GetGenericName( root )
                # read classname for root, and return as a class object:
                jbsclass = self.ImportClass( root_generic+'.class' )
                # create object for this root:
                jbs = jbsclass( root, self.rcfile, rcbase=self.rcbase )
                # get dictionairy with job variables:
                vars = jbs.GetVariables(element)
                # loop:
                for key in vars.keys() :
                    # current:
                    val = vars[key]
                    # combined key:
                    envkey = '%s.%s' % (root_generic,key)
                    # not defined yet?
                    if envkey not in envkeys :
                        # add definition:
                        if type(val) == str :
                            self.Append( 'env["%s"] = "%s"' % (envkey,val) )
                        else :
                            self.Append( 'env["%s"] = %s' % (envkey,val) )
                        #endif
                        # add key:
                        envkeys.append( envkey )
                    #endif
                #endfor
            #endif
            # extend:
            root = (root+'.'+element).lstrip('.')
        #endfor # elements
        
        # extra enviroment variables:
        if len(self.env) > 0 :
            # loop:
            for envkey in self.env.keys() :
                # not defined yet?
                if envkey not in envkeys :
                    # extract value:
                    val = self.env[envkey]
                    # add definition:
                    if type(val) == str :
                        self.Append( 'env["%s"] = "%s"' % (envkey,val) )
                    else :
                        self.Append( 'env["%s"] = %s' % (envkey,val) )
                    #endif
                    # add key:
                    envkeys.append( envkey )
                #endif
            #endfor
        #endif
        
        # finish:
        self.Append( '' )
        self.Append( '' )
        
    #enddef AddVariables

    # *
    
    def AddTasks( self ) :
    
        """
        Add task command to script lines, that consists of a class import
        and initialization of an object of this class.
        The class name and the arguments for the initialization are
        defined in the rcfile settings::
            
            <name>.task.class     :  mymod.MyTask
            <name>.task.args      :  msg='Do something'
            
        This will insert the following lines in the job script::
            
            tskclass = utopya.ImportClass( "mymod.MyTask" )
            tsk = tskclass(msg='Do something')
        
        The arguments could include the templates '%{workdir}', '%{name}' and '%{rcfile}'
        to insert the workdirectory, name of the job, or the settings file::

            appl.task.args      :  rcfile='%{rcfile}', rcbase='applx'
        
        If the class name is left empty, nothing is inserted
        and no arguments need to be specified.
        
        Instead of a single ``task``, it is also possible to provide
        a list of tasks::

            ! tasks:
            appl.tasks                      :  wakeup work sleep

            ! task:
            appl.task.wakeup.class          :  utopya.UtopyaJobTask
            appl.task.wakeup.args           :  msg='Wake up!'

            ! task:
            appl.task.work.class            :  utopya.UtopyaJobTask
            appl.task.work.args             :  msg='Work ...'

            ! task:
            appl.task.sleep.class           :  utopya.UtopyaJobTask
            appl.task.sleep.args            :  msg='... and go to sleep'

        """
        
        # modules:
        import os
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # virtual?
        virtual = self.GetSetting( name_generic+'.virtual', totype='bool', default=False )
        
        # only tasks if not virtual ...
        if not virtual :
        
            # info:
            self.Append( '#' )
            self.Append( '# *** task(s)' )
            self.Append( '#' )
            self.Append( '' )

            # list with task names:
            tasknames = self.GetSetting( name_generic+'.tasks', default='None' ).split()
            # loop:
            for taskname in tasknames :

                # taskname part in rcfile keys:
                if taskname == 'None' :
                    taskkey = 'task'
                else :
                    taskkey = 'task.%s' % taskname
                #endif

                # work directory:
                wdir = self._GetWorkdir()

                # module and class name as defined in rcfile:
                cls = self.GetSetting( name_generic+('.%s.class' % taskkey) )

                # arguments for initialization:
                args = self.GetSetting( name_generic+('.%s.args' % taskkey) )
                # replace templates:
                args = args.replace( '%{workdir}', wdir )
                args = args.replace( '%{name}', self.name )
                args = args.replace( '%{rcfile}', self.rcfile )

                # task lines:
                self.Append( '# task class:' )
                self.Append( 'tskclass = utopya.ImportClass( "%s" )' % cls )
                self.Append( '# create task object and initialize, which does the actual work:' )
                self.Append( 'tsk = tskclass( %s )' % args )
                self.Append( '' )

            #endfor # tasks

            # add space:
            self.Append( '' )
            
        #endif # not virtual

    #enddef AddTasks
    
    # *
    
    def AddNextJob( self ) :
    
        """
        Add lines that create and submit the next job step if necessary.
        """
        
        # modules:
        import os
        
        # name of next job:
        nextname,nexttest = self.GetNextName()
        # next job to be created ?
        if nextname != 'None' :

            # get generic name:
            nextname_generic = self._GetGenericName( nextname )

            # module and class name as defined in rcfile:
            cls = self.GetSetting( nextname_generic+'.class' )
            # rcfile to be used for initialization:
            cls_rcfile = self.GetSetting( nextname_generic+'.rcfile', default=self.rcfile )
            
            # info ...
            self.Append( '#' )
            self.Append( '# *** next job' )
            self.Append( '#' )
            self.Append( '' )
            
            # should test be performed ?
            if nexttest == 'None' :

                # add lines to change to workdir:
                self._AddWorkdirLines( nextname_generic, nextname, '' )
                # add lines:
                self.Append( '# job step class:' )
                self.Append( 'jbsclass = utopya.ImportClass( "%s" )' % cls )
                self.Append( '# create job step object:' )
                self.Append( 'jbs = jbsclass( "%s", "%s", rcbase="%s", env=env )' % (nextname,cls_rcfile,self.rcbase) )
                self.Append( '# write job file and start:' )
                self.Append( 'jbs.Run()' )
                self.Append( '' )

            else :
                
                # if test should be done for 'run.iter0001', 
                # then settings are defined for root 'run' ; get this root:
                testroot,testelement = self._SplitName( nexttest )
                
                # get generic version:
                testroot_generic = self._GetGenericName( testroot )
                
                # class name for the level that should perform the check:
                testcls = self.GetSetting( testroot_generic+'.class' )
                
                # next name if iteration is finished:
                finishname,finishtest = self.GetNextName(finish=True)
                # defined ?
                if finishname != 'None' :
                    # generic version:
                    finishname_generic = self._GetGenericName( finishname )
                    # corresponding class:
                    finishcls = self.GetSetting( finishname_generic+'.class' )
                #endif
            
                # add lines to import continuation check class and create object:
                self.Append( '# UtopyaJobTree class that should perform the check:' )
                self.Append( 'jbsclass = utopya.ImportClass( "%s" )' % testcls )
                self.Append( '# create check object:' )
                self.Append( 'jbs = jbsclass( "%s", "%s", rcbase="%s" )' % (testroot,self.rcfile,self.rcbase) )
                self.Append( '# evaluate check given latest element;' )
                self.Append( '# return a character instruction key and an informative message:' )
                self.Append( 'order,msg = jbs.CheckContinuation( "%s" )' % testelement )
                self.Append( '' )
                self.Append( '# info ...' )
                self.Append( 'jbs.logger.info( "continuation order : %s" % order )' )
                self.Append( 'jbs.logger.info( "           message : %s" % msg )' )
                self.Append( '# next job ?' )
                self.Append( 'if order == "continue" :' )
                self.Append( '' )
                # add lines to change to workdir:
                self._AddWorkdirLines( nextname_generic, nextname, '    ' )
                # add job lines:
                self.Append( '    # job step class:' )
                self.Append( '    jbsclass = utopya.ImportClass( "%s" )' % cls )
                self.Append( '    # create job step object:' )
                self.Append( '    jbs = jbsclass( "%s", "%s", rcbase="%s", env=env )' % (nextname,self.rcfile,self.rcbase) )
                self.Append( '    # write job file and start:' )
                self.Append( '    jbs.Run()' )
                self.Append( '' )
                self.Append( 'elif order == "finish" :' )
                self.Append( '' )
                if finishname == 'None' :
                    self.Append( '    # iteration finished, no next job defined ... ' )
                    self.Append( '    pass' )
                else :
                    # add lines to change to workdir:
                    self._AddWorkdirLines( finishname_generic, finishname, '    ' )
                    # add job lines:
                    self.Append( '    # job step class:' )
                    self.Append( '    jbsclass = utopya.ImportClass( "%s" )' % finishcls )
                    self.Append( '    # create job step object:' )
                    self.Append( '    jbs = jbsclass( "%s", "%s", rcbase="%s", env=env )' % (finishname,self.rcfile,self.rcbase) )
                    self.Append( '    # write job file and start:' )
                    self.Append( '    jbs.Run()' )
                #endif
                self.Append( '' )
                self.Append( 'else :' )
                self.Append( '' )
                self.Append( '    # unknown ...' )
                self.Append( '    jbs.logger.error( "unsupported order \'%s\'" % order )' )
                self.Append( '    raise ValueError' )
                self.Append( '' )
                self.Append( '#endif' )

            #endif

            # finish:
            self.Append( '' )
            self.Append( '' )

        #endif  # next job defined
        
    #enddef AddNextJob
    
    # *
    
    def _AddWorkdirLines( self, nextname_generic, nextname, indent ) :
    
        """
        Add lines to job script to create and change to a work directory
        before creating the next job.
        Argument 'indent' should include the white space at the start of the lines.
        """
        
        # workdir:
        wdir = self.GetSetting( nextname_generic+'.workdir', default='' )
        # workdir specified ?
        if len(wdir) > 0 :
            # special: if the special template is in the work dir,
            # insert the next name with '.' replaced by '/':
            wdir = wdir.replace( '__NAME2PATH__', nextname.replace('.','/') )
            # add lines:
            self.Append( indent+'# workdir:' )
            self.Append( indent+'wdir = "%s"' % wdir )
            self.Append( indent+'# create if necessary:' )
            self.Append( indent+'if not os.path.isdir( wdir ) : os.makedirs( wdir )' )
            self.Append( indent+'# goto:' )
            self.Append( indent+'os.chdir( wdir )' )
            self.Append( '' )
        #endif
        
    #enddef _AddWorkdirLines

    # *
    
    def AddFooter( self ) :
    
        """
        Add closing lines to script content.
        """
        
        # end script:
        self.Append( '#' )
        self.Append( '# *** end' )
        self.Append( '#' )
        
    #enddef AddFooter
    
    # *
    
    def GetFileName( self, name ) :
    
        """
        Return name of job file to be written.
        Here use the job name and add extension '.jb'.
        """
        
        # modules:
        import os

        # target file:
        #jbfile = '%s.jb' % name.replace('.','_')
        jbfile = '%s.jb' % name
        
        # ok
        return jbfile
        
    #enddef GetFileName
        
    # *
    
    def WriteAndSubmit( self ) :
    
        """
        Write script content to file with provided name, and submit the created file.
        Creation and submission is performed by an object derived from
        the :py:class:`.JobScript` class. The name of this class is defined 
        in the settings, as well as the working directory (empty for current)::
            
            <name>.script.class         :  utopya.JobScriptForeground
            <name>.workdir              :  /work/appl/run
        """
        
        # modules:
        import os
            
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # import a module class derived from 'JobScript' :
        jbs_cls = self.ImportClass( name_generic+'.script.class' )
        # job script object:
        jbs = jbs_cls()

        # workdir:
        wdir = self.GetSetting( name_generic+'.workdir', default='' )
        # defined?
        if len(wdir) > 0 :
            # special: if the special template is in the work dir,
            # insert the next name with '.' replaced by '/':
            wdir = wdir.replace( '__NAME2PATH__', self.name.replace('.','/') )
            # create if necessary:
            if not os.path.isdir(wdir) : os.makedirs( wdir )
            # info ...
            self.logger.info( '  change to "%s" ..' % wdir )
            # change:
            os.chdir( wdir )
        #endif
        
        ## target name for pre-processed settings:
        #rcfile = '%s.rc' % self.name
        ## info ...
        #self.logger.info( '  write preprocessed settings to %s ...' % rcfile )
        ## write:
        #self.WriteSettings( rcfile )

        # target file:
        filename = self.GetFileName( self.name )
        # info ...
        self.logger.info( '  write "%s" ...' % filename )
        # write:
        jbs.Create( filename, self.lines )

        # full path for info message:
        fullname = os.path.join( os.getcwd(), filename )
        # add message to logfile;
        self._LogMessage( 'created %s' % fullname )
    
        # submit:
        jbs.Submit( filename )
        
    #enddef WriteAndSubmit
    
    # *
    
    def Run( self ) :
    
        """
        Create and submit the job for the named item in the job step.
        
        The content of the job file is filled using calls to methods:

        * :py:meth:`AddHeader`
        * :py:meth:`AddOptions`
        * :py:meth:`AddModules`
        * :py:meth:`AddEnvModules`
        * :py:meth:`AddLines`
        * :py:meth:`AddCwd`
        * :py:meth:`AddVariables`
        * :py:meth:`AddTasks`
        * :py:meth:`AddNextJob`
        * :py:meth:`AddFooter`
            
        The file is written and submitted by a call to the method::
            
        * :py:meth:`WriteAndSubmit`
            
        """
        
        # info ..
        self.logger.info( 'create job "%s" ...' % self.name )

        # init lines:
        self.AddHeader()
        # add batch job options:
        self.AddOptions()
        
        # import tools etc:
        self.AddModules()

        # load environment modules:
        self.AddEnvModules()
        
        # insert user-defined lines:
        self.AddLines()
        
        # change to workdir:
        self.AddCwd()

        # define job variables:
        self.AddVariables()
        # add task lines:
        self.AddTasks()
        
        # add lines to create and submit next job:
        self.AddNextJob()

        # finish lines:
        self.AddFooter()

        # write, and submit:
        self.WriteAndSubmit()
        
    #enddef Run
    
    # *
    
    def Start( self ) :
        
        """
        Create object for the first job along tree that is not marked as virtual,
        and call its :py:meth:`Run` method.
        """
        
        # get first non-virtual name, might be current:
        xname = self.GetFirstName( self.name )
        # get generic name:
        xname_generic = self._GetGenericName( xname )
        # read classname for root, and return as a class object:
        jbsclass = self.ImportClass( xname_generic+'.class' )
        # create object for this root:
        jbs = jbsclass( xname, self.rcfile, rcbase=self.rcbase )
        # run:
        jbs.Run()
        
    #enddef Start
    
    # *
    
    def GetNextElement( self, element, indent='' ) :
    
        """
        Returns the next element in a list (iteration?) of job steps.
        For derived classes such as 'UtopyaJobTree' it is sufficient to re-define only
        this method.
        
        If no next sub-element is available, the value 'None' is returned.
        If the requested element is 'None', the name of the first sub-element
        is returned (if present).
        """
        
        # no sub-elements:
        nextelement = 'None'
        
        # ok
        return nextelement
        
    #enddef GetNextElement

    # *
    
    def GetNextName( self, finish=False, indent='' ) :
    
        """
        Return information on the next job name in a chain, including information 
        on performing a test on continuation of the chain (if necesary).
        
        Two values are returned: the next job name,
        and the name of the job step that should decide on continuation
        (or 'None' if not needed).

        The name of the next job is read from the job tree definition in the rcfile.
        The next job could also be specified directly using an adhoc setting,
        which is useful to skip a part of the tree::

           <name>.jump-to     :  nextname

        """

        ## testing ...
        #print indent+'xxx UtopyaJobStep/GetNextName: for name "%s"' % (self.name)
        
        # ADHOC: optionally jump from one job to the other defined by:
        #   <name>.jump-to      :  nextname
        # Use destination 'None' to terminate job list.
        # read destination (if defined):
        jump_to = self.GetSetting( self.name+'.jump-to', default='' )
        # defined?
        if len(jump_to) > 0 :
            # info ...
            self.logger.warning( '  JUMP TO NEXT JOB "%s" ...' % jump_to )
            # return as next name, no continuation test:
            return jump_to,'None'
        #endif

        # get first sub-element if present, or 'None' otherwise:
        nextelement = self.GetNextElement( 'None', indent=indent+'  ' )
        ctest = self._WithContinuationTest()
        ## testing ...
        #print indent+'  x UtopyaJobStep/GetNextName: next element "%s", test %s' % (nextelement,ctest)                

        # no next branch ? then one level back:
        if nextelement == 'None' :
            ## testing ...
            #print indent+'  x UtopyaJobStep/GetNextName: no next branch, trace back ...'
            # current name:
            xname = self.name
            # current finish flag:
            xfinish = finish
            # loop unit next name is found:
            while True :
                ## testing ...
                #print indent+'  x UtopyaJobStep/GetNextName: try for xname "%s" ...' % xname
                # split current name in root and element:
                root,element = self._SplitName( xname )
                # already at lowest level ?
                if len(root) == 0 :
                    # end of chain:
                    nextname = 'None'
                    # no continuation test needed either:
                    nexttest = 'None'
                    ## testing ...
                    #print indent+'  x UtopyaJobStep/GetNextName: end of chain, return next name "%s", test "%s"' % (nextname,nexttest)
                    # quit loop:
                    break
                else :
                    # get generic name:
                    root_generic = self._GetGenericName( root )
                    # read classname for root, and return as a class object:
                    jbsclass = self.ImportClass( root_generic+'.class' )
                    # create object for this root:
                    jbs = jbsclass( root, self.rcfile, rcbase=self.rcbase )
                    # get next sub-element if present, or 'None' otherwise:
                    nextelement = jbs.GetNextElement( element, indent=indent+'  ' )
                    ctest = jbs._WithContinuationTest()
                    # test needed and requested next name after test indicates finish?
                    if ctest and xfinish :
                        # reset result to no next sub-elment:
                        nextelement = 'None'
                        # reset finish flag for higher levels:
                        xfinish = False
                        ## testing ...
                        #print indent+'  x UtopyaJobStep/GetNextName: reset nextelement to "%s"' % nextelement
                    #endif
                    # defined ?
                    if nextelement != 'None' :
                        # combine:
                        nextname = root+'.'+nextelement
                        # need to perform a continuation test ?
                        if ctest :
                            nexttest = xname
                        else :
                            nexttest = 'None'
                        #endif
                        ## testing ...
                        #print indent+'  x UtopyaJobStep/GetNextName: next name "%s", test "%s"' % (nextname,nexttest)
                        # quit loop:
                        break
                    #endif
                #endif
                # try again with shorter name:
                xname = root
            #endwhile # trace back loop
        else :
            # combine:
            nextname = self.name+'.'+nextelement
            # no testing ..
            nexttest = 'None'
            ## testing ...
            #print indent+'  x UtopyaJobStep/GetNextName: next name "%s", test %s' % (nextname,nexttest)                
        #endif

        # eventually jump to the first tip of a branch:
        nextname = self.GetFirstName( nextname, indent=indent+'  ' )

        # ok
        return nextname,nexttest

    #enddef GetNextName
    
    # *
    
    def GetFirstName( self, name, indent='' ) :
        
        """
        Return first name along tree that is not marked as "virtual".
        """
        
        # init name along tree:
        xname = name
        # loop:
        while True :
            # end ?
            if xname == 'None' : break
            # get generic name:
            xname_generic = self._GetGenericName( xname, indent='  ' )
            # virtual job ? then continue searching:
            virtual = self.GetSetting( xname_generic+'.virtual', totype='bool', default=False )
            if not virtual : break
            # read classname for root, and return as a class object:
            jbsclass = self.ImportClass( xname_generic+'.class' )
            # create object for this root:
            jbs = jbsclass( xname, self.rcfile, rcbase=self.rcbase )
            # get first sub-element if present, or 'None' otherwise:
            nextelement = jbs.GetNextElement( 'None', indent=indent+'  ' )
            # leave if tip of branch:
            if nextelement == 'None' : break
            # extend name and walk on:
            xname = xname+'.'+nextelement
        #endwhile
        ## testing ...
        #print indent+'  x UtopyaJobStep/GetNextName: name "%s" nextname "%s" xname "%s"' % (self.name,nextname,xname)                
        
        # ok:
        return xname
        
    #enddef GetFirstName
        
    # *
    
    def _WithContinuationTest( self ) :
        
        """
        Return bool flag to enable a continuation test.
        The flag is usually False, but set to True if the element is part of a sequence 
        from which continuation has to be deceided; see the 'CheckContinuation' method.
        """
        
        # by default no test:
        return False
        
    #enddef _WithContinuationTest
    
    # *
    
    def CheckContinuation( self, element ) :
    
        """
        Check status given job name and its generic represenation.
        The job name could be used to derive the iteration step number
        and to read output files; the generic name is usefull to read settings.
        
        Returns two str values:
         * order  : one of 'continue', 'finish', 'error'
         * msg    : informative message to explain the order
        
        Here always return 'finish' since this is a single job.
        """
        
        # ok
        return 'finish','no job chain, single step only'
        
    #enddef CheckContinuation
    
#endclass UtopyaJobStep


# ======================================================================
# ===
# === UtopyaJobTree
# ===
# ======================================================================


class UtopyaJobTree( UtopyaJobStep ) :

    """
    Class to create and submit a chain of job scripts
    defined in the rcfile as a tree with branches.
    
    Example of a tree defined for name 'appl'::
    
       appl
           .build
           .init
                .emis
                .obs
                    .point
                    .sat
                    .valid
           .run
               .fwd
               .dep
               .grd
               .opt
           .done
              
    This will create the job chain (see below for the meaning of "v")::
    
        appl.jb                   (v)
        appl.build.jb
        appl.init.jb              (v)
        appl.init.emis.jb
        appl.init.obs.jb          (v)
        appl.init.obs.point.jb
        appl.init.obs.sat.jb
        appl.init.obs.valid.jb
        appl.run.jb               (v)
        appl.run.fwd.jb
        appl.run.dep.jb
        appl.run.grd.jb
        appl.run.opt.jb
        appl.done.jb
        
    The tree is defined by lists of element names.
    The definitions in the settings should first define the
    elements of the main trunk.
    Example for settings for a job named 'appl'::

        ! class to create a job tree:
        appl.class                     :  utopya.UtopyaJobTree
        
        ! list of sub-elements:
        appl.elements                  :  build init run done
        
    .. _virtual-jobs:
        
    A job is also created for the trunk, in this case "appl".
    This is useful in case resources (memory, cpu's) should
    be allocated once for all jobs in a sub-tree;
    to achieve this, define resources for the trunk,
    and submit the elements to the foreground.
    The "trunk" job is skipped if it is declared to be "virtual";
    in the above example, the jobs that can be skipped in this way
    are marked with a "v".
    A trunk is declared virtual by an optional rc setting,
    which is False by default if not defined::
        
        ! virtual main job?
        appl.virtual                   :  True
        
    For each element in the list it is necessary to define the
    class that should be used to create it.
    For the "appl.build" job, this could be simply the "UtopyaJobStep" class
    since no sub-jobs are necessary::

        ! job step class for this branch:
        appl.build.class                     :  utopya.UtopyaJobStep
        
    For the "appl.init" job however, sub-jobs are defined.
    Use for this the "UtopyaJobTree" class again, and define a list with the sub elements::

        ! UtopyaJobTree class for this branch:
        appl.init.class                      :  utopya.UtopyaJobTree
        appl.init.elements                   :  emis obs
        
    The elements are combined with the 'parent' elements and form
    together the full job name, for example "appl.init.emis".

    For all names (non-virtual) in the tree, define the class that should be used 
    to create and submit individual jobs.
    If the jobs are to be submitted to a queue, specify job options too.
    Example for the 'appl.build' step::

        ! submit to LSF queue:
        appl.build.script.class                 :  utopya.JobScriptLSF

        ! batch options:
        appl.build.batch.lsf.format             :  batch.lsf
        appl.build.batch.lsf.options            :  name output error
        appl.build.batch.lsf.option.name        :  J %(env:name)
        appl.build.batch.lsf.option.output      :  oo %(name).out
        appl.build.batch.lsf.option.error       :  eo %(name).err

    The actual work is again performed by an object derived from
    the :py:class:`.UtopyaJobTask` class, for which proper
    initialization arguments should be specified::
            
        appl.build.task.class     :  mymod.MyTask
        appl.build.task.args      :  msg='Do something'

    While testing the job tree it is sometimes useful to skip a number
    of sub-jobs.
    This could be specified by a 'jump-to' specification.
    If this is present for a certain job name, the value should
    be the name of the next job that should be created.
    For example, to have build new scripts and executables
    but skip the initialization steps, include in the settings::
    
        appl.build.jump-to      :  appl.run

    """

    # *
    
    def GetNextElement( self, element, indent='' ) :
    
        """
        Returns next element from the list to which the requested element belongs.
        For example, for the name "init.obs" the list of elements is defined by::
        
            init.obs.elements   :  point sat valid
            
        In this example, the next element after "point" is "sat".
        
        If the requested element is 'None', return the first element ("point").
        Otherwise, the requested element should be in the defined list,
        and either the next element is returned, or 'None' if the
        last element was requested.
        """
        
        ## testing ...
        #print indent+'yyy UtopyaJobTree/GetNextElement of "%s" after "%s"' % (self.name,element)

        # get generic name:
        name_generic = self._GetGenericName( self.name )

        # read list:
        listkey = name_generic+'.elements'
        elements = self.GetSetting( listkey ).split()
        # no name yet? then first element of list should be returned:
        if element == 'None' :
            # empty branch ...
            if len(elements) == 0 :
                self.logger( 'empty list defined by "%s"' % listkey )
                raise Exception
            else :
                nextelement = elements[0]
                ## testing ..
                #print indent+'  y requsted empty name, returned first element "%s"' % nextelement
            #endif
        else :
            # check ..
            if element not in elements :
                self.logger.error( 'element "%s" not in list "%s" defined by rc key "%s"' % (element,elements,listkey) )
                raise Exception
            #endif
            # index:
            ind = elements.index(element)
            # last ?
            if ind == len(elements)-1 :
                # no next one ...
                nextelement = 'None'
                ## testing ...
                #print indent+'  y end of list, return "%s"' % nextelement
            else :
                # next:
                nextelement = elements[ind+1]
                ## testing ...
                #print indent+'  y selected next in list : "%s"' % nextelement
            #endif
        #endif
        
        # ok
        return nextelement

    #enddef GetNextElement

#endclass UtopyaJobTree


# ======================================================================
# ===
# === UtopyaJobIteration
# ===
# ======================================================================


class UtopyaJobIteration( UtopyaJobTree ) :

    """
    Class to create and submit a chain of job scripts
    that are defined as iteration steps.
    
    Example of a tree defined for name 'appl'::
    
       appl
           .build
           .init
                .emis
                .obs
                    .point
                    .sat
                    .valid
           .run
               .iter-0001
                         .fwd
                         .dep
                         .grd
                         .opt
           .run
               .iter-0002
                         .fwd
                         .dep
                         .grd
                         .opt
                      :
           .done
              
    This will create the job chain:: 
    
        appl.jb
        appl_build.jb
        appl_init.jb                 (v)
        appl_init_emis.jb
        appl_init_obs.jb             (v)
        appl_init_obs_point.jb
        appl_init_obs_sat.jb
        appl_init_obs_valid.jb
        appl_run.jb                  (v)
        appl_run_iter-0001.jb        (v)
        appl_run_iter-0001_fwd.jb
        appl_run_iter-0001_dep.jb
        appl_run_iter-0001_grd.jb
        appl_run_iter-0001_opt.jb
        appl_run_iter-0002.jb        (v)
        appl_run_iter-0002_fwd.jb
        appl_run_iter-0002_dep.jb
        appl_run_iter-0002_grd.jb
        appl_run_iter-0002_opt.jb
                         :
        done.jb
        
    The jobs marked with "(v)" are virtual and actually not created;
    see the description of :ref:`virtual jobs <virtual-jobs>`.

    The iteration list is defined by the following properties:
        
    * a generic name used to read settings;
    * a format to create an element name (iteration) from an integer step number;
    * the initial step number;
    * a maximum posible step number, used to perform a loop over possible names.
        
    Example rcfile settings for name 'appl.run' that define the iteration steps::

        ! job iteration class:
        appl.run.class                                   :  utopya.UtopyaJobIteration
        
        ! generic for step name used in settings:
        appl.run.generic                                 :  iter-NNNN
        ! formatting rule for actual step names given
        ! an integer step numper; 
        ! syntax should follow str.format() rules ; 
        ! here 4 digits with zero padding:
        appl.run.step_format                             :  iter-{step:0>4}
        ! initial step numbers:
        appl.run.step_start                              :  1
        ! maximum possible number for defined format:
        appl.run.step_max                                :  9999
        ! optional step size, default 1:
        appl.run.step_size                               :  1

    For the above example, each iteration consists of 4 sub-jobs.
    Define these using the generic name; eventual mark the iteration job as virtual::
        
        ! sub list:
        appl.run.iter-NNNN.class                        :  utopya.UtopyaJobTree
        appl.run.iter-NNNN.virtual                      :  True
        appl.run.iter-NNNN.elements                     :  fwd dep grd opt

    The method :py:meth:`CheckContinuation` from the base class is re-implemented
    to decide if a next iteration step should be performed or that 
    the loop should be terminated, given a step number.
    For derived classes that implement a new iteration loop,
    it might be sufficent to only re-define this method too in order
    to terminate the loop at the right step.
 
   """
    
    def GetStepNumbers( self, element ) :
    
        """
        Extract step number from element name, and also return maximum number.
        This method is used by 'GetNextElement' and 'CheckContinuation'
        to translate a job name to a number and decide on the next step.
        
        Return value is a three element tupple with integers: 
           step,step_max,step_size
        
        Current implementation uses the step range setting in the rcfile
        (from 'step_start' to 'step_max') to perform a loop over possible
        step numbers. For each step, the 'step_format' is evaluated and
        compared to the provided element; if a match is found, the step is known.
        This brute-force test should in future be replaced by a more
        elegant method reading the number given the format.
        """
        
        # get generic version:
        name_generic = self._GetGenericName( self.name )

        # read iteration settings:
        step_format = self.GetSetting( name_generic+'.step_format' )
        step_start  = self.GetSetting( name_generic+'.step_start', 'int' )
        step_max    = self.GetSetting( name_generic+'.step_max', 'int' )
        step_size   = self.GetSetting( name_generic+'.step_size', 'int', default=1 )
        
        # get current step number ;
        # until we found how to read the number from the name 
        # and the format, just search ...
        found = False
        for step in range(step_start,step_max+1,step_size) :
            step_name = step_format.format(step=step)
            #print indent+'  y try step %i "%s" == "%s"' % (step,genericX,nameX)
            found = step_name == element
            if found : break
        #endfor
        if not found :
            self.logger.error( 'could not find step number in range %i .. %i in element "%s" for format "%s"' % (step_start,step_max,element,step_format) )
            raise Exception
        #endif
        
        # ok
        return step,step_max,step_size
        
    #enddef GetStepNumbers

    # *
    
    def GetVariables( self, element ) :
    
        """
        Return dictionairy with job variables,
        for this class the iteration step::
            
            { '__step__' : 4 }
        """
        
        # get current step number from element:
        step,step_max,step_size = self.GetStepNumbers( element )
        # fill:
        var = { '__step__' : step }
        
        # ok
        return var
        
    #enddef GetVariables
    
    # *
    
#    #def _GetGenericElement( self, element ) :
#    def _GetGenericElement( self, root_generic, element ) :
#        
#        """
#        Return generic name of iteration elements.
#        """
#        
#        # read iteration settings:
#        #element_generic = self.GetSetting( self.name+'.generic' )
#        element_generic = self.GetSetting( root_generic+'.generic' )
#                
#        # ok
#        return element_generic
#        
#    #endddef _GetGenericElement

    # *
    
    def GetNextElement( self, element, indent='' ) :
    
        """
        Return info on next iteration step after provided job name.
        Two values are returned: the next name, 
        and a bool flag to enable a continuation test.
        In this class the later is always True, to have the iteration loop
        terminated at some point defined by the 'CheckContinuation' method.
        
        As example, if the requested element is equal to 'iter-0004',
        then the returned tupple is ('iter-0005',True)
            
        If the requested element is 'None', then the name that corresponds
        to the first iteration step is returned based on the rcfile settings.
        """
        
        ## testing ...
        #print indent+'yyy UtopyaJobIteration/GetNextElement for "%s" after "%s"' % (self.name,element)
        
        # get generic name:
        name_generic = self._GetGenericName( self.name )
        
        # no name yet? then first iteration step:
        if element == 'None' :

            # first element:
            nextstep = self.GetSetting( name_generic+'.step_start', 'int' )
            ## testing ...
            #print indent+'  y first step number %i' % nextstep

        else :

            # get current step number from element:
            step,step_max,step_size = self.GetStepNumbers( element )
            ## testing ...
            #print indent+'  y extracted step number %i' % step
            # next:
            nextstep = step + step_size

        #endif
        
        # read iteration settings:
        step_format = self.GetSetting( name_generic+'.step_format' )
        # form element name from new step number:
        nextelement = step_format.format(step=nextstep )
        ## testing ...
        #print indent+'  y defined next element "%s"' % (nextelement)

        # ok
        return nextelement

    #enddef GetNextElement
    
    # *
    
    def _WithContinuationTest( self ) :
        
        """
        Return bool flag to enable a continuation test.
        For iteration class the value is True.
        """
        
        # enable test:
        return True
        
    #enddef _WithContinuationTest
    
    # *
    
    def CheckContinuation( self, element ) :
    
        """
        Check iteration status given job name.
        The job name could be used to derive the iteration step number
        to decide whether a maximum is reached,
        or it could be used to read output files to decide whether
        convergence is reached.
        
        Returns two str values:
         * order  : one of 'continue', 'finish', 'error'
         * msg    : informative message to explain the order
        
        In this implementation, the loop is finished if the maximum
        step number is reached, defined by the 'step_max' value in the settings.
        The current step number and the maximum are obtained from 
        a call to the :py:meth:`GetStepNumbers` method.
        """
        
        # get current step number and allowed maximum for element:
        step,step_max,step_size = self.GetStepNumbers( element )

        # dummy ...
        if step >= step_max :
            order = 'finish'
            msg = 'reached maximum number of iterations, finish'
        else :
            order = 'continue'
            msg = 'no iteration thresholds reached yet, continue'
        #endif
        
        # ok
        return order,msg
        
    #enddef CheckContinuation 

#endclass UtopyaJobIteration



# ======================================================================
# ===
# === UtopyaJobIteration with continuation check by file
# ===
# ======================================================================


class UtopyaJobIteration_CheckFile( UtopyaJobIteration ) :

    """
    UTOPyA JobIteration class with the :py:meth:`CheckContinuation` method
    defined to read instructions from an input file.
    """

    def CheckContinuation( self, element ) :
    
        """
        Check iteration status by reading a text file named::
        
           <name>.<element>.msg
           
        The text file should consist of 2 lines that are read and
        provided as return values:
        
         * order  : one of 'continue', 'finish', 'error'
         * msg    : informative message to explain the order
        """
        
        # modules:
        import os
        
        # message file:
        fname = '%s.%s.msg' % (self.name,element)
        # check ..
        if not os.path.isfile(fname) :
            self.logger.error( 'continuation message file not found: %s' % fname )
            self.logger.error( '  current directory : %s' % os.getcwd() )
            raise Exception
        #endif
        
        # open:
        f = open( fname, 'r' )
        # read:
        order = f.readline().strip()
        msg = f.readline().strip()
        # close:
        f.close()
        
        # ok
        return order,msg
        
    #enddef CheckContinuation

#endclass UtopyaJobIteration_CheckFile



# ======================================================================
# ===
# === UtopyaJobTask
# ===
# ======================================================================


class UtopyaJobTask( object ) :
    
    """
    Dummy class for illustration of the :py:class:`UtopyaJobStep` class
    and its derivatives.
    
    In this implementation, an optional str message could be passed
    on initialization, which is printed if pressent.
    """
    
    def __init__( self, msg=None ) :
        
        """
        Initialize application.
        If message 'msg' is supplied it is printed.
        """
        
        # print info ...
        if msg is not None : print( msg )
        
    #enddef __init__
    
#endclass UtopyaJobTask

# *

class UtopyaJobTaskSubmit( UtopyaJobTask ) :
    
    """
    Job task to create and submit a job(tree).
    
    The arguments specify the name of the jobstep (jobtree) settings
    and the name of the file with these settings.
    For example, with ``name`` equal to 'appl' the first lines
    of the setting could be::
    
        ! single job:
        appl.class        :  utopya.UtopyaJobStep

        ! task:
        appl.task.class   :  utopya.UtopyaJobTask
        appl.task.args    :  msg='Perform application task.'

    An instance of the from :py:class:`UtopyaJobStep` derived class will be created,
    and after initialization, the ``Run`` method of the instance is called.
    
    The optional message is passed to the initialization of the parent,
    which will print it if defined.
    """
    
    def __init__( self, name, rcfile, msg=None ) :
        
        """
        Create and submit job(tree).
        """
        
        # modules:
        import os
        
        # tools:
        import utopya
        
        # parent init:
        UtopyaJobTask.__init__( self, msg=msg )
        
        # info ...
        print( 'submit job task:' )
        print( '  name    : %s' % name )
        print( '  rcfile  : %s' % rcfile )
        
        # base ojbect with rcfile settings:
        urc = utopya_rc.UtopyaRc( rcfile=rcfile )
        
        # target for expanded rcfile:
        rcfile_out = '%s.rc' % name
        # info ...
        print( 'write rcfile: %s' % rcfile_out )
        # testing ...
        urc.rcf.WriteFile( rcfile_out )
        
        # read optional workdir:
        workdir = urc.GetSetting( '%s.workdir' % name, default='None' )
        # defined?
        if workdir != 'None' : 
            # expand:
            workdir = workdir.replace( '__NAME2PATH__', name.replace('.','/') )
            # create?
            if not os.path.isdir(workdir) : os.makedirs( workdir )
            # info ...
            print( 'change to workdir: %s' % workdir )
            # change:
            os.chdir( workdir )
        #endif

        # info ...
        print( 'create and submit ...' )        
        # job step class:
        jbsclass = urc.ImportClass( name+'.class' )
        # create job step object:
        jbs = jbsclass( name, rcfile )
        # write job file and start:
        jbs.Run()

    #enddef __init__
    
#endclass UtopyaJobTaskSubmit

# *

class UtopyaJobTaskRun( UtopyaJobTask ) :
    
    """
    Job task to run executable.
    
    The argument specifies the command::
    
        'appl.x --flag=1 input.txt'
    """
    
    def __init__( self, command='' ) :
        
        """
        Start executable.
        """
        
        # modules:
        import subprocess
        
        # parent init:
        UtopyaJobTask.__init__( self )
        
        # defined?
        if len(command) > 0 :
            # convert to list:
            cmnd = command.split()
            # run:
            subprocess.check_call( cmnd )
        #endif

    #enddef __init__
    
#endclass UtopyaJobTaskRun

# *

class UtopyaRcCommand( utopya_rc.UtopyaRc ) :
    
    """
    Job task to run executable.
    
    Command to run is read from settings::
    
        [<rcbase>.]command    :  mpirun -n 4 prog.x
    
    Optionally define to run the command in a shell,
    this could be used to change resourcs that can only be set using shell options::
    
        [<rcbase>.]shell      :  True
        [<rcbase>.]command    :  ulimit -s unlimited ; mpirun -n 4 prog.x
    
    """
    
    def __init__( self, rcfile, rcbase='', env={}) :

        """
        Run command.
        """
        
        # tools:
        import utopya_tools

        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # command string:
        command = self.GetSetting( 'command' )
        in_shell = self.GetSetting( 'shell', totype='bool', default=False )
        
        # switch:
        if in_shell :
            # info ..
            self.logger.info( 'run shell command:' )
            for cmd in command.split(';') :
                self.logger.info( '  %s' % cmd )
            #endfor
            # run:
            utopya_tools.Call( command, shell=True )
        else :
            # info ..
            self.logger.info( 'run command:' )
            self.logger.info( '  %s' % command )
            # run:
            utopya_tools.Call( command.split() )
        #endif
        
    #enddef __init__
    
#endclass UtopyaRcCommand
        
        

# ======================================================================
# ===
# === End
# ===
# ======================================================================


