#! /usr/bin/env python3

"""
***************************
``utopya_jobscript`` module
***************************

Introduction
============

A job script is a (usually) small text file to do something important.
A simple example is::

    #! /usr/bin/env python
    
    # do something:
    print( "boe!" )
    
Such a script could be run in the foreground, in which case the user has
to wait for the job to be finished before control is regained.
Standard output ("boe!") and eventually standard error are printed directly
and can be watched by the user.

When the script is run in background, control is given back to the user
while the script remains running. The standard output and error should
be redirected to files.

Alternatively, the script could be submitted to a batch system.
In this case, options to identify the job, to specify destination of
standard output and error, and to request resoures (memory, cpu's),
could be inserted in the header of the script::

    #! /usr/bin/env python

    #BSUB -J myjob
    #BSUB -oo myjob.out
    #BSUB -eo myjob.err

    # do something:
    print( "boe!" )

The classes provided by this module facilitate creation of job scripts
for which the run destination (foreground, background,  batch system)
is flexible.

Jobscripts to run in foreground
===============================

To run a job in the foreground, use one of the following classes:

* :py:class:`UtopyaJobScriptForeground`
* :py:class:`UtopyaJobScriptRedirect`, which redirects standard output (and error) to files.

Jobscripts to run in background
===============================

To run in background, use the :py:class:`UtopyaJobScriptBackground` class.

Jobscripts to be submitted to a batch system
============================================

High performace clusters with a high number of processors
and many users logged in at the same time are always equiped
with a batch system to handle jobs.
Batch jobs are submitted to a queue, and the batch system
empties the queue by assigning jobs to first available processors.

Special commands are required to submit jobs to the queue,
list the currently submitted and running jobs, and eventually
remove jobs from the queue.

Batch job files typically have special comments in the top
to tell the batch system the destination queue, the name of
output/error files to be used, required memory and maximum
run time, etc.

Which batch system is available usually depends on the machine vendor
and/or administrator.
Each type of batch system has its own job handling command
and format for the batch options in the top of the job file.
For each type, a seperate class needs to be defined to handle
creation and submission.
A base class :py:class:`UtopyaJobScriptBatch` is provided from
which batch type specific class could be derived;
see it's documentation for the methods to be re-defined.

The following specific batch systems are already supported:

* For LSF, which uses the 'bsub' command to submit, use the :py:class:`UtopyaJobScriptBatchLSF` class.
* For SLURM, which uses the 'sbatch' command to submit, use the :py:class:`UtopyaJobScriptBatchSlurm` class.
* For PBS, which uses the 'qsub' command to submit, use the :py:class:`UtopyaJobScriptBatchPBS` class.
* For the IBM LoadLeveler queue, use the :py:class:`UtopyaJobScriptBatchLoadLeveler` class.

If queue systems are fully occupied, testing creation of job files could suffer
from long waiting times. To avoid wasting of precious development time, the
special :py:class:`UtopyaJobScriptBatchTest` class is provided.
This will create jobs using fake job options, and run the script in foreground
while redirecting standard output and error.


Class hierarchy
===============

The classes provided by this module have been derived with the following hierchy:

* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaJobScript`

    * :py:class:`.UtopyaJobScriptForeground`
    * :py:class:`.UtopyaJobScriptRedirect`
    * :py:class:`.UtopyaJobScriptBackground`
    * :py:class:`.UtopyaJobScriptBatch`

      * :py:class:`.UtopyaJobScriptBatchTest`
      * :py:class:`.UtopyaJobScriptBatchLSF`
      * :py:class:`.UtopyaJobScriptBatchSlurm`
      * :py:class:`.UtopyaJobScriptBatchPBS`
      * :py:class:`.UtopyaJobScriptBatchLoadLeveler`

Classes
=======

"""

# ======================================================================
# ===
# === modules
# ===
# ======================================================================

# base class:
from utopya_base import UtopyaBase


# ======================================================================
# ===
# === JobScript
# ===
# ======================================================================

class UtopyaJobScript( UtopyaBase ) :

    """
    Base class for an object that can be used to create and submit
    a job script.
    
    Derived classes probably need to re-define the :py:meth:`Submit` methode only.

    Example of usage::
        
        # init object:
        jb = UtopyaJobScript()
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # submit created script 
        # (here dummy method that issues an error):
        jb.Submit( 'myjob.jb' )
        
    """
    
    def __init__( self ) :
    
        """
        Initialize job script.        
        Defines an attribute 'logger' to issue messages.
        """
        
        # init from parent:
        UtopyaBase.__init__( self )
        
    #enddef __init__
    
    # *
    
    def Create( self, jbfile, lines ) :
    
        """
        Create job file with a name provided by `jbfile`,
        and content in the `lines` argument.
        
        The 'lines' should either be a list of 'str' objects,
        or a single 'str' object, but in both cases newline
        characters should be included.
        
        The created job file is made executable to allow
        execution from a command line.
        """
        
        # modules:
        import os
        import stat
        
        # directory:
        jbdir = os.path.dirname( jbfile )
        # create if not present yet:
        if len(jbdir) > 0 :
            if not os.path.isdir(jbdir) : os.makedirs( jbdir )
        #endif
    
        # create new file:
        fid = open( jbfile, 'w' )
        # write:
        if type(lines) == list :
            fid.writelines( lines )
        else :
            fid.write( lines )
        #endif
        # close:
        fid.close()                                   

        # make executable:
        os.chmod( jbfile, stat.S_IRWXU + stat.S_IRGRP + stat.S_IXGRP  + stat.S_IROTH + stat.S_IXOTH )
        
    #enddef Create
    
    # *
    
    def GetOptionsRc( self, rcfile, name='batch', rcbase='', env={} ) :
        
        """
        Template for derived classes that need to include batch job
        options in a job file. This version returns an empty line.
        """
        
        # ok
        return ''
        
    #enddef GetOptionsRc
    
    # *
    
    def Submit( self, jbfile ) :
    
        """
        Template for job submision methode by derived classes.
        """
        
        # not yet ...
        self.logger.error( 'not implemented for base class' )
        
    #enddef Submit
    
#endclass UtopyaJobScript


# ======================================================================
# ===
# === jobs running in foreground
# ===
# ======================================================================


class UtopyaJobScriptForeground( UtopyaJobScript ) :

    """
    Class to create a job that runs in foreground.

    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptForeground()
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # run script in foreground:
        jb.Submit( 'myjob.jb' )

    """
    
    def Submit( self, jbfile ) :
    
        """
        Run job file in foreground.
        Standard output and error are not redirected.
        """
        
        # modules:
        import os
        import subprocess
        
        # info ...
        self.logger.info( 'run "%s" in foreground ...' % jbfile )
        
        # command:
        if os.path.isabs( jbfile ) :
            command = [ jbfile ]
        else :
            command = [ os.path.join( os.curdir, jbfile ) ]
        #endif
        
        # run:
        subprocess.check_call( command )
        
    #enddef Submit

#endclass UtopyaJobScriptForeground


# ======================================================================
# ===
# === jobs running in foreground with redirection of output and error
# ===
# ======================================================================


class UtopyaJobScriptRedirect( UtopyaJobScript ) :

    """
    Class to create job for running in foreground,
    while re-directing std.output and std.error to files.

    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptRedirect()
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # run script in foreground, redirect output:
        jb.Submit( 'myjob.jb' )

    """
    
    def Submit( self, jbfile ) :
    
        """
        Run job file in foreground, but with standard output and error
        redirected to files.
        The output and error files have the same name as the job file but extensions
        '.out' and '.err' respectively.
        """
        
        # modules:
        import os
        import subprocess
        
        # info ...
        self.logger.info( 'run "%s" in foreground (redirect output and error) ...' % jbfile )
        
        # split in base name and extension:
        bname,ext = os.path.splitext( jbfile )
        
        # open output and error:
        f_stdout = open( bname+'.out', 'w' )
        f_stderr = open( bname+'.err', 'w' )
        
        # command:
        if os.path.isabs( jbfile ) :
            command = [ jbfile ]
        else :
            command = [ os.path.join(os.curdir,jbfile) ]
        #endif
        
        # run:
        subprocess.check_call( command, stdout=f_stdout, stderr=f_stderr )
        
    #enddef Submit

#endclass UtopyaJobScriptRedirect


# ======================================================================
# ===
# === jobs running in background
# ===
# ======================================================================


class UtopyaJobScriptBackground( UtopyaJobScript ) :

    """
    Class to create job for running in background.

    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptBackground()
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # run script in background:
        jb.Submit( 'myjob.jb' )

    """
    
    def Submit( self, jbfile ) :
    
        """
        Run job file in background.
        Standard output and error are written to files
        with the same name as the job file but extensions
        '.out' and '.err' respectively.
        """
        
        # modules:
        import os
        import subprocess
        
        # info ...
        self.logger.info( 'run "%s" in background ...' % jbfile )
        
        # split in base name and extension:
        bname,ext = os.path.splitext( jbfile )
        
        # open output and error:
        f_stdout = open( bname+'.out', 'w' )
        f_stderr = open( bname+'.err', 'w' )
        
        # command:
        if os.path.isabs( jbfile ) :
            command = [ jbfile ]
        else :
            command = [ os.path.join(os.curdir,jbfile) ]
        #endif
        
        # run:
        subprocess.Popen( command, stdout=f_stdout, stderr=f_stderr )
        
    #enddef Submit

#endclass UtopyaJobScriptBackground


# ======================================================================
# ===
# === UtopyaJobScriptBatch
# ===
# ======================================================================

class UtopyaJobScriptBatch( UtopyaJobScript ) :

    """
    Base class to create and submit a batch job script.    
    The base class itself does not support a particular system,
    but only provides generic methods that are shared by the 
    specific classes.
    Derived classes probably need to re-define the following methods only:
    
    * :py:meth:`GetOptionsRc` that should call the same method
      from the parent (this) class with proper arguments;
        
    * :py:meth:`Submit` to submit a job script to batch queue.
    
    Example of usage of the base class, derived classes should be
    used in the same way::
        
        # init object:
        jb = utopya.UtopyaJobScriptBatch()
        
        # obtain line with job options from rcfile:
        options = jb.GetOptionsRc( 'UtopyaJobScriptBatch.rc', rcbase='appl', \\
                                        env={'name':'myjob'} )
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( options )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # submit:
        jb.Submit( 'myjob.jb' )
        
    """
    
    def GetOptions( self, keys, options, 
                       comment='#', prefix='BATCH', 
                       arg='-', assign=' ' ) :
    
        """
        Return str line (with newline breaks) with job options.
        For example for a load-leveler job this could look like::

            '#@ job_name = myjob\\n#@ output = myjob.out\\n#@ error = myjob.err\\n#@ queue'

        In an actual job script this is then expanded into seperate lines::
            
            #@ job_name = myjob
            #@ output   = myjob.out                     
            #@ error    = myjob.err                    
            #@ queue

        With a different queue system, the job options should probably look different.
        For an LSF job for example, the options could look like::
        
            #BSUB -J myjob
            #BSUB -o myjob.out
            #BSUB -e myjob.err
         
        Although formatting is different, the options in general consist
        of a lines starting with unique comment pragma, followed by a flag,
        and usually a value.
        Flags might be preceeded by a '-' sign; flag and value might be seperated
        by whitespace or a '=' sign.
        In a general notation the format of a job option is::
            
            <comment><prefix> <arg><flag>[<assign><value>]
            
        The formatting to be used is defined by the optional keyword arguments 
        passed to this method.
        For the LoadLeveler job for example the formatting is defined by::
            
          comment = '#'
          prefix  = '@'
          arg     = ''
          assign  = '='
          
        which will give option lines::
          
          #@ <flag>[=<value>]

        For the LSF system these keywords are necessary::

          comment = '#'
          prefix  = 'BSUB'
          arg     = '-'
          assign  = ' '
          
        which will give the lines::
          
          #BSUB -<flag>[ <value>]
        
        The flag/value pairs should be passed to this method by two arguments:
        a list of 'keys' that defines the order of the flag/value pairs,
        and a dictionairy 'options' with flag/value tupples for each of the 
        elements in the 'keys' list.
        For example for the LoadLeveler job, the arguments could be:::

          keys = ['name','out','err','end']
          options = { 'name'  : ('job_name','myjob'),
                      'out'   : ('output','myjob.out'),
                      'err'   : ('error','myjob.err'),
                      'end'   : ('queue',None) }

        and for the LSF job::

          keys = ['name','out','err']
          options = { 'name' : ('J','myjob'),
                      'out'  : ('o','myjob.out'),
                      'err'  : ('e','myjob.err') }

        """
        
        # write lines:
        lines = ''
        for key in keys :
            flag,value = options[key]
            if value is None :
                line = '%s%s %s%s' % (comment,prefix,arg,flag)
            else :
                line = '%s%s %s%s%s%s' % (comment,prefix,arg,flag,assign,value)
            #endif
            lines = lines+line+'\n'
        #endfor
        
        # ok
        return lines
        
    #enddef GetOptions
    
    # *
    
    def _RcGet( self, rcf, name, rcbase='' ) :
    
        """
        Read setting for 'name' from rcfile.
        Optional 'rcbase' is used as prefix for the name.
        Default values might be defined as if 'rcbase' is set to '*'::
        
          *.<name>               :  <default>
          [<rcbase>.]<name>      :  <value>
        """
        
        # actual key:
        akey = ('%s.%s' % (rcbase,name)).lstrip('.')
        # present ?
        if rcf.has_key(akey) :
            # read:
            value = rcf.get( akey )
        else :
            # default key:
            dkey = '*.%s' % name
            # read default:
            if rcf.has_key(dkey) :
                value = rcf.get( dkey )
            else :
                self.logger.error( 'neither default or actual keys defined in "%s":' % (rcf.filename) )
                self.logger.error( '  %s' % dkey )
                self.logger.error( '  %s' % akey )
                raise Exception
            #endif
        #endif
        
        # ok
        return value
        
    #enddef _RcGet
        
    # *
    
    def GetOptionsRc( self, rcfile, name='batch', rcbase='', env={} ) :
    
        """
        Read batch options and formatting rules from rcfile settings,
        and pass these to the :py:meth:`GetOptions` method.
        The result is similar, a 'str' line with batch options
        to be written to a job file, including newline characters.
        
        The 'name' is the first part of the rcfile keys that
        should be used.

        The 'rcbase' is an optional prefix for the name.
        This allows a rcfile to have multiple batch job
        definitions, each with a different 'rcbase' to be used for different
        job files. 

        Example rcfile settings for rcbase 'appl' and name 'batch'::
        
            ! which keywords:
            appl.batch.options           :  jobname output error nodes memory queue
        
            ! values:
            appl.batch.option.jobname    :  n %(env:name)
            appl.batch.option.output     :  o %(jobname).out
            appl.batch.option.error      :  e %(jobname).err
            appl.batch.option.nodes      :  r nodes:4
            appl.batch.option.memory     :  r memory:2G
            appl.batch.option.queue      :  q
        
            ! batch job option format keyword:
            appl.batch.format            :  myformat

            ! Define format of batch options, e.g.:
            !   #BATCH -flag value
            ! If whitespace is essential, enclose value by single quotes:
            myformat.comment       :  #
            myformat.prefix        :  BATCH
            myformat.arg           :  -
            myformat.assign        :  ' '
            ! format of template to subsitute (environment) values:
            myformat.template      :  %(key)
            myformat.envtemplate   :  %(env:key)
        
        The first setting is a list of option keys, here 'jobname', 'output', etc.
        
        For each option key, a flag/value pair needs to be defined;
        in the example above, for key 'jobname' the flag/value pair 
        'J %(env:name)' is defined.
        Typically, the flag could consist of just a cryptic
        letter ('J') as required by the batch system,
        while the option key is longer and more descriptive.
        
        The first part of a flag/value is the flag, seperated by whitespace
        from the value(s) which form the remainder of the line.
        Multiple options could therefore share the same flag;
        this is useful when the job scheduler requires the same flag
        to be used for different settings, for example to set resources.
        
        The value could be defined literaly, such as "nodes:4" is defined
        above as the value for the number of nodes.
        
        The value could also contain templates for substitution of a value
        from a flag/value pair assigned to another key.
        In the example the value assigned to "jobname" is for example also
        used for the output and error files.
        The format of the template is defined here as '%(key)', where
        only the presence of the word "key" is required.
        A loop over all keys ("jobname", "output", etc) will be performed
        and the word "key" in the template is replaced by the current key;
        if the result is found in a value, it is decided that the template
        is used and it is replaced by the corresponding value.
        
        A special substitution is defined for variables from the so-called environment.
        In this case, the template should contain the word "env:key".
        The environment is an optional dictionairy that is passed by the calling
        method and contains specific values at run time; a loop over the environment
        keys is performed to search for matching templates. For example, the following
        environment contains the job name to be used::
            
            env = { 'name' : 'myjob' }
                
        With this environment, the above example will result in the
        following job options (actually a single str with newlines)::
        
            #BATCH -n myjob
            #BATCH -o myjob.out
            #BATCH -e myjob.err
            #BATCH -r nodes:4
            #BATCH -r memory:2G
            #BATCH -q

        """
        
        # tools:
        import rc
        
        # init result:
        lines = ''
        
        # read settings:
        rcf = rc.RcFile( rcfile )
        
        # format name:
        format = self._RcGet( rcf, name+'.format', rcbase=rcbase )
        
        # read batch job settings, remove quotes:
        comment     = rcf.get( format+'.comment'     ).replace("'","")
        prefix      = rcf.get( format+'.prefix'      ).replace("'","")
        arg         = rcf.get( format+'.arg'         ).replace("'","")
        assign      = rcf.get( format+'.assign'      ).replace("'","")
        template    = rcf.get( format+'.template'    ).replace("'","")
        envtemplate = rcf.get( format+'.envtemplate' ).replace("'","")
        
        # option keywords:
        keys = self._RcGet( rcf, name+'.options', rcbase=rcbase ).split()
        # collect flag/value pairs:
        opts = {}
        for key in keys :
            # corresponding flag/value pair:
            flagvalue = self._RcGet( rcf, name+'.option.'+key, rcbase=rcbase )
            # split at first whitespace if present:
            if ' ' in flagvalue :
                # split:
                flag,value = flagvalue.split(' ',1)
                # cleanup:
                flag  = flag.strip()
                value = value.strip()
                # replace environment values if necessary;
                # loop over environment keys:
                for kk in env.keys() :
                    # corresponding value:
                    vv = env[kk]
                    # search pattern:
                    tmpl = envtemplate.replace('key',kk)
                    # present ?
                    if tmpl in value :
                        # replace:
                        value = value.replace(tmpl,vv)
                    #endif
                #endfor
            else :
                # no value, set flag only:
                flag,value = flagvalue,None
            #endif
            # store:
            opts[key] = (flag,value)
        #endfor

        # replace option templates:
        for key in keys :
            # corresponding flag/value pair:
            flag,value = opts[key]
            # skip empty:
            if value is None : continue
            # loop over option keys:
            for kk in keys :
                # prevent infinite loop ...
                if kk == key : continue
                # corresponding flag/value pair:
                ff,vv = opts[kk]
                # skip empty:
                if vv is None : continue
                # search pattern:
                tmpl = template.replace('key',kk)
                # present ?
                if tmpl in value :
                    # replace:
                    value = value.replace(tmpl,vv)
                    # store:
                    opts[key] = (flag,value)
                #endif
            #endfor
        #endfor
        
        # format lines:
        lines = self.GetOptions( keys, opts, comment=comment, prefix=prefix, 
                                   arg=arg, assign=assign )
        
        # ok
        return lines
        
    #enddef GetOptionsRc
    
    # *
    
    def Submit( self, jbfile ) :
        
        """
        Template for submit method.
        Derived classes should re-implement this method to
        submit the named 'jbfile' to a batch queue.
        """
        
        # info ..
        self.logger.error( 'This is a template method, should not be called.' )
        
    #enddef Submit
    
#endclass UtopyaJobScriptBatch


# ======================================================================
# ===
# === test batch jobs
# ===
# ======================================================================


class UtopyaJobScriptBatchTest( UtopyaJobScriptBatch ) :

    """
    Class to create test job files that are not submitted to a queue
    but simply run in foreground. Useful for testing job script creation.

    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptBatchTest()

        # obtain line with job options from rcfile:
        options = jb.GetOptionsRc( 'settings-UtopyaJobScriptBatchTest.rc', rcbase='appl', \\
                                        env={'name':'myjob'} )
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( options )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # submit:
        jb.Submit( 'myjob.jb' )
        
    """
    
    def GetOptionsRc( self, rcfile, rcbase='', env={} ) :
    
        """
        Return str line (with newline characters) with job options based on rcfile settings.
        
        The rcfile settings should start with '[<rcbase>.]batch.test',
        where the rcbase might be empty or equal to '*' for default settings.
        
        Example settings for rcbase 'appl'::
                    
            ! job format for this application:
            appl.batch.test.format            :  test_format
            ! which keywords:
            appl.batch.test.options           :  name output error
            ! flags and values:
            appl.batch.test.option.name       :  J myjob
            appl.batch.test.option.output     :  oo %(name).out
            appl.batch.test.option.error      :  eo %(name).err
    
            ! Define format of batch options, e.g.:
            !   #TEST -flag value
            test_format.comment       :  #
            test_format.prefix        :  TEST
            test_format.arg           :  '-'
            test_format.assign        :  ' '
            test_format.template      :  %(key)
            test_format.envtemplate   :  %(env:key)

        This will return the following job options as a str with newline characters::
    
            #TEST -J myjob
            #TEST -oo myjob.out
            #TEST -eo myjob.err
            
        """
        
        # call parent:
        lines = UtopyaJobScriptBatch.GetOptionsRc( self, rcfile, 'batch.test', rcbase=rcbase, env=env )
        
        # ok
        return lines
        
    #enddef GetOptionsRc

    # *
    
    def Submit( self, jbfile ) :
    
        """
        Test routine that runs the batch job file in background.

        Standard output and error are written to files
        with the same name as the job file but extensions
        '.out' and '.err' respectively.
        """
        
        # modules:
        import os
        import subprocess
        
        # info ...
        self.logger.info( 'batch job test: run "%s" in background ...' % jbfile )
        
        # split in base name and extension:
        bname,ext = os.path.splitext( jbfile )
        
        # open output and error:
        f_stdout = open( bname+'.out', 'w' )
        f_stderr = open( bname+'.err', 'w' )
        
        # command:
        if os.path.isabs( jbfile ) :
            command = [ jbfile ]
        else :
            command = [ os.path.join(os.curdir,jbfile) ]
        #endif
        
        # run:
        subprocess.Popen( command, stdout=f_stdout, stderr=f_stderr )
        
    #enddef Submit

#endclass UtopyaJobScriptBatchTest


# ======================================================================
# ===
# === jobs submitted to LSF
# ===
# ======================================================================


class UtopyaJobScriptBatchLSF( UtopyaJobScriptBatch ) :

    """
    Class to create job for submission to LSF batch system.
    
    Example of job options::
    
      #BSUB -J myjob
      #BSUB -oo myjob.out
      #BSUB -eo myjob.err
    
    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptBatchLSF()

        # obtain line with job options from rcfile:
        options = jb.GetOptionsRc( 'UtopyaJobScriptBatchLSF.rc', rcbase='appl', \\
                                env={'name':'myjob'} )
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( options )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # submit:
        jb.Submit( 'myjob.jb' )
    
    See also man pages of batch job commands:
            
    * :manpage:`bsub`
    * :manpage:`bjobs`
    * :manpage:`bkill`

    """
    
    def GetOptionsRc( self, rcfile, rcbase='', env={} ) :
    
        """
        Return str line (with newline characters) with job options based on rcfile settings.
        
        The rcfile settings should start with '[<rcbase>.]batch.lsf',
        where the rcbase might be empty or equal to '*' for default settings.
        
        Example settings for rcbase 'appl'::
                    
            ! job format for this application:
            appl.batch.lsf.format            :  lsf_format
            ! which keywords:
            appl.batch.lsf.options           :  name output error
            ! values:
            appl.batch.lsf.option.name       :  J myjob
            appl.batch.lsf.option.output     :  oo %(name).out
            appl.batch.lsf.option.error      :  eo %(name).err
    
            ! Define format of batch options, e.g.:
            !   #BSUB -flag value
            lsf_format.comment       :  #
            lsf_format.prefix        :  BSUB
            lsf_format.arg           :  '-'
            lsf_format.assign        :  ' '
            lsf_format.template      :  %(key)
            lsf_format.envtemplate   :  %(env:key)

        This will return the following job options as a str with newline characters::
    
            #BSUB -J myjob
            #BSUB -oo myjob.out
            #BSUB -eo myjob.err

        """
        
        # call parent:
        lines = UtopyaJobScriptBatch.GetOptionsRc( self, rcfile, 'batch.lsf', rcbase=rcbase, env=env )
        
        # ok
        return lines
        
    #enddef GetOptionsRc
    
    # *
    
    def Submit( self, jbfile ) :
    
        """
        Submit job file.
        Information on job id and commands to follow and cancel the job
        are written to a file with the same name but extension '.info' .
        """
        
        # modules:
        import sys
        import os
        import subprocess        

        # basename for scripts etc is name of rcfile minus extension:
        bname,ext = os.path.splitext(jbfile)    

        # setup command line, should read jbfile from standard input:
        #    bsub < jbfile
        command = 'bsub'
        # input lines:
        f = open( jbfile, 'r' )
        jblines = f.read()
        f.close()
        # info ...
        self.logger.info( 'command: %s' % command )
        self.logger.info( 'stdin  : %s' % jbfile )

        # init submission info file:
        infotext = []
        infotext.append(  '\n' )

        # start command, trap errors (command does not exist etc):
        try:
            # handle standard input, output, and error later on by 'communicate' :
            p = subprocess.Popen( command.split(), stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        except :
            self.logger.error( sys.exc_info()[1] )
            self.logger.error( 'from command : %s' % command )
            raise Exception
        #endtry
        # send job lines to standard input, obtain standard output and error in return:
        stdout,stderr = p.communicate( input=jblines.encode('utf-8') )
        # convert:
        stdout = stdout.decode('utf-8').strip()
        stderr = stderr.decode('utf-8').strip()
        # error ?
        if len(stderr) > 0 :
            self.logger.error( stderr )
            self.logger.error( 'from command : %s' % command )
            raise Exception
        #endif
        # add to help info message:
        infotext = infotext + stdout.split('\n')
        infotext = infotext + stderr.split('\n')
        # standard output is:
        #   <jobname> <jobnr>
        # extract job id from first line:
        try :
            job_id = stdout.split()[1].strip('<>')
        except :
            self.logger.error( 'could not extract job id from stdout:' )
            self.logger.error( stdout )
            raise
        #endtry
    
        # help text:
        infotext.append(  '\n' )
        infotext.append(  'To manage LSF jobs:\n' )
        infotext.append(  '\n' )
        infotext.append(  '  bsub < %s   # submit to queue\n' % jbfile )
        infotext.append(  '  bjobs                # list your current jobs\n' )
        infotext.append(  '  bjobs %s             # list this job\n' % job_id )
        infotext.append(  '  bkill %s             # kill this job\n' % job_id )
        infotext.append(  '\n' )

        # write to loging system:
        for line in infotext :
            self.logger.info( line.rstrip() )
        #endfor

        # log file for submission:
        job_info = bname+'.info'
        # write to file:
        f = open( job_info, 'w' )
        f.writelines(infotext)
        f.close()

    #enddef Submit

#endclass UtopyaJobScriptBatchLSF


# ======================================================================
# ===
# === jobs submitted to PBS
# ===
# ======================================================================


class UtopyaJobScriptBatchPBS( UtopyaJobScriptBatch ) :

    """
    Class to create job for submission to PBS batch system.
    
    Example of job options::
    
      #PBS -N myjob
      #PBS -o myjob.out
      #PBS -e myjob.err
    
    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptBatchPBS()

        # obtain line with job options from rcfile:
        options = jb.GetOptionsRc( 'UtopyaJobScriptBatchPBS.rc', rcbase='appl', \\
                                env={'name':'myjob'} )
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( options )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # submit:
        jb.Submit( 'myjob.jb' )
    
    See also man pages of batch job commands:
            
    * :manpage:`qsub`
    * :manpage:`qscan`
    * :manpage:`qdel`
    * :manpage:`qstat`

    """
    
    def GetOptionsRc( self, rcfile, rcbase='', env={} ) :
    
        """
        Return str line (with newline characters) with job options based on rcfile settings.
        
        The rcfile settings should start with '[<rcbase>.]batch.pbs',
        where the rcbase might be empty or equal to '*' for default settings.
        
        Example settings for rcbase 'appl'::
                    
            ! job format for this application:
            appl.batch.pbs.format            :  pbs_format
            ! which keywords:
            appl.batch.pbs.options           :  name output error
            ! values:
            appl.batch.pbs.option.name       :  J myjob
            appl.batch.pbs.option.output     :  o %(name).out
            appl.batch.pbs.option.error      :  e %(name).err
    
            ! Define format of batch options, e.g.:
            !   #PBS -flag value
            pbs_format.comment       :  #
            pbs_format.prefix        :  PBS
            pbs_format.arg           :  '-'
            pbs_format.assign        :  ' '
            pbs_format.template      :  %(key)
            pbs_format.envtemplate   :  %(env:key)

        This will return the following job options as a str with newline characters::
    
            #PBS -J myjob
            #PBS -o myjob.out
            #PBS -e myjob.err

        """
        
        # call parent:
        lines = UtopyaJobScriptBatch.GetOptionsRc( self, rcfile, 'batch.pbs', rcbase=rcbase, env=env )
        
        # ok
        return lines
        
    #enddef GetOptionsRc
    
    # *
    
    def Submit( self, jbfile ) :
    
        """
        Submit job file.
        Information on job id and commands to follow and cancel the job
        are written to a file with the same name but extension '.info' .
        """
        
        # modules:
        import sys
        import os
        import subprocess        

        # basename for scripts etc is name of rcfile minus extension:
        bname,ext = os.path.splitext(jbfile)    

        # setup command line, should read jbfile from standard input:
        #    qsub jbfile
        command = 'qsub'
        # info ...
        self.logger.info( 'command    : %s' % command )
        self.logger.info( 'jobs script: %s' % jbfile )

        # init submission info file:
        infotext = []
        infotext.append(  '\n' )

        # start command, trap errors (command does not exist etc):
        try:
            # handle standard input, output, and error later on by 'communicate' :
            p = subprocess.Popen( [command,jbfile], 
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        except :
            self.logger.error( sys.exc_info()[1] )
            self.logger.error( 'from command : %s' % command )
            raise Exception
        #endtry
        # send job lines to standard input, obtain standard output and error in return:
        stdout,stderr = p.communicate()
        # add to help info message:
        infotext = infotext + stdout.split('\n')
        infotext = infotext + stderr.split('\n')
        # standard output is:
        #   <jobnr>.<jbqueue>
        # extract job id from first line:
        job_id = stdout.strip()
    
        # help text:
        infotext.append(  '\n' )
        infotext.append(  'To manage PBS jobs:\n' )
        infotext.append(  '\n' )
        infotext.append(  '  qsub %s   # submit to queue\n' % jbfile )
        infotext.append(  '  qscan [-u $USER]     # list your current jobs\n' )
        infotext.append(  '  qscan %s             # list this job\n' % job_id )
        infotext.append(  '  qdel %s              # kill this job\n' % job_id )
        infotext.append(  '\n' )

        # write to loging system:
        for line in infotext :
            self.logger.info( line.rstrip() )
        #endfor

        # log file for submission:
        job_info = bname+'.info'
        # write to file:
        f = open( job_info, 'w' )
        f.writelines(infotext)
        f.close()

    #enddef Submit

#endclass UtopyaJobScriptBatchPBS
    
    
# ======================================================================
# ===
# === jobs submitted to Slurm
# ===
# ======================================================================


class UtopyaJobScriptBatchSlurm( UtopyaJobScriptBatch ) :

    """
    Class to create job for submission to SLURM batch system.
    
    Example of job options::
    
      #SBATCH --job-name=myjob
      #SBATCH -output=myjob.out
      #SBATCH -error=myjob.err
    
    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptBatchSlurm()

        # obtain line with job options from rcfile:
        options = jb.GetOptionsRc( 'UtopyaJobScriptBatchSlurm.rc', rcbase='appl', \\
                                env={'name':'myjob'} )
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( options )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # submit:
        jb.Submit( 'myjob.jb' )
    
    See also man pages of batch job commands:
            
    * :manpage:`sbatch`
    * :manpage:`squeue`
    * :manpage:`scancel`

    """
    
    def GetOptionsRc( self, rcfile, rcbase='', env={} ) :
    
        """
        Return str line (with newline characters) with job options based on rcfile settings.
        
        The rcfile settings should start with '[<rcbase>.]batch.slurm',
        where the rcbase might be empty or equal to '*' for default settings.
        
        Example settings for rcbase 'appl'::
                    
            ! job format for this application:
            appl.batch.slurm.format            :  slurm_format
            ! which keywords:
            appl.batch.slurm.options           :  name output error
            ! values:
            appl.batch.slurm.option.name       :  job-name myjob
            appl.batch.slurm.option.output     :  output %(name).out
            appl.batch.slurm.option.error      :  error %(name).err
    
            ! Define format of batch options, e.g.:
            !   #SBATCH --flag=value
            slurm_format.comment       :  #
            slurm_format.prefix        :  SBATCH
            slurm_format.arg           :  '--'
            slurm_format.assign        :  '='
            slurm_format.template      :  %(key)
            slurm_format.envtemplate   :  %(env:key)

        This will return the following job options as a str with newline characters::
    
            #SBATCH --job-name=myjob
            #SBATCH --output=myjob.out
            #SBATCH --error=myjob.err

        """
        
        # call parent:
        lines = UtopyaJobScriptBatch.GetOptionsRc( self, rcfile, 'batch.slurm', rcbase=rcbase, env=env )
        
        # ok
        return lines
        
    #enddef GetOptionsRc
    
    # *
    
    def Submit( self, jbfile ) :
    
        """
        Submit job file.
        Information on job id and commands to follow and cancel the job
        are written to a file with the same name but extension '.info' .
        """
        
        # modules:
        import sys
        import os
        import subprocess        

        # basename for scripts etc is name of rcfile minus extension:
        bname,ext = os.path.splitext(jbfile)    

        # setup command line, should read jbfile from standard input:
        #    sbatch jbfile
        command = 'sbatch %s' % jbfile
        # info ...
        self.logger.info( 'command: %s' % command )

        # init submission info file:
        infotext = []
        infotext.append(  '\n' )

        # start command, trap errors (command does not exist etc):
        try:
            # handle standard input, output, and error later on by 'communicate' :
            p = subprocess.Popen( command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        except :
            self.logger.error( sys.exc_info()[1] )
            self.logger.error( 'from command : %s' % command )
            raise Exception
        #endtry
        # send job lines to standard input, obtain standard output and error in return:
        stdout,stderr = p.communicate()
        # convert:
        stdout = str(stdout)
        stderr = str(stderr)
        # add to help info message:
        infotext = infotext + stdout.split('\n')
        infotext = infotext + stderr.split('\n')
        # standard output differs per version, these formats were found:
        #   Submitted batch job 123456
        #   MYJOBNAME <123456>
        # extract job id from line:
        try :
            # switch:
            if stdout.startswith('Submitted batch job') :
                # extract from 4th key:
                job_id = stdout.split()[3]
            else :
                # split at whitepace, extract 2nd field,
                # remove "<" and ">" chars:
                job_id = stdout.split()[1].strip('<>')
            #endif
        except :
            self.logger.error( 'could not extract job id from stdout:' )
            self.logger.error( stdout )
            raise
        #endtry
    
        # help text:
        infotext.append(  '\n' )
        infotext.append(  'To manage Slurm jobs:\n' )
        infotext.append(  '\n' )
        infotext.append(  '  sbatch %s   # submit to queue\n' % jbfile )
        infotext.append(  '  squeue               # list your current jobs\n' )
        infotext.append(  '  scancel %s           # kill this job\n' % job_id )
        infotext.append(  '\n' )

        # write to loging system:
        for line in infotext :
            self.logger.info( line.rstrip() )
        #endfor

        # log file for submission:
        job_info = bname+'.info'
        # write to file:
        f = open( job_info, 'w' )
        f.writelines(infotext)
        f.close()

    #enddef Submit

#endclass UtopyaJobScriptBatchSlurm
    
    
# ======================================================================
# ===
# === jobs submitted to LoadLeveler queue
# ===
# ======================================================================


class UtopyaJobScriptBatchLoadLeveler( UtopyaJobScriptBatch ) :

    """
    Class to create job for submission to LoadLeveler batch system.
    
    Example of job options::
        
        #@ name = myjob
        #@ output = myjob.out
        #@ error = myjob.err
        #@ queue
    
    Example of usage::
        
        # init object:
        jb = utopya.UtopyaJobScriptBatchLoadLeveler()

        # obtain line with job options from rcfile:
        options = jb.GetOptionsRc( 'UtopyaJobScriptBatchLoadLeveler.rc', \\
                                        rcbase='appl', \\
                                        env={'name':'myjob'} )
        
        # fill script lines:
        lines = []
        lines.append( '#! /usr/bin/env python\\n' )
        lines.append( '\\n' )
        lines.append( options )
        lines.append( '\\n' )
        lines.append( '# do something:\\n' )
        lines.append( 'print( "boe!" )\\n' )
        lines.append( '\\n' )
        
        # write:
        jb.Create( 'myjob.jb', lines )
        
        # submit:
        jb.Submit( 'myjob.jb' )
        
    See also man pages of batch job commands:
            
    * :manpage:`llsubmit`
    * :manpage:`llq`
    * :manpage:`llcancel`

    """
    
    def GetOptionsRc( self, rcfile, rcbase='', env={} ) :
    
        """
        Return str line (with newline characters) with job options based on rcfile settings.
        
        The rcfile settings should start with '[<rcbase>.]batch.loadleveler',
        where the rcbase might be empty or equal to '*' for default settings.
        
        Example settings for rcbase 'appl'::
            
            ! job format for this application:
            appl.batch.loadleveler.format            :  loadleveler_format
            ! which keywords:
            appl.batch.loadleveler.options           :  name output error queue
            ! values:
            appl.batch.loadleveler.option.name       :  name myjob
            appl.batch.loadleveler.option.output     :  output %(name).out
            appl.batch.loadleveler.option.error      :  error %(name).err
            appl.batch.loadleveler.option.queue      :  queue
    
            ! Define format of batch options, e.g.:
            !   #@ key = value
            loadleveler_format.comment       :  #
            loadleveler_format.prefix        :  @
            loadleveler_format.arg           :  ''
            loadleveler_format.assign        :  ' = '
            loadleveler_format.template      :  %(key)
            loadleveler_format.envtemplate   :  %(env:key)

        This will return the following job options as a str with newline characters::

            #@ name = myjob
            #@ output = myjob.out
            #@ error = myjob.err
            #@ queue

        """
        
        # call parent:
        lines = UtopyaJobScriptBatch.GetOptionsRc( self, rcfile, 'batch.loadleveler', rcbase=rcbase, env=env )
        
        # ok
        return lines
        
    #enddef
    
    # *
    
    def Submit( self, jbfile ) :
    
        """
        Submit job file.
        Information on job id and commands to follow and cancel the job
        are written to a file with the same name but extension '.info' .
        """
        
        # modules:
        import sys
        import os
        import subprocess

        # basename for scripts etc is name of rcfile minus extension:
        bname,ext = os.path.splitext(jbfile)    

        # setup command line, last argument is script:
        command = 'llsubmit %s' % jbfile
        # info ...
        self.logger.info( 'command: %s' % command )

        # init submission info file:
        infotext = []
        infotext.append(  '\n' )

        # call submit command, trap errors:
        try:
            # submit; redirect errors to standard output:
            p = subprocess.Popen( command.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
        except :
            self.logger.error( sys.exc_info()[1] )
            self.logger.error( 'from command : %s' % command )
            raise Exception
        #endtry
        # extract:
        outlines = p.stdout.readlines()
        # add to help info message:
        infotext = infotext + outlines
        # extract job id from last line:
        #   llsubmit: The job "c1a0303.4290133" with 3 job steps has been submitted.
        firstwords = 'llsubmit: The job'
        lastline = outlines[-1]
        if lastline.startswith(firstwords) :
            job_id = lastline.lstrip(firstwords).split()[0].replace('"','')
        else :
            job_id = '<job-id>'
        #endif
    
        # help text:
        infotext.append(  '\n' )
        infotext.append(  'To manage LoadLeveler jobs:\n' )
        infotext.append(  '\n' )
        infotext.append(  '  llq [-u ${USER}]         # list [your] current jobs\n' )
        infotext.append(  '  llq %s             # list this job\n' % job_id )
        infotext.append(  '  llcancel %s        # kill this job\n' % job_id )
        infotext.append(  '\n' )

        # write to loging system:
        for line in infotext :
            self.logger.info( line.rstrip() )
        #endfor

        # log file for submission:
        job_info = bname+'.info'
        # write to file:
        f = open( job_info, 'w' )
        f.writelines(infotext)
        f.close()

    #enddef Submit

#endclass UtopyaJobScriptBatchLoadLeveler
    
    
# ======================================================================
# ===
# === end
# ===
# ======================================================================
    

