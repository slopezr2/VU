
#-------------------------------------------------
# help
#-------------------------------------------------

"""
*************************
``utopya_archive`` module
*************************

The :py:mod:`utopya_archive` module provides classes to collect
and store output files.
All classes are accesible through the top level :py:mod:`utopya` module.

* The :py:class:`.UtopyaArchive` class creates an instance of a top level
  object that, based on rcfile settings, performs all archiving tasks
  requested in the settings.
  
The configuration of the objects is done using rcfile settings.


Class hierarchy
===============

The classes provided by this module have been derived with the following hierchy:

* :py:class:`.UtopyaRc`

  * :py:class:`.UtopyaArchive`

Classes
=======

"""


#-------------------------------------------------
# modules
#-------------------------------------------------

import utopya_rc


#-------------------------------------------------
# builder
#-------------------------------------------------

class UtopyaArchive( utopya_rc.UtopyaRc ) :

    """
    Base class for objects that collect and store output files.
    
    The initialization arguments look like::
    
        UtopyaArchive( 'settings.rc', 'model.archive', 
                          env={ 'WORKDIR' : '/scratch/yourname/work' } )
                          
    The first argument defines the name of the rcfile with settings,
    and the second is the base name of the settings that will be used.
    The optional 'env' dictionairy could be used to pass environment variables
    that are expanded in the rcfile.
    
    The first configuration in the settings should be the definition of a 
    collection list. This is just a list of keywords with nicknames for
    the file collections to be archived::
    
        ! list of output collections to be made:
        model.archive.collections       :  scripts output
    
    where 'model.archive' is the 'rcbase' passed as argument.
    In this example, two collections are to be made, one containing the settings
    for the experiment, and one containing the output.
    For each of the collections, the details should be specified by keywords
    that start with the 'rcbase' extended with the collection name.
    An example for the 'scripts' collection is::

        ! define collection details:
        !~ base directory:
        model.archive.scripts.basedir   :  ${WORKDIR}
        !~ target file:
        model.archive.scripts.target    :  /archive/model/scripts.tar.gz
        !~ filename filter:
        model.archive.scripts.files     :  py/*.py
        
    These settings are:
    
    * The 'basedir' specifies the location from where files should be collected for archiving.
    * The 'target' specifies the name of the collection file to be created.
      Typically this is a '.tar.gz' file that collects and compresses a set of files.
      The target could be a remote location, with the description supported by
      the :py:mod:`gss` module.
    * The 'files' settings specifies a space-seperated list of filename filters,
      for example to collect all python scripts.

    The target specifications could have different extensions that define 
    the collecting and compression program. 
    Currenlty supported extensions are::
    
        .tar .tar.gz .tgz .zip

    """
    
    def __init__( self, rcfile, rcbase='', env={} ) :
    
        """
        Initialize archiver object.
        """
        
        # modules:
        import os
        
        # init base object:
        utopya_rc.UtopyaRc.__init__( self, rcfile=rcfile, rcbase=rcbase, env=env )
        
        # info ...
        self.logger.info( 'archive ...' )
        self.logger.info( '  settings from : %s' % rcfile )
        self.logger.info( '  base key      : %s' % rcbase )
        if len(env) > 0 :
            self.logger.info( '  environment:' )
            for key in env.keys() :
                self.logger.info( '    %s = %s' % (key,env[key]) )
            #endfor
        #endif
        
        # read collection list:
        collections = self.GetSetting( 'collections' ).split()
        # info ...
        self.logger.info( '  loop over collections ...' )
        # loop:
        for collection in collections :

            # info ...
            self.logger.info( '    collection "%s" ...' % collection )

            # base directory to collect output:
            basedir = self.GetSetting( collection+'.basedir' )
            # change base dir?
            if len(basedir) > 0 :
                # info ...
                self.logger.info( '      change to %s ...' % basedir )
                # save current:
                owd = os.getcwd()
                # change:
                os.chdir( basedir)
            #endif # basedir defined?
            
            # target file:
            target = self.GetSetting( collection+'.target' )
            # source files:
            files = self.GetSetting( collection+'.files' )
            # create target:
            self.CreateAndStore( target, files, indent='      ' )

            # change base dir?
            if len(basedir) > 0 :
                # back:
                os.chdir( owd )
            #endif # basedir defined?
            
        #endfor # collections
        
    #enddef __init__
    
    # *
    
    def CreateAndStore( self, target, files, indent='' ) :
    
        """
        Create target archive with specified files.

        The 'target' filename might be at a remote location;
        the format of the target should then be supported by the
        :py:mod:`gss` module.
        In that case, first a local target file will be created,
        which is after completion stored to the remote location.
        
        The 'files' argument is a string with filename patterns,
        which are passed to the :py:meth:`Create` method.
        """
        
        # modules:
        import os
        import sys
        
        # gss path?
        if ':' in target :
            # info ...
            self.logger.info( indent+'target is at remote location, create local target first ...' )
            # temporary local file:
            tfile = 'tmp.%i.%s' % (os.getpid(),os.path.basename(target))
            # need to store this file later on:
            store = True
        else :
            # target is already local...
            tfile = target
            # no need to store:
            store = False
        #endif
        
        # create local target:
        self.Create( tfile, files, indent=indent )
        
        # store if ncessary:
        if store :
            # info ...
            self.logger.info( indent+'store %s ..' % target )
            # tools:
            import gss
            # copy:
            gss.Copy( tfile, target, verbose=True )
            # cleanup:
            os.remove( tfile )
        #endif
        
    #enddef CreateAndStore
    
    # *
    
    def Create( self, target, files, indent='' ) :
    
        """
        Create target archive with specified files.
        
        Arguments:

        * 'target' : name of archive to be created, for example 'settings.tar.gz'
        * 'files' : line with filename matching patterns, for example 'rc data/\*.txt'
        """
        
        # modules:
        import sys
        import os
        import subprocess
        
        # tools:
        import utopya_tools

        # info ...
        self.logger.info( indent+'create %s ..' % target )
        self.logger.info( indent+'files to be included:' )
        for f in files.split() : self.logger.info( indent+'  %s' % f )
        
        # target directory:
        tdir = os.path.dirname( target )
        # create?
        if (len(tdir) > 0) and (not os.path.isdir(tdir)) :
            # info ...
            self.logger.info( indent+'  create %s directory ...' % tdir )
            # create:
            os.makedirs( tdir )
        #endif

        # run shell commands (needed to expand filename patterns):
        if target.endswith('.tar') :
            # use 'tar' as archiver command, verbose:
            command = 'tar c -v -f %s %s' % (target,files)
        #
        elif target.endswith('.tar.gz') or target.endswith('.tgz') :
            # use 'tar' as archiver command, verbose, gzip the result:
            command = 'tar zc -v -f %s %s' % (target,files)
        #
        elif target.endswith('.zip') :
            # use 'zip' as archiver command, 
            # always recursive to ensure collection of directories:
            command = 'zip -r %s %s' % (target,files)
        #
        else :
            self.logger.error( 'unsupported extension in target file: %s' % target )
            raise Exception
        #endif

        # info ...
        self.logger.info( indent+'command: %s' % command )
        # call subprocess, trap errors:
        self.Call( command, shell=True, indent=indent+'  ' )
        # info ...
        self.logger.info( indent+'(note that informative messages might have beeen displayed as errors ...)' )
        
    #enddef Create

#endclass UtopyaArchive


#-------------------------------------------------
# end
#-------------------------------------------------

