
#-------------------------------------------------
# help
#-------------------------------------------------

"""
***********************
``utopya_tools`` module
***********************

The ``utopya_tools`` module provides access to a number
of general purpose methods and classes.

Overview
========

The following methods and classes are provided:
  
* A method is provided to import classes from a module using a path-like name:
  
  * :py:meth:`ImportClass`
  
* Two methods are provided to create directories:

  * :py:meth:`CreateFilePath`
  * :py:meth:`CreateDirs`

* The following tools facilitate (re)writting of text files:
  
  * :py:meth:`WriteTextFile`
  * :py:meth:`UpdateTextFile`
  
* A class is provided to call a subprocess with more informative exception handling:
    
  * :py:class:`Call`


Methods and classes
===================
"""

#-------------------------------------------------
# routines
#-------------------------------------------------

def ImportClass( moduleclass, logger=None ) :

    """
    Import class from a module. Specify the `moduleclass` as a string::

        [path/]module.classname

    The path to the module is optional; if not specified, the module
    should be available on the search path.
    Example::

        mycls = ImportClass( '/path/to/utopya/vtrunk/py/utopya.UtopyaBase' )
        
    Optionally provide a 'logger' instance derived from the :py:class:`logging.Logger` 
    class to issue messages.

    """
    
    # modules:
    import os
    import sys
    import time
    import logging
    import imp
    #import importlib
    
    # get logger object if not provided:
    if logger is None : 
        logging.basicConfig()
        logger = logging.getLogger()
    #endif

    # check ..
    if '.' not in moduleclass :
        logger.error( 'Argument `moduleclass` should have form: <modulename>.<classname>' )
        raise ValueError
    #endif

    # split at last dot:
    modulename,classname = moduleclass.rsplit('.',1)

    # module name might include path:
    modulepath,modulename = os.path.split( modulename )
    # extend path ?
    if len(modulepath) > 0 : sys.path.insert(0,modulepath)
    
    # import this module with the specified class:
    try :
        #~ this sometimes fails if module file was just copied to a newly created directory ...
        #mod = __import__( modulename, fromlist=[classname] )
        #~ this seems to work for the moment ...
        fp,pathname,description = imp.find_module( modulename, sys.path )
        mod = imp.load_module( modulename, fp, pathname, description )
    except ImportError :
        logger.error( 'Could not import module "%s" ...' % (modulename) )
        logger.error( 'Search path:' )
        for p in sys.path :
            logger.error( '  %s' % p )
            modulefile = os.path.join(p,modulename+'.py')
            if os.path.isfile( modulefile ) : logger.error( '    FOUND: %s' % modulefile )
        #endfor
        logger.error( 'Error from importing:' )
        raise
    except :
        raise
    #endtry

    # create object, prepare for errors:
    try:
        cls = mod.__dict__[classname]
    except KeyError:
        logger.error( "Class '%s' not defined in module '%s'" % (classname,modulename) )
        sys.exit()
    except:
        logger.error( "Unknown error importing '%s' class from '%s' module" % (classname,modulename) )
        sys.exit()
    #endtry

    # remove extra path if necessary:
    if len(modulepath) > 0 : sys.path.pop(0)

    # ok
    return cls

#enddef ImportClass


# ***


def CreateFilePath( filename, logger=None ) :

    """
    Create directory in provided 'filename' if not present yet.
    Optionally provide a 'logger' instance derived from the 
    :py:class:`logging.Logger` class to issue messages.
    """

    # modules:
    import os

    # directory name:
    dirname = os.path.dirname(filename)
    # defined ?
    if len(dirname) > 0 :
        # not present yet ?
        if not os.path.isdir(dirname) :
            os.makedirs( dirname )
        #endif
    #endif

#enddef CreateFilePath

# *

def CreateDirs( dirname, renew=False, logger=None ) :

    """
    Create directory (including subdirs) if not present yet. 
    If renew flag is set, remove existing version first.
    Optionally provide a 'logger' instance derived from the
    :py:class:`logging.Logger` class to issue messages.
    """

    # external:
    import os
    import shutil
    import logging

    # get logger object if not provided:
    if logger is None : logger = logging.getLogger()

    # already present ?
    if os.path.isdir( dirname ) :
        # remove existing directory ?
        if renew :
            # info ...
            logger.info( 'remove existing %s ...' % dirname )
            # use the chainsaw:
            shutil.rmtree( dirname )
        #endif
    #endif

    # not present (anymore) ?
    if not os.path.isdir( dirname ) :
        # info ...
        logger.info( 'create new directory %s ...' % dirname )
        # create directory including parents if necesary:
        os.makedirs(dirname)
    #endif

#enddef CreateDirs


# ***


def WriteTextFile( filename, text, logger=None ) :

    """
    Write a text file.

    Arguments:

    * filename   : target file name
    * text       : list of strings, a line should end with the newline character '\\\\n'
    
    Optional arguments:
      
    * logger : instance of :py:class:`logging.Logger` class.
    """

    # external:
    import os
    import logging

    # get logger object if not provided:
    if logger is None : logger = logging.getLogger()

    ## info ...
    #logger.debug( '    write %s ...' % filename )

    # create path if necessary:
    CreateFilePath( filename, logger=logger )

    # write new text:
    f = open( filename, 'w' )
    f.writelines( text )
    f.close()

#enddef WriteTextFile


# ***


def UpdateTextFile( filename, newtext, logger=None ) :

    """
    Replace a file by a new text if the later differs
    from the current content.

    Arguments:

    * filename   : target file name
    * newtext    : list of strings, a line should end with the newline character '\\\\n'
    
    Optional arguments:
      
    * logger : instance of :py:class:`logging.Logger` class.
    """

    # external:
    import os
    import logging

    # get logger object if not provided:
    if logger is None : logger = logging.getLogger()

    ## info ...
    #logger.info( '  update %s ...' % filename )
    # file exists alread?
    if os.path.exists(filename) :
        # read current content:
        f = open( filename, 'r' )
        oldtext = f.readlines()
        f.close()
        # differences ?
        rewrite = newtext != oldtext
        ## info ...
        #for iline in range(len(oldtext)) :
        #    if iline < len(newtext) :
        #        if oldtext[iline] != newtext[iline] :
        #            logger.info( '    first different lines:' )
        #            logger.info( '      old: %s' % oldtext[iline] )
        #            logger.info( '      new: %s' % newtext[iline] )
        #            break
        #        #endif
        #    #endif
        ##endfor
    else :
        # no file yet, always rewrite:
        rewrite = True
    #endif

    # write file ?
    if rewrite :
        # for info message:
        stat = 'replace'
        # write new text:
        f = open( filename, 'w' )
        f.writelines( newtext )
        f.close()
    else :
        # for info message:
        stat = 'keep'
    #endif

#enddef UpdateTextFile
    
# *

def ReadNcVariable( filename, varname, logger=None ) :

    """
    Return data and attributes from netcdf variable.
    
    Arguments:
      filename   :  input file (netcdf)
      varname    :  (path to) variable
    
    Return values:
      data     : array with data values
      attrs    : dictionairy with the variable attributes
    """

    # modules:
    import os
    import netCDF4
    import logging

    # get logger object if not provided:
    if logger is None : logger = logging.getLogger()

    # check ..
    if not os.path.isfile(filename) :
        logger.error( 'file not found: %s' % filename )
        raise Exception
    #endif

    # open:
    ncid = netCDF4.Dataset( filename, 'r' )

    # main group:
    gid = ncid
    vname = varname
    vpath = '/'
    # groups?
    while '/' in vname :
        # split:
        gname,vname = vname.split('/',1)
        # check ...
        if gname not in gid.groups.keys() :
            logger.error( 'group "%s" not found in: %s' % (gname,filename) )
            raise Exception
        #endif
        # new group:
        gid = gid.groups[gname]
        # extend:
        vpath = vpath+gname+'/'
    #endwhile
    # check ...
    if vname not in gid.variables.keys() :
        logger.error( 'variable "%s" not found in group "%s" in: %s' % (vname,vpath,filename) )
        raise Exception
    #endif
    # select:
    varid = gid.variables[vname]
    # read:
    data = varid[:]
    # attributes:
    attrs = {}
    for key in varid.ncattrs() :
        attrs[key] = varid.getncattr(key)
    #endfor

    # close:
    ncid.close()
    
    # ok
    return data,attrs

#enddef ReadNcVariable


#-------------------------------------------------
# subprocess
#-------------------------------------------------

class Call( object ) :

    """
    Class to call a subprocess and trap exceptions.
    
    This routine differs from those provided by the standard 
    :py:mod:`subprocess` module in the following ways.
    
    * In case the command could not be executed a :py:exc:`CallingError` is raised.
      The original command as well as the error message are included in the exception.
    * If the command returns with a non-zero exit status a :py:exc:`StatusError` is raised.
    * Eventually, send standard output and standard error to the supplied logger.
      The standard output is written with loglevel `logging.INFO` if verbose=True,
      and otherwise with loglevel 'logging.DEBUG'.
      The standard error is written with level `logging.ERROR`.
      All lines written by the logger are preceded by the optional indent.
      With the 'err2out' flag enable, std.error is written to std.output,
      which is sometimes useful to maintain the order in which messages are issued.
    
    The command and optional keyword arguments are passed to the
    :py:class:`subprocess.Popen` class.
    The opional input (list of str objects) will be passed to
    its :py:meth:`subprocess.Popen.communicate` method. 
    
    Returns an object with attributes:

    * retcode   (integer)
    * stdout    (list of strings)
    * stderr    (list of strings)

    Example of usage::

        # modules:
        import logging
        
        # get logger object:
        logger = logging.getLogger()

        # list file, log result outside call:
        try :
            p = Call( ['/bin/ls','-l'] )
        except CallingError as err :
            logger.error( str(err) )
            raise Exception
        except StatusError as err :
            for line in err.stderr : logging.error(line)
            logger.error( str(err) )
            raise Exception
        except Exception as err :
            logger.error( 'unknown error : %s' % err )
            raise Exception
        #endtry
        # display output:
        for line in p.stdout : logger.info(line)
        for line in p.stderr : logger.error(line)
        
        # idem, but exception handling and logging is done in the call:
        p = Call( ['/bin/ls','-l'], logger=logger )

    
    """

    def __init__( self, command, input=None, logger=None, verbose=True, indent='', err2out=False, **kwargs ) :

        # external:
        import subprocess
        
        # init outputs:
        self.stdout = []
        self.stderr = []
        self.returncode = -999
        
        # error destination:
        stderr_dest = subprocess.PIPE
        if err2out : stderr_dest = subprocess.STDOUT
        
        # call suprocess, catch calling errors:
        try :
            p = subprocess.Popen( command,
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=stderr_dest,
                                  **kwargs )
        except Exception as err :
            if logger is not None : logger.error( indent+str(err) )
            raise CallingError( command, err )
        #endtry
        
        # also input provided?
        if input is not None :

            # send input if necessary ;
            # wait for process to terminate, receive standard output and error:
            stdout_str,stderr_str = p.communicate( input=input )

            # remove final end-of-lines:
            stdout_str = stdout_str.strip('\n')
            stderr_str = stderr_str.strip('\n')

            # store lists:
            if len(stdout_str) > 0 :
                self.stdout = stdout_str.split('\n')
            else :
                self.stdout = []
            #endif
            if len(stderr_str) > 0 :
                self.stderr = stderr_str.split('\n')
            else :
                self.stderr = []
            #endif

            # write standard output and error:
            if logger is not None :
                # output messages:
                if verbose :
                    for line in self.stdout : logger.info( indent+line )
                else :
                    for line in self.stdout : logger.debug( indent+line )
                #endif
                # output errors:
                for line in self.stderr : logger.error( indent+line )
            else :
                # output messages:
                for line in self.stdout : print( indent+line )
                # output errors:
                for line in self.stderr : print( indent+line )
            #endif
            
        else :

            # wait for process to terminate:
            while p.poll() == None :
                # read next line from standard output:
                line = p.stdout.readline()
                # filled?
                if len(line) > 0 :
                    # convert:
                    line = line.strip().decode('utf-8')
                    # show:
                    if logger is not None :
                        # output message:
                        if verbose :
                            logger.info( indent+line )
                        else :
                            logger.debug( indent+line )
                        #endif
                    else :
                        # output message:
                        print( indent+line )
                    #endif
                    # store:
                    self.stdout.append( line )
                #endif
            #endwhile
            
            # read remaing lines:
            for line in p.stdout.readlines() :
                # convert:
                line = line.strip().decode('utf-8')
                # show:
                if logger is not None :
                    # send message to logger:
                    if verbose :
                        logger.info( indent+line )
                    else :
                        logger.debug( indent+line )
                    #endif
                else :
                    # print:
                    print( indent+line )
                #endif
                # store:
                self.stdout.append( line )
            #endfor
            
            # error lines?
            if p.stderr is not None :
                # read error lines:
                for line in p.stderr.readlines() :
                    # convert:
                    line = line.strip().decode('utf-8')
                    # output message:
                    if logger is not None :
                        logger.error( indent+line )
                    else :
                        print( indent+line )
                    #endif
                    # store:
                    self.stderr.append( line )
                #endfor
            #endif

        #endif # input provided?

        # store return code:
        self.returncode = p.returncode        
        # error ?
        if self.returncode != 0 :
            # return status error:
            raise StatusError( command, self.returncode, self.stdout, self.stderr )
        #endif

    #enddef __init__

#endclass Call

# *

class CallingError( Exception ) :

    """
    This exception is raised when an attempt to call a command fails.
    The following attributes are defined:

    * command   : str object or list of str objects with the command that failed
    * strerror  : error message raised by command
    
    Use the 'str()' method to obtain a nicely formatted message
    including the command and the error message.
    
    Example of usage::
      
        try :
            ...
        except CallingError as err :
            print( str(err) )        
    """

    def __init__( self, command, strerror ) :
        """
        Store exception attributes.
        """
        # store:
        self.command = command
        self.strerror = strerror
    #enddef __init__
    
    def __str__( self ) :
        """
        Return formatted erorr message.
        """
        # error message:
        return 'Error from calling "%s" : %s' % (self.command,self.strerror)
    #enddef __str__
    
#endclass CallingError

# *

class StatusError( Exception ) :

    """
    This exception is raised when a called command returns a non-zero exit status.
    The following attributes are defined:

    * command   : str object (list of) with original command
    * returncode : integer return status
    * stdout     : str list
    * stderr     : str list
    
    Use the 'str()' method to obtain a nicely formatted message
    including the command and the return code.
    
    Example of usage::
      
        try :
            ...
        except StatusError as err :
            print( str(err) )
    """

    def __init__( self, command, returncode, stdout, stderr ) :
        """
        Store exception attributes.
        """
        # store:
        self.command = command
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
    #enddef __init__
    
    def __str__( self ) :
        """
        Format error message.
        """
        # error message:
        return 'Call "%s" returned non-zero status %i' % (str(self.command),self.returncode)
    #enddef __str__
    
#endclass StatusError

#-------------------------------------------------
# end
#-------------------------------------------------

