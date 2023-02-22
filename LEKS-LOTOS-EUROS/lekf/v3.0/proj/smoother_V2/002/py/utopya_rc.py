
#-------------------------------------------------
# help
#-------------------------------------------------

"""
********************
``utopya_rc`` module
********************

The ``utopya_rc`` module holds the base class
for objects that have settings defined in a rcfile.

Class hierarchy
===============

The classes are defined according to following hierarchy:

* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaRc`

Classes
=======


"""

#-------------------------------------------------
# modules
#-------------------------------------------------

from utopya_base import UtopyaBase


#-------------------------------------------------
# rcfile class
#-------------------------------------------------

class UtopyaRc( UtopyaBase ) :

    """
    Base class for UTOPyA objects configured by a text file called the 'rcfile'.
    See the documentation of the :py:mod:`rc` module for details
    on the rcfile formating and preprocessing on reading.

    The settings are read on initialization if arguments 'rcfile' and 'rcbase' are provided,
    or later on by a call to :py:meth:`InitRc`.
    
    Use the :py:meth:`GetSetting` method to obtain the value of a setting.
    For example, if the rcfile 'test.rc' has content::
      
        ! sample setting:
        data.value    :  10
        
    this could be retrieved using::
      
        uto = UtopyaRc( 'test.rc' )
        n = uto.GetSetting( 'data.value', 'int' )
        
    The optional 'rcbase' defines the first part of the settings,
    and can be used to define groups of similar settings.
    For example, if the rcfile contains settings for 'Model-A' and 'Model-B'::
    
        ! sample setting:
        Model-A.data.value    :  10
        Model-B.data.value    :  20

    then the settings for 'Model-A' are selected using::
              
        uto = UtopyaRc( 'test.rc', rcbase='Model-A' )
        n = uto.GetSetting( 'data.value', 'int' )

    Settings that apply to all bases could be defined by a template::
      
        ! sample setting:
        *.data.value    :  10
    
    The following rcfile settings should be present
    to initialize the logging sytem::
        
       ! Enable debug messages ?
       ! If True, the logging level is 'DEBUG', otherwise 'INFO' :
       [<rcbase>.]logging.verbose       :  True
       
       ! logging format, empty for default:
       [<rcbase>.]logging.format        :  %(message)s

    In summary, the following methods are provided with this class:
        
    * :py:meth:`InitRc` : read rcfile if not done on initialization;
    * :py:meth:`GetSetting` : extract value of setting;
    * :py:meth:`ImportClass` : import class specified by an rcfile setting.
        
    """
    
    def __init__( self, rcfile=None, rcbase='', env={} ) :
    
        """
        Initialize UTOPyA RC class.
        """
        
        # init base class:
        UtopyaBase.__init__( self )
        
        # initialize rcfile if provided
        if rcfile is None :
            # check ...
            if rcbase != '' :
                self.logger( 'rcbase provided without rcfile' )
                raise Exception
            #endif
            # not yet ..
            self.rcf    = None
            self.rcfile = None
            self.rcbase = None
        else :
            # define rcfile object, read settings:
            self.InitRc( rcfile, rcbase, env=env )
        #endif

    #enddef __init__
    
    # *
    
    def InitRc( self, rcfile, rcbase='', env={} ) :
    
        """
        Initialize rcfile if not done on initialization.
        This also initializes the logging based on the rcfile settings.
        """
        
        # tools:
        import rc
        
        # define rcfile object, read settings:
        self.rcf = rc.RcFile( rcfile, env=env )
        # store arguments:
        self.rcfile = rcfile
        self.rcbase = rcbase

        # setup logging given rc params:
        self._InitLogging()
        
    #enddef
        
    
    # *
    
    def __del__( self ) :
    
        """
        Destructor for UTOPyA RC class.
        """
        
        # delete rcfile object:
        if hasattr(self,'rcf') : del self.rcf
        
        # call base desctuctor:
        UtopyaBase.__del__( self )
        
    #enddef __del__
    
    # *
    
    def GetSetting( self, name, totype='', **kwargs ) :
    
        """
        Return value of the named entry in the rcfile::
            
            <name>       :  <value>
            
        If at initialization an optional 'rcbase' was provided,
        the entry that is searched for is prefixed by this base,
        or is a generic version::
            
            <rcbase>.<name>    :  <value>
            *.<name>           :  <value>
  
        Optional keyword arguments are pased to the :py:meth:`.RcFile.get` method.
        For example, a keyword argument 'default' could be
        used to set the return value if none of the keys is found.

        Note that double dots '..' in a name are replaced by a single '.',
        this is useful when the provided name is combined from some
        underlying parts and some of these start or end with a dot already.
        """
        
        # check ...
        if self.rcf is None :
            self.logger.error( 'no rcfile defined for object' )
            raise Exception
        #endif
        
        # form key:
        if self.rcbase is None :
            xkey = name
        else :
            xkey = self.rcbase+'.'+name
        #endif

        # remove double and leading seperations:
        xkey = xkey.replace('..','.').lstrip('.')
        
        # add all generic versions: 'a.b.c.c', '*.b.c.d', '*.c.d', ...
        keys = [xkey,'*.'+xkey]
        while '.' in xkey :
            bkey,xkey = xkey.split('.',1)
            keys.append( '*.'+xkey )
        #endwhile

        # flag:
        found = False
        # loop:
        for key in keys :
            # present ?
            if key in self.rcf.keys() :
                # read from rcfile:
                value = self.rcf.get( key, totype=totype, **kwargs )
                # set flag:
                found = True
            #endif
            # leave ?
            if found : break
        #endif
        # not found ?
        if not found :
            # default specified ?
            if 'default' in kwargs.keys() :
                value = kwargs['default']
            else :
                self.logger.error( 'none of the following key(s) found in: "%s"' % self.rcfile )
                for key in keys :
                    self.logger.error( '  %s' % key )
                #endfor
                raise Exception
            #endif
        #endif
        
        # ok:
        return value
        
    #enddef GetSetting
    
    # *
    
    def WriteSettings( self, filename ) :
    
        """
        Write evaluated settings to file.
        """
        
        # check ...
        if self.rcf is None :
            self.logger.error( 'no rcfile defined for object' )
            raise Exception
        #endif
        
        # write:
        import os
        self.rcf.WriteFile( filename )
        
    #enddef WriteSettings

    # *
    
    def _InitLogging( self ) :
    
        """
        Internal routine to setup logging properties from rcfile settings.
        """
        
        # modules:
        import logging
        
        # verbose messages ?
        verbose = self.GetSetting( 'logging.verbose', totype='bool', default=False )
        # set level keyword:
        if verbose :
            level = 'debug'
        else :
            level = 'info'
        #endif

        # format description:
        format = self.GetSetting( 'logging.format', totype='str', default='[%(levelname)-8s] %(message)s' )
        # remove quotes:
        if format.startswith("'") :
            if not format.endswith("'") :
                self.logger.error( "format not properly enclosed by quotes: %s" % format )
                raise Exception
            #endif
        #endif
        # remove quotes:
        if format.startswith('"') :
            if not format.endswith('"') :
                self.logger.error( 'format not properly enclosed by double quotes: %s' % format )
                raise Exception
            #endif
        #endif
        # reset empty to default:
        if len(format) == 0 : format = None

        # set properties:
        self.LoggerSet( format=format, level=level )

    #enddef _InitLogging

    # *
    
    def ImportClass( self, name ) :

        """
        Import class from a module, 
        where the names are specified in the rcfile.
        
        For example, to import 'MyClass' from the module 'MyMod'
        defined in the file '/path/to/MyMod.py', use the following setting::
            
            [<rcbase>.]name      :  /path/to/MyMod.MyClass
        """
        
        # modules:
        import os
        import sys
        
        # tools:
        import utopya_tools
        
        # name of model module and class:
        moduleclass = self.GetSetting( name )
        
        # use tool function, provide current logger:
        cls = utopya_tools.ImportClass( moduleclass, logger=self.logger )
        
        # ok
        return cls
        
    #enddef ImportClass

#endclass UtopyaRc



#-------------------------------------------------
# end
#-------------------------------------------------


