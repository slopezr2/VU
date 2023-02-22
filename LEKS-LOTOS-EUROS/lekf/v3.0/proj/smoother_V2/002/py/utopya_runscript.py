
#-------------------------------------------------
# help
#-------------------------------------------------

"""
***************************
``utopya_runscript`` module
***************************

The ``utopya_runscript`` module provides access to UTOPyA's run script classes.
Note that all entities could be imported directly from the main ``utopya`` module too.

Two classes are provided:
  
* The :py:class:`UtopyaRunScript` class could be used for a standard script without 
  specific arguments.
* The :py:class:`UtopyaRunScriptRc` class requires the name of an rcfile with settings
  as command line argument, and settings could be read directly from this class.

Class hierarchy
===============

The classes are defined according to following hierarchy:

* :py:class:`.UtopyaBase`

  * :py:class:`.UtopyaRunScript`

* :py:class:`.UtopyaRc`

  * :py:class:`.UtopyaRunScriptRc`

Classes
=======

"""

#-------------------------------------------------
# modules
#-------------------------------------------------

# tools:
from utopya_base import UtopyaBase
from utopya_rc   import UtopyaRc


#-------------------------------------------------
# Runscript
#-------------------------------------------------

class UtopyaRunScript( UtopyaBase ) :

    """
    Base class for script object, to be used in user defined scripts 
    to quickly perform some common tasks:
    
    * setup logging to put out informative messages;
    * define some default arguments;
    * parse the arguments.
    
    Example script::
        
        #! /usr/bin/env python
    
        # modules:
        import utopya
        
        # init script:
        utos = utopya.UtopyaRunScript()

        # setup default command line arguments:
        utos.ArgumentsSetup( description='UTOPyA test script' )

        # evaluate known arguments, store the other ;
        # might show help text and exit:
        args,xargs = utos.ArgumentsParse()

        # get logger object to shout something:
        logger = utos.LoggerGet()
        
        # do something:
        logger.info( 'boe!')
    
    """
    
    def __init__( self ) :
    
        """
        Performs the initialization of the ancestor class,
        in this mainly the setup of the logging system.
        """
        
        # setup base class:
        UtopyaBase.__init__( self )
        
    #enddef __init__
        
    # *
    
    def ArgumentsSetup( self, description=None ) :
    
        """
        Defines data object 'parser' of class :py:class:`argparse.ArgumentParser` .
        The parser is populated with default arguments for:
          
        * verbose option to shout not only info but also debug messages.
        
        """
        
        # modules:
        import argparse
        
        #: setup argument parser:
        self.parser = argparse.ArgumentParser( description=description )
        
        # messages:
        self.parser.add_argument( '-v', '--verbose', help='Show debug messages.', 
                                   action='store_true', dest='verbose' )
        
    #enddef ArgumentsSetup
    
    # *
    
    def ArgumentsAdd( self, *args, **kwargs ) :

        """
        Add a new command line argument.
        
        The arguments provided to this routine are passed directly
        to the :py:meth:`argparse.ArgumentParser.add_argument` routine of the 'parser' 
        object that was defined by the 'SetupArguments' routine.
        See the documentation of this underlying  method for useful arguments.
        
        Examples::
        
            # argument with value:
            ArgumentsAdd( '-n', '--name', help='name used for something',
                                   action='store', dest='name' )
            # flag:
            ArgumentsAdd( '-f', '--flag', help='enable something', 
                                   action='store_true', dest='flag' )

        After the call to :meth:`ArgumentsParse`, the arguments are available through::
        
            # parse arguments, return known and extra:
            args,xargs = ArgumentsParse()
            # show:
            print( 'name = ', args.name )
            print( 'flag = ', args.flag )
        
        """
        
        # pass settings:
        self.parser.add_argument( *args, **kwargs )
        
    #enddef ArgumentsAdd
    
    # *
    
    def ArgumentsParse( self ) :
    
        """
        Parse known arguments.
        
        The optional verbose setting is evaluated, and if necessary, the logging level is changed.
        
        In case the arguments included the '-h' or '--help' flag,
        an automatically created help tekst is displayed and the program is terminated.
        Otherwise, the method returns with the following arguments:
            
        * args    : namespace with defined arguments and values
        * xargs   : extra arguments
        """
        
        # modules:
        import logging

        # evaluate known arguments, store the other:
        args,xargs = self.parser.parse_known_args()
        
        # Explict request for debugging info ?
        # Note there is no way to reset the logging ...
        if args.verbose : self.LoggerSet( level='debug' )

        # ok
        return args,xargs
        
    #enddef ArgumentsParse
        
#endclass UtopyaRunScript
        

#-------------------------------------------------
# Runscript with rcfile settings
#-------------------------------------------------

class UtopyaRunScriptRc( UtopyaRc ) :

    """
    Base class for script object, to be used in user defined scripts 
    that need a rcfile argument. 
    The following tasks are performed by the class:
    
    * setup logging to put out informative messages;
    * define default arguments including rcfile name;
    * parse the arguments, and read the settings from the rcfile.
    
    Example script::
        
        #! /usr/bin/env python
    
        # modules:
        import utopya
        
        # init script, pass default rcfile settings:
        utos = utopya.UtopyaRunScriptRc()

        # setup default command line arguments:
        utos.ArgumentsSetup( description='UTOPyA test script',
                              rcfile='test.rc', rcbase='appl' )

        # evaluate known arguments, store the other ;
        # might show help text and exit:
        args,xargs = utos.ArgumentsParse()

        # get logger object to shout something:
        logger = utos.LoggerGet()
        
        # do something:
        logger.info( 'boe!')
    
    """
    
    def __init__( self ) :
    
        """
        Performs the initialization of the ancestor class,
        in this mainly the setup of the logging system.
        """
        
        # setup base class, currently without rcfile settings,
        # these will be added after parsing command line arguments:
        UtopyaRc.__init__( self )
        
    #enddef __init__
        
    # *
    
    def ArgumentsSetup( self, description=None, 
                          rcfile=None, rcbase='' ) :
    
        """
        Defines data object 'parser' of class :py:class:`argparse.ArgumentParser` .
        The parser is populated with default arguments for:
          
        * verbose option to shout not only info but also debug messages;            
        * name of the rcfile with settings (default to optional argument 'rcfile');
        * base key for relevant rcfile settings (default to optional argument 'rcbase').
        """
        
        # modules:
        import argparse
        
        # setup argument parser:
        self.parser = argparse.ArgumentParser( description=description )
        
        # messages:
        self.parser.add_argument( '-v', '--verbose', help='Show debug messages.', 
                                   action='store_true', dest='verbose' )

        # settings file:
        self.parser.add_argument( 'rcfile', help='Settings file (default "%s")' % rcfile, 
                                     nargs='?', default=rcfile )
        # prefix of settings:
        self.parser.add_argument( '-r', '--rcbase', help='Base of rcfile keys (default "%s")' % rcbase,
                                   action='store', dest='rcbase', default=rcbase )
        
    #enddef ArgumentsSetup
    
    # *
    
    def ArgumentsAdd( self, *args, **kwargs ) :

        """
        Add a new command line argument.
        
        The arguments provided to this routine are passed directly
        to the :py:meth:`argparse.ArgumentParser.add_argument` routine of the 'parser' 
        object that was defined by the 'SetupArguments' routine.
        See the documentation of this underlying  method for useful arguments.
        
        Examples::
        
            # argument with value:
            ArgumentsAdd( '-n', '--name', help='name used for something',
                                   action='store', dest='name' )
            # flag:
            ArgumentsAdd( '-f', '--flag', help='enable something', 
                                   action='store_true', dest='flag' )

        After the call to :meth:`ArgumentsParse`, the arguments are available through::
        
            # parse arguments, return known and extra:
            args,xargs = ArgumentsParse()
            # show:
            print( 'name = ', args.name )
            print( 'flag = ', args.flag )
        
        """
        
        # pass settings:
        self.parser.add_argument( *args, **kwargs )
        
    #enddef ArgumentsAdd
        
    # *
    
    def ArgumentsParse( self, env={} ) :
    
        """
        Parse known command line arguments.
        
        In case the arguments included the '-h' or '--help' flag,
        an automatically created help tekst is displayed and the program is terminated.
        Otherwise, the method returns with the following arguments:
            
        * args    : namespace with defined arguments and values
        * xargs   : extra arguments
        
        The 'rcfile' and 'rcbase' passed as command line arguments,
        and eventually the optional 'env' dictionairy passed to this method,
        are used to initalize an attribute 'rcf' of the :py:class:`.RcFile` class 
        from which settings could be read.
        
        Necessary settings are those that define the logging settings,
        see the description of the :py:class:`.UtopyaRc` class.
        
        """
        
        # modules:
        import logging

        # evaluate known arguments, store the other:
        args,xargs = self.parser.parse_known_args()
        
        # Explict request for debugging info ?
        # Note there is no way to reset the logging ...
        if args.verbose : self.LoggerSet( level='debug' )
        
        # init rcfile settings:
        self.InitRc( args.rcfile, rcbase=args.rcbase, env=env ) 
        
        # ok
        return args,xargs
        
    #enddef ArgumentsParse
        
#endclass UtopyaRunScriptRc


#-------------------------------------------------
# end
#-------------------------------------------------

