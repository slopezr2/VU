
#-------------------------------------------------
# help
#-------------------------------------------------

"""
**********************
``utopya_base`` module
**********************

The ``utopya_base`` module provides access to UTOPyA's base class.

"""

#-------------------------------------------------
# base class
#-------------------------------------------------

class UtopyaBase( object ) :

    """
    Base class for UTOPyA objects.
    
    The main purpose is to define an attribute :py:attr:`logger` 
    of class :py:class:`logging.Logger`, 
    based on the idea that every object needs to issue messages sometimes, 
    either for information or to describe an error situation.
    Use this attribute to issue messages to the logging system::
    
        self.logger.info( 'Important message.' )
        self.logger.error( 'Important issue.' )
        
    The name of the defined logger is 'UTOPyA'.
    Applications using UTOPyA objects might use this name to configure the 
    logging messages according to the applications needs; no idea how actually ...
    
    The following methods are provided:
    
    * :py:meth:`LoggerGet` to return the 'logger' instance;
    * :py:meth:`LoggerSet` and :py:meth:`LoggerReset` to (temporary) configure the logging level and formatting;
    * :py:meth:`LoggerTest` to test the logging system.
    """
    
    def __init__( self ) :
    
        """
        Initialize UTOPyA base class.
        """
        
        # modules:
        import logging
        import sys
        
        # tools:
        import rc
        
        # defaults for root logger:
        logging.basicConfig()

        # define logger attribute.
        self.logger = logging.Logger('UTOPyA')
        
        # log only info and higher:
        self.logger.setLevel(logging.INFO)
        # standard format and level:
        logformat = '%(asctime)s [%(levelname)-8s] %(message)s'
        loglevel = logging.DEBUG
        # create standard output handler:
        self._stdout_handler = logging.StreamHandler(sys.stdout)
        # set level to debug or info:
        self._stdout_handler.setLevel(loglevel)
        # set format:
        self._stdout_handler.setFormatter(logging.Formatter(logformat))
        # add:
        self.logger.addHandler( self._stdout_handler )

        # maintain stack of logging properties:
        self._logstack = [ { 'format' : logformat, 'level' : loglevel } ]

    #enddef __init__
    
    # *
    
    def __del__( self ) :
    
        """
        Destructor for UTOPyA base class.
        """
        
        # cleanup logger:
        if hasattr(self,'logger') : del self.logger
        
    #enddef __del__
    
    # *
    
    def LoggerGet( self ) :
    
        """
        Return logger instance.
        """
        
        # ok
        return self.logger
        
    #enddef LoggerGet
    
    # *
    
    def LoggerTest( self ) :
    
        """
        Test logging system by printing messages for different logging levels.
        """

        # messages ...
        self.logger.info( 'test of logging:' )
        self.logger.debug   ( '  debug    message ...' )
        self.logger.info    ( '  info     message ...' )
        self.logger.warning ( '  warning  message ...' )
        self.logger.error   ( '  error    message ...' )
        self.logger.critical( '  critical message ...' )
        
    #enddef LoggerTest

    # *

    def LoggerSet( self, level='info', format=None, indent=None ) :

        """
        Temporary change logging configuratation to only issue messages 
        from the specified level onwards, with eventually a new format
        or extra indent.
        
        Each call to this routine should be followed later on
        by a call to :py:meth:`LoggerReset`.
        
        The 'level' argument is either an integer constants from the logging
        module (``logging.DEBUG``) or an equivalent character string.
        
        Example of format::
        
            '%(asctime)s [%(levelname)-8s] %(message)s'
        
        See also the details in the documentation of the :py:mod:`logging` module:
        
        * `Logging Levels <https://docs.python.org/2/library/logging.html#levels>`_  (keywords for levels)
        * `LogRecord attributes <https://docs.python.org/2/library/logging.html#logrecord-attributes>`_ (formatting templates)
        
        """

        # modules:
        import logging
        
        # copy previous format if no alternative is provided:
        if format is None : format = self._logstack[-1]['format']
        
        # add extra indent ?
        if indent is not None :
            # insert indent before message:
            format = format.replace( '%(message)', indent+'%(message)' )
        #endif

        # reset format:
        self._stdout_handler.setFormatter(logging.Formatter(format))
        
        # integer level number:
        if type(level) == int :
            loglevel = level
        elif type(level) == str :
            if level.lower() == 'info' :
                loglevel = logging.INFO
            elif level.lower() == 'debug' :
                loglevel = logging.DEBUG
            else :
                self.logger.error( 'unsupported level keyword "%s"' % level )
                raise Exception
            #endif
        else :
            self.logger.error( 'usupported type "%s" for level argument "%s"' % (type(level),level) )
            raise Exception
        #endif
        # reset level:
        self.logger.setLevel(loglevel)

        # add property to stack:
        self._logstack.append( { 'format' : format, 'level' : loglevel } )

    #enddef LoggerSet

    # *

    def LoggerReset( self ) :

        """
        Reset loging parameters to those used prior to the
        preceeding call to :py:meth:`LoggerSet` .
        """

        # modules:
        import logging

        # stack should not be empty:
        if len(self._logstack) == 0 :
            self.logger.error( 'ResetLogging called with empty stack; corresponding call to SetLogging is missing ?' )
            raise Exception
        #endif

        # remove top:
        self._logstack.pop()
        # extract top:
        top = self._logstack[-1]
        format = top['format']
        level  = top['level']
        # reset:
        self._stdout_handler.setFormatter(logging.Formatter(format))
        self._stdout_handler.setLevel(level)

    #enddef LoggerReset
    
    # *

    def Call( self, command, **kwargs ) :

        """
        Class to call a subprocess and trap exceptions.
        Output and error messages are written to logging system.
        Optional keyword arguments are passed to the :py:class:`Call` class
        of the :py:mod:`utopya_tools` module;
        an instance of this class is returned.
        """
        
        # tools:
        import utopya_tools

        # call supprocess, trap errors:
        p = utopya_tools.Call( command, logger=self.logger, **kwargs )
        
        # ok
        return p
        
    #enddef Call

#endclass UtopyaBase


#-------------------------------------------------
# end
#-------------------------------------------------


