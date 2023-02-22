"""
C3PO - CF Convention Compliance Python Objects

Tools module with general purpose routines.
  
"""


def HistoryLine( msg, t=None, tstamp=None ) :

    """
    Create history message line with time stamp:

      Mon 1 Jan 01:01:01 2001: message

    The time stamp label might be explicitly provided as argument.
    If not provided, it is created from datetime 't' if present,
    or otherwise from the current time.
    """

    # modules:
    import datetime

    # set time stamp if not specified:
    if tstamp == None :
        # fill time stamp, by default current time:
        if t == None : t = datetime.datetime.today()
        # convert to "Mon 1 Jan 01:01:01 2001"
        tstamp = t.strftime('%a %b %d %H:%M:%S %Y' )
    #endif

    # fill event:
    return '%s: %s' % (tstamp,msg)

#enddef  HistoryLine
