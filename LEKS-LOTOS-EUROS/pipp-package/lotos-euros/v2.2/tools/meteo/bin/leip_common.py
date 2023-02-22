#! /usr/bin/env python3

"""
Leip common stuff.
"""


########################################################################
###
### logger
###
########################################################################

# modules:
import sys
import logging

# setup logging if not done yet:
logging.basicConfig( stream=sys.stdout, level=logging.INFO,
                       format='[%(levelname)-8s] %(message)s' )

# get root logger instance:
logger = logging.getLogger()


########################################################################
###
### exceptions
###
########################################################################

class LeipError( Exception ) :
    """
    Raised after Leip routine encounters a problem.
    """
#endclass


########################################################################
###
### end
###
########################################################################
