
"""
GeneS - Generate Sources

Common tools and objects.
"""


# ------------------------------------------------
# logging
# ------------------------------------------------

# modules:
import sys
import logging

# setup logging if not done yet:
logging.basicConfig( stream=sys.stdout, level=logging.INFO,
                       format='[%(levelname)-8s] %(message)s' )

# get root logger instance:
logger = logging.getLogger()


# ------------------------------------------------
# end
# ------------------------------------------------
