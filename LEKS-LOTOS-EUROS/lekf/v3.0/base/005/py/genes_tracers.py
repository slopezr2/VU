#! /usr/bin/env python

# ------------------------------------------------
# help
# ------------------------------------------------

"""
GeneS - Generate Sources

GeneSTracers - class to store info on tracers

"""

# ------------------------------------------------
# common
# ------------------------------------------------

# get common stuff:
#   - logger instance
from genes_common import *


# ------------------------------------------------
# routines
# ------------------------------------------------

def split_formula( formula ) :

    """
    Split formula 'C 2 + H 6 + O' into dictionairy with 
    per atom the number of occurances.
    Atoms C, S, and N are always added.
    Atom name 'R' represents rest of molecule.
    """
    # init result:
    res = {}
    # split at '+ sign:
    parts = formula.split('+')
    # loop:
    for part in parts :
        # split in factors:
        factors = part.split()
        # check ...
        if len(factors) == 1 :
            # single atom
            number = 1
            atom = factors[0]
        elif len(factors) == 2 :
            # extract atom and number:
            if factors[0].isalpha() and factors[1].isdigit() :
                number = int(factors[1])
                atom = factors[0]
            elif factors[1].isalpha() and factors[0].isdigit() :
                number = int(factors[0])
                atom = factors[1]
            else :
                logger.error( 'could not split factors "%s" in number and atoms' % part )
                raise Exception
            #endif
        else :
            logger.error( 'too many factors in "%s" of formula "%s"' % (part,formula) )
            raise Exception
        #endif
        # store:
        if atom not in res.keys() : res[atom] = 0
        res[atom] = res[atom] + number
    #endfor
    # add some defaults:
    for atom in ['C','N','S'] :
        if atom not in res.keys() : res[atom] = 0
    #endfor

    # ok
    return res

#enddef

# ------------------------------------------------
# class
# ------------------------------------------------

class GeneSTracers( object ) :

    """
    Class to store info on tracers.
    """
    
    def __init__( self, txtfile, groups, props, skip_tracers ) :

        """
        Read species from table.

        Example of input file:

          # Species names for the CBM4 core mechanism.
          name    , description                      , units,  formula, aerdens,  aerhygro, groups, properties
          #                                                             kg/m3       1
          NO      , Nitric oxide                     , ppb  , N + O   ,        ,          , cbm4  ,
          NO2     , Nitrogen dioxide                 , ppb  , N + O 2 ,        ,          , cbm4  ,
          SO4a_c  , Sulfate aerosol coarse mode      , ug/m3, S + O 4 , 1.77e3 , 0.8      , sia   , aerosol coarse_mode
          
        Provide list of groups to which a tracer belongs to enable certain tracers.
        Provide list of selected properties to create groups of tracers which
        share some features.
        
        Arguments:
          txtfile      : name of input file
          groups       : dictionairy with fields:
                           'supported' : list with supported group names
                           'selected'  : list with selected group names
          props        : dictionairy with fields:
                           'supported' : list with supported property names
          skip_tracers : list with tracer names (lower case) that should be disabled
                            regardless of the selected groups.
        
        Result is an object with fields:
          nspec       : total number of tracers
          name        : list with the names
          desc        : list with the descriptions
          units       : list with the tracer units
          aerdens     : list with aerosol densities
          aerhygro    : list with aerosol hygroscopicities
          formula     : list with formula dictionairies, e.g.:
                            { 'N' : 1, 'O' : 2, 'C' : 0, 'S' : 0 }
                        atoms 'C', 'N', and 'S' are always present
          groups      : list with group lists
          props       : list with property lists
          enabled     : list with booleans
        
        """

        # external:
        import os

        # file parameters:
        fcomment = '#'

        # init output:
        self.name     = []
        self.desc     = []
        self.units    = []
        self.aerdens  = []
        self.aerhygro = []
        self.formula  = []
        self.groups   = []
        self.props    = []
        self.enabled  = []

        # found ?
        if not os.path.exists(txtfile) :
            logger.error( 'input file not found:' )
            logger.error( '  %s' % txtfile )
            raise Exception
        #endif

        # open file:
        f = open( txtfile, mode='r' )
        
        # no header yet:
        headers = None
        # line counter:
        iline = 0

        # read lines until eof
        while True:

            # read line:
            line = f.readline()
            iline = iline + 1

            # zero length ? then eof (empty line has always a cr/lf):
            if len(line) == 0 : break
            
            # cleanup:
            line = line.strip()

            # skip comment:
            if line.startswith(fcomment) : continue
            # skip empty lines:
            if len(line) == 0 : continue

            # debug ...
            #logger.debug( line )

            # split at comma's :
            fields = line.split(',')
            
            # header line ?
            if headers == None :
                headers = fields
                nfield = len(headers)
                continue
            #endif
            
            # check ...
            if len(fields) != nfield :
                logger.error( 'number of fields %i does not match with number of headers %i' % (len(fields),nfield) )
                logger.error( '(line %i of file %s)' % (iline,txtfile) )
                raise Exception
            #endif
            
            # loop over fields:
            for ifield in range(nfield) :
                header = headers[ifield].strip()
                field = fields[ifield].strip()
                if header == 'name' :
                    tracer_name = field
                    self.name.append( tracer_name )
                elif header == 'description' :
                    self.desc.append( field )
                elif header == 'formula' :
                    formula = split_formula(field)
                    self.formula.append( formula )
                elif header == 'units' :
                    self.units.append( field )
                elif header == 'aerdens' :
                    # replace missing values by dummy:
                    if len(field) == 0 : field = '-999.9e3'
                    # store:
                    self.aerdens.append( field )
                elif header == 'aerhygro' :
                    # replace missing values by dummy:
                    if len(field) == 0 : field = '-999.9'
                    # store:
                    self.aerhygro.append( field )
                elif header == 'groups' :
                    self.groups.append( field.split() )
                elif header == 'properties' :
                    self.props.append( field.split() )
                else :
                    logger.error( 'unsupported header "%s" in file "%s"' % (header,txtfile) )
                    raise Exception
                #endif
            #endif

            # check if all groups are supported ...
            for group in self.groups[-1] :
                # supported ?
                if group not in groups['supported'] :
                    logger.error( 'tracer grouperty "%s" not one of supported : %s' % (group,str(groups['supported'])) )
                    raise Exception
                #endif
            #endfor

            # check if all properties are supported ...
            for prop in self.props[-1] :
                # supported ?
                if prop not in props['supported'] :
                    logger.error( 'tracer property "%s" not one of supported : %s' % (prop,str(props['supported'])) )
                    raise Exception
                #endif
            #endfor
            
            # * enabled ?

            # not enabled by default ...
            enabled = False

            # loop over groups set for current tracer;
            # if group is selected, then enable this tracer:
            for group in self.groups[-1] :
                # match ?
                if group in groups['selected'] :
                    enabled = True
                    break
                #endif
            #endfor
            
            # explicitly disabled ?
            if tracer_name.lower() in skip_tracers :
                enabled = False
            #endif
            
            # store:
            self.enabled.append( enabled )
            
            # *

        #endwhile

        # post process:
        self.nspec = len(self.name)

        # close:
        f.close()

    #enddef  # __init__
    
        
#endclass   # GeneSTracers


# ------------------------------------------------
# end
# ------------------------------------------------

