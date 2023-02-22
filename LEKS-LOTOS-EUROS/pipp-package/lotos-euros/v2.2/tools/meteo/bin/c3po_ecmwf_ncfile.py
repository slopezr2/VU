#! /usr/bin/env python3

#######################################################################
###
### main
###
#######################################################################

"""
C3PO - CF Convention Compliance Python Objects

Tools to convert NetCDF files produced by standard software available
at ECMWF, i.e. mars2netcdf or cdo, towards full CF-Convention compliance.
The standard tools do not always provide complied results.


**************************************************
DEPRICATED: convert directly from GRIB in future !
**************************************************


Usage:

    # open existing file for read/write:
    ecfile = EcGridFile( 'myfile.nc' )

    # change variable names and attributes following parameter tables:
    ecfile.ComplyVariables()
    
    # common changes to time variable:
    ecfile.ComplyTimeVariable()
    
    # improve hybride level desciption, 
    # add the surface pressure fields obtained from a log-normal-surface pressure file:
    ecfile.ComplyModelLevels( lnspfile )
    
    # gather non-zero lon/lat values into a single coordinate:
    ecfile.Gather()
    
    # convert from accumulated to temporal avarge fields:
    ecfile.DeAccum( 'time=00/12, step=03/06/09/12' )

    # close:
    ecfile.close()

HISTORY

    2013-01, Arjo Segers

"""


#######################################################################
###
### ECMWF gridded file
###
#######################################################################



# modules:
import c3po_file

# define clas:
class EcGridFile( c3po_file.GridFile ) :
    
    """
    Access NetCDF file produced at ECMWF.
    """
    
    def ComplyVariables( self, indent='' ) :
        
        """
        Rename variables if necessary, add appropriate attribues.
        
        Some variables are only identified by their parameter number;
        the corresponding table number is in the attributes:

          float var152(time, lev, lat, lon) ;
                  var152:table = 128 ;
        
        The name corresponding to parameter (152) in a talbe (128) is
        read from the tables included in the GRIB_API library;
        the variable will be renamed to this.
        From the GRIB_API tables also attributes like units, long names,
        and CF standard names are sometimes available;
        if present, these are added as attributes.

        Optionally specify indent for logging messages.
        """
        
        # modules:
        import logging
        
        # tools:
        import ecparamtable
        
        # no parameter table read yet:
        partab = None
        
        # flag:
        changed = False
        
        # info ...
        logging.debug( indent+'comply variables ...' )
        # loop over variables:
        for varname in self.ncid.variables.keys() :
            # info ...
            logging.debug( indent+'  found variable "%s" ...' % varname )
            # remove extensions '_2' etc for variable originating from same grib code
            # (but in different table ...):
            if '_' in varname :
                varname = varname.split('_')[0]
                logging.debug( indent+'    reduced name to "%s" ...' % varname )
            #endif
            # is this one of the variables we are looking for?
            if varname.startswith('var') and varname[3:].isalnum() :
                # info ...
                logging.debug( indent+'    try to rename ...' )
                # extact parameter number:
                parid = int( varname[3:] )
                # current variable:
                varid = self.ncid.variables[varname]
                # should have an attribute with the table number:
                attrname = 'table'
                if attrname not in varid.ncattrs() :
                    # info ...
                    logging.debug( indent+'    no "%s" attribute found, skip ...' % attrname )
                    # next:
                    continue
                #endif
                # extract table number:
                tabid = int( varid.getncattr(attrname) )
                # no parameter table yet ? then need to read:
                if partab == None :
                    # info ...
                    logging.debug( indent+'    read parameter tables ...' )
                    # read tables, use environment for search path:
                    partab = ecparamtable.GRIB_API_Tables( keys=['paramId','shortName','name','units','cfName'] )
                #endif
                # table should be available:
                if tabid not in partab.keys() :
                    # info ...
                    logging.debug( indent+'    no local table "%i" available.' % tabid )
                    # next:
                    continue
                #endif
                # local table:
                ltable = partab[tabid]
                # parameter should be available:
                if parid not in ltable.keys() :
                    # info ...
                    logging.debug( indent+'    no parameter "%i" found in local table "%i".' % (parid,tabid) )
                    # next:
                    continue
                #endif
                # extract parameter info:
                pari = ltable[parid]
                # rename variable to shortname if possible:
                key = 'shortName'
                if key in pari.keys() :
                    # new name:
                    newname = pari[key]
                    # adhoc fix, cf-checker does not like variables starting with a digit:
                    if newname == '2t' : newname = 'T2M'
                    if newname == '2d' : newname = 'D2M'
                    # switch:
                    if newname == '~' :
                        # info ...
                        logging.debug( indent+'    do not rename to "%s", experimental product ...' % newname )
                    else :
                        # info ...
                        logging.debug( indent+'    rename to "%s" ...' % newname )
                        # rename:
                        self.ncid.renameVariable( varname, newname )
                        # reset flag:
                        changed = True
                    #endif
                else :
                    # info ...
                    logging.debug( indent+'    no "%s" found in parameter table, keep current name' % key )
                #endif
                # map parameter info to attributes:
                amap = {}
                amap['units'  ] = 'units'
                amap['name'   ] = 'long_name'
                amap['cfName' ] = 'standard_name'
                # loop over parameter info keys:
                for key in amap.keys() :
                    # defined ?
                    if key in pari.keys() :
                        # add attribute if not defined yet:
                        ch = self.DefaultAttribute( varid, amap[key], pari[key], indent=indent+'    ' )
                        changed = changed or ch
                    #endif
                #endfor  # attribute mapping
                # add attribute with parameter id:
                ch = self.DefaultAttribute( varid, 'param', parid, indent=indent+'    ' )
                changed = changed or ch
            #endif  # variable 'varNNN'
        #endfor  # variables in file
        
        # add line to history:
        if changed : self.ExtendHistory( 'c3po: renamed variables and/or added attributes using GRIB_API tables' )

    #enddef  # ComplyVariables
    
    
    # ***
    
    
    def ComplyTimeVariable( self, indent='' ) :
        
        """
        Comply the 'time' variable to some conventions:
        - add 'standard_name' attribute with value 'time';
        - add 'axis' attribute with value 'T' .
        
        Optionally specify indent for logging messages.
        """
        
        # modules:
        import logging
        
        # loop over variables:
        for varname in self.ncid.variables.keys() :
            # time variable ?
            if varname == 'time' :
                # info ...
                logging.debug( indent+'comply time variable ...' )
                # current variable:
                varid = self.ncid.variables[varname]
                # add attributes if not present yet:
                self.DefaultAttribute( varid, 'standard_name', 'time', indent=indent+'  ' )
                self.DefaultAttribute( varid, 'axis'         , 'T'   , indent=indent+'  ' )
                # add line to history:
                self.ExtendHistory( 'c3po: added attributes to time variable' )
            #endif  # time variable
        #endfor  # variables in file
        

    #enddef  # ComplyTimeVariable


    # ***
    
    def ComplyModelLevels( self, lnspfile, indent='', cleanup=True ) :

        """
        Change model level definition following CF-conventions for the 
        vertical coordinate.

        The NetCDF files converted from GRIB using 'cdo -f nc copy'
        is not completely conform these conventions because:
         - hybride coefficients 'hy[ab][mi]' are defined with dimensions
           named 'nhy[mi]'
         - no surface pressue field is present for the hybride coordinate

        Dimensions and/or variables cannot be removed from a file.
        Therefore, the open file is first closed and a temporary copy is created
        with extension '_ml'. The required variables/attributes are then copied
        from there into a newly created version of the current file.
        Set optional argument 'cleanup' to 'False' to keep the temporary copy.

        The new file that is created is basically a copy of the original
        model level file, with the following changes:
         - the hybride coeff at the mid of the layer are renamed to
           'hy[ab]' and have dimension 'lev' ;
         - boundary variables 'lev_bnds' and 'hy[ab]_bnds' of shape (nlev,2) 
           are defined for the hybride coeff at the layer interfaces;
         - a surface pressure field is extraced from a file with
           the log-normal-surface pressure, this is usually the quantity
           that is available for model level fields.
        """

        # modules:
        import os
        import logging
        import numpy
        import netCDF4
        import datetime

        # info ...
        logging.debug( indent+'comply model level definition ...' )
        
        # name of temporary copy:
        mlfile = self.filename.replace('.nc','_ml.nc')
        # create copy, return its dimension id;
        # create new empty file with original name, copy global attributes:
        mlf = self.CreateTemporaryCopy( mlfile )

        # check ...
        if not os.path.isfile(lnspfile) :
            logging.error( 'log-normal-surface-pressure l file not found : %s' % lnspfile )
            raise Exception
        #endif
        # open:
        lnspf = netCDF4.Dataset( lnspfile, 'r' )

        # extend history:
        self.ExtendHistory( 'c3po: added surface pressure from %s, re-formatted level definition' % os.path.basename(lnspfile) )

        # loop over dimensions:
        for dimname in mlf.dimensions.keys() :
            # skip, change, or copy ?
            if dimname in ['mlev','nhym'] :
                # used for hybride coeff over all original layers;
                # will be replaced by the subset 'lev' :
                continue
            elif dimname in ['ilev','nhyi'] :
                # used for hybride coeff at original layer interfaaces;
                # will be replaced by boundary varaibles of shape (nlev,2),
                # so define a dimension for the number of vertices:
                self.createDimension( 'bnds', 2 )
            else :
                # size:
                if mlf.dimensions[dimname].isunlimited() :
                    n = 0
                else :
                    n = len(mlf.dimensions[dimname])
                #endif
                # create:
                self.createDimension( dimname, n )
            #endif
        #endfor

        # level indices, 1-based:
        levs = mlf.variables['lev'][:]
        # count:
        nlev = len(levs)

        # init values, might become defined:
        varid_lev = None
        with_lev_bnds = False
        
        # init list of data variables,
        # need to copy these later on after strange errors ...
        xvarnames = []

        # loop over variables:
        for varname in mlf.variables.keys() :
            # variable object:
            varid = mlf.variables.get(varname)
            #
            # special ?
            #
            #~ level indices:
            if varname == 'lev' :
                # info ...
                logging.debug( indent+'  create new version of variable "%s" ...' % varname )
                # create new variable;
                # originally contains pressures using standard surface pressure,
                # but this is not in agreement with the standard units '1' for
                # variable 'atmosphere_hybrid_sigma_pressure_coordinate' :
                newvarid = self.ncid.createVariable( varname, 'i4', varid.dimensions )
                # copy attributes:
                for attrname in varid.ncattrs() :
                    # current value:
                    value = varid.getncattr(attrname)
                    # add to new file:
                    newvarid.setncattr( attrname, value )
                #endfor
                # add units, probably replaced the depricated value 'level' :
                newvarid.setncattr( 'units', '1' )
                # replace formula attributes:
                newvarid.setncattr( 'formula'      , 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)' )
                newvarid.setncattr( 'formula_terms', 'ap: hya b: hyb ps: ps' )
                # replace standard name:
                newvarid.setncattr( 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate' )
                # set values as numbers 1,..,nlev
                newvarid[:] = range(1,nlev+1)
                # store for later usage:
                varid_lev = newvarid
                #
            #~ hybride coefficients for full levels:
            elif varname in ['hyam','hybm'] :
                # info ...
                logging.debug( indent+'  create new version of variable "%s" ...' % varname )
                # create variable:
                newvarid = self.ncid.createVariable( varname[0:3], varid.dtype, ('lev',) )
                # copy attributes:
                for attrname in varid.ncattrs() :
                    # current value:
                    value = varid.getncattr(attrname)
                    # add to new file:
                    newvarid.setncattr( attrname, value )
                #endfor
                ## add attributes for boundaries, named 'hya_bnds' or 'hyb_bnds' :
                #newvarid.setncattr( 'bounds', varname[0:3]+'_bnds' )
                # extract subset from all layers:
                values = []
                for ilev in range(nlev) :
                    k = int(levs[ilev])
                    newvarid[ilev] = varid[k-1]
                #endfor
                #
            #~ hybride coefficients for half levels:
            elif varname in ['hyai','hybi'] :
                # info ...
                logging.debug( indent+'  create new version of variable "%s" ...' % varname )
                # create new boundary variables 'hya_bnds' or 'hyb_bnds' :
                newvarid = self.ncid.createVariable( varname[0:3]+'_bnds', 
                                                  varid.dtype, ('lev','bnds') )
                # copy attributes:
                for attrname in varid.ncattrs() :
                    # current value:
                    value = varid.getncattr(attrname)
                    # add to new file:
                    newvarid.setncattr( attrname, value )
                #endfor
                # extract subset from all layers:
                for ilev in range(nlev) :
                    k = int(levs[ilev])  # 1..L
                    newvarid[ilev,0] = varid[k-1]
                    newvarid[ilev,1] = varid[k]
                #endfor
                # reset flag:
                with_lev_bnds = True
                #
            #~ standard variable:
            else :
                # store name to copy later on:
                xvarnames.append( varname )
            #endif
            
        #endfor  # variables:

        # add level boundaries ?
        if with_lev_bnds :
            # new variable:
            varname = 'lev_bnds'
            # info ...
            logging.debug( indent+'  create new variable "%s" ...' % varname )
            ## adhoc fix ..
            #self.ncid.close()
            #self.ncid = netCDF4.Dataset( self.filename, mode='r+' )
            # add level boundaries variable:
            newvarid = self.ncid.createVariable( varname, 'i4', ('lev','bnds') )
            # add attributes:
            newvarid.setncattr( 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate' )
            newvarid.setncattr( 'units', '1' )
            # according to the CF-compliance-checker this is not a coordinate variable,
            # and therefore the 'formula_terms' attribute is not allowed ;
            # however, it is useful to have it, so here a different name is used:
            newvarid.setncattr( '_formula'      , 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)' )
            newvarid.setncattr( '_formula_terms', 'ap: hya_bnds b: hyb_bnds ps: ps' )
            # fill dummy values:
            values = numpy.zeros((nlev,2),int)
            values[:,0] = numpy.arange(0,nlev  )
            values[:,1] = numpy.arange(1,nlev+1)
            # write:
            newvarid[:] = values
            # add reference to bounds variable::
            varid_lev.setncattr( 'bounds', varname )
        #endif  # with level boundaries
        
        # add surface pressure field:
        if True :
            # new variable:
            varname = 'ps'
            # info ...
            logging.debug( indent+'  create new variable "%s" ...' % varname )
            # get lnsp-variable object:
            lnspvname = None
            xnames = ['lnsp','LNSP']
            for xname in xnames :
                if xname in lnspf.variables.keys() :
                    lnspvname = xname
                    break
                #endif
            #endfor
            if lnspvname == None :
                logging.error( 'no lnsp variable found under one of these names : %s' % xnames )
                raise Exception
            #endif
            varid = lnspf.variables[lnspvname]
            # extract dimension names:
            dims = varid.dimensions
            # the 'lev' dimesion is 1 here, remove if present;
            # since dims is a tupple there is no remove procedure ...
            newdims = ()
            for dim in dims :
                if dim != 'lev' : newdims = newdims+(dim,)
            #endfor
            # create variable:
            newvarid = self.ncid.createVariable( varname, varid.dtype, newdims )
            # copy attributes:
            for attrname in varid.ncattrs() :
                # extract value:
                avalue = varid.getncattr(attrname)
                # store:
                newvarid.setncattr( attrname, avalue )
            #endfor
            # (re)new attributes:
            newvarid.setncattr( 'units'        , 'Pa'                   )
            newvarid.setncattr( 'param'        , 134                    )
            newvarid.setncattr( 'long_name'    , 'surface air pressure' )
            newvarid.setncattr( 'standard_name', 'surface_air_pressure' )
            # copy values, convert from log-normal to surface-pressure: 
            newvarid[:] = numpy.exp( varid[:] )
        #endif  # ps variable

        # loop over data variables:
        for varname in xvarnames :
            # variable object:
            varid = mlf.variables.get(varname)
            # info ...
            logging.debug( indent+'  copy "%s" ...' % varname )
            # create variable:
            newvarid = self.ncid.createVariable( varname, varid.dtype, varid.dimensions )
            # copy attributes:
            for attrname in varid.ncattrs() :
                # current value:
                value = varid.getncattr(attrname)
                # add to new file:
                newvarid.setncattr( attrname, value )
            #endfor
            # copy values; use small chunks to save memory:
            rank = len(varid.shape)
            if rank == 1 :
                newvarid[:] = varid[:]
            else :
                nrec = varid.shape[0]
                for irec in range(nrec) :
                    newvarid[irec,:] = varid[irec,:]
                #endfor
            #endif
        #endfor  # variables

        # close temporary copy:
        mlf.close()
        # close surface pressure input file:
        lnspf.close()

        # remove temporary file ?
        if cleanup : os.remove( mlfile )

    #enddef   # ComplyModelLevels


    # ***

    
    def ComplyModelLevels2( self, sp_ncfile, lnsp_ncfile_cdo, indent='', cleanup=True ) :

        """
        Change model level definition following CF-conventions for the 
        vertical coordinate for file created by 'grib_to_netcdf'.
        
        The files created by 'grib_to_netcdf' lack the following:
         - no surface pressure field yet;
         - attributes of 'lev' are incompleet;
         - no 'levi' coordinate for level interfaces;
         - no hybride coordinate variables.
         
         As input, two file names should be supplied:
         - A file 'sp_ncfile' for a netcdf file with a surface pressure field,
           created with 'grib_to_netcdf' on a model level lnsp file
           and conversion from lnsp to sp using ncap2.
         - A file 'lnsp_ncfile_cdo' which is the result of a conversion
           from grib to netcdf using cdo, for example using:
              cdo -f nc copy lnsp.gb lnsp.nc
           This conversion will and suitable hybride coefficents 'hy[ab][mi]'.  
        
         The following changes are applied:
         - copy 'sp' from 'sp_ncfile':
           - convert from 'lnsp'
           - keep packing
         - fix 'lev' coordinate:
           - add attributes
           - copy 'hy[ab]m' from 'lnsp_ncfile_cdo'
         - create 'levi' coordinate:
           - create 'levi' dimension
           - create 'levi' coordinate
           - copy 'hy[ab]i' from 'lnsp_ncfile_cdo'

        """

        # modules:
        import os
        import logging
        import numpy
        import netCDF4

        # info ...
        logging.debug( indent+'comply model level definition ...' )
        
        # check ...
        if not os.path.isfile(sp_ncfile) :
            logging.error( 'surface-pressure file not found : %s' % sp_ncfile )
            raise Exception
        #endif
        # open:
        spfile = netCDF4.Dataset( sp_ncfile, 'r' )
        # variable:
        sp_vname = 'sp'
        # check ...
        if sp_vname not in spfile.variables.keys() :
            logging.error( 'no variable "%s" found in "%s"' % (sp_vname,sp_ncfile) )
            raise Exception
        #endif
        # create copy in current file, no data yet:
        varid = self.createVariableCopy( sp_vname, spfile.variables[sp_vname], data=True )
        # close:
        spfile.close()
        
        # extend history:
        self.ExtendHistory( 'c3po: added surface pressure from %s' % os.path.basename(sp_ncfile) )
        
        # *

        # check ...
        if not os.path.isfile(lnsp_ncfile_cdo) :
            logging.error( 'surface-pressure file not found : %s' % lnsp_ncfile_cdo )
            raise Exception
        #endif
        # open:
        spfile = netCDF4.Dataset( lnsp_ncfile_cdo, 'r' )
        
        # name of level coordinate in input file:
        vname_in = 'lev'
        # check ...
        if vname_in not in spfile.variables.keys() :
            logging.error( 'no variable "%s" found in "%s"' % (vname_in,lnsp_ncfile_cdo) )
            raise Exception
        #endif
        # input variable:
        varid_in = spfile.variables[vname_in]
        
        # name of level coordinate in output file:
        vname = 'level'
        # check ...
        if vname not in self.ncid.variables.keys() :
            logging.error( 'no variable "%s" found in "%s"' % (vname,self.filename) )
            raise Exception
        #endif
        # target variable:
        varid = self.ncid.variables[vname]
        # loop over input attributes:
        for aname in varid_in.ncattrs() :
            # original value:
            value = varid_in.getncattr( aname )
            # reset some:
            if aname == 'standard_name' : value = 'atmosphere_hybrid_sigma_pressure_coordinate'
            if aname == 'units'         : value = '1'
            if aname == 'formula'       : value = 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)'
            if aname == 'formula_terms' : value = 'ap: hyam b: hybm ps: %s' % sp_vname
            # copy:
            varid.setncattr( aname, value )
        #endfor # attributes
        # size:
        nlevel = len(self.ncid.dimensions[vname])
        # values:
        levels = varid[:]
        # copy coefficients:
        for cname in ['hyam','hybm'] :
            # copy including data, use new dimension name:
            self.createVariableCopy( cname, spfile.variables[cname], 
                                       dimensions=(vname),
                                       data=True, selection=levels-1 )
        #endfor
        
        # create dimension for level intervaces:
        hvname = 'hlevel'
        # size:
        nhlevel = nlevel + 1
        # create values:
        hlevels = numpy.zeros( (nhlevel), dtype=varid.dtype )
        if levels[0] < levels[-1] :
            hlevels[0] = levels[0] - 1
            hlevels[1:] = levels
        else :
            hlevels[0:nlevel-1] = levels
            hlevels[nlevel] = levels[-1] - 1
        #endif
        # create new dimension:
        self.ncid.createDimension( hvname, nhlevel )
        # create coordinate variable:
        varidi = self.ncid.createVariable( hvname, varid.dtype, (hvname) )
        # fill data:
        varidi[:] = hlevels
        # loop over attributes of full level coordinate:
        for aname in varid.ncattrs() :
            # original value:
            value = varid.getncattr( aname )
            # reset some:
            if aname == 'long_name'     : value = 'hybrid level at layer interfaces'
            if aname == 'formula_terms' : value = 'ap: hyai b: hybi ps: %s' % sp_vname
            # copy:
            varidi.setncattr( aname, value )
        #endfor # attributes
        # copy coefficients:
        for cname in ['hyai','hybi'] :
            # copy including data, use new dimension name:
            self.createVariableCopy( cname, spfile.variables[cname], 
                                       dimensions=(hvname),
                                       data=True, selection=hlevels )
        #endfor

        # close:
        spfile.close()
        
        # extend history:
        self.ExtendHistory( 'c3po: added level coordinates from %s' % os.path.basename(lnsp_ncfile_cdo) )

    #enddef   # ComplyModelLevels2
    
    
    # ***


    def Gather( self, indent='', cleanup=True ) :

        """
        Gather non zero grid cells.
        Now only for 2D arrays.

        Dimensions and/or variables cannot be removed from a file.
        Therefore, the open file is first closed and a temporary copy is created
        with extension '_ml'. The required variables/attributes are then copied
        from there into a newly created version of the current file.
        Set optional argument 'cleanup' to 'False' to keep the temporary copy.
        """

        # modules:
        import os
        import logging
        import numpy

        # info ...
        logging.debug( indent+'gather non-zero grid cells ...' )
        
        # name of temporary copy:
        gridfile = self.filename.replace('.nc','_gridded.nc')
        # create copy, return its dimension id;
        # create new empty file with original name, copy global attributes:
        gridf = self.CreateTemporaryCopy( gridfile )

        # dimension names:
        dimnames = gridf.dimensions.keys()
        # collect variable names:
        varnames = []
        for vname in gridf.variables.keys() :
            # special ?
            if vname in dimnames : continue
            if vname.endswith('_bounds') : continue
            # add:
            varnames.append( vname )
        #endfor

        # gathering dimension names:
        if 'lon' in gridf.dimensions.keys() :
            londimname = 'lon'
            latdimname = 'lat'
        else :
            londimname = 'longitude'
            latdimname = 'latitude'
        #endif
        # shape:
        nlon = len(gridf.dimensions[londimname])
        nlat = len(gridf.dimensions[latdimname])

        # flag to check if a gridcell has any data:
        hasdata = numpy.zeros((nlat,nlon),bool)
        # update for all variables:
        for varname in varnames :
            # variable:
            varid = gridf.variables[varname]
            # check, expecting (time,lat,lon) :
            if varid.ndim != 3 :
                logging.error( 'could only gather 3D fields, found "%s" with shape %s' % (varname,varid.shape) )
                raise Exception
            #endif
            # check ...
            if varid.dimensions[-2:] != (latdimname,londimname) :
                logging.error( 'expected last dimensions to be "%s", found "%s"' % ((latdimname,londimname),varid.dimensions[-2:]) )
                raise Exception
            #endif
            # read field:
            pat = varid[:]
            # where data ?
            hdata = pat > 0.0
            # dimensions:
            nt,ny,nx = pat.shape
            # uptate flags:
            for it in range(nt) : hasdata = hasdata | hdata[it,:,:]
        #endfor

        # indices with data:
        jj,ii = numpy.where( hasdata )

        # point dimension:
        pointname = 'point'
        npoint = len(ii)
        # zero based 1D indices in 0,..,nlat*nlon-1
        point = jj * nlon + ii

        # copy dimensions:
        for dimname in gridf.dimensions.keys() :
            # original:
            dimid = gridf.dimensions[dimname]
            # length:
            dimlen = len(dimid)
            if dimid.isunlimited() : dimlen = 0
            # create:
            self.createDimension( dimname, dimlen )
        #endfor
        # extra:
        self.createDimension( pointname, npoint )
        # create point index variable:
        varid = self.ncid.createVariable( pointname, 'i4', (pointname,) )
        varid.setncattr( 'long_name', 'zero-based indices in 1D array' )
        varid.setncattr( 'compress', '%s %s' % (latdimname,londimname) )
        varid.setncattr( 'description', 'original zero-based indices: %s_index = point/len(%s), %s_index = point mod len(%s)' % (latdimname,londimname,londimname,londimname) )
        varid[:] = point
        # loop over variables:
        for varname in gridf.variables.keys() :
            # original:
            varid = gridf.variables[varname]
            # data variable ?
            if varname in varnames :
                # new dimansions:
                dims = (varid.dimensions[0],pointname)
                # create copy:
                gvarid = self.ncid.createVariable( varname, varid.dtype, dims )
                # attributes:
                for aname in varid.ncattrs() : gvarid.setncattr( aname, varid.getncattr(aname) )
                # extract data:
                gg = varid[:]
                # dimensions:
                nt,ny,nx = gg.shape
                # fill:
                for it in range(nt) : gvarid[it,:] = gg[it,jj,ii]
            else :
                # create copy:
                gvarid = self.ncid.createVariable( varname, varid.dtype, varid.dimensions )
                # attributes:
                for aname in varid.ncattrs() : gvarid.setncattr( aname, varid.getncattr(aname) )
                # fill:
                gvarid[:] = varid[:]
            #endif
        #endfor

        # extend history:
        self.ExtendHistory( 'c3po: gathered non-zero cells' )

        # close temporary copy:
        gridf.close()

        # remove temporary file ?
        if cleanup : os.remove( gridfile )

    #enddef   # Gather


    # ***
    
    
    def DeAccum( self, time_steps, indent='' ) :

        """
        Convert from accumulated to temporal average fields.
        Provide the time/step values as passed to the mars job:
          [time = ]00/12, step=03/06/09/12
        """

        # modules:
        import logging
        import netCDF4
        import numpy
        import datetime
        
        ## testing ...
        #import shutil
        ## create backup:
        #shutil.copy( self.filename, self.filename+'.bak' )

        # info ...
        logging.debug( indent+'convert from accumulated values to time averages ...' )

        # split time/steps keywords:
        timelist,stepdef = time_steps.replace(' ','').split(',')
        # extract times:
        times = map(int,timelist.split('/'))
        # extract step values:
        if not stepdef.startswith('step=') :
            logger.error( 'second field of "%s" should be of the form "step = 03/06/09/12"' % time_steps )
            raise Exception
        #endif
        steplist = stepdef.split('=')[1].lower()
        if '/to/' in steplist :
            fields = steplist.split('/')
            if len(fields) == 5 :   #  3/to/120/by/3
                if (fields[1] != 'to') or (fields[3] != 'by') :
                    logger.error( 'unsupported step list "%s"' % steplist )
                #endif
                step_from = int(fields[0])
                step_to   = int(fields[2])
                step_by   = int(fields[4])
                steps = range(step_from,step_to+1,step_by)
            else :
                logger.error( 'unexpected number of fields in step list "%s"' % steplist )
                raise Exception
            #endif
        else :
            # list with values:
            steps = map(int,steplist.split('/'))
        #endif
        # info ...
        logging.debug( indent+'  time settings   : %s' % time_steps )
        logging.debug( indent+'  extracted steps : %s' % steps )

        # dimension names:
        dimnames = self.ncid.dimensions.keys()
        # collect variable names:
        varnames = []
        for vname in self.ncid.variables.keys() :
            # special ?
            if vname in dimnames : continue
            if vname.endswith('_bounds') : continue
            if vname in ['longitude','latitude'] : continue
            if vname in ['ulat_nlon','ulat_i0'] : continue
            if vname.endswith('_bnds') : continue
            # add:
            varnames.append( vname )
        #endfor
        
        # extract time values:
        tname = 'time'
        varid_time = self.ncid.variables[tname]
        time_units = varid_time.getncattr('units')
        time_vals = netCDF4.num2date( self.ncid.variables[tname][:], time_units )

        # loop over variables:
        for varname in varnames :
            # info ...
            logging.debug( indent+'  convert variable "%s" ...' % varname )
            # variable id:
            varid = self.ncid.variables[varname]
            # extract data:
            fields = varid[:]
            # info ...
            logging.debug( indent+'    value range : (%e,%e)' % (fields.min(),fields.max()) )
            # guess sign:
            vmin,vmax = None,None
            if fields.min() >= 0.0 : vmin = 0.0
            if fields.max() <= 0.0 : vmax = 0.0
            # init new time bounds data:
            time_bounds = []
            # dims:
            ndim = len(fields.shape)
            # loop over time records:
            irec = 0
            for time in times :
                ## info ...
                #logging.debug( indent+'    time = %s' % time )
                # init previous field, start with zero:
                step_prev = 0
                if ndim == 2 :
                    field_prev = fields[0,    :].copy() * 0.0
                elif ndim == 3 :
                    field_prev = fields[0,  :,:].copy() * 0.0
                elif ndim == 4 :
                    field_prev = fields[0,:,:,:].copy() * 0.0
                else :
                    logging.error( 'could not init field for ndim=%i' % ndim )
                    raise Exception
                #endif
                # loop over steps:
                for step in steps :
                    ## info ...
                    #logging.debug( indent+'      step = %s' % step )
                    # length of averaging interval, convert from hours to seconds:
                    dtsec = ( step - step_prev ) * 3600
                    # add time bounds:
                    time_bounds.append( [time_vals[irec]-datetime.timedelta(0,dtsec),time_vals[irec]] )
                    # current (accumulated) field:
                    if ndim == 2 :
                        field = fields[irec,    :].copy()
                    elif ndim == 3 :
                        field = fields[irec,  :,:].copy()
                    elif ndim == 4 :
                        field = fields[irec,:,:,:].copy()
                    else :
                        logging.error( 'could not copy field for ndim=%i' % ndim )
                        raise Exception
                    #endif
                    # replace record:
                    if ndim == 2 :
                        fields[irec,    :] = ( field - field_prev )/dtsec
                    elif ndim == 3 :
                        fields[irec,  :,:] = ( field - field_prev )/dtsec
                    elif ndim == 4 :
                        fields[irec,:,:,:] = ( field - field_prev )/dtsec
                    else :
                        logging.error( 'could not reset field for ndim=%i' % ndim )
                        raise Exception
                    #endif
                    # store current values as next previous:
                    step_prev = step
                    field_prev = field.copy()
                    # next record:
                    irec = irec + 1
                #endfor
            #endfor
            # clip if necessary:
            if (vmin is not None) or (vmax is not None) :
                # info ...
                logging.debug( indent+'    cliping interval: [%s,%s]' % (vmin,vmax) )
                # truncate:
                fields = fields.clip(vmin,vmax)
            #endif
            # packed?
            if 'add_offset' in varid.ncattrs() :
                # info on packed data type:
                info = numpy.iinfo(varid.dtype)
                # new offset and scale factor ;
                # enlarge factor with 1% to avoid values juste outside range after packing:
                offset = 0.5*( fields.max() + fields.min() )
                scale  = 1.01*( fields.max() - fields.min() )/( info.max - info.min )
                # reset attributes:
                varid.setncattr( 'add_offset', offset )
                varid.setncattr( 'scale_factor', scale )
            #endif
            # write, packed automatically if packing attributes are present:
            varid[:] = fields
            # add new time bounds ?
            tbname = '%s_bounds' % tname
            if tbname not in self.ncid.variables.keys() :
                # add dimension for coordinate bounds if necessary:
                bname = 'bnds'
                if bname not in self.ncid.dimensions.keys() : self.createDimension( bname, 2 )
                # create new variable for time bounds:
                varid_time_bounds = self.ncid.createVariable( tbname, varid_time.dtype, (tname,bname) )
                # fill; round to integers here, otherwise netcdf lib truncates:
                varid_time_bounds[:] = numpy.round( netCDF4.date2num( numpy.array(time_bounds), time_units ) )
                varid_time_bounds.setncattr( 'units', time_units )
                varid_time_bounds.setncattr( 'calendar', varid_time.getncattr('calendar') )
                # add reference:
                varid_time.setncattr( 'bounds', tbname )
            #endif
            # change attributes:
            aname = 'units'
            if aname in varid.ncattrs() :
                avalue = varid.getncattr(aname)
                if avalue == 's' :
                    newvalue = '1'
                elif avalue.endswith(' s') :
                    newvalue = avalue[0:-2]
                else :
                    newvalue = avalue+' s**-1'
                #endif
                varid.setncattr( aname, newvalue )
            #endif
        #endfor  # data variables

        # extend history:
        self.ExtendHistory( 'c3po: converted from accumulated to time average' )

    #enddef   # DeAccum



#endclass   # EcGridFile
            


#######################################################################
###
### test
###
#######################################################################

if __name__ == "__main__" :
    
    # modules
    import os
    import shutil
    import logging
    
    # setup logging:
    logging.basicConfig( level=logging.DEBUG )
    
    # data directory:
    datadir = '/scratch/segers/tmp/test-leip/data'
    
    # which file?
    infile = os.path.join( datadir, 'eurext_lnsp_20060701_fc00.nc' )
    # copy:
    outfile = os.path.join( datadir, 'lnsp.nc' )
    
    # copy:
    if os.path.isfile(outfile) : os.remove( outfile )
    shutil.copy( infile, outfile )
    
    # open:
    ecfile = EcGridFile( outfile )
    # rename variables and add approriate attributes
    # using parameter tables passed with GRIB_API:
    ecfile.ComplyVariables()
    # done:
    ecfile.close()

#endif   # main


#######################################################################
###
### end
###
#######################################################################


    
    

