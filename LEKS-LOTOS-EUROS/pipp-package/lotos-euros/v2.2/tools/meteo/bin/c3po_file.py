#! /usr/bin/env python3

#######################################################################
###
### main
###
#######################################################################

"""
C3PO - CF Convention Compliance Python Objects

CLASSES
  
    Define classes for NetCDF files with meteo:
    
        netCDF4.Dataset   -   base class to access NetCDF files

            GridFile          -   general file with gridded data
            
                ReducedGridFile   -  special case with reduced grid representation

    The classes have routines defined that convert the content towards 
    full CF-Convention compliance.

HISTORY

    2013-01, Arjo Segers

"""


#######################################################################
###
### modules
###
#######################################################################

# modules:
import netCDF4


#######################################################################
###
### base class for gridded file
###
#######################################################################

# define clas:
class GridFile( object ) :
    
    """
    Access NetCDF file for conversion towards CF-compliance.

    Usage:
    
        (see routine specific help for details)

        # open existing file for read/write:
        ncfile = GridFile( 'myfile.nc' )
        
        # add global attributes defined in dictionairy:
        ncfile.SetGlobalAttributes( { 'name' : value, ... } )
        
        # add or change attribute, return boolean flag:
        is_added = ncfile.DefaultAttribute( varid, 'standard_name', 'pressure' )

        # extend history:
        ncfile.ExtendHistory( 'Added CF standard name attributes.' )

        # remove dimension:
        ncfile.RemoveDimension( 'time' )

        # add bound variables to coordinate variable:
        ncfile.AddBounds( 'longitude' )
        
        # define time dimension(s) following conventions in GRIB files:
        ncfile.DefineTime( tref, t1, t2, 'Forecast product valid at reference time + P1 (P1>0)' )

        # close:
        ncfile.close()

    Expert usage:    

        # open file:
        ncfile = GridFile( 'myfile.nc' )
        # close, and create a temporary copy:
        tmpfile = ncfile.CreateTemporaryCopy( 'myfile_tmp.nc' )
        # close:
        tmpfile.close()

    """
    
    
    # ***
    
    # define names used by other routines:
    def GetTimeName  ( self ) : return 'time'
    def GetBoundsName( self ) : return 'nv'
    def GetBoundsExt ( self ) : return '_bnds'
    
    
    # ***
    
    def __init__( self, filename, mode='r+', format=None ) :
        
        """
        Open NetCDF file. Default mode is to open for read/write ; 
        use mode='w' to (re)create an empty file.
        """
        
        # open file according to acces mode:
        self.ncid = netCDF4.Dataset( filename, mode=mode, format=format )
        
        # store filename:
        self.filename = filename
        
    #enddef
    
    
    # ***
    
    
    def close( self ) :

        """
        Close file.
        """
        
        # close:
        self.ncid.close()

    #enddef close
    
    
    # ***
    
    def createDimension( self, name, length ) :
    
        """
        Check if dimension with name is present, and if it has the same lengt;
        otherwise create new dimension.
        """
        
        # modules:
        import logging
        
        # check on existence first ...
        if name in self.ncid.dimensions.keys() :
            # current lengt:
            length_curr = len(self.ncid.dimensions[name])
            # chec
            if length_curr != length :
                logging.error( 'file "%s" has already dimenstion "%s" with length %i instead of %i' % (self.filename,name,length_curr,length) )
                raise Exception
            #endif
        else :
            # new:
            self.ncid.createDimension(  name, length )
        #endif
        
    #enddef
    
    
    # ***
    
    
    def CreateTemporaryCopy( self, tmpfilename ) :
    
        """
        Close file, rename, and open again for input.
        Create a new file with the original name for writing,
        copy already the global attributes.
        Useful for copying parts from the current file.
        Retruns the netCDF4.Dataset object of the opened temporary copy.
        """
        
        # modules:
        import os
        import netCDF4
    
        # keep own name:
        filename = self.filename
        # close current file:
        self.close()
        # rename:
        os.rename( filename, tmpfilename )
        # open copy:
        tempf = netCDF4.Dataset( tmpfilename, 'r' )

        # new file:
        GridFile.__init__( self, filename, mode='w', format=tempf.file_format )
        
        # copy global attributes:
        for attrname in tempf.ncattrs() :
            # current value:
            value = tempf.getncattr(attrname)
            # add to new file:
            self.ncid.setncattr( attrname, value )
        #endfor
        
        # ok
        return tempf
        
    #enddef  # CreateTemporaryCopy


    # ***
    
    
    def SetGlobalAttributes( self, attrs ) :
    
        """
        Add global attributes specified by dictionairy.
        """
        
        # loop:
        for key in attrs.keys() :
            # add:
            self.ncid.setncattr( key, attrs[key] )
        #endfor
        
    #enddef  # SetGlobalAttributes
    
    # *
    
    def GetGlobalAttributes( self ) :
    
        """
        Return global attributes.
        """
        
        # init result:
        attrs = {}
        # loop:
        for key in self.ncid.ncattrs() :
            # extract:
            attrs[key] = self.ncid.getncattr( key )
        #endfor
        
        # ok
        return attrs
        
    #enddef  # SetGlobalAttributes
    
    # *
    
    def GetGlobalAttribute( self, name ) :
    
        """
        Return global attribute with given name.
        """
        
        # modules:
        import logging
        
        # check ...
        if name not in self.ncid.ncattrs() :
            logging.error( 'attribute "%s" not found in %s' % (name,self.filename) )
            raise Exception
        #endif
        
        # extract:
        return self.ncid.getncattr(name)
        
    #enddef  # GetGlobalAttribute


    # ***
    
    
    def DefaultAttribute( self, varid, attrname, attrvalue, indent='' ) :
    
        """
        Add new attribute with specified value if it does not exist yet.
        
        Return value:
            boolean flag which is True if the the attribute has been added.
        """
        
        # modules:
        import logging

        # default:
        changed = False
        # not defined yet ?
        if attrname not in varid.ncattrs() :
            # info ...
            logging.debug( indent+'add attribute "%s" = "%s"' % (attrname,attrvalue) )
            # add:
            varid.setncattr( attrname, attrvalue )
            # set flag:
            changed = True
        #endif
        
        # ok
        return changed
        
    #enddef  # DefaultAttribute
    

    # ***


    def ExtendHistory( self, msg, t=None, tstamp=None, old=False ) :

        """
        Add new message to history global attribute.

        Each line in the history starts with a time stamp.
        If no timestamp is provided, it is formed from the datetime
        object 't', which defaults to the current time.

        New messages are prepend in the history, unless old=True
        in which case the message is placed at the end,
        for example to desecribe the histroy before conversion to netcdf.
        """

        # modules:
        import datetime
        
        # tools:
        import c3po_tools
        
        # get history line (time label and message):
        line = c3po_tools.HistoryLine( msg, t=t, tstamp=tstamp )

        # fill history, already existing?
        if 'history' in self.ncid.ncattrs() :
            # get current:
            history = self.ncid.getncattr('history')
            # old or new ?
            if old :
                history = '%s\n%s' % (history,line)
            else :
                history = '%s\n%s' % (line,history)
            #endif
        else :
            # copy:
            history = line
        #endif
        # write:
        self.ncid.setncattr( 'history', history )

    #enddef  ExtendHistory
    
    
    # ***
    
    
    def RemoveDimension( self, name, indent='', cleanup=True, allow_missing=False ) :

        """
        Remove named dimension and related coordinate variables;
        the dimension to be removed should have length 1.

        Dimensions and/or variables cannot be removed from a file.
        Therefore, the open file is first closed and a temporary copy is created
        with extension '_with_<name>'. The required variables/attributes are then copied
        from there into a newly created version of the current file,
        Set optional argument 'cleanup' to 'False' to keep the temporary copy.

        Specifiy optional whitespace to indent logging messages.
        """
        
        # modules:
        import os
        import logging
        
        # info ...
        logging.debug( indent+'remove singleton dimension %s ...' % name )
        
        # check if dimension is present:
        if name not in self.ncid.dimensions.keys() :
            if allow_missing :
                # just leave ...
                return
            else :
                # error ...
                logging.error( 'dimension "%s" not found in "%s"' % (name,self.filename) )
                raise Exception
            #endif
        #endif
        
        # get length:
        length = len(self.ncid.dimensions[name])
        ## should be singleton dimension:
        #if length != 1 :
        #    logging.error( 'could only remove dimension if it has length 1 :' )
        #    logging.error( 'input file : %s' % self.filename )
        #    logging.error( 'dimension  : %s (%i)' % (name,length) )
        #    raise Exception
        ##endif

        # name of temporary copy:
        dimfile = self.filename.replace( '.nc', '_with_%s.nc' % name )
        # create copy, return its dimension id, and create new empty file:
        dimf = self.CreateTemporaryCopy( dimfile )

        # copy dimensions:
        for dimname in dimf.dimensions.keys() :
            # skip dimension to be removed:
            if dimname == name : continue
            # original:
            dimid = dimf.dimensions[dimname]
            # length:
            dimlen = len(dimid)
            if dimid.isunlimited() : dimlen = 0
            # create:
            self.ncid.createDimension( dimname, dimlen )
        #endfor

        # loop over variables:
        for varname in dimf.variables.keys() :
            # skip coordinate variable:
            if varname == name : continue
            # original:
            varid = dimf.variables[varname]
            # skip variables with only removed dimension:
            if (name in varid.dimensions) and (len(varid.dimensions) == 1) : continue
            # now only length 1 allowed ...
            if (name in varid.dimensions) and (length != 1) :
                logging.error( 'could only remove dimension from multi-dimensional data if it has length 1 :' )
                logging.error( 'input file : %s' % self.filename )
                logging.error( 'dimension  : %s (%i)' % (name,length) )
                logging.error( 'variable   : %s %s' % (varname,str(varid.dimensions)) )
                raise Exception
            #endif
            # copy dimension list, skip dimension if necesary:
            newdims = ()
            for dim in varid.dimensions :
                if dim != name : newdims = newdims+(dim,)
            #endif
            # create copy:
            newvarid = self.ncid.createVariable( varname, varid.dtype, newdims )
            # attributes:
            for aname in varid.ncattrs() : newvarid.setncattr( aname, varid.getncattr(aname) )
            # copy data:
            newvarid[:] = varid[:]
        #endfor

        # add line to history:
        self.ExtendHistory( 'c3po: removed singleton dimension %s' % (name) )

        # close temporary copy:
        dimf.close()

        # remove temporary file ?
        if cleanup : os.remove( dimfile )

    #enddef   # RemoveDimension


    # ***
    

    def AddCoordinateVariable( self, name, values, attr={}, dtype=None, indent='' ) :

        """
        Add coordinate variable with data.
        
        Optional arguments:
          attr      : optional attributes
          indent    : specifiy whitespace to indent logging messages
        """

        # modules:
        import numpy
        import logging
        
        # info ...
        logging.debug( indent+'add variable "%s" ...' % name )
        
        # should be coordinate variable, thus dimension with same name:
        if name not in self.ncid.dimensions.keys() :
            logging.error( 'dimension "%s" not found in "%s"' % (name,self.filename) )
            raise Exception
        #endif
        
        # data type:
        if dtype == None : dtype = values.dtype
        
        # create:
        varid = self.ncid.createVariable( name, dtype, (name,) )
        # fill data:
        varid[:] = values

        # attributes:
        for key in attr.keys() : varid.setncattr( key, attr[key] )

    #enddef  # AddCoordinateVariable


    # ***
    

    def AddBounds( self, name, bounds=None, indent='', extend_history=True ) :

        """
        Add bounds for a coordinate variable.
        Guess from center values.
        
        Optional arguments:
          bounds    : (n,2) array with the values to be used instead of guessing
          indent    : specifiy whitespace to indent logging messages
        """

        # modules:
        import numpy
        import logging
        
        # info ...
        logging.debug( indent+'add boundary arrays for variable "%s" ...' % name )
        
        # check if variable exists:
        if name not in self.ncid.variables.keys() :
            logging.warning( 'varible "%s" not found in "%s"' % (name,self.filename) )
            return
        #endif
        
        # should be coordinate variable, thus dimension with same name:
        if name not in self.ncid.dimensions.keys() :
            logging.error( 'dimension "%s" not found in "%s"' % (name,self.filename) )
            raise Exception
        #endif

        # add boundary dimension if necessary:
        bname = self.GetBoundsName()
        nv = 2
        # add if not present yet:
        self.createDimension( bname, nv )

        # variable:
        varid = self.ncid.variables[name]
        
        # guess boundary values ?
        if bounds == None :
        
            # center values:
            vals = varid[:]
            # length:
            nval = len(vals)
            # boundaries between center values:        
            bvals = numpy.interp( numpy.arange(nval-1)+0.5, numpy.arange(nval), vals )

            # boundary range:
            lower = vals[0] - (bvals[0]-vals[0])
            upper = vals[-1] + (vals[-1]-bvals[-1])

            # has standard name ?
            if 'standard_name' in varid.ncattrs() :
                # extract:
                standard_name = varid.getncattr('standard_name')
                # switch:
                if standard_name == 'latitude' :
                    # truncate latitudes to [-90,90] :
                    lower = min( max( -90.0, lower ), 90.0 )
                    upper = min( max( -90.0, upper ), 90.0 )
                #endif
            #endif

            # 2d bounds array:
            bounds = numpy.zeros((nval,nv),vals.dtype)
            # fill lower bounds:
            bounds[     0  ,0] = lower
            bounds[1:nval  ,0] = bvals
            # fill upper bounds:
            bounds[0:nval-1,1] = bvals
            bounds[  nval-1,1] = upper

        #endif # create bounds

        # name of new variable:
        name_bounds = name+'_bounds'
        # create new variable:
        varid_bounds = self.ncid.createVariable( name_bounds, bounds.dtype, (name,bname) )
        # loop over original attributes:
        for aname in varid.ncattrs() :
            # skip some:
            if aname in ['axis','calendar','positive','formula','formula_terms'] : continue
            # copy:
            varid_bounds.setncattr( aname, varid.getncattr(aname) )
        #endfor
        # fill:
        varid_bounds[:] = bounds
        # add attribute to original variable:
        varid.setncattr( 'bounds', name_bounds )

        # add line to history:
        if extend_history :
            self.ExtendHistory( 'c3po: added boundaries for coordinate variable %s' % (name) )
        #endif

    #enddef  # AddBounds
    
    
    # *
    
    def AddCellArea( self, lonname, latname ) :
    
        """
        Compute cell area field given lon and lat bounds;
        dims and names are supposed to be:
           <lonname> = 123
           <lonname>_bounds(<lonname>,nv)
        Return values:
          name of cell area variable created, use in attribute:
             cell_measures = "area: cell_area"
          area field
        """
        
        # modules:
        import numpy
        
        # tools:
        import grid
        
        # dims:
        nlon = len(self.ncid.dimensions[lonname])
        nlat = len(self.ncid.dimensions[latname])
        # bounds:
        lon_bounds = self.ncid.variables[lonname+'_bounds'][:]
        lat_bounds = self.ncid.variables[latname+'_bounds'][:]
        
        # create field:
        area = numpy.zeros((nlat,nlon),float)
        # loop over lat bands:
        for j in range(nlat) :
            # use bounds for first lon, others are same:
            area[j,:] = grid.ll_area( lon_bounds[0,0], lon_bounds[0,1],
                                      lat_bounds[j,0], lat_bounds[j,1]  )
        #endfor
        
        # target variable:
        vname = 'cell_area'
        # write:
        self.AddVariable( vname, (latname,lonname), area,
                            attr={ 'standard_name' : 'cell_area', 'units' : 'm2' } )
        
        # ok
        return vname,area
        
    #enddef  AddCellArea
        
    

    # ***
    

    def AddVariable( self, name, dimnames, values, dtype='f4', attr={}, indent='' ) :

        """
        Add variable with data.
        
        Optional arguments:
          dtype     : 'f4' | 'i4' | ...
                      if not defined, value.dtype is used
          attr      : optional attributes
          indent    : specifiy whitespace to indent logging messages
        """

        # modules:
        import logging
        
        # info ...
        logging.debug( indent+'add variable "%s" ...' % name )
        
        # data type:
        if dtype == None : dtype = values.dtype
        
        # create::
        varid = self.ncid.createVariable( name, dtype, dimnames )
        # fill data:
        if values != None :
            if dtype == 'S1' :
                varid[:] = netCDF4.stringtochar(values)
            else :
                varid[:] = values
            #endif
        #endif

        # attributes:
        for key in attr.keys() : varid.setncattr( key, attr[key] )

    #enddef  # AddVariable

    
    def DefineTime( self, tref, t1, t2, longname ) :
        
        """
        Define time dimension(s) and time coordinate variables given lists of values
        
        The meaning of the time variables is given by a character string 'longname'.
        Possible values are those in grib table '5' selected following the 'timeRangeIndicator'
        in a grib file sections.
        Currently supported:
            'Forecast product valid at reference time + P1 (P1>0)'

        Arguments:
            tref       :  reference time values
            t1         :  first time values
            t2         :  second time values
            longname   :  character string describing meaning of (tref,t1,t2)
          
        """
        
        # modules:
        import logging
        import numpy
        
        # time length:
        nt = len(tref)
        
        # check lengths ...
        if (len(t1) != nt) or (len(t2) != nt) :
            logging.error( 'length of time arrays should be the same: tref (%i), t1 (%i), t2 (%i)' % (nt,len(t1),len(t2)) )
            raise Exception
        #endif
        
        # check if times are increasing:
        for k in range(2,nt) :
            if tref[k] < tref[k-1] :
                logging.error( 'reference times should be increasing:' )
                for i in range(nt) : logging.error( '  tref[%i] = %s' % (i,tref[i]) )
                raise Exception
            #endif
            if t1[k] < t1[k-1] :
                logging.error( 'start times should be increasing:' )
                for i in range(nt) : logging.error( '  t1[%i] = %s' % (i,t1[i]) )
                raise Exception
            #endif
            if tref[k] < tref[k-1] :
                logging.error( 'end times should be increasing:' )
                for i in range(nt) : logging.error( '  t2[%i] = %s' % (i,t2[i]) )
                raise Exception
            #endif
        #endfor

        # check if t2 is after t1 :
        for k in range(2,nt) :
            if t2[k] < t1[k] :
                logging.error( 'end of time interval should be after start:' )
                for i in range(nt) : logging.error( '  (t1,t2)[%i] = (%s,%s)' % (i,t1[i],t2[i]) )
                raise Exception
            #endif
        #endfor

        # name used for time dimension:
        tname = self.GetTimeName()
        # define dimension:
        self.createDimension( tname, nt )

        # store in long format:
        tdtype = 'f8'
        # store as seconds since first reference time:
        tstepunits = 'seconds'
        tunits = '%s since %s' % (tstepunits,tref[0])

        # switch:
        if longname == 'Forecast product valid at reference time + P1 (P1>0)' :

            # check that t2 == t1 ..
            for k in range(nt) :
                if t2[k] != t1[k] :
                    logging.error( 'product should be instantaneous:' )
                    for i in range(nt) : logging.error( '  (t1,t2)[%i] = (%s,%s)' % (i,t1[i],t2[i]) )
                    raise Exception
                #endif
            #endfor

            # define reference time:
            varid = self.ncid.createVariable( '%s_ref' % tname, tdtype, (tname,) )
            # annote:
            varid.setncattr( 'standard_name', 'time' )
            varid.setncattr( 'longname', 'reference time of forecast' )
            varid.setncattr( 'units', tunits )
            #varid.setncattr( 'calendar', 'standard' )  # don't include this in reftime, otherwise cf checker will complain
            # fill:
            varid[:] = numpy.round( netCDF4.date2num( tref, tunits ) )

            # time step P1
            varid = self.ncid.createVariable( '%s_P1' % tname, tdtype, (tname,) )
            # annote:
            varid.setncattr( 'standard_name', 'time' )
            varid.setncattr( 'longname', 'time step of forecast' )
            varid.setncattr( 'units', tstepunits )
            # compute:
            P1 = []
            for k in range(nt) : P1.append( (t1[k]-tref[k]).total_seconds() )
            # fill:
            varid[:] = P1

            # actual time:
            varid = self.ncid.createVariable( tname, tdtype, (tname,) )
            # annote:
            varid.setncattr( 'standard_name', 'time' )
            varid.setncattr( 'longname', longname )
            varid.setncattr( 'units', tunits )
            varid.setncattr( 'calendar', 'standard' )
            # fill:
            varid[:] = numpy.round( netCDF4.date2num( t1, tunits ) )
            
        #elif ... accumulated field ...
        #
        #    add attribute:
        #      cell_methods = 'time: sum'
        #
    
        else :
            
            logging.error( 'unsupported timing longname : %s' % longname )
            raise Exception
            
        #endif  # switch over longnames
       
    #enddef   # DefineTime
    
    # *
    
    def createVariableCopy( self, name, varid, dimensions=None,
                              data=False, selection=None ) :
    
        """
        Create new variable with provided name,
        and data type, dimensions, and attributes of the input variable. 
        A changed set of dimension names is optional.
        Eventually copy data too, eventually with the given selection.
        Returns new variable object.
        """
        
        # new dimensions:
        if dimensions is None : dimensions = varid.dimensions

        # create new variable:
        v = self.ncid.createVariable( name, varid.dtype, dimensions )
        # copy attributes:
        for key in varid.ncattrs() :
            v.setncattr( key, varid.getncattr(key) )
        #endfor
        
        # copy data?
        if data : 
            # copy:
            if selection is None :
                v[:] = varid[:]
            else :
                v[:] = varid[selection]
            #endif
        #endif
        
        # ok
        return v
        
    #enddef createVariableCopy


#endclass   # GridFile
    



#######################################################################
###
### base class for reduced gridded file
###
#######################################################################
            

# define class:
class ReducedGridFile( GridFile ) :
    
    """
    Write/read reduced grid data to NetCDF file.
    
    Definition following CF convention document:
      http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/ch05s03.html

    Example of NetCDF header:
 
      dimensions:
        latdim =   160 ;    // number of latitude bands
        londim =   320 ;    // number of longitudes in unreduced grid
        rgrid  = 35718 ;    // number of grid points
      variables:
        // zero based index that specifies location in regular (latdim,londim) grid :
        //   j = rgrid(n) / londim
        //   i = rgrid(n) - londim * j
        int rgrid(rgrid);
          rgrid:compress    = "latdim londim";
        // longitude for each point:
        float lon(rgrid) ;
          lon:standard_name = "longitude" ;
          lon:long_name     = "longitude" ;
          lon:units         = "degrees_east" ;
        float lat(rgrid) ;
          lat:standard_name = "latitude" ;
          lat:long_name     = "latitude" ;
          lat:units         = "degrees_north" ;

    Example of variable defined on this grid:

        float PS(rgrid) ;
          PS:long_name = "surface pressure" ;
          PS:units = "Pa" ;
          PS:coordinates = "lon lat" ;
          
    Extra facilities are present compared to the base class "GridFile" :
    
        # create and fill grid variables passed in object 'rgg'
        # of 'grid.rgg.ReducedGrid' class:
        ncid.DefineGrid( rgg )
        
        # create new empty gridded variable:
        ncid.CreateGridVariable( 'name', *args )
        
        # fill (slice of) data in named variable:
        ncid.AddValues( 'name', values, itime=2, ilevel=3 )
    
    """
    
    # * names:
    
    def GetRGridName ( self ) : return 'rgrid'
    def GetLonDimName( self ) : return 'ulon'
    def GetLatDimName( self ) : return 'ulat'
    def GetLonName   ( self ) : return 'longitude'
    def GetLatName   ( self ) : return 'latitude'
    
    # *
    
    def DefineGrid( self, rgg ) :
        
        """
        Create dimensions and variables that define the grid.
        Arguments:
          rgg   :  ReducedGaussianGrid object.
        """
        
        # add dimensions:
        self.createDimension( self.GetLatDimName(), rgg.nlat )
        self.createDimension( self.GetLonDimName(), rgg.unlon )
        self.createDimension( self.GetRGridName(), rgg.npoint )
        self.createDimension( self.GetBoundsName(), 2 )
        
        # ~ index
        
        # add index variable:
        varid = self.ncid.createVariable( self.GetRGridName(), 'i4', (self.GetRGridName(),) )
        varid.setncattr( 'long_name', 'reduced grid index' )
        varid.setncattr( 'compress', '%s %s' % (self.GetLatDimName(),self.GetLonDimName()) )
        varid.setncattr( 'description', 'zero based index that specifies location in regular (%s,%s) grid :  j = rgrid(n) / %s ; i = rgrid(n) - %s * j' % (self.GetLatDimName(),self.GetLonDimName(),self.GetLonDimName(),self.GetLonDimName()) )
        varid[:] = rgg.irgrid
        
        # ~ centers
        
        # longitude values for each point:
        varid = self.ncid.createVariable( self.GetLonName(), 'f4', (self.GetRGridName(),) )
        varid.setncattr( 'standard_name', 'longitude' )
        varid.setncattr( 'long_name', 'longitude' )
        varid.setncattr( 'units', 'degrees_east' )
        varid.setncattr( 'bounds', self.GetLonName()+self.GetBoundsExt() )
        varid[:] = rgg.longitude

        # latitude values for each point:
        varid = self.ncid.createVariable( self.GetLatName(), 'f4', (self.GetRGridName(),) )
        varid.setncattr( 'standard_name', 'latitude' )
        varid.setncattr( 'long_name', 'latitude' )
        varid.setncattr( 'units', 'degrees_north' )
        varid.setncattr( 'bounds', self.GetLatName()+self.GetBoundsExt() )
        varid[:] = rgg.latitude
        
        # ~ boundary values
        
        # longitude bondary values for each point:
        varid = self.ncid.createVariable( self.GetLonName()+self.GetBoundsExt(), 'f4', 
                                              (self.GetRGridName(),self.GetBoundsName()) )
        varid.setncattr( 'long_name', 'longitude cell boundaries' )
        varid.setncattr( 'units', 'degrees_east' )
        varid[:] = rgg.longitude_bnds

        # latitude boundary values for each point:
        varid = self.ncid.createVariable( self.GetLatName()+self.GetBoundsExt(), 'f4', 
                                              (self.GetRGridName(),self.GetBoundsName()) )
        varid.setncattr( 'long_name', 'latitude cell boundaries' )
        varid.setncattr( 'units', 'degrees_north' )
        varid[:] = rgg.latitude_bnds
        
        # ~ latitude bands (extra info)
        
        # latitude per band:
        varid = self.ncid.createVariable( self.GetLatDimName(), 'f4', (self.GetLatDimName(),) )
        varid.setncattr( 'standard_name', 'latitude' )
        varid.setncattr( 'long_name', 'latitude' )
        varid.setncattr( 'units', 'degrees_north' )
        varid.setncattr( 'bounds', self.GetLatDimName()+self.GetBoundsExt() )
        varid[:] = rgg.band_lats
        
        # latitude boundaries per band:
        varid = self.ncid.createVariable( self.GetLatDimName()+self.GetBoundsExt(), 'f4', 
                                              (self.GetLatDimName(),self.GetBoundsName()) )
        varid.setncattr( 'long_name', 'latitude cell boundaries' )
        varid.setncattr( 'units', 'degrees_north' )
        varid[:] = rgg.band_lats_bnds
        
        # number of longitudes per row:
        varid = self.ncid.createVariable( self.GetLatDimName()+'_nlon', 'i2', (self.GetLatDimName(),) )
        varid.setncattr( 'long_name', 'number of longitudes on latitude band' )
        varid.setncattr( 'units', '1' )
        varid[:] = rgg.band_nlon
        
        ## first longitude on row:
        #varid = self.ncid.createVariable( self.GetLatDimName()+'_'+'lon0', 'f4', (self.GetLatDimName(),) )
        #varid.setncattr( 'standard_name', 'longitude' )
        #varid.setncattr( 'long_name', 'longitude of first point on latitude band' )
        #varid.setncattr( 'units', 'degrees_east' )
        #varid[:] = rgg.band_lon0
        
        # index in reduced grid:
        varid = self.ncid.createVariable( self.GetLatDimName()+'_i0', 'i4', (self.GetLatDimName(),) )
        varid.setncattr( 'long_name', 'zero-based index in reduced grid of first point on latitude band' )
        varid.setncattr( 'units', '1' )
        varid[:] = rgg.band_i0
        
    #enddef  # DefineGrid
    
    # *
    
    def GetGridDefinition( self ) :
        
        """
        Read grid variable and return  'grid_rgg.ReducedGrid' ojbect with defintion.
        """
        
        # modules:
        import logging
        
        # tools:
        import grid_rgg
        
        # dimensions:
        vname = self.GetLatDimName()
        nlat = len(self.ncid.dimensions[vname])
        
        # latitude per band:
        vname = self.GetLatDimName()
        band_lats = self.ncid.variables[vname][:]
        
        # latitude boundaries per band:
        vname = self.GetLatDimName()+self.GetBoundsExt()
        band_lats_bnds = self.ncid.variables[vname][:]
        
        # longitude values for each point:
        vname = self.GetLonName()
        longitude = self.ncid.variables[vname][:]
        
        # longitude bondary values for each point:
        vname = self.GetLonName()+self.GetBoundsExt()
        longitude_bnds = self.ncid.variables[vname][:]
        
        # number of longitudes per row:
        vname = self.GetLatDimName()+'_nlon'
        band_nlon = self.ncid.variables[vname][:]
        
        # create lists with for each latitude band the longitude arrays:
        band_lons = []
        band_lons_bnds = []
        # init base index:
        i0 = 0
        # loop over lat bands:
        for j in range(nlat) :
            # number of longitudes in this band:
            nlon = band_nlon[j]
            # index range:
            i1 = i0 + nlon
            # copy:
            band_lons.append( longitude[i0:i1] )
            band_lons_bnds.append( longitude_bnds[i0:i1,:] )
            # shift:
            i0 = i1
        #endfor
        # check ...
        if i0 != len(longitude) :
            logging.error( 'after loop i0 = %i while len(longitude) = %i' % (i0,len(longitude)) )
            raise Exception
        #endif
                    
        # define grid definition:
        rgg = grid_rgg.ReducedGrid( band_lats, band_lons,
                                    band_lats_bnds=None, band_lons_bnds=None )
                                    
        # ok
        return rgg
        
    #enddef  # GetGridDefinition
    
    # *
    
    def CreateGridVariable( self, varname, datatype='f4', hastime=False, haslevel=False, attrs={} ) :
    
        """
        Create grided variable. 
        Extra dimensions are leading to the grid.
        Returns variable id.
        """
        
        # dimension list:
        dims = (self.GetRGridName(),)
        #if haslevel : dims = (self.GetLevelName(),)+dims
        if hastime  : dims = (self.GetTimeName(),)+dims
        
        # create:
        varid = self.ncid.createVariable( varname, datatype, dims )
        # add description of auxilary coordinate variables:
        varid.setncattr( 'coordinates', '%s %s' % (self.GetLonName(),self.GetLatName()) )
        # add other attributes:
        for key,value in attrs.iteritems() : varid.setncattr( key, value )
        
    #enddef  # CreateGridVariable
        
    # *
    
    def AddValues( self, shortname, values, itime=None, ilevel=None ) :
    
        """
        Add data.
        """
        
        # time record ?
        if itime == None :
            # level ?
            if ilevel == None :
                # write constant 2D field:
                self.ncid.variables[shortname][:] = values
            else :
                # write level of 3D field:
                self.ncid.variables[shortname][ilevel,:] = values
            #endif
        else :
            # level ?
            if ilevel == None :
                # write time record of 2D field:
                self.ncid.variables[shortname][itime,:] = values
            else :
                # write time record of 3D field:
                self.ncid.variables[shortname][itime,ilevel,:] = values
            #endif
        #endif
        
    #enddef  # AddValues
        
    
#endclass ReducedGridFile





#######################################################################
###
### extract domain
###
#######################################################################

def ExtractDomain( inputfile, outputfile, domain,
                        verbose=False, indent='' ) :

    """
    Extract domain (west,easth,south,north) from input file
    and store in output file.
    """
    
    # modules:
    import logging
    
    # tools:
    import c3po_file
    
    # info ...
    if verbose : logging.info( indent+'open input file %s ..' % inputfile )
    
    # open input file:
    ncid_in = c3po_file.ReducedGridFile( inputfile, 'r', format='NETCDF3_CLASSIC' )
    
    # check ...
    key = 'gridType'
    if key not in ncid_in.ncid.ncattrs() :
        logging.error( 'input file "%s" should have global attribute "%s"' % (inputfile,key) )
        raise Exception
    #endif
    # extract:
    gridType = ncid_in.ncid.getncattr(key)
    # switch:
    if gridType == 'reduced_gg' :
    
        # info ...
        if verbose : logging.info( indent+'read reduced grid defnition ...' )
        # get grid defition from variables:
        rgg = ncid_in.GetGridDefinition()

        # info ...
        if verbose : logging.info( indent+'select domain ...' )
        # select domain:
        rg,inds = rgg.GetRegion( domain )
    
        # info ...
        if verbose : logging.info( indent+'create %s ...' % outputfile )
        # create new file:
        ncid = c3po_file.ReducedGridFile( outputfile, 'w', format='NETCDF3_CLASSIC' )

        # copy global attributes:
        for key in ncid_in.ncid.ncattrs() :
            ncid.ncid.setncattr( key, ncid_in.ncid.getncattr(key) )
        #endif

        # add dimensions etc
        ncid.DefineGrid( rg )
        # name of grid dimension:
        gridname = ncid.GetRGridName()
        
        # copy other dimensions:
        for dimname in ncid_in.ncid.dimensions.keys() :
            # not available yet ?
            if dimname not in ncid.ncid.dimensions.keys() :
                # create:
                ncid.createDimension( dimname, len(ncid_in.ncid.dimensions[dimname]) )
            #endif
        #endfor

        # loop over variables:
        for vname in ncid_in.ncid.variables.keys() :
    
            # already present ? probably a grid variable:
            if vname in ncid.ncid.variables.keys() : continue
    
            # info ...
            if verbose : logging.info( indent+'  write %s ...' % vname )

            # current:
            varidg = ncid_in.ncid.variables[vname]
            dimnames = varidg.dimensions

            # create new variable, copy datatype, dimensions, and attributes:
            varid = ncid.createVariableCopy( vname, varidg )

            # gridded ?
            if gridname in dimnames :
                # switch:
                if len(dimnames) == 1 :
                    # copy selection:
                    varid[:] = varidg[inds]
                elif (len(dimnames) == 2) and (dimnames[1] == gridname) :
                    # copy selection:
                    varid[:,:] = varidg[:,inds]
                elif (len(dimnames) == 3) and (dimnames[2] == gridname) :
                    # copy selection, for some reason does not work for 1 level:
                    if varidg.shape[1] == 1 :
                        varid[:,0,:] = varidg[:,0,inds]
                    else :
                        varid[:,:,:] = varidg[:,:,inds]
                    #endif
                else :
                    # not yet ...
                    logging.error( 'unsupported dimensions "%s" for local selection' % str(dimnames) )
                    raise Exception
                #endif  # shape
            else :
                # copy data:
                varid[:] = varidg[:]
            #endif  # gridded ?

        #endfor  # variables

        # extend history:
        ncid.ExtendHistory( 'c3po: extracted region %s from %s' % (str(domain),inputfile) )
        # close target file:
        ncid.close()

    else :

        # not yet ...
        logging.error( 'unsupported gridType "%s"' % gridType )
        raise Exception

    #endif  # gridType

    # close source file:
    ncid_in.close()
    
#enddef ExtractDomain


#######################################################################
###
### test
###
#######################################################################

# main program ?
if __name__ == "__main__" :

    # modules:
    #import os
    #import sys
    
    ## local version:
    #leipdir = '/home/ms/nl/nl5/work/le/tools/leip/v1.2'
    
    ## local:
    #sys.path.insert( 0, os.path.join(leipdir,'bin') )
    
    # info ...
    print( 'extract local domain ...' )
    
    # split grib file in basename and extenstion:
    bname,ext = os.path.splitext(gribfilename)
    # new file:
    ncfilename = '%s.nc' % bname

    # sub region:
    domain = -15.0,35.0,35.0,70.0
    # target file:
    ncfilename_domain = '%s_domain.nc' % bname
    
    # run:
    ExtractDomain( ncfilename, ncfilename_domain, domain )

#endif  # main program


#######################################################################
###
### end
###
#######################################################################


    
    

