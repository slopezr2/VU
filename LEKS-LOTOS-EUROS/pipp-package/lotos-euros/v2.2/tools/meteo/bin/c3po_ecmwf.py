"""
C3PO - CF Convention Compliance Python Objects

ECMWF data specific conversions.

Main module that provides access to other modules:

  tables     -   Read grib definition and other tables.
  
  gribfile   -   Tools to convert GRIB files to NetCDF
                 using the Python interface in GRIB_API.

  gridfile   -   Tools to comply NetCDF files created by ECMWF
                 standard softare, i.e. mars2netcdf or cdo
                 towards CF standard.
  
"""

# sub modules:
#import c3po_ecmwf_tables   as tables
#import c3po_ecmwf_gribfile as gribfile
import c3po_ecmwf_ncfile   as ncfile
