OpenLE meteo tools
==================

Scripts to extract OpenLE meteorological data from the ECMWF archive,
and to transfer the files to your local computer.

An ECMWF member state account is required.


Extract data from ECMWF archive
-------------------------------

The best (fastest) way to extract meteorological data from the ECMWF archive
is to retrieve it from the MARS archive on the member state server 'ecgate'
(https://software.ecmwf.int/wiki/display/UDOC/ecgate).
After copying this 'tools/meteo' directory to ecgate,
edit and run the following script to create data files:
  
  ./bin/OpenLE-meteo-ecmwf
  
This will first retrieve files from MARS in grib format,
and then convert to netcdf. The standard grib_to_netcdf tool at ecgate
does not provide completely CF-compliant files, and therefore small
additional modifications are needed.

Variables that are originally accumulated from the start of a forecast
are converted to averages over an interval. For example, precipitation 
fields are originally stored in 'm water equivalent' since start,
but converted to 'm/s' average during a (3 hourly) time interval.

The data volume is reduced significantly by combining vertical layers.
How to combine levels is defined in the header of the script;
for OpenLE we created files with about 20 layers in the troposphere.
To combine the levels a Fortran program is compiled that is present
in the 'src' directory. Compilation settings are defined in:
  rc/openle-meteo.rc
  
Extracting data could take a long time! Especially fields on ecmwf model levels
(temperature, wind) are expensive to extract. The ecgate server does not
allow jobs to run in foreground for longer than 30 minutes; longer jobs should
be committed to the queue system. Use the following script to create a
sequence of jobs that will be submitted one by one to the queue system:
  ./bin/OpenLE-meteo-ecmwf-jobs
Edit this file to set the time range and select the long/short jobs.
To see which jobs are currently running on your account, use:
  squeue -u $USER
The job-id number in the first column could be used to cancel a running job:
  scancel <JOBID>


Transfer meteo files from ecgate to local computer
--------------------------------------------------

Meteo files produced at ecgate could be transfered to a local computer
using the script:

  ./bin/OpenLE-meteo-get
  
The script uses the 'ecaccess' tools to scan a source directory on ecgate
and get the files to the local machine if not present yet
(https://software.ecmwf.int/wiki/display/ECAC/).
Edit the settings in the top of this script to change locations etc.

