LOTOS-EUROS model
=================


Quick Start
-----------

LOTOS-EUROS can be run by taking the following steps. This is
just to give a first impression or a quick reminder; for
details, take a look at the User Guide.

1. Go to the required version directory:

      cd lotos-euros/v2.2

2. Decide which patch you want to run; these are sub-dirs of the base directory:

      base/000
      base/001
             :

3. Best is to not change the (frozen!) base files.
   Therefore, create you own "project" directory to store your changes.
   For clarity, include the patch number in the directory name, 
   so that you know that your changes are a modification on top of this patch:
   
      mkdir -p proj/myproj/001/rc
      
   Copy the main settings to this project directory:
   
      cp base/001/rc/lotos-euros.rc proj/myproj/001/rc

4. As test, just try to setup and start a run with these settings:

      ./base/001/bin/setup-le  proj/myproj/001/rc/lotos-euros.rc
      
   It will probably fail, for example because the settings are not correct for 
   the specific machine that you are using.
   Use a text editor to change the settings where necessary.
   Common changes are probably related to the machine that you use,
   for example the compiler name and the locations of the NetCDF library.
   You might create new machine specific settings as copy of already available machine settings:
   
       cp base/001/rc/machine-tno-hpc3.rc  proj/myproj/001/rc/machine-myinstitute-mymachine.rc
       
   Ensure that these settings are included in:
       proj/myproj/001/rc/lotos-euros.rc
   using lines similar to:

      ! include settings:
      !#include base/${my.le.patch}/rc/machine-tno-hpc3.rc
      #include proj/myproj/${my.le.patch}/rc/machine-myinstitute-mymachine.rc
      
5. Depending on the settings, the model will run in foreground or be submitted to a queue.
   Read the messages displayed.





