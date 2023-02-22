LOTOS-EUROS model
=================


Quick Start
-----------

LOTOS-EUROS can be run by taking the following steps. This is
just to give a first impression or a quick reminder; for
details, take a look at the User Guide.

1. Go to the required version directory:

      cd lotos-euros/v2.2
     
2. If not present yet, create a link to 'setup_le' script from
   the base/patch version you want to run. 
   For base/patch v2.2.000, this would be:
   
      ln  -s  base/000/bin/setup_le .

3. Copy a template rc-file and give it a suitable name:

      cp  base/000/rc/lotos-euros.rc  lotos-euros-test.rc
      
   Modify it using a text editor, or at least browse through it to see
   if the settings for your run are ok.
   
4. Setup a run-directory and compile an executable using the setup script:

      ./setup_le  lotos-euros-test.rc
      
5. Follow the instructions written by the setup script.
   Probably, the instructions tell you to go to the run directory
   and start the submit script:

      ./submit_le lotos-euros-test.x lotos-euros-test.rc
      
   Alternative: supply the '-s' or '--submit' option 
   to the 'setup_le' script to submit automatically after 
   the setup is finished.
   
   The 'submit_le' script accepts options to let the job run
   in the background or in a queue system; 
   for more info, call the script with argument '--help' .



