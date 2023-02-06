# AnanseSeurat 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* updated code formatting to adhere to the bestpractices() package
* updated Readme file
* added a smaller test file for the example of DEG_function to test within the allowed timeslot
* changed print statements to message(), warning() or stop() statements
* removed 'outputdir' variable from Maelstrom_Motif2TF, since it was not used and data is
  now returned to the S4 objects as assays (and no longer as output files)
* changed output_dir in all the export functions to be a mandatory argument instead of taking 
  the working dir as a default.
* added apache license to include patent info

