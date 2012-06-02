The "transcript-align.sh" script is intended for use under
Unix/Linux systems but should also work under Cygwin.

To use the transcript-align.sh script you must do the following:

1. In the directory with the script, place a goby.jar and a
   goby.properties file that specifies the appropriate locations
   for the aligners (such as bwa) that you intend to use.

2. Edit the top of the transcript-align.sh script to specify the
   reference you want to use for the alignment
   (should be a .fa.gz file), the organism and version for that
   reference, and the reads file you want to align
   (should be a .compact-reads file). Also make sure you have the
   aligner specified and that you have specified the tag for this
   alignment (a unique string with no spaces that is meaningful
   to you).

3. Execute the script.

If all went well, a new directory named "result-TAG" should have
been created that contains the alignment.
