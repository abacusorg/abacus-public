This directory contains the Direct IO library used by singlestep and convolution
(as well as some analysis routines).  Direct IO is used to read files from disk
without filesystem or OS caching.

Some of the classes here are not strictly Direct IO related and are instead
general file/IO utils.  I think it's a good idea to collect these in one
place, though.
