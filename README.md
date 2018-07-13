# gaitextract-corrected
Prasanna Sritharan, July 2018

Updates and corrections to GaitExtract171 by Tim Dorn. This is a continual work in progress, updated as needed by me. It may be superseded by more up to date tools in the near future.

Tim Dorn's original implementation of GaitExtract can be found [here](https://simtk.org/projects/c3dtoolbox).

GaitExtract is a Matlab package originally written in 2010 by Tim Dorn ([website](http://timzone.net/)) for his PhD to extract C3D data taken at the Biomotion Laboratory, Mechanical Engineering, University of Melbourne, and also the Australian Institute of Sport indoor track, for subsequent use in OpenSim. Despite being a few years old, GaitExtract still does the job well, but it is temperamental. Feel free to use, but I would suggest users consider similar implementations from other researchers that utilise BTK or more recent versions of C3Dserver.

It utilises a legacy version of C3Dserver to extract C3D point and analog datastreams, and writes OpenSim TRC, MOT and STO tab-delimited text files using the laboratory and model configurations defined in the file LoadLabels.m. It previously could write OpenSim XML setup files, however I have disabled this it included old Matlab p-code files which are no longer supported.
