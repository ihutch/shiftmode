This GIT repo contains the code used to calculate the unstable
dispersion relation of plasma electron holes in two dimensions it may
be used freely provided scholarly research acknowledges the author and
cites the paper "Transverse instability of electron phase-space holes
in multi-dimensional Maxwellian plasmas" by I H Hutchinson.

Copyright I H Hutchinson 2018. No warranty of fitness for any task
whatsoever is given or implied.

Build the code on a linux system by entering for example

$ make dFtdWs

and run it by

$ ./dFtdWs

The executables so formed will generate versions of the figures in the
paper mentioned. 

To change parameters in the executables, edit the main programs.
dFtdWs.f90, fcontko.f90, kpsiarraya.f90, omarray.f90, omsolve.f90 etc.
The parameters declaring array sizes in shiftmode.f90 may need to be
changed to get high resolution results. The main programs have not
been tidied up; so beware.

The graphics library will be generated in ~/src/accis/ which needs to
be creatable and writeable. Also libX11 needs to be linkable via -lX11
which often requires the development version of the xorg package. If
difficulties with graphics are encountered or a non screen displaying
version is desired, see the documentation in ~/src/accis/. Output will
be saved in the form of encapsulated postscript files plot0001.ps etc.
