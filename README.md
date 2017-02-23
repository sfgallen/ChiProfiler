# ChiProfiler

A series of MATLAB functions that utilizes TopoToolbox to analyze river longitudinal profiles

ChiProfiler is a series of Matlab functions that utilize TopoToolbox version 2 (Schwanghart and Scherler, 2014) to conduct river profile analysis using the chi or integral method (Perron and Royden, 2013). All users need to do is download TopoToolbox (https://topotoolbox.wordpress.com/) and it is easy to run ChiProfiler in Matlab.

If you have any questions or comments please contact the author:

Sean F. Gallen

sean.gallen[at]erdw.ethz.ch

## References

Perron, J.T., Royden, L.: An integral approach to bedrock river profile analysis, Earth Surf. Processes Landforms, 38, 570-576, 2013. http://dx.doi.org/10.1002/esp.3302

Schwanghart, W., Scherler, D.: Short Communication: TopoToolbox 2 – MATLAB-based software for topographic analysis and modeling in Earth surface sciences, Earth Surf. Dynam., 2, 1-7, 2014. http://dx.doi.org/10.5194/esurf-2-1-2014

## Requirements

TMatlab 2011b or higher, the MATLAB Image Processing Toolbox, and TopoToolbox.

## Getting started

1) Download TopoToolbox
2) Download ChiProfiler and same in the your TopoToolbox folder
3) Open MATLAB and set paths to TopoToolbox and ChiProfiler with the MATLAB command 'addpath'

Enter following code into the command line:

        addpath C:\path\to\where\you\installed\TopoToolbox-2
        addpath C:\path\to\where\you\installed\TopoToolbox-2\utilities
        addpath C:\path\to\where\you\installed\TopoToolbox-2\topoapp
        addpath C:\path\to\where\you\installed\TopoToolbox-2\DEMdata
        addpath C:\path\to\where\you\installed\TopoToolbox-2\ChiProfiler
