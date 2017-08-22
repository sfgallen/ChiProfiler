# ChiProfiler

A series of MATLAB functions that utilizes TopoToolbox to analyze river longitudinal profiles

ChiProfiler is a series of Matlab functions that utilize TopoToolbox version 2 (Schwanghart and Scherler, 2014) to conduct river profile analysis using the chi or integral method (Perron and Royden, 2013). ChiProfiler also allows users to generate maps of river network metrics, such as the normalized steepness index (ksn) and the integral quantity chi (Wobus et al., 2006; Willett et al., 2014). ChiProfiler was developed by Sean Gallen and is shortly described and applied in Gallen and Wegmann (2017). Please cite Gallen and Wegmann (2017) if you use these codes for scientific research. All users need to do is download TopoToolbox (https://topotoolbox.wordpress.com/) and it is easy to run ChiProfiler in Matlab.

If you have any questions or comments please contact the author:

Sean F. Gallen

sean.gallen[at]erdw.ethz.ch

## References

Gallen, S.F., Wegmann, K.W.: River profile response to normal fault growth and linkage: An example from the Hellenic forearc of south-central Crete, Greece, Earth Surf. Dynam., 2017, http://www.earth-surf-dynam.net/5/161/2017/.

Perron, J.T., Royden, L.: An integral approach to bedrock river profile analysis, Earth Surf. Processes Landforms, 38, 570-576, 2013. http://dx.doi.org/10.1002/esp.3302

Schwanghart, W., Scherler, D.: Short Communication: TopoToolbox 2 – MATLAB-based software for topographic analysis and modeling in Earth surface sciences, Earth Surf. Dynam., 2, 1-7, 2014. http://dx.doi.org/10.5194/esurf-2-1-2014

Willett, S.D., McCoy, S.W., Perron, J.T., Goren, L., Chen, C.-Y.: Dynamic Reorganization of River Basins, Science, 343, 2014. http://dx.doi.org/10.1126/science.1248765

Wobus, C., Whipple, K.X., Kirby, E., Snyder, N., Johnson, J., Spyropolou, K., Crosby, B., Sheehan, D.: Tectonics from topography: Procedures, promise, and pitfalls, Geological Society of America Special Papers, 398, 55-74, 2006. http://dx.doi.org/10.1130/2006.2398(04)


## Requirements

Matlab 2014b or higher, the MATLAB Image Processing Toolbox, and TopoToolbox.

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
        
ChiProfiler also comes with a user guide that can help you get started.

## Notes

Some people have reported errors associated with the river profile analysis GUI for specific versions of Matlab (2015a, I think). This apparently has to do with a bug in the Matlab command 'gtext' that is called to put information on the ksn regressions on the figure. If you run into this problem, try updating the version of Matlab that you are using.
