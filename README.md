### AWEToolchain 
This repository consists of the Main codes of the Toolchain from downloading MERRA wind data, processing the data for the optimization model and estimating the optimal power output of one power cycle of an AWE system

#PowerOptimization 
Contains
- Main_ClusterRun: RUN this code as main, here you set the city, and size of the kite. It runs through the different clusters of the wind data in order to speed up the optimization
-  _OCPmodel: Optimal Control Problem (OCP) formulation including the model of the AWE system
- init_Drag_wind: is called by the Main in order to run the homotopy through the different OCP
- Collocation: collocation code for discretizing the OCP
- parameters: includes the parameters for the aerodynamic model and the  different sizes of the AWE system
- _aero_poly: aerodynamic model of the kite including the wind data as polynomial
- _fcts: includes function as creating the Lagrange polynomial of the downloaded wind data
- sol_ClusterRun: plots the results 
- 
#MERRA_Download
Contains:
- merra_download: RUN this code for download. Choose the cite via latitude and longitude coordinates,  the time span and the variables you are interested in as well as the required altitudes of the chosen variables. 
- _downloadwind: functions for the whole download: NOTE!! You need your own access to the NASA big data site in order to run the automated download. See more info on: https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget

- clustering: generates from  the east / north wind a main wind speed and a deviation so that this can be aligned with the x direction of the kite model; further it clusters the wind speeds to 7 different groups, such that the optimization can run in a kind of homotopy through similar wind shears.
- merra2npy: reads the *.nc4 files including the wind data and saves into an *.npy, which is better readable with python.
