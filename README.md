# Matlab-scripts-for-Calcium
Scripts for Analyse of Calcium oscilations
Analyse_Calcium_Data.m For determination of peaks of all waves and turn them into binary data. In detail, algorithm determined base line of wave by 9 fixed points and smooth with “Loess” method with 30% point windows around each point. Point to point difference between base line and raw data was determined by formula Delta F = calcium traces – base line
AnalyseDuree_ONOFF_Raster.m To calculate total amount of frames in seconds where wave in “On” mode defined as “Total time in On”.
Analyse_FFT.m Fast Fourier Transformation (FFT) algorithm with computation of power spectrum (Mourao et al., 2014a; Mourao et al., 2014b; Uhlen, 2004).
Script_multiGraph.m to obtain figures of individual wave traces 
Analyse_Ensemble_Calcium_Data.m Analysis of coactive networks described by (Miller et al., 2014). Four inputs file are required: a zip file from ImageJ, for Regions of Interest (ROI), the csv of raster ON-OFF, one projected tiff for the tissue and csv for coordinate value of the ROI center. 
Analyse_Calcium_Data_Network.m For Analysis of functional networks was described by (Rubinov and Sporns, 2010) Three inputs file are used: the csv of F/F0 waves, one projected tiff for the tissue and csv for coordinate value of the ROI center. 
All change in existing code and generation of new code done by Pierre Fontanaud
