## Use Case: Monitoring Tropical Forest Recovery Capacity Using **RADAR** satellite data

__Contact persons:__ 
* Milutin Milenković, milutin.milenkovic(at)wur.nl, Wageningen University & Research (WUR), The Netherlands
* Raymond Oonk, raymond.oonk(at)surf.nl, SURF (Samenwerkende Universitaire Reken Faciliteiten), The Netherlands
## Scope
This use case utilizes time series of Sentinel-1 (S1) backscatter intensity (Sigmma0) images to monitor tropical forest recovery capacity over the Amazon basin. Sentinel-1 is a C-band, synthetic aperture RADAR imaging mission consisting of two satellites (S1A and S1B) operating together since 2016 and acquiring images globally, with two polarizations channels (VV and VH), a spatial resolution of 20 x 5 m<sup>2</sup>  (in the interferometric wide-swath, IW, mode), and a revisiting time of maximum 6 days. A C-SCALE data provider, Earth Observation Data Center (EODC) GmbH, has provided the S1 backscatter intensity datacube over the Amazon basin for this use case (Wagner et al. 2021). For each location (pixel) of the datacube, a 5-year-long S1 backscatter intensity time series (2017–2021) has been analyzed to detect the magnitude (_Mag_) and recovery time (_t<sub>total</sub>_) of the signal disturbance.  The analysis was based on an algorithm that combines the original and smoothed time series and the statistical properties of the first- and last-year or the S1 data to detect the signal disturbances. 

![Alt text](figures/ts_figure.png?raw=true "XXXXX")

This use case aims to (a) derive Amazon-wide disturbance magnitude and recovery time maps from the Sentinel-1 image time series and (b) analyze the functional relation between derived magnitude and disturbance time.  
## Implementation

The use case has been tested and deployed on the SURF computing infrastructure, i.e. the Spider cluster, and the analysis-ready data (the Amazon-wide S1 datacube) provided by EODC. Access to the datacube has been ensured through a license agreement between WUR and EODC. This GitHub repository contains an algorithm (Python code) developed by WUR to detect the disturbance magnitude and recovery time of the S1 backscatter intensity time series. 

The use case workflow takes the S1 datacube as the input and applies the disturbance detection algorithm, providing the disturbance magnitude and recovery time per each pixel, i.e. x, y location. As the input datacube is a file-based image collection split in the 300 x 300 km<sup>2</sup>  Equi7Grid tiles (Bauer-Marschallinger et al. 2014), the output raster maps with disturbance magnitude values (in dB) and the recovery periods (in days) are also split in the 300 x 300 km<sup>2</sup>  tiles.

### Examples

The folder Notebooks contains examples on how to: access Sentinel-1 data, make a data cube, plot time series and maps, and analyze the data mostly based on yeoda and xarray packages as well as some custom-defined helper functions. Below is the list of Python notebooks with short explanations.

* plot_ts_yeoda.ipynb - Plot time seres from a specifed location using the yeoda modole to query the EODC S1 datacube
* test_ts_yeoda_processing.ipynb - Take a S1 time series from a specified location and calculate the disturbance magnitude and disturbance time
* tempReducer_yeoda.ipynb - Reduce the S1 data cube along the temporal axis to a raster
* checkFileSizes.ipynb - Check the file size of each image in the data folders and write those image files with an reading error

    * annualStat_Equi7tile.ipynb - Calculate annual statistic per time seres (for each pixel) within a single Equi7grid tile considering all images
    * annualStat_Equi7tile_AscOnly.ipynb - Calculate annual statistic per time seres (for each pixel) within a single Equi7grid tile considering only ascending images
    * annualStat_Equi7tile_DescOnly.ipynb - Calculate annual statistic per time seres (for each pixel) within a single Equi7grid tile considering only descending images
    * process_Equi7tile.ipynb - Apply the time series analysis per each x, y location (pixel) within a single Equi7grid tile
    * debug_process_Equi7tile.ipynb - Simple debug setup for the process_Equi7tile function
    * debug_Equi7tile_DescOnly.ipynb - Debug setup for the processing of descending images
    * debug_row_col_ts.ipynb - debug setup for a single chank within a Equi7Tile
    * output_and_forestAge.ipynb - Analysis of the output (S1 time series features) with forest age map
    * outputChanks_and_forestAge_inSituData.ipynb - Analysis of the output (S1 time series features) with forest age map exstended to in-situ data analysis and recovery time statistics
    * spatialUnit_analysis.ipynb - Analyzing optimal spatial unit for the space-for-time analysis
    * read_and_inspect_output.ipynb - Analyze single chank and derive annual statistics maps
    * outputChunks_explore.ipynb - Analyze single chanks and deal with Equi7Grid
    * outputChunks_geoWombat.ipynb - Test reading chanks with geoWombat
    * mosaic_output_with_stac.ipynb - Test reading chanks with stack 
    * s1_age_and_IntactForest.ipynb -  Analysis of the output (S1 time series features) with forest age map using the old growth (intact) forest as the refernce
    * outExplore_plotDistance.ipynb -  Analysis of the output (S1 time series features) as the feature space difference from the old growth vector
    * outExplore_mapFeatDistance.ipynb - Plot a map with magnitute detected in disturbed S1 time seres
    
    
    
    
    
    



## References

Wagner, W.; Bauer-Marschallinger, B.; Navacchi, C.; Reuß, F.; Cao, S.; Reimer, C.; Schramm, M.; Briese, C. A Sentinel-1 Backscatter Datacube for Global Land Monitoring Applications. Remote Sens. 2021, 13, 4622. https://doi.org/10.3390/rs13224622

Bauer-Marschallinger, B.; Sabel D.; agner W. Optimisation of global grids for high-resolution remote sensing data, Computers & Geosciences, 2014, Volume 72,  Pages 84-93, ISSN 0098-3004. https://doi.org/10.1016/j.cageo.2014.07.005.
