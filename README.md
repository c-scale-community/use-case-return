# Use Case:
_Title_: Monitoring Tropical Forest Recovery Capacity Using **RADAR** satellite data

_Contact persons_: 
* Milutin Milenković, milutin.milenkovic(at)wur.nl, Wageningen University & Research (WUR), The Netherlands
* Raymond Oonk, raymond.oonk(at)surf.nl, SURF (Samenwerkende Universitaire Reken Faciliteiten), The Netherlands
## Scope
This use case utilizes time series of Sentinel-1 (S1) backscatter intensity (Sigmma0) images to monitor tropical forest recovery capacity over the Amazon basin. Sentinel-1 is a C-band, synthetic aperture RADAR imaging mission consisting of two satellites (S1A and S1B) operating together since 2016 and acquiring images globally, with two polarizations channels (VV and VH), a spatial resolution of 20 x 5 m<sup>2</sup>  (in the interferometric wide-swath, IW, mode), and a revisiting time of maximum 6 days. A C-SCALE data provider, Earth Observation Data Center (EODC) GmbH, has provided the S1 backscatter intensity datacube over the Amazon basin for this use case (Wagner et al. 2021). For each location (pixel) of the datacube, a 5-year-long S1 backscatter intensity time series (2017–2021) has been analyzed to detect the magnitude and recovery time of the signal disturbance.  The analysis was based on an algorithm that combines the original and smoothed time series and the statistical properties of the first- and last-year or the S1 data to detect the signal disturbances. 

![Alt text](figures/ts_figure.png?raw=true "XXXXX")

This use case aims to (a) derive Amazon-wide disturbance magnitude and recovery time maps from the Sentinel-1 image time series and (b) analyze the functional relation between derived magnitude and disturbance time.  
## Implementation

The use case has been tested and deployed on the SURF computing infrastructure, i.e. the Spider cluster, and the analysis-ready data (the Amazon-wide S1 datacube) provided by EODC. Access to the datacube has been ensured through a license agreement between WUR and EODC. This GitHub repository contains an algorithm (Python code) developed by WUR to detect the disturbance magnitude and recovery time of the S1 backscatter intensity time series. 

The use case workflow takes the S1 datacube as the input and applies the disturbance detection algorithm, providing the disturbance magnitude and recovery time per each pixel, i.e. x, y location. As the input datacube is a file-based image collection split in the 300 x 300 km<sup>2</sup>  Equi7Grid tiles (Bauer-Marschallinger et al. 2014), the output raster maps with disturbance magnitude values (in dB) and the recovery periods (in days) are also split in the 300 x 300 km<sup>2</sup>  tiles.

## References

Wagner, W.; Bauer-Marschallinger, B.; Navacchi, C.; Reuß, F.; Cao, S.; Reimer, C.; Schramm, M.; Briese, C. A Sentinel-1 Backscatter Datacube for Global Land Monitoring Applications. Remote Sens. 2021, 13, 4622. https://doi.org/10.3390/rs13224622

Bernhard Bauer-Marschallinger, Daniel Sabel, Wolfgang Wagner, Optimisation of global grids for high-resolution remote sensing data, Computers & Geosciences, Volume 72, 2014, Pages 84-93, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2014.07.005.
