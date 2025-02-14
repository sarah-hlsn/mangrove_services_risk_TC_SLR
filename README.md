# Mangroves and their services are at risk from tropical cyclones and sea level rise under climate change

These scripts reproduce the main results of the manuscript:

 **Mangroves and their services are at risk from tropical cyclones and sea level rise under climate change.**
 
Sarah Hülsen (1,2), Laure E. Dee (3,4), Chahan, M. Kropf (1,2), Simona Meiler (1,2,5), and David N. Bresch (1,2)

Publication status: under review ([preprint](https://doi.org/10.21203/rs.3.rs-5346064/v1))

(1) Institute for Environmental Decisions, ETH Zurich, Switzerland

(2) Federal Office of Meteorology and Climatology MeteoSwiss, Switzerland

(3) Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, USA

(4) Nature-Based Solutions Initiative, University of Colorado – Boulder, Boulder, CO, USA

(5) Department of Civil and Environmental Engineering, Stanford University, CA, USA

Contact: [Sarah Hülsen](https://wcr.ethz.ch/people/person-detail.sarah-huelsen.html) ([sarah.huelsen@usys.ethz.ch](sarah.huelsen@usys.ethz.ch))

## Content:
This directory contains all the necessary code to replicate the analysis of the above mentioned manuscript. The files should be run in the order listed here:

`01_mangrove_data_prep.ipynb`: Notebook to calculate mangrove centroids and area from GMW global mangrove extents for 2020, extract RSLR data and merge with global mangrove centroids.

`02_calc_frequency_change.py`: Script to calculate TC frequency changes and merge with mangrove centroids.

`03_risk_index_analysis.ipynb`: Notebook to calculate combined RSLR and TC risk index, merge index with mangrove conservation priority areas (including ecosystem service data) from Dabalà et al., and calculate all further results. Includes tables shown in manuscript.

`04_figures.ipynb`: Notebook containing code for all main figures of the manuscript. Data to replicate the figures shown is available from Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14869916.svg)](https://doi.org/10.5281/zenodo.14869916)

## Data sources:
Mangroves: Bunting, P., Rosenqvist, A., Lucas, R., Rebelo, L.-M., Hilarides, L., Thomas, N., Hardy, A., Itoh, T., Shimada, M., & Finlayson, C. (2018). The Global Mangrove Watch—A New 2010 Global Baseline of Mangrove Extent. Remote Sensing, 10(10), 1669. https://doi.org/10.3390/rs10101669

Mangrove conservation priority areas: Dabalà, A., Dahdouh-Guebas, F., Dunn, D. C., Everett, J. D., Lovelock, C. E., Hanson, J. O., Buenafe, K. C. V., Neubert, S., & Richardson, A. J. (2023). Priority areas to protect mangroves and maximise ecosystem services. Nature Communications, 14(1), Article 1. https://doi.org/10.1038/s41467-023-41333-3

Sea level rise data: Garner, G. G., Hermans, T., Kopp, R. E., Slangen, A. B. A., Edwards, T. L., Levermann, A., Nowicki, S., Palmer, M. D., Smith, C., Fox-Kemper, B., Hewitt, H. T., Xiao, C., Aðalgeirsdóttir, G., Drijfhout, S. S., Golledge, N. R., Hemer, M., Krinner, G., Mix, A., Notz, D., ... Pearson, B. (2021). IPCC AR6 Sea Level Projections (Version 20210809). Zenodo. https://doi.org/10.5281/zenodo.5914710

The tropical cyclone tracks were kindly provided by Kerry Emanuel, and are based on the following sources: 
Emanuel, K., Ravela, S., Vivant, E., & Risi, C. (2006). A Statistical Deterministic Approach to Hurricane Risk Assessment. Bulletin of the American Meteorological Society, 87(3), 299– 314. https://doi.org/10.1175/BAMS-87-3-299

Emanuel, K. (2008). The Hurricane-Climate Connection. Bulletin of the American Meteorological Society, 89, ES10–ES20. https://doi.org/DOI:10.1175/BAMS-89-5-Emanuel

## Requirements
Requires:
- Python 3.9+ environment (best to use conda for CLIMADA repository)
- *CLIMADA* repository version 3.3.3+: https://github.com/CLIMADA-project/climada_python

## ETH cluster
Computationally demanding calculations were run on the [Euler cluster of ETH Zurich](https://scicomp.ethz.ch/wiki/Euler).

Documentation for CLIMADA is available on Read the Docs:
* [online (recommended)](https://climada-python.readthedocs.io/en/stable/)
* [PDF file](https://buildmedia.readthedocs.org/media/pdf/climada-python/stable/climada-python.pdf)

## How to cite this source:
Sarah Hülsen. (2025). sarah-hlsn/mangrove_services_risk_TC_SLR: Code for "Mangroves and their services are at risk from tropical cyclones and sea level rise under climate change" (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.14870133 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14870133.svg)](https://doi.org/10.5281/zenodo.14870133)
