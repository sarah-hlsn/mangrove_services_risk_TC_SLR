# Mangroves and their services are at risk from tropical cyclones and sea level rise under climate change

These scripts reproduce the main results of the manuscript:

 **Mangroves and their services are at risk from tropical cyclones and sea level rise under climate change.**
 
Sarah Hülsen (1,2), Laure E. Dee (3,4), Chahan, M. Kropf (1,2), Simona Meiler (1,2,5), and David N. Bresch (1,2)

Publication status: under review

(1) Institute for Environmental Decisions, ETH Zurich, Switzerland

(2) Federal Office of Meteorology and Climatology MeteoSwiss, Switzerland

(3) Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, USA

(4) Nature-Based Solutions Initiative, University of Colorado – Boulder, Boulder, CO, USA

(5) Department of Civil and Environmental Engineering, Stanford University, CA, USA

Contact: Sarah Hülsen ([sarah.huelsen@usys.ethz.ch](sarah.huelsen@usys.ethz.ch))

## Content:
`extract_mangrove_centroids.ipynb`: Script to calculate mangrove centroids and area from GMW global mangrove extents for 2020
`extract_priority_areas.ipynb`: Script to merge mangrove conservation priority areas (including ecosystem service data) from Dabalà et al. and merge with mangrove centroids
`extract_rslr.ipynb`: Script to extract RSLR data and merge with global mangrove centroids
`calc_frequency_change.py`: Script to calculate TC frequency changes and merge with mangrove centroids

## Requirements
Requires:
- Python 3.9+ environment (best to use conda for CLIMADA repository)
- *CLIMADA* repository version 3.3.3+: https://github.com/CLIMADA-project/climada_python

## ETH cluster
Computationally demanding calculations were run on the [Euler cluster of ETH Zurich](https://scicomp.ethz.ch/wiki/Euler).

Documentation for CLIMADA is available on Read the Docs:
* [online (recommended)](https://climada-python.readthedocs.io/en/stable/)
* [PDF file](https://buildmedia.readthedocs.org/media/pdf/climada-python/stable/climada-python.pdf)