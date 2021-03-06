# SapsuckRz
Repo for Aphid Suction trap project.

We aim to investigate how climate and landscape (including crop cover) affect aphid emergence and peak abundace across a large extent of the United States over a period of 11 years (2005-2015). We hope to understand potential drivers of pest and non-pest aphid species abundance and how this might change due to antropogenic climate change.

## Data we are using

We will be combining aphid suction trap data from 46 sites across the U.S. with climate and landscape data

Aphid abuance data set available from [KBS LTER database](http://lter.kbs.msu.edu/datatables/122)

weather data pulled from the [prism_weather_data_pull.R](https://github.com/ReproducibleQM/sapsuckRz/blob/master/prism_weather_data_pull.R) scirpt can be found on figureshare: [here](https://figshare.com/articles/weather_data_csv/5012747)
The analysis script pulls this data into R from figshare.

landcover data pulled from [cdl_landcover_data_pull.R](https://github.com/ReproducibleQM/sapsuckRz/blob/master/cdl_landcover_data_pull.R) via the USDA Cropland Data Layer. The intermediate data product is stored directly on this repository: [here](https://github.com/ReproducibleQM/sapsuckRz/blob/master/Data/cdl_landcover_data.csv).

## EDA- Exploratory Data Analysis
The folder with EDA was generated during the process of working through the data and potential projects/vinettes we hoped to move forward with. As these are exploratory, the code isn't super clean; however, we thought it useful to highlight the scientific process and any code we utilized prior to the final analysis.

## Follow our progress
[google drive folder](https://drive.google.com/drive/folders/0B7EmIF4p0bakV01yanpHSHFYYzA?usp=sharing) for non-code analysis

## More about the project
This project is a group project the [Reproducible Quantitative methods course](https://cbahlai.github.io/rqm-template/) offered at MSU by [Christie Bahlai](https://sites.google.com/site/cbahlai/).
