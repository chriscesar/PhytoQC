# PhytoQC
Import and analyse phytoplankton QC data

Data are stored in a mixture of .xls and .xlsx formats.

Formats of individual files is consistent

Approx 250 files need to be imported into a single data set for processing

##Script naming protocol
R script names follow the basic structure
* Prefix 00 - *data import and formatting*
* Prefix 10 - analysis/investigations of **non labswap** data
* Prefix 20 - analysis/investigations of **labswap** data

##To-do list
* Consider generating the report using Rmarkdown to create report as html. This will allow inclusion of 'explorable' 3D NMDS
  * see example here: https://youtu.be/6dyjG8p8U0Q?si=nfBQI41miV3rHm7u&t=924
  