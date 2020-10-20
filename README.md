# R_Rarefaction 1.0
R code for resampling and thus normalizing of count data to the minimum number of counts across a set of samples (e.g. sedaDNA sequence/pollen taxa counts per sample along a sediment core).

![Barplot comparing original to resampled data](https://github.com/StefanKruse/R_Rarefaction/blob/master/output/replace_resampled_speciesnumber_Sampleeffort1710_aggregated_comparisonplot.png)

Figure 1. Rarefaction is needed to bring original count data to a common ground for downstream analyses. In the figure you see the uneven reduction of total number of species present in the samples when resampled to the minimum counts across all considered samples.

### Version history
- 20.10.2020 updated with Birks & Line (1992) sampling, minor fixes
- 27.06.2019 initial upload

### Authors
- Stefan Kruse - stefan.kruse@awi.de

## Containing files
1. master file including the script and explanations "rarefaction_coredata.r"
1. folder "data" for data input
	1. contains data from Niemeyer et al. (2016)
1. folder "output" for output
	1. contains output that is produced by the running the code
	1. filename contains first sampling method, and the minimum number that to which the data is resampled

## References
Birks, H. J. B., & Line, J. M. (1992). The use of Rarefaction Analysis for Estimating Palynological Richness from Quaternary Pollen-Analytical Data. The Holocene, 2(1), 1–10. doi:10.1177/095968369200200101
Niemeyer, Bastian; Epp, Laura Saskia; Stoof-Leichsenring, Kathleen Rosmarie; Pestryakova, Ludmila A; Herzschuh, Ulrike (2016): DNA records from surface sediments of lacustrine lakes along a forest tundra - tundra transect at the Taymyr peninsula, northern Siberia. PANGAEA, https://doi.org/10.1594/PANGAEA.860628 - In supplement to: Niemeyer, Bastian; Epp, Laura Saskia; Stoof-Leichsenring, Kathleen Rosmarie; Pestryakova, Luidmila A; Herzschuh, Ulrike (2017): A comparison of sedimentary DNA and pollen from lake sediments in recording vegetation composition at the Siberian treeline. Molecular Ecology Resources, 17(6), e46-e62, https://doi.org/10.1111/1755-0998.12689
