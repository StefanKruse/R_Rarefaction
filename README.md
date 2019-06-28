# R_Rarefaction 1.0
R code for resampling and thus normalizing of count data to the minimum number of counts across a set of samples (e.g. sedaDNA sequence/pollen taxa counts per sample along a sediment core).

![alt text](https://github.com/StefanKruse/R_Rarefaction/blob/master/output/resampled_speciesnumber_Sampleeffort1710_aggregated_comparisonplot.png)
Figure 1. Rafefaction is needed to bring raw count data to a common ground for downstream analyses. In the figure you see the uneven reduction of total number of species present in the samples when resampled to the minimum counts across all considered samples.

### Version history:
- 27.06.2019 initial upload

### Authors:
- Stefan Kruse - stefan.kruse@awi.de

## Containing files:
1. master file including the script and explanations "rarefaction_coredata.r"
2. folder "data" for data input
2.1. contains data from Niemeyer, Bastian; Epp, Laura Saskia; Stoof-Leichsenring, Kathleen Rosmarie; Pestryakova, Ludmila A; Herzschuh, Ulrike (2016): DNA records from surface sediments of lacustrine lakes along a forest tundra - tundra transect at the Taymyr peninsula, northern Siberia. PANGAEA, https://doi.org/10.1594/PANGAEA.860628, In supplement to: Niemeyer, B et al. (submitted): Recording vegetation composition at the Siberian boreal treeline: A comparison between sedimentary DNA and pollen. Molecular Ecology
3. folder "output" for output
3.1. contains output that is produced by the running the code
