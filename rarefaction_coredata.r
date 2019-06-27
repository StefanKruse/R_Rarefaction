# Script for rarefaction by simple resampling of count data.
# ... application possible for Pollen, aDNA, eDNA, etc. 
# ... 27.06.2019 - v1.0 by Stefan Kruse (stefan.kruse@awi.de)

# ... here I use data from Niemeyer, Bastian; Epp, Laura Saskia; Stoof-Leichsenring, Kathleen Rosmarie; Pestryakova, Ludmila A; Herzschuh, Ulrike (2016): DNA records from surface sediments of lacustrine lakes along a forest tundra - tundra transect at the Taymyr peninsula, northern Siberia. PANGAEA, https://doi.org/10.1594/PANGAEA.860628, In supplement to: Niemeyer, B et al. (submitted): Recording vegetation composition at the Siberian boreal treeline: A comparison between sedimentary DNA and pollen. Molecular Ecology


### prerequisites
# folder set up and data input
	options(stringsAsFactors=FALSE)
	setwd(choose.dir())# select main folder containing folders named exactly "data" and "output"

	specseq=read.csv2("data/DNA_data.CSV")# data table was exported from xlsx here csv2: sep=";", dec=","
	str(specseq)# check data is of required format, count data as integer/double, taxa names as strings
	specseq[is.na(specseq)]=0# make sure additional rows or NA values are converted to zeros

# define columns that contain raw count data
	names(specseq)
	colstart=6
	colend=36
	names(specseq)[colstart:colend]# check
	COLUMNNAMESAREYEARS=FALSE

# species must have unique sample names and need to be merged with the family (for technical reasons)
	SPECIESNAMECOLUMN=1# position of the species assignments column in the data frame
	names(specseq)[SPECIESNAMECOLUMN]="scientific_name"# change species name column to "scientific_name"
	
	FAMILYNAMEPRESENT=FALSE
	FAMILYNAMECOLUMN=2# position of the family assignments column in the data frame
	if(FAMILYNAMEPRESENT)
	{
		names(specseq)[FAMILYNAMECOLUMN]="family_name"# change family name column to "family_name"
	} else
	{
		specseq$family_name="NoFamilyName"# add empty family name column
	}
	
	specseq_final_name=paste(specseq$family_name, make.unique(specseq$scientific_name))# make sure to have individual names for each species/taxa entry





### rarification
# resample loop for each sample/year present in the data table
	(nsampleff=sort(apply(specseq[,c(colstart:colend)], 2, sum))[1])# determine min. read counts for rarefaction, here automatic procedure to find the minimum within the data set
	resamplingnumber=100# set here the number of resamplings, standard==100
	genrare=list()
	yr=NULL# in case of column/sample names are numbers corresponding to depths, an X was given in front of the number on data table import, in the loop this X is deleted and the year is assigned as the name of the entry in the filled loop
	missing=0# this counter is needed to reorganize the list in case samples without any reads are present. These will be skipped (should be not the case otherwise the minimal counts would be zero what will produce in an error so please check before running the script if you have such samples in your data set end exclude them prior analyses).
	
	for(yrcoli in colstart:colend)
	{
	  print(yrcoli)
	  
	  allspec=specseq_final_name[which(specseq[,yrcoli]>0)]
	  allspec_counts=specseq[which(specseq[,yrcoli]>0),yrcoli]
	  
	  if(length(allspec)>0)
	  {
		if(COLUMNNAMESAREYEARS)
		{# convert to number
			yr=c(yr, as.numeric(gsub(",", ".", gsub("X","",names(specseq)[yrcoli]))))
		} else
		{# take the name as it is
			yr=c(yr, names(specseq)[yrcoli])
		}
		
		sampleeffort=list()
		for(nsampleeffi in nsampleff)# would also be possible to process many differrent sampling efforts with the same loop e.g. replace here in the loop call nsampleeff by c(nsampleeff, 100, 1000, ...) but then the script needs to be adapted to take them in the analyses into account
		{
		  repeatsample=list()
		  for(repi in 1:resamplingnumber)
		  {
			repeatsample[[repi]]=sample(allspec,nsampleeffi,replace=TRUE, prob=allspec_counts/sum(allspec_counts))# weighted resampling
		  }
		  sampleeffort[[which(nsampleff==nsampleeffi)]]=repeatsample
		}
		genrare[[yrcoli-missing-(colstart-1)]]=list(allspec,sampleeffort)
	  } else
	  {
		missing=missing+1
	  }
	}
	names(genrare)=yr

	# how to access the data in a resampled dataset (list called genrare)
		# ... show how many taxa in the first newly generated data set were drawn
		length(unique(genrare[[1]][[2]][[1]][[1]]))
		# ... genrare[[1...number of columns]]
		# ... genrare[[1...number of columns]][[1]] == potential unique taxa names in the corresponding sample
		# ... genrare[[1...number of columns]][[2]] == resampled data
		# ... genrare[[1...number of columns]][[2]][[1]] == resampled data of nsampleff, usually only 1 present
		# ... genrare[[1...number of columns]][[2]][[1]][[1:...resamplingnumber]] == resampled data (individual taxa) as determined earlier





# processing of the resampled data set
	# count total species/family number and reads of individual families
	famorig=specseq$family_name
	familylevels=names(rev(sort(table(famorig))))

	totspec=NULL# number of species per sample
	totfam=NULL# number of species per family per sample
	for(li in 1:length(genrare))
	{
	  for(li2 in 1:length(genrare[[li]][[2]]))
	  {
		print(paste0(li," - ",li2))
		
		spectot=NULL
		spectot4fam=NULL
		for(repi in 1:100)
		{
		  pei=unique(genrare[[li]][[2]][[li2]][[repi]])
		  spectot=c(spectot,length(pei))
		  spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split=" ",pei),function(x)return(x[1]))),levels=familylevels)))
		}
		totspec=rbind(totspec, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),Nspecies=spectot))
		totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
	  }
	}
	str(totspec)
	str(totfam)
	
	# simple plots of the processed data
	# modify column names ifnot years that can be coerced to numbers
	if(!COLUMNNAMESAREYEARS)
	{
		totspec$T=factor(totspec$T, levels=names(genrare))# ordered levels to contain original sorting of columns
		totfam$T=factor(totfam$T, levels=names(genrare))# ordered levels to contain original sorting of columns
	}
	png(paste0("output/resampled_totalspecies_Sampleeffort",nsampleff,"_plot.png"), width=480,height=480)
	par(mar=c(8,4,3,1),las=2)
	with(totspec,plot(Nspecies~T, main="number of species per sample"))
	dev.off()

	# save processed data
	write.csv2(totspec, paste0("output/resampled_totalspecies_Sampleeffort",nsampleff,".csv"), row.names=FALSE)	
	write.csv2(totfam, paste0("output/resampled_totalfamilies_Sampleeffort",nsampleff,".csv"), row.names=FALSE)	
	
	

	# aggregate data on species level
	famorig=specseq$scientific_name
	familylevels=names(rev(sort(table(famorig))))
	
	totfam=NULL
	for(li in 1:length(genrare))
	{
	  for(li2 in 1:length(genrare[[li]][[2]]))
	  {
		print(paste0(li," - ",li2))

		spectot4fam=NULL
		for(repi in 1:100)
		{
		  pei=genrare[[li]][[2]][[li2]][[repi]]
		  spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split=" ",pei),function(x)return(paste(x[-1],collapse = " ")))),levels=familylevels)))
		}
		totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
	  }
	}
	str(totfam)

	# modify column names ifnot years that can be coerced to numbers
	if(!COLUMNNAMESAREYEARS)
	{
		totfam$T=factor(totfam$T, levels=names(genrare))# ordered levels to contain original sorting of columns
	}

	# calculate mean values for each species/taxa
	speciesfamiliesdf_totfam=NULL
	pdf(paste0("output/resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated.pdf"))
	par(mar=c(8,4,3,1),las=2)
	for(fami in names(totfam)[3:dim(totfam)[2]])
	{
		plot(totfam[,fami]~totfam$T,col=rainbow(length(names(totfam)[3:dim(totfam)[2]]),s=0.6)[which(names(totfam)[3:dim(totfam)[2]]==fami)],type="n",lwd=2,main=fami,ylab="Sequence counts",xlab="Depth (m)")
		  
		mn=aggregate(totfam[,fami],list(as.numeric(totfam$T)),mean)
		ti=mn$Group.1
		mn=mn$x
		sd=aggregate(totfam[,fami],list(as.numeric(totfam$T)),sd)$x*1.96
		
		polygon(y=c(mn+sd,rev(mn-sd)), x=c(ti,rev(ti)), col="gray", border=NA)
		lines(mn~ti,col="steelblue2",lwd=3)
		  
		if(!COLUMNNAMESAREYEARS)
		{
			ti=levels(totfam$T)
		}
		speciesfamiliesdf_totfam=rbind(speciesfamiliesdf_totfam, data.frame(Species=fami,TBP=ti,Mean=mn,CI95=sd))
	}
	dev.off()

	str(speciesfamiliesdf_totfam)
	
	# save processed data
	write.csv2(speciesfamiliesdf_totfam, paste0("output/resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated.csv"), row.names=FALSE)




### post processing	
# reformat for pca input ... rows are samples, cols are species/taxa
	ordidf=NULL
	for(fi in levels(factor(speciesfamiliesdf_totfam$Species)))
	{
	  numberi=speciesfamiliesdf_totfam[speciesfamiliesdf_totfam$Species==fi, ]$Mean
	  ordidf=rbind(ordidf, data.frame(Spec=fi,t(numberi)))
	}
	names(ordidf)[2:dim(ordidf)[2]]=speciesfamiliesdf_totfam[speciesfamiliesdf_totfam$Species==fi, ]$TBP
	str(ordidf)
	
	# remove species/taxa column
	row.names(ordidf)=ordidf$Spec
	ordidf=ordidf[,-1]
	str(ordidf)
	head(ordidf)
	
	# exclude samples or species that have no records
	rowsumsnotzero=which(apply(ordidf,1,sum)>0)
	colsumsnotzero=which(apply(ordidf,2,sum)>0)
	ordidf=ordidf[rowsumsnotzero,colsumsnotzero]

	# export data
	write.csv2(t(ordidf), paste0("output/resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated_pcainput.csv"))

# comparison of original and resampled data set
	# prepare the data
	ordiorigdf=specseq[colstart:colend]
	row.names(ordiorigdf)=make.unique(specseq$scientific_name)
	rowsumsnotzero=which(apply(ordiorigdf,1,sum)>0)
	colsumsnotzero=which(apply(ordiorigdf,2,sum)>0)
	ordiorigdf=ordiorigdf[rowsumsnotzero,colsumsnotzero]

	# species count (counts >= 1)
	png(paste0("output/resampled_speciesnumber_Sampleeffort",nsampleff,"_aggregated_comparisonplot.png"), width=480,height=480)
		par(mar=c(8,4,3,1),las=2)
		barplot(apply(ordiorigdf,2,function(x)length(which(x>=1))), col="tomato", border=FALSE)
		barplot(apply(ordidf,2,function(x)length(which(x>=1))), add=TRUE, col="skyblue", border=FALSE)
		legend("topright", c("original","resampled"), col=c("tomato","skyblue"), pch=15, pt.cex=1.5, title="Species number >= 1 count")
	dev.off()
	
	# ordination
	pca_original=prcomp(sqrt(sqrt(t(ordiorigdf))))
	pca_resampled=prcomp(sqrt(sqrt(t(ordidf))))
	
	png(paste0("output/resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated_pca_comparisonplot.png"), width=960,height=480)
		par(mfrow=c(1,2))
		biplot(pca_original, main="original data")
		biplot(pca_resampled, main="rarefied data")
	dev.off()

