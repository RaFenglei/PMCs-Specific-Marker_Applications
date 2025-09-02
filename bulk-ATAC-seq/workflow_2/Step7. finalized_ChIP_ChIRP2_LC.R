###Finalized script for submission
###Lauren C. March 27 2023


library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
options(scipen=999)

#See script scripttoget50bpbintable.sbatch for how to make the txt file with the 50bp bins.
###see file bedintersectcommands_finalized.txt for the commands for bedtools intersect .

######Part 0: getting avg vers of the 50 bp bin files that are outputted from Multi Big Wig Summary.###
# managing the need to get average of Rep1-2:
# 1.	BamCoverage on all individual files with CPM norm. Output as bigwig. 
# 2.	MultibigwigSummary, with all individual files (so rep1 and 2 separately), with bins = 1000, outputting the raw table with --outRawCounts 
# 3.	Then I can put the resulting table from MBWS into R, and manually get mean values between reps 1-2 and make the plots with either individual replicates or as mean values


#import MBWS and then do plot of all indiv and take mean of the rep1-2s, and directly plot data as samples and as chr
option1table<-read.table("Average_CPM_bins50bp_MBWSummary_Option1.txt", header = F, sep = "\t")
head(option1table)
colnames(option1table)<-c("chr", "start", "end", "ctrinput", "ctr_rep1", "ctr_rep2", "l1kdinput", "l1kd_rep1", "l1kd_rep2")
#drop extra chr
option1table<-option1table[!(grepl("alt", option1table$chr)),]; option1table<-option1table[!(grepl("rand", option1table$chr)),]
option1table<-option1table[!(grepl("unk", option1table$chr)),]; option1table<-option1table[!(grepl("chrM", option1table$chr)),]
option1table<-option1table[!(grepl("chrY", option1table$chr)),];

# handling of 0s disrupting the boxplot display when u compare indiv vs means
#use new method for means where we set 0s to NA, then when u take mean it will just be the remaining value.
option1tablewithmeansNAs<-option1table
option1tablewithmeansNAs<-option1tablewithmeansNAs %>% mutate(across(c("ctr_rep1", "ctr_rep2", "l1kd_rep1", "l1kd_rep2"), ~ na_if(., 0)))
head(option1tablewithmeansNAs)
option1tablewithmeansNAs$l1kdr1and2mean<-rowMeans(option1tablewithmeansNAs[,c("l1kd_rep1","l1kd_rep2")], na.rm = T)
option1tablewithmeansNAs$ctr1and2mean<-rowMeans(option1tablewithmeansNAs[,c("ctr_rep1", "ctr_rep2")], na.rm = T)
option1tablewithmeansNAs<-option1tablewithmeansNAs[,c("chr", "start", "end", "ctrinput", "l1kdinput", "l1kdr1and2mean", "ctr1and2mean"),]

write.table(option1tablewithmeansNAs, "option1tablewithmeansNAs.txt", quote = F, sep = "\t", row.names = F)
read.table("option1tablewithmeansNAs.txt", sep = "\t", header = T)
#also manually made a version of this without header, for bedtools use.

#also record table with indiv reps with the 0 replaced by NAs
option1tableNAs<-option1table
option1tableNAs<-option1tableNAs %>% mutate(across(c("ctr_rep1", "ctr_rep2", "l1kd_rep1", "l1kd_rep2"), ~ na_if(., 0)))
head(option1tableNAs)
option1tableNAs<-option1tableNAs[,c("chr", "start", "end", "ctrinput", "l1kdinput", "l1kdr1and2mean", "ctr1and2mean"),]
write.table(option1tableNAs, "option1tablewithreps_NAs.txt", quote = F, sep = "\t", row.names = F)
write.table(option1tableNAs, "option1tablewithreps_NAs_noHeader.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#also manually made a version of this without header, for bedtools use.

###end of making 50bp cpm tables.####




####Part 1: Select strict peaks in L1KD and Ctr-ChIP conditions.######
##Need to assess replicability between the L1KD and Ctr-ChIP replicates.
##get bed files for ctr chip and l1kd chip reps1-2
#and then do each comparison from both sides. so rep1 int rep2, keep originals in A; and then rep2 int rep1, keep originals in new A.

 
###need bed intersect commands for these few files.

#now load in the filtered beds:
ctrchiprep1<-read.table("ctr-ChIP-rep1_BFpeaksinrep2.bed", header=F, sep = "\t")
head(ctrchiprep1)
#the 4th col is the transformed -log10(qvalue) 
colnames(ctrchiprep1)<-c("chr", "start", "end", "qvaltrans", "summit", "strand")

ctrchiprep2<-read.table("ctr-ChIP-rep2_BFpeaksinrep1.bed", header = F, sep = "\t" )
colnames(ctrchiprep2)<-c("chr", "start", "end", "qvaltrans", "summit", "strand")
l1kdrep1<-read.table("L1KD-ChIP-rep1_BFpeaksinrep2.bed", header = F, sep = "\t" )
colnames(l1kdrep1)<-c("chr", "start", "end", "qvaltrans", "summit", "strand")
l1kdrep2<-read.table("L1KD-ChIP-rep2_BFpeaksinrep1.bed", header = F, sep = "\t" )
colnames(l1kdrep2)<-c("chr", "start", "end", "qvaltrans", "summit", "strand")


##can merge these, then get min start of pair and max start of pair, to create a 'region'.
l1kdboth<-cbind(l1kdrep1, l1kdrep2)
colnames(l1kdboth)<-c("chr",  "rep1start", "rep1end", "rep1qvaltrans", "rep1summit", "strand", "chr2",  "rep2start", "rep2end", "rep2qvaltrans", "rep2summit", "rep2strand")
ctrchipboth<-cbind(ctrchiprep1, ctrchiprep2)
colnames(ctrchipboth)<-c("chr",  "rep1start", "rep1end", "rep1qvaltrans", "rep1summit", "strand", "chr2",  "rep2start", "rep2end", "rep2qvaltrans", "rep2summit", "rep2strand")

l1kdboth<-l1kdboth[,c("chr",  "rep1start", "rep1end", "rep1qvaltrans", "rep1summit", "strand", "rep2start", "rep2end", "rep2qvaltrans", "rep2summit")]
ctrchipboth<-ctrchipboth[,c("chr",  "rep1start", "rep1end", "rep1qvaltrans", "rep1summit", "strand", "rep2start", "rep2end", "rep2qvaltrans", "rep2summit")]

l1kdboth$minstart<-apply(l1kdboth[,c("rep1start","rep2start")], 1, FUN = min)
l1kdboth$maxend<-apply(l1kdboth[,c("rep1end","rep2end")], 1, FUN = max)
ctrchipboth$minstart<-apply(ctrchipboth[,c("rep1start","rep2start")], 1, FUN = min)
ctrchipboth$maxend<-apply(ctrchipboth[,c("rep1end","rep2end")], 1, FUN = max)
head(l1kdboth)

##next load in the 50bp bin processed multi bigwig summary file.
cpmtable<-read.table("Average_CPM_bins50bp_MBWSummary_Option1.txt", header = F, sep = "\t")
colnames(cpmtable)<-c("chr", "start", "end", "ctrinput", "ctr_rep1", "ctr_rep2", "l1kdinput", "l1kd_rep1", "l1kd_rep2")
#drop extra chromosomes.
cpmtable<-cpmtable[!(grepl("alt", cpmtable$chr)),]; cpmtable<-cpmtable[!(grepl("rand", cpmtable$chr)),]
cpmtable<-cpmtable[!(grepl("unk", cpmtable$chr)),]; cpmtable<-cpmtable[!(grepl("chrM", cpmtable$chr)),]
cpmtable<-cpmtable[!(grepl("chrY", cpmtable$chr)),];

l1kdboth<-l1kdboth[!(grepl("alt", l1kdboth$chr)),]; l1kdboth<-l1kdboth[!(grepl("rand", l1kdboth$chr)),]
l1kdboth<-l1kdboth[!(grepl("unk", l1kdboth$chr)),]; l1kdboth<-l1kdboth[!(grepl("chrM", l1kdboth$chr)),]
l1kdboth<-l1kdboth[!(grepl("chrY", l1kdboth$chr)),];

ctrchipboth<-ctrchipboth[!(grepl("alt", ctrchipboth$chr)),]; ctrchipboth<-ctrchipboth[!(grepl("rand", ctrchipboth$chr)),]
ctrchipboth<-ctrchipboth[!(grepl("unk", ctrchipboth$chr)),]; ctrchipboth<-ctrchipboth[!(grepl("chrM", ctrchipboth$chr)),]
ctrchipboth<-ctrchipboth[!(grepl("chrY", ctrchipboth$chr)),];



##need loop to compare 50bp bins between reps in the chosen regions.
#set up l1kd table
write.table(c("chr\trep1start\trep1end\trep1qvaltrans\trep1summit\tstrand\trep2start\trep2end\trep2qvaltrans\trep2summit\tminstart\tmaxend\tnumBins\tl1kdinputMean\tl1kdrep1mean\tl1kdrep2mean\tctrchipinputMean\tctrchiprep1mean\tctrchiprep2mean"), 
            file="l1kd_peaks_withCPMmeans.txt", append=F, sep="\t", col.names=F, quote = F, row.names = F)

for(i in 1:nrow(l1kdboth)){
  print(i)
  currentchr<-l1kdboth[i,1]
  minstart<-l1kdboth[i,11]
  maxend<-l1kdboth[i,12]
  subsetofcpm<-subset(cpmtable, cpmtable$chr == currentchr)
  subsetofcpm<-subset(subsetofcpm, subsetofcpm$start > minstart & subsetofcpm$end < maxend ) #might lose ends of peak but ok 
  #get input mean, rep1 mean, and rep 2 mean, and number of bins used
  subrows<-nrow(subsetofcpm)
  sub_inputmean<-mean(na.omit(subsetofcpm$l1kdinput))
  sub_rep1mean<-mean(na.omit(subsetofcpm$l1kd_rep1))
  sub_rep2mean<-mean(na.omit(subsetofcpm$l1kd_rep2))
  sub_ctrinputmean<-mean(na.omit(subsetofcpm$ctrinput))
  sub_ctrrep1mean<-mean(na.omit(subsetofcpm$ctr_rep1))
  sub_ctrrep2mean<-mean(na.omit(subsetofcpm$ctr_rep2))

  write.table(c(l1kdboth[i,],subrows,sub_inputmean, sub_rep1mean, sub_rep2mean, sub_ctrinputmean, sub_ctrrep1mean, sub_ctrrep2mean), 
              file="l1kd_peaks_withCPMmeans_original.txt", append=T, sep="\t", col.names=F, quote = F, row.names = F)
  
}

#ctr-chip table
write.table(c("chr\trep1start\trep1end\trep1qvaltrans\trep1summit\tstrand\trep2start\trep2end\trep2qvaltrans\trep2summit\tminstart\tmaxend\tnumBins\tl1kdinputMean\tl1kdrep1mean\tl1kdrep2mean\tctrchipinputMean\tctrchiprep1mean\tctrchiprep2mean"), 
            file="ctrchip_peaks_withCPMmeans.txt", append=F, sep="\t", col.names=F, quote = F, row.names = F)

###repeat this loop for the ctr list
for(i in 1:nrow(ctrchipboth)){
  print(i)
  currentchr<-ctrchipboth[i,1]
  minstart<-ctrchipboth[i,11]
  maxend<-ctrchipboth[i,12]
  subsetofcpm<-subset(cpmtable, cpmtable$chr == currentchr)
  subsetofcpm<-subset(subsetofcpm, subsetofcpm$start > minstart & subsetofcpm$end < maxend ) #might lose ends of peak but ok 
  #get input mean, rep1 mean, and rep 3 mean, and number of bins used
  subrows<-nrow(subsetofcpm)
  sub_inputmean<-mean(na.omit(subsetofcpm$l1kdinput))
  sub_rep1mean<-mean(na.omit(subsetofcpm$l1kd_rep1))
  sub_rep2mean<-mean(na.omit(subsetofcpm$l1kd_rep2))
  sub_ctrinputmean<-mean(na.omit(subsetofcpm$ctrinput))
  sub_ctrrep1mean<-mean(na.omit(subsetofcpm$ctr_rep1))
  sub_ctrrep2mean<-mean(na.omit(subsetofcpm$ctr_rep2))

  write.table(c(ctrchipboth[i,],subrows,sub_inputmean, sub_rep1mean, sub_rep2mean, sub_ctrinputmean, sub_ctrrep1mean, sub_ctrrep2mean), 
              file="ctrchip_peaks_withCPMmeans_original.txt", append=T, sep="\t", col.names=F, quote = F, row.names = F)
  
}


####reload tables of regions and assess replicability####
chrorder<-c("chr1",  "chr2", "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", 
          "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

###get stringent vers. When comparing two values, the lower one should be at least 0.8x the higher one. 
l1kdpeakbins_strict<-read.table("l1kd_peaks_withCPMmeans_original.txt", header = T, sep = "\t")
head(l1kdpeakbins_strict)

l1kdpeakbins_strict$l1kd_rep1to2diff<-abs(l1kdpeakbins_strict$l1kdrep2mean - l1kdpeakbins_strict$l1kdrep1mean)
l1kdpeakbins_strict$min80percentofhigherrep<-"fail"
l1kdpeakbins_strict$l1kdrep1to2max<-apply(l1kdpeakbins_strict[,c("l1kdrep2mean", "l1kdrep1mean")], 1, function(x) max(x))
l1kdpeakbins_strict$half_l1kdrep1to2max<-0.2 * l1kdpeakbins_strict$l1kdrep1to2max
l1kdpeakbins_strict[l1kdpeakbins_strict$l1kd_rep1to2diff < l1kdpeakbins_strict$half_l1kdrep1to2max,]$min80percentofhigherrep<-"pass"

##and check if they're both higher than input 
passingl1kdpeaksstrict<-subset(l1kdpeakbins_strict, l1kdpeakbins_strict$min80percentofhigherrep == "pass")
#save the file 
write.table(passingl1kdpeaksstrict, file="passingl1kdpeaks_80strict.txt", sep = "\t", col.names = T, row.names = F, quote = F)


####Ctr-Chip Peak list filtering#####
###repeat this process for ctr chip peak list.
ctrchippeakbins<-read.table("ctrchip_peaks_withCPMmeans_original.txt", header = T, sep = "\t")

####get strict peak list for Ctr ChIP 
ctrchippeakbins_strict<-read.table("ctrchip_peaks_withCPMmeans_original.txt", header = T, sep = "\t")
head(ctrchippeakbins_strict)

ctrchippeakbins_strict$ctrchip_rep1to2diff<-abs(ctrchippeakbins_strict$ctrchiprep2mean - ctrchippeakbins_strict$ctrchiprep1mean)
ctrchippeakbins_strict$min80percentofhigherrep<-"fail"
ctrchippeakbins_strict$ctrchiprep1to2max<-apply(ctrchippeakbins_strict[,c("ctrchiprep2mean", "ctrchiprep1mean")], 1, function(x) max(x))
ctrchippeakbins_strict$half_ctrchiprep1to2max<-0.2 * ctrchippeakbins_strict$ctrchiprep1to2max
#if the abs val of diff is more than half the max ?
#this is like saying, the difference between rep1-2 has to be less than 20% of the max value, to be close enough
#because then it would be with 80%
ctrchippeakbins_strict[ctrchippeakbins_strict$ctrchip_rep1to2diff < ctrchippeakbins_strict$half_ctrchiprep1to2max,]$min80percentofhigherrep<-"pass"

#save file.
passingctrchippeaksstrict<-subset(ctrchippeakbins_strict, ctrchippeakbins_strict$min80percentofhigherrep == "pass")
write.table(passingctrchippeaksstrict, file="passingctrchippeaks_80strict.txt", sep = "\t", col.names = T, row.names = F, quote = F)



#####separate the L1kd filtered list and the ctr chip filtered lists in overlapping and non overlapping.#####
##assign ids
passingctrchippeaksstrict$strictpeakid<-seq(1:nrow(passingctrchippeaksstrict))
passingl1kdpeaksstrict$strictpeakid<-seq(1:nrow(passingl1kdpeaksstrict))


#get sub set of cols for each so we can do bed intersect 
##use min and max range values 
passingctrchippeaksstrictforbed<-passingctrchippeaksstrict[,c("chr", "minstart", "maxend", "strictpeakid")]
passingl1kdpeaksstrictforbed<-passingl1kdpeaksstrict[,c("chr", "minstart", "maxend", "strictpeakid")]
write.table(passingctrchippeaksstrictforbed, file = "strict_ctrchip_peaks_bed.bed", sep = "\t", quote=F, row.names = F, col.names = F)
write.table(passingl1kdpeaksstrictforbed, file = "strict_l1kdchip_peaks_bed.bed", sep = "\t", quote=F, row.names = F, col.names = F)

#####end of selecting strict peaks. #####



#####Final set of plots: use the strict ctr-ChIP list even if it overlaps with the L1KD######
##reload and work on graphs.###
##load the all cpm table since that has to be plotted
#this is the all cpm reported as means 
allcpmwithmeans<-read.table("option1tablewithmeansNAs.txt", header = T)

#this is our final set of strict regions in the two conditions 
allstrict<-read.table("final50bpbincpmtable_allstrictregionscombined.bed")
colnames(allstrict)<-colnames(allcpmwithmeans)

#this is the all bin info 
allcpm<-read.table("Average_CPM_bins50bp_MBWSummary_Option1.txt")
colnames(allcpm)<-c("chr", "start", "end", "ctr-chip-input", "ctr-chip-rep1", "ctr-chip-rep2", "l1kd-input", "l1kd-rep1", "l1kd-rep2")


###drop unused chr from these files 
#drop extra chr
allcpm<-allcpm[!(grepl("alt", allcpm$chr)),]; allcpm<-allcpm[!(grepl("rand", allcpm$chr)),]
allcpm<-allcpm[!(grepl("unk", allcpm$chr)),]; allcpm<-allcpm[!(grepl("chrM", allcpm$chr)),]
allcpm<-allcpm[!(grepl("chrY", allcpm$chr)),];

allstrict<-allstrict[!(grepl("alt", allstrict$chr)),]; allstrict<-allstrict[!(grepl("rand", allstrict$chr)),]
allstrict<-allstrict[!(grepl("unk", allstrict$chr)),]; allstrict<-allstrict[!(grepl("chrM", allstrict$chr)),]
allstrict<-allstrict[!(grepl("chrY", allstrict$chr)),];

allcpmwithmeans<-allcpmwithmeans[!(grepl("alt", allcpmwithmeans$chr)),]; allcpmwithmeans<-allcpmwithmeans[!(grepl("rand", allcpmwithmeans$chr)),]
allcpmwithmeans<-allcpmwithmeans[!(grepl("unk", allcpmwithmeans$chr)),]; allcpmwithmeans<-allcpmwithmeans[!(grepl("chrM", allcpmwithmeans$chr)),]
allcpmwithmeans<-allcpmwithmeans[!(grepl("chrY", allcpmwithmeans$chr)),];



####Plot 1:  boxplot of K3K27me3 peaks by Lad, nad type, using strict peak lists####
allstrict<-unique(allstrict)

ladsstrict<-read.table("final50bpbincpmtable_allstrictregionscombined_ladsonly.bed")
nadsstrict<-read.table("final50bpbincpmtable_allstrictregionscombined_nadsonly.bed")
ladsnadsstrict<-read.table("final50bpbincpmtable_allstrictregionscombined_ladsnadsoverlaps.bed")
noladsnadsstrict<-read.table("final50bpbincpmtable_allstrictregionscombined_NoLadsNads.bed")
#colnames
colnames(ladsstrict)<-colnames(allstrict)
colnames(nadsstrict)<-colnames(allstrict)
colnames(ladsnadsstrict)<-colnames(allstrict)
colnames(noladsnadsstrict)<-colnames(allstrict)

ladsstrict<-unique(ladsstrict)
nadsstrict<-unique(nadsstrict)
ladsnadsstrict<-unique(ladsnadsstrict)
noladsnadsstrict<-unique(noladsnadsstrict)

##order is all h3k27me3, then lads, coladnads, nads, then nonladnads
allstrictplot1<-allstrict
allstrictplot1$sites<-"All H3K27me3"
ladsstrictplot1<-ladsstrict
ladsstrictplot1$sites<-"LADs H3K27me3"
ladsnadsstrictplot1<-ladsnadsstrict
ladsnadsstrictplot1$sites<-"co-LADs/NADs H3K27me3"
nadsstrictplot1<-nadsstrict
nadsstrictplot1$sites<-"NADs H3K27me3"
noladsnadsstrictplot1<-noladsnadsstrict
noladsnadsstrictplot1$sites<-"non-LADs/NADs H3K27me3"

#rbind
allinfoplot1<-rbind(allstrictplot1, ladsstrictplot1, ladsnadsstrictplot1, nadsstrictplot1, noladsnadsstrictplot1 )
sitesorder<-c("All H3K27me3", "LADs H3K27me3", "co-LADs/NADs H3K27me3",   "NADs H3K27me3",  "non-LADs/NADs H3K27me3" )

allinfoplot1$chr_F<-factor(allinfoplot1$chr, levels=chrorder)
allinfoplot1$sites_F<-factor(allinfoplot1$sites, levels=sitesorder)

##quick plot of the number of peaks per chromosomes 
pexp1<-ggplot(allinfoplot1, aes(chr_F)) + geom_bar(position = "dodge") + facet_grid(sites_F~.) + 
      ggtitle("Strict ChIP Peak regions, as 50bp bins") + #scale_fill_manual(values = colopts) +
  xlab("chromosome") + ylab("count") + theme_bw() + labs(fill = "Sample") + scale_x_discrete(drop=F) +#coord_cartesian(xlim = c(-8,-1)) +
    theme(axis.title = element_text(size=12), legend.position = c("bottom"), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 6,  hjust = 0.5),
        axis.text.y = element_text(hjust = 1, size = 8))
pdf("strictPeakList_chrBreakdown_LADsNADs.pdf", width = 8, height = 10)
pexp1
dev.off()

###Sub plots to go with plot 1: 
#plot 1. a figure like above per chr, but show ctr counts and L1KD counts side by side? 
#plot 2. a figure (do not split the counts per chromosome), and plot it like the actual Plot#1 (instead of log2CPM, plot counts#). 
#need to go back to strict lists per condition and merge to the 4 lad region lists, then get #bins per those as needed.

###load all bed files from bedtools intersect in.. then rbind the regions with info on l1kd or ctr, then with site region.

ctr_allh3k27me3_50bp<-read.table("50bpbintable_strict_ctr.bed")
ctr_allh3k27me3_50bp<-unique(ctr_allh3k27me3_50bp)

l1kd_allh3k27me3_50bp<-read.table("50bpbintable_strict_l1kd.bed")
l1kd_allh3k27me3_50bp<-unique(l1kd_allh3k27me3_50bp)

ctr_nads_50bp<-read.table("50bpbintable_strict_ctr_nadsonly.bed")
ctr_nads_50bp<-unique(ctr_nads_50bp)
ctr_lads_50bp<-read.table("50bpbintable_strict_ctr_lads_only.bed")
ctr_lads_50bp<-unique(ctr_lads_50bp)
ctr_ladsnads_50bp<-read.table("50bpbintable_strict_ctr_ladsnadsoverlaps_only.bed")
ctr_ladsnads_50bp<-unique(ctr_ladsnads_50bp)
ctr_noladsnads_50bp<-read.table("50bpbintable_strict_ctr_nonladnads_only.bed")
ctr_noladsnads_50bp<-unique(ctr_noladsnads_50bp)

l1kd_nads_50bp<-read.table("50bpbintable_strict_l1kd_nads_only.bed")
l1kd_nads_50bp<-unique(l1kd_nads_50bp)
l1kd_lads_50bp<-read.table("50bpbintable_strict_l1kd_lads_only.bed")
l1kd_lads_50bp<-unique(l1kd_lads_50bp)
l1kd_ladsnads_50bp<-read.table("50bpbintable_strict_l1kd_ladsnadsoverlaps_only.bed")
l1kd_ladsnads_50bp<-unique(l1kd_ladsnads_50bp)
l1kd_noladsnads_50bp<-read.table("50bpbintable_strict_l1kd_nonladnads_only.bed")
l1kd_noladsnads_50bp<-unique(l1kd_noladsnads_50bp)


##add sites info and sample info 
l1kd_allh3k27me3_50bp$site<-"All H3K27me3"
ctr_allh3k27me3_50bp$site<-"All H3K27me3"
ctr_nads_50bp$site<-"H3K27me3 + NADs"
ctr_lads_50bp$site<-"H3K27me3 + LADs"
ctr_ladsnads_50bp$site<-"H3K27me3 + co-LADs/NADs"
ctr_noladsnads_50bp$site<-"H3K27me3 + non-LADs/NADs"
l1kd_nads_50bp$site<-"H3K27me3 + NADs"
l1kd_lads_50bp$site<-"H3K27me3 + LADs"
l1kd_ladsnads_50bp$site<-"H3K27me3 + co-LADs/NADs"
l1kd_noladsnads_50bp$site<-"H3K27me3 + non-LADs/NADs"

l1kd_allh3k27me3_50bp$sample<-"L1KD"
ctr_allh3k27me3_50bp$sample<-"Ctr-ChIP"
ctr_nads_50bp$sample<-"Ctr-ChIP"
ctr_lads_50bp$sample<-"Ctr-ChIP"
ctr_ladsnads_50bp$sample<-"Ctr-ChIP"
ctr_noladsnads_50bp$sample<-"Ctr-ChIP"
l1kd_nads_50bp$sample<-"L1KD"
l1kd_lads_50bp$sample<-"L1KD"
l1kd_ladsnads_50bp$sample<-"L1KD"
l1kd_noladsnads_50bp$sample<-"L1KD"


##rbind all
allindiv50bpbins<-rbind(ctr_allh3k27me3_50bp, ctr_nads_50bp,ctr_lads_50bp,ctr_ladsnads_50bp,ctr_noladsnads_50bp,l1kd_allh3k27me3_50bp, l1kd_nads_50bp, l1kd_lads_50bp, l1kd_ladsnads_50bp, l1kd_noladsnads_50bp)
colnames(allindiv50bpbins)<-c("chr","start","end", "ctrinput", "l1kdinput", "l1kdr1and2mean", "ctr1and2mean","sites", "sample")

###plot
allindiv50bpbins$chr_F<-factor(allindiv50bpbins$chr, levels=chrorder)
ladnadonlysiteorder<-c("All H3K27me3", "H3K27me3 + LADs","H3K27me3 + co-LADs/NADs", "H3K27me3 + NADs", "H3K27me3 + non-LADs/NADs")
allindiv50bpbins$sites_F<-factor(allindiv50bpbins$sites, levels=ladnadonlysiteorder)

##quick plot of the number of peaks per chromosomes 
pexp1.1<-ggplot(allindiv50bpbins, aes(chr_F, fill = sample)) + geom_bar(position = position_dodge2(preserve = "single")) + facet_grid(sites_F~.) + 
      ggtitle("Strict ChIP Peak regions, as count of 50bp bins") + #scale_fill_manual(values = colopts) +
  xlab("chromosome") + ylab("count of bins") + theme_bw() + labs(fill = "Sample") + scale_x_discrete(drop=F) +#coord_cartesian(xlim = c(-8,-1)) +
    theme(axis.title = element_text(size=12), legend.position = c("bottom"), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), strip.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 6,  hjust = 0.5), 
        axis.text.y = element_text(hjust = 1, size = 8))
pdf("strictPeakLists_chrBreakdown_individual_CTR-ChIP_and_L1KD_LADsNADs.pdf", width = 8, height = 10)
pexp1.1
dev.off()

##boxplot showing the regions split between L1KD vs ctr across the 4 lad nad regions####
pexp1.2<-ggplot(allindiv50bpbins, aes(sites_F, fill = sample)) + geom_bar(position = position_dodge2(preserve = "single")) + #facet_grid(sites_F~.) + 
      ggtitle("Strict ChIP Peak regions, as count of 50bp bins") + #scale_fill_manual(values = colopts) +
  xlab("Region Group") + ylab("count of bins") + theme_bw() + labs(fill = "Sample") + scale_x_discrete(drop=F) +#coord_cartesian(xlim = c(-8,-1)) +
    theme(axis.title = element_text(size=12), legend.position = c("bottom"), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5), #legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), strip.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 6,  hjust = 1, angle = 45), 
        axis.text.y = element_text(hjust = 1, size = 8))
pdf("strictPeakLists_individual_CTR-ChIP_and_L1KD_LADsNADs.pdf", width =5, height = 5)
pexp1.2
dev.off()

##give lists of numbers on the plots
persamplesite<-allindiv50bpbins %>% dplyr::count(sample, sites)
persamplesitechr<-allindiv50bpbins %>% dplyr::count(sample, sites, chr)
persamplesitechrwide<-spread(persamplesitechr, sample, n)

####end of extra graphs showing diff between LADsNADS for the two conditions####
##continue on to make plot 1

#drop the inputs
allinfoplot1<-allinfoplot1[,c("chr","start","end", "sites", "ctrinput","l1kdinput","l1kdr1and2mean", "ctr1and2mean"  )]

#now make long, then graph as samples and as means, and as chr for both
allinfoplot1$chr_F<-NULL
allinfoplot1$sites_F<-NULL

allinfoplot1long<-melt(allinfoplot1, id.vars=c("chr", "start", "end", "sites"))
colnames(allinfoplot1long)<-c("chr", "start", "end", "sites", "sample", "CPM")

allinfoplot1long$chr_F<-factor(allinfoplot1long$chr, levels=chrorder)
allinfoplot1long$sites_F<-factor(allinfoplot1long$sites, levels=sitesorder)

allinfoplot1long<-subset(allinfoplot1long, allinfoplot1long$sample %in% c("l1kdr1and2mean", "ctr1and2mean"))
allinfoplot1long$sample_F<-factor(allinfoplot1long$sample, levels= c( "ctr1and2mean", "l1kdr1and2mean"))
allinfoplot1long<-allinfoplot1long[allinfoplot1long$CPM > 0,]
allinfoplot1long<-allinfoplot1long[!(is.na(allinfoplot1long$CPM)),]
allinfoplot1long$CPM<-as.numeric(allinfoplot1long$CPM)
allinfoplot1long$log2CPM<-log2(allinfoplot1long$CPM)

colopts<-c("ctr1and2mean" = "gray", "l1kdr1and2mean" = "turquoise")

pdf("plot1_strictpeaks_ladsnadsinfo.pdf", width = 5, height = 5)
plotchip.1<-ggplot(allinfoplot1long, aes(x = sites_F, y = log2CPM, fill = sample_F)) + geom_boxplot(position="dodge", outlier.shape = NA) +
      ggtitle(paste("ChIP: log2 CPM per 50bp interval for peak region subsets\n(outliers not shown)")) + scale_fill_manual(values = colopts) +
  xlab("Sample") + ylab("log2 CPM") + theme_bw() + coord_cartesian(ylim = c(-6,1)) +
    theme(axis.title = element_text(size=12), legend.position = c("none"), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5, size = 8), #legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 7,  hjust = 1, angle = 45),  axis.text.y = element_text(hjust = 1, size = 8))
plotchip.1
dev.off()

###need t tests for these 
#I assume it is just between ctr and l1kd for each site group since we don't show chr here?
plot1h3k27me3subset<-allinfoplot1long[allinfoplot1long$sites == "All H3K27me3",]
plot1h3k27me3subset.tt<-t.test(plot1h3k27me3subset[plot1h3k27me3subset$sample == "l1kdr1and2mean",]$CPM, 
       plot1h3k27me3subset[plot1h3k27me3subset$sample == "ctr1and2mean",]$CPM, alternative = "two.sided", paired=F)

plot1ladsubset<-allinfoplot1long[allinfoplot1long$sites == "LADs H3K27me3",]
plot1ladsubset.tt<-t.test(plot1ladsubset[plot1ladsubset$sample == "l1kdr1and2mean",]$CPM, 
       plot1ladsubset[plot1ladsubset$sample == "ctr1and2mean",]$CPM, alternative = "two.sided", paired=F)

plot1coladnadsubset<-allinfoplot1long[allinfoplot1long$sites == "co-LADs/NADs H3K27me3",]
plot1coladnadsubset.tt<-t.test(plot1coladnadsubset[plot1coladnadsubset$sample == "l1kdr1and2mean",]$CPM, 
       plot1coladnadsubset[plot1coladnadsubset$sample == "ctr1and2mean",]$CPM, alternative = "two.sided", paired=F)

plot1nadsubset<-allinfoplot1long[allinfoplot1long$sites == "NADs H3K27me3",]
plot1nadsubset.tt<-t.test(plot1nadsubset[plot1nadsubset$sample == "l1kdr1and2mean",]$CPM, 
       plot1nadsubset[plot1nadsubset$sample == "ctr1and2mean",]$CPM, alternative = "two.sided", paired=F)

plot1nonladnadsubset<-allinfoplot1long[allinfoplot1long$sites == "non-LADs/NADs H3K27me3",]
plot1nonladnadsubset.tt<-t.test(plot1nonladnadsubset[plot1nonladnadsubset$sample == "l1kdr1and2mean",]$CPM, 
       plot1nonladnadsubset[plot1nonladnadsubset$sample == "ctr1and2mean",]$CPM, alternative = "two.sided", paired=F)


print(c("t-test for Ctr vs L1KD, for CPM per 50bp interval from strict peaks region lists\nchr; tstat; pvalue; df"))
print(c("H3K27me3", plot1h3k27me3subset.tt$statistic, plot1h3k27me3subset.tt$p.value, plot1h3k27me3subset.tt$parameter))
print(c("LADs H3K27me3", plot1ladsubset.tt$statistic, plot1ladsubset.tt$p.value, plot1ladsubset.tt$parameter))
print(c("co-LADs/NADs H3K27me3", plot1coladnadsubset.tt$statistic, plot1coladnadsubset.tt$p.value, plot1coladnadsubset.tt$parameter))
print(c("NADs H3K27me3", plot1nadsubset.tt$statistic, plot1nadsubset.tt$p.value, plot1nadsubset.tt$parameter))
print(c("non-LADs/NADs H3K27me3", plot1nonladnadsubset.tt$statistic, plot1nonladnadsubset.tt$p.value, plot1nonladnadsubset.tt$parameter))

#for t tests add the comparison just among the ctr, per LADs/NADs categories (--> ctr-LADs/ctr-coLN/ctr-NADs/ctr-nonLN each to the ctr-All). 

plot1h3k27me3subsetCtr<-allinfoplot1long[allinfoplot1long$sample == "ctr1and2mean",]
plot1h3k27me3subsetCtr_LADs.tt<-t.test(plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "All H3K27me3",]$CPM, 
       plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "LADs H3K27me3",]$CPM, alternative = "two.sided", paired=F)
 plot1h3k27me3subsetCtr_NADs.tt<-t.test(plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "All H3K27me3",]$CPM, 
       plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "NADs H3K27me3",]$CPM, alternative = "two.sided", paired=F)
 plot1h3k27me3subsetCtr_coLADNADs.tt<-t.test(plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "All H3K27me3",]$CPM, 
       plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "co-LADs/NADs H3K27me3",]$CPM, alternative = "two.sided", paired=F)
 plot1h3k27me3subsetCtr_noLADNADs.tt<-t.test(plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "All H3K27me3",]$CPM, 
       plot1h3k27me3subsetCtr[plot1h3k27me3subsetCtr$sites == "non-LADs/NADs H3K27me3",]$CPM, alternative = "two.sided", paired=F)

print(c("LADs H3K27me3", plot1h3k27me3subsetCtr_LADs.tt$statistic, plot1h3k27me3subsetCtr_LADs.tt$p.value, plot1h3k27me3subsetCtr_LADs.tt$parameter))
print(c("co-LADs/NADs H3K27me3", plot1h3k27me3subsetCtr_coLADNADs.tt$statistic, plot1h3k27me3subsetCtr_coLADNADs.tt$p.value, plot1h3k27me3subsetCtr_coLADNADs.tt$parameter))
print(c("NADs H3K27me3", plot1h3k27me3subsetCtr_NADs.tt$statistic, plot1h3k27me3subsetCtr_NADs.tt$p.value, plot1h3k27me3subsetCtr_NADs.tt$parameter))
print(c("non-LADs/NADs H3K27me3", plot1h3k27me3subsetCtr_noLADNADs.tt$statistic, plot1h3k27me3subsetCtr_noLADNADs.tt$p.value, plot1h3k27me3subsetCtr_noLADNADs.tt$parameter))


 
#####plot 2, all strict only###
#drop inputs
allstrictforplot2<-allstrict
allstrictforplot2<-allstrictforplot2[,c("chr","start","end", "l1kdr1and2mean", "ctr1and2mean")]
allstrictforplot2long<-melt(allstrictforplot2, id.vars=c("chr", "start", "end"))
colnames(allstrictforplot2long)<-c("chr", "start", "end", "sample", "CPM")

allstrictforplot2long$chr_F<-factor(allstrictforplot2long$chr, levels=chrorder)

allstrictforplot2long$sample_F<-factor(allstrictforplot2long$sample, levels= c( "ctr1and2mean", "l1kdr1and2mean"))
allstrictforplot2long<-allstrictforplot2long[allstrictforplot2long$CPM > 0,]
allstrictforplot2long<-allstrictforplot2long[!(is.na(allstrictforplot2long$CPM)),]
allstrictforplot2long$CPM<-as.numeric(allstrictforplot2long$CPM)
allstrictforplot2long$log2CPM<-log2(allstrictforplot2long$CPM)

pdf("plot2_strictpeaks_h3k27me3_perChr.pdf", width = 9, height = 5)
plotchip.2<-ggplot(allstrictforplot2long, aes(x = chr_F, y = log2CPM, fill = sample_F)) + geom_boxplot(position="dodge", outlier.shape = NA) +
      ggtitle(paste("ChIP: log2 CPM per 50bp interval from strict H3K27me3 region list\n(outliers not shown)")) + scale_fill_manual(values = colopts) +
  xlab("Chromosome") + ylab("log2CPM") + theme_bw() + labs(fill = "Sample") + scale_x_discrete(drop=F) + coord_cartesian(ylim = c(-6,2)) +
    theme(axis.title = element_text(size=12), legend.position = c("bottom"), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5, size = 9), #legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 8,  hjust = 1, angle = 45),
        axis.text.y = element_text(hjust = 1, size = 8))
plotchip.2
dev.off()


#t tests
write.table(c("t-test outs for log2 CPM per 50bp interval from all strict k27me3 peaks region list\nchr; tstat; pvalue; df"), file = "plot2_ttestandStats_allstrict.txt", append=F, quote = F)
for(i in chrorder){
  subs<-subset(allstrictforplot2long, allstrictforplot2long$chr == i)
  
  ttestres<-t.test(subs[subs$sample == "l1kdr1and2mean",]$CPM, subs[subs$sample == "ctr1and2mean",]$CPM, paired = F, alternative = "t")
  print(c(i, ttestres$statistic, ttestres$p.value, ttestres$parameter))
  write.table(paste(i, ttestres$statistic, ttestres$p.value, ttestres$parameter), file = "plot2_ttestandStats_allstrict.txt", append=T, quote = F)

}

t.test(allstrictforplot2long[allstrictforplot2long$sample == "l1kdr1and2mean",]$CPM, allstrictforplot2long[allstrictforplot2long$sample == "ctr1and2mean",]$CPM, paired = F, alternative = "t", na.action=na.omit)
#t = -40.608, df = 119370, p-value < 0.00000000000000022

##also report mean and sd of the values per samp per chr.
write.table(c("\nstats for the samples: mean of L1KD, SD of L1KD, mean of ctrChip, SD of ctrChip\n"), file = "plot2_ttestandStats_allstrict.txt", append=T, quote = F)
for(i in chrorder){
  subs<-subset(allstrictforplot2long, allstrictforplot2long$chr == i)
  meanvall1kd<-mean(na.omit(subs[subs$sample == "l1kdr1and2mean",]$CPM))
  meanvalctr<-mean(na.omit(subs[subs$sample == "ctr1and2mean",]$CPM))
  sdctr<-sd(na.omit(subs[subs$sample == "ctr1and2mean",]$CPM))
  sdl1kd<-sd(na.omit(subs[subs$sample == "l1kdr1and2mean",]$CPM))
  write.table(paste(i, meanvall1kd, sdl1kd, meanvalctr, sdctr), file = "plot2_ttestandStats_allstrict.txt", append=T, quote = F)
} 

###end plot 2


####plot 3 is all reads OR all h3k27me3 with rep1 and rep2 beside each other to show the reps agree####
#need to select from the all non mean cpm table with filtered peaks...
##also include t test
#drop extra chr

allstrict<-unique(allstrict)
allstrict$chrstartend<-paste(allstrict$chr, allstrict$start, allstrict$end, sep = "_")

allcpm_copy<-allcpm
allcpm_copy$chrstartend<-paste(allcpm_copy$chr, allcpm_copy$start, allcpm_copy$end, sep = "_")

allcpm_copy<-subset(allcpm_copy, allcpm_copy$chrstartend %in% allstrict$chrstartend)
allcpm_copy<-allcpm_copy[,c("chr","start","end", "ctr-chip-rep1", "ctr-chip-rep2",  "l1kd-rep1", "l1kd-rep2")]

allcpm_copylong<-melt(allcpm_copy, id.vars=c("chr", "start", "end"))
colnames(allcpm_copylong)<-c("chr", "start", "end", "sample", "CPM")

allcpm_copylong$chr_F<-factor(allcpm_copylong$chr, levels=chrorder)
allcpm_copylong$sample_F<-factor(allcpm_copylong$sample, levels= c("ctr-chip-rep1", "ctr-chip-rep2","l1kd-rep1","l1kd-rep2" ))
allcpm_copylong<-allcpm_copylong[allcpm_copylong$CPM > 0,]
allcpm_copylong<-allcpm_copylong[!(is.na(allcpm_copylong$CPM)),]
allcpm_copylong$CPM<-as.numeric(allcpm_copylong$CPM)
allcpm_copylong$log2CPM<-log2(allcpm_copylong$CPM)


allcpm_copylong$condition<-"L1KD"
allcpm_copylong[allcpm_copylong$sample %in% c("ctr-chip-rep1", "ctr-chip-rep2"),]$condition<-"Ctr-ChIP"

colopts<-c("Ctr-ChIP" = "gray", "L1KD" = "turquoise")

pdf("plot3_strictpeaks_individualReps.pdf", width = 5, height = 5) #geom_bar(position = position_dodge2(preserve = "single"))
plotchip.3<-ggplot(allcpm_copylong, aes(x = sample_F, y = log2CPM, fill = condition)) + geom_boxplot(outlier.shape = NA) +
      ggtitle(paste("ChIP: log2 CPM per 50bp interval for peak region subsets\n(outliers not shown)")) + scale_fill_manual(values = colopts) +
  xlab("Sample") + ylab("log2 CPM") + theme_bw() + coord_cartesian(ylim = c(-6,1)) +
    theme(axis.title = element_text(size=12), legend.position = c("none"), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5, size = 8), #legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(size = 7,  hjust = 1, angle = 45),  axis.text.y = element_text(hjust = 1, size = 8))
plotchip.3
dev.off()

##do t tests
ttestres<-t.test(allcpm_copylong[allcpm_copylong$sample == "ctr-chip-rep1",]$CPM, allcpm_copylong[allcpm_copylong$sample == "ctr-chip-rep2",]$CPM, paired = F, alternative = "t")
print(c( ttestres$statistic, ttestres$p.value, ttestres$parameter))

ttestres<-t.test(allcpm_copylong[allcpm_copylong$sample == "ctr-chip-rep1",]$CPM, allcpm_copylong[allcpm_copylong$sample == "l1kd-rep1",]$CPM, paired = F, alternative = "t")
print(c( ttestres$statistic, ttestres$p.value, ttestres$parameter))

ttestres<-t.test(allcpm_copylong[allcpm_copylong$sample == "ctr-chip-rep1",]$CPM, allcpm_copylong[allcpm_copylong$sample == "l1kd-rep2",]$CPM, paired = F, alternative = "t")
print(c( ttestres$statistic, ttestres$p.value, ttestres$parameter))

ttestres<-t.test(allcpm_copylong[allcpm_copylong$sample == "ctr-chip-rep2",]$CPM, allcpm_copylong[allcpm_copylong$sample == "l1kd-rep1",]$CPM, paired = F, alternative = "t")
print(c( ttestres$statistic, ttestres$p.value, ttestres$parameter))

ttestres<-t.test(allcpm_copylong[allcpm_copylong$sample == "ctr-chip-rep2",]$CPM, allcpm_copylong[allcpm_copylong$sample == "l1kd-rep2",]$CPM, paired = F, alternative = "t")
print(c( ttestres$statistic, ttestres$p.value, ttestres$parameter))

ttestres<-t.test(allcpm_copylong[allcpm_copylong$sample == "l1kd-rep1",]$CPM, allcpm_copylong[allcpm_copylong$sample == "l1kd-rep2",]$CPM, paired = F, alternative = "t")
print(c( ttestres$statistic, ttestres$p.value, ttestres$parameter))


###end plot 3.


##end of script.