#Necessary libraries
library(PopGenome)
library(pegas)
library(vcfR)
library(dplyr)
library(ggplot2)
library(tidyverse)

#Set working directory
setwd("~/Chilense/")


#Split VCF into scaffolds; this helps load in data
VCF_split_into_scaffolds("VCFNew/NoIndelsNoMissing.sorted.recode.vcf", 
                         "scaffoldVCFs")

#Corresponding GFF scaffolds
GFF_split_into_scaffolds("GFF/NoIndel.gff","scaffoldGFFs")

#Create genome object
genome <- readData("scaffoldVCFs/", format = "VCF",
                   gffpath = "scaffoldGFFs/", 
                   include.unknown = TRUE)


#Create population variables populated with list of individuals
pop1963 <- c("1963_t1", 
             "1963_t3",
             "1963_t5",
             "1963_t7",
             "1963_t9",
             "1963_t1.2", 
             "1963_t3.2",
             "1963_t5.2",
             "1963_t7.2",
             "1963_t9.2")

pop2931 <- c("2931_t2",
             "2931_t3",
             "2931_t4",
             "2931_t5",
             "2931_t6",
             "2931_t2.2",
             "2931_t3.2",
             "2931_t4.2",
             "2931_t5.2",
             "2931_t6.2")

pop2932 <- c("2932_1",
             "2932_12",
             "2932_20",
             "2932_22",
             "2932_8",
             "2932_1.2",
             "2932_12.2",
             "2932_20.2",
             "2932_22.2",
             "2932_8.2")

pop3111 <- c("3111_t10",
             "3111_t15",
             "3111_t3",
             "3111_t5",
             "3111_t9",
             "3111_t10.2",
             "3111_t15.2",
             "3111_t3.2",
             "3111_t5.2",
             "3111_t9.2")

pop4107 <- c("4107_3",
             "4107_6",
             "4107_9",
             "4107_t11",
             "4107_t5",
             "4107_3.2",
             "4107_6.2",
             "4107_9.2",
             "4107_t11.2",
             "4107_t5.2")

pop4117A <- c("4117A_1",
              "4117A_10",
              "4117A_15",
              "4117A_4",
              "4117A_5",
              "4117A_1.2",
              "4117A_10.2",
              "4117A_15.2",
              "4117A_4.2",
              "4117A_5.2")

pop4330 <- c("4330_t1",
             "4330_t12",
             "4330_t4",
             "4330_t6",
             "4330_t9",
             "4330_t1.2",
             "4330_t12.2",
             "4330_t4.2",
             "4330_t6.2",
             "4330_t9.2")

#Set populations
genome <- set.populations(genome, diploid = TRUE, list(c("4330_t1",
                                                         "4330_t12",
                                                         "4330_t4",
                                                         "4330_t6",
                                                         "4330_t9",
                                                         "4330_t1.2",
                                                         "4330_t12.2",
                                                         "4330_t4.2",
                                                         "4330_t6.2",
                                                         "4330_t9.2"),
                                                       c("4117A_1",
                                                         "4117A_10",
                                                         "4117A_15",
                                                         "4117A_4",
                                                         "4117A_5",
                                                         "4117A_1.2",
                                                         "4117A_10.2",
                                                         "4117A_15.2",
                                                         "4117A_4.2",
                                                         "4117A_5.2"),
                                                       c("4107_3",
                                                         "4107_6",
                                                         "4107_9",
                                                         "4107_t11",
                                                         "4107_t5",
                                                         "4107_3.2",
                                                         "4107_6.2",
                                                         "4107_9.2",
                                                         "4107_t11.2",
                                                         "4107_t5.2"),
                                                       c("3111_t10",
                                                         "3111_t15",
                                                         "3111_t3",
                                                         "3111_t5",
                                                         "3111_t9",
                                                         "3111_t10.2",
                                                         "3111_t15.2",
                                                         "3111_t3.2",
                                                         "3111_t5.2",
                                                         "3111_t9.2"),
                                                       c("2932_1",
                                                         "2932_12",
                                                         "2932_20",
                                                         "2932_22",
                                                         "2932_8",
                                                         "2932_1.2",
                                                         "2932_12.2",
                                                         "2932_20.2",
                                                         "2932_22.2",
                                                         "2932_8.2"),
                                                       c("2931_t2",
                                                         "2931_t3",
                                                         "2931_t4",
                                                         "2931_t5",
                                                         "2931_t6",
                                                         "2931_t2.2",
                                                         "2931_t3.2",
                                                         "2931_t4.2",
                                                         "2931_t5.2",
                                                         "2931_t6.2"),
                                                       c("1963_t1", 
                                                         "1963_t3",
                                                         "1963_t5",
                                                         "1963_t7",
                                                         "1963_t9",
                                                         "1963_t1.2", 
                                                         "1963_t3.2",
                                                         "1963_t5.2",
                                                         "1963_t7.2",
                                                         "1963_t9.2")) )


#Split genes on individual chromosomes into exons and introns
genome.splitChr1Exons <- split_data_into_GFF_features(genome,
                                                      "GFF/NoIndel.gff",
                                                      "Scaffold_12631_Chr1",
                                                      "exon")

genome.splitChr6Exons <- split_data_into_GFF_features(genome,
                                                      "GFF/NoIndel.gff",
                                                      "Scaffold_8756_Chr6",
                                                      "exon")
genome.splitChr5Exons <- split_data_into_GFF_features(genome,
                                                      "GFF/NoIndel.gff",
                                                      "Scaffold_12633_Chr5",
                                                      "exon")

genome.splitChr7Exons <- split_data_into_GFF_features(genome,
                                                      "GFF/NoIndel.gff",
                                                      "Scaffold_12636_Chr7",
                                                      "exon")

genome.splitChr9Exons <- split_data_into_GFF_features(genome,
                                                      "GFF/NoIndel.gff",
                                                      "Scaffold_12638_Chr9",
                                                      "exon")

genome.splitChr12Exons <- split_data_into_GFF_features(genome,
                                                       "GFF/NoIndel.gff",
                                                       "Scaffold_12629_Chr12",
                                                       "exon")

genome.splitChr11Exons <- split_data_into_GFF_features(genome,
                                                       "GFF/NoIndel.gff",
                                                       "Scaffold_12632_Chr11",
                                                       "gene")


genome.splitChr10Exons <- split_data_into_GFF_features(genome,
                                                       "GFF/NoIndel.gff",
                                                       "Scaffold_12637_Chr10",
                                                       "gene")


genome.splitChr8Exons <- split_data_into_GFF_features(genome,
                                                      "GFF/NoIndel.gff",
                                                      "Scaffold_12635_Chr8",
                                                      "gene")

genome.splitChr1Introns <- split_data_into_GFF_features(genome,
                                                        "GFF/NoIndel.gff",
                                                        "Scaffold_12631_Chr1",
                                                        "intron")

genome.splitChr6Introns <- split_data_into_GFF_features(genome,
                                                        "GFF/NoIndel.gff",
                                                        "Scaffold_8756_Chr6",
                                                        "intron")
genome.splitChr5Introns <- split_data_into_GFF_features(genome,
                                                        "GFF/NoIndel.gff",
                                                        "Scaffold_12633_Chr5",
                                                        "intron")


genome.splitChr7Introns <- split_data_into_GFF_features(genome,
                                                        "GFF/NoIndel.gff",
                                                        "Scaffold_12636_Chr7",
                                                        "intron")

genome.splitChr9Introns <- split_data_into_GFF_features(genome,
                                                        "GFF/NoIndel.gff",
                                                        "Scaffold_12638_Chr9",
                                                        "intron")

genome.splitChr12Introns <- split_data_into_GFF_features(genome,
                                                         "GFF/NoIndel.gff",
                                                         "Scaffold_12629_Chr12",
                                                         "intron")

genome.splitChr11Introns <- split_data_into_GFF_features(genome,
                                                         "GFF/NoIndel.gff",
                                                         "Scaffold_12632_Chr11",
                                                         "intron")


genome.splitChr10Introns <- split_data_into_GFF_features(genome,
                                                         "GFF/NoIndel.gff",
                                                         "Scaffold_12637_Chr10",
                                                         "intron")


genome.splitChr8Introns <- split_data_into_GFF_features(genome,
                                                        "GFF/NoIndel.gff",
                                                        "Scaffold_12635_Chr8",
                                                        "intron")

#"False" genes needed to initialize Gene.matrix for analysis, flag to remove downstream
row.names.remove <- c("Chr8 6693760 - 6693761",
                      "Chr11 58704695 - 58704696",
                      "Chr10 68424263 - 68424264")
row.names.remove.introns <- c("Chr5 3723258 - 3723259")

#Calculate Neutrality Statistics and save as table to extract row names (chromosome positions) from   
genome.splitChr1Exons <- neutrality.stats(genome.splitChr1Exons)
neut1 <- get.neutrality(genome.splitChr1Exons)[[1]]

genome.splitChr8Exons <- neutrality.stats(genome.splitChr8Exons)
neut8 <- get.neutrality(genome.splitChr8Exons)[[2]]

genome.splitChr11Exons <- neutrality.stats(genome.splitChr11Exons)
neut11 <- get.neutrality(genome.splitChr11Exons)[[2]]

genome.splitChr10Exons <- neutrality.stats(genome.splitChr10Exons)
neut10 <- get.neutrality(genome.splitChr10Exons)[[2]]

genome.splitChr6Exons <- neutrality.stats(genome.splitChr6Exons)
neut6 <- get.neutrality(genome.splitChr6Exons)[[2]]

genome.splitChr5Exons <- neutrality.stats(genome.splitChr5Exons)
neut5 <- get.neutrality(genome.splitChr5Exons)[[2]]

genome.splitChr7Exons <- neutrality.stats(genome.splitChr7Exons)
neut7 <- get.neutrality(genome.splitChr7Exons)[[1]]

genome.splitChr9Exons <- neutrality.stats(genome.splitChr9Exons)
neut9 <- get.neutrality(genome.splitChr9Exons)[[1]]

genome.splitChr12Exons <- neutrality.stats(genome.splitChr12Exons)
neut12 <- get.neutrality(genome.splitChr12Exons)[[1]]

genome.splitChr1Introns <- neutrality.stats(genome.splitChr1Introns)
neut1Introns <- get.neutrality(genome.splitChr1Introns)[[1]]

genome.splitChr6Introns <- neutrality.stats(genome.splitChr6Introns)
neut6Introns <- get.neutrality(genome.splitChr6Introns)[[2]]

genome.splitChr5Introns <- neutrality.stats(genome.splitChr5Introns)
neut5Introns <- get.neutrality(genome.splitChr5Introns)[[2]]

genome.splitChr7Introns <- neutrality.stats(genome.splitChr7Introns)
neut7Introns <- get.neutrality(genome.splitChr7Introns)[[1]]

genome.splitChr9Introns <- neutrality.stats(genome.splitChr9Introns)
neut9Introns <- get.neutrality(genome.splitChr9Introns)[[1]]

genome.splitChr12Introns <- neutrality.stats(genome.splitChr12Introns)
neut12Introns <- get.neutrality(genome.splitChr12Introns)[[1]]

#Create data frames for each chromosome where rows are chromosome positions and columns are populations

DataFrameListChromosomes <- list() #List to hold dataframes for each chromosome

ListTajimasChr12Pops <- list()
for(i in 1:7){
  ListTajimasChr12Pops[[i]] <- data.frame(genome.splitChr12Exons@Tajima.D[1:nrow(
    genome.splitChr12Exons@Tajima.D),i],
    row.names = paste("Chr12",rownames(neut12), sep = " "))
  
}

Chr12TajimaDataFrame <- dplyr::bind_cols(ListTajimasChr12Pops)
colnames(Chr12TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[1]] <- Chr12TajimaDataFrame

ListTajimasChr11Pops <- list()
for(i in 1:7){
  ListTajimasChr11Pops[[i]] <- data.frame(genome.splitChr11Exons@Tajima.D[1:nrow(
    genome.splitChr11Exons@Tajima.D),i],
    row.names = paste("Chr11",rownames(neut11), sep = " "))
  
}

Chr11TajimaDataFrame <- dplyr::bind_cols(ListTajimasChr11Pops)
colnames(Chr11TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[2]] <- Chr11TajimaDataFrame

ListTajimasChrPops <- list()
for(i in 1:7){
  ListTajimasChrPops[[i]] <- data.frame(genome.splitChr10Exons@Tajima.D[1:nrow(
    genome.splitChr10Exons@Tajima.D),i],
    row.names = paste("Chr10",rownames(neut10), sep = " "))
  
}

Chr10TajimaDataFrame <- dplyr::bind_cols(ListTajimasChrPops)
colnames(Chr10TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[3]] <- Chr10TajimaDataFrame

ListTajimasChrPops <- list()
for(i in 1:7){
  ListTajimasChrPops[[i]] <- data.frame(genome.splitChr1Exons@Tajima.D[1:nrow(
    genome.splitChr1Exons@Tajima.D),i],
    row.names = paste("Chr1",rownames(neut1), sep = " "))
  
}

Chr1TajimaDataFrame <- dplyr::bind_cols(ListTajimasChrPops)
colnames(Chr1TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[4]] <- Chr1TajimaDataFrame

ListTajimasChrPops <- list()
for(i in 1:7){
  ListTajimasChrPops[[i]] <- data.frame(genome.splitChr5Exons@Tajima.D[1:nrow(
    genome.splitChr5Exons@Tajima.D),i],
    row.names = paste("Chr5",rownames(neut5), sep = " "))
  
}

Chr5TajimaDataFrame <- dplyr::bind_cols(ListTajimasChrPops)
colnames(Chr5TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[5]] <- Chr5TajimaDataFrame

ListTajimasChrPops <- list()
for(i in 1:7){
  ListTajimasChrPops[[i]] <- data.frame(genome.splitChr6Exons@Tajima.D[1:nrow(
    genome.splitChr6Exons@Tajima.D),i],
    row.names = paste("Chr6",rownames(neut6), sep = " "))
  
}

Chr6TajimaDataFrame <- dplyr::bind_cols(ListTajimasChrPops)
colnames(Chr6TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[6]] <- Chr6TajimaDataFrame

ListTajimasChrPops <- list()
for(i in 1:7){
  ListTajimasChrPops[[i]] <- data.frame(genome.splitChr7Exons@Tajima.D[1:nrow(
    genome.splitChr7Exons@Tajima.D),i],
    row.names = paste("Chr7",rownames(neut7), sep = " "))
  
}

Chr7TajimaDataFrame <- dplyr::bind_cols(ListTajimasChrPops)
colnames(Chr7TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[7]] <- Chr7TajimaDataFrame

ListTajimasChrPops <- list()
for(i in 1:7){
  ListTajimasChrPops[[i]] <- data.frame(genome.splitChr8Exons@Tajima.D[1:nrow(
    genome.splitChr8Exons@Tajima.D),i],
    row.names = paste("Chr8",rownames(neut8), sep = " "))
  
}

Chr8TajimaDataFrame <- dplyr::bind_cols(ListTajimasChrPops)
colnames(Chr8TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[8]] <- Chr8TajimaDataFrame

ListTajimasChrPops <- list()
for(i in 1:7){
  ListTajimasChrPops[[i]] <- data.frame(genome.splitChr9Exons@Tajima.D[1:nrow(
    genome.splitChr9Exons@Tajima.D),i],
    row.names = paste("Chr9",rownames(neut9), sep = " "))
  
}

Chr9TajimaDataFrame <- dplyr::bind_cols(ListTajimasChrPops)
colnames(Chr9TajimaDataFrame) <- c(1:7)
DataFrameListChromosomes[[9]] <- Chr9TajimaDataFrame
TajimaExonMaster <- dplyr::bind_rows(DataFrameListChromosomes)
colnames(TajimaExonMaster) <- c("4330", "4117A", "4107", "3111", "2932", "2931", "1963")
TajimaExonMaster <- TajimaExonMaster[!(row.names(TajimaExonMaster)
                                             %in% row.names.remove), ]
##Convert to Tidyverse compatible data

TajimaExonMaster <- tibble::rownames_to_column(TajimaExonMaster, var = "position" )
TajimaExonMaster <- TajimaExonMaster %>%
  pivot_longer(c("4330", "4117A", "4107", "3111", "2932", "2931", "1963"), 
               names_to = "population", values_to = "tajimasD")

#Plot

ExonTajimasDPlot <- ggplot(TajimaExonMaster, aes(x = population, y = tajimasD)) +
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")
print(ExonTajimasDPlot + ggtitle("Exon Tajima's D"))

#Same for Introns

DataFrameListChromosomesIntrons <- list() #List to hold dataframes for each chromosome

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr12Introns@Tajima.D[1:nrow(
    genome.splitChr12Introns@Tajima.D),i],
    row.names = paste("Chr12",rownames(neut12Introns), sep = " "))
  
}

Chr12IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr12IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[1]] <- Chr12IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr11Introns@Tajima.D[1:nrow(
    genome.splitChr11Introns@Tajima.D),i],
    row.names = paste("Chr11",rownames(neut11Introns), sep = " "))
  
}

Chr11IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr11IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[2]] <- Chr11IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr10Introns@Tajima.D[1:nrow(
    genome.splitChr10Introns@Tajima.D),i],
    row.names = paste("Chr10",rownames(neut10Introns), sep = " "))
  
}

Chr10IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr10IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[3]] <- Chr10IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr1Introns@Tajima.D[1:nrow(
    genome.splitChr1Introns@Tajima.D),i],
    row.names = paste("Chr1",rownames(neut1Introns), sep = " "))
  
}

Chr1IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr1IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[4]] <- Chr1IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr5Introns@Tajima.D[1:nrow(
    genome.splitChr5Introns@Tajima.D),i],
    row.names = paste("Chr5",rownames(neut5Introns), sep = " "))
  
}

Chr5IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr5IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[5]] <- Chr5IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr6Introns@Tajima.D[1:nrow(
    genome.splitChr6Introns@Tajima.D),i],
    row.names = paste("Chr6",rownames(neut6Introns), sep = " "))
  
}

Chr6IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr6IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[6]] <- Chr6IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr7Introns@Tajima.D[1:nrow(
    genome.splitChr7Introns@Tajima.D),i],
    row.names = paste("Chr7",rownames(neut7Introns), sep = " "))
  
}

Chr7IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr7IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[7]] <- Chr7IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr8Introns@Tajima.D[1:nrow(
    genome.splitChr8Introns@Tajima.D),i],
    row.names = paste("Chr8",rownames(neut8Introns), sep = " "))
  
}

Chr8IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr8IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[8]] <- Chr8IntronTajimaDataFrame

ListTajimasChrIntronPops <- list()
for(i in 1:7){
  ListTajimasChrIntronPops[[i]] <- data.frame(genome.splitChr9Introns@Tajima.D[1:nrow(
    genome.splitChr9Introns@Tajima.D),i],
    row.names = paste("Chr9",rownames(neut9Introns), sep = " "))
  
}

Chr9IntronTajimaDataFrame <- dplyr::bind_cols(ListTajimasChrIntronPops)
colnames(Chr9IntronTajimaDataFrame) <- c(1:7)
DataFrameListChromosomesIntrons[[9]] <- Chr9IntronTajimaDataFrame

TajimaIntronMaster <- dplyr::bind_rows(DataFrameListChromosomesIntrons)

colnames(TajimaIntronMaster) <- c("4330", "4117A", "4107", "3111", "2932", "2931", "1963")
TajimaIntronMaster <- TajimaIntronMaster[!(row.names(TajimaIntronMaster)
                                       %in% row.names.remove.introns), ]
##Convert to Tidyverse compatible data

TajimaIntronMaster <- tibble::rownames_to_column(TajimaIntronMaster, var = "position" )
TajimaIntronMaster <- TajimaIntronMaster %>%
  pivot_longer(c("4330", "4117A", "4107", "3111", "2932", "2931", "1963"), 
               names_to = "population", values_to = "tajimasD")

#Plot

IntronTajimasDPlot <- ggplot(TajimaIntronMaster, aes(x = population, y = tajimasD)) +
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")
print(IntronTajimasDPlot +  ggtitle("Intron Tajima's D"))

########################################################
#Nucleotide Diversity Analysis
########################################################

genome.splitChr1Exons <- diversity.stats(genome.splitChr1Exons)
divers1 <- get.diversity(genome.splitChr1Exons)[[1]]

genome.splitChr8Exons <- diversity.stats(genome.splitChr8Exons)
divers8 <- get.diversity(genome.splitChr8Exons)[[2]]

genome.splitChr11Exons <- diversity.stats(genome.splitChr11Exons)
divers11 <- get.diversity(genome.splitChr11Exons)[[2]]

genome.splitChr10Exons <- diversity.stats(genome.splitChr10Exons)
divers10 <- get.diversity(genome.splitChr10Exons)[[2]]

genome.splitChr6Exons <- diversity.stats(genome.splitChr6Exons)
divers6 <- get.diversity(genome.splitChr6Exons)[[2]]

genome.splitChr5Exons <- diversity.stats(genome.splitChr5Exons)
divers5 <- get.diversity(genome.splitChr5Exons)[[2]]

genome.splitChr7Exons <- diversity.stats(genome.splitChr7Exons)
divers7 <- get.diversity(genome.splitChr7Exons)[[1]]

genome.splitChr9Exons <- diversity.stats(genome.splitChr9Exons)
divers9 <- get.diversity(genome.splitChr9Exons)[[1]]

genome.splitChr12Exons <- diversity.stats(genome.splitChr12Exons)
divers12 <- get.diversity(genome.splitChr12Exons)[[1]]

genome.splitChr1Introns <- diversity.stats(genome.splitChr1Introns)
divers1Introns <- get.diversity(genome.splitChr1Introns)[[1]]

genome.splitChr6Introns <- diversity.stats(genome.splitChr6Introns)
divers6Introns <- get.diversity(genome.splitChr6Introns)[[2]]

genome.splitChr5Introns <- diversity.stats(genome.splitChr5Introns)
divers5Introns <- get.diversity(genome.splitChr5Introns)[[2]]

genome.splitChr7Introns <- diversity.stats(genome.splitChr7Introns)
divers7Introns <- get.diversity(genome.splitChr7Introns)[[1]]

genome.splitChr9Introns <- diversity.stats(genome.splitChr9Introns)
divers9Introns <- get.diversity(genome.splitChr9Introns)[[1]]

genome.splitChr12Introns <- diversity.stats(genome.splitChr12Introns)
divers12Introns <- get.diversity(genome.splitChr12Introns)[[1]]

#####Calculate Stats

DataFrameListChromDivers <- list() #List to hold dataframes for each chromosome

ListDiversityChrPops <- list()

#Nucleotide diversity must be divided by loci length; this function gets length for each loci
namedPositionIntervals12 <- sapply(rownames(divers12), function(x) eval(parse(text = x)))

positionIntervals12 <- namedPositionIntervals12

names(positionIntervals12) <- NULL

a <- genome.splitChr12Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr12Exons@nuc.diversity.within),i]/positionIntervals12,
    row.names = paste("Chr12",rownames(divers12), sep = " "))
  
}


Chr12DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr12DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[1]] <- Chr12DiversityDataFrame

ListDiversityChrPops <- list()

namedPositionIntervals11 <- sapply(rownames(divers11), function(x) eval(parse(text = x)))

positionIntervals11 <- namedPositionIntervals11

names(positionIntervals11) <- NULL

a <- genome.splitChr11Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr11Exons@nuc.diversity.within),i]/positionIntervals11,
    row.names = paste("Chr11",rownames(divers11), sep = " "))
  
}

Chr11DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr11DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[2]] <- Chr11DiversityDataFrame

ListDiversityChrPops <- list()

namedPositionIntervals10 <- sapply(rownames(divers10), function(x) eval(parse(text = x)))

positionIntervals10 <- namedPositionIntervals10

names(positionIntervals10) <- NULL

a <- genome.splitChr10Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr10Exons@nuc.diversity.within),i]/positionIntervals10,
    row.names = paste("Chr10",rownames(divers10), sep = " "))
  
}

Chr10DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr10DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[3]] <- Chr10DiversityDataFrame

ListDiversityChrPops <- list()

namedPositionIntervals1 <- sapply(rownames(divers1), function(x) eval(parse(text = x)))

positionIntervals1 <- namedPositionIntervals1

names(positionIntervals1) <- NULL

a <- genome.splitChr1Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr1Exons@nuc.diversity.within),i]/positionIntervals1,
    row.names = paste("Chr1",rownames(divers1), sep = " "))
  
}


Chr1DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr1DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[4]] <- Chr1DiversityDataFrame

ListDiversityChrPops <- list()

namedPositionIntervals5 <- sapply(rownames(divers5), function(x) eval(parse(text = x)))

positionIntervals5 <- namedPositionIntervals5

names(positionIntervals5) <- NULL

a <- genome.splitChr5Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr5Exons@nuc.diversity.within),i]/positionIntervals5,
    row.names = paste("Chr5",rownames(divers5), sep = " "))
  
}


Chr5DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr5DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[5]] <- Chr5DiversityDataFrame

ListDiversityChrPops <- list()
namedPositionIntervals6 <- sapply(rownames(divers6), function(x) eval(parse(text = x)))
positionIntervals6 <- namedPositionIntervals6
names(positionIntervals6) <- NULL

a <- genome.splitChr6Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr6Exons@nuc.diversity.within),i]/positionIntervals6,
    row.names = paste("Chr6",rownames(divers6), sep = " "))
  
}


Chr6DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr6DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[6]] <- Chr6DiversityDataFrame

ListDiversityChrPops <- list()
namedPositionIntervals7 <- sapply(rownames(divers7), function(x) eval(parse(text = x)))
positionIntervals7 <- namedPositionIntervals7
names(positionIntervals7) <- NULL

a <- genome.splitChr7Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr7Exons@nuc.diversity.within),i]/positionIntervals7,
    row.names = paste("Chr7",rownames(divers7), sep = " "))
  
}

Chr7DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr7DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[7]] <- Chr7DiversityDataFrame

ListDiversityChrPops <- list()
namedPositionIntervals8 <- sapply(rownames(divers8), function(x) eval(parse(text = x)))
positionIntervals8 <- namedPositionIntervals8
names(positionIntervals8) <- NULL

a <- genome.splitChr8Exons@nuc.diversity.within


for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr8Exons@nuc.diversity.within),i]/positionIntervals8,
    row.names = paste("Chr8",rownames(divers8), sep = " "))
  
}

Chr8DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr8DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[8]] <- Chr8DiversityDataFrame

ListDiversityChrPops <- list()
namedPositionIntervals9 <- sapply(rownames(divers9), function(x) eval(parse(text = x)))
positionIntervals9 <- namedPositionIntervals9
names(positionIntervals9) <- NULL

a <- genome.splitChr9Exons@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPops[[i]] <- data.frame(a[1:nrow(
    genome.splitChr9Exons@nuc.diversity.within),i]/positionIntervals9,
    row.names = paste("Chr9",rownames(divers9), sep = " "))
  
}

Chr9DiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPops)
colnames(Chr9DiversityDataFrame) <- c(1:7)
DataFrameListChromDivers[[9]] <- Chr9DiversityDataFrame

DiversityExonMaster <- dplyr::bind_rows(DataFrameListChromDivers)
colnames(DiversityExonMaster) <- c("4330", "4117A", "4107", "3111", "2932", "2931", "1963")
DiversityExonMaster <- DiversityExonMaster[!(row.names(DiversityExonMaster)
                                             %in% row.names.remove), ]
##Convert to Tidyverse compatible data


DiversityExonMaster <- tibble::rownames_to_column(DiversityExonMaster, var = "position" )
DiversityExonMaster <- DiversityExonMaster %>%
  pivot_longer(c("4330", "4117A", "4107", "3111", "2932", "2931", "1963"), 
               names_to = "population", values_to = "pi")

DiversityExonMaster <- mutate(DiversityExonMaster, pi = pi * -1)

#Plot

ExonNucleotideDiversityPlot <- ggplot(DiversityExonMaster, 
                                      aes(x = population, y = pi)) +
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")
print(ExonNucleotideDiversityPlot + ggtitle("Exon Nucleotide Diversity"))

#Same for Introns

DataFrameListChromDiversIntrons <- list() #List to hold dataframes for each chromosome

ListDiversityChrPopsIntrons <- list()

namedPositionIntervalsIntrons12 <- sapply(rownames(divers12Introns), function(x) eval(parse(text = x)))

positionIntervalsIntrons12 <- namedPositionIntervalsIntrons12

names(positionIntervalsIntrons12) <- NULL

a <- genome.splitChr12Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr12Introns@nuc.diversity.within),i]/positionIntervalsIntrons12,
    row.names = paste("Chr12",rownames(divers12Introns), sep = " "))
  
}


Chr12IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr12IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[1]] <- Chr12IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()

namedPositionIntervalsIntrons11 <- sapply(rownames(divers11Introns), function(x) eval(parse(text = x)))

positionIntervalsIntrons11 <- namedPositionIntervalsIntrons11

names(positionIntervalsIntrons11) <- NULL

a <- genome.splitChr11Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr11Introns@nuc.diversity.within),i]/positionIntervalsIntrons11,
    row.names = paste("Chr11",rownames(divers11Introns), sep = " "))
  
}

Chr11IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr11IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[2]] <- Chr11IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()

namedPositionIntervalsIntrons10 <- sapply(rownames(divers10Introns), function(x) eval(parse(text = x)))

positionIntervalsIntrons10 <- namedPositionIntervalsIntrons10

names(positionIntervalsIntrons10) <- NULL

a <- genome.splitChr10Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr10Introns@nuc.diversity.within),i]/positionIntervalsIntrons10,
    row.names = paste("Chr10",rownames(divers10Introns), sep = " "))
  
}

Chr10IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr10IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[3]] <- Chr10IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()

namedPositionIntervalsIntrons1 <- sapply(rownames(divers1Introns), function(x) eval(parse(text = x)))

positionIntervalsIntrons1 <- namedPositionIntervalsIntrons1

names(positionIntervalsIntrons1) <- NULL

a <- genome.splitChr1Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr1Introns@nuc.diversity.within),i]/positionIntervalsIntrons1,
    row.names = paste("Chr1",rownames(divers1Introns), sep = " "))
  
}


Chr1IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr1IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[4]] <- Chr1IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()

namedPositionIntervalsIntrons5 <- sapply(rownames(divers5Introns), function(x) eval(parse(text = x)))

positionIntervalsIntrons5 <- namedPositionIntervalsIntrons5

names(positionIntervalsIntrons5) <- NULL

a <- genome.splitChr5Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr5Introns@nuc.diversity.within),i]/positionIntervalsIntrons5,
    row.names = paste("Chr5",rownames(divers5Introns), sep = " "))
  
}


Chr5IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr5IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[5]] <- Chr5IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()
namedPositionIntervalsIntrons6 <- sapply(rownames(divers6Introns), function(x) eval(parse(text = x)))
positionIntervalsIntrons6 <- namedPositionIntervalsIntrons6
names(positionIntervalsIntrons6) <- NULL

a <- genome.splitChr6Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr6Introns@nuc.diversity.within),i]/positionIntervalsIntrons6,
    row.names = paste("Chr6",rownames(divers6Introns), sep = " "))
  
}


Chr6IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr6IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[6]] <- Chr6IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()
namedPositionIntervalsIntrons7 <- sapply(rownames(divers7Introns), function(x) eval(parse(text = x)))
positionIntervalsIntrons7 <- namedPositionIntervalsIntrons7
names(positionIntervalsIntrons7) <- NULL

a <- genome.splitChr7Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr7Introns@nuc.diversity.within),i]/positionIntervalsIntrons7,
    row.names = paste("Chr7",rownames(divers7Introns), sep = " "))
  
}

Chr7IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr7IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[7]] <- Chr7IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()
namedPositionIntervalsIntrons8 <- sapply(rownames(divers8Introns), function(x) eval(parse(text = x)))
positionIntervalsIntrons8 <- namedPositionIntervalsIntrons8
names(positionIntervalsIntrons8) <- NULL

a <- genome.splitChr8Introns@nuc.diversity.within


for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr8Introns@nuc.diversity.within),i]/positionIntervalsIntrons8,
    row.names = paste("Chr8",rownames(divers8Introns), sep = " "))
  
}

Chr8IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr8IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[8]] <- Chr8IntronDiversityDataFrame

ListDiversityChrPopsIntrons <- list()
namedPositionIntervalsIntrons9 <- sapply(rownames(divers9Introns), function(x) eval(parse(text = x)))
positionIntervalsIntrons9 <- namedPositionIntervalsIntrons9
names(positionIntervalsIntrons9) <- NULL

a <- genome.splitChr9Introns@nuc.diversity.within

for(i in 1:7){
  ListDiversityChrPopsIntrons[[i]] <- data.frame(a[1:nrow(
    genome.splitChr9Introns@nuc.diversity.within),i]/positionIntervalsIntrons9,
    row.names = paste("Chr9",rownames(divers9Introns), sep = " "))
  
}

Chr9IntronDiversityDataFrame <- dplyr::bind_cols(ListDiversityChrPopsIntrons)
colnames(Chr9IntronDiversityDataFrame) <- c(1:7)
DataFrameListChromDiversIntrons[[9]] <- Chr9IntronDiversityDataFrame

DiversityIntronMaster <- dplyr::bind_rows(DataFrameListChromDiversIntrons)
colnames(DiversityIntronMaster) <- c("4330", "4117A", "4107", "3111", "2932", "2931", "1963")
DiversityIntronMaster <- DiversityIntronMaster[!(row.names(DiversityIntronMaster)
                                             %in% row.names.remove.introns), ]
##Convert to Tidyverse compatible data


DiversityIntronMaster <- tibble::rownames_to_column(DiversityIntronMaster, var = "position" )
DiversityIntronMaster <- DiversityIntronMaster %>%
  pivot_longer(c("4330", "4117A", "4107", "3111", "2932", "2931", "1963"), 
               names_to = "population", values_to = "pi")

DiversityIntronMaster <- mutate(DiversityIntronMaster, pi = pi * -1)

#Plot

IntronNucleotideDiversityPlot <- ggplot(DiversityIntronMaster, 
                                        aes(x = population, y = pi)) +
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")
print(IntronNucleotideDiversityPlot + ggtitle("Intron Nucleotide Diversity"))


###################################################
#Fst
###################################################

genome.splitChr1Exons <- F_ST.stats(genome.splitChr1Exons)
genome.splitChr1Exons@nuc.F_ST.vs.all

genome.splitChr8Exons <- F_ST.stats(genome.splitChr8Exons)
F_ST8 <- get.F_ST(genome.splitChr8Exons)
genome.splitChr8Exons@nuc.F_ST.vs.all

genome.splitChr11Exons <- F_ST.stats(genome.splitChr11Exons)
F_ST11 <- get.F_ST(genome.splitChr11Exons)[[2]]
genome.splitChr11Exons@nuc.F_ST.vs.all

genome.splitChr10Exons <- F_ST.stats(genome.splitChr10Exons)
F_ST10 <- get.F_ST(genome.splitChr10Exons)[[1]]
genome.splitChr10Exons@nuc.F_ST.vs.all

genome.splitChr6Exons <- F_ST.stats(genome.splitChr6Exons)
F_ST6 <- get.F_ST(genome.splitChr6Exons)[[2]]
genome.splitChr6Exons@nuc.F_ST.vs.all

genome.splitChr5Exons <- F_ST.stats(genome.splitChr5Exons)
F_ST5 <- get.F_ST(genome.splitChr5Exons)[[2]]
genome.splitChr5Exons@nuc.F_ST.vs.all

genome.splitChr7Exons <- F_ST.stats(genome.splitChr7Exons)
F_ST7 <- get.F_ST(genome.splitChr7Exons)[[1]]
genome.splitChr7Exons@nuc.F_ST.vs.all

genome.splitChr9Exons <- F_ST.stats(genome.splitChr9Exons)
F_ST9 <- get.F_ST(genome.splitChr9Exons)[[2]]
genome.splitChr9Exons@nuc.F_ST.vs.all


genome.splitChr12Exons <- F_ST.stats(genome.splitChr12Exons)
F_ST12 <- get.F_ST(genome.splitChr12Exons)[[1]]

genome.splitChr1Introns <- F_ST.stats(genome.splitChr1Introns)
F_ST1Introns <- get.F_ST(genome.splitChr1Introns)[[1]]

genome.splitChr6Introns <- F_ST.stats(genome.splitChr6Introns)
F_ST6Introns <- get.F_ST(genome.splitChr6Introns)[[2]]

genome.splitChr5Introns <- F_ST.stats(genome.splitChr5Introns)
F_ST5Introns <- get.F_ST(genome.splitChr5Introns)[[2]]

genome.splitChr7Introns <- F_ST.stats(genome.splitChr7Introns, FAST = TRUE)
F_ST7Introns <- get.F_ST(genome.splitChr7Introns)[[1]]

genome.splitChr9Introns <- F_ST.stats(genome.splitChr9Introns)
F_ST9Introns <- get.F_ST(genome.splitChr9Introns)[[1]]

genome.splitChr12Introns <- F_ST.stats(genome.splitChr12Introns)
F_ST12Introns <- get.F_ST(genome.splitChr12Introns)[[1]]

DataFrameListChromosomes <- list() #List to hold dataframes for each chromosome

#Calculate Fst per chromosome per locus

ListFSTChr1Pops <- list()
for(i in 1:7){
  ListFSTChr1Pops[[i]] <- data.frame(genome.splitChr1Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr1Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr1",c(genome.splitChr1Exons@region.names), sep = " "))
  
}

Chr1FSTDataFrame <- dplyr::bind_cols(ListFSTChr1Pops)
colnames(Chr1FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[1]] <- Chr1FSTDataFrame

ListFSTChr12Pops <- list()
for(i in 1:7){
  ListFSTChr12Pops[[i]] <- data.frame(genome.splitChr12Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr12Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr12",c(genome.splitChr12Exons@region.names), sep = " "))
  
}

Chr12FSTDataFrame <- dplyr::bind_cols(ListFSTChr12Pops)
colnames(Chr12FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[2]] <- Chr12FSTDataFrame

ListFSTChr11Pops <- list()
for(i in 1:7){
  ListFSTChr11Pops[[i]] <- data.frame(genome.splitChr11Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr11Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr11",c(genome.splitChr11Exons@region.names), sep = " "))
  
}

Chr11FSTDataFrame <- dplyr::bind_cols(ListFSTChr11Pops)
colnames(Chr11FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[3]] <- Chr11FSTDataFrame

ListFSTChr10Pops <- list()
for(i in 1:7){
  ListFSTChr10Pops[[i]] <- data.frame(genome.splitChr10Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr10Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr10",c(genome.splitChr10Exons@region.names), sep = " "))
  
}

Chr10FSTDataFrame <- dplyr::bind_cols(ListFSTChr10Pops)
colnames(Chr10FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[4]] <- Chr10FSTDataFrame

ListFSTChr5Pops <- list()
for(i in 1:7){
  ListFSTChr5Pops[[i]] <- data.frame(genome.splitChr5Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr5Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr5",c(genome.splitChr5Exons@region.names), sep = " "))
  
}

Chr5FSTDataFrame <- dplyr::bind_cols(ListFSTChr5Pops)
colnames(Chr5FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[5]] <- Chr5FSTDataFrame

ListFSTChr6Pops <- list()
for(i in 1:7){
  ListFSTChr6Pops[[i]] <- data.frame(genome.splitChr6Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr6Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr6",c(genome.splitChr6Exons@region.names), sep = " "))
  
}

Chr6FSTDataFrame <- dplyr::bind_cols(ListFSTChr6Pops)
colnames(Chr6FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[6]] <- Chr6FSTDataFrame

ListFSTChr7Pops <- list()
for(i in 1:7){
  ListFSTChr7Pops[[i]] <- data.frame(genome.splitChr7Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr7Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr7",c(genome.splitChr7Exons@region.names), sep = " "))
  
}

Chr7FSTDataFrame <- dplyr::bind_cols(ListFSTChr7Pops)
colnames(Chr7FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[7]] <- Chr7FSTDataFrame

ListFSTChr8Pops <- list()
for(i in 1:7){
  ListFSTChr8Pops[[i]] <- data.frame(genome.splitChr8Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr8Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr8",c(genome.splitChr8Exons@region.names), sep = " "))
  
}

Chr8FSTDataFrame <- dplyr::bind_cols(ListFSTChr8Pops)
colnames(Chr8FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[8]] <- Chr8FSTDataFrame

ListFSTChr9Pops <- list()
for(i in 1:7){
  ListFSTChr9Pops[[i]] <- data.frame(genome.splitChr9Exons@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr9Exons@nuc.F_ST.vs.all),i],
    row.names = paste("Chr9",c(genome.splitChr9Exons@region.names), sep = " "))
  
}

Chr9FSTDataFrame <- dplyr::bind_cols(ListFSTChr9Pops)
colnames(Chr9FSTDataFrame) <- c(1:7)
DataFrameListChromosomes[[9]] <- Chr9FSTDataFrame

FSTExonMaster <- dplyr::bind_rows(DataFrameListChromosomes)
colnames(FSTExonMaster) <- c("4330", "4117A", "4107", "3111", "2932", "2931", "1963")
FSTExonMaster <- FSTExonMaster[!(row.names(FSTExonMaster)
                                             %in% row.names.remove), ]
##Convert to Tidyverse compatible data


FSTExonMaster <- tibble::rownames_to_column(FSTExonMaster, var = "position" )
FSTExonMaster <- FSTExonMaster %>%
  pivot_longer(c("4330", "4117A", "4107", "3111", "2932", "2931", "1963"), 
               names_to = "population", values_to = "FST")

#Plot

ExonNucleotideFSTPlot <- ggplot(FSTExonMaster, 
                                aes(x = population, y = FST)) +
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")
print(ExonNucleotideFSTPlot + ggtitle("Exon Nucleotide FST"))

#Same for introns

ListFSTChr1IntronPops <- list()
for(i in 1:7){
  ListFSTChr1IntronPops[[i]] <- data.frame(genome.splitChr1Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr1Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr1",c(genome.splitChr1Introns@region.names), sep = " "))
  
}

Chr1FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr1IntronPops)
colnames(Chr1FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[1]] <- Chr1FSTIntronsDataFrame

ListFSTChr12IntronPops <- list()
for(i in 1:7){
  ListFSTChr12IntronPops[[i]] <- data.frame(genome.splitChr12Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr12Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr12",c(genome.splitChr12Introns@region.names), sep = " "))
  
}

Chr12FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr12IntronPops)
colnames(Chr12FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[2]] <- Chr12FSTIntronsDataFrame

ListFSTChr11IntronPops <- list()
for(i in 1:7){
  ListFSTChr11IntronPops[[i]] <- data.frame(genome.splitChr11Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr11Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr11",c(genome.splitChr11Introns@region.names), sep = " "))
  
}

Chr11FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr11IntronPops)
colnames(Chr11FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[3]] <- Chr11FSTIntronsDataFrame

ListFSTChr10IntronPops <- list()
for(i in 1:7){
  ListFSTChr10IntronPops[[i]] <- data.frame(genome.splitChr10Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr10Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr10",c(genome.splitChr10Introns@region.names), sep = " "))
  
}

Chr10FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr10IntronPops)
colnames(Chr10FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[4]] <- Chr10FSTIntronsDataFrame

ListFSTChr5IntronPops <- list()
for(i in 1:7){
  ListFSTChr5IntronPops[[i]] <- data.frame(genome.splitChr5Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr5Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr5",c(genome.splitChr5Introns@region.names), sep = " "))
  
}

Chr5FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr5IntronPops)
colnames(Chr5FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[5]] <- Chr5FSTIntronsDataFrame

ListFSTChr6IntronPops <- list()
for(i in 1:7){
  ListFSTChr6IntronPops[[i]] <- data.frame(genome.splitChr6Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr6Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr6",c(genome.splitChr6Introns@region.names), sep = " "))
  
}

Chr6FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr6IntronPops)
colnames(Chr6FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[6]] <- Chr6FSTIntronsDataFrame

ListFSTChr7IntronPops <- list()
for(i in 1:7){
  ListFSTChr7IntronPops[[i]] <- data.frame(genome.splitChr7Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr7Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr7",c(genome.splitChr7Introns@region.names), sep = " "))
  
}

Chr7FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr7IntronPops)
colnames(Chr7FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[7]] <- Chr7FSTIntronsDataFrame

ListFSTChr8IntronPops <- list()
for(i in 1:7){
  ListFSTChr8IntronPops[[i]] <- data.frame(genome.splitChr8Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr8Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr8",c(genome.splitChr8Introns@region.names), sep = " "))
  
}

Chr8FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr8IntronPops)
colnames(Chr8FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[8]] <- Chr8FSTIntronsDataFrame

ListFSTChr9IntronPops <- list()
for(i in 1:7){
  ListFSTChr9IntronPops[[i]] <- data.frame(genome.splitChr9Introns@nuc.F_ST.vs.all[1:nrow(
    genome.splitChr9Introns@nuc.F_ST.vs.all),i],
    row.names = paste("Chr9",c(genome.splitChr9Introns@region.names), sep = " "))
  
}

Chr9FSTIntronsDataFrame <- dplyr::bind_cols(ListFSTChr9IntronPops)
colnames(Chr9FSTIntronsDataFrame) <- c(1:7)
DataFrameListChromosomes[[9]] <- Chr9FSTIntronsDataFrame

FSTIntronMaster <- dplyr::bind_rows(DataFrameListChromosomes)
colnames(FSTIntronMaster) <- c("4330", "4117A", "4107", "3111", "2932", "2931", "1963")
FSTIntronMaster <- FSTIntronMaster[!(row.names(FSTIntronMaster)
                                     %in% row.names.remove.introns), ]
##Convert to Tidyverse compatible data


FSTIntronMaster <- tibble::rownames_to_column(FSTIntronMaster, var = "position" )
FSTIntronMaster <- FSTIntronMaster %>%
  pivot_longer(c("4330", "4117A", "4107", "3111", "2932", "2931", "1963"), 
               names_to = "population", values_to = "FST")

#Plot

IntronNucleotideFSTPlot <- ggplot(FSTIntronMaster, 
                                  aes(x = population, y = FST)) +
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")
print(IntronNucleotideFSTPlot + ggtitle("Intron Nucleotide FST"))

#################################
#Pull interesting loci for analysis
#################################

sigFSTintrons <- subset(FSTIntronMaster, FSTIntronMaster$FST > .6)
conservedFSTintrons <- subset(FSTIntronMaster, FSTIntronMaster$FST < .1)

sigFSTexons<- subset(FSTExonMaster, FSTExonMaster$FST >.6)
conservedFSTexons <- subset(FSTExonMaster, FSTExonMaster$FST < .1)

sigTajimaExonsNeg <- subset(TajimaExonMaster, TajimaExonMaster$tajimasD < -1.5 )
sigTajimaExonsPos <- subset(TajimaExonMaster, TajimaExonMaster$tajimasD > 2)

sigTajimaIntronsNeg <- subset(TajimaIntronMaster, TajimaIntronMaster$tajimasD < -1.5)
sigTajimaIntronsPos <- subset(TajimaIntronMaster, TajimaIntronMaster$tajimasD > 2)


sigTajimaIntrons <- rbind(sigTajimaIntronsNeg, sigTajimaIntronsPos)
sigTajimaExons <- rbind(sigTajimaExonsNeg, sigTajimaExonsPos)

write.table(sigTajimaExons, file = "sigTajimaExons.txt", col.names = FALSE,
            row.names = FALSE)

write.table(sigTajimaIntrons, file = "sigTajimaIntrons.txt", col.names = FALSE, 
            row.names = FALSE)

write.table(sigFSTexons, file = "sigFSTExons.txt", col.names = FALSE, 
            row.names = FALSE)

write.table(sigFSTintrons, file = "sigFSTIntrons.txt", col.names = FALSE, 
            row.names = FALSE)

Exon_Pi_Tajima_DF <- full_join(TajimaExonMaster,DiversityExonMaster)
PiXTajimaExon <- ggplot(Exon_Pi_Tajima_DF,
                    aes(x = pi, y = tajimasD)) +
  geom_point()

Intron_Pi_Tajima_DF <- full_join(TajimaIntronMaster,DiversityIntronMaster)
PiXTajimaIntron <- ggplot(Intron_Pi_Tajima_DF,
                    aes(x = pi, y = tajimasD)) +
  geom_point()

ExonMeanTajima <- TajimaExonMaster %>% drop_na() %>%
  group_by(population) %>%
  summarise(tajimasD = mean(tajimasD))

IntronMeanTajima <- TajimaIntronMaster %>% drop_na() %>%
  group_by(population) %>%
  summarise(tajimasD = mean(tajimasD))

ExonMeanFST <- FSTExonMaster %>% drop_na() %>%
  group_by(population) %>%
  summarise(FST = mean(FST))

IntronMeanFST <- FSTIntronMaster %>% drop_na() %>%
  group_by(population) %>%
  summarise(FST = mean(FST))

ExonMeanDiversity <- DiversityExonMaster %>% drop_na() %>%
  group_by(population) %>%
  summarise(pi = mean(pi))

IntronMeanDiversity <- DiversityIntronMaster %>% drop_na() %>%
  group_by(population) %>%
  summarise(pi = mean(pi))


