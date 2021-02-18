#!/usr/bin/env Rscript

#title: "MEP-LINCS Preprocessing"
#author: "Mark Dane"
library(MEMA)#merge, annotate and normalize functions
library(tidyverse)
library(RUVnormalize)
library(ruv)
library(rrscale)
suppressPackageStartupMessages(library(optparse))
#library(synapser)


# Get the command line arguments and options
# returns a list with options and args elements
getCommandLineArgs <- function(){
  parser <- OptionParser(usage = "%prog [options] INPUT_PATH OUTPUT_PATH STUDY_NAME K")
  arguments <- parse_args(parser, positional_arguments = 4)
}

###update MEMA function to synapser package
getBarcodes <- function(studyName, synId = "syn10846457"){
  synLogin()
  synObj <- synGet(synId)
  barcodes <- read_tsv(synObj$path, col_types = cols(
    StudyName = col_character(),
    Barcode = col_character()
  )) %>%
    filter(StudyName==str_to_lower(studyName)) %>%
    select(Barcode) %>%
    stringr::str_split(",") %>%
    unlist()
  return(barcodes)
}

#Debug add signif reduction
writeLevelStudyData <- function(x, level){
  # if(!dir.exists(paste0(output_path,"study/",studyName))) dir.create(paste0(output_path,"study/",studyName), recursive = TRUE)
  # if(!dir.exists(paste0(output_path,"study/",studyName,"/Annotated"))) dir.create(paste0(output_path,"study/",studyName,"/Annotated"), recursive = TRUE)
  if(!dir.exists(paste0(output_path,"Annotated/GT1"))) dir.create(paste0(output_path,"Annotated/GT1"), recursive = TRUE)
  write_csv(x, paste0(output_path,studyName,"_level_",level,".csv"))
}

###End of functions

#Specify the command line options
if(!interactive()){
  cl <- getCommandLineArgs()
  input_path <- cl$args[1]
  output_path <- cl$args[2]
  studyName <- cl$args[3]
  k <- cl$args[4]
} else {
  input_path <- "~/Documents/A549_low_serum/data/"
  output_path <- "~/Documents/A549_low_serum/data/"
  studyName <- "a549_low_serum"
  k <- 0
}

#get level 2 files
#Get barcodes based on Synapse data
barcodes <-c(paste0("LI8X0113", 2:9), "LI8X01140")

level2_data <- lapply(barcodes, function(barcode){
  sd <- read_csv(paste0(input_path,barcode,"_level_2.csv"),
                 col_types = cols())
}) %>%
  bind_rows() %>%
  filter(!ECMp=="None") %>%
  select(-matches("Cells_")) %>%
  mutate(BW = paste(Barcode, Well, sep = "_")) %>%
  data.table()
colnames(level2_data) <- make.names(colnames(level2_data))

#Hack in change of no vehicle in some control wells
level2_data$Drug[level2_data$Drug == "air"] <- "DMSO"

rrscaleColumn <- function(x, zeros = 1){
  xRRValues <- as.matrix(x) %>%
    rrscale(zeros = zeros) 
  x_RR <- xRRValues$RR %>%
    pull
  return(x_RR)
}

level2_data <- level2_data %>%
  mutate(across(matches("SpotCellCount|Intensity"), rrscaleColumn, .names = "{col}_RR"))

#RUVLoess normalize all signals
if(!k==0){
  message(paste("Normalizing", studyName,"\n"))
  signalsMinimalMetadata <- level2_data %>%
    select(matches("SpotCellCount|_RR$|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^Drug1Conc$|^ArrayRow$|^ArrayColumn$|^CellLine$"),
           -matches("_SE$")) %>%
    data.table()
  nDT <- norm_RR_RUVResiduals(signalsMinimalMetadata, k)
  nDT$NormMethod <- "RR_RUV_Residuals"
  level2_data$k <- k
  level3_data <- merge(level2_data, nDT, by = c("BW","PrintSpot"))
  rm(level2_data)
} else {
  level3_data <- level2_data %>%
    mutate(NormMethod = "none",
           k = k)
}

#Add QA flags to the data
level3_data <- QASpotLevelData(level3_data, lowSpotCellCountThreshold=5,
                        lowRegionCellCountThreshold = 0.4,
                        lowWellQAThreshold = .7)

writeLevelStudyData(level3_data, level = 3)

level4_numeric_data <- level3_data %>%
  select(-contains("Spot")) %>%
  select(-matches("xi|ImageID")) %>%
  group_by(MEP_Drug) %>%
  summarise_if(is.numeric, median)  %>%
  ungroup() %>%
  mutate_if(is.numeric, signif, digits=4)

level4_character_data <- level3_data %>%
  select(-matches("ImageID|xi|recordSet|BW|Barcode|PlateID|Well|ID|ConcUnit")) %>%
  group_by(MEP_Drug) %>%
  summarise_if(is.character, unique) %>%
  ungroup()

level4_logical_data <- level3_data %>%
  select(-matches("QA_LowSpotCellCount")) %>%
  group_by(MEP_Drug) %>%
  summarise_if(is.logical, unique) %>%
  ungroup()

level4_data <- full_join(level4_character_data, level4_numeric_data, by = c("MEP_Drug")) %>%
  right_join(level4_logical_data, by = c("MEP_Drug")) %>%
  mutate(QA_LowReplicateCount = Spot_PA_ReplicateCount < 3)

#Add in the barcodes for each MEP_Drug
#level4_data <- addBarcodes(dt3 = level3_data, dt4 = level4_data)

writeLevelStudyData(level4_data, 4)
