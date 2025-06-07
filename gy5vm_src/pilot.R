
library(foreign)
library(survival)
library(nnet)                  # for multinomial regression
library(plyr)                  # for fast reassembling of data frames

ptm <- proc.time()

set.seed(10501)

filepath_code <- "C:\\HomeDir\\P3RL\\R_code"
filepath_read_link <- "E:\\P3RL\\linkage"
filepath_read_SHCS <- "E:\\P3RL\\SHCS"
filepath_read_NICER <- "E:\\P3RL\\NICER"
filepath_read_R <- "E:\\P3RL\\R"

tumor_restriction <- "anus"                # "none"=keeping all patients regardless of cancer incidence, "any"=single category for all cancers, also possible: "other NADCs", "skin melanoma", "skin non-melanoma"
info_source <- "both"                       # which source of cancer information to use: "linkage" (cancer registries), "SHCS", or "both"
multiplier <- 1000
round_digits <- 2

start_point <- "reg"                 # start of risk period, either HIV diagnostic ("diag") or registration ("reg")
lower_cutoff_year <- 1980            # left-truncation threshold (should be 1988 or 1996 when computing SIRs, otherwise 1980 or before)
nicer_classification <- FALSE        # whether to label tumors based on NICER's categories (e.g. for SIRs), otherwise using inhouse classification
year_restriction <- TRUE             # whether to restrict tumors to those occuring after entry into the registry and within years spanned by cancer registries - possibly turn off if counting linkage cases
impute_KSA <- FALSE                  # whether to impute NICER classification of KSA cases reported in SHCS but not in cancer registries
KSA_cat <- TRUE                      # whether to have a category for KSA in NICER classification (necessary to remove KSA from other categories)

if(!year_restriction)
  warning("no year restriction!")

#NADC_vect <- c()
NADC_vect <- c("anus","lung","prostate","liver")    # any NADC not in this list will be classified as 'other NADCs' (not applicable to NICER classification currently), give empty vector for analyzing all NADCs
#NADC_vect <- c("anus","lung","head and neck","prostate","liver","HL","bladder","breast","colon and rectum","kidney","oesophagus","skin melanoma","stomach","testis")

nacc_topo <- read.csv(file=paste(filepath_read_SHCS,"\\csv\\NA-ACCORD_mapping.csv",sep=''),sep=",", stringsAsFactors=FALSE)[-c(335,336),]
colnames(nacc_topo) <- c("V1","V2","V3","V4")
nacc_skin_morpho <- read.csv(file=paste(filepath_read_link,"\\skin_cancers_morpho.csv",sep=''),sep=";", stringsAsFactors=FALSE,header=FALSE)
leuk_morpho <- read.csv(file=paste(filepath_read_link,"\\leukaemia_morpho.csv",sep=''),sep=";", stringsAsFactors=FALSE,header=FALSE)$V1

KSA_morpho <- c(9140)
NHL_morpho <- c(9705,9695,9691,9698,9690,9716,9596,9590,9709,9714,9687,9717,9680,9684,9671,9673,9699,9689,9702,9679,9675,9700,9719,9591,9728,9727,9729,9701,9670,9708,9678,9827,9823,9761)
HL_morpho <- c(9661,9654,9659,9655,9651,9653,9652,9664,9665,9667,9663,9650,9662)
BNS_topo <- c(470:476,478,700,701,709:725,728,729)                   # brain and nervous system
anus_topo <- c(210,211,212)
cervix_topo <- c(530,531,538,539)
head_and_neck_topo <- as.numeric(subset(nacc_topo,V4%in%c(22,35,40,42,44,45,46,47,48,50))$V1) # some numbers have 1 or 2 digits
liver_topo <- c(220)
lung_topo <- c(340:343,348,349)
prostate_topo <- c(610,619)          # shows up as "C61" in SHCS dataset
#skin_topo <- c(440:449,760,764,765)           
skin_topo <- c(430:449,760,764,765)      # warning: differences between ICD-0-3 (i.e. used by cancer registries) and ICD-10 (used by SHCS and NICER)
bladder_topo <- 670:679
breast_topo <- 500:509
colon_rectum_topo <- 180:209
kidney_topo <- 640:649
oesophagus_topo <- 150:159
soft_tissues_topo <- c(470:479,490:499)
stomach_topo <- 160:169
testis_topo <- 620:629

#skin cancer subtypes
melanoma_morpho <- subset(nacc_skin_morpho,V4==64)$V1          # used with ICD-0-3
nonmelanoma_morpho <- subset(nacc_skin_morpho,V4==65)$V1

africa_region_id <- c(17,11,18,14,15)            # all of Africa
#africa_region_id <- c(17,11,18,14)              # Sub-Saharan Africa
europe_region_id <- c(154,39,155,151)
other_region_id <- c(34,61,29,5,145,53,13,21,35,30,54,57,143)
if(length(unique(union(union(africa_region_id,europe_region_id),other_region_id)))!=22)
 print("problem with region IDs")

as.Year <- function(date)
 as.numeric(substr(date,1,4))

as.Month <- function(date)
 as.numeric(substr(date,6,7))

as.YearMonth <- function(date)
 substr(date,1,7)

month_delta <- function(d2,d1)                                      # number of months seperating two dates
 12*(as.Year(d2)-as.Year(d1))+as.Month(d2)-as.Month(d1)

source(paste(filepath_code,"\\preprocess_SHCS.R",sep=''))
source(paste(filepath_code,"\\preprocess_linkage.R",sep=''))
source(paste(filepath_code,"\\merge.R",sep=''))

print(paste("tumor:",tumor_restriction))
print(paste("source of cancer info:",info_source))
print(proc.time()-ptm)
