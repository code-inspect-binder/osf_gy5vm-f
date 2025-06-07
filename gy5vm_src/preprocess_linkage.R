

data_l <- read.dta(paste(filepath_read_link,"\\HIV_BS_GE_ZH_Cancer.dta",sep=""))          # linkage data
data_l <- data_l[c("shcsid","registry","linked","cr_dinc","hiv_sex","cr_sex","hiv_dob_y","hiv_c_dob_y","hiv_dod_y","hiv_c_dod_y","hiv_dreg","cr_topo","cr_morpho","hiv_dstop")]

data_l$registry <- factor(data_l$registry)
data_l$cr_topo[data_l$cr_topo==""] <- NA
data_l$cr_topo <- factor(data_l$cr_topo)

# tumors from cancer registries
data_l$morpho <- NA
data_l$morpho[data_l$linked==1] <- "other/unknown"
data_l$morpho[data_l$cr_morpho%in%KSA_morpho] <- "KSA"
data_l$morpho[data_l$cr_morpho%in%NHL_morpho] <- "NHL"
data_l$morpho[data_l$cr_morpho%in%HL_morpho] <- "HL"
data_l$morpho[data_l$cr_morpho%in%leuk_morpho] <- "leuk"
data_l$topo <- NA
data_l$topo[data_l$linked==1] <- "other/unknown"
cr_topo_vect <- as.numeric(substr(data_l$cr_topo,start=2,stop=4))
data_l$topo[cr_topo_vect%in%BNS_topo] <- "BNS"
data_l$topo[cr_topo_vect%in%anus_topo] <- "anus"
data_l$topo[cr_topo_vect%in%cervix_topo] <- "cervix"
data_l$topo[cr_topo_vect%in%head_and_neck_topo] <- "head and neck"
data_l$topo[cr_topo_vect%in%liver_topo] <- "liver"
data_l$topo[cr_topo_vect%in%lung_topo] <- "lung"
data_l$topo[cr_topo_vect%in%prostate_topo] <- "prostate"
data_l$topo[cr_topo_vect%in%skin_topo] <- "skin"
data_l$topo[cr_topo_vect%in%bladder_topo] <- "bladder"
data_l$topo[cr_topo_vect%in%breast_topo] <- "breast"
data_l$topo[cr_topo_vect%in%colon_rectum_topo] <- "colon and rectum"
data_l$topo[cr_topo_vect%in%kidney_topo] <- "kidney"
data_l$topo[cr_topo_vect%in%oesophagus_topo] <- "oesophagus"
#data_l$topo[cr_topo_vect%in%soft_tissues_topo] <- "soft tissues"
data_l$topo[cr_topo_vect%in%stomach_topo] <- "stomach"
data_l$topo[cr_topo_vect%in%testis_topo] <- "testis"

data_l$skin_morpho <- NA
data_l$skin_morpho[data_l$topo=="skin" & data_l$cr_morpho%in%melanoma_morpho] <- "melanoma"           # ICD-0-3 definition
data_l$skin_morpho[data_l$topo=="skin" & !data_l$cr_morpho%in%melanoma_morpho] <- "non-melanoma"    
data_l$source <- "linkage"

cr_topo_vect2 <- substr(data_l$cr_topo,start=2,stop=3)
data_l$cr_nicer <- NA
data_l$cr_nicer[data_l$linked==1] <- "Other Sites"
data_l$cr_nicer[cr_topo_vect2=="21"] <- "Anus & Anal Canal"
data_l$cr_nicer[cr_topo_vect2=="67"] <- "Bladder"
data_l$cr_nicer[cr_topo_vect2%in%c("40","41")] <- "Bones, Joints, Cartilage"
data_l$cr_nicer[cr_topo_vect2%in%c("70","71","72")] <- "Brain & Central Nerves"
data_l$cr_nicer[cr_topo_vect2=="50"] <- "Breast"
data_l$cr_nicer[cr_topo_vect2=="53"] <- "Cervix Uteri"
data_l$cr_nicer[cr_topo_vect2%in%c("18","19","20")] <- "Colon, Rectum"
data_l$cr_nicer[cr_topo_vect2%in%c("54","55")] <- "Corpus Uteri & NOS"
data_l$cr_nicer[cr_topo_vect2=="69"] <- "Eye"
data_l$cr_nicer[cr_topo_vect2%in%c("23","24")] <- "Gallbladder & Extrahepatic Bile Tract"
data_l$cr_nicer[cr_topo_vect2=="64"] <- "Kidney"
data_l$cr_nicer[cr_topo_vect2=="32"] <- "Larynx"
data_l$cr_nicer[cr_topo_vect2=="22"] <- "Liver & Intrahepatic Bile Ducts"
data_l$cr_nicer[cr_topo_vect2%in%c("33","34")] <- "Lung, Bronchus, Trachea"
data_l$cr_nicer[cr_topo_vect2=="15"] <- "Oesophagus"
data_l$cr_nicer[as.numeric(cr_topo_vect2)%in%(0:14)] <- "Oral Cavity & Pharynx"
data_l$cr_nicer[cr_topo_vect2%in%c("65","66","68")] <- "Other Urinary Organs"
data_l$cr_nicer[cr_topo_vect2=="56"] <- "Ovary"
data_l$cr_nicer[cr_topo_vect2=="25"] <- "Pancreas"
data_l$cr_nicer[cr_topo_vect2=="61"] <- "Prostate"
data_l$cr_nicer[cr_topo_vect%in%c(384,450)] <- "Pleura"
data_l$cr_nicer[cr_topo_vect2=="44" & data_l$cr_morpho%in%melanoma_morpho] <- "Skin Melanoma"        # ICD-0-3 definition
data_l$cr_nicer[cr_topo_vect2=="17"] <- "Small Intestine"
data_l$cr_nicer[cr_topo_vect2%in%c("47","49")] <- "Soft Tissues"
data_l$cr_nicer[cr_topo_vect2=="16"] <- "Stomach"
data_l$cr_nicer[cr_topo_vect2=="62"] <- "Testis"
data_l$cr_nicer[cr_topo_vect2=="73"] <- "Thyroid"

# overrides topo classifications
data_l$cr_nicer[data_l$morpho=="HL"] <- "Hodgkin Lymphoma"
data_l$cr_nicer[data_l$morpho=="NHL"] <- "Non Hodgkin Lymphoma"
if(KSA_cat)
  data_l$cr_nicer[data_l$morpho=="KSA"] <- "Kaposi Sarcoma"

data_l1 <- subset(data_l,linked==1)
data_l0 <- subset(data_l,linked==0)
 
# tumors reported in SHCS
neoplasms <- c("NEO","NE2","NE3","NE4")
data_l1_SHCS <- with(data_dis,data.frame(shcsid=shcsid,disease=disease,cr_dinc=newdate,comments=comments))
data_l1_SHCS <- subset(data_l1_SHCS,shcsid%in%data_l$shcsid)                                                 #only considering patients for which linkage was attempted
data_l1_SHCS <- subset(data_l1_SHCS,disease%in%c(neoplasms,"KSA","NHL","MHO","ICC","LOB"))
#data_l1_SHCS <- subset(data_l1_SHCS,disease%in%c(neoplasms,"KSA","NHL","MHO","ICC","CDC","LOB"))
pat_to_dis <- match(data_l1_SHCS$shcsid,data_pat$shcsid)
admi_to_dis <- match(data_l1_SHCS$shcsid,data_admi$shcsid)
data_l1_SHCS$registry <- data_pat$registry[pat_to_dis]
data_l1_SHCS$linked <- 1
data_l1_SHCS$hiv_sex <- data_pat$sex[pat_to_dis]
data_l1_SHCS$hiv_sex <- factor(ifelse(data_l1_SHCS$hiv_sex==2,"Female","Male"),levels=levels(data_l$hiv_sex))
data_l1_SHCS$cr_sex <- NA
data_l1_SHCS$hiv_dob_y <- data_pat$born[pat_to_dis]
data_l1_SHCS$hiv_c_dob_y <- NA
data_l1_SHCS$hiv_dod_y <- as.Year(data_admi$exitdate[admi_to_dis])
data_l1_SHCS$hiv_c_dod_y <- NA
data_l1_SHCS$hiv_dreg <- data_pat$regdate[pat_to_dis]
data_l1_SHCS$cr_topo <- NA
data_l1_SHCS$cr_morpho <- NA
data_l1_SHCS$hiv_dstop <- data_admi$stopdate[admi_to_dis]
data_l1_SHCS$source <- "SHCS"
#data_l1_SHCS$skin_cat <- NA

# extracting site codes from the awful mess that is the comments variable
data_l1_SHCS$comments[data_l1_SHCS$comments==""] <- NA
data_l1_SHCS$comments_num <- substr(gsub(".","",data_l1_SHCS$comments,fixed=TRUE),start=2,stop=4)         
data_l1_SHCS$comments_num <- gsub("(","0",data_l1_SHCS$comments_num,fixed=TRUE)
data_l1_SHCS$comments_num <- gsub(" ","0",data_l1_SHCS$comments_num,fixed=TRUE)

data_l1_SHCS <- transform(data_l1_SHCS,comments_num=ifelse(nchar(comments_num)==2,paste(comments_num,"0",sep=''),comments_num))

#id_temp <- nchar(data_l1_SHCS$comments_num)==2 & !is.na(data_l1_SHCS$comments_num)
#data_l1_SHCS$comments_num[id_temp] <- paste(data_l1_SHCS$comments_num[id_temp],"0",sep='')
#rm(id_temp)
#data_l1_SHCS$comments_num <- suppressWarnings(as.numeric(data_l1_SHCS$comments_num))

data_l1_SHCS <- subset(data_l1_SHCS,registry%in%c("BS","GE","ZH"))             # restricting to Basel, Geneva, Zurich
data_l1_SHCS$registry <- factor(data_l1_SHCS$registry,levels=levels(data_l$registry))
data_l1_SHCS$morpho <- NA
data_l1_SHCS$morpho[data_l1_SHCS$disease=="KSA"] <- "KSA"
data_l1_SHCS$morpho[data_l1_SHCS$disease=="NHL"] <- "NHL"                     # some of these may be PBL, there's no way to know
data_l1_SHCS$morpho[data_l1_SHCS$disease=="MHO"] <- "HL"                    
data_l1_SHCS$topo[data_l1_SHCS$disease=="ICC"] <- "cervix"               # leaving out disease=CDC (cervical carcinoma)
data_l1_SHCS$topo[data_l1_SHCS$disease=="LOB"] <- "BNS"                      
data_l1_SHCS$morpho[data_l1_SHCS$disease=="LOB"] <- "NHL"                    # (questionable) assumption here: the PBL are all NHL
data_l1_SHCS$morpho[data_l1_SHCS$disease%in%c(neoplasms,"ICC")] <- "other/unknown"

#data_l1_SHCS$morpho[data_l1_SHCS$disease%in%c(neoplasms,"ICC","CDC")] <- "other/unknown"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%c(neoplasms,"KSA","NHL","MHO")] <- "other/unknown"
data_l1_SHCS$morpho[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%leuk_morpho] <- "leuk"             # leukaemia=morphology
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%anus_topo] <- "anus"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%head_and_neck_topo] <- "head and neck"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%liver_topo] <- "liver"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%lung_topo] <- "lung"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%prostate_topo] <- "prostate"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%skin_topo] <- "skin"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%BNS_topo] <- "BNS"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%bladder_topo] <- "bladder"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%breast_topo] <- "breast"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%colon_rectum_topo] <- "colon and rectum"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%kidney_topo] <- "kidney"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%oesophagus_topo] <- "oesophagus"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%stomach_topo] <- "stomach"
data_l1_SHCS$topo[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%testis_topo] <- "testis"
data_l1_SHCS$skin_morpho <- NA
data_l1_SHCS$skin_morpho[data_l1_SHCS$comments_num%in%430:439] <- "melanoma"              # ICD-10 definition
data_l1_SHCS$skin_morpho[data_l1_SHCS$comments_num%in%440:449] <- "non-melanoma"

data_l1_SHCS$cr_nicer <- NA
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms]<- "Other Sites"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(210:219)] <- "Anus & Anal Canal"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(670:679)] <- "Bladder"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(400:419)] <- "Bones, Joints, Cartilage"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(700:729)] <- "Brain & Central Nerves"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(500:509)] <- "Breast"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(180:209)] <- "Colon, Rectum"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(540:559)] <- "Corpus Uteri & NOS"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(690:699)] <- "Eye"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(230:249)] <- "Gallbladder & Extrahepatic Bile Tract"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(640:649)] <- "Kidney"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(320:329)] <- "Larynx"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(220:229)] <- "Liver & Intrahepatic Bile Ducts"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(330:349)] <- "Lung, Bronchus, Trachea"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(150:159)] <- "Oesophagus"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(0:149)] <- "Oral Cavity & Pharynx"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%c(650:669,680:689)] <- "Other Urinary Organs"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(560:569)] <- "Ovary"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(250:259)] <- "Pancreas"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%c(384,450)] <- "Pleura"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(610:619)] <- "Prostate"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(430:439)] <- "Skin Melanoma"            # C44 is non-melanoma (ICD-10)
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(170:179)] <- "Small Intestine"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%c(470:479,490:499)] <- "Soft Tissues"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(160:169)] <- "Stomach"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(620:629)] <- "Testis"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease%in%neoplasms & data_l1_SHCS$comments_num%in%(730:739)] <- "Thyroid"

data_l1_SHCS$cr_nicer[data_l1_SHCS$disease=="KSA"] <- "Kaposi Sarcoma"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease=="MHO"] <- "Hodgkin Lymphoma"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease=="NHL"] <- "Non Hodgkin Lymphoma"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease=="ICC"] <- "Cervix Uteri"
data_l1_SHCS$cr_nicer[data_l1_SHCS$disease=="LOB"] <- "Other Sites"               # unclear where to assign this

if(impute_KSA)
{
  load(file=paste(filepath_read_R,"\\KSA_imputations.RData",sep=''))
  data_l1_SHCS$cr_nicer[match(SHCS_imp$shcsid,data_l1_SHCS$shcsid)] <- SHCS_imp$cr_nicer
}

if(!nicer_classification)
{data_l1_SHCS <- subset(data_l1_SHCS,!is.na(morpho) | !is.na(topo))} else
{data_l1_SHCS <- subset(data_l1_SHCS,!is.na(cr_nicer))}

data_l1 <- switch(info_source,"both"=rbind(data_l1,data_l1_SHCS[names(data_l1)]),"linkage"=data_l1,"SHCS"=data_l1_SHCS[names(data_l1)])
data_l0 <- subset(data_l0,!shcsid%in%data_l1$shcsid)                             # removing unlinked patients that also appear as linked

# breaking cancers down into chosen categories: PBL -> morpho -> topo
data_l1$topo[data_l1$topo=="skin" & data_l1$skin_morpho=="melanoma"] <- "skin melanoma"
data_l1$topo[data_l1$topo=="skin"] <- "skin non-melanoma"
data_l1$topo[!data_l1$topo%in%c("BNS","cervix",NADC_vect)] <- "other NADCs"            # category for misc NADCs
data_l1$cr_cat <- with(data_l1,ifelse(topo=="BNS" & morpho=="NHL","PBL",ifelse(!morpho%in%c("leuk","other/unknown"),morpho,topo)))    # leukaemia classified as 'other NADC'
if(!"HL"%in%NADC_vect)
  data_l1$cr_cat[data_l1$cr_cat=="HL"] <- "other NADCs"
#data_l1$topo[data_l1$topo=="other NADCs"] <- "other"
data_l0$cr_cat <- NA
data_l1_res <- data_l1
data_l0_res <- data_l0

if(nicer_classification)
 data_l1$cr_cat <- data_l1$cr_nicer

data_l1 <- data_l1[order(data_l1$source,data_l1$cr_dinc),]
 
if(!(tumor_restriction%in%c("none","any")))
{
 data_l1_res <- subset(data_l1,cr_cat==tumor_restriction)
 data_l0_res <- rbind(data_l0,subset(data_l1,cr_cat!=tumor_restriction & !(shcsid%in%data_l1_res$shcsid)))   # unlinking patients that only have other types of tumor
} else if(tumor_restriction=="none")
{
 data_l1_res <- data_l1
}

data_l0_res <- data_l0_res[sample(1:nrow(data_l0_res),replace=FALSE),]            # randomizing order of unlinked patients
data_l0_res <- data_l0_res[!duplicated(data_l0_res$shcsid),]                      # removing patients appearing several times (random choice)
data_l0_res$linked <- 0

data_l1_res <- data_l1_res[!duplicated(data_l1_res$shcsid),]                               # only keeping earliest reported tumor of given type, and linkage info preferred to SHCS info

# breakpoints: cancer incidence, start of haart
X1 <- subset(data_l,linked==1,select=c("shcsid","cr_dinc"))
names(X1)[2] <- "date"
X2 <- subset(data_dis,shcsid%in%data_l$shcsid & disease%in%c(neoplasms,"KSA","NHL","MHO","ICC","LOB"),select=c("shcsid","newdate"))
names(X2)[2] <- "date"
X3 <- subset(data_tail,!is.na(haart_start_date),select=c("shcsid","haart_start_date"))
names(X3)[2] <- "date"
X <- rbind(X1,X2,X3)
X <- X[order(X$date),]
X <- subset(X,!duplicated(paste(X$shcsid,X$date)))
breakpoints <- split(X$date,X$shcsid)
rm(X,X1,X2,X3)
