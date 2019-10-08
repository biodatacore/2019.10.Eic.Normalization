library(matrixStats) #install.packages('matrixStats')
library(R.utils) #install.packages('R.utils')
library(utils)
library(plyr) #install.packages('plyr')
library(reshape2) #install.packages('reshape2')
library(dplyr) #install.packages('dplyr')
library(tidyr) #install.packages('tidyr')
library(methods)
library(data.table) #install.packages('data.table')
library(readxl) #install.packages('readxl')
library(stringr)
library(magrittr)

#source("https://bioconductor.org/biocLite.R")
#biocLite("vsn")
library(vsn)


# Import Data -------------------------------------------------------------

#setwd('~/Dropbox (Partners HealthCare)/Susan-Mir Shared/Eicosanoids-FS-CODE/eicdata/FINRISK')
#setwd("~/Knowledge/SHC")

mainDir <- '~/Dropbox (Partners HealthCare)/Susan-Mir Shared/Eicosanoids-FS-CODE/eicdata/FINRISK'
eicdatafile <- '171025_FR02_Analysis_Parsed.xlsx'
is_file <- '2018.09.03_On-Plate_PooledPlasma.csv'
clinfile <- 'FR02.xlsx'
runorderfile <- '180223_FR02_Sample_Key_Run_Order_AK_edit_wplate.xlsx'


fr02 <- readxl::read_xlsx(file.path(mainDir, eicdatafile), sheet = 'Eicosanoids_Docosanoids') %>% 
  as.data.frame() %>% 
  subset(select = -c(Charge,ID,Adducts)) %>% 
  mutate(`m/z` = formatC(`m/z`, digits = 5, format = "f")) %>% 
  rename(mz = `m/z`) %>% 
  mutate(RT = formatC(RT, digits = 4, format = "f")) %>% 
  select(sort(names(.))) %>% 
  mutate(mzid = paste0(mz, '_', RT)) %>%
  select(mzid, everything()) %>% 
  arrange(mzid) %>% 
  subset(select = -c(mz, RT)) 
is <- read.csv(file.path(mainDir, is_file)) %>% 
  as.data.frame() %>% 
  subset(select = -c(featID)) %>% 
  mutate(mz = formatC(mz, digits = 5, format = "f")) %>% 
  mutate(RT = formatC(RT, digits = 4, format = "f")) %>% 
  select(sort(names(.))) %>% 
  mutate(mzid = paste0(mz, '_', RT)) %>%
  select(mzid, everything()) %>% 
  arrange(mzid) %>% 
  subset(select = -c(mz, RT)) %>% 
  `names<-`(gsub('[X]', '', names(.))) %>% 
  mutate_all(funs(replace(., is.na(.), 0)))
clin <- readxl::read_xlsx(file.path(mainDir, clinfile)) %>% 
  as.data.frame() %>% 
  subset(select = c(Sample_ID, MEN))
idmatch <- readxl::read_xlsx(file.path(mainDir, runorderfile)) %>% 
  as.data.frame() %>% 
  mutate(Plate = sub('.*_', '', `Matrix plate name`)) %>%
  mutate(well_corrected = str_pad(well, 2, pad = "0")) %>% 
  mutate(Plate_well = paste0(Plate, "_", well_corrected)) %>% 
  subset(select = c(Sample_ID, Plate, Plate_well, EIC_ID)) %>% 
  merge(clin, by = 'Sample_ID')

rm(clin)

# Need to fix certain ID names in the internal standards df

#names_to_keep <- setdiff(names(fr02),names(is))
#ind <- match(setdiff(names(is),names(fr02)), names(is))
#names(is)[ind] <- names_to_keep

#is <- is[names(fr02)] # Reorder IS columns to match FR02 cols

# Need to account here for the Relo IDs

#names(fr02) <- idmatch$Plate_well[match(names(fr02),idmatch$EIC_ID)]

names(is) <- idmatch$EIC_ID[match(names(is),idmatch$Plate_well)] # Change to Relo values
names(is)[1] <- 'mzid' # Rename first column

is <- is[1:(ncol(is)-2)] # Last two columns don't have match ids yet

# Format plate_well in pw correctly

# add col names and turn into row names

fr02 <- data.frame(fr02[,-1], row.names = fr02[,1])
is <- data.frame(is[,-1], row.names = is[,1])

# platenum <- gsub('_.*$', '', pw$Plate_Position) %>% as.numeric()
# platenum <- sprintf("%02d", platenum)
# 
# wellnum <- gsub('.*_', '', pw$Plate_Position) %>% as.numeric()
# wellnum <- sprintf("%02d", wellnum)
# 
# pw$Plate_Position <- paste0(platenum, "_", wellnum)

# # Eicosanoids
# mzid_info <- fhs[,c(1:4)]
# 
# fhs <- select(fhs, one_of(intersect(a, b)))

# Remove X from colnames

names(fr02) <- gsub('[X]', '', names(fr02))

batch <- sub("_[^_]+$", "", names(fr02))
# Gender

idmatch <- idmatch[order(idmatch$EIC_ID), ]

#fr02 <- fr02[idmatch$EIC_ID]
#is <- is[idmatch$EIC_ID]

# Checking that the names are all sorted correctly (should be numeric)
all.equal(names(fr02), idmatch$Plate_well)
all.equal(names(is), idmatch$Plate_well)

gr <- idmatch$MEN

# FINRISK also has no variables with all zeros, so ignore the following code
#allzero <- which(apply(is, 1, function(x) all(x == 0))) # remove variables with all zero values
#fr02 <- fr02[-allzero, ]


fr02_clean <- fr02[-which(rowMeans(isZero(fr02)) > 0.25), ]
is_clean <- is[-which(rowMeans(isZero(is)) > 0.05), ]

# Replace 0s with 1/4 * minimum
fr02_mins <- apply(fr02_clean, 1, function(x) min(x[x>0]*0.25))
is_mins <- apply(is_clean, 1, function(x) min(x[x>0]*0.25))

replace_with_min <- function(x) {
  ifelse(x==0,fr02_mins,x)
}


# tmp <- apply(test, 1, function(x) ifelse(x==0,fr02_mins,x)) %>% t() %>% as.data.frame()
# 
# tmp <- apply(test, 1, function(x) ifelse(x==0,fr02_mins,x)) %>% 
#   log() %>% 
#   cbind(batch[1:20]) %>% 
#   as.data.frame()
# 
# tmp <- setDT(tmp)[, lapply(.SD, function(x) (x - mean(x)))/sd(x), by = batch] # Subtract batch means

dat <- as.data.frame(fr02)

# InternalStandards
# is <- is[-16, ] # remove the special one: fpHILIC
# is_mzids <- is[, 3]
# is <- select(is, one_of(intersect(a, b)))
# is <- cbind(is_mzids, is)
# is <- is[order(is$is_mzids), ]
# is$is_mzids <- paste0("mzid_", is$is_mzids) # add 'mzid_' to numbers
# row.names(is) <- is[, 1]
# is <- is[, -1]
is <- as.data.frame(t(is))

# combine Eicosanoids and is
fr02t <- cbind(is, fr02)

# group (factor of interest)
# gr[which(gr == "Unknown")] <- "M" # there is one sample with "Unknown", just value it with "M".
# gr <- as.character(gr)
# grp <- gr
# gr[gr == "M"] <- 1
# gr[gr == "F"] <- 0
# gr <- as.matrix(as.numeric(gr))


###################################################################################################
## Normalization Methods                                                                         ##
###################################################################################################

# MeanCenter ####
centr.a <- MeanCenter(dat, batch)

# MedianCenter ####
centr.b <- MedianCenter(dat, batch)

# Quantile ####
centr.q <- Quantile(dat)

# ComBat ####
centr.cb <- combat(dat, batch)

# Quantile + ComBat ####
centr.qcb <- combat(as.data.frame(t(centr.q)), batch)


# RUV ####
ctl <- 1:ncol(is)
k <- 5

## ruv2
fit2 <- RUV2(fr02t, gr, ctl, k)
centr.ruv2 <- as.data.frame(t(dat - (fit2$W)%*%(fit2$alpha))) # corrected = Y - W*alpha

## ruv4
#k <- getK(shct, gr, ctl) # get k of unwanted variation factors
fit4 <- RUV4(fr02t, gr, ctl, k)
centr.ruv4 <- as.data.frame(t(dat - (fit4$W)%*%(fit4$alpha)))

# EigenMS ####
prot.info <- cbind(data.frame(colnames(dat)), data.frame(colnames(dat)))
rv <- eig_norm1(m = as.data.frame(t(dat)), treatment = as.factor(batch), prot.info)
rv2 <- eig_norm2(rv)
centr.EigenMS <- as.data.frame(rv2$norm_m)

## Working on getting the bottom 3 methods to work atm

# # PQN ####
# datx <- as.data.frame(t(dat))
# datx$name <- row.names(datx)
# 
# datp <- new("metaXpara") # create new object of 'metaCpara'
# sfile <- "samplelist.txt"
# rawPeaks(datp) <- datx # value object with actual dataset
# sampleListFile(datp) <- sfile # make sure the sample column in samplelist.txt is exactly same as rownames of datx, samplelist.txt is a tab delimited txt file.
# datp <- reSetPeaksData(datp)
# 
# centr.pqn <- normalize(datp, method = "pqn", valueID="value")
# centr.pqn <- getPeaksTable(centr.pqn) # get normalized dataset 
# row.names(centr.pqn) <- centr.pqn[, 1]
# centr.pqn <- as.data.frame(t(centr.pqn[, -c(1:4)]))
# 
# # SVR ####
# centr.svr <- normalize(datp, method = "svr", valueID="value") 
# centr.svr <- getPeaksTable(centr.svr)
# row.names(centr.svr) <- centr.svr[, 1]
# centr.svr <- as.data.frame(t(centr.svr[, -c(1:4)]))
# 
# # VSN ####
# centr.vsn <- normalize(datp, method = "vsn", valueID="value")
# centr.vsn <- getPeaksTable(centr.vsn)
# row.names(centr.vsn) <- centr.vsn[, 1]
# centr.vsn <- as.data.frame(t(centr.vsn[, -c(1:4)]))

outputDir <- '~/Dropbox (Partners HealthCare)/2018 Applied Bioinformatics Work/Personal folders/Andy Kim/NormalizationYY/'
subDir <- 'FINRISK_normalization'
dir.create(file.path(outputDir, subDir), showWarnings = FALSE)

write.csv(dat, 'FINRISK_normalization/FR02_eic_data.csv')
write.csv(centr.a, 'NormalizationAK/FR02_centr_a.csv')
write.csv(centr.b, 'NormalizationAK/FR02_centr_b.csv')
write.csv(centr.q, 'NormalizationAK/FR02_centr_q.csv')
write.csv(centr.cb, 'NormalizationAK/FR02_centr_cb.csv')
write.csv(centr.qcb, 'NormalizationAK/FR02_centr_qcb.csv')
write.csv(centr.ruv2, 'NormalizationAK/FR02_centr_ruv2.csv')
write.csv(centr.ruv4, 'NormalizationAK/FR02_centr_ruv4.csv')
write.csv(centr.EigenMS, 'NormalizationAK/FR02_centr_EigenMS.csv')

rm(centr.a)
rm(centr.b)
rm(centr.cb)
rm(centr.q)
rm(centr.qcb)
rm(centr.EigenMS)
rm(centr.ruv2)
rm(centr.ruv4)


# visualizing finrisk data completeness -----------------------------------

hist(nrow(fr02) - colSums(isZero(fr02[,-1])), main = 'Number of Present Metabolites')
