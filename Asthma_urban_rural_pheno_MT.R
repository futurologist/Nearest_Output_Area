library(dplyr)
library(data.table)

#############################################################################################################
########################### INPUT: ##########################################################################
#############################################################################################################
# The goal is first to create a Master Table with:
# IID, Urban_brth, Urban, ethncity, Asthma total,  Asthma <40, Asthma age grps, genotypic array and principal components
# rows: all samples from UKBiobank
# Then filter the table to a table with the same columns but
# rows: only from bgen sample 
# Then extract tables with rows: only European, only form bgen file and only either urban or rural 
# Asthma_v0 with NA entries are not explicitly excluded, as SAIGE can handle NA phenotype entries.
# However, the result is such that no sample has NA in Asthma_v0, possibly due to 
# including only samples with genetically imputed data, hence samples who reported on Asthma_v0

# Asthma file with age groups, generated from an appropriate script
# Asthma:  IID | Asthma_v0 | Asthma_age_6_v0 | Asthma_age_4_v0 | Asthma_age_5_v0 | Asthma_v1 | Asthma_age_6_v1 | Asthma_age_4_v1 | Asthma_age_5_v1 | ...
# ur_ethn: Urban_birht classification at time of birth | Urban classification at time of assessment | info of full sample of European origin
# Covars is a table of Covariates like Age, Sex, PCs and genotypic array

urban_rural_ethnicity_f <- "E:\\Genomics\\Resources\\Output_data_ph2\\urb_rur_ethn.txt"
Asthma_with_age_grps_f <- "E:\\Genomics\\Resources\\Output_data_ph2\\asthma_age_groups.txt"
covar_f<-"E:\\Genomics\\Resources\\Output_data_ph2\\demogr_geno.txt"

sample_f <- "E:\\Genomics\\Resources\\UKB_sample_file_newest\\ukb_bgen_sample_new.txt"

output_folder <- "E:\\Genomics\\Resources\\Asthma_ur\\Pheno_Cov\\"

#############################################################################################################
############################################  EXECTUE:        ###############################################

##### Import and pre-process urban-rural ethnicity file:
gc()
ur_ethn <- fread(urban_rural_ethnicity_f)

##### Process Asthma file: Create Asthma_under40
##### Remark: 0 in one column of asthma iplies 0 for all columns from visit (i.e. 0 is control)
##### Remark: NA in one column could mean missing asthma data, whenever all other oclumns from visit are NA 
##### Remark: or a case, then there is 1 in the Asthma_v* and 1 in one ASthma_age_*, NA else (in same visit)
Asthma <- fread(Asthma_with_age_grps_f)
names(Asthma)[1] <- 'IID'

Asthma[, Asthma_under40_v0:= ifelse(is.na(Asthma_Age_4_v0) & is.na(Asthma_Age_6_v0), NA, 
                                    ifelse( Asthma_Age_4_v0 == 1 | Asthma_Age_6_v0 == 1, 
                                            1, 0))] 
#View(Asthma[ , .(IID, Asthma_v0, Asthma_under40_v0, Asthma_Age_6_v0, Asthma_Age_4_v0, Asthma_Age_5_v0)])
cols <- c("IID", "Asthma_v0", "Asthma_under40_v0", colnames(Asthma)[c(3:(ncol(Asthma)-1))])
Asthma <- Asthma[,..cols]

gc()

##### Process Covariates file:
Covars <- fread(covar_f)
names(Covars)[1] <- 'IID'
Covars <- select(Covars, IID, Genot_array, starts_with('PC'))

gc()



##########################################################################################################
##### Create a master table: 
##### rows == bgen sample, 
##### cols == all Asthma | all visits | all age groups | full ethnicity | full urban rural | covariates

Asthma_ur_ethn_cov <- as.data.table(left_join(ur_ethn, Asthma, by='IID'))
 rm(Asthma)
 rm(ur_ethn)
 gc()
Asthma_ur_ethn_cov <- as.data.table(left_join(Asthma_ur_ethn_cov, Covars, by='IID'))
 rm(Covars)
 gc()
#####
##########################################################################################################

##### Process sample file: sort and remove the first 15 rows, which are dummy samples

sample <- fread(sample_f)
sample <- sample[order(sample)]
names(sample) <- 'IID'
sample <- sample[IID > 0]

##### Process the master table: intersect with the sample
##### rows = bgen sample

Asthma_ur_ethn_cov <- as.data.table(left_join(sample, Asthma_ur_ethn_cov, by='IID'))
gc()

##########################################################################################################
##### Extract from master table: 

##### rows = bgen sample, rural, European 
##### cols = asthma all | Asthma age grps | covariates
cols <- c("IID", "Asthma_v0", "Asthma_under40_v0", "Asthma_Age_6_v0", "Asthma_Age_4_v0", "Asthma_Age_5_v0", "Genot_array", names(Asthma_ur_ethn_cov)[grepl('PC', names(Asthma_ur_ethn_cov))])
Asthma_rr_eur_c <- Asthma_ur_ethn_cov[(European == 1) & (Urban_birth == 0) & (Urban == 0), ..cols]  
Asthma_rr_eur_c <- data.table(FID = Asthma_rr_eur_c[,IID], Asthma_rr_eur_c)
gc()

N_r <- nrow(Asthma_rr_eur_c[!is.na(Asthma_v0), ]) 

write.table(Asthma_rr_eur_c, paste(output_folder, "Asthma_rr_eur.txt",sep=""), append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)

View(fread(paste(output_folder, "Asthma_rr_eur.txt",sep="")))


##### rows = bgen sample, urban, European 
##### cols = asthma all | Asthma age grps | covariates
cols <- c("IID", "Asthma_v0", "Asthma_under40_v0", "Asthma_Age_6_v0", "Asthma_Age_4_v0", "Asthma_Age_5_v0", "Genot_array", names(Asthma_ur_ethn_cov)[grepl('PC', names(Asthma_ur_ethn_cov))])
Asthma_uu_eur_c <- Asthma_ur_ethn_cov[(European == 1) & (Urban_birth == 1) & (Urban == 1), ..cols]  
Asthma_uu_eur_c <- data.table(FID = Asthma_uu_eur_c[,IID], Asthma_uu_eur_c)
gc()

N_u <- nrow(Asthma_uu_eur_c[!is.na(Asthma_v0), ]) 

write.table(Asthma_uu_eur_c, paste(output_folder, "Asthma_uu_eur.txt", sep=""), append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)

View(fread(paste(output_folder, "Asthma_uu_eur.txt",sep="")))

##### rows = bgen sample, urban, European 
##### cols = asthma all | Asthma age grps | covariates


##########################################################################################################

