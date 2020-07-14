library(dplyr)
library(data.table)

############################################ FUNCTIONS #######################################

# Asthma_in:  IID | Asthma_v0 | Asthma_age_6_v0 | Asthma_age_4_v0 | Asthma_age_5_v0 | Asthma_v1 | Asthma_age_6_v1 | Asthma_age_4_v1 | Asthma_age_5_v1 | ...
create_asthma_under40 <- function(Asthma_in){
  Asthma_1 <- filter(Asthma_in, !is.na(Asthma_in$Asthma_v0))
  Asthma_1 <- mutate(Asthma_1, 
                     Asthma_under40_v0 = ifelse(is.na(Asthma_1$Asthma_Age_4_v0) & is.na(Asthma_1$Asthma_Age_6_v0), NA, 
                                                ifelse( Asthma_1$Asthma_Age_4_v0 == 1 | Asthma_1$Asthma_Age_6_v0 == 1, 
                                                        1, 0)))
  return(as.data.table(select(Asthma_1, IID, Asthma_under40_v0)))
}

join_phenotypes <- function(Asthma_in, Urban_birth_in, Urban_in, IID_w_in){
  Asthma_1 <- create_asthma_under40(Asthma_in)
  Urban_1 <- select(Urban_in, 1, Urban)
  Urban_1 <- inner_join(select(Urban_birth_in, 1, Urban_birth), Urban_1, by='IID')
  Asthma_1 <- inner_join(Asthma_1, Urban_1, by='IID')
  return(as.data.table(inner_join(IID_w_in, Asthma_1, by='IID')))
}

setup_asthma_u_or_r <- function(Asthma_in, Urban_birth_in, Urban_in, IID_w_in, u_or_r){
  asthma_uu <- join_phenotypes(Asthma_in, Urban_birth_in, Urban_in, IID_w_in)
  asthma_uu <- filter(asthma_uu, asthma_uu$Urban_birth == u_or_r & asthma_uu$Urban == u_or_r)
  return(as.data.table(select(asthma_uu, 1, 2)))
}

#############################################################################################################

# Asthma with 3 age groups, 0-18, 18-40, 40-100
# Urban_birth, Urban, 
# IID_Eur (Europeans only from the bgen sample!)
# Exclude all samples with asthma diagnosed in 40-100 and all NA, NA, NA
# include all non NA samples with asthma in 

# Asthma file with age groups, generated from a nappropriate script
# Asthma_in:  IID | Asthma_v0 | Asthma_age_6_v0 | Asthma_age_4_v0 | Asthma_age_5_v0 | Asthma_v1 | Asthma_age_6_v1 | Asthma_age_4_v1 | Asthma_age_5_v1 | ...
# Urban classification at time of assessment
# Urban_birht classification at time of birth
#IID_w of all people from the bgen sample that are of European origin
#Covars is a table of Covariates like Age, Sex, PCs and genotypic array

Asthma <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\Output_data_ph2\\asthma_age_groups.txt")
Urban <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\output_data\\UKB_urban_2011.txt")
Urban_birth <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\output_data\\UKB_urban_1991.txt")
IID_w <- fread("C:\\MY_FOLDERS\\Pain_GWAS\\PCA_results\\Ethnicity_outputs\\IDs_UKB_Whites_3_PCs.txt")
Covars <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\Output_data_ph2\\demogr_geno.txt")
Covars <- select(Covars, 1, Genot_array, starts_with('PC'))

names(Asthma)[1] <- 'IID'
names(IID_w)[1] <- 'IID'
names(Covars)[1] <- 'IID'
Covars <- left_join(IID_w, Covars, by='IID')

Asthma_uu <- setup_asthma_u_or_r(Asthma, Urban_birth, Urban, IID_w, 1)
Asthma_rr <- setup_asthma_u_or_r(Asthma, Urban_birth, Urban, IID_w, 0)
Asthma_uu <- rename(Asthma_uu,  Asthma_under40_urban_v0 = Asthma_under40_v0)
Asthma_rr <- rename(Asthma_rr,  Asthma_under40_rural_v0 = Asthma_under40_v0)
Pheno_Covar_table <- right_join(Asthma_rr, Covars, by='IID')
Pheno_Covar_table <- right_join(Asthma_uu, Pheno_Covar_table, by='IID')
Pheno_Covar_table <- cbind(select(Pheno_Covar_table, 1), Pheno_Covar_table)
names(Pheno_Covar_table)[1] <- 'FID'

write.table(Pheno_Covar_table, "C:\\MY_FOLDERS\\Asthma_and_Pain\\Pheno_and_Covariates\\Pheno_Cov_asthma_urb.txt", append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
View(fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\Pheno_and_Covariates\\Pheno_Cov_asthma_urb.txt"))
