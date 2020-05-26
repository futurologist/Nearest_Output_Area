library(dplyr)
library(data.table)

#########################################################################################################
####

relabel <- function(table, fields, arrays_length, instances, labels){
  headers <- names(table)
  #n_headers <- length(headers)
  n_fields <- length(fields)
  h <- 2
  for(i in 1:n_fields){
    for(v in 1:instances[i]){
      if(arrays_length[i] != 1){
        for(l in 1:arrays_length[i]){
          if(instances[i] == 1){
            headers[h] <- paste(labels[i], "_", as.character(l), sep = "")
          }
          else{
            headers[h] <- paste(labels[i], "_", as.character(l),"_v", as.character(v-1), sep = "")
          }
          h <- h+1
        }
      } 
      else{
        if(instances[i] == 1){
          headers[h] <- labels[i]
        }
        else{
          headers[h] <- paste(labels[i], "_v", as.character(v-1), sep = "")
        }
        h <- h+1
      }
    }
  }
  return(headers)
}

#########################################################################################################

extract_UKB_subtable <- function(data_base, list_of_codes){
  headers <- names(data_base)
  pos <- c(1)
  for(code in list_of_codes){
    pos <- c(pos, grep(code, headers))
  } 
  subt <- data_base[, ..pos]
  subt <- subt[order(subt[,'f.eid']), ]
  return(subt)
}


# For a set of fields that need to be extracted from the central UKB data table
# list a vector of strings, each string in the format "f.[Field_number].", for each field.
# Skip the first column, id column "f.eid", it is authomatically included

################################### YOUR INPUT GOES HERE: ###############################################

codes <- c(130, 129, 20074, 20075, 20118, 
           22700, 22702, 22704, 20115)
labels <- c('Brth_plc_east', 'Brth_plc_north', 'Home_loc_asses_east', 'Home_loc_asses_north', 'Urban_Rural', 
            'Date_1st_at_loc', 'home_loc_east', 'home_loc_north', 'Country_brth')


list_of_codes <- paste('f.', codes, sep="")

filepath_in <- "C:\\MY_FOLDERS\\Asthma_and_Pain\\UKB_raw_data\\ukb41654.r.tmp.tab" #string
filepath_out <- "C:\\MY_FOLDERS\\Asthma_and_Pain\\UKB_raw_data\\ukb_loc_data.txt"

#########################################################################################################
##################################### Execute Relabeling of the headers #################################

data_base <- fread(filepath_in)
ukb_data <- extract_UKB_subtable(data_base, list_of_codes)
#write.table(subtable, filepath_out, append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
rm(data_base)

####
#########################################################################################################

names(ukb_data)[1] <- 'IID'

fields <- list_of_codes
rm(list_of_codes)
array_length <- rep (1, length(fields))
array_length[which(grepl('f.227', fields))] <- 15
instances   <-   rep (1, length(fields))
instances[which(grepl('f.129', fields))] <- 3
instances[which(grepl('f.130', fields))] <- 3
instances[which(grepl('f.20074', fields))] <- 3
instances[which(grepl('f.20075', fields))] <- 3

names(ukb_data) <- relabel(ukb_data, fields, array_length, instances, labels)

write.table(ukb_data, filepath_out, append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)

###################################### Urban/Rural for current location ###########################################
ukb_data_1 <- fread(filepath_out)
ukb_brt_loc <- as.data.table(filter(ukb_data, !is.na(ukb_data$Brth_plc_east_v0) &  ukb_data$Brth_plc_east_v0 !=-1))
ukb_brt_loc <- select(ukb_brt_loc, IID, Brth_plc_east_v0, Brth_plc_north_v0)

ukb_curr_loc <- select(ukb_data, 1, 8, 11, 14)
ukb_curr_loc <- mutate(ukb_curr_loc, Urban = ifelse(ukb_curr_loc$Urban_Rural %in% c(5, 11, 12), 1, 
                                                    ifelse(ukb_curr_loc$Urban_Rural == 9, -1, 0)))

###################################### Urban/Rural 1991 England and Wales: ##################################################

en_ur_1991 <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\England_urbrur_classn_1991\\England_urbrur_classn_1991.csv")
wl_ur_1991 <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\Wales_urbrur_classn_1991\\wales_urbrur_classn_1991.csv")

### England 1991:

en_1991 <- select(en_ur_1991,x,y,cat)

### Wales 1991:

wl_1991 <- select(wl_ur_1991,x,y,cat)

### merge England and Wales 1991 urban information 

# England/Wales: Urban --> Rural <=> 1 --> 6

ew_1991 <- bind_rows(en_1991, wl_1991)

###########################  Urban/Rural 2001 Scotland:    ###########################

sc_oa_cntr_01 <- fread('C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\Scotland_2001\\Scotland_oa_2001_gen3\\scotland_oa_2001_gen3.csv')
sc_oa_cntr_01 <- select(sc_oa_cntr_01, c(2,3,4))

sc_ur_oa <- fread('C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\Scotland_2001\\oa2001_urban_rural_2003_2004.csv')
names(sc_ur_oa)[1] <- 'name'

# Scotalnd: Urban --> Rural <=> 1 --> 6 or 8
sc_ur <- left_join(sc_oa_cntr_01, sc_ur_oa, by='name')
sc_ur <- select(sc_ur, x, y, UR8_2003_2004)
sc_ur <- as.data.table(sc_ur)
names(sc_ur)[3] <- c('cat') 

###########################  Urban/Rural for England, Wales and Scotland: ###########################

ews <- bind_rows(ew_1991, sc_ur)

###########################  Nearest Neighbour Urban/Rural for England, Wales and Scotland: ##################

index_of_nrst_nbr <- function(point_x, point_y, set_in){
  d_2 <- function(p_x, p_y){ return((p_x-point_x)^2 + (p_y-point_y)^2) }
  set <- mutate(set_in, distance_2 = d_2(set_in$x, set_in$y))
  i = which.min(set$distance_2)
  return(i)
}

add_col_of_idxs <- function(ukb_in, wards){
  ukb <- mutate(ukb_in, index = 0)
  nnb <- function(r){ return(index_of_nrst_nbr(as.double(r[2]), as.double(r[3]), wards)) }
  ukb[,'index'] <- as.data.table(apply(ukb, 1, nnb))
  return(ukb)
}

add_cols_of_coord_and_urb <- function(ukb_in, wards){
  M <- as.matrix(wards)
  nnb <- function(i,j,M){ return(M[i,j]) }
  ukb <- mutate(ukb_in, ward_east = nnb(index,1,M), 
                ward_north = nnb(index,2,M), 
                category=nnb(index,3,M),
                distance = sqrt(((Brth_plc_east_v0) - ward_east)^2 +
                                ((Brth_plc_north_v0) - ward_north)^2)
                #distance = sqrt(((Home_loc_asses_east_v0) - as.double(ward_east))^2 +
                #                ((Home_loc_asses_north_v0) - as.double(ward_north))^2)
  )
  return(ukb)
}

start_time <- Sys.time()
ukb_ur <- add_col_of_idxs(ukb_brt_loc, ews)
end_time <- Sys.time()
end_time - start_time

ukb_urb <- add_cols_of_coord_and_urb(ukb_ur, ews)

write.table(ukb_urb, "C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\output_data\\UKB_urban_1991.txt", append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)

###########################  Nearest Neighbour Urban/Rural for England, Wales and Scotland in 2011: ##################

# REPLICATION:

###################################### Urban/Rural 2011 England and Wales: ##################################################

en_ur_2011 <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\England_oa_ru_classn_2011\\england_oa_ru_classn_2011.csv")
wl_ur_2011 <- fread("C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\Wales_oa_ru_classn_2011\\wales_oa_ru_classn_2011.csv")

en_2011 <- select(en_ur_2011,x,y,ruc11)

### Wales 1991:

wl_2011 <- select(wl_ur_2011,x,y,ruc11)

cat <- distinct(en_2011, ruc11)
names(en_2011)[3] <- 'cat'
names(wl_2011)[3] <- 'cat'


### merge England and Wales 1991 urban information 

# England/Wales: Urban --> Rural <=> 1 --> 6

ew_2011 <- bind_rows(en_2011, wl_2011)

sc_ur[,'cat'] <- as.character(sc_ur$cat)

ews_11 <- bind_rows(ew_2011, sc_ur)

start_time <- Sys.time()
ukb_ur_11 <- add_col_of_idxs(filter(ukb_curr_loc, !is.na(ukb_curr_loc$Home_loc_asses_east_v0)), ews_11)
end_time <- Sys.time()
end_time - start_time

ukb_urb_11 <- add_cols_of_coord_and_urb(ukb_ur_11, ews_11)

write.table(ukb_urb_11, "C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\output_data\\UKB_urban_2011.txt", append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)

