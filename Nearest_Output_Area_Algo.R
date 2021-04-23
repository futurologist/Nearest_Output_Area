library(dplyr)
library(data.table)

#########################################################################################################

# a function that relables the headers of the table initially extracted from the UKBiobank
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
################################## functions for UKB data extraction: ###################################


# extracts a relevant subtable of the raw UKB table, used in generate_location_table 
extract_UKB_subtable <- function(data_base, list_of_codes){
  headers <- names(data_base)
  pos <- c(1)
  for(code in list_of_codes){
    pos <- c(pos, grep(code, headers))
  } 
  subt <- data_base[, ..pos]
  subt <- subt[order(subt[, 1]), ]
  subt
}

# bulids and labels the complete table for home location at birth, 
# at visits 0,1,2 and date and home location history
generate_location_table <- function(filepath_in){
  codes <- c(130, 129, 20074, 20075, 20118, 
             22700, 22702, 22704, 20115)
  labels <- c('Brth_plc_east', 'Brth_plc_north', 'Home_loc_assess_east', 'Home_loc_assess_north', 'Urban_Rural', 
              'Date_1st_at_loc', 'home_loc_east', 'home_loc_north', 'Country_brth')
  fields <- paste('f.', codes, sep="")
  data_base <- fread(filepath_in)
  data_base <- extract_UKB_subtable(data_base, fields)
  names(data_base)[1] <- 'IID'
  array_length <- rep (1, length(fields))
  array_length[which(grepl('f.227', fields))] <- 15
  instances   <-   rep (1, length(fields))
  instances[which(grepl('f.129', fields))] <- 3
  instances[which(grepl('f.130', fields))] <- 3
  instances[which(grepl('f.20074', fields))] <- 3
  instances[which(grepl('f.20075', fields))] <- 3
  names(data_base) <- relabel(data_base, fields, array_length, instances, labels)
  return(data_base)
}

# extracts north-east coordinates of home location at birth, visit 0, 
extract_birth_location <- function(ukb_table){
  table_out <- select(ukb_table, IID, Brth_plc_east_v0, Brth_plc_north_v0)
  table_out <- filter(table_out, !is.na(ukb_table$Brth_plc_east_v0) &  ukb_table$Brth_plc_east_v0 !=-1)
  return(as.data.table(table_out))
}

# extracts north-east coordinates of home location at assessment, visit 0, 
# exctracts the category of urbanization 
# classifies urban rural locations base on category of urbanization
extract_assess_location <- function(ukb_table){
  ukb_curr_loc <- select(ukb_table, IID, Home_loc_assess_east_v0, Home_loc_assess_north_v0, Urban_Rural)
  ukb_curr_loc <- filter(ukb_curr_loc, !is.na(ukb_curr_loc$Home_loc_assess_east_v0))
  ukb_curr_loc <- mutate(ukb_curr_loc, 
                         Urban = ifelse(is.na(ukb_curr_loc$Urban_Rural) | ukb_curr_loc$Urban_Rural == 9, NA, 
                                        ifelse(ukb_curr_loc$Urban_Rural %in% c(5, 11, 12), 1, 0)))
  
  return(as.data.table(ukb_curr_loc))
}

###################################### Urban/Rural 1991 England and Wales: ##################################

# excracts and puts together a table with
# north-east coordinates of inhabited locations and their urban-rural classification in terms of categories 
# the inhabited locations are 1991 wards for England and Wales, and 2001 ouput areas of Scotland  
make_urban_rural_GB <- function(Eng_1991_ur, Wls_1991_ur, Scot_2001_oa, Scot_2003_ur){
  en_1991 <- fread(Eng_1991_ur)
  wl_1991 <- fread(Wls_1991_ur)
  en_1991 <- select(en_1991,x,y,cat)
  wl_1991 <- select(wl_1991,x,y,cat)
  ew_1991 <- bind_rows(en_1991, wl_1991)
  rm(en_1991)
  rm(wl_1991)
  sc_ur <- fread(Scot_2001_oa)
  sc_ur <- select(sc_ur, c(2,3,4))
  sc_ur_03 <- fread(Scot_2003_ur)
  names(sc_ur_03)[1] <- 'name'
  sc_ur <- left_join(sc_ur, sc_ur_03, by='name')
  sc_ur <- select(sc_ur, x, y, UR8_2003_2004)
  names(sc_ur)[3] <- c('cat') 
  return(as.data.table(bind_rows(ew_1991, sc_ur)))
}

######################## functions for Nearest Neighbour matching for England, Wales and Scotland: ########

# for a home-location with given north-east coordinates, 
# finds the index of the nearest ward from England_Whalse_Scotland table of 
# ward coordinates and their urban rural cathegory 
index_of_nrst_nbr <- function(point_x, point_y, set_in){
  d_2 <- function(p_x, p_y){ return((p_x-point_x)^2 + (p_y-point_y)^2) }
  set <- mutate(set_in, distance_2 = d_2(set_in$x, set_in$y))
  i = which.min(set$distance_2)
  return(i)
}

# by using the function index_of_nrst_nbr(),
# adds a column of indexes where for each row, each index entry is the index of the 
# ward from the England_Whalse_Scotland table of 
# ward coordinates, nearest to the home-location listed in that row 
add_col_of_idxs <- function(ukb_in, wards){
  ukb <- mutate(ukb_in, index = 0)
  nnb <- function(r){ return(index_of_nrst_nbr(as.double(r[2]), as.double(r[3]), wards)) }
  ukb[,'index'] <- as.data.table(apply(ukb, 1, nnb))
  return(as.data.table(ukb))
}

# takes the output from the funciton add_col_of_idxs(ukb_in, wards)
# and adds and calculates the entries of columns
# containing the coordiantes of the nearest ward, 
# the urban-rural category of the nearest ward,
# and the distance to the nearest ward,
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
  return(as.data.table(ukb))
}

# executes one after the other the functions 
# add_col_of_idxs() and then
# add_cols_of_coord_and_urb()
# to rpoduce a final table with necessary data, ready for urban rural 
# binary classification
assign_nearest_ward <- function(ukb_brt_loc, ews){
  start_time <- Sys.time()
  ukb_ur <- add_col_of_idxs(ukb_brt_loc, ews)
  end_time <- Sys.time()
  print(end_time - start_time)
  ukb_ur <- add_cols_of_coord_and_urb(ukb_ur, ews)
  return(ukb_ur)
}

#########################################################################################################

################################### YOUR INPUT GOES HERE: ###############################################

filepath_in <- 'E:\\Genomics\\Resources\\Output_data_ph1\\ukb41654.r.tmp.tab'  # "C:\\MY_FOLDERS\\Asthma_and_Pain\\UKB_raw_data\\ukb41654.r.tmp.tab" #string
filepath_out <- 'E:\\Genomics\\Resources\\Output_data_ph2\\UKB_urban_1991.txt' # 'E:\\Genomics\\Resources\\Output_data_ph2\\ukb_loc_data.txt' #"C:\\MY_FOLDERS\\Asthma_and_Pain\\UKB_raw_data\\ukb_loc_data.txt"

# bulids and labels the complete table for home location at birth, 
# at visits 0,1,2 and date and home location history
ukb_data <- generate_location_table(filepath_in)

# extracts north-east coordinates of home location at birth, visit 0, 
ukb_birth_loc <- extract_birth_location(ukb_data)

# extracts north-east coordinates of home location at assessment, visit 0, 
# exctracts the category of urbanization 
# classifies urban rural locations base on category of urbanization
ukb_assess_loc <- extract_assess_location(ukb_data)


filepath_Eng_1991_ur <- "C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\England_urbrur_classn_1991\\england_urbrur_classn_1991.csv"
filepath_Wls_1991_ur <- "C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\Wales_urbrur_classn_1991\\wales_urbrur_classn_1991.csv"
filepath_Scot_2001_oa <- 'C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\Scotland_2001\\Scotland_oa_2001_gen3\\scotland_oa_2001_gen3.csv'
filepath_Scot_2003_ur <- 'C:\\MY_FOLDERS\\Asthma_and_Pain\\urban_rural\\Scotland_2001\\oa2001_urban_rural_2003_2004.csv'

# excracts and puts together a table with
# north-east coordinates of inhabited locations and their urban-rural classification in terms of categories 
# the inhabited locations are 1991 wards for England and Wales, and 2001 ouput areas of Scotland  
Eng_Wls_Scot <- make_urban_rural_GB(filepath_Eng_1991_ur, 
                                filepath_Wls_1991_ur, 
                                filepath_Scot_2001_oa, 
                                filepath_Scot_2003_ur)

# takes the table with coordinates of birth home location from UKB 
# and adds and calculates the entries of the new columns,
# containing the coordiantes of the nearest ward, 
# the urban-rural category of the nearest ward,
# and the distance to the nearest ward
ukb_ur_at_brth <- assign_nearest_ward(ukb_birth_loc, Eng_Wls_Scot)

# decide about the definition of urban vs rural. If places of category 1 are considered urban whilw everything else is rural, choose TRUE; 
# if places of categories 1 and 2 are considered urban, choose FALSE
only_cat_1_is_urban <- TRUE

if(only_cat_1_is_urban){  
  ukb_ur_at_brth_1<- mutate(ukb_ur_at_brth, Urban_birth = ifelse(ukb_ur_at_brth$category == 1, 1, 0))
} else {
  ukb_ur_at_brth_1<- mutate(ukb_ur_at_brth, Urban_birth = ifelse(ukb_ur_at_brth$category == 1 | ukb_ur_at_brth$category == 2, 1, 0))
}
  
write.table(ukb_ur_at_brth, filepath_out, 
            append = FALSE, sep = "\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
