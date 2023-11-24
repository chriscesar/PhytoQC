## 00_dataImport.R
# import data from >250 workbooks into long format

# admin ####
# load packages #
ld_pkgs <- c("tidyverse","purrr","readxl")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)


# Set the path to the parent directory containing Excel files
path_to_files <- "C:/Users/ccesar/OneDrive - Defra/Desktop/MyFiles-notBackedUp/Tasks/Phyto QC/QA Project"

# Get a list of all Excel files in the directory and sub-directories
excel_files <-
  list.files(
    path_to_files,
    pattern = ".xls$|.xlsx$", #files ending in Excel file types
    full.names = TRUE, #returns full filepath
    recursive = TRUE # look in subfolders
  )

# Function to read specific cells from a workbook and associate with filename
read_excel_data <- function(file_path) {
  
  ## Extract info from Sample_Details worksheet ##

  # Read cells A1:B13 from the worksheet "Sample_Details"
  sample_details <- readxl::read_excel(file_path,
                                       sheet = "Sample_Details",
                                       range = "A1:B13",
                                       col_names = FALSE)
  
  # Assign column names to sample_details
  colnames(sample_details) <- c("Variable", "Value")
  
  # Handle potentially blank cells and convert them to NA
  sample_details$Value <-
    ifelse(sample_details$Value == "", NA, sample_details$Value)
  
  # Read additional values from Input_Data sheet (cells C1 and C2)
  water_sample_vol <- readxl::read_excel(
    file_path,
    sheet = "Input_Data",
    range = "C1",
    col_names = FALSE,
    col_types = "text"
    )
  
  sub_sample_vol <- readxl::read_excel(
    file_path,
    sheet = "Input_Data",
    range = "C2",
    col_names = FALSE,
    col_types = "text"
    )
  
  # Add rows for WaterSampleVol_ml and SubSampleVol_ml to sample_details
  sample_details <- rbind(
    sample_details,
    c("WaterSampleVol_ml", as.character(water_sample_vol[[1]])),
    c("SubSampleVol_ml", as.character(sub_sample_vol[[1]]))
    )
  
  # Rename rows to maintain correct order
  sample_details[, 1] <- c(
    "SD01_AnaylsisLab",
    "SD02_LabSwap",
    "SD03_OriginalAnalyst",
    "SD04_AnalysisDateOrig",
    "SD05_NameOfSurvey_WFD",
    "SD06_SampleDate",
    "SD07_EAOldSiteCode",
    "SD08_EAWIMSCode",
    "SD09_InternalSampleID",
    "SD10_AuditAnalyst",
    "SD11_AuditDateBaseRec",
    "SD12_AuditDateRepSub",
    "SD13_Comments",
    "SD14_WaterSampleVol_ml",
    "SD15_SubSampleVol_ml"
    )
  
  # Convert sample_details to wide format
  sample_details <-
    spread(sample_details, key = Variable, value = Value)
  
  # Extract file name with extension from the file path
  file_name <- basename(file_path)
  
  # Add the file_name to the sample_details data frame
  sample_details <- cbind("SD00FileName" = file_name, sample_details)
  
  ## Extract info from Input_Data worksheet ##
  
  # Read cells B6:I238 from the worksheet "Input_Data" as character values
  input_data <- readxl::read_excel(
    file_path,
    sheet = "Input_Data",
    range = "B6:I238",
    col_names = FALSE,
    col_types = "text"
    )
  
  # Assign custom column names to input_data
  colnames(input_data) <- c(
    "Taxon",
    "Qualifier",
    "Original_dens",
    "Baseplate_dens",
    "Replicate_dens",
    "Original_prop",
    "Baseplate_prop",
    "Replicate_prop"
    )
  
  # Concatenate "Taxon" and "Qualifier" into a new variable "Tax_Qual"
  input_data$Tax_Qual <- ifelse(
    !is.na(input_data$Qualifier),
    paste(input_data$Taxon, input_data$Qualifier, sep = "_"),
    input_data$Taxon
    )
  
  # Remove "Taxon" and "Qualifier" variables
  input_data <- input_data %>%
    select(
      Tax_Qual,
      Original_dens,
      Baseplate_dens,
      Replicate_dens,
      Original_prop,
      Baseplate_prop,
      Replicate_prop
    ) %>%
    # Remove rows where all variables are NA or zero
    filter(
      !(
        is.na(Original_dens) &
          is.na(Baseplate_dens) & is.na(Replicate_dens) &
          is.na(Original_prop) &
          is.na(Baseplate_prop) & is.na(Replicate_prop)
      ) &
        !(
          Original_dens == 0 & Baseplate_dens == 0 & Replicate_dens == 0 &
            Original_prop == 0 &
            Baseplate_prop == 0 & Replicate_prop == 0
        )
    ) %>%
    
    ##convert to long
    pivot_longer(-Tax_Qual,
                 names_to = c("AnalysisType", ".value"),
                 names_sep = "_") %>%
    #remove NA or 0 values
    filter(!(is.na(dens)) & !(dens == 0))
  
  ## Extract info from QA_Summary worksheet ##
  
  # Read the "QA_Summary" table from cells A1:B5
  qa_summary <- readxl::read_excel(file_path,
                                   sheet = "QA_Summary",
                                   range = "A1:B5",
                                   col_names = TRUE)
  
  # Rename specific row values in qa_summary
  qa_summary[, 1] <-
    c("QA01_TotAbund",
      "QA02_DomTax",
      "QA03_SharedTax",
      "QA04_Final")
  
  # Convert qa_summary to wide format
  qa_summary <- spread(qa_summary, key = Test, value = Outcome)
  
  # Replicate sample_details and qa_summary rows to match the number of rows in input_data
  replicated_sample_details <- sample_details[rep(seq_len(nrow(sample_details)), nrow(input_data)), ]
  replicated_qa_summary <- qa_summary[rep(seq_len(nrow(qa_summary)), nrow(input_data)), ]
  
  # Bind replicated sample_details and qa_summary with input_data
  merged_data <- as_tibble(cbind(replicated_sample_details, input_data, replicated_qa_summary))
  
  
  # Return a named list with filename, data frames, and renamed QA_Summary table
  return(
    merged_data = merged_data
  )
  }

# Use purrr::map() to apply the modified function to all files and read the data
start.time <- Sys.time() # start timer

# extract data as a list
extracted_data_list <- purrr::map(excel_files, read_excel_data)

# Combine list elements into a single data frame using dplyr::bind_rows()
extracted_data <- dplyr::bind_rows(extracted_data_list)

### save data
## all data
write.csv(extracted_data,
          file = "data/out/extracted_data_ALL.csv",
          row.names = FALSE)

## Lab Swap
extracted_dataLS <- extracted_data %>% 
  dplyr::filter(.,SD02_LabSwap == "Yes")

write.csv(extracted_dataLS,
          file = "data/out/extracted_data_LabSwap.csv",
          row.names = FALSE)
rm(extracted_dataLS)

## Non-Lab Swap
extracted_dataNLS <- extracted_data %>% 
  dplyr::filter(.,SD02_LabSwap == "No") %>% 
  dplyr::filter(., !is.na(SD00FileName))

write.csv(extracted_dataNLS,
          file = "data/out/extracted_data_NonLabSwap.csv",
          row.names = FALSE)
rm(extracted_dataNLS)

end.time <- Sys.time() #stop timer
time.taken <- round(end.time - start.time,2)
time.taken

### tidy up ####
rm(extracted_data, extracted_data_list,
   end.time,
   excel_files,
   path_to_files,
   start.time,
   time.taken,
   read_excel_data)

detach("package:readxl", unload = TRUE)
detach("package:tidyverse", unload = TRUE)
detach("package:purrr", unload = TRUE)

### TO DO:
# investigate worksheet errors with lab
# where do NA values arise? - suspect linked to 'partially duplicated' "WB QA 150 TCO007P Nov 2022.xlsx"
# check that all files are imported into the data