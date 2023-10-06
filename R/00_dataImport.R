## 00_dataImport.R
# import data from >250 workbooks into long format

# load packages #
ld_pkgs <- c("tidyverse","purrr","readxl")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)


# Set the path to the directory containing your Excel files
path_to_files <- "path/to/data"

# Get a list of all Excel files in the directory
excel_files <-
  list.files(path_to_files, pattern = ".xls$|.xlsx$", full.names = TRUE)

# Function to read specific cells from a workbook and associate with filename
read_excel_data <- function(file_path) {

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
    "01_AnaylsisLab",
    "02_LabSwap",
    "03_OriginalAnalyst",
    "04_AnalysisDate_Orig",
    "05_NameOfSurvey_WFD",
    "06_SampleDate",
    "07_EAOldSiteCode",
    "08_EAWIMSCode",
    "09_InternalSampleID",
    "10_AuditAnalyst*",
    "11_AuditDate_BaseRec",
    "12_AuditDate_RepSub",
    "13_Comments",
    "14_WaterSampleVol_ml",
    "15_SubSampleVol_ml"
  )
  
  # Convert sample_details to wide format
  sample_details <-
    spread(sample_details, key = Variable, value = Value)
  
  # Extract file name with extension from the file path
  file_name <- basename(file_path)
  
  # Add the file_name to the sample_details data frame
  sample_details <- cbind("00FileName" = file_name, sample_details)
  
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
  
  # Read the "QA_Summary" table from cells A1:B5
  qa_summary <- readxl::read_excel(file_path,
                                   sheet = "QA_Summary",
                                   range = "A1:B5",
                                   col_names = TRUE)
  
  # Rename specific row values in qa_summary
  qa_summary[, 1] <-
    c("QA1_TotAbund", "QA2_DomTax", "QA3_SharedTax", "QA4_Final")
  
  # Convert qa_summary to wide format
  qa_summary <- spread(qa_summary, key = Test, value = Outcome)
  
  # Return a named list with filename, data frames, and renamed QA_Summary table
  return(
    list(
      sample_details = sample_details,
      input_data = input_data,
      qa_summary = qa_summary
    )
  )
}

# Use purrr::map() to apply the modified function to all files and read the data
start.time <- Sys.time()
extracted_data <- purrr::map(excel_files, read_excel_data)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

#### to do:
# convert sample_details and qa_summary dfs to WIDE format.
# append file name as a variable to the sample_details object
# join sample_details and qa_summary dfs to the input_data df
