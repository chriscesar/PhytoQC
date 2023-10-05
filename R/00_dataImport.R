## 00_dataImport.R
# import data from >250 workbooks into long format

# load packages #
ld_pkgs <- c("tidyverse","purrr","readxl")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)


# Set the path to the directory containing your Excel files
path_to_files <- "path/to/files"

# Get a list of all Excel files in the directory
excel_files <- list.files(path_to_files, pattern = ".xls$|.xlsx$", full.names = TRUE)

# Function to read specific cells from a workbook and associate with filename
read_excel_data <- function(file_path) {
  # Get the filename with the extension
  file_name <- basename(file_path)
  
  # Read cells A1:B13 from the worksheet "Sample_Details"
  sample_details <- readxl::read_excel(file_path, sheet = "Sample_Details", range = "A1:B13", col_names = FALSE)
  
  # Assign column names to sample_details
  colnames(sample_details) <- c("Variable", "Value")
  
  # Add "Type" column to sample_details with value "SampleDetails"
  sample_details$Type <- "SampleDetails"
  
  # Read cells B6:I238 from the worksheet "Input_Data"
  input_data <- readxl::read_excel(file_path, sheet = "Input_Data", range = "B6:I238", col_names = FALSE)
  
  # Assign custom column names to input_data
  colnames(input_data) <- c("Taxon", "Qualifier", "densOriginal", "densBaseplate", "densReplicate", "propOriginal", "propBaseplate", "propReplicate")
  
  # Concatenate "Taxon" and "Qualifier" into a new variable "Tax_Qual"
  input_data$Tax_Qual <- ifelse(!is.na(input_data$Qualifier), paste(input_data$Taxon, input_data$Qualifier, sep = "_"), input_data$Taxon)
  
  # Add "Type" column
  input_data$Type <- "TaxonAbundance"
  
  # Remove "Taxon" and "Qualifier" variables
  input_data <- input_data %>%
    select(Tax_Qual, densOriginal, densBaseplate, densReplicate, propOriginal, propBaseplate, propReplicate, Type) %>%
    # Remove rows where all three variables are NA or zero
    filter(!(is.na(densOriginal) & is.na(densBaseplate) & is.na(densReplicate)) &
             !(densOriginal == 0 & densBaseplate == 0 & densReplicate == 0))
  
  # Read the "QA_Summary" table from cells A1:B5
  qa_summary <- readxl::read_excel(file_path, sheet = "QA_Summary", range = "A1:B5", col_names = TRUE)
  
  # Add "Type" column
  qa_summary$Type <- "QASummary"
  
  # Return a named list with filename, data frames, and QA_Summary table
  return(list(filename = file_name, sample_details = sample_details, input_data = input_data, qa_summary = qa_summary))
}

# Use purrr::map() to apply the function to all files and read the data
extracted_data <- map(excel_files, read_excel_data)
