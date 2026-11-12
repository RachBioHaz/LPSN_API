
#' Utilities and functions used in the LPSN Taxonomic correction walkthrough
#'
#' These functions were created by R Gore in April 2026 with some AI help. 
#' 
#' Use these as needed, however these functions are not being actively maintained,  
#' so I highly recommend you use a subset of your taxonomy file in a test run
#' to validate the function is returning the correct  and manually
#' checking the results in LPSN. A list of taxa that failed to return a result, 
#' is also provided to enable manual searches.
#' 
#' How long does it take? 
#' The whole process is shorter if you prepare your data, there will be fewer manual searches 
#' needed if the data is prepared properly. 
#' 
#' Are there records that cannot be retrieved?
#' The LPSN API currently does not return Candidatus records, these taxa names can be manually searched 
#' using the website and then selecting the parent for each record, working backward to create a lineage.
#' 
#' The LPSN API will not retrieve superseded synonyms taxa like Firmicutes. If starting from lower levels,
#' like genus or family, the code works backwards to create a lineage and will not encounter an error. 
#' If the search starts at the superseded synonym, it will reject as not found and a manual search is 
#' required. 
#' 
#' 1. Clean your data.
#' 
#' 1.1 Remove underscores and any orphaned characters/numbers that result.
#' AGRF often append letter to taxonomy levels to indicate which 
#' database the information came from, e.g. Enterococcus A avium, Enterococcus AB avium, 
#' Enterococcus AD avium B. Sometimes databases have spurious names that will need to be fixed 
#' before the LPSN API call or after the run has finished, e.g. Candidatus Koribacter Koribacter. 
#' 
#' 1.2 Any solo words should be removed from the species column, e.g. soil, bacterium, rectale, etc.
#' In AGRF files, both the genus and species names populate the species field. Check what is the norm 
#' for the specific database you're using.
#' 
#' 1.3 Some databases are not able to resolve which closely related genera a sequence 
#' belongs to and will provide a concatenated name, e.g. genus
#' Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium,  this will be ignored if left as is,
#' and the search will begin at the next highest taxonomic level. You can choose to manually 
#' change the name before running it or deal with it afterwards  
#'  
#' 2. Place holder names from the various databases will be ignored, e.g. Genus CAG0056 
#' isn't a real name, so it will be ignored and the next highest valid taxonomic level is used.
#' You may prefer to delete these before the run to reduce the number of API calls that are made
#' 
#' 3.To shorten the run, use a list of unique names rather that the same names multiple times in the 
#' search list, this will result in fewer API calls and a quicker run
#' 
#' After cleaning the data, the API run took just over 4 hours, followed by the manual validation 
#' of a subset of results search of rejected taxa


# Remove Taxonomy Level Prefixes ##############################################

# This function cleans a taxonomy data frame by stripping the standard prefixes
# (e.g., "d__", "p__") used to indicate lineage level.
#
# @param tax A data frame or tibble containing taxonomy rank columns.
#
# @return A modified data frame (tibble) with prefixes removed from available
#   taxonomic columns.
# @export
#
# @examples
# \dontrun{
# df_clean <- remove_taxa_prefix(df)
# }
#

remove_taxa_prefix <- function(tax) {
  library(dplyr)
  library(stringr)
  
  
  # Define expected taxonomy columns and their prefixes
  taxa_prefixes <- c(
    Domain  = "^d__",
    Kingdom = "^k__",
    Phylum  = "^p__",
    Class   = "^c__",
    Order   = "^o__",
    Family  = "^f__",
    Genus   = "^g__",
    Species = "^s__"
  )
  
  
  
  # Convert to data frame
  tax <- base::as.data.frame(tax)
  
  # Identify taxonomy columns that are present
  existing_taxa <- base::intersect(
    names(taxa_prefixes),
    base::names(tax)
  )
  
  # Apply prefix removal only to existing columns
  if (base::length(existing_taxa) > 0) {
    
    tax <- tax |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(existing_taxa),
          ~ stringr::str_replace(
            .x,
            taxa_prefixes[dplyr::cur_column()],
            "")
          )
        )
    } else {
    
    base::warning(
      "No expected taxonomy columns were found. Returning original data.")
      }
  
  # Print final column names
  base::message("Final column names:")
  base::print(base::names(tax))
  
  return(tax)
}



# # Usage example
# # Remove taxonomy prefixes
#
# tax_clean <- remove_taxa_prefix(tax_df)
#
# ### will also display final column names
#
# # View cleaned taxonomy table
# tax_clean



#  Removing database flags appended to taxa lineage ############################ 

#' Clean Orphaned Letters and Report Changes
#' 
#' This function, clean_orphaned_letters(), is designed to remove isolated, 
#' short uppercase letter codes (e.g., "A", "BC", "XYZ") from a 
#' specified column of a data frame, while preserving a record of 
#' exactly which rows were changed and how they were modified.
#'
#' @param df The dataframe to process.
#' @param target_col The name of the column to clean (e.g., "Species").
#' @param max_letters The maximum number of uppercase letters to remove (default 3).
#'
#' @return A list containing the cleaned dataframe and a summary of changes.
clean_orphaned_letters <- function(df, target_col, max_letters = 3) {
  
  # Create the regex pattern for 1 to max_letters uppercase letters
  # Regex components:
  # (^| ) specifies the start of the string OR a space before the letters.
  # [A-Z]{1,3} Matches: one, two or 3 UPPERCASE letters
  # ( |$) Requires the letters to be followed by either a space or end of string
  
  pattern <- base::sprintf("(^| )[A-Z]{1,%d}( |$)", max_letters)

  
  # Identify rows that WILL be changed 
  affected_indices <- base::which(base::grepl(pattern, df[[target_col]]))
  
  # Store original values for the report
  original_values <- df[[target_col]][affected_indices]
  
  # Perform the cleaning on the specific column
  # Step 1: Replace patterns with space
  df[[target_col]] <- base::gsub(pattern, " ", df[[target_col]])
  
  # Step 2: Collapse double spaces and trim
  df[[target_col]] <- base::trimws(base::gsub("\\s+", " ", df[[target_col]]))
  
  # Store the new values for the report
  cleaned_values <- df[[target_col]][affected_indices]
  
  # Create a comparison table for the user to check
  report <- base::data.frame(
    Row_Index = affected_indices,
    Original  = original_values,
    Cleaned   = cleaned_values,
    stringsAsFactors = FALSE
  )
  
  base::message(base::sprintf("Cleaning complete. %d rows modified.", base::nrow(report)))
  
  # Return both the dataframe and the report
  base::return(base::list(
    data = df,
    report = report
  ))
}


## usage example
# # Run the function
# results <- clean_orphaned_letters(df, "Species", max_letters = 3)
# 
# # Access the cleaned data
# df_trimmed <- results$data
# 
# # Check the changes (shows only the rows that were modified)
# View(results$report)
# 
# # Quick look at the top of the report in the console
# head(results$report)

# Detection example
  # | Text              | Detected? |
  # | ----------------- | --------- |
  # | `"Bacteria A"`    | Yes       |
  # | `"A Bacteria"`    | Yes       |
  # | `"Bacteria AB"`   | Yes       |
  # | `"Bacteria ABC"`  | Yes       |
  # | `"Bacteria ABCD"` | No        |
  # | `"Bacteria abc"`  | No        |
  # 


# Identify Suspicious Taxonomy Names ######################################### 

#' 
#' This function, find_suspicious_taxa(), identifies potentially problematic
#' taxonomy names within a specified dataframe column. Taxa are flagged based
#' on user-defined quality-control criteria:
#'
#' 1. Presence of numerical characters
#' 2. Fewer than a specified number of alphabetic letters
#' 3. Multiple words within a taxonomy name (optional)
#' Please NOTE: taxa can be prefixed with Candidatus, 
#' so the list returned should be checked before taxa names are changed
#'
#' The function is intended as a taxonomy quality-control step after initial
#' cleaning. It does not modify the dataframe. Instead, it returns a report
#' containing suspicious taxonomy names and the criteria responsible for
#' flagging them.
#'
#' @param df The dataframe containing taxonomy information.
#'
#' @param target_col The taxonomy column to inspect
#' (e.g., "Phylum", "Genus", "Species").
#'
#' @param min_letters Minimum number of alphabetic characters required for a
#' taxon name to avoid being flagged (default = 5).
#'
#' @param check_multi_word Logical. Should taxa containing multiple words be
#' flagged? Default = FALSE.
#'
#' @param export_csv Logical. Should the report be exported as a CSV file?
#' Default = FALSE.
#'
#' @param file Filename for CSV output if export_csv = TRUE.
#'
#' @return A list containing:
#'
#' \describe{
#' \item{taxa}{A character vector of suspicious taxonomy names.}
#'
#' \item{report}{A dataframe containing:
#' \itemize{
#' \item Taxon
#' \item Has_Numbers
#' \item Letter_Count
#' \item Word_Count
#' \item Has_Multiple_Words
#' }}
#' }
#'
#' @examples
#'
#' genus_check <- find_suspicious_taxa(
#'   df = df_trimmed,
#'   target_col = "Genus",
#'   min_letters = 6,
#'   check_multi_word = TRUE,
#'   export_csv = TRUE,
#'   file = "genus_suspicious_taxa.csv"
#' )
#'
#' genus_check$report
#' 
#' 
#' #' phylum_check <- find_suspicious_taxa(
#'     df = tax_df,
#'     target_col = "Phylum",
#'     check_multi_word = TRUE
#' )
#'
#' phylum_check$report
#'
#'
#' @importFrom stringr str_detect str_count str_remove_all
#' @importFrom dplyr arrange desc
#' @export


#' @importFrom stringr str_detect str_count str_remove_all str_squish
#' @importFrom dplyr arrange desc
#' @export


find_suspicious_taxa <- function(
    df,
    target_col,
    min_letters = 5,
    check_multi_word = FALSE,
    export_csv = FALSE,
    file = NULL
) {
  
  # Validate inputs

  if (!base::is.data.frame(df)) {
    stop("'df' must be a dataframe.")
  }
  
  
  if (!target_col %in% base::names(df)) {
    stop(base::sprintf("Column '%s' not found.", target_col))
  }
  
  
  # Extract unique taxa

  taxa <- unique(df[[target_col]])
  
  taxa <- taxa[!is.na(taxa)]
  
  
  # Calculate QC features

  has_numbers <- stringr::str_detect(taxa, "\\d")
  
  
  letter_count <- stringr::str_count(
    stringr::str_remove_all(taxa, "[^A-Za-z]")
  )
  
  
  word_count <- stringr::str_count(
    stringr::str_squish(taxa), "\\S+"
  )
  
  
  has_multi_word <- word_count > 1
  
  

  # Apply QC rules

  suspicious <- 
    has_numbers |
    letter_count < min_letters
  
  
  # Add multi-word rule if requested
  
  if (isTRUE(check_multi_word)) {
    
    suspicious <- suspicious |
      has_multi_word
    
  }
  

  # Create report

  
  report <- data.frame(
    Taxon = taxa[suspicious],
    Has_Numbers = has_numbers[suspicious],
    Letter_Count = letter_count[suspicious],
    Word_Count = word_count[suspicious],
    Has_Multiple_Words = has_multi_word[suspicious],
    stringsAsFactors = FALSE
  )
  
  
  
  if (nrow(report) > 0) {
    
    report <- report |>
      dplyr::arrange(
        dplyr::desc(Has_Numbers),
        dplyr::desc(Has_Multiple_Words),
        Letter_Count,
        Taxon
      )
    
  }
  

  # Optional export

  
  if (isTRUE(export_csv)) {
    
    if (is.null(file)) {
      
      file <- paste0(target_col, "_suspicious_taxa.csv")
      
    }
    
    utils::write.csv(
      report,
      file = file,
      row.names = FALSE
    )
    
  }
  

  # output message

  
  message(
    sprintf(
      "%d suspicious unique taxa identified in '%s'.",
      nrow(report),
      target_col)
  )
  
  
  
  return(
    list(
      taxa = report$Taxon,
      report = report
    )
  )
  
}


# Usage example
#
# genus_check <- find_suspicious_taxa(
#   df = df_trimmed,
#   target_col = "Genus",
#   min_letters = 6,
#   check_multi_word = TRUE,  # set to FALSE for sp and maybe genus
#   export_csv = TRUE,
#   file = "genus_suspicious_taxa.csv"
# )





# Clean Taxonomy Names with Change Log ######################################## 


#' This function recodes taxonomy names within one or multiple dataframe
#' columns using user-defined replacement mappings. It accepts either a single
#' named character vector for one column or a named list of vectors for
#' multiple columns.
#'
#' The function returns the updated dataframe and a log showing which
#' taxonomy names were changed.
#'
#' @param df The dataframe containing taxonomy information.
#'
#' @param mappings A named character vector for a single column, or a named
#' list of named character vectors for multiple columns.
#'
#' For a single column:
#' \code{c("Old_name" = "New_name")}
#'
#' For multiple columns:
#' \code{list(Column_name = c("Old_name" = "New_name"))}
#'
#' @param column Character string giving the target column when supplying
#' a single mapping vector. Ignored when \code{mappings} is supplied as a
#' named list.
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{The dataframe containing updated taxonomy names.}
#'
#'   \item{report}{A dataframe containing:
#'   \itemize{
#'     \item Column: taxonomy column modified
#'     \item Original: original taxonomy name
#'     \item Updated: replacement taxonomy name
#'   }}
#' }
#'
#' @examples
#'
#' # Single column example
#'
#' phylum_mapping <- c(
#'   "Bacteroidetes" = "Bacteroidota",
#'   "Firmicutes" = "Bacillota"
#' )
#'
#' result <- clean_taxonomy_with_log(
#'   df = tax_df,
#'   mappings = phylum_mapping,
#'   column = "Phylum"
#' )
#'
#'
#' # Multiple column example
#'
#' taxonomy_mappings <- list(
#'   Phylum = c(
#'     "Bacteroidetes" = "Bacteroidota"
#'   ),
#'   Genus = c(
#'     "Lactobacillus" = "Lacticaseibacillus"
#'   )
#' )
#'
#' result <- clean_taxonomy_with_log(
#'   df = tax_df,
#'   mappings = taxonomy_mappings
#' )
#'
#' @importFrom dplyr recode bind_rows distinct
#' @export

clean_taxonomy_with_log <- function(
    df,
    mappings,
    column = NULL
) {
  
  # Validate dataframe
  
  if (!base::is.data.frame(df)) {
    
    stop("'df' must be a dataframe.")
    
  }
  
  
  # Validate mappings
  
  if (base::is.atomic(mappings)) {
    
    if (!base::is.character(mappings) ||
        base::is.null(base::names(mappings))) {
      
      stop(
        "'mappings' must be either a named character vector ",
        "or a named list of named character vectors."
      )
      
    }
    
    if (base::is.null(column)) {
      
      stop(
        "Argument 'column' must be supplied when using ",
        "a single mapping vector."
      )
      
    }
    
    mappings <- list(
      .single_column = mappings
    )
    
    single_column <- TRUE
    
  } else if (base::is.list(mappings)) {
    
    if (base::is.null(base::names(mappings))) {
      
      stop(
        "'mappings' must be a named list."
      )
      
    }
    
    invalid_mapping <- base::vapply(
      
      mappings,
      
      function(x) {
        
        !base::is.character(x) ||
          base::is.null(base::names(x))
        
      },
      
      logical(1)
      
    )
    
    if (base::any(invalid_mapping)) {
      
      stop(
        "Each element of 'mappings' must be a ",
        "named character vector."
      )
      
    }
    
    single_column <- FALSE
    
  } else {
    
    stop(
      "'mappings' must be either a named character vector ",
      "or a named list of named character vectors."
    )
    
  }
  
  
  # Store output and change log
  
  out <- df
  
  change_log <- list()
  
  
  # Loop through columns
  
  for (col_name in base::names(mappings)) {
    
    # Handle single-column input
    
    if (single_column) {
      
      target_col <- column
      
    } else {
      
      target_col <- col_name
      
    }
    
    
    # Check column exists
    
    if (!target_col %in% base::names(out)) {
      
      warning(
        base::sprintf(
          "Column '%s' not found. Skipping.",
          target_col
        )
      )
      
      next
      
    }
    
    
    # Store original values
    
    old_vals <- out[[target_col]]
    # Apply replacements
    
    out[[target_col]] <- dplyr::recode(
      out[[target_col]],
      !!!mappings[[col_name]]
    )
    
    
    # Identify changes
    
    new_vals <- out[[target_col]]
    
    changed <- !base::is.na(old_vals) &
      old_vals != new_vals
    
    
    if (base::any(changed)) {
      
      change_log[[target_col]] <- base::data.frame(
        
        Column = target_col,
        
        Original = old_vals[changed],
        
        Updated = new_vals[changed],
        
        stringsAsFactors = FALSE
        
      )
      
    }
    
    
    base::message(
      base::sprintf(
        "Column [%s]: %d names updated.",
        target_col,
        base::sum(changed)
      )
    )
    
  }
  
  
  # Combine reports
  
  if (base::length(change_log) > 0) {
    
    report <- dplyr::bind_rows(change_log) |>
      dplyr::distinct()
    
  } else {
    
    report <- base::data.frame(
      
      Column = character(),
      
      Original = character(),
      
      Updated = character(),
      
      stringsAsFactors = FALSE
      
    )
    
  }
  
  
  # Return
  
  return(
    
    list(
      
      data = out,
      
      report = report
      
    )
    
  )
  
}


# ### usage examples
# ### Single column to update
# 
# result <- clean_taxonomy_with_log(
#   df = df_trimmed,
#   mappings = genus_map,
#   column = "Genus"
# )
# 
# genus_clean <- result$data
# 
# result$report
# 
# ### Multiple columns to update
# 
# taxonomy_maps <- list(
#   
#   Phylum = c(
#     "Bacteroidetes" = "Bacteroidota"
#   ),
#   
#   Genus = c(
#     "Lactobacillus" = "Lacticaseibacillus",
#     "Clostridium" = "Clostridioides"
#   )
#   
# )
# 
# 
# result <- clean_taxonomy_with_log(
#   df = df_trimmed,
#   mappings = taxonomy_maps
# )

# Compare Two Vectors Side-by-Side ############################################################################### 

#' Compare Two Vectors Side-by-Side
#'
#' This is useful for checking if your "before and after" lists match 
#' up after a filtering step or name changes
#' Takes two vectors of potentially different lengths, displays their 
#' length ratio, and aligns them into a single data frame for easy comparison.
#'
#' @param vec1 The first (original) vector.
#' @param vec2 The second (updated) vector.
#'
#' @return A data frame with two columns: 'Original' and 'Updated'. 
#'   The shorter vector is padded with NAs to match the longer one.
#' @export
#' 
compare_vectors <- function(vec1, vec2) {
  
  # Print the length comparison in X:Y format
  base::cat(base::length(vec1), ":", base::length(vec2), sep = "", "\n\n")
  
  # Equalise lengths to the maximum found to prevent recycling errors
  max_len <- base::max(base::length(vec1), base::length(vec2))
  
  v1_aligned <- vec1
  base::length(v1_aligned) <- max_len
  
  v2_aligned <- vec2
  base::length(v2_aligned) <- max_len
  
  # Create the data frame to inspect the vectors
  df <- base::data.frame(
    Original = v1_aligned,
    Updated = v2_aligned,
    stringsAsFactors = FALSE
  )
  
  return(df)
}


# Example usage
# phylum_comparison <- compare_vectors(phylum_list, phylum_list_uptd)
# 
## To view the new list
# print(phylum_comparison)


# Process Taxonomy Lookup and Metadata Appending  ############################################################################### 


#' application of lineage retrieved from LPSN API to original taxonomy file
#' 
#' 
#' Prerequisites:
#' This function requires the 'openxlsx' package for reading master Excel files.
#' install.packages("openxlsx")
#'
#' Input Requirements:
#' 1. Master Data: Must contain a 'search_term' column used as the lookup key.
#' 2. Taxonomy Data: A CSV where the first column is the ASV/OTU ID.

#' Process Taxonomy Lookup and Metadata Appending
#'
#' This R function takes a taxonomy file and cross-references it against a master 
#' reference file that is the result of the API calls to the LPSN database. It searches through 
#' specified taxonomic ranks (e.g., Species, Genus) until it finds a match in the master list, 
#' then appends the corresponding data.  It will produce a file with the ASV_IDs,
#' the original taxonomy, the field that found a match in the LPSN database and the lineage 
#' returned from LPSN 
#'
#' @param master_path String. Path to the master reference file (.csv, .xlsx, or .xls). 
#'   Must contain a 'search_term' column.
#' @param taxonomy_path String. Path to the input taxonomy CSV file.
#' @param output_path String. Path where the resulting CSV will be saved.
#' @param check_ranks Character vector. The column names in the taxonomy file to 
#'   search through, ordered by priority (usually specific to general).
#' @param append_cols Character vector. The columns from the master file to 
#'   append to the results.
#'
#' @return Writes a CSV file to the specified \code{output_path} and prints a success message.
#' @export
#'
#' @examples
#' \dontrun{
#' ### Define file paths
#' master_file <- "LPSN_result_Master.xlsx"
#' tax_file    <- "orifinal_taxonomy.csv"
#' result_file <- "Final_Taxonomy_Matched.csv"
#' 
#' ### Define Custom Ranks/Columns (Optional)
#' ### LPSN returns kingdom but you do not need to append it to your file
#' # You can provide these in any case (lower, UPPER, Sentence)
#' # ranks in original taxonomy to check
#' my_check_ranks <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")
#' # columns from LPSN results to append to your file
#' my_append_cols <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# 
#' #Run the Function ---
#' process_taxonomy_lookup(
#'   master_path = master_file,
#'   taxonomy_path = tax_file,
#'   output_path = result_file,
#'   check_ranks = my_check_ranks, 
#'   append_cols = my_append_cols  
#' )
#' }
#' 
process_taxonomy_lookup <- function(master_path, 
                                    taxonomy_path, 
                                    output_path, 
                                    check_ranks = c("Species", "Genus", 
                                                    "Family", "Order", 
                                                    "Class", "Phylum", 
                                                    "Kingdom"),
                                    append_cols = c("domain", "kingdom", 
                                                    "phylum", "class", 
                                                    "order", "family", 
                                                    "genus", "species")) {
  
  # Load Master Data based on file extension
  # Detects file extension to use the correct loading function
  ext <- base::tolower(tools::file_ext(master_path))
  
  if (ext == "csv") {
    master <- utils::read.csv(master_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    master <- openxlsx::read.xlsx(master_path)
  } else {
    base::stop("Master file must be .csv, .xlsx, or .xls")
  }
  
  # Load Taxonomy Data
  # Loads the target file; assumes the first column is the row identifier
  taxonomy <- utils::read.csv(taxonomy_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Column Mapping: case-insensitive
  # This section ensures that even if columns are named "genus" vs "Genus", they are found.
  
  # Map check_ranks to actual taxonomy column names
  actual_tax_cols <- base::colnames(taxonomy)
  mapped_check_ranks <- base::character()
  for (r in check_ranks) {
    idx <- base::grep(base::paste0("^", r, "$"), actual_tax_cols, ignore.case = TRUE)
    if (base::length(idx) > 0) mapped_check_ranks <- base::c(mapped_check_ranks, actual_tax_cols[idx])
  }
  
  # Map append_cols to actual master column names
  actual_master_cols <- base::colnames(master)
  mapped_append_cols <- base::character()
  for (c in append_cols) {
    idx <- base::grep(base::paste0("^", c, "$"), actual_master_cols, ignore.case = TRUE)
    if (base::length(idx) > 0) mapped_append_cols <- base::c(mapped_append_cols, actual_master_cols[idx])
  }
  
  # Initialize Result Columns
  # Create a column to track which rank/term was successfully matched
  taxonomy$matched_search_term <- "none"
  for (col in mapped_append_cols) {
    taxonomy[[base::paste0("master_", col)]] <- NA_character_
  }
  
  # Row-wise Processing
  # Iterate through every row of the taxonomy file
  # If a match isn't found at the Species level, it moves to Genus, then Family,
  # and so on. This maximizes the amount of data recovered.
  for (i in 1:base::nrow(taxonomy)) {
    # Check each rank (Species, then Genus, etc.) sequentially
    for (rank in mapped_check_ranks) {
      
      current_value <- base::as.character(taxonomy[i, rank])
      
      # Skip if the rank is empty or explicitly labeled 'unclassified'
      if (base::is.na(current_value) || current_value == "" || base::tolower(current_value) == "unclassified") {
        next
      }
      
      # Perform case-insensitive search against the master 'search_term' column
      match_idx <- base::which(base::tolower(master$search_term) == base::tolower(current_value))
      
      if (base::length(match_idx) > 0) {
        first_match <- match_idx[1] # Explicitly take the first hit
        taxonomy$matched_search_term[i] <- current_value
        
        # Copy values from Master to Taxonomy for all requested columns
        for (m_col in mapped_append_cols) {
          taxonomy[i, base::paste0("master_", m_col)] <- master[first_match, m_col]
        }
        # Once a match is found for a row, stop checking broader ranks
        break
      }
    }
  }
  
  # Export 
  utils::write.csv(taxonomy, output_path, row.names = TRUE)
  base::cat(base::sprintf("Done. Results in: %s\n", output_path))
}





# Log your session details to a .csv ################################################################################

log_session <- function(folder_path) {
  # Create filename with date: "session_log_20240522.csv"
  date_stamp <- format(Sys.Date(), "%Y%m%d")
  full_path <- file.path(folder_path, paste0("session_log_", date_stamp, ".csv"))
  
  # Collect Platform and Package info
  platform <- sessioninfo::platform_info()
  plat_string <- paste0("R ", platform$version, " on ", platform$os)
  
  full_log <- sessioninfo::package_info() |>
    as.data.frame() |>
    dplyr::select(package, version = loadedversion, source) |>
    dplyr::mutate(
      run_timestamp = Sys.time(),
      platform = plat_string
    )
  
  # Append to CSV (header only if file is new)
  write.table(
    full_log, 
    file = full_path, 
    sep = ",", 
    append = TRUE, 
    row.names = FALSE, 
    col.names = !file.exists(full_path)
  )
  
  message("Log updated: ", full_path)
}


# usage example
### LPSN_api$log_session("./filepath")





# LPSN search by lineage rank - with error logs ########


lpsn_by_rank_search <- function(
    taxa_input,
    csv_path,
    workbook_path,
    log_path,
    batch_size = 100,
    pause_min = 1,
    pause_max = 3,
    use_timestamp = TRUE) {
  
  
  # load libraries

  library(LPSN)
  library(progress)
  library(openxlsx)
  
  

  # 1. Record process start time and check login details
  
  start_time <-
    base::Sys.time()
  
  
  base::message("")
  base::message("Checking LPSN credentials...")
  base::message("")
  
  
  user <-
    base::Sys.getenv("DSMZ_API_USER")
  
  
  password <-
    base::Sys.getenv("DSMZ_API_PASSWORD")
  
  
  if (
    user == "" ||
    password == ""
  ) {
    
    base::stop(
      "LPSN credentials not found. ",
      "Please set DSMZ_API_USER and DSMZ_API_PASSWORD ",
      "in your .Renviron file."
    )
  }
  
  
  lpsn <-
    base::tryCatch(
      
      LPSN::open_lpsn(
        user,
        password
      ),
      
      error = function(e) {
        
        base::stop(
          "LPSN login failed. ",
          "Please check that DSMZ_API_USER and DSMZ_API_PASSWORD ",
          "are correct and have not expired.\n",
          "LPSN error: ",
          base::conditionMessage(e)
        )
      }
    )
  
  
  base::message(
    "LPSN login successful. Credentials verified."
  )
  
  
  # 2. Validate taxa_input list

  
  if (
    !base::is.list(taxa_input)
  ) {
    
    base::stop(
      "`taxa_input` must be a named list.",
      "\nExample:",
      "\ntaxa_input <- list(",
      "\n  phylum = phylum_list,",
      "\n  class = class_list,",
      "\n  genus = genus_list,",
      "\n  species = species_list",
      "\n)"
    )
  }
  
  
  if (
    base::is.null(
      base::names(taxa_input)
    ) ||
    base::any(
      base::names(taxa_input) == ""
    ) ||
    base::any(
      base::is.na(
        base::names(taxa_input)
      )
    )
  ) {
    
    base::stop(
      "`taxa_input` must have names corresponding to taxonomic ranks."
    )
  }
  
  

  # 3. Taxonomic rank order and validation

  
  RANK_ORDER <-
    c(
      "domain",
      "kingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )
  
  
  supplied_ranks <-
    base::tolower(
      base::names(taxa_input)
    )
  
  
  invalid_ranks <-
    base::setdiff(
      supplied_ranks,
      RANK_ORDER
    )
  
  
  if (
    base::length(invalid_ranks) > 0
  ) {
    
    base::stop(
      "Invalid taxonomic rank(s) supplied in `taxa_input`: ",
      base::paste(
        invalid_ranks,
        collapse = ", "
      ),
      "\nValid ranks are: ",
      base::paste(
        RANK_ORDER,
        collapse = ", "
      )
    )
  }
  
  
  base::names(taxa_input) <-
    supplied_ranks
  
  
  # 4. Create output directories

  if (
    !base::dir.exists(log_path)
  ) {
    
    base::dir.create(
      log_path,
      recursive = TRUE
    )
  }
  
  
  if (
    !base::dir.exists(csv_path)
  ) {
    
    base::dir.create(
      csv_path,
      recursive = TRUE
    )
  }
  
  
  # 5. De-duplication and pre-processing

  
  original_taxa_count <- 0
  
  
  processed_taxa_input <- base::list()
  
  
  duplicate_log_lines <- character()
  
  
  for (
    supplied_rank in base::names(taxa_input)
  ) {
    
    
    taxa_list <- taxa_input[[supplied_rank]]
    
    
    original_count <- base::length(taxa_list)
    
    
    original_taxa_count <-
      original_taxa_count +
      original_count
    
    
    taxa_list <- base::as.character(taxa_list)
    
    
    taxa_list <- base::trimws(taxa_list)
    
    
    taxa_list <- base::gsub("[[:space:]]+",
                            " ",
                            taxa_list)
    
    
    taxa_list <- taxa_list[!base::is.na(taxa_list) &
                             taxa_list != ""]
    
    
    unique_taxa <- base::unique(taxa_list)
    
    
    duplicates <- base::setdiff(taxa_list, unique_taxa)
    
    
    if (
      base::length(duplicates) > 0
    ) {
      
      duplicate_log_lines <- base::c(duplicate_log_lines,
                                     base::paste0("[", supplied_rank, "]"),
                                     base::sort(base::unique(duplicates)),
                                     "")
    } else {
      
      duplicate_log_lines <- base::c(duplicate_log_lines,
                                     base::paste0("[", supplied_rank, "]"),
                                     "No duplicates found.",
                                     "")
    }
    
    processed_taxa_input[[supplied_rank]] <- unique_taxa
    
  }
  
  taxa_input <- processed_taxa_input
  
  total_unique_taxa <- base::sum(base::vapply(
    taxa_input,
    base::length,
    integer(1)
  )
  )
  
  
  base::writeLines(
    duplicate_log_lines,
    base::file.path(
      log_path,
      "deduplicated_taxa.txt"
    )
  )
  
  
  base::cat(
    "\nTotal input taxa: ",
    original_taxa_count,
    "\n",
    sep = ""
  )
  
  
  base::cat(
    "\nTotal unique taxa after preprocessing: ",
    total_unique_taxa,
    "\n",
    sep = ""
  )
  
  

  # 6. Set up cache and workbook path

  
  api_cache <-
    base::list()
  
  
  if (
    use_timestamp
  ) {
    
    ts <-
      base::format(
        base::Sys.time(),
        "%Y%m%d_%H%M"
      )
    
    
    workbook_path <-
      base::gsub(
        "\\.xlsx$",
        base::paste0(
          "_",
          ts,
          ".xlsx"
        ),
        workbook_path
      )
  }
  
  

  # 7. LPSN fields retained in lineage CSV files

  
  export_fields <-
    c(
      "id",
      "category",
      "full_name",
      "authority",
      "proposed_as",
      "lpsn_parent_id",
      "lpsn_correct_name_id",
      "lpsn_taxonomic_status",
      "nomenclatural_status"
    )
  
  

  # 8. Helper: safe_rbind()

  
  safe_rbind <- function(df_list) {
    
    
    if (
      base::length(df_list) == 0
    ) {
      
      return(NULL)
    }
    
    
    all_cols <-
      base::unique(
        base::unlist(
          base::lapply(
            df_list,
            base::names
          )
        )
      )
    
    
    res_list <-
      base::lapply(
        df_list,
        function(df) {
          
          missing_cols <-
            base::setdiff(
              all_cols,
              base::names(df)
            )
          
          
          for (
            col in missing_cols
          ) {
            
            df[[col]] <-
              NA
          }
          
          
          df[
            ,
            all_cols,
            drop = FALSE
          ]
        }
      )
    
    
    base::do.call(
      base::rbind,
      res_list
    )
  }
  
  

  # 9. Helper: retry_api_call()

  
  retry_api_call <- function(
    expr,
    retries = 3,
    wait = 5) {
    
    
    for (
      i in base::seq_len(retries)
    ) {
      
      
      res <-
        base::tryCatch(
          
          expr(),
          
          error = function(e) {
            
            NULL
          }
        )
      
      
      if (
        !base::is.null(res)
      ) {
        
        return(res)
      }
      
      
      if (
        i < retries
      ) {
        
        base::Sys.sleep(
          wait
        )
      }
    }
    
    
    NULL
  }
  
  

  # 10. Helper: get_taxonomic_status()

  
  get_taxonomic_status <- function(df) {
    
    
    if (
      base::is.null(df) ||
      base::nrow(df) == 0 ||
      !"lpsn_taxonomic_status" %in%
      base::names(df)
    ) {
      
      return(NA_character_)
    }
    
    
    status <-
      base::tolower(
        base::trimws(
          base::as.character(
            df$lpsn_taxonomic_status[1]
          )
        )
      )
    
    
    if (
      base::length(status) == 0 ||
      base::is.na(status) ||
      status == ""
    ) {
      
      return(NA_character_)
    }
    
    
    status
  }
  
  

  # 11. Helper: is_placeholder_record()

  
  is_placeholder_record <- function(df) {
    
    
    if (
      base::is.null(df) ||
      base::nrow(df) == 0
    ) {
      
      return(FALSE)
    }
    
    
    proposed <-
      ""
    
    
    nomenclatural <-
      ""
    
    
    if (
      "proposed_as" %in%
      base::names(df)
    ) {
      
      proposed <-
        base::tolower(
          base::trimws(
            base::as.character(
              df$proposed_as[1]
            )
          )
        )
      
      
      if (
        base::is.na(proposed)
      ) {
        
        proposed <-
          ""
      }
    }
    
    
    if (
      "nomenclatural_status" %in%
      base::names(df)
    ) {
      
      nomenclatural <-
        base::tolower(
          base::trimws(
            base::as.character(
              df$nomenclatural_status[1]
            )
          )
        )
      
      
      if (
        base::is.na(nomenclatural)
      ) {
        
        nomenclatural <-
          ""
      }
    }
    
    
    base::grepl(
      "lpsn placeholder",
      proposed,
      fixed = TRUE
    ) ||
      base::grepl(
        "artificial parent node of taxa incertae sedis",
        nomenclatural,
        fixed = TRUE
      )
  }
  
  

  # 12. Helper: classify_lineage_record()
  #
  # AMENDED: uses startsWith() instead of exact equality
  # so that status strings with trailing qualifiers, e.g.
  #   "correct name (and explicitly recommended for medical use)"
  #   "preferred name (not correct name)"
  # are still correctly classified.

  
  classify_lineage_record <- function(df) {
    
    
    if (
      base::is.null(df) ||
      base::nrow(df) == 0
    ) {
      
      return(NA_character_)
    }
    
    
    if (
      is_placeholder_record(df)
    ) {
      
      return("LPSN Placeholder")
    }
    
    
    status <-
      get_taxonomic_status(
        df
      )
    
    
    if (
      !base::is.na(status) &&
      base::startsWith(status, "correct name")
    ) {
      
      return("Correct Name")
    }
    
    
    if (
      !base::is.na(status) &&
      base::startsWith(status, "pro-correct name")
    ) {
      
      return("Pro-correct Name")
    }
    
    
    if (
      !base::is.na(status) &&
      base::startsWith(status, "preferred name")
    ) {
      
      return("Preferred Name")
    }
    
    
    NA_character_
  }
  
  

  # 13. Helper: lineage_status_priority()

  
  lineage_status_priority <- function(df) {
    
    
    status <-
      classify_lineage_record(
        df
      )
    
    
    if (
      !base::is.na(status) &&
      status == "Correct Name"
    ) {
      
      return(1)
    }
    
    
    if (
      !base::is.na(status) &&
      status == "Pro-correct Name"
    ) {
      
      return(2)
    }
    
    
    if (
      !base::is.na(status) &&
      status == "Preferred Name"
    ) {
      
      return(3)
    }
    
    
    if (
      !base::is.na(status) &&
      status == "LPSN Placeholder"
    ) {
      
      return(5)
    }
    
    
    999
  }
  
  

  # 14. Helper: get_correct_name_id()

  
  get_correct_name_id <- function(df) {
    
    
    if (
      base::is.null(df) ||
      base::nrow(df) == 0 ||
      !"lpsn_correct_name_id" %in%
      base::names(df)
    ) {
      
      return(NULL)
    }
    
    
    id <-
      df$lpsn_correct_name_id[1]
    
    
    if (
      base::length(id) == 0 ||
      base::is.na(id) ||
      id == 0
    ) {
      
      return(NULL)
    }
    
    
    if (
      "id" %in%
      base::names(df) &&
      !base::is.na(df$id[1]) &&
      base::as.character(id) ==
      base::as.character(df$id[1])
    ) {
      
      return(NULL)
    }
    
    
    id
  }
  
  

  # 15. Helper: normalize_taxon_name()

  
  normalize_taxon_name <- function(x) {
    
    
    if (
      base::length(x) == 0 ||
      base::is.na(x)
    ) {
      
      return(NA_character_)
    }
    
    
    x <-
      base::trimws(
        base::as.character(x)
      )
    
    
    x <-
      base::gsub(
        "[[:space:]]+",
        " ",
        x
      )
    
    
    base::tolower(x)
  }
  
  

  # 16. Helper: fetch_record_by_id()

  
  fetch_record_by_id <- function(
    target_id,
    target_category = NULL) {
    
    
    if (
      base::length(target_id) == 0 ||
      base::is.na(target_id) ||
      target_id == 0
    ) {
      
      return(NULL)
    }
    
    
    id_key <-
      base::as.character(
        target_id
      )
    
    
    if (
      id_key %in%
      base::names(api_cache)
    ) {
      
      rec_df <-
        api_cache[[id_key]]
      
    } else {
      
      base::Sys.sleep(
        stats::runif(
          1,
          pause_min,
          pause_max
        )
      )
      
      
      rec_res <-
        retry_api_call(
          function() {
            
            LPSN::fetch(
              lpsn,
              ids = target_id
            )
          }
        )
      
      
      if (
        base::is.null(rec_res)
      ) {
        
        return(NULL)
      }
      
      
      fetch_count <-
        if (
          !base::is.null(rec_res$count) &&
          base::length(rec_res$count) > 0
        ) {
          
          rec_res$count
          
        } else {
          
          NA_integer_
        }
      
      
      if (
        !base::is.na(fetch_count) &&
        fetch_count == 0
      ) {
        
        return(NULL)
      }
      
      
      rec_df <-
        base::as.data.frame(
          rec_res
        )
      
      
      if (
        base::nrow(rec_df) == 0
      ) {
        
        return(NULL)
      }
      
      
      rec_df <-
        rec_df[
          1,
          ,
          drop = FALSE
        ]
      
      
      api_cache[[id_key]] <-
        rec_df
    }
    
    
    if (
      base::is.null(rec_df) ||
      base::nrow(rec_df) == 0
    ) {
      
      return(NULL)
    }
    
    
    if (
      !base::is.null(target_category) &&
      "category" %in%
      base::names(rec_df)
    ) {
      
      if (
        base::is.na(rec_df$category[1]) ||
        base::tolower(
          base::as.character(
            rec_df$category[1]
          )
        ) !=
        base::tolower(
          target_category
        )
      ) {
        
        return(NULL)
      }
    }
    
    
    record_type <-
      classify_lineage_record(
        rec_df
      )
    
    
    if (
      !base::is.na(record_type)
    ) {
      
      return(
        list(
          record = rec_df,
          record_type = record_type,
          method = "Correct-name ID"
        )
      )
    }
    
    
    NULL
  }
  
  

  # 17. Helper: select_best_lineage_record()
  #
  # AMENDED:
  #   - candidatus parameter and matching removed.
  #   - If exact-name candidates exist but none have an
  #     acceptable taxonomic status (e.g. status =
  #     "synonym (and no standing)"), attempt to redirect
  #     via lpsn_correct_name_id before giving up.

  
  select_best_lineage_record <- function(
    candidate_df,
    target_name,
    target_category) {
    
    
    if (
      base::is.null(candidate_df) ||
      base::nrow(candidate_df) == 0
    ) {
      
      return(
        list(
          record = NULL,
          record_type = NULL,
          method = NULL,
          search_record = NULL
        )
      )
    }
    
    
    df <-
      candidate_df
    
    
    if (
      !"category" %in%
      base::names(df)
    ) {
      
      return(
        list(
          record = NULL,
          record_type = NULL,
          method = NULL,
          search_record = NULL
        )
      )
    }
    
    
    df <-
      df[
        !base::is.na(df$category) &
          base::tolower(
            base::as.character(
              df$category
            )
          ) ==
          base::tolower(
            target_category
          ),
        ,
        drop = FALSE
      ]
    
    
    if (
      base::nrow(df) == 0
    ) {
      
      return(
        list(
          record = NULL,
          record_type = NULL,
          method = NULL,
          search_record = NULL
        )
      )
    }
    
    
    if (
      !"full_name" %in%
      base::names(df)
    ) {
      
      return(
        list(
          record = NULL,
          record_type = NULL,
          method = NULL,
          search_record = NULL
        )
      )
    }
    
    
    normalized_target <-
      normalize_taxon_name(
        target_name
      )
    
    
    normalized_names <-
      base::vapply(
        df$full_name,
        normalize_taxon_name,
        character(1)
      )
    
    
    exact_idx <-
      base::which(
        !base::is.na(normalized_names) &
          normalized_names ==
          normalized_target
      )
    
    
    if (
      base::length(exact_idx) == 0
    ) {
      
      return(
        list(
          record = NULL,
          record_type = NULL,
          method = NULL,
          search_record = NULL
        )
      )
    }
    
    
    exact_df <-
      df[
        exact_idx,
        ,
        drop = FALSE
      ]
    
    
    priorities <-
      base::vapply(
        base::seq_len(
          base::nrow(exact_df)
        ),
        function(i) {
          
          lineage_status_priority(
            exact_df[
              i,
              ,
              drop = FALSE
            ]
          )
        },
        numeric(1)
      )
    
    
    acceptable <-
      base::which(
        priorities < 999
      )
    
    

    # No acceptable status among exact-name matches.
    # Try redirecting each candidate via its
    # lpsn_correct_name_id before giving up entirely.
    # This handles the case where the searched name is
    # itself a synonym, e.g.:
    #   "Prochlorophyceae" -> status "synonym (and no
    #   standing)" -> correct name "Chroococcophyceae"

    
    if (
      base::length(acceptable) == 0
    ) {
      
      
      for (
        i in base::seq_len(
          base::nrow(exact_df)
        )
      ) {
        
        
        candidate_for_redirect <-
          exact_df[
            i,
            ,
            drop = FALSE
          ]
        
        
        redirected_id <-
          get_correct_name_id(
            candidate_for_redirect
          )
        
        
        if (
          base::is.null(redirected_id)
        ) {
          
          next
        }
        
        
        redirected <-
          fetch_record_by_id(
            target_id = redirected_id,
            target_category = target_category
          )
        
        
        if (
          !base::is.null(redirected)
        ) {
          
          return(
            list(
              record = redirected$record,
              record_type = redirected$record_type,
              method = "Synonym Redirect",
              search_record = candidate_for_redirect
            )
          )
        }
      }
      
      
      return(
        list(
          record = NULL,
          record_type = NULL,
          method = NULL,
          search_record = exact_df[
            1,
            ,
            drop = FALSE
          ]
        )
      )
    }
    
    
    best_idx <-
      acceptable[
        base::which.min(
          priorities[
            acceptable
          ]
        )
      ]
    
    
    best_record <-
      exact_df[
        best_idx,
        ,
        drop = FALSE
      ]
    
    
    best_type <-
      classify_lineage_record(
        best_record
      )
    
    
    list(
      record = best_record,
      record_type = best_type,
      method = "Direct Search",
      search_record = best_record
    )
  }
  
  

  # 18. Helper: perform_lpsn_search()
  #
  # AMENDED: candidatus parameter and branching removed, as at
  # this juncture - no API call available for Candidatus taxa.
  # Search order is always: Advanced -> Flexible.

  
  perform_lpsn_search <- function(
    search_name,
    search_category,
    process_label = NULL) {
    
    
    base::Sys.sleep(
      stats::runif(
        1,
        pause_min,
        pause_max
      )
    )
    
    

    # Search 1: Advanced search

    
    raw_res <-
      retry_api_call(
        function() {
          
          LPSN::retrieve(
            lpsn,
            search = "advanced",
            `taxon-name` = search_name,
            category = search_category
          )
        }
      )
    
    
    if (
      !base::is.null(raw_res)
    ) {
      
      raw_df <-
        base::as.data.frame(
          raw_res
        )
      
      
      if (
        base::nrow(raw_df) > 0
      ) {
        
        resolved <-
          select_best_lineage_record(
            candidate_df = raw_df,
            target_name = search_name,
            target_category = search_category
          )
        
        
        if (
          !base::is.null(
            resolved$record
          )
        ) {
          
          resolved$process <-
            process_label
          
          
          resolved$search_method <-
            "Advanced Search"
          
          
          return(resolved)
        }
      }
    }
    
    

    # Search 2: Flexible search

    
    base::Sys.sleep(
      stats::runif(
        1,
        pause_min,
        pause_max
      )
    )
    
    
    flexible_res <-
      retry_api_call(
        function() {
          
          LPSN::retrieve(
            lpsn,
            search = "flexible",
            query = list(
              full_name = search_name,
              category = search_category
            )
          )
        }
      )
    
    
    if (
      !base::is.null(flexible_res)
    ) {
      
      flexible_df <-
        base::as.data.frame(
          flexible_res
        )
      
      
      if (
        base::nrow(flexible_df) > 0
      ) {
        
        resolved <-
          select_best_lineage_record(
            candidate_df = flexible_df,
            target_name = search_name,
            target_category = search_category
          )
        
        
        if (
          !base::is.null(
            resolved$record
          )
        ) {
          
          resolved$process <-
            process_label
          
          
          resolved$search_method <-
            "Flexible Search"
          
          
          return(resolved)
        }
      }
    }
    
    
    NULL
  }
  
  

  # 19. Helper: process_single_taxon()
  #
  # Searches the requested taxonomic level and then walks
  # ONLY upward through lpsn_parent_id.
  #
  # AMENDED: candidatus argument removed from the search
  # call; lineage_fallback_summary now reflects whether
  # the starting taxon was resolved via a synonym redirect.

  
  process_single_taxon <- function(
    taxon,
    search_rank) {
    
    

    # Direct search at requested rank

    
    result <-
      perform_lpsn_search(
        search_name = taxon,
        search_category = search_rank,
        process_label = "Direct Search"
      )
    
    

    # Starting taxon not found
  
    
    if (
      base::is.null(result)
    ) {
      
      return(
        list(
          status = "NOT FOUND",
          search_term = taxon,
          search_rank = search_rank,
          record = NULL,
          search_record = NULL,
          lineage_list = list(),
          lineage_status = "Starting taxon not found",
          match_type = NA_character_,
          corrected_name = NA_character_,
          nomenclatural_status = NA_character_,
          lineage_record_type_summary = NA_character_,
          lineage_fallback_summary = NA_character_,
          search_method = "Not Found"
        )
      )
    }
    
    

    # Accepted starting record

    
    first_hit <-
      result$record
    
    
    first_record_type <-
      result$record_type
    
    
    current_id <-
      first_hit$id[1]
    
    
    corrected_name <-
      if (
        "full_name" %in%
        base::names(first_hit) &&
        base::length(first_hit$full_name) > 0
      ) {
        
        first_hit$full_name[1]
        
      } else {
        
        NA_character_
      }
    
    
    search_record <-
      if (
        !base::is.null(
          result$search_record
        )
      ) {
        
        result$search_record
        
      } else {
        
        first_hit
      }
    
    
    match_type <-
      result$method
    
    
    search_rank_index <-
      base::match(
        search_rank,
        RANK_ORDER
      )
    
    
    if (
      base::is.na(search_rank_index)
    ) {
      
      base::stop(
        base::sprintf(
          "Unable to determine valid taxonomic rank for '%s'.",
          search_rank
        )
      )
    }
    
    

    # Walk upward through parent IDs

    
    lineage_list <-
      base::list()
    
    
    lineage_types <-
      character()
    
    
    last_rank_index <-
      search_rank_index
    
    
    while (
      !base::is.null(current_id) &&
      !base::is.na(current_id) &&
      current_id != 0
    ) {
      
      
      rec_df <-
        NULL
      
      
      id_key <-
        base::as.character(
          current_id
        )
      
      

      # Retrieve from cache or LPSN

      
      if (
        id_key %in%
        base::names(api_cache)
      ) {
        
        rec_df <-
          api_cache[[id_key]]
        
      } else {
        
        base::Sys.sleep(
          stats::runif(
            1,
            pause_min,
            pause_max
          )
        )
        
        
        rec_res <-
          retry_api_call(
            function() {
              
              LPSN::fetch(
                lpsn,
                ids = current_id
              )
            }
          )
        
        
        if (
          base::is.null(rec_res)
        ) {
          
          break
        }
        
        
        fetch_count <-
          if (
            !base::is.null(rec_res$count) &&
            base::length(rec_res$count) > 0
          ) {
            
            rec_res$count
            
          } else {
            
            NA_integer_
          }
        
        
        if (
          !base::is.na(fetch_count) &&
          fetch_count == 0
        ) {
          
          break
        }
        
        
        rec_df <-
          base::as.data.frame(
            rec_res
          )
        
        
        if (
          base::nrow(rec_df) == 0
        ) {
          
          break
        }
        
        
        rec_df <-
          rec_df[
            1,
            ,
            drop = FALSE
          ]
        
        
        if (
          "id" %in%
          base::names(rec_df) &&
          base::length(rec_df$id) > 0 &&
          !base::is.na(rec_df$id[1])
        ) {
          
          cache_id <-
            base::as.character(
              rec_df$id[1]
            )
          
          
          api_cache[[cache_id]] <-
            rec_df
        }
      }
      
      

      # Validate record

      
      if (
        base::is.null(rec_df) ||
        base::nrow(rec_df) == 0
      ) {
        
        break
      }
      
      
      if (
        !"category" %in%
        base::names(rec_df) ||
        base::length(rec_df$category) == 0 ||
        base::is.na(rec_df$category[1])
      ) {
        
        break
      }
      
      
      current_rank <-
        base::as.character(
          rec_df$category[1]
        )
      
      
      current_rank_index <-
        base::match(
          current_rank,
          RANK_ORDER
        )
      
      
      if (
        base::is.na(current_rank_index)
      ) {
        
        break
      }
      
      

      # Never move below requested starting rank

      
      if (
        current_rank_index >
        search_rank_index
      ) {
        
        break
      }
      
      

      # Never move sideways or downward

      
      if (
        base::length(lineage_list) > 0 &&
        current_rank_index >= last_rank_index
      ) {
        
        break
      }
      
      

      # Classify lineage record

      
      lineage_record_type <-
        classify_lineage_record(
          rec_df
        )
      
      

      # Redirect synonyms/non-accepted records to their
      # correct-name ID when available.

      
      if (
        base::is.na(lineage_record_type)
      ) {
        
        redirected_id <-
          get_correct_name_id(
            rec_df
          )
        
        
        if (
          !base::is.null(redirected_id)
        ) {
          
          redirected_record <-
            fetch_record_by_id(
              target_id = redirected_id,
              target_category = current_rank
            )
          
          
          if (
            !base::is.null(
              redirected_record
            )
          ) {
            
            current_id <-
              redirected_id
            
            
            next
          }
        }
        
        
        break
      }
      
      

      # Add accepted lineage record
      #
      # IMPORTANT:
      # Keep this assignment on one line.

      
      if (
        current_rank != "division"
      ) {
        
        lineage_record <-
          rec_df[
            ,
            base::intersect(
              export_fields,
              base::names(rec_df)
            ),
            drop = FALSE
          ]
        
        
        lineage_record$lineage_record_type <-
          lineage_record_type
        
        
        lineage_list[[base::length(lineage_list) + 1]] <- lineage_record
        
        
        lineage_types <-
          base::c(
            lineage_types,
            lineage_record_type
          )
        
        
        last_rank_index <-
          current_rank_index
      }
      
      

      # Move ONLY to parent

      
      if (
        "lpsn_parent_id" %in%
        base::names(rec_df) &&
        base::length(
          rec_df$lpsn_parent_id
        ) > 0 &&
        !base::is.na(
          rec_df$lpsn_parent_id[1]
        ) &&
        rec_df$lpsn_parent_id[1] != 0
      ) {
        
        current_id <-
          rec_df$lpsn_parent_id[1]
        
      } else {
        
        current_id <-
          NULL
      }
    }
    
    

    # 20. Determine lineage outcome

    
    if (
      base::length(lineage_list) == 0
    ) {
      
      lineage_status <-
        "Starting taxon resolved but no lineage records retained"
      
    } else {
      
      resolved_rank_indices <-
        base::vapply(
          lineage_list,
          function(x) {
            
            base::match(
              x$category[1],
              RANK_ORDER
            )
          },
          numeric(1)
        )
      
      
      resolved_rank_indices <-
        resolved_rank_indices[
          !base::is.na(
            resolved_rank_indices
          )
        ]
      
      
      if (
        base::length(resolved_rank_indices) == 0
      ) {
        
        lineage_status <-
          "Starting taxon resolved but no recognised lineage ranks retained"
        
      } else if (
        base::min(
          resolved_rank_indices
        ) == 1
      ) {
        
        lineage_status <-
          "Complete to domain"
        
      } else {
        
        lineage_status <-
          "Partial lineage / orphan"
      }
    }
    
    

    # 21. Summarise lineage record types

    
    if (
      base::length(lineage_types) > 0
    ) {
      
      lineage_type_order <-
        c(
          "Correct Name",
          "Pro-correct Name",
          "Preferred Name",
          "LPSN Placeholder"
        )
      
      
      lineage_types <-
        base::unique(
          lineage_types
        )
      
      
      lineage_types <-
        lineage_types[
          lineage_types %in%
            lineage_type_order
        ]
      
      
      lineage_types <-
        lineage_types[
          base::order(
            base::match(
              lineage_types,
              lineage_type_order
            )
          )
        ]
      
      
      lineage_record_type_summary <-
        base::paste(
          lineage_types,
          collapse = " + "
        )
      
    } else {
      
      lineage_record_type_summary <-
        NA_character_
    }
    
    

    # 22. Summarise fallback
    #
    # AMENDED: candidatus branch removed. Now reflects
    # whether the starting taxon was resolved via a
    # synonym redirect (lpsn_correct_name_id), which is
    # the only fallback mechanism remaining at this stage.

    
    lineage_fallback_summary <-
      if (
        base::identical(
          match_type,
          "Synonym Redirect"
        )
      ) {
        
        "Synonym Redirect"
        
      } else {
        
        "No"
      }
    
    

    # 23. Return processed taxon

    
    list(
      status = "FOUND",
      search_term = taxon,
      search_rank = search_rank,
      record = first_hit,
      search_record = search_record,
      lineage_list = lineage_list,
      lineage_status = lineage_status,
      match_type = match_type,
      corrected_name = corrected_name,
      nomenclatural_status =
        if (
          "nomenclatural_status" %in%
          base::names(search_record) &&
          !base::is.na(
            search_record$nomenclatural_status[1]
          )
        ) {
          
          search_record$nomenclatural_status[1]
          
        } else {
          
          "Unknown"
        },
      lineage_record_type_summary =
        lineage_record_type_summary,
      lineage_fallback_summary =
        lineage_fallback_summary,
      search_method = result$search_method
    )
  }
  
  

  # 24. Create processing table

  
  processing_list <-
    base::lapply(
      base::names(taxa_input),
      function(rank_name) {
        
        taxa_values <-
          taxa_input[[rank_name]]
        
        
        if (
          base::length(taxa_values) == 0
        ) {
          
          return(NULL)
        }
        
        
        base::data.frame(
          search_term = taxa_values,
          search_rank = rank_name,
          stringsAsFactors = FALSE
        )
      }
    )
  
  
  processing_list <-
    processing_list[
      !base::vapply(
        processing_list,
        base::is.null,
        logical(1)
      )
    ]
  
  
  if (
    base::length(processing_list) == 0
  ) {
    
    base::stop(
      "No taxa remain after preprocessing."
    )
  }
  
  
  processing_table <-
    base::do.call(
      base::rbind,
      processing_list
    )
  
  
  total_rows <-
    base::nrow(
      processing_table
    )
  
  

  # 25. Search start time

  
  search_start_time <-
    base::Sys.time()
  
  
  base::cat(
    "\nLPSN search started: ",
    base::format(
      search_start_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    "\n\n",
    sep = ""
  )
  
  

  # 26. Progress bar

  
  pb <-
    progress::progress_bar$new(
      format =
        "[:bar] :percent | :current/:total | :taxon ........ :status | Search Avg: :avg | Finish: :finish",
      total = total_rows,
      clear = FALSE,
      width = 120,
      show_after = 0
    )
  
  

  # 27. Main processing loop

  
  batch_number <-
    0
  
  
  num_batches <-
    base::ceiling(
      total_rows / batch_size
    )
  
  
  for (
    batch_start in base::seq(
      1,
      total_rows,
      by = batch_size
    )
  ) {
    
    
    batch_number <-
      batch_number + 1
    
    
    batch_end <-
      base::min(
        batch_start + batch_size - 1,
        total_rows
      )
    
    
    current_batch <-
      processing_table[
        batch_start:batch_end,
        ,
        drop = FALSE
      ]
    
    
    batch_workbook_rows <-
      list()
    
    
    batch_not_found <-
      character()
    
    
    batch_partial_lineages <-
      character()
    
    

    # Process each taxon

    
    for (
      row_index in base::seq_len(
        base::nrow(current_batch)
      )
    ) {
      
      
      taxon <-
        current_batch$search_term[
          row_index
        ]
      
      
      search_rank <-
        current_batch$search_rank[
          row_index
        ]
      
      

      # Process taxon

      
      result <-
        process_single_taxon(
          taxon = taxon,
          search_rank = search_rank
        )
      
      

      # Determine final status shown in progress bar

      
      if (
        result$status == "NOT FOUND"
      ) {
        
        display_status <-
          "NOT FOUND"
        
      } else if (
        result$lineage_status ==
        "Partial lineage / orphan"
      ) {
        
        display_status <-
          "PARTIAL"
        
      } else {
        
        display_status <-
          "SUCCESS"
      }
      
      

      # Calculate search average and estimated completion

      
      completed_taxa <-
        batch_start +
        row_index -
        1
      
      
      elapsed_seconds <-
        base::as.numeric(
          base::difftime(
            base::Sys.time(),
            search_start_time,
            units = "secs"
          )
        )
      
      
      if (
        completed_taxa > 0
      ) {
        
        average_seconds_per_search <-
          elapsed_seconds /
          completed_taxa
        
      } else {
        
        average_seconds_per_search <-
          NA_real_
      }
      
      
      remaining_taxa <-
        total_rows -
        completed_taxa
      
      
      if (
        !base::is.na(
          average_seconds_per_search
        ) &&
        remaining_taxa > 0
      ) {
        
        estimated_completion <-
          base::Sys.time() +
          (
            average_seconds_per_search *
              remaining_taxa
          )
        
      } else {
        
        estimated_completion <-
          base::Sys.time()
      }
      
      

      # Format Search Avg

      
      if (
        base::is.na(
          average_seconds_per_search
        )
      ) {
        
        average_display <-
          "calculating..."
        
      } else if (
        average_seconds_per_search < 60
      ) {
        
        average_display <-
          base::sprintf(
            "%.1fs",
            average_seconds_per_search
          )
        
      } else {
        
        average_minutes <-
          base::floor(
            average_seconds_per_search / 60
          )
        
        
        average_seconds_remainder <-
          base::round(
            average_seconds_per_search -
              (
                average_minutes *
                  60
              )
          )
        
        
        average_display <-
          base::sprintf(
            "%dm %02ds",
            average_minutes,
            average_seconds_remainder
          )
      }
      
      

      # Format estimated completion

      
      finish_display <-
        base::format(
          estimated_completion,
          "%Y-%m-%d %H:%M:%S"
        )
      
      

      # Update the SINGLE progress bar

      
      pb$tick(
        tokens =
          base::list(
            taxon =
              taxon,
            status =
              display_status,
            avg =
              average_display,
            finish =
              finish_display
          )
      )
      
      

      # Not found

      
      if (
        result$status == "NOT FOUND"
      ) {
        
        batch_not_found <-
          base::c(
            batch_not_found,
            base::paste0(
              search_rank,
              "\t",
              taxon
            )
          )
        
        
        next
      }
      
      

      # Partial lineage

      
      if (
        result$lineage_status ==
        "Partial lineage / orphan"
      ) {
        
        batch_partial_lineages <-
          base::c(
            batch_partial_lineages,
            base::paste0(
              search_rank,
              "\t",
              taxon
            )
          )
      }
      
      

      # Build lineage data

      
      lineage_data <-
        if (
          base::length(
            result$lineage_list
          ) > 0
        ) {
          
          safe_rbind(
            result$lineage_list
          )
          
        } else {
          
          NULL
        }
      
      

      # Export individual lineage CSV

      
      if (
        !base::is.null(lineage_data) &&
        base::nrow(lineage_data) > 0
      ) {
        
        
        safe_name <-
          base::gsub(
            "[^[:alnum:]]",
            "_",
            taxon
          )
        
        
        safe_rank <-
          base::gsub(
            "[^[:alnum:]]",
            "_",
            search_rank
          )
        
        
        utils::write.csv(
          lineage_data,
          base::file.path(
            csv_path,
            base::paste0(
              "Batch_",
              batch_number,
              "_",
              safe_rank,
              "_",
              safe_name,
              ".csv"
            )
          ),
          row.names = FALSE
        )
      }
      
      

      # Build Excel row

      
      wb_row <-
        stats::setNames(
          base::as.list(
            base::rep(
              NA,
              base::length(
                RANK_ORDER
              )
            )
          ),
          RANK_ORDER
        )
      
      
      wb_row$search_term <-
        taxon
      
      
      wb_row$search_rank <-
        search_rank
      
      
      wb_row$match_type <-
        result$match_type
      
      
      wb_row$search_method <-
        result$search_method
      
      
      wb_row$nomenclatural_status <-
        result$nomenclatural_status
      
      
      wb_row$corrected_name <-
        result$corrected_name
      
      
      wb_row$lineage_record_type <-
        result$lineage_record_type_summary
      
      
      wb_row$lineage_fallback <-
        result$lineage_fallback_summary
      
      
      wb_row$lineage_status <-
        result$lineage_status
      
      

      # Populate only resolved ranks

      
      if (
        !base::is.null(lineage_data) &&
        base::nrow(lineage_data) > 0
      ) {
        
        
        for (
          i in base::seq_len(
            base::nrow(lineage_data)
          )
        ) {
          
          
          rank <-
            lineage_data$category[
              i
            ]
          
          
          if (
            !base::is.na(rank) &&
            rank %in%
            RANK_ORDER
          ) {
            
            
            rank_index <-
              base::match(
                rank,
                RANK_ORDER
              )
            
            
            search_rank_index <-
              base::match(
                search_rank,
                RANK_ORDER
              )
            
            
            if (
              rank_index <=
              search_rank_index
            ) {
              
              wb_row[[rank]] <-
                lineage_data$full_name[
                  i
                ]
            }
          }
        }
      }
      
      
      batch_workbook_rows[[taxon]] <-
        base::as.data.frame(
          wb_row,
          stringsAsFactors = FALSE
        )
    }
    
    

    # 28. Write batch results to Excel

    
    if (
      base::length(
        batch_workbook_rows
      ) > 0
    ) {
      
      
      final_batch_df <-
        safe_rbind(
          batch_workbook_rows
        )
      
      
      sheet_name <-
        base::paste0(
          "Batch_",
          batch_number
        )
      
      
      wb <-
        if (
          base::file.exists(
            workbook_path
          )
        ) {
          
          openxlsx::loadWorkbook(
            workbook_path
          )
          
        } else {
          
          openxlsx::createWorkbook()
        }
      
      
      if (
        sheet_name %in%
        base::names(wb)
      ) {
        
        openxlsx::removeWorksheet(
          wb,
          sheet_name
        )
      }
      
      
      openxlsx::addWorksheet(
        wb,
        sheetName = sheet_name
      )
      
      
      cols <-
        base::c(
          "search_term",
          "search_rank",
          "match_type",
          "search_method",
          "nomenclatural_status",
          "corrected_name",
          "lineage_record_type",
          "lineage_fallback",
          "lineage_status",
          RANK_ORDER
        )
      
      
      for (
        c in base::setdiff(
          cols,
          base::names(
            final_batch_df
          )
        )
      ) {
        
        final_batch_df[[c]] <-
          NA
      }
      
      
      openxlsx::writeData(
        wb,
        sheet = sheet_name,
        x =
          final_batch_df[
            ,
            cols,
            drop = FALSE
          ]
      )
      
      
      openxlsx::saveWorkbook(
        wb,
        workbook_path,
        overwrite = TRUE
      )
    }
    
    

    # 29. Write not-found log

    
    not_found_path <-
      base::file.path(
        log_path,
        base::paste0(
          "taxa_not_found_batch_",
          batch_number,
          ".txt"
        )
      )
    
    
    if (
      base::length(
        batch_not_found
      ) > 0
    ) {
      
      base::writeLines(
        base::sort(
          batch_not_found
        ),
        not_found_path
      )
      
    } else {
      
      base::writeLines(
        "No taxa failed to resolve at the requested starting rank.",
        not_found_path
      )
    }
    
    

    # 30. Write partial lineage log

    
    partial_path <-
      base::file.path(
        log_path,
        base::paste0(
          "partial_lineages_batch_",
          batch_number,
          ".txt"
        )
      )
    
    
    if (
      base::length(
        batch_partial_lineages
      ) > 0
    ) {
      
      base::writeLines(
        base::sort(
          batch_partial_lineages
        ),
        partial_path
      )
      
    } else {
      
      base::writeLines(
        "No partial/orphan lineages encountered.",
        partial_path
      )
    }
  }
  
  

  # 31. Completion

  
  finish_time <-
    base::Sys.time()
  
  
  duration <-
    base::difftime(
      finish_time,
      search_start_time,
      units = "auto"
    )
  
  
  base::cat(
    "\n\nLPSN Run Complete at: ",
    base::format(
      finish_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    "\n",
    sep = ""
  )
  
  
  base::cat(
    "Total time elapsed: ",
    base::format(
      duration
    ),
    "\n",
    sep = ""
  )
  
  
  base::cat(
    "Total unique taxa processed: ",
    total_unique_taxa,
    "\n",
    sep = ""
  )
  
  
  base::cat(
    "Total batches: ",
    num_batches,
    "\n",
    sep = ""
  )
  
  
  invisible(
    list(
      taxa_input = taxa_input,
      total_unique_taxa = total_unique_taxa,
      total_batches = num_batches,
      workbook_path = workbook_path,
      csv_path = csv_path,
      log_path = log_path
    )
  )
}



#  LPSN search with no lineage rank - with error logs ########



lpsn_no_rank_search <- function(
    taxa_names,
    csv_path,
    workbook_path,
    log_path,
    batch_size = 100,
    pause_min = 1,
    pause_max = 3,
    use_timestamp = TRUE) {
  
  

  # load libraries

  
  library(LPSN)
  library(progress)
  library(openxlsx)
  
  

  # 1. Record process start time and check login details

  
  start_time <- base::Sys.time()
  
  
  base::message("")
  base::message("Checking LPSN credentials...")
  base::message("")
  
  
  user     <- base::Sys.getenv("DSMZ_API_USER")
  password <- base::Sys.getenv("DSMZ_API_PASSWORD")
  
  
  if (user == "" || password == "") {
    
    base::stop(
      "LPSN credentials not found. ",
      "Please set DSMZ_API_USER and DSMZ_API_PASSWORD ",
      "in your .Renviron file."
    )
  }
  
  
  lpsn <-
    base::tryCatch(
      LPSN::open_lpsn(user, password),
      error = function(e) {
        base::stop(
          "LPSN login failed. ",
          "Please check that DSMZ_API_USER and DSMZ_API_PASSWORD ",
          "are correct and have not expired.\n",
          "LPSN error: ",
          base::conditionMessage(e)
        )
      }
    )
  
  
  base::message("LPSN login successful. Credentials verified.")
  
  

  # 2. Validate taxa_names

  
  if (
    !base::is.vector(taxa_names) ||
    base::is.list(taxa_names)
  ) {
    
    base::stop(
      "`taxa_names` must be a simple character vector of taxon names ",
      "(no ranks). Example:\n",
      "taxa_names <- c(\"Bacillus cereus\", \"Lactobacillaceae\", \"Escherichia coli\")"
    )
  }
  
  

  # 3. Taxonomic ranks (used only to bound/validate the
  #    walk once a starting rank is discovered)

  
  RANK_ORDER <-
    c(
      "domain", "kingdom", "phylum", "class",
      "order", "family", "genus", "species"
    )
  
  

  # 4. Create output directories

  
  if (!base::dir.exists(log_path)) {
    base::dir.create(log_path, recursive = TRUE)
  }
  
  if (!base::dir.exists(csv_path)) {
    base::dir.create(csv_path, recursive = TRUE)
  }
  
  

  # 5. De-duplication and preprocessing (flat list, no ranks)

  
  original_taxa_count <- base::length(taxa_names)
  
  
  taxa_names <- base::as.character(taxa_names)
  taxa_names <- base::trimws(taxa_names)
  taxa_names <- base::gsub("[[:space:]]+", " ", taxa_names)
  taxa_names <- taxa_names[!base::is.na(taxa_names) & taxa_names != ""]
  
  
  unique_taxa <- base::unique(taxa_names)
  
  
  duplicates <- base::setdiff(taxa_names, unique_taxa)
  
  
  duplicate_log_lines <-
    if (base::length(duplicates) > 0) {
      base::c("[unranked]", base::sort(base::unique(duplicates)), "")
    } else {
      base::c("[unranked]", "No duplicates found.", "")
    }
  
  
  base::writeLines(
    duplicate_log_lines,
    base::file.path(log_path, "deduplicated_taxa.txt")
  )
  
  
  taxa_names <- unique_taxa
  
  
  total_unique_taxa <- base::length(taxa_names)
  
  
  base::cat("\nTotal input taxa: ", original_taxa_count, "\n", sep = "")
  base::cat("\nTotal unique taxa after preprocessing: ", total_unique_taxa, "\n", sep = "")
  
  

  # 6. Set up cache and workbook path

  
  api_cache <- base::list()
  
  
  if (use_timestamp) {
    
    ts <- base::format(base::Sys.time(), "%Y%m%d_%H%M")
    
    workbook_path <-
      base::gsub("\\.xlsx$", base::paste0("_", ts, ".xlsx"), workbook_path)
  }
  
  

  # 7. LPSN fields retained in lineage CSV files

  
  export_fields <-
    c(
      "id", "category", "full_name", "authority", "proposed_as",
      "lpsn_parent_id", "lpsn_correct_name_id",
      "lpsn_taxonomic_status", "nomenclatural_status"
    )
  
  

  # 8. Helper: safe_rbind()

  
  safe_rbind <- function(df_list) {
    
    if (base::length(df_list) == 0) return(NULL)
    
    all_cols <- base::unique(base::unlist(base::lapply(df_list, base::names)))
    
    res_list <- base::lapply(df_list, function(df) {
      
      missing_cols <- base::setdiff(all_cols, base::names(df))
      
      for (col in missing_cols) {
        df[[col]] <- NA
      }
      
      df[, all_cols, drop = FALSE]
    })
    
    base::do.call(base::rbind, res_list)
  }
  
  

  # 9. Helper: retry_api_call()

  
  retry_api_call <- function(expr, retries = 3, wait = 5) {
    
    for (i in base::seq_len(retries)) {
      
      res <- base::tryCatch(expr(), error = function(e) NULL)
      
      if (!base::is.null(res)) return(res)
      
      if (i < retries) base::Sys.sleep(wait)
    }
    
    NULL
  }
  
  

  # 10. Helper: get_taxonomic_status()

  
  get_taxonomic_status <- function(df) {
    
    if (base::is.null(df) || base::nrow(df) == 0 ||
        !"lpsn_taxonomic_status" %in% base::names(df)) {
      return(NA_character_)
    }
    
    status <- base::tolower(base::trimws(base::as.character(df$lpsn_taxonomic_status[1])))
    
    if (base::length(status) == 0 || base::is.na(status) || status == "") {
      return(NA_character_)
    }
    
    status
  }
  
  

  # 11. Helper: is_placeholder_record()

  
  is_placeholder_record <- function(df) {
    
    if (base::is.null(df) || base::nrow(df) == 0) return(FALSE)
    
    proposed <- ""
    nomenclatural <- ""
    
    if ("proposed_as" %in% base::names(df)) {
      proposed <- base::tolower(base::trimws(base::as.character(df$proposed_as[1])))
      if (base::is.na(proposed)) proposed <- ""
    }
    
    if ("nomenclatural_status" %in% base::names(df)) {
      nomenclatural <- base::tolower(base::trimws(base::as.character(df$nomenclatural_status[1])))
      if (base::is.na(nomenclatural)) nomenclatural <- ""
    }
    
    base::grepl("lpsn placeholder", proposed, fixed = TRUE) ||
      base::grepl("artificial parent node of taxa incertae sedis", nomenclatural, fixed = TRUE)
  }
  
  

  # 12. Helper: classify_lineage_record()
  #
  # AMENDED: uses startsWith() instead of exact equality
  # so that status strings with trailing qualifiers, e.g.
  #   "correct name (and explicitly recommended for medical use)"
  #   "preferred name (not correct name)"
  # are still correctly classified.

  
  classify_lineage_record <- function(df) {
    
    if (base::is.null(df) || base::nrow(df) == 0) return(NA_character_)
    
    if (is_placeholder_record(df)) return("LPSN Placeholder")
    
    status <- get_taxonomic_status(df)
    
    if (!base::is.na(status) && base::startsWith(status, "correct name")) {
      return("Correct Name")
    }
    
    if (!base::is.na(status) && base::startsWith(status, "pro-correct name")) {
      return("Pro-correct Name")
    }
    
    if (!base::is.na(status) && base::startsWith(status, "preferred name")) {
      return("Preferred Name")
    }
    
    NA_character_
  }
  
  

  # 13. Helper: lineage_status_priority()

  
  lineage_status_priority <- function(df) {
    
    status <- classify_lineage_record(df)
    
    if (!base::is.na(status) && status == "Correct Name")     return(1)
    if (!base::is.na(status) && status == "Pro-correct Name") return(2)
    if (!base::is.na(status) && status == "Preferred Name")   return(3)
    if (!base::is.na(status) && status == "LPSN Placeholder") return(5)
    
    999
  }
  
  

  # 14. Helper: get_correct_name_id()

  
  get_correct_name_id <- function(df) {
    
    if (base::is.null(df) || base::nrow(df) == 0 ||
        !"lpsn_correct_name_id" %in% base::names(df)) {
      return(NULL)
    }
    
    id <- df$lpsn_correct_name_id[1]
    
    if (base::length(id) == 0 || base::is.na(id) || id == 0) return(NULL)
    
    if ("id" %in% base::names(df) && !base::is.na(df$id[1]) &&
        base::as.character(id) == base::as.character(df$id[1])) {
      return(NULL)
    }
    
    id
  }
  
  

  # 15. Helper: normalize_taxon_name()

  
  normalize_taxon_name <- function(x) {
    
    if (base::length(x) == 0 || base::is.na(x)) return(NA_character_)
    
    x <- base::trimws(base::as.character(x))
    x <- base::gsub("[[:space:]]+", " ", x)
    
    base::tolower(x)
  }
  
  

  # 16. Helper: fetch_record_by_id()

  
  fetch_record_by_id <- function(target_id, target_category = NULL) {
    
    if (base::length(target_id) == 0 || base::is.na(target_id) || target_id == 0) {
      return(NULL)
    }
    
    id_key <- base::as.character(target_id)
    
    if (id_key %in% base::names(api_cache)) {
      
      rec_df <- api_cache[[id_key]]
      
    } else {
      
      base::Sys.sleep(stats::runif(1, pause_min, pause_max))
      
      rec_res <- retry_api_call(function() LPSN::fetch(lpsn, ids = target_id))
      
      if (base::is.null(rec_res)) return(NULL)
      
      fetch_count <-
        if (!base::is.null(rec_res$count) && base::length(rec_res$count) > 0) {
          rec_res$count
        } else {
          NA_integer_
        }
      
      if (!base::is.na(fetch_count) && fetch_count == 0) return(NULL)
      
      rec_df <- base::as.data.frame(rec_res)
      
      if (base::nrow(rec_df) == 0) return(NULL)
      
      rec_df <- rec_df[1, , drop = FALSE]
      
      api_cache[[id_key]] <- rec_df
    }
    
    if (base::is.null(rec_df) || base::nrow(rec_df) == 0) return(NULL)
    
    if (!base::is.null(target_category) && "category" %in% base::names(rec_df)) {
      
      if (base::is.na(rec_df$category[1]) ||
          base::tolower(base::as.character(rec_df$category[1])) != base::tolower(target_category)) {
        return(NULL)
      }
    }
    
    record_type <- classify_lineage_record(rec_df)
    
    if (!base::is.na(record_type)) {
      return(list(record = rec_df, record_type = record_type, method = "Correct-name ID"))
    }
    
    NULL
  }
  
  

  # 17. Helper: select_best_lineage_record_any_rank()
  #
  # RANK-AGNOSTIC VERSION.
  #
  # AMENDED: if exact-name candidates exist but none have
  # an acceptable taxonomic status (e.g. status =
  # "synonym (and no standing)"), attempt to redirect via
  # lpsn_correct_name_id before giving up. The redirect is
  # tried first at the synonym's own rank, then retried
  # with no rank restriction, since the correct name may
  # sit at a different rank than the synonym itself.

  
  select_best_lineage_record_any_rank <- function(candidate_df, target_name) {
    
    if (base::is.null(candidate_df) || base::nrow(candidate_df) == 0) {
      return(list(record = NULL, record_type = NULL, method = NULL,
                  search_record = NULL, resolved_rank = NULL))
    }
    
    df <- candidate_df
    
    if (!"full_name" %in% base::names(df) || !"category" %in% base::names(df)) {
      return(list(record = NULL, record_type = NULL, method = NULL,
                  search_record = NULL, resolved_rank = NULL))
    }
    
    normalized_target <- normalize_taxon_name(target_name)
    normalized_names   <- base::vapply(df$full_name, normalize_taxon_name, character(1))
    
    exact_idx <- base::which(!base::is.na(normalized_names) & normalized_names == normalized_target)
    
    if (base::length(exact_idx) == 0) {
      return(list(record = NULL, record_type = NULL, method = NULL,
                  search_record = NULL, resolved_rank = NULL))
    }
    
    exact_df <- df[exact_idx, , drop = FALSE]
    
    priorities <-
      base::vapply(
        base::seq_len(base::nrow(exact_df)),
        function(i) lineage_status_priority(exact_df[i, , drop = FALSE]),
        numeric(1)
      )
    
    acceptable <- base::which(priorities < 999)
    
    
    # ....................................................
    # No acceptable status among exact-name matches.
    # Try redirecting each candidate via its
    # lpsn_correct_name_id before giving up entirely.
    # This handles the case where the searched name IS an
    # exact match, but is itself a synonym, e.g.:
    #   "Prochlorophyceae" -> status "synonym (and no
    #   standing)" -> correct name "Chroococcophyceae"
    # ....................................................
    
    if (base::length(acceptable) == 0) {
      
      for (i in base::seq_len(base::nrow(exact_df))) {
        
        candidate_for_redirect <- exact_df[i, , drop = FALSE]
        
        redirected_id <- get_correct_name_id(candidate_for_redirect)
        
        if (base::is.null(redirected_id)) next
        
        own_category <- base::as.character(candidate_for_redirect$category[1])
        
        # Try same-rank redirect first (most common case)
        redirected <- fetch_record_by_id(target_id = redirected_id, target_category = own_category)
        
        # If that fails, retry without category restriction,
        # in case the correct name sits at a different rank
        if (base::is.null(redirected)) {
          redirected <- fetch_record_by_id(target_id = redirected_id, target_category = NULL)
        }
        
        if (!base::is.null(redirected)) {
          
          resolved_rank <-
            base::tolower(base::as.character(redirected$record$category[1]))
          
          return(
            list(
              record        = redirected$record,
              record_type   = redirected$record_type,
              method        = "Synonym Redirect",
              search_record = candidate_for_redirect,
              resolved_rank = resolved_rank
            )
          )
        }
      }
      
      return(list(record = NULL, record_type = NULL, method = NULL,
                  search_record = exact_df[1, , drop = FALSE], resolved_rank = NULL))
    }
    
    
    best_idx    <- acceptable[base::which.min(priorities[acceptable])]
    best_record <- exact_df[best_idx, , drop = FALSE]
    best_type   <- classify_lineage_record(best_record)
    best_rank   <- base::tolower(base::as.character(best_record$category[1]))
    
    list(
      record        = best_record,
      record_type   = best_type,
      method        = "Direct Search",
      search_record = best_record,
      resolved_rank = best_rank
    )
  }
  
  

  # 18. Helper: perform_lpsn_search_any_rank()
  #
  # RANK-AGNOSTIC VERSION. No `category` argument is sent
  # to LPSN::retrieve(). Search order: Advanced -> Flexible.

  
  perform_lpsn_search_any_rank <- function(search_name, process_label = NULL) {
    
    base::Sys.sleep(stats::runif(1, pause_min, pause_max))
    
    # Search 1: Advanced search (no category)
    raw_res <-
      retry_api_call(function() {
        LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = search_name)
      })
    
    if (!base::is.null(raw_res)) {
      
      raw_df <- base::as.data.frame(raw_res)
      
      if (base::nrow(raw_df) > 0) {
        
        resolved <- select_best_lineage_record_any_rank(raw_df, search_name)
        
        if (!base::is.null(resolved$record)) {
          resolved$process       <- process_label
          resolved$search_method <- "Advanced Search"
          return(resolved)
        }
      }
    }
    
    # Search 2: Flexible search (no category)
    base::Sys.sleep(stats::runif(1, pause_min, pause_max))
    
    flexible_res <-
      retry_api_call(function() {
        LPSN::retrieve(lpsn, search = "flexible", query = list(full_name = search_name))
      })
    
    if (!base::is.null(flexible_res)) {
      
      flexible_df <- base::as.data.frame(flexible_res)
      
      if (base::nrow(flexible_df) > 0) {
        
        resolved <- select_best_lineage_record_any_rank(flexible_df, search_name)
        
        if (!base::is.null(resolved$record)) {
          resolved$process       <- process_label
          resolved$search_method <- "Flexible Search"
          return(resolved)
        }
      }
    }
    
    NULL
  }
  
  

  # 19. Helper: process_single_taxon_any_rank()
  #
  # Searches with NO rank specified. The starting rank is
  # taken from resolved$resolved_rank (the category of the
  # winning exact-name match, whether reached directly or
  # via a synonym redirect). From that point onward, the
  # upward walk through lpsn_parent_id is IDENTICAL to the
  # ranked version.
  #
  # AMENDED: lineage_fallback_summary now reflects whether
  # the starting taxon was resolved via a synonym redirect.

  
  process_single_taxon_any_rank <- function(taxon) {
    
    result <- perform_lpsn_search_any_rank(search_name = taxon, process_label = "Direct Search")
    
    if (base::is.null(result) || base::is.null(result$resolved_rank)) {
      
      return(
        list(
          status = "NOT FOUND",
          search_term = taxon,
          search_rank = NA_character_,
          record = NULL,
          search_record = NULL,
          lineage_list = list(),
          lineage_status = "Starting taxon not found",
          match_type = NA_character_,
          corrected_name = NA_character_,
          nomenclatural_status = NA_character_,
          lineage_record_type_summary = NA_character_,
          lineage_fallback_summary = NA_character_,
          search_method = "Not Found"
        )
      )
    }
    
    
    # ....................................................
    # Rank discovered from the search result
    # ....................................................
    
    search_rank <- result$resolved_rank
    
    
    search_rank_index <- base::match(search_rank, RANK_ORDER)
    
    if (base::is.na(search_rank_index)) {
      
      return(
        list(
          status = "NOT FOUND",
          search_term = taxon,
          search_rank = search_rank,
          record = NULL,
          search_record = NULL,
          lineage_list = list(),
          lineage_status = "Resolved rank not recognised",
          match_type = NA_character_,
          corrected_name = NA_character_,
          nomenclatural_status = NA_character_,
          lineage_record_type_summary = NA_character_,
          lineage_fallback_summary = NA_character_,
          search_method = result$search_method
        )
      )
    }
    
    
    first_hit  <- result$record
    current_id <- first_hit$id[1]
    
    corrected_name <-
      if ("full_name" %in% base::names(first_hit) && base::length(first_hit$full_name) > 0) {
        first_hit$full_name[1]
      } else {
        NA_character_
      }
    
    search_record <- if (!base::is.null(result$search_record)) result$search_record else first_hit
    
    match_type <- result$method
    
    
    # ....................................................
    # Walk upward through parent IDs
    # [IDENTICAL LOGIC to the ranked version]
    # ....................................................
    
    lineage_list  <- base::list()
    lineage_types <- character()
    last_rank_index <- search_rank_index
    
    while (!base::is.null(current_id) && !base::is.na(current_id) && current_id != 0) {
      
      rec_df <- NULL
      id_key <- base::as.character(current_id)
      
      if (id_key %in% base::names(api_cache)) {
        
        rec_df <- api_cache[[id_key]]
        
      } else {
        
        base::Sys.sleep(stats::runif(1, pause_min, pause_max))
        
        rec_res <- retry_api_call(function() LPSN::fetch(lpsn, ids = current_id))
        
        if (base::is.null(rec_res)) break
        
        fetch_count <-
          if (!base::is.null(rec_res$count) && base::length(rec_res$count) > 0) {
            rec_res$count
          } else {
            NA_integer_
          }
        
        if (!base::is.na(fetch_count) && fetch_count == 0) break
        
        rec_df <- base::as.data.frame(rec_res)
        
        if (base::nrow(rec_df) == 0) break
        
        rec_df <- rec_df[1, , drop = FALSE]
        
        if ("id" %in% base::names(rec_df) && base::length(rec_df$id) > 0 && !base::is.na(rec_df$id[1])) {
          cache_id <- base::as.character(rec_df$id[1])
          api_cache[[cache_id]] <- rec_df
        }
      }
      
      if (base::is.null(rec_df) || base::nrow(rec_df) == 0) break
      
      if (!"category" %in% base::names(rec_df) ||
          base::length(rec_df$category) == 0 ||
          base::is.na(rec_df$category[1])) {
        break
      }
      
      current_rank       <- base::as.character(rec_df$category[1])
      current_rank_index <- base::match(current_rank, RANK_ORDER)
      
      if (base::is.na(current_rank_index)) break
      
      if (current_rank_index > search_rank_index) break
      
      if (base::length(lineage_list) > 0 && current_rank_index >= last_rank_index) break
      
      lineage_record_type <- classify_lineage_record(rec_df)
      
      if (base::is.na(lineage_record_type)) {
        
        redirected_id <- get_correct_name_id(rec_df)
        
        if (!base::is.null(redirected_id)) {
          
          redirected_record <- fetch_record_by_id(target_id = redirected_id, target_category = current_rank)
          
          if (!base::is.null(redirected_record)) {
            current_id <- redirected_id
            next
          }
        }
        
        break
      }
      
      if (current_rank != "division") {
        
        lineage_record <- rec_df[, base::intersect(export_fields, base::names(rec_df)), drop = FALSE]
        lineage_record$lineage_record_type <- lineage_record_type
        
        lineage_list[[base::length(lineage_list) + 1]] <- lineage_record
        lineage_types <- base::c(lineage_types, lineage_record_type)
        
        last_rank_index <- current_rank_index
      }
      
      if ("lpsn_parent_id" %in% base::names(rec_df) &&
          base::length(rec_df$lpsn_parent_id) > 0 &&
          !base::is.na(rec_df$lpsn_parent_id[1]) &&
          rec_df$lpsn_parent_id[1] != 0) {
        
        current_id <- rec_df$lpsn_parent_id[1]
        
      } else {
        
        current_id <- NULL
      }
    }
    
    
    # ....................................................
    # Determine lineage outcome
    # ....................................................
    
    if (base::length(lineage_list) == 0) {
      
      lineage_status <- "Starting taxon resolved but no lineage records retained"
      
    } else {
      
      resolved_rank_indices <-
        base::vapply(lineage_list, function(x) base::match(x$category[1], RANK_ORDER), numeric(1))
      
      resolved_rank_indices <- resolved_rank_indices[!base::is.na(resolved_rank_indices)]
      
      if (base::length(resolved_rank_indices) == 0) {
        lineage_status <- "Starting taxon resolved but no recognised lineage ranks retained"
      } else if (base::min(resolved_rank_indices) == 1) {
        lineage_status <- "Complete to domain"
      } else {
        lineage_status <- "Partial lineage / orphan"
      }
    }
    
    
    # ....................................................
    # Summarise lineage record types
    # ....................................................
    
    if (base::length(lineage_types) > 0) {
      
      lineage_type_order <- c("Correct Name", "Pro-correct Name", "Preferred Name", "LPSN Placeholder")
      
      lineage_types <- base::unique(lineage_types)
      lineage_types <- lineage_types[lineage_types %in% lineage_type_order]
      lineage_types <- lineage_types[base::order(base::match(lineage_types, lineage_type_order))]
      
      lineage_record_type_summary <- base::paste(lineage_types, collapse = " + ")
      
    } else {
      
      lineage_record_type_summary <- NA_character_
    }
    
    
    # ....................................................
    # Summarise fallback
    #
    # AMENDED: reflects whether the starting taxon was
    # resolved via a synonym redirect.
    # ....................................................
    
    lineage_fallback_summary <-
      if (base::identical(match_type, "Synonym Redirect")) {
        "Synonym Redirect"
      } else {
        "No"
      }
    
    
    # ....................................................
    # Return processed taxon
    # ....................................................
    
    list(
      status = "FOUND",
      search_term = taxon,
      search_rank = search_rank,
      record = first_hit,
      search_record = search_record,
      lineage_list = lineage_list,
      lineage_status = lineage_status,
      match_type = match_type,
      corrected_name = corrected_name,
      nomenclatural_status =
        if ("nomenclatural_status" %in% base::names(search_record) &&
            !base::is.na(search_record$nomenclatural_status[1])) {
          search_record$nomenclatural_status[1]
        } else {
          "Unknown"
        },
      lineage_record_type_summary = lineage_record_type_summary,
      lineage_fallback_summary = lineage_fallback_summary,
      search_method = result$search_method
    )
  }
  
  

  # 20. Create processing table (flat, no rank column input)

  
  processing_table <- base::data.frame(search_term = taxa_names, stringsAsFactors = FALSE)
  
  total_rows <- base::nrow(processing_table)
  
  if (total_rows == 0) {
    base::stop("No taxa remain after preprocessing.")
  }
  
  

  # 21. Search start time

  
  search_start_time <- base::Sys.time()
  
  base::cat(
    "\nLPSN search started: ",
    base::format(search_start_time, "%Y-%m-%d %H:%M:%S"),
    "\n\n", sep = ""
  )
  
  

  # 22. Progress bar

  
  pb <-
    progress::progress_bar$new(
      format = "[:bar] :percent | :current/:total | :taxon ........ :status | Search Avg: :avg | Finish: :finish",
      total = total_rows,
      clear = FALSE,
      width = 120,
      show_after = 0
    )
  
  

  # 23. Main processing loop

  
  batch_number <- 0
  num_batches  <- base::ceiling(total_rows / batch_size)
  
  for (batch_start in base::seq(1, total_rows, by = batch_size)) {
    
    batch_number <- batch_number + 1
    batch_end    <- base::min(batch_start + batch_size - 1, total_rows)
    
    current_batch <- processing_table[batch_start:batch_end, , drop = FALSE]
    
    batch_workbook_rows   <- list()
    batch_not_found        <- character()
    batch_partial_lineages <- character()
    
    for (row_index in base::seq_len(base::nrow(current_batch))) {
      
      taxon <- current_batch$search_term[row_index]
      
      result <- process_single_taxon_any_rank(taxon = taxon)
      
      display_status <-
        if (result$status == "NOT FOUND") {
          "NOT FOUND"
        } else if (result$lineage_status == "Partial lineage / orphan") {
          "PARTIAL"
        } else {
          "SUCCESS"
        }
      
      completed_taxa  <- batch_start + row_index - 1
      elapsed_seconds <- base::as.numeric(base::difftime(base::Sys.time(), search_start_time, units = "secs"))
      
      average_seconds_per_search <-
        if (completed_taxa > 0) elapsed_seconds / completed_taxa else NA_real_
      
      remaining_taxa <- total_rows - completed_taxa
      
      estimated_completion <-
        if (!base::is.na(average_seconds_per_search) && remaining_taxa > 0) {
          base::Sys.time() + (average_seconds_per_search * remaining_taxa)
        } else {
          base::Sys.time()
        }
      
      average_display <-
        if (base::is.na(average_seconds_per_search)) {
          "calculating..."
        } else if (average_seconds_per_search < 60) {
          base::sprintf("%.1fs", average_seconds_per_search)
        } else {
          average_minutes <- base::floor(average_seconds_per_search / 60)
          average_seconds_remainder <- base::round(average_seconds_per_search - (average_minutes * 60))
          base::sprintf("%dm %02ds", average_minutes, average_seconds_remainder)
        }
      
      finish_display <- base::format(estimated_completion, "%Y-%m-%d %H:%M:%S")
      
      pb$tick(tokens = base::list(taxon = taxon, status = display_status, avg = average_display, finish = finish_display))
      
      
      if (result$status == "NOT FOUND") {
        batch_not_found <- base::c(batch_not_found, base::paste0("[unranked]\t", taxon))
        next
      }
      
      
      if (result$lineage_status == "Partial lineage / orphan") {
        batch_partial_lineages <- base::c(batch_partial_lineages, base::paste0(result$search_rank, "\t", taxon))
      }
      
      
      lineage_data <-
        if (base::length(result$lineage_list) > 0) safe_rbind(result$lineage_list) else NULL
      
      
      if (!base::is.null(lineage_data) && base::nrow(lineage_data) > 0) {
        
        safe_name <- base::gsub("[^[:alnum:]]", "_", taxon)
        safe_rank <- base::gsub("[^[:alnum:]]", "_", result$search_rank)
        
        utils::write.csv(
          lineage_data,
          base::file.path(csv_path, base::paste0("Batch_", batch_number, "_", safe_rank, "_", safe_name, ".csv")),
          row.names = FALSE
        )
      }
      
      
      wb_row <- stats::setNames(base::as.list(base::rep(NA, base::length(RANK_ORDER))), RANK_ORDER)
      
      wb_row$search_term    <- taxon
      wb_row$search_rank    <- result$search_rank
      wb_row$match_type     <- result$match_type
      wb_row$search_method  <- result$search_method
      wb_row$nomenclatural_status <- result$nomenclatural_status
      wb_row$corrected_name <- result$corrected_name
      wb_row$lineage_record_type <- result$lineage_record_type_summary
      wb_row$lineage_fallback    <- result$lineage_fallback_summary
      wb_row$lineage_status      <- result$lineage_status
      
      if (!base::is.null(lineage_data) && base::nrow(lineage_data) > 0) {
        
        for (i in base::seq_len(base::nrow(lineage_data))) {
          
          rank <- lineage_data$category[i]
          
          if (!base::is.na(rank) && rank %in% RANK_ORDER) {
            
            rank_index <- base::match(rank, RANK_ORDER)
            search_rank_index <- base::match(result$search_rank, RANK_ORDER)
            
            if (rank_index <= search_rank_index) {
              wb_row[[rank]] <- lineage_data$full_name[i]
            }
          }
        }
      }
      
      batch_workbook_rows[[taxon]] <- base::as.data.frame(wb_row, stringsAsFactors = FALSE)
    }
    
    
  
    # Write batch results to Excel
  
    
    if (base::length(batch_workbook_rows) > 0) {
      
      final_batch_df <- safe_rbind(batch_workbook_rows)
      sheet_name <- base::paste0("Batch_", batch_number)
      
      wb <-
        if (base::file.exists(workbook_path)) {
          openxlsx::loadWorkbook(workbook_path)
        } else {
          openxlsx::createWorkbook()
        }
      
      if (sheet_name %in% base::names(wb)) {
        openxlsx::removeWorksheet(wb, sheet_name)
      }
      
      openxlsx::addWorksheet(wb, sheetName = sheet_name)
      
      cols <-
        base::c(
          "search_term", "search_rank", "match_type", "search_method",
          "nomenclatural_status", "corrected_name", "lineage_record_type",
          "lineage_fallback", "lineage_status", RANK_ORDER
        )
      
      for (c in base::setdiff(cols, base::names(final_batch_df))) {
        final_batch_df[[c]] <- NA
      }
      
      openxlsx::writeData(wb, sheet = sheet_name, x = final_batch_df[, cols, drop = FALSE])
      openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
    }
    
    
  
    # Write not-found / partial lineage logs
  
    
    not_found_path <- base::file.path(log_path, base::paste0("taxa_not_found_batch_", batch_number, ".txt"))
    
    if (base::length(batch_not_found) > 0) {
      base::writeLines(base::sort(batch_not_found), not_found_path)
    } else {
      base::writeLines("No taxa failed to resolve.", not_found_path)
    }
    
    partial_path <- base::file.path(log_path, base::paste0("partial_lineages_batch_", batch_number, ".txt"))
    
    if (base::length(batch_partial_lineages) > 0) {
      base::writeLines(base::sort(batch_partial_lineages), partial_path)
    } else {
      base::writeLines("No partial/orphan lineages encountered.", partial_path)
    }
  }
  
  

  # 24. Completion

  
  finish_time <- base::Sys.time()
  duration <- base::difftime(finish_time, search_start_time, units = "auto")
  
  base::cat("\n\nLPSN Run Complete at: ", base::format(finish_time, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  base::cat("Total time elapsed: ", base::format(duration), "\n", sep = "")
  base::cat("Total unique taxa processed: ", total_unique_taxa, "\n", sep = "")
  base::cat("Total batches: ", num_batches, "\n", sep = "")
  
  invisible(
    list(
      taxa_names = taxa_names,
      total_unique_taxa = total_unique_taxa,
      total_batches = num_batches,
      workbook_path = workbook_path,
      csv_path = csv_path,
      log_path = log_path
    )
  )
}


# # usage example
# # Output directory structure
# CSV_DIR <- "./XXXXX/csv_lineages"
# WB_PATH <- "./XXXXX/LPSN_result_Master.xlsx"
# LOG_DIR <- "./XXXXX/logs"
# 
# taxa_names <- c("Bacillus cereus", "Clostridium botulinum", "Lactobacillaceae", "Escherichia coli")
# 
# result <- lpsn_levels_minimal_err_no_rank(
#   taxa_names    = taxa_names,
#   csv_path      = CSV_DIR,
#   workbook_path = WB_PATH,
#   log_path      = LOG_DIR
# )





# Validation by direct record number retrieval ################################


# Function: fetch_lpsn_records()
#
# Fetches one or more LPSN records by numeric ID and
# reports which IDs (if any) failed to return a result.


direct_record_id_fetch <- function(record_ids, lpsn_connection = NULL) {
  
  # load library
  library(LPSN)
  
  

  # 1. Validate input

  
  if (
    base::missing(record_ids) ||
    base::length(record_ids) == 0
  ) {
    
    base::stop("`record_ids` must be a non-empty vector of numeric LPSN IDs.")
  }
  
  
  record_ids <- base::as.numeric(record_ids)
  
  
  if (
    base::any(base::is.na(record_ids))
  ) {
    
    base::stop("`record_ids` contains values that could not be converted to numeric.")
  }
  
  
  record_ids <- base::unique(record_ids)
  
  

  # 2. Authenticate 
  # (unless an existing connection was passed in, e.g. to reuse a session across calls)

  
  lpsn <- lpsn_connection
  
  
  if (base::is.null(lpsn)) {
    
    user     <- base::Sys.getenv("DSMZ_API_USER")
    password <- base::Sys.getenv("DSMZ_API_PASSWORD")
    
    
    if (user == "" || password == "") {
      
      base::stop("LPSN credentials not found. Set DSMZ_API_USER and DSMZ_API_PASSWORD in .Renviron.")
    }
    
    
    lpsn <- LPSN::open_lpsn(user, password)
  }
  
  

  # 3. Fetch records

  
  result <- base::tryCatch(
    LPSN::fetch(lpsn, ids = record_ids),
    error = function(e) {
      
      NULL
    }
  )
  
  

  # 4. Handle a completely failed / empty request

  
  if (
    base::is.null(result) ||
    (!base::is.null(result$count) && result$count == 0)
  ) {
    
    return(
      base::list(
        found   = base::data.frame(),
        missing = record_ids
      )
    )
  }
  
  

  # 5. Convert to data frame

  
  record_df <- base::as.data.frame(result)
  
  
  if (base::nrow(record_df) == 0) {
    
    return(
      base::list(
        found   = base::data.frame(),
        missing = record_ids
      )
    )
  }
  
  

  # 6. (Optional) keep only the fields of interest

  
  export_fields <- base::c(
    "id", "category", "full_name", "authority", "proposed_as",
    "lpsn_parent_id", "lpsn_correct_name_id",
    "lpsn_taxonomic_status", "nomenclatural_status"
  )
  
  
  record_df <- record_df[, base::intersect(export_fields, base::names(record_df)), drop = FALSE]
  
  

  # 7. Identify which requested IDs failed to return

  
  returned_ids <- base::as.numeric(record_df$id)
  
  
  missing_ids <- base::setdiff(record_ids, returned_ids)
  
  

  # 8. Return both found records and missing IDs

  
  base::list(
    found   = record_df,
    missing = missing_ids
  )
}

# # Usage example
# 
# result <- fetch_lpsn_records(c(917, 884, 25877, 29408, 9750))
# 
# # Records successfully retrieved
# print(result$found)
# 
# # IDs that failed to return a record
# print(result$missing)



# LPSN Citation ###############################################################################

# When using LPSN database (including the Application Programming Interface or API), cite:

# Freese, H.M., Meier-Kolthoff, J.P., Sardà Carbasse, J., Afolayan, A.O., Göker, M. (2026). TYGS and LPSN in 2025: a Global Core Biodata Resource for genome-based classification and nomenclature of prokaryotes within DSMZ Digital Diversity. Nucleic Acids Research, 54, D884–D891;
# DOI: 10.1093/nar/gkaf1110

# When citing LPSN always mention the date on which you accessed the website directly or 
# through the API.


## citing the LPSN package

# citation(package = "LPSN")

#' To cite package ‘LPSN’ in publications use:
#'   
#'   Goeker M (2022). _LPSN: API Client for LPSN_. R package version 0.14.1/r892,
#' <https://R-Forge.R-project.org/projects/opm/>.
#' 
#' A BibTeX entry for LaTeX users is
#' 
#' @Manual{,
#'   title = {LPSN: API Client for LPSN},
#'   author = {Markus Goeker},
#'   year = {2022},
#'   note = {R package version 0.14.1/r892},
#'   url = {https://R-Forge.R-project.org/projects/opm/},
#' }
#' 
#' ATTENTION: This citation information has been auto-generated from the package DESCRIPTION file and may need manual editing, see ‘help("citation")’.
