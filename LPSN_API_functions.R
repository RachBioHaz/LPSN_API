
#' Utilities and functions used in the LPSN Taxonomic correction walkthrough
#'
#' These functions were created by R Gore in April 2026 with some AI help. 
#' 
#' Use these as needed, however these functions are not being actively maintained,  
#' so I highly recommend you use a subset of your taxonomy file in a test run
#' to validate the function is returning the correct,  Manually
#' checking the results in LPSN. A list of taxa that failed to return a result, 
#' is also provided to enable manual searches.
#' 
#' How long does it take? 
#' The whole process is shorter if you prepare your data, there will be fewer manual searches 
#' needed if the data is prepared properly. 
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

#'
#' Remove Taxonomy Level Prefixes
#'
#' This function cleans a taxonomy data frame by stripping the standard prefixes 
#' (e.g., "d__", "p__") often found in reference databases like SILVA or Greengenes.
#'
#' @param tax A data frame or tibble where columns 1 through 7 represent 
#'   taxonomic ranks (Domain to Species) containing prefixes.
#'
#' @return A modified data frame (tibble) with the prefixes removed from the 
#'   first seven taxonomic columns.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'df' has entries like "p__Firmicutes"
#' df_clean <- remove_taxa_prefix(df)
#' }
#' 

remove_taxa_prefix <- function(tax) {
  # Ensure the stringr library is available for the replacement function
  # In a package context, use stringr::str_replace directly
  
  tax |> 
    # 1. Convert to data frame so dplyr::mutate can process it
    base::as.data.frame() |>
    # 2. Strip prefixes using column names
    dplyr::mutate(
      # We refer to the columns directly. 
      # If your columns aren't named exactly like this, 
      # we can use across() for a more dynamic approach.
      Domain  = stringr::str_replace(Domain,  "^d__", ""),
      Phylum  = stringr::str_replace(Phylum,  "^p__", ""),
      Class   = stringr::str_replace(Class,   "^c__", ""),
      Order   = stringr::str_replace(Order,   "^o__", ""),
      Family  = stringr::str_replace(Family,  "^f__", ""),
      Genus   = stringr::str_replace(Genus,   "^g__", ""),
      Species = stringr::str_replace(Species, "^s__", "")
    )
}






#' Update Taxonomic Names with Logging
#'
#' Replaces old taxonomic names with new ones based on a provided mapping list.
#' It processes the data frame column by column and prints a summary of how 
#' many changes were made in each rank.
#' 
#' This can be used to bulk change taxa names as part of data preparation for
#' the LPSN API function
#'
#' @param df A data frame containing taxonomic data.
#' @param mappings A named list of named vectors. The list names must match 
#'   column names in \code{df}. The vectors should be in the format 
#'   \code{c("old_name" = "new_name")}.
#'
#'# Define the mapping list
#'# list each column and the values to change in that column
#'# Format: column_name =C("Old Name" = "New Name", ....) 
#'# Only changes the names in the lists
#'
#'my_mappings <- list(
#'  Phylum = c(
#'  "Bacteroidetes" = "Bacteroidota", 
#'  "Firmicutes"    = "Bacillota"),
#'  Genus = c(
#'  "Lactobacillus" = "Lacticaseibacillus",
#'  "Clostridium"   = "Clostridioides")
#'  )
#'
#'
#' @return A data frame with updated taxonomic names and a printed log 
#'   of changes per column.
#' @export
#'
clean_taxonomy_with_log <- function(df, mappings) {
  out <- df
  
  for (col_name in base::names(mappings)) {
    if (col_name %in% base::colnames(out)) {
      
      # 1. Capture the current column values to compare later
      old_vals <- out[[col_name]]
      
      # 2. Perform the recoding
      # !!col_name unquotes the variable name, !!!mappings splices the list of changes
      out <- out |> 
        dplyr::mutate(
          !!col_name := dplyr::recode(.data[[col_name]], !!!mappings[[col_name]])
        )
      
      # 3. Calculate and print the number of changes
      new_vals <- out[[col_name]]
      n_changes <- base::sum(old_vals != new_vals, na.rm = TRUE)
      
      base::message(base::paste0("Column [", col_name, "]: ", n_changes, " names updated."))
    }
  }
  return(out)
}




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
  
  # 1. Print the length comparison in X:Y format
  base::cat(base::length(vec1), ":", base::length(vec2), sep = "", "\n\n")
  
  # 2. Equalise lengths to the maximum found to prevent recycling errors
  max_len <- base::max(base::length(vec1), base::length(vec2))
  
  v1_aligned <- vec1
  base::length(v1_aligned) <- max_len
  
  v2_aligned <- vec2
  base::length(v2_aligned) <- max_len
  
  # 3. Create the data frame for visual inspection
  df <- base::data.frame(
    Original = v1_aligned,
    Updated = v2_aligned,
    stringsAsFactors = FALSE
  )
  
  return(df)
}

#### Example usage
# phylum_comparison <- compare_vectors(phylum_list, phylum_list_uptd)
# 
## To view the new list
# print(phylum_comparison)





################ Using the LSPN API
#
#
### LPSN citation
# When using LPSN (including the Application Programming Interface or API), cite:
#   
# Freese, H.M., Meier-Kolthoff, J.P., Sardà Carbasse, J., Afolayan, A.O., Göker, M. (2026). TYGS and LPSN in 2025: a Global Core Biodata Resource for genome-based classification and nomenclature of prokaryotes within DSMZ Digital Diversity. Nucleic Acids Research, 54, D884–D891; 
# DOI: 10.1093/nar/gkaf1110
# 
# When citing LPSN always mention the date on which you accessed the website or the API.


#### API calls to LPSN - works in batches
## does not check for candidatus - use when you want confirmed names only

process_lpsn_batches <- function(taxa_list, csv_path, workbook_path, log_path, user, password, 
                                        batch_size = 20, pause_min = 1, pause_max = 3, use_timestamp = TRUE) {
  
  start_time <- base::Sys.time()
  
  # --- 1. Deduplication & Logging ---
  if (!base::dir.exists(log_path)) base::dir.create(log_path, recursive = TRUE)
  
  original_count <- base::length(taxa_list)
  unique_taxa <- base::unique(base::as.character(taxa_list))
  unique_taxa <- unique_taxa[!base::is.na(unique_taxa) & unique_taxa != ""]
  duplicates <- base::setdiff(taxa_list, unique_taxa)
  
  log_content <- if (base::length(duplicates) > 0) { 
    base::sort(base::unique(duplicates)) 
  } else { 
    "No duplicates found in the provided taxa list." 
  }
  base::writeLines(log_content, base::file.path(log_path, "deduplicated_taxa.txt"))
  
  taxa_list <- unique_taxa
  new_count <- base::length(taxa_list)
  base::cat(base::sprintf("Input: %d taxa. Unique: %d. Deduplication log saved to 'deduplicated_taxa.txt'.\n", original_count, new_count))
  
  # --- 2. Setup Environments & Cache ---
  api_cache <- base::list()
  
  if (use_timestamp) {
    ts <- base::format(base::Sys.time(), "%Y%m%d_%H%M")
    workbook_path <- base::gsub("\\.xlsx$", base::paste0("_", ts, ".xlsx"), workbook_path)
  }
  
  if (!base::dir.exists(csv_path)) base::dir.create(csv_path, recursive = TRUE)
  
  RANK_ORDER <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  HIGHER_RANKS <- c("phylum", "class", "order", "family", "genus", "kingdom", "domain")
  export_fields <- c("id", "category", "full_name", "lpsn_parent_id", "lpsn_correct_name_id", "lpsn_taxonomic_status")
  
  # Helper: Merges dataframes with different columns
  safe_rbind <- function(df_list) {
    if (base::length(df_list) == 0) return(NULL)
    all_cols <- base::unique(base::unlist(base::lapply(df_list, base::names)))
    res_list <- base::lapply(df_list, function(df) {
      for (col in base::setdiff(all_cols, base::names(df))) df[[col]] <- NA
      df[all_cols]
    })
    base::do.call(base::rbind, res_list)
  }
  
  # New Helper: Retry logic for API calls
  retry_api_call <- function(expr, retries = 3, wait = 5) {
    for (i in 1:retries) {
      res <- base::tryCatch(expr, error = function(e) NULL)
      if (!base::is.null(res)) return(res)
      if (i < retries) {
        base::message(base::sprintf("\n  API timeout or error. Retrying (%d/%d) in %ds...", i, retries, wait))
        base::Sys.sleep(wait)
      }
    }
    return(NULL)
  }
  
  find_direct_hit <- function(query_res, target_name, target_category) {
    df_res <- base::as.data.frame(query_res)
    if (base::nrow(df_res) == 0) return(NULL)
    hit_idx <- base::which(base::tolower(df_res$full_name) == base::tolower(target_name) & df_res$category == target_category)
    if (base::length(hit_idx) > 0) return(df_res[hit_idx, , drop = FALSE])
    return(NULL)
  }
  
  lpsn <- LPSN::open_lpsn(user, password)
  num_batches <- base::ceiling(new_count / batch_size)
  pb <- progress::progress_bar$new(format = " [:bar] :percent Batch :batch/:total_batches ETA: :eta", total = new_count, width = 80)
  
  # --- 3. Main Loop ---
  for (b in 1:num_batches) {
    start_idx <- ((b - 1) * batch_size) + 1
    end_idx <- base::min(b * batch_size, new_count)
    current_batch_taxa <- taxa_list[start_idx:end_idx]
    
    batch_workbook_rows <- list()
    batch_not_found <- character()
    
    for (taxon in current_batch_taxa) {
      pb$tick(tokens = base::list(batch = b, total_batches = num_batches))
      
      match_type <- "Full Match"
      res <- NULL
      is_binomial <- base::grepl(" ", taxon)
      
      if (is_binomial) {
        base::Sys.sleep(stats::runif(1, pause_min, pause_max))
        raw_res <- retry_api_call(LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = taxon, category = "species"))
        res <- find_direct_hit(raw_res, taxon, "species")
        
        if (base::is.null(res)) {
          genus_only <- if (base::grepl("^Candidatus", taxon, ignore.case = TRUE)) base::sub("^([^ ]+ +[^ ]+).*", "\\1", taxon) else base::sub(" .*", "", taxon)
          match_type <- "Genus Fallback"
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          raw_res <- retry_api_call(LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = genus_only, category = "genus"))
          res <- find_direct_hit(raw_res, genus_only, "genus")
        }
      } else {
        for (rank_to_check in HIGHER_RANKS) {
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          raw_res <- retry_api_call(LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = taxon, category = rank_to_check))
          res <- find_direct_hit(raw_res, taxon, rank_to_check)
          if (!base::is.null(res)) {
            match_type <- base::paste0("Full Match (", tools::toTitleCase(rank_to_check), ")")
            break
          }
        }
      }
      
      if (base::is.null(res)) {
        batch_not_found <- base::c(batch_not_found, base::as.character(taxon))
        next
      }
      
      first_hit <- base::as.data.frame(res)[1, , drop = FALSE]
      is_synonym <- !base::is.null(first_hit$lpsn_correct_name_id) && 
        !base::is.na(first_hit$lpsn_correct_name_id) && 
        first_hit$lpsn_correct_name_id != 0 && 
        first_hit$lpsn_correct_name_id != first_hit$id
      
      current_id <- if (is_synonym) first_hit$lpsn_correct_name_id else first_hit$id
      lineage_list <- list()
      corrected_name <- NA 
      
      while (!base::is.null(current_id) && !base::is.na(current_id) && current_id != 0) {
        id_key <- base::as.character(current_id)
        if (id_key %in% base::names(api_cache)) {
          rec_df <- api_cache[[id_key]]
        } else {
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          rec_res <- retry_api_call(LPSN::fetch(lpsn, ids = current_id))
          if (base::is.null(rec_res)) break
          rec_df <- base::as.data.frame(rec_res)[1, , drop = FALSE]
          api_cache[[id_key]] <- rec_df
        }
        
        if (rec_df$category != "division") {
          lineage_list[[base::length(lineage_list) + 1]] <- rec_df[, base::intersect(export_fields, base::names(rec_df)), drop = FALSE]
        }
        if (base::length(lineage_list) == 1) corrected_name <- rec_df$full_name
        current_id <- if (base::length(rec_df$lpsn_parent_id) > 0 && rec_df$lpsn_parent_id != 0) rec_df$lpsn_parent_id else NULL
      }
      
      if (base::length(lineage_list) > 0) {
        lineage_data <- safe_rbind(lineage_list)
        safe_name <- base::gsub("[^[:alnum:]]", "_", taxon)
        utils::write.csv(lineage_data, base::file.path(csv_path, base::paste0("Batch_", b, "_", safe_name, ".csv")), row.names = FALSE)
        
        wb_row <- stats::setNames(base::as.list(base::rep(NA, base::length(RANK_ORDER))), RANK_ORDER)
        wb_row$search_term <- taxon
        wb_row$match_type <- if (is_synonym) base::paste0(match_type, " (Synonym Redirect)") else match_type
        wb_row$nomenclatural_status <- if (!base::is.null(first_hit$lpsn_taxonomic_status)) first_hit$lpsn_taxonomic_status else "Unknown"
        wb_row$corrected_name <- corrected_name 
        
        for (i in 1:base::nrow(lineage_data)) {
          rank <- lineage_data$category[i]
          if (!base::is.null(rank) && !base::is.na(rank) && rank %in% RANK_ORDER) wb_row[[rank]] <- lineage_data$full_name[i]
        }
        batch_workbook_rows[[taxon]] <- base::as.data.frame(wb_row)
      }
    }
    
    if (base::length(batch_workbook_rows) > 0) {
      final_batch_df <- safe_rbind(batch_workbook_rows)
      sheet_name <- base::paste0("Batch_", b)
      wb <- if (base::file.exists(workbook_path)) openxlsx::loadWorkbook(workbook_path) else openxlsx::createWorkbook()
      if (sheet_name %in% base::names(wb)) openxlsx::removeWorksheet(wb, sheet_name)
      openxlsx::addWorksheet(wb, sheetName = sheet_name)
      cols <- base::c("search_term", "match_type", "nomenclatural_status", "corrected_name", RANK_ORDER)
      for (c in base::setdiff(cols, base::names(final_batch_df))) final_batch_df[[c]] <- NA
      openxlsx::writeData(wb, sheet = sheet_name, x = final_batch_df[, cols])
      openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
    }
    base::writeLines(base::sort(base::as.character(base::unlist(batch_not_found))), base::file.path(log_path, base::paste0("taxa_not_found_batch_", b, ".txt")))
  }
  
  finish_time <- base::Sys.time()
  duration <- base::difftime(finish_time, start_time, units = "auto")
  base::cat(base::sprintf("\nLPSN Run Complete at: %s\n", base::format(finish_time, "%Y-%m-%d %H:%M:%S")))
  base::cat(base::sprintf("Total time elapsed: %s\n", base::format(duration)))
}

# # Your cleaned taxa list 
# my_taxa <- c("Firmicutes",..., "Clostridium")

# # Output directory structure
# CSV_DIR <- "./LPSN/results/csv_lineages"
# WB_PATH <- "./LPSN/results/LPSN_result_Master.xlsx"
# LOG_DIR <- "./LPSN/results/logs"
# 
# # run API function
# # set batch size as needed. Small batches is a bit slower as it opens and saves to excel more often
# 
# process_lpsn_batches(
#   taxa_list     = my_taxa_clean,
#   csv_path      = CSV_DIR,
#   workbook_path = WB_PATH,
#   log_path      = LOG_DIR,
#   user          = Sys.getenv("DSMZ_API_USER"),
#   password      = Sys.getenv("DSMZ_API_PASSWORD"),
#   batch_size    = 50,    # Processes 50 unique taxa per Excel sheet
#   pause_min     = 0.5,   # Minimum delay between API calls (seconds)
#   pause_max     = 1.5,   # Maximum delay between API calls (seconds)
#   use_timestamp = TRUE   # Appends date/time to the Excel filename
# )



##############
############ LPSN API call - searches for candidatus no quotes



#' Process LPSN Taxonomy Batches
#'
#' Connects to the LPSN API to retrieve validated taxonomic information for a list of taxa. 
#' It handles deduplication, hierarchical fallback searching, synonym redirection, 
#' and exports results into both per-taxon CSVs and a master Excel workbook.
#' 
#' --- Prerequisites ---
#' Ensure you have your LPSN credentials set in your environment (.Renviron)
#' DSMZ_API_USER="your_email@example.com"
#' DSMZ_API_PASSWORD="your_password"
#'
#' @param taxa_list Character vector. List of taxonomic names to search.
#' @param csv_path String. Directory path to save individual lineage CSV files.
#' @param workbook_path String. Path to the master .xlsx file.
#' @param log_path String. Directory path to save deduplication and "not found" logs.
#' @param user String. LPSN API username.
#' @param password String. LPSN API password.
#' @param batch_size Integer. Number of taxa to process per Excel worksheet (default: 20).
#' @param pause_min Numeric. Minimum seconds to wait between API calls (default: 1).
#' @param pause_max Numeric. Maximum seconds to wait between API calls (default: 3).
#' @param use_timestamp Logical. If TRUE, appends a timestamp to the workbook filename.
#'
#' @return No return value. Outputs files to the specified paths and prints progress to console.
#' @export


##### this is old 1
# process_lpsn_cached_batches <- function(taxa_list, csv_path, workbook_path, 
#                                              log_path, user, password, batch_size = 20, 
#                                              pause_min = 1, pause_max = 3, use_timestamp = TRUE) {
#   start_time <- base::Sys.time()
#   
#   # --- 1. Deduplication & Logging ---
#   # Ensures we don't pay for the same API call twice.
#   if (!base::dir.exists(log_path)) base::dir.create(log_path, recursive = TRUE)
#   original_count <- base::length(taxa_list)
#   
#   unique_taxa <- base::as.character(taxa_list)
#   unique_taxa <- base::trimws(unique_taxa) # Remove leading/trailing spaces
#   unique_taxa <- base::gsub("\\s+", " ", unique_taxa) # Collapse internal double spaces
#   unique_taxa <- base::unique(unique_taxa)
#   unique_taxa <- unique_taxa[!base::is.na(unique_taxa) & unique_taxa != ""]
#   
#   duplicates <- base::setdiff(taxa_list, unique_taxa)
#   log_content <- if (base::length(duplicates) > 0) base::sort(base::unique(duplicates)) else "No duplicates found."
#   base::writeLines(log_content, base::file.path(log_path, "deduplicated_taxa.txt"))
#   
#   taxa_list <- unique_taxa
#   new_count <- base::length(taxa_list)
#   base::cat(base::sprintf("Input: %d taxa. Unique/Cleaned: %d.\n", original_count, new_count))
#   
#   # --- 2. Setup Environments & Helpers ---
#   api_cache <- base::list()
#   if (use_timestamp) {
#     ts <- base::format(base::Sys.time(), "%Y%m%d_%H%M")
#     workbook_path <- base::gsub("\\.xlsx$", base::paste0("_", ts, ".xlsx"), workbook_path)
#   }
#   if (!base::dir.exists(csv_path)) base::dir.create(csv_path, recursive = TRUE)
#   
#   RANK_ORDER <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
#   HIGHER_RANKS <- c("phylum", "class", "order", "family", "genus", "kingdom", "domain")
#   export_fields <- c("id", "category", "full_name", 
#                      "lpsn_parent_id", "lpsn_correct_name_id", 
#                      "lpsn_taxonomic_status")
#   
#   safe_rbind <- function(df_list) {
#     if (base::length(df_list) == 0) return(NULL)
#     all_cols <- base::unique(base::unlist(base::lapply(df_list, base::names)))
#     res_list <- base::lapply(df_list, function(df) {
#       for (col in base::setdiff(all_cols, base::names(df))) df[[col]] <- NA
#       df[all_cols]
#     })
#     base::do.call(base::rbind, res_list)
#   }
#   
#   retry_api_call <- function(expr, retries = 3, wait = 5) {
#     for (i in 1:retries) {
#       res <- base::tryCatch(expr, error = function(e) NULL)
#       if (!base::is.null(res)) return(res)
#       if (i < retries) {
#         base::message(base::sprintf("\n API error. Retry %d/%d...", i, retries))
#         base::Sys.sleep(wait)
#       }
#     }
#     return(NULL)
#   }
#   
#   find_direct_hit <- function(query_res, target_name, target_category) {
#     df_res <- base::as.data.frame(query_res)
#     if (base::nrow(df_res) == 0) return(NULL)
#     hit_idx <- base::which(base::tolower(df_res$full_name) == base::tolower(target_name) & df_res$category == target_category)
#     if (base::length(hit_idx) > 0) return(df_res[hit_idx, , drop = FALSE])
#     return(NULL)
#   }
#   
#   # Candidatus retry wrapper for Species/Genus level only
#   search_with_candidatus_retry <- function(lpsn_obj, name, rank) {
#     base::Sys.sleep(stats::runif(1, pause_min, pause_max))
#     raw <- retry_api_call(LPSN::retrieve(lpsn_obj, search = "advanced", 
#                                          `taxon-name` = name, 
#                                          category = rank))
#     hit <- find_direct_hit(raw, name, rank)
#     
#     if (base::is.null(hit) && !base::grepl("^Candidatus", name, ignore.case = TRUE)) {
#       cand_name <- base::paste("Candidatus", name)
#       base::message(base::sprintf("\n '%s' not found. Retrying as '%s'...", name, cand_name))
#       base::Sys.sleep(stats::runif(1, pause_min, pause_max))
#       raw_cand <- retry_api_call(LPSN::retrieve(lpsn_obj, 
#                                                 search = "advanced", 
#                                                 `taxon-name` = cand_name, 
#                                                 category = rank))
#       hit <- find_direct_hit(raw_cand, cand_name, rank)
#       if (!base::is.null(hit)) match_type <<- base::paste0(match_type, " (Candidatus)")
#     }
#     return(hit)
#   }
#   
#   lpsn <- LPSN::open_lpsn(user, password)
#   num_batches <- base::ceiling(new_count / batch_size)
#   pb <- progress::progress_bar$new(format = " [:bar] :percent Batch :batch/:total_batches ETA: :eta", total = new_count, width = 80)
#   
#   # --- 3. Main Loop ---
#   for (b in 1:num_batches) {
#     start_idx <- ((b - 1) * batch_size) + 1
#     end_idx <- base::min(b * batch_size, new_count)
#     current_batch_taxa <- taxa_list[start_idx:end_idx]
#     batch_workbook_rows <- list()
#     batch_not_found <- character()
#     
#     for (taxon in current_batch_taxa) {
#       pb$tick(tokens = base::list(batch = b, total_batches = num_batches))
#       match_type <- "Full Match"
#       res <- NULL
#       is_binomial <- base::grepl(" ", taxon)
#       
#       if (is_binomial) {
#         res <- search_with_candidatus_retry(lpsn, taxon, "species")
#         if (base::is.null(res)) {
#           genus_only <- if (base::grepl("^Candidatus", taxon, ignore.case = TRUE)) {
#             base::sub("^([^ ]+ +[^ ]+).*", "\\1", taxon)
#           } else {
#             base::sub(" .*", "", taxon)
#           }
#           match_type <- "Genus Fallback"
#           res <- search_with_candidatus_retry(lpsn, genus_only, "genus")
#         }
#       } else {
#         # Loop through Higher Ranks: Standard search ONLY (no Candidatus retry)
#         for (rank_to_check in HIGHER_RANKS) {
#           base::Sys.sleep(stats::runif(1, pause_min, pause_max))
#           raw_res <- retry_api_call(LPSN::retrieve(lpsn, 
#                                                    search = "advanced", 
#                                                    `taxon-name` = taxon, 
#                                                    category = rank_to_check))
#           res <- find_direct_hit(raw_res, taxon, rank_to_check)
#           if (!base::is.null(res)) {
#             match_type <- base::paste0("Full Match (", tools::toTitleCase(rank_to_check), ")")
#             break
#           }
#         }
#       }
#       
#       if (base::is.null(res)) {
#         batch_not_found <- base::c(batch_not_found, base::as.character(taxon))
#         next
#       }
#       
#       # --- 4. Lineage Retrieval & Synonym Resolution ---
#       first_hit <- base::as.data.frame(res)[1, , drop = FALSE]
#       is_synonym <- !base::is.null(first_hit$lpsn_correct_name_id) && 
#         !base::is.na(first_hit$lpsn_correct_name_id) && 
#         first_hit$lpsn_correct_name_id != 0 && 
#         first_hit$lpsn_correct_name_id != first_hit$id
#       
#       current_id <- if (is_synonym) first_hit$lpsn_correct_name_id else first_hit$id
#       lineage_list <- list()
#       corrected_name <- NA
#       
#       while (!base::is.null(current_id) && !base::is.na(current_id) && current_id != 0) {
#         id_key <- base::as.character(current_id)
#         if (id_key %in% base::names(api_cache)) {
#           rec_df <- api_cache[[id_key]]
#         } else {
#           base::Sys.sleep(stats::runif(1, pause_min, pause_max))
#           rec_res <- retry_api_call(LPSN::fetch(lpsn, ids = current_id))
#           if (base::is.null(rec_res)) break
#           rec_df <- base::as.data.frame(rec_res)[1, , drop = FALSE]
#           api_cache[[id_key]] <- rec_df
#         }
#         
#         # Safety check for "length zero" error
#         if (base::length(rec_df$category) > 0 && !base::is.na(rec_df$category) && rec_df$category != "division") {
#           lineage_list[[base::length(lineage_list) + 1]] <- rec_df[, base::intersect(export_fields, base::names(rec_df)), drop = FALSE]
#         }
#         
#         if (base::length(lineage_list) == 1) corrected_name <- rec_df$full_name
#         current_id <- if (base::length(rec_df$lpsn_parent_id) > 0 && rec_df$lpsn_parent_id != 0) rec_df$lpsn_parent_id else NULL
#       }
#       
#       # --- 5. Export Logic (Individual CSVs) ---
#       if (base::length(lineage_list) > 0) {
#         lineage_data <- safe_rbind(lineage_list)
#         safe_name <- base::gsub("[^[:alnum:]]", "_", taxon)
#         utils::write.csv(lineage_data, base::file.path(csv_path, base::paste0("Batch_", b, "_", safe_name, ".csv")), row.names = FALSE)
#         
#         wb_row <- stats::setNames(base::as.list(base::rep(NA, base::length(RANK_ORDER))), RANK_ORDER)
#         wb_row$search_term <- taxon
#         wb_row$match_type <- if (is_synonym) base::paste0(match_type, " (Synonym Redirect)") else match_type
#         wb_row$nomenclatural_status <- if (base::length(first_hit$lpsn_taxonomic_status) > 0) first_hit$lpsn_taxonomic_status else "Unknown"
#         wb_row$corrected_name <- corrected_name
#         
#         for (i in 1:base::nrow(lineage_data)) {
#           rank <- lineage_data$category[i]
#           if (!base::is.null(rank) && !base::is.na(rank) && rank %in% RANK_ORDER) wb_row[[rank]] <- lineage_data$full_name[i]
#         }
#         batch_workbook_rows[[taxon]] <- base::as.data.frame(wb_row)
#       }
#     }
#     
#     # --- 6. Save Excel Batch with tryCatch ---
#     if (base::length(batch_workbook_rows) > 0) {
#       final_batch_df <- safe_rbind(batch_workbook_rows)
#       sheet_name <- base::paste0("Batch_", b)
#       
#       base::tryCatch({
#         wb <- if (base::file.exists(workbook_path)) openxlsx::loadWorkbook(workbook_path) else openxlsx::createWorkbook()
#         if (sheet_name %in% base::names(wb)) openxlsx::removeWorksheet(wb, sheet_name)
#         openxlsx::addWorksheet(wb, sheetName = sheet_name)
#         
#         cols <- base::c("search_term", "match_type", "nomenclatural_status", "corrected_name", RANK_ORDER)
#         for (col_name in base::setdiff(cols, base::names(final_batch_df))) final_batch_df[[col_name]] <- NA
#         
#         openxlsx::writeData(wb, sheet = sheet_name, x = final_batch_df[, cols])
#         openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
#       }, error = function(e) {
#         base::message(base::sprintf("\n[ERROR] Could not save Batch %d to Excel. Is '%s' open?", b, base::basename(workbook_path)))
#       })
#     }
#     
#     base::writeLines(base::sort(base::as.character(base::unlist(batch_not_found))), 
#                      base::file.path(log_path, base::paste0("taxa_not_found_batch_", b, ".txt")))
#   }
#   
#   duration <- base::difftime(base::Sys.time(), start_time, units = "auto")
#   base::cat(base::sprintf("\nRun Complete. Total time: %s\n", base::format(duration)))
# }


##### check this works- replaced old 1
#' Process LPSN Taxonomy Batches
#'
#' Connects to the LPSN API to retrieve validated taxonomic information for a list of taxa. 
#' It handles deduplication, hierarchical fallback searching, synonym redirection, 
#' and exports results into both per-taxon CSVs and a master Excel workbook.
#'
#' @param taxa_list Character vector. List of taxonomic names to search.
#' @param csv_path String. Directory path to save individual lineage CSV files.
#' @param workbook_path String. Path to the master .xlsx file.
#' @param log_path String. Directory path to save deduplication and "not found" logs.
#' @param user String. LPSN API username.
#' @param password String. LPSN API password.
#' @param batch_size Integer. Number of taxa to process per Excel worksheet (default: 20).
#' @param pause_min Numeric. Minimum seconds to wait between API calls (default: 1).
#' @param pause_max Numeric. Maximum seconds to wait between API calls (default: 3).
#' @param use_timestamp Logical. If TRUE, appends a timestamp to the workbook filename.
#'
#' @return No return value. Outputs files to the specified paths and prints progress to console.
#' @export
process_lpsn_batches <- function(taxa_list, csv_path, workbook_path, log_path, user, password, 
                                 batch_size = 20, pause_min = 1, pause_max = 3, use_timestamp = TRUE) {
  
  start_time <- base::Sys.time()
  
  # --- 1. Deduplication & Logging ---
  # Ensures we don't pay for the same API call twice.
  if (!base::dir.exists(log_path)) base::dir.create(log_path, recursive = TRUE)
  
  original_count <- base::length(taxa_list)
  unique_taxa <- base::unique(base::as.character(taxa_list))
  unique_taxa <- unique_taxa[!base::is.na(unique_taxa) & unique_taxa != ""]
  
  duplicates <- base::setdiff(taxa_list, unique_taxa)
  log_content <- if (base::length(duplicates) > 0) {
    base::sort(base::unique(duplicates))
  } else {
    "No duplicates found."
  }
  base::writeLines(log_content, base::file.path(log_path, "deduplicated_taxa.txt"))
  
  taxa_list <- unique_taxa
  new_count <- base::length(taxa_list)
  base::cat(base::sprintf("Input: %d taxa. Unique: %d.\n", original_count, new_count))
  
  # --- 2. Setup Environments & Internal Helpers ---
  api_cache <- base::list() # Minimizes redundant API fetches for parent ranks
  
  if (use_timestamp) {
    ts <- base::format(base::Sys.time(), "%Y%m%d_%H%M")
    workbook_path <- base::gsub("\\.xlsx$", base::paste0("_", ts, ".xlsx"), workbook_path)
  }
  
  if (!base::dir.exists(csv_path)) base::dir.create(csv_path, recursive = TRUE)
  
  RANK_ORDER <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  HIGHER_RANKS <- c("phylum", "class", "order", "family", "genus", "kingdom", "domain")
  export_fields <- c("id", "category", "full_name", "lpsn_parent_id", "lpsn_correct_name_id", "lpsn_taxonomic_status")
  
  # Helper: Handles column mismatches when combining dataframes
  safe_rbind <- function(df_list) {
    if (base::length(df_list) == 0) return(NULL)
    all_cols <- base::unique(base::unlist(base::lapply(df_list, base::names)))
    res_list <- base::lapply(df_list, function(df) {
      for (col in base::setdiff(all_cols, base::names(df))) df[[col]] <- NA
      df[all_cols]
    })
    base::do.call(base::rbind, res_list)
  }
  
  # Helper: Manages API timeouts with retries
  retry_api_call <- function(expr, retries = 3, wait = 5) {
    for (i in 1:retries) {
      res <- base::tryCatch(expr, error = function(e) NULL)
      if (!base::is.null(res)) return(res)
      if (i < retries) {
        base::message(base::sprintf("\n API timeout. Retrying (%d/%d)...", i, retries))
        base::Sys.sleep(wait)
      }
    }
    return(NULL)
  }
  
  # Helper: Verifies that the API result actually matches the query string and rank
  find_direct_hit <- function(query_res, target_name, target_category) {
    df_res <- base::as.data.frame(query_res)
    if (base::nrow(df_res) == 0) return(NULL)
    hit_idx <- base::which(base::tolower(df_res$full_name) == base::tolower(target_name) & df_res$category == target_category)
    if (base::length(hit_idx) > 0) return(df_res[hit_idx, , drop = FALSE])
    return(NULL)
  }
  
  # Authenticate with LPSN
  lpsn <- LPSN::open_lpsn(user, password)
  num_batches <- base::ceiling(new_count / batch_size)
  pb <- progress::progress_bar$new(format = " [:bar] :percent Batch :batch/:total_batches ETA: :eta", total = new_count, width = 80)
  
  # --- 3. Main Loop: Processing Taxa ---
  for (b in 1:num_batches) {
    start_idx <- ((b - 1) * batch_size) + 1
    end_idx <- base::min(b * batch_size, new_count)
    current_batch_taxa <- taxa_list[start_idx:end_idx]
    
    batch_workbook_rows <- list()
    batch_not_found <- character()
    
    for (taxon in current_batch_taxa) {
      pb$tick(tokens = base::list(batch = b, total_batches = num_batches))
      
      match_type <- "Full Match"
      res <- NULL
      is_binomial <- base::grepl(" ", taxon)
      
      # A. Primary Search (Species or Higher Rank)
      if (is_binomial) {
        base::Sys.sleep(stats::runif(1, pause_min, pause_max))
        raw_res <- retry_api_call(LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = taxon, category = "species"))
        res <- find_direct_hit(raw_res, taxon, "species")
        
        # B. Fallback to Genus if Species is not found
        if (base::is.null(res)) {
          genus_only <- if (base::grepl("^Candidatus", taxon, ignore.case = TRUE)) base::sub("^([^ ]+ +[^ ]+).*", "\\1", taxon) else base::sub(" .*", "", taxon)
          match_type <- "Genus Fallback"
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          raw_res <- retry_api_call(LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = genus_only, category = "genus"))
          res <- find_direct_hit(raw_res, genus_only, "genus")
        }
      } else {
        # Search higher ranks if input is a single word
        for (rank_to_check in HIGHER_RANKS) {
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          raw_res <- retry_api_call(LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = taxon, category = rank_to_check))
          res <- find_direct_hit(raw_res, taxon, rank_to_check)
          if (!base::is.null(res)) {
            match_type <- base::paste0("Full Match (", tools::toTitleCase(rank_to_check), ")")
            break
          }
        }
      }
      
      if (base::is.null(res)) {
        batch_not_found <- base::c(batch_not_found, base::as.character(taxon))
        next
      }
      
      # C. Handle Synonyms and Correct Names
      first_hit <- base::as.data.frame(res)[1, , drop = FALSE]
      is_synonym <- !base::is.null(first_hit$lpsn_correct_name_id) && !base::is.na(first_hit$lpsn_correct_name_id) && first_hit$lpsn_correct_name_id != 0 && first_hit$lpsn_correct_name_id != first_hit$id
      current_id <- if (is_synonym) first_hit$lpsn_correct_name_id else first_hit$id
      
      # D. Recursive Lineage Fetching
      # Climbs the taxonomic tree using lpsn_parent_id
      lineage_list <- list()
      corrected_name <- NA
      while (!base::is.null(current_id) && !base::is.na(current_id) && current_id != 0) {
        id_key <- base::as.character(current_id)
        if (id_key %in% base::names(api_cache)) {
          rec_df <- api_cache[[id_key]]
        } else {
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          rec_res <- retry_api_call(LPSN::fetch(lpsn, ids = current_id))
          if (base::is.null(rec_res)) break
          rec_df <- base::as.data.frame(rec_res)[1, , drop = FALSE]
          api_cache[[id_key]] <- rec_df
        }
        
        if (rec_df$category != "division") {
          lineage_list[[base::length(lineage_list) + 1]] <- rec_df[, base::intersect(export_fields, base::names(rec_df)), drop = FALSE]
        }
        if (base::length(lineage_list) == 1) corrected_name <- rec_df$full_name
        current_id <- if (base::length(rec_df$lpsn_parent_id) > 0 && rec_df$lpsn_parent_id != 0) rec_df$lpsn_parent_id else NULL
      }
      
      # E. Data Formatting and CSV Export
      if (base::length(lineage_list) > 0) {
        lineage_data <- safe_rbind(lineage_list)
        safe_name <- base::gsub("[^[:alnum:]]", "_", taxon)
        utils::write.csv(lineage_data, base::file.path(csv_path, base::paste0("Batch_", b, "_", safe_name, ".csv")), row.names = FALSE)
        
        wb_row <- stats::setNames(base::as.list(base::rep(NA, base::length(RANK_ORDER))), RANK_ORDER)
        wb_row$search_term <- taxon
        wb_row$match_type <- if (is_synonym) base::paste0(match_type, " (Synonym Redirect)") else match_type
        wb_row$nomenclatural_status <- if (!base::is.null(first_hit$lpsn_taxonomic_status)) first_hit$lpsn_taxonomic_status else "Unknown"
        wb_row$corrected_name <- corrected_name
        
        for (i in 1:base::nrow(lineage_data)) {
          rank <- lineage_data$category[i]
          if (!base::is.null(rank) && !base::is.na(rank) && rank %in% RANK_ORDER) wb_row[[rank]] <- lineage_data$full_name[i]
        }
        batch_workbook_rows[[taxon]] <- base::as.data.frame(wb_row)
      }
    }
    
    # --- 4. Finalize Batch Excel Save ---
    if (base::length(batch_workbook_rows) > 0) {
      final_batch_df <- safe_rbind(batch_workbook_rows)
      sheet_name <- base::paste0("Batch_", b)
      wb <- if (base::file.exists(workbook_path)) openxlsx::loadWorkbook(workbook_path) else openxlsx::createWorkbook()
      
      if (sheet_name %in% base::names(wb)) openxlsx::removeWorksheet(wb, sheet_name)
      openxlsx::addWorksheet(wb, sheetName = sheet_name)
      
      cols <- base::c("search_term", "match_type", "nomenclatural_status", "corrected_name", RANK_ORDER)
      for (c in base::setdiff(cols, base::names(final_batch_df))) final_batch_df[[c]] <- NA
      openxlsx::writeData(wb, sheet = sheet_name, x = final_batch_df[, cols])
      openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
    }
    
    # Save the "Not Found" log for this batch
    base::writeLines(base::sort(base::as.character(base::unlist(batch_not_found))), base::file.path(log_path, base::paste0("taxa_not_found_batch_", b, ".txt")))
  }
  
  finish_time <- base::Sys.time()
  base::cat(base::sprintf("\nLPSN Run Complete at: %s\n", base::format(finish_time, "%Y-%m-%d %H:%M:%S")))
}


###example

# # Your cleaned taxa list (can include synonyms, candidatus, or higher ranks)
# my_taxa <- c("Firmicutes",..., "Clostridium")

# # Output directory structure
# CSV_DIR <- "./LPSN/results/csv_lineages"
# WB_PATH <- "./LPSN/results/LPSN_result_Master.xlsx"
# LOG_DIR <- "./LPSN/results/logs"
# 
# # run API function
# # set batch size as needed. Small batches is a bit slower as it opens and saves to excel more often
# 
# process_lpsn_cached_batches(
#   taxa_list     = my_taxa_clean,
#   csv_path      = CSV_DIR,
#   workbook_path = WB_PATH,
#   log_path      = LOG_DIR,
#   user          = Sys.getenv("DSMZ_API_USER"),
#   password      = Sys.getenv("DSMZ_API_PASSWORD"),
#   batch_size    = 50,    # Processes 50 unique taxa per Excel sheet
#   pause_min     = 0.5,   # Minimum delay between API calls (seconds)
#   pause_max     = 1.5,   # Maximum delay between API calls (seconds)
#   use_timestamp = TRUE   # Appends date/time to the Excel filename
# )



###### application of lineage retrieved from LPSN API to taxonomy file

#' --- Prerequisites ---
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
  
  # 1. Load Master Data based on file extension
  # Detects file extension to use the correct loading function
  ext <- base::tolower(tools::file_ext(master_path))
  
  if (ext == "csv") {
    master <- utils::read.csv(master_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    master <- openxlsx::read.xlsx(master_path)
  } else {
    base::stop("Master file must be .csv, .xlsx, or .xls")
  }
  
  # 2. Load Taxonomy Data
  # Loads the target file; assumes the first column is the row identifier
  taxonomy <- utils::read.csv(taxonomy_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  # --- Case-Insensitive Column Mapping ---
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
  
  # --- 3. Initialize Result Columns ---
  # Create a column to track which rank/term was successfully matched
  taxonomy$matched_search_term <- "none"
  for (col in mapped_append_cols) {
    taxonomy[[base::paste0("master_", col)]] <- NA_character_
  }
  
  # --- 4. Row-wise Processing ---
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
  
  # --- 5. Export ---
  utils::write.csv(taxonomy, output_path, row.names = TRUE)
  base::cat(base::sprintf("Done. Results in: %s\n", output_path))
}







##########################  
##### modified function to test - candidatus now in "" - check if function returns the pro-correct name for candidatus

# --- Prerequisites ---
# Ensure you have your LPSN credentials set in your environment (.Renviron)
# DSMZ_API_USER="your_email@example.com"
# DSMZ_API_PASSWORD="your_password"

#' Process Taxa via LPSN API in Batches
#'
#' @param taxa_list Character vector of taxonomic names.
#' @param csv_path Directory for individual lineage CSVs.
#' @param workbook_path Path to the output .xlsx file.
#' @param log_path Directory for error/deduplication logs.
#' @param user DSMZ API Username.
#' @param password DSMZ API Password.
#' @param batch_size Number of taxa per batch (default 20).
#' @param pause_min Min pause between API calls (default 1).
#' @param pause_max Max pause between API calls (default 3).
#' @param use_timestamp Append timestamp to workbook name (default TRUE).

process_lpsn_cached_batches_full <- function(taxa_list, csv_path, workbook_path, log_path, user, password, batch_size = 20, pause_min = 1, pause_max = 3, use_timestamp = TRUE) {
  start_time <- base::Sys.time()
  
  # --- 1. Deduplication & Name Cleaning ---
  if (!base::dir.exists(log_path)) base::dir.create(log_path, recursive = TRUE)
  original_count <- base::length(taxa_list)
  
  unique_taxa <- base::as.character(taxa_list)
  unique_taxa <- base::trimws(unique_taxa) # Strip whitespace
  unique_taxa <- base::gsub("\\s+", " ", unique_taxa) # Collapse double spaces
  unique_taxa <- base::unique(unique_taxa)
  unique_taxa <- unique_taxa[!base::is.na(unique_taxa) & unique_taxa != ""]
  
  duplicates <- base::setdiff(taxa_list, unique_taxa)
  log_content <- if (base::length(duplicates) > 0) base::sort(base::unique(duplicates)) else "No duplicates found."
  base::writeLines(log_content, base::file.path(log_path, "deduplicated_taxa.txt"))
  
  taxa_list <- unique_taxa
  new_count <- base::length(taxa_list)
  base::cat(base::sprintf("Input: %d taxa. Unique/Cleaned: %d.\n", original_count, new_count))
  
  # --- 2. Setup Environments & Helpers ---
  api_cache <- base::list()
  if (use_timestamp) {
    ts <- base::format(base::Sys.time(), "%Y%m%d_%H%M")
    workbook_path <- base::gsub("\\.xlsx$", base::paste0("_", ts, ".xlsx"), workbook_path)
  }
  if (!base::dir.exists(csv_path)) base::dir.create(csv_path, recursive = TRUE)
  
  RANK_ORDER <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  HIGHER_RANKS <- c("phylum", "class", "order", "family", "genus", "kingdom", "domain")
  export_fields <- c("id", "category", "full_name", "lpsn_parent_id", "lpsn_correct_name_id", "lpsn_taxonomic_status")
  
  safe_rbind <- function(df_list) {
    if (base::length(df_list) == 0) return(NULL)
    all_cols <- base::unique(base::unlist(base::lapply(df_list, base::names)))
    res_list <- base::lapply(df_list, function(df) {
      for (col in base::setdiff(all_cols, base::names(df))) df[[col]] <- NA
      df[all_cols]
    })
    base::do.call(base::rbind, res_list)
  }
  
  retry_api_call <- function(expr, retries = 3, wait = 5) {
    for (i in 1:retries) {
      res <- base::tryCatch(expr, error = function(e) NULL)
      if (!base::is.null(res)) return(res)
      if (i < retries) {
        base::message(base::sprintf("\n API error. Retry %d/%d...", i, retries))
        base::Sys.sleep(wait)
      }
    }
    return(NULL)
  }
  
  find_direct_hit <- function(query_res, target_name, target_category) {
    df_res <- base::as.data.frame(query_res)
    if (base::nrow(df_res) == 0) return(NULL)
    # Note: Comparison is done against raw name (without quotes)
    hit_idx <- base::which(base::tolower(df_res$full_name) == base::tolower(target_name) & df_res$category == target_category)
    if (base::length(hit_idx) > 0) return(df_res[hit_idx, , drop = FALSE])
    return(NULL)
  }
  
  # Wrapper for Species/Genus level only: Uses QUOTED search terms and retry logic
  search_with_candidatus_retry <- function(lpsn_obj, name, rank) {
    # Wrap name in double quotes for API: "Name"
    quoted_name <- base::sprintf('"%s"', name)
    
    base::Sys.sleep(stats::runif(1, pause_min, pause_max))
    raw <- retry_api_call(LPSN::retrieve(lpsn_obj, search = "advanced", `taxon-name` = quoted_name, category = rank))
    hit <- find_direct_hit(raw, name, rank)
    
    if (base::is.null(hit) && !base::grepl("^Candidatus", name, ignore.case = TRUE)) {
      cand_name <- base::paste("Candidatus", name)
      quoted_cand_name <- base::sprintf('"%s"', cand_name)
      
      base::message(base::sprintf("\n '%s' not found. Retrying as '%s'...", name, quoted_cand_name))
      base::Sys.sleep(stats::runif(1, pause_min, pause_max))
      raw_cand <- retry_api_call(LPSN::retrieve(lpsn_obj, search = "advanced", `taxon-name` = quoted_cand_name, category = rank))
      hit <- find_direct_hit(raw_cand, cand_name, rank)
      if (!base::is.null(hit)) match_type <<- base::paste0(match_type, " (Candidatus)")
    }
    return(hit)
  }
  
  lpsn <- LPSN::open_lpsn(user, password)
  num_batches <- base::ceiling(new_count / batch_size)
  pb <- progress::progress_bar$new(format = " [:bar] :percent Batch :batch/:total_batches ETA: :eta", total = new_count, width = 80)
  
  # --- 3. Main Loop ---
  for (b in 1:num_batches) {
    start_idx <- ((b - 1) * batch_size) + 1
    end_idx <- base::min(b * batch_size, new_count)
    current_batch_taxa <- taxa_list[start_idx:end_idx]
    batch_workbook_rows <- list()
    batch_not_found <- character()
    
    for (taxon in current_batch_taxa) {
      pb$tick(tokens = base::list(batch = b, total_batches = num_batches))
      match_type <- "Full Match"
      res <- NULL
      is_binomial <- base::grepl(" ", taxon)
      
      if (is_binomial) {
        res <- search_with_candidatus_retry(lpsn, taxon, "species")
        if (base::is.null(res)) {
          genus_only <- if (base::grepl("^Candidatus", taxon, ignore.case = TRUE)) {
            base::sub("^([^ ]+ +[^ ]+).*", "\\1", taxon)
          } else {
            base::sub(" .*", "", taxon)
          }
          match_type <- "Genus Fallback"
          res <- search_with_candidatus_retry(lpsn, genus_only, "genus")
        }
      } else {
        # Higher Ranks: Standard search using QUOTES (no Candidatus retry)
        for (rank_to_check in HIGHER_RANKS) {
          quoted_taxon <- base::sprintf('"%s"', taxon)
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          raw_res <- retry_api_call(LPSN::retrieve(lpsn, search = "advanced", `taxon-name` = quoted_taxon, category = rank_to_check))
          res <- find_direct_hit(raw_res, taxon, rank_to_check)
          if (!base::is.null(res)) {
            match_type <- base::paste0("Full Match (", tools::toTitleCase(rank_to_check), ")")
            break
          }
        }
      }
      
      if (base::is.null(res)) {
        batch_not_found <- base::c(batch_not_found, base::as.character(taxon))
        next
      }
      
      # --- 4. Lineage Retrieval & Status Handling ---
      first_hit <- base::as.data.frame(res)[1, , drop = FALSE]
      
      # Handle Candidatus statuses like "pro-correct name"
      nom_status <- if (base::length(first_hit$lpsn_taxonomic_status) > 0) {
        first_hit$lpsn_taxonomic_status
      } else {
        "Unknown"
      }
      
      is_synonym <- !base::is.null(first_hit$lpsn_correct_name_id) && 
        !base::is.na(first_hit$lpsn_correct_name_id) && 
        first_hit$lpsn_correct_name_id != 0 && 
        first_hit$lpsn_correct_name_id != first_hit$id
      
      current_id <- if (is_synonym) first_hit$lpsn_correct_name_id else first_hit$id
      lineage_list <- list()
      corrected_name <- NA
      
      while (!base::is.null(current_id) && !base::is.na(current_id) && current_id != 0) {
        id_key <- base::as.character(current_id)
        if (id_key %in% base::names(api_cache)) {
          rec_df <- api_cache[[id_key]]
        } else {
          base::Sys.sleep(stats::runif(1, pause_min, pause_max))
          rec_res <- retry_api_call(LPSN::fetch(lpsn, ids = current_id))
          if (base::is.null(rec_res)) break
          rec_df <- base::as.data.frame(rec_res)[1, , drop = FALSE]
          api_cache[[id_key]] <- rec_df
        }
        
        # Safety Check for zero-length logical
        if (base::length(rec_df$category) > 0 && !base::is.na(rec_df$category) && rec_df$category != "division") {
          lineage_list[[base::length(lineage_list) + 1]] <- rec_df[, base::intersect(export_fields, base::names(rec_df)), drop = FALSE]
        }
        
        if (base::length(lineage_list) == 1) corrected_name <- rec_df$full_name
        current_id <- if (base::length(rec_df$lpsn_parent_id) > 0 && rec_df$lpsn_parent_id != 0) rec_df$lpsn_parent_id else NULL
      }
      
      # --- 5. Export Logic (CSVs) ---
      if (base::length(lineage_list) > 0) {
        lineage_data <- safe_rbind(lineage_list)
        safe_name <- base::gsub("[^[:alnum:]]", "_", taxon)
        utils::write.csv(lineage_data, base::file.path(csv_path, base::paste0("Batch_", b, "_", safe_name, ".csv")), row.names = FALSE)
        
        wb_row <- stats::setNames(base::as.list(base::rep(NA, base::length(RANK_ORDER))), RANK_ORDER)
        wb_row$search_term <- taxon
        wb_row$match_type <- if (is_synonym) base::paste0(match_type, " (Synonym Redirect)") else match_type
        wb_row$nomenclatural_status <- nom_status
        wb_row$corrected_name <- corrected_name
        
        for (i in 1:base::nrow(lineage_data)) {
          rank <- lineage_data$category[i]
          if (!base::is.null(rank) && !base::is.na(rank) && rank %in% RANK_ORDER) wb_row[[rank]] <- lineage_data$full_name[i]
        }
        batch_workbook_rows[[taxon]] <- base::as.data.frame(wb_row)
      }
    }
    
    # --- 6. Save Excel with tryCatch protection ---
    if (base::length(batch_workbook_rows) > 0) {
      final_batch_df <- safe_rbind(batch_workbook_rows)
      sheet_name <- base::paste0("Batch_", b)
      
      base::tryCatch({
        wb <- if (base::file.exists(workbook_path)) openxlsx::loadWorkbook(workbook_path) else openxlsx::createWorkbook()
        if (sheet_name %in% base::names(wb)) openxlsx::removeWorksheet(wb, sheet_name)
        openxlsx::addWorksheet(wb, sheetName = sheet_name)
        
        cols <- base::c("search_term", "match_type", "nomenclatural_status", "corrected_name", RANK_ORDER)
        for (col_name in base::setdiff(cols, base::names(final_batch_df))) final_batch_df[[col_name]] <- NA
        
        openxlsx::writeData(wb, sheet = sheet_name, x = final_batch_df[, cols])
        openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
      }, error = function(e) {
        base::message(base::sprintf("\n[ERROR] Could not save Batch %d to Excel. Is '%s' open in another program?", b, base::basename(workbook_path)))
      })
    }
    
    base::writeLines(base::sort(base::as.character(base::unlist(batch_not_found))), 
                     base::file.path(log_path, base::paste0("taxa_not_found_batch_", b, ".txt")))
  }
  
  duration <- base::difftime(base::Sys.time(), start_time, units = "auto")
  base::cat(base::sprintf("\nRun Complete. Total time: %s\n", base::format(duration)))
}

  
###example

# # Your cleaned taxa list (can include synonyms, candidatus, or higher ranks)
# my_taxa <- c("Firmicutes",..., "Clostridium")

# # Output directory structure
# CSV_DIR <- "./LPSN/results/csv_lineages"
# WB_PATH <- "./LPSN/results/LPSN_result_Master.xlsx"
# LOG_DIR <- "./LPSN/results/logs"
# 
# # run API function
# # set batch size as needed. Small batches is a bit slower as it opens and saves to excel more often
# 
# 
# process_lpsn_cached_batches_full(
#   taxa_list     = my_taxa_clean,
#   csv_path      = CSV_DIR,
#   workbook_path = WB_PATH,
#   log_path      = LOG_DIR,
#   user          = Sys.getenv("DSMZ_API_USER"),
#   password      = Sys.getenv("DSMZ_API_PASSWORD"),
#   batch_size    = 50,    # Processes 50 unique taxa per Excel sheet
#   pause_min     = 0.5,   # Minimum delay between API calls (seconds)
#   pause_max     = 1.5,   # Maximum delay between API calls (seconds)
#   use_timestamp = TRUE   # Appends date/time to the Excel filename
# )



