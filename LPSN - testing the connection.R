# Test LPSN package to check the functions work  ##############################

# Examples provided by LPSN package documentation

## Registration for LPSN is required but free and easy to accomplish.
## In real applications username and password could of course also be stored
## within the R code itself, or read from a file.

#### See code sections in test for different search and retrieval options


## Load Credentials and create the LPSN access object  ####

credentials <- Sys.getenv(c("DSMZ_API_USER", "DSMZ_API_PASSWORD"))

if (all(nzchar(credentials))) {
  ## create the LPSN access object
  lpsn <- open_lpsn(credentials[[1L]], credentials[[2L]])
  print(lpsn)
  # it would be frustrating if the object was already expired on creation
  stopifnot(
    inherits(lpsn, "lpsn_access"),
    !summary(lpsn)[c("expired", "refresh_expired")]
  )

  ## Using LPSN IDs (record numbers) to fetch records ##############################

  ### Example 1: (id1) each ID given as separate argument

  id1 <- fetch(lpsn, 520424, 4948, 17724)
  print(id1)
  stopifnot(
    inherits(id1, "lpsn_result"),
    summary(id1)[["count"]] == 3L
  )
  # conversion to data frame to view results
  id1d <- as.data.frame(id1)
  head(id1d)
  stopifnot(is.data.frame(id1d), nrow(id1d) == 3L)

  ### Example 2: (id2) all IDs in vector

  id2 <- fetch(lpsn, c(520424, 4948, 17724))
  stopifnot(identical(id1, id2))

  ### Example 3: (id3) all IDs in a list

  id3 <- fetch(lpsn, list(520424, 4948, 17724))
  stopifnot(identical(id1, id3))

  ## Run flexible search:  ########################################################

  ## Flexible search is flexible because one can query for
  ## all key-value pairs available in some API entry, but
  ## the values are matched exactly,
  ## Flexible search also allows for nested queries

  ### Example 4: (fs1) each part of the query given as separate argument

  fs1 <- request(lpsn, search = "flexible", monomial = "Flexithrix")
  print(fs1)
  stopifnot(
    inherits(fs1, "lpsn_result"),
    summary(fs1)[["count"]] >= 2L
  )

  # conversion to data frame is possible but may not contain all
  # data if the API outcome was paginated; see 'retrieve' below

  fs1d <- as.data.frame(fs1)
  head(fs1d)
  stopifnot(is.data.frame(fs1d), dim(fs1d) > 0L)


  ### Example 5: (fs2) all parts of the query given in vector

  fs2 <- request(lpsn, c(monomial = "Flexithrix"), "flexible")
  stopifnot(identical(fs1, fs2))


  ### Example 6: (fs3) all parts of the query given in list

  fs3 <- request(lpsn, list(monomial = "Flexithrix"), "flexible")
  stopifnot(identical(fs1, fs3))



  ## Run advanced search:  ########################################################

  ## Advanced search uses the equivalent of a predefined set of
  ## keys (or combinations thereof) to query. The values are not
  ## matched exactly only; case is not distinguished and substrings
  ## are recognized.
  ## Advanced search does not permit nested queries

  ### Example 7: (as1) each part of the query given as separate argument

  as1 <- request(lpsn, search = "advanced", `taxon-name` = "Flexithrix")
  print(as1)
  stopifnot(
    inherits(as1, "lpsn_result"),
    summary(as1)[["count"]] >= 4L
  )

  # conversion to data frame is possible but may here not contain all
  # data if the API outcome was paginated; see 'retrieve' below

  as1d <- as.data.frame(as1)
  head(as1d)
  stopifnot(is.data.frame(as1d), dim(as1d) > 0L)

  ### Example 8: (as2) all parts of the query given in vector

  as2 <- request(lpsn, c(`taxon-name` = "Flexithrix"), "advanced")
  stopifnot(identical(as1, as2))


  ### Example 9: (as3) all parts of the query given in list

  as3 <- request(lpsn, list(`taxon-name` = "Flexithrix"), "advanced")
  stopifnot(identical(as1, as3))


  ## Run search + fetch in one step:   #############################################

  ### Example 10: (fsf) via flexible search

  fsf <- retrieve(lpsn, search = "flexible", monomial = "Flexithrix")
  stopifnot(
    inherits(fsf, "records"),
    length(fsf) == summary(fs1)[["count"]]
  )
  # Convert to a data frame
  fsfd <- as.data.frame(fsf)
  stopifnot(is.data.frame(fsfd), nrow(fsfd) == length(fsf))
  print(dim(fsfd))


  ### Example 11: (asf) via advanced search

  asf <- retrieve(lpsn, search = "advanced", `taxon-name` = "Flexithrix")
  stopifnot(length(asf) == summary(as1)[["count"]],
    identical(class(asf), class(fsf)))
  # convert to data frame possible
  asfd <- as.data.frame(asf)
  stopifnot(is.data.frame(asfd), nrow(asfd) == length(asf))
  print(dim(asfd))

  ### Example 12: (asfh) use a handler
  asfh <- list()
  retrieve(lpsn, search = "advanced", `taxon-name` = "Flexithrix",
    handler = function(x) asfh <<- c(asfh, x))
  stopifnot(length(asfh) == length(asf))
  # there are of course better ways to use a handler

  ### Example 13: (nil) if nothing is found
  nil <- retrieve(lpsn, search = "flexible", monomial = "unknown")
  stopifnot(length(nil) == 0L,
    identical(class(nil), class(fsf)))
  # convert to data frame possible
  nild <- as.data.frame(nil)
  stopifnot(is.data.frame(nild), nrow(nild) == length(nil))
  print(dim(nild))

  ### Example 14: Lots of records retuned
  bac <- retrieve(lpsn, list(monomial = "Bacillus"),
    sleep = 0.1) # just for the sake of testing
  stopifnot(length(bac) > 400L,
    identical(class(bac), class(fsf)))
  # onversion to data frame possible
  bacdf <- as.data.frame(bac)
  stopifnot(is.data.frame(bacdf),
    nrow(bacdf) == length(bac))

  ## Refresh API call ####
  refresh(lpsn, TRUE)
  # this is also done internally and automatically
  # in some situations when apparently needed

} else {

  warning("username or password missing, cannot run examples")

}

# Test output

# lpsn_access                             TRUE
# dsmz_keycloak                           TRUE
# expired                                 FALSE
# refresh_expired                         FALSE
# class                                   lpsn_result < dsmz_result
# parts                                   count, next, previous, results
# count                                   3
# next                                    FALSE
# previous                                FALSE
# results                                 TRUE
# class                                   lpsn_result < dsmz_result
# parts                                   count, next, previous, results
# count                                   2
# next                                    FALSE
# previous                                FALSE
# results                                 TRUE
# class                                   lpsn_result < dsmz_result
# parts                                   count, next, previous, results
# count                                   4
# next                                    FALSE
# previous                                FALSE
# results                                 TRUE
# [1]  2 27
# [1]  4 27
# [1] 0 0
# lpsn_access                             TRUE
# dsmz_keycloak                           TRUE
# expired                                 FALSE
# refresh_expired                         FALSE



# Search (id1) using ids returns following variables.
# Convert to data frame.

# res_data <- as.data.frame(id1)
# names(res_data)
#
# [1] "id"                    "category"              "monomial"
# [4] "authority"             "full_name"             "proposed_as"
# [7] "lpsn_address"          "lpsn_version"          "is_legitimate"
# [10] "ijsem_list_doi"        "lpsn_parent_id"        "ijsem_list_kind"
# [13] "ijsem_list_text"       "publication_doi"       "publication_kind"
# [16] "publication_pmid"      "publication_text"      "validly_published"
# [19] "lpsn_correct_name_id"  "nomenclatural_status"  "is_spelling_corrected"
# [22] "lpsn_taxonomic_status" "nomenclatural_type_id"
#
# t(head(res_data)[1, ])
# 1
# id                    "520424" ## id for taxonomic name
# category              "genus"  ## Taxonomic rank
# monomial              "Abditibacterium"  ## single name
# authority             "Tahon et al. 2018"
# full_name             "Abditibacterium"
# proposed_as           "gen. nov."
# lpsn_address          "https://doi.org/10.83108/rn.520424"
# lpsn_version          "1"
# is_legitimate         "TRUE"
# ijsem_list_doi        "10.1099/ijsem.0.003071"
# lpsn_parent_id        "4812"
# ijsem_list_kind       "Validation List"
# ijsem_list_text       "Oren A, Garrity GM. Validation list no. 184. ..."
# publication_kind      "effective"
# publication_pmid      "29475572"
# publication_text      "Tahon G, Tytgat B, Lebbe L, Carlier A, Willems A. ...."
# validly_published     "ICNP"
# lpsn_correct_name_id  "520424"
# nomenclatural_status  "validly published under the ICNP"
# is_spelling_corrected "FALSE"
# lpsn_taxonomic_status "correct name"
# nomenclatural_type_id "797965"
