#---------------------------
#The purpose of the scripts below is generate a summary table from a spreadsheet of specimens that are being considered for accession to KU. The script also generates a new version of the original spreadsheet that is clear of common formatting errors and includes additional columns to indicate the protected status of each individual specimen. In addition to providing numbers of specimens in several categories (e.g., frogs, snakes) the resulting summary table identifies which specimens are domestic versus international and which species are protected under either CITES or US ESA. Protected taxa are identified by querying the accession spreadsheet against databases of taxa that are protected either by CITES or the US FWS's Endangered Species Act.
#---------------------------
library(readr)
library(dplyr)
library(stringr)
library(taxize)
library(purrr)
library(tibble)
library(tidyverse)

rm(list=ls())
#----------------------------
#UPLOAD LISTS IF PROTECTED TAXA
#----------------------------
#Obtain current list of CITES protected species from https://checklist.cites.org/#/en. Download by clicking the "Download" link and selecting CSV option. 

cites <- read.csv(
  "Index_of_CITES_Species_2025-12-16 13_47.csv",
  stringsAsFactors = FALSE
)

# restrict list of taxa to cases where Class is Reptilia or Amphibia
cites_herps <- cites %>% 
  filter(Class %in% c("Reptilia", "Amphibia"))

#Obtain current list of FWS ESA protected species. You can obtain the most recent data for all vertebrates from the US by going to https://ecos.fws.gov/ecp/report/species-listings-by-tax-group-totals clicking on the "All Vertebrate Species" link and selecting the "CSV" download option. You will also need to download a second list that includes foreign species that are protected by the USFWS by going to https://ecos.fws.gov/ecp/species-reports and clicking on the Foreign Species link. Unlike the US list, the list of foreign species is comprehensive of all taxonomic groups and can be downloaded by clicking the CSV link on the page.

# read the file containing FWS species from US into R
fws <- read.csv(
  "species-listings-by-tax-group-report_251217.csv",
  stringsAsFactors = FALSE
)

# restrict vertebrate lest to just reptiles and amphibians
fws_herps <- fws %>%
  filter(str_detect(tolower(Group), "reptile|amphib"))

#split the column with binomial names in the downloaded spreadsheet into separate columns for genus and species
fws_herps <- fws_herps %>%
  mutate(
    Scientific_Name_clean = str_squish(`Scientific.Name`),
    genus   = word(Scientific_Name_clean, 1),
    species = word(Scientific_Name_clean, 2),
    binomial = paste(genus, species)
  )

# read the file containing FWS international species into R
fws_international <- read_csv("species-listings-foreign-report_251217.csv", show_col_types = FALSE)

# process and qc of FWS international list
sci_col_int <- names(fws_international)[str_detect(tolower(names(fws_international)), "scientific")][1]
status_col_int <- names(fws_international)[str_detect(tolower(names(fws_international)), "status")][1]

if (is.na(sci_col_int)) stop("Couldn't find a Scientific Name column in fws_international (expected name containing 'scientific').")
if (is.na(status_col_int)) stop("Couldn't find a Status column in fws_international (expected name containing 'status').")

fws_int_parsed <- fws_international %>%
  mutate(
    Scientific_Name_clean = str_squish(as.character(.data[[sci_col_int]])),
    genus   = word(Scientific_Name_clean, 1),
    species = word(Scientific_Name_clean, 2)
  )

#----------------------------
#UPLOAD SPECIMEN LIST FOR ACCESSION
#----------------------------
#read in specimen data from UMC
umc <- read.csv(
  "UMC_Herpetology_Combined_NoBlanks_251216.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

#fix formatting error involving extra spaces in the Order column
umc$Order <- trimws(gsub("\\s+", " ", umc$Order))
unique(umc$Order)

#ad column that includes possible alternative species names that might vary depending on recent taxonomic revisions
umc <- umc %>%
  mutate(
    alt_genus = case_when(
      Genus == "Rana"    ~ "Lithobates",
      Genus == "Bufo"    ~ "Anaxyrus",
      Genus == "Eumeces" ~ "Plestiodon",
      TRUE               ~ NA_character_
    )
  )

#----------------------------
#IDENTIFY TAXA THAT ARE SHARED BETWEEN SPECIMEN LIST AND PROTECTED LISTS AND APPEND INFORMATION TO SPECIMEN TABLE
#----------------------------

#############################
# Identify CITES protected taxa
#############################

# helper to clean text consistently
clean_txt <- function(x) str_trim(str_to_lower(as.character(x)))

#clean up content in priority column
priority <- function(x) {
  x2 <- str_to_upper(str_trim(x))
  case_when(
    str_detect(x2, "APPENDIX\\s*I\\b|\\bI\\b")   ~ 3L,
    str_detect(x2, "APPENDIX\\s*II\\b|\\bII\\b") ~ 2L,
    str_detect(x2, "APPENDIX\\s*III\\b|\\bIII\\b")~ 1L,
    TRUE ~ 0L
  )
}

# Generate look-up table that includes protected genera and their CITES status. This table is generated by identifying rows that specifically reference a protected genus.
cites_genus_lookup <- cites %>%
  mutate(
    rank_std  = str_to_upper(str_trim(RankName)),
    genus_std = str_to_lower(str_trim(Genus))
  ) %>%
  filter(rank_std == "GENUS", !is.na(genus_std), genus_std != "") %>%
  mutate(listing_priority = priority(CurrentListing)) %>%
  arrange(desc(listing_priority)) %>%
  distinct(genus_std, .keep_all = TRUE) %>%
  transmute(genus_std, cites_CurrentListing = CurrentListing)

umc_out <- umc %>%
  mutate(
    genus_std     = tolower(trimws(Genus)),
    alt_genus_std = tolower(trimws(alt_genus))
  ) %>%
  left_join(
    cites_genus_lookup,
    by = c("genus_std" = "genus_std")
  ) %>%
  left_join(
    cites_genus_lookup,
    by = c("alt_genus_std" = "genus_std"),
    suffix = c("", "_alt")
  ) %>%
  mutate(
    genus_CITES_listing = case_when(
      !is.na(cites_CurrentListing)     ~ cites_CurrentListing,
      !is.na(cites_CurrentListing_alt) ~ cites_CurrentListing_alt,
      TRUE                             ~ "not listed"
    )
  ) %>%
  select(-genus_std, -alt_genus_std,
         -cites_CurrentListing, -cites_CurrentListing_alt)


# Generate look-up table that includes protected species and their CITES status. This table is generated by identifying rows based on Family and Species names to avoid potential issues with generic revision.
cites_fs_lookup <- cites %>%
  mutate(
    family_std  = str_to_lower(str_trim(Family)),
    species_std = str_to_lower(str_trim(Species))
  ) %>%
  filter(!is.na(family_std), family_std != "", !is.na(species_std), species_std != "") %>%
  mutate(listing_priority = priority(CurrentListing)) %>%
  arrange(desc(listing_priority)) %>%
  distinct(family_std, species_std, .keep_all = TRUE) %>%
  transmute(family_std, species_std, cites_species_listing = CurrentListing)

umc_out <- umc_out %>%
  mutate(
    family_std  = str_to_lower(str_trim(Family)),
    species_std = str_to_lower(str_trim(Species))
  ) %>%
  left_join(cites_fs_lookup, by = c("family_std", "species_std")) %>%
  mutate(
    species_CITES_listing = if_else(is.na(cites_species_listing), "not listed", cites_species_listing)
  ) %>%
  select(-family_std, -species_std, -cites_species_listing)

#############################
# Identify USFWS protected taxa from the US
#############################
# find the status column in fws_herps
status_col_us <- names(fws_herps)[str_detect(tolower(names(fws_herps)), "status")][1]
if (is.na(status_col_us)) stop("Couldn't find a status column in fws_herps (expected 'status').")

fws_us_lookup <- fws_herps %>%
  mutate(
    fws_genus_std   = str_to_lower(str_trim(genus)),
    fws_species_std = str_to_lower(str_trim(species)),
    fws_status_us   = str_squish(as.character(.data[[status_col_us]]))
  ) %>%
  filter(fws_genus_std != "", fws_species_std != "") %>%
  distinct(fws_genus_std, fws_species_std, .keep_all = TRUE) %>%
  select(fws_genus_std, fws_species_std, fws_status_us)

umc_out2 <- umc_out %>%
  mutate(
    genus_std     = str_to_lower(str_trim(Genus)),
    alt_genus_std = str_to_lower(str_trim(alt_genus)),
    species_std   = str_to_lower(str_trim(Species))
  ) %>%
  left_join(fws_us_lookup,
            by = c("genus_std" = "fws_genus_std", "species_std" = "fws_species_std")) %>%
  left_join(fws_us_lookup,
            by = c("alt_genus_std" = "fws_genus_std", "species_std" = "fws_species_std"),
            suffix = c("", "_alt")) %>%
  mutate(
    fws_status_us = coalesce(fws_status_us, fws_status_us_alt, "not listed")
  ) %>%
  select(-genus_std, -alt_genus_std, -species_std, -fws_status_us_alt)

#############################
# Identify USFWS international protected taxa
#############################

fws_int_lookup <- fws_int_parsed %>%
  mutate(
    fws_genus_std   = str_to_lower(str_trim(genus)),
    fws_species_std = str_to_lower(str_trim(species)),
    fws_status_international = str_squish(as.character(.data[[status_col_int]]))
  ) %>%
  filter(fws_genus_std != "", fws_species_std != "") %>%
  distinct(fws_genus_std, fws_species_std, .keep_all = TRUE) %>%
  select(fws_genus_std, fws_species_std, fws_status_international)

umc_out3 <- umc_out2 %>%
  mutate(
    genus_std     = str_to_lower(str_trim(Genus)),
    alt_genus_std = str_to_lower(str_trim(alt_genus)),
    species_std   = str_to_lower(str_trim(Species))
  ) %>%
  left_join(fws_int_lookup,
            by = c("genus_std" = "fws_genus_std", "species_std" = "fws_species_std")) %>%
  left_join(fws_int_lookup,
            by = c("alt_genus_std" = "fws_genus_std", "species_std" = "fws_species_std"),
            suffix = c("", "_alt")) %>%
  mutate(
    fws_status_international =
      coalesce(fws_status_international, fws_status_international_alt, "not listed")
  ) %>%
  select(-genus_std, -alt_genus_std, -species_std, -fws_status_international_alt)


# normalize Order text for consistent matching
umc_out3 <- umc_out3 %>%
  mutate(
    Order_clean   = str_to_lower(str_squish(as.character(Order))),
    Country_clean = str_to_upper(str_squish(as.character(Country)))
  )

# identify CITES and FWS columns by name
cites_cols <- grep("cites", names(umc_out3), ignore.case = TRUE, value = TRUE)
fws_cols   <- grep("fws",   names(umc_out3), ignore.case = TRUE, value = TRUE)

# helper: counts for one subset
summarize_group <- function(df) {
  if (nrow(df) == 0) {
    return(tibble(
      total_records = 0L,
      usa_records   = 0L,
      nonusa_records= 0L,
      cites_listed  = 0L,
      fws_listed    = 0L
    ))
  }
  
  # CITES listed = any CITES-named column has value other than "not listed" (case-insensitive) and not blank/NA
  cites_listed_vec <- if (length(cites_cols) == 0) {
    rep(FALSE, nrow(df))
  } else {
    mat <- df[, cites_cols, drop = FALSE] %>% mutate(across(everything(), as.character))
    apply(mat, 1, function(r) any(!is.na(r) & str_squish(str_to_lower(r)) != "" & str_squish(str_to_lower(r)) != "not listed"))
  }
  
  # FWS listed = any FWS-named column has value other than "not listed" and not blank/NA
  fws_listed_vec <- if (length(fws_cols) == 0) {
    rep(FALSE, nrow(df))
  } else {
    mat <- df[, fws_cols, drop = FALSE] %>% mutate(across(everything(), as.character))
    apply(mat, 1, function(r) any(!is.na(r) & str_squish(str_to_lower(r)) != "" & str_squish(str_to_lower(r)) != "not listed"))
  }
  
  tibble(
    total_records  = nrow(df),
    usa_records    = sum(df$Country_clean == "USA", na.rm = TRUE),
    nonusa_records = sum(df$Country_clean != "USA" & !is.na(df$Country_clean) & df$Country_clean != "", na.rm = TRUE),
    cites_listed   = sum(cites_listed_vec),
    fws_listed     = sum(fws_listed_vec)
  )
}

# define groups (order-based filters)
groups <- list(
  frogs       = umc_out3 %>% filter(Order_clean == "anura"),
  salamanders = umc_out3 %>% filter(Order_clean == "caudata"),
  caecilians  = umc_out3 %>% filter(Order_clean == "gymnophiona"),
  amphibians  = umc_out3 %>% filter(Order_clean %in% c("anura","caudata","gymnophiona")),
  lizards     = umc_out3 %>% filter(str_detect(Order_clean, "lacertilia")),
  snakes      = umc_out3 %>% filter(str_detect(Order_clean, "serpentes")),
  turtles     = umc_out3 %>% filter(str_detect(Order_clean, "testudines")),
  crocs       = umc_out3 %>% filter(str_detect(Order_clean, "crocodilia"))
)

# build the first 8 rows
sum_tbl <- bind_rows(lapply(groups, summarize_group), .id = "group")

# reptiles total = lizards + snakes + turtles + crocs
reptiles_total <- sum_tbl %>%
  filter(group %in% c("lizards","snakes","turtles","crocs")) %>%
  summarise(across(where(is.numeric), base::sum)) %>%
  mutate(group = "reptiles")

# totals = across all rows in original
totals <- summarize_group(umc_out3) %>% mutate(group = "totals")

# final 10-row table
sum_tbl <- bind_rows(
  sum_tbl %>% filter(group %in% c("frogs","salamanders","caecilians","amphibians",
                                  "lizards","snakes","turtles","crocs")),
  reptiles_total,
  totals
)

sum_tbl

write.csv(
  sum_tbl,
  "sum_tbl.csv",
  row.names = FALSE
)

write.csv(
  umc_out3,
  "umc_out3.csv",
  row.names = FALSE
)

#OPTIONAL STEPS TO TAXONOMICALLY DEFINE ORGANISMS IN FWS LISTS
#------------
#lines below are to ad taxonomic information about family and order, but this probabaly isn't necessary
#function to obtain data on order and family and append to table
get_order_family_itis <- function(sciname) {
  res <- tryCatch({
    cl <- taxize::classification(sciname, db = "itis")[[1]]
    if (is.null(cl) || nrow(cl) == 0)
      return(tibble(order = NA_character_, family = NA_character_))
    
    tibble(
      order  = cl$name[cl$rank == "order"][1],
      family = cl$name[cl$rank == "family"][1]
    )
  }, error = function(e) {
    tibble(order = NA_character_, family = NA_character_)
  })
  res
}

# obtain order and family data by querying file against the taxonomic database in the R function "taxize". Note that this depends on taxize being up to date and the species binomials in your original list having matches to this database. This function takes a while and probably isn't even necessary.
tax_tbl <- fws_herps %>%
  distinct(binomial) %>%
  mutate(tax = map(binomial, get_order_family_itis)) %>%
  unnest(tax)

fws_out <- fws_herps %>%
  left_join(tax_tbl, by = "binomial") %>%
  # drop helper columns
  select(-Scientific_Name_clean, -binomial)
#------------