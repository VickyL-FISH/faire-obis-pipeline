# =============================================================================
# DFO Northern Shelf Bioregion eDNA — OBIS Submission Pipeline
# =============================================================================
#
# Purpose:   Standardise and aggregate DFO NSB eDNA occurrence data across
#            survey years (2022–2024) for submission to the Ocean Biodiversity
#            Information System (OBIS). Produces three Darwin Core extension
#            files per assay: occurrence, DNADerivedData, and eMoF.
#
# Authors:   Vicky Lee, Emily Rubidge, Kristen Westfall
# Affiliation: Fisheries and Oceans Canada (DFO), Pacific Science Branch
# Contact:   Emily.Rubidge@dfo-mpo.gc.ca
#
# Data DOI:  [add upon publication]
# Version:   1.0.0
# Date:      2026-05-26
#
# Outputs (written to the project root, i.e. alongside the .Rproj file):
#   OBIS_NSB_<assay>_occurrence.csv
#   OBIS_NSB_<assay>_dnaderiveddata.csv
#   OBIS_NSB_<assay>_emof.csv
#
# Dependencies:
#   here, openxlsx, tidyr, dplyr, lubridate, purrr, worrms,
#   stringr, tibble, obistools, Hmisc
#
# Notes:
#   - This script uses the 'here' package for portable path resolution. Open
#     the project via its .Rproj file and all paths will resolve correctly
#     regardless of the user's local directory structure.
#   - WoRMS and match_taxa lookups are cached locally to avoid redundant API
#     calls across runs. Delete cache files when processing a new project or
#     if taxonomic names have been updated.
#   - When match_taxa() prompts for manual disambiguation, enter '1' to accept
#     the first (best) match, or the appropriate integer for an alternative.
#
# License:   Open Government Licence – Canada
#            https://open.canada.ca/en/open-government-licence-canada
# =============================================================================

# --- 1. Dependencies ----------------------------------------------------------

library(here)        # portable path resolution relative to .Rproj root
library(openxlsx)
library(tidyr)
library(dplyr)
library(lubridate)
library(purrr)
library(worrms)
library(stringr)
library(tibble)
library(obistools)
library(Hmisc)

# --- 2. Configuration ---------------------------------------------------------

YEARS  <- c(2022L, 2023L, 2024L)
ASSAYS <- c("12s-mifish", "16s-mcinnes", "co1-geller")

# Assay-level metadata: primer sequences and reference database identifiers.
# Primer sequences follow the 5'→3' convention.
ASSAY_META <- list(
  "12s-mifish" = list(
    targetTaxonomicAssay    = "chordata",
    target_gene             = "12S rRNA (SSU mitochondria)",
    otu_db                  = "MIDORI_UNIQ_NUC_GB256_srRNA_BLAST",
    pcr_primer_forward      = "CGTGCCAGCCACCGCGGTT",
    pcr_primer_reverse      = "CATAGTGGGGTATCTAATCCCAGTTTG",
    pcr_primer_name_forward = "248F",
    pcr_primer_name_reverse = "jgHCO2198"
  ),
  "16s-mcinnes" = list(
    targetTaxonomicAssay    = "chordata",
    target_gene             = "16S rRNA (LSU mitochondria)",
    otu_db                  = "MIDORI_UNIQ_NUC_GB256_lrRNA_BLAST",
    pcr_primer_forward      = "AGCGYAATCACTTGTCTYTTAA",
    pcr_primer_reverse      = "CRBGGTCGCCCCAACCRAA",
    pcr_primer_name_forward = "16SFish_forward",
    pcr_primer_name_reverse = "16SFish_reverse"
  ),
  "co1-geller" = list(
    targetTaxonomicAssay    = "metazoa",
    target_gene             = "COI",
    otu_db                  = "MIDORI_UNIQ_NUC_GB256_CO1_BLAST",
    pcr_primer_forward      = "GGWACWGGWTGAACWGTWTAYCCYCC",
    pcr_primer_reverse      = "TAIACYTCIGGRTGICCRAARAAYCA",
    pcr_primer_name_forward = "mlCOIint",
    pcr_primer_name_reverse = "jgHCO2198"
  )
)

# Measurement type URIs for eMoF fields.
# Fields not listed here default to the FAIRe checklist URI.
EMOF_URI_MAP <- c(
  protocol_ID     = "https://github.com/GOTeDNA-OBON",
  protocolVersion = "https://github.com/GOTeDNA-OBON",
  LClabel         = "https://github.com/GOTeDNA-OBON",
  ownerContact    = "https://github.com/GOTeDNA-OBON",
  samplingStation = "https://github.com/GOTeDNA-OBON",
  filtrationType  = "https://github.com/GOTeDNA-OBON",
  totalDNAconc    = "https://github.com/GOTeDNA-OBON",
  unitsDNAconc    = "https://github.com/GOTeDNA-OBON",
  dateFiltration  = "https://github.com/GOTeDNA-OBON",
  timeFiltration  = "https://github.com/GOTeDNA-OBON",
  volumeFiltered  = "https://github.com/GOTeDNA-OBON",
  depthWaterTemp  = "https://github.com/GOTeDNA-OBON",
  occurrenceID    = "https://manual.obis.org/darwin_core.html"
)

# Cache file paths — stored at project root, excluded from version control.
# Add both filenames to .gitignore if using Git.
WORMS_CACHE_PATH <- here("worms_cache.rds")
TAXA_CACHE_PATH  <- here("taxa_match_cache.rds")

# --- 3. File-reading helpers --------------------------------------------------

#' Read and reshape a taxaOccurrence spreadsheet to long format.
#'
#' @param path Path to the taxaOccurrence .xlsx file.
#' @return A tibble with columns: samp_name, seq_id, scientificName, genus,
#'   organismQuantity.
clean_taxaOccu <- function(path) {
  df <- read.xlsx(path, startRow = 3)
  
  # Remove higher taxonomy columns already present in sample metadata
  df <- df[, !colnames(df) %in% c("domain", "phylum", "class", "order", "family")]
  
  # Fill missing scientificName from genus (genus-level identifications)
  df$scientificName <- ifelse(
    is.na(df$scientificName),
    df$genus,
    df$scientificName
  )
  
  df_long <- pivot_longer(
    df,
    cols      = -c(seq_id, scientificName, genus),
    names_to  = "samp_name",
    values_to = "organismQuantity"
  )
  
  dplyr::select(df_long, samp_name, everything())
}

#' Extract sequence ID and representative DNA sequence from a taxaRaw file.
#'
#' @param path Path to the taxaRaw .xlsx file.
#' @return A tibble with columns: seq_id, dna_sequence.
clean_taxaRaw <- function(path) {
  df <- read.xlsx(path, startRow = 3)
  df[, c("seq_id", "dna_sequence")]
}

#' Read sample metadata, selecting and coercing relevant columns.
#'
#' Column types are explicitly coerced to ensure consistent binding across
#' survey years where input encoding may differ (e.g. numeric stored as
#' character in one year's spreadsheet).
#'
#' @param path Path to the master metadata .xlsx file.
#' @return A tibble of sample-level metadata.
clean_metadata <- function(path) {
  df <- read.xlsx(
    path,
    startRow    = 3,
    sheet       = "sampleMetadata",
    detectDates = TRUE
  )
  
  df <- df[, c(
    "samp_name", "region", "samp_category", "station_name", "eventDate",
    "decimalLongitude", "decimalLatitude", "filter_material", "filter_name",
    "prepped_samp_store_sol", "prepped_samp_store_temp", "nucl_acid_ext_kit",
    "nucl_acid_ext_modify", "env_broad_scale", "env_local_scale", "env_medium",
    "samp_collect_device", "filter_passive_active_0_1", "samp_size",
    "samp_size_unit", "minimumDepthInMeters", "maximumDepthInMeters",
    "samp_mat_process", "geo_loc_name"
  )]
  
  df <- df %>%
    mutate(
      prepped_samp_store_sol    = as.character(prepped_samp_store_sol),
      prepped_samp_store_temp   = as.numeric(prepped_samp_store_temp),
      filter_passive_active_0_1 = as.numeric(filter_passive_active_0_1),
      samp_size                 = as.numeric(samp_size),
      minimumDepthInMeters      = as.numeric(minimumDepthInMeters),
      maximumDepthInMeters      = as.numeric(maximumDepthInMeters)
    )
  
  return(df)
}

#' Read sequencing run metadata for a given sample set.
#'
#' @param path Path to the master metadata .xlsx file.
#' @return A tibble with columns: samp_name, assay_name, seq_run_id.
clean_run_meta <- function(path) {
  df <- read.xlsx(path, startRow = 3, sheet = "experimentRunMetadata")
  df[, c("samp_name", "assay_name", "seq_run_id")]
}

# --- 4. WoRMS taxon-ID cache --------------------------------------------------

#' Load the WoRMS ID cache from disk, or return an empty list.
load_worms_cache <- function() {
  if (file.exists(WORMS_CACHE_PATH)) readRDS(WORMS_CACHE_PATH) else list()
}

#' Persist the WoRMS ID cache to disk.
save_worms_cache <- function(cache) {
  saveRDS(cache, WORMS_CACHE_PATH)
}

#' Retrieve the WoRMS AphiaID for a scientific name, using the local cache.
#'
#' Returns NA_integer_ when the name is not found or the API call fails.
#'
#' @param name    Scientific name string.
#' @param cache   Named list (name → AphiaID) loaded via load_worms_cache().
#' @return Integer AphiaID, or NA_integer_.
lookup_taxon_id <- function(name, cache) {
  if (name %in% names(cache)) return(cache[[name]])
  
  id <- tryCatch(
    wm_name2id(name, marine_only = FALSE),
    error = function(e) NA_integer_
  )
  
  if (is.null(id) || length(id) == 0) NA_integer_ else as.integer(id[[1L]])
}

#' Retrieve full WoRMS classification and LSID for an AphiaID.
#'
#' @param id Integer AphiaID.
#' @return Single-row tibble with taxonID, scientificNameID, kingdom, phylum,
#'   class, order, family, genus.
lookup_taxonomy <- function(id) {
  classif <- tryCatch(wm_classification(id), error = function(e) NULL)
  rec     <- tryCatch(wm_record(id),         error = function(e) NULL)
  
  # Helper: extract a rank value from the classification data frame
  rank_val <- function(rank_name) {
    if (is.null(classif) ||
        !all(c("rank", "scientificname") %in% names(classif))) {
      return(NA_character_)
    }
    v <- classif$scientificname[classif$rank == rank_name]
    if (length(v) == 0L) NA_character_ else as.character(v[[1L]])
  }
  
  tibble(
    taxonID          = id,
    scientificNameID = if (!is.null(rec) && "lsid" %in% names(rec))
      rec$lsid else NA_character_,
    kingdom = rank_val("Kingdom"),
    phylum  = rank_val("Phylum"),
    class   = rank_val("Class"),
    order   = rank_val("Order"),
    family  = rank_val("Family"),
    genus   = rank_val("Genus")
  )
}

# --- 5. match_taxa cache ------------------------------------------------------

#' Load the obistools::match_taxa() cache from disk, or return an empty list.
load_taxa_cache <- function() {
  if (file.exists(TAXA_CACHE_PATH)) readRDS(TAXA_CACHE_PATH) else list()
}

#' Persist the obistools::match_taxa() cache to disk.
save_taxa_cache <- function(cache) {
  saveRDS(cache, TAXA_CACHE_PATH)
}

# --- 6. Core processing function: one year × one assay -----------------------

#' Process a single year–assay combination.
#'
#' Reads occurrence, sequence, sample metadata, and run metadata files;
#' joins them; constructs Darwin Core identifiers; and performs WoRMS
#' taxon-ID lookups. Returns a list with element \code{occurrence} (a
#' tibble), or NULL if any required input file is absent.
#'
#' Where a scientificName maps to multiple seq_ids within the same sample
#' (i.e. duplicate occurrenceIDs), only the row with the highest
#' organismQuantity is retained. This is conservative: it selects the most
#' abundant representative sequence while preserving one unambiguous record
#' per occurrence.
#'
#' @param year  Integer survey year (e.g. 2022L).
#' @param assay Character assay identifier (e.g. "12s-mifish").
#' @return Named list with element \code{occurrence}, or NULL.
process_one <- function(year, assay) {
  folder <- here("dfo-msea-nsb", paste0("dfo-msea-nsb-", year))
  prefix <- paste0("dfo-msea-nsb-", year)
  
  occ_path  <- file.path(folder, paste0("taxaOccurrence_", prefix, "_", assay, ".xlsx"))
  raw_path  <- file.path(folder, paste0("taxaRaw_",        prefix, "_", assay, ".xlsx"))
  meta_path <- file.path(folder, paste0(prefix, ".xlsx"))
  
  if (!file.exists(occ_path) || !file.exists(raw_path) || !file.exists(meta_path)) {
    message(sprintf("  Skipping %d / %s — one or more input files not found.", year, assay))
    return(NULL)
  }
  
  message(sprintf("  Processing %d / %s ...", year, assay))
  
  # --- Join occurrence, sequence, and metadata tables
  OBIS_prep <- clean_taxaOccu(occ_path)
  OBIS_prep <- left_join(OBIS_prep, clean_taxaRaw(raw_path),   by = "seq_id")
  OBIS_prep <- left_join(OBIS_prep, clean_metadata(meta_path), by = "samp_name")
  OBIS_prep <- left_join(
    OBIS_prep,
    clean_run_meta(meta_path),
    by           = "samp_name",
    relationship = "many-to-many"
  )
  
  # --- Construct Darwin Core identifiers
  # eventID: uniquely identifies a sampling event (one filter/sample)
  # occurrenceID: uniquely identifies one taxon in one sample for one assay.
  # A row suffix is appended to guard against the rare case where the same
  # scientificName appears under multiple seq_ids in the same sample.
  OBIS_prep <- OBIS_prep %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
    mutate(
      eventID      = paste0("DFO_MSEA_NSB-", samp_name),
      occurrenceID = paste0(eventID, "-", assay_name, "-", scientificName),
      row_suffix   = row_number(),
      occurrenceID = paste0(occurrenceID, "-", row_suffix)
    ) %>%
    select(-row_suffix)
  
  am <- ASSAY_META[[assay]]
  
  # --- Build occurrence table with all Darwin Core + eDNA fields
  occurrence <- OBIS_prep %>%
    mutate(
      materialSampleID        = as.character(samp_name),
      protocol_ID             = 1L,
      protocolVersion         = 1L,
      project_contact         = "Emily.Rubidge@dfo-mpo.gc.ca",
      bibliographicCitation   = "Not published",
      basisOfRecord           = "MaterialSample",
      locationID              = station_name,
      technical_rep_id        = "1",
      recordedBy              = "Emily Rubidge, Kristen Westfall",
      nameAccordingTo         = "WoRMS",
      country                 = "Canada",
      datasetID               = "DFO-MSEA-NSB",
      size_frac               = "0.22",
      samp_store_sol          = prepped_samp_store_sol,
      samp_store_temp         = prepped_samp_store_temp,
      organismQuantityType    = "DNA sequence reads",
      occurrenceStatus        = case_when(
        organismQuantity >  0 ~ "present",
        organismQuantity == 0 ~ "absent",
        TRUE                  ~ NA_character_
      ),
      language                = "en",
      month                   = month(eventDate),
      year                    = year(eventDate),
      project_name            = "Northern Shelf Bioregion eDNA",
      targetTaxonomicAssay    = am$targetTaxonomicAssay,
      assay_type              = "metabarcoding",
      target_gene             = am$target_gene,
      target_subfragment      = assay_name,
      seq_meth                = "ILLUMINA NovaSeq6000|NovaSeqX Plus",
      seq_kit                 = "SP (2 x 150 bp)|25B Reagent Kit (2 × 150 bp)",
      otu_seq_comp_appr       = "blastn v2.12",
      otu_db                  = am$otu_db,
      pcr_primer_forward      = am$pcr_primer_forward,
      pcr_primer_reverse      = am$pcr_primer_reverse,
      pcr_primer_name_forward = am$pcr_primer_name_forward,
      pcr_primer_name_reverse = am$pcr_primer_name_reverse,
      checkls_ver             = "1.0.2",
      sampleSizeUnit          = samp_size_unit,
      samp_collec_device      = samp_collect_device,
      project_id              = "DFO-MSEA-NSB",
      pcr_0_1                 = "1"
    ) %>%
    filter(
      !is.na(minimumDepthInMeters),
      !is.na(maximumDepthInMeters),
      !is.na(target_gene),
      !is.na(scientificName)
    )
  
  # --- WoRMS AphiaID lookup (with local cache)
  cache      <- load_worms_cache()
  species    <- unique(na.omit(occurrence$scientificName))
  new_species <- setdiff(species, names(cache))
  
  if (length(new_species) > 0L) {
    message(sprintf(
      "    Looking up %d new taxon name(s) in WoRMS ...",
      length(new_species)
    ))
    new_ids <- setNames(
      map_int(new_species, ~ lookup_taxon_id(.x, cache)),
      new_species
    )
    cache <- c(cache, as.list(new_ids))
    save_worms_cache(cache)
  }
  
  species_lookup <- tibble(
    scientificName = species,
    taxonID        = map_int(species, function(nm) {
      v <- cache[[nm]]
      if (is.null(v)) NA_integer_ else as.integer(v)
    })
  )
  
  occurrence <- occurrence %>%
    select(-any_of(c(
      "taxonID", "scientificNameID", "kingdom", "phylum",
      "class", "order", "family", "genus", "aphiaID"
    ))) %>%
    left_join(species_lookup, by = "scientificName") %>%
    mutate(
      worms_match = !is.na(taxonID),
      aphiaID     = taxonID
    )
  
  # --- Full classification lookup for matched taxa
  taxonomy_lookup <- map_dfr(
    discard(unique(occurrence$taxonID), is.na),
    lookup_taxonomy
  )
  
  occurrence <- occurrence %>%
    select(-any_of(c(
      "kingdom", "phylum", "class", "order", "family", "genus", "scientificNameID"
    ))) %>%
    left_join(taxonomy_lookup, by = "taxonID") %>%
    filter(!is.na(organismQuantity))
  
  # Deduplicate: retain the highest-read-count representative per occurrenceID
  occurrence <- occurrence %>%
    group_by(occurrenceID) %>%
    slice_max(organismQuantity, n = 1L, with_ties = FALSE) %>%
    ungroup()
  
  list(occurrence = occurrence)
}

# --- 7. Run across all year × assay combinations -----------------------------

message("\n=== DFO NSB eDNA OBIS Pipeline — starting run ===\n")

results <- list()

for (assay in ASSAYS) {
  message(sprintf("Assay: %s", assay))
  results[[assay]] <- list()
  
  for (yr in YEARS) {
    res <- process_one(yr, assay)
    if (!is.null(res)) {
      results[[assay]][[as.character(yr)]] <- res$occurrence
    }
  }
}

# --- 8. Aggregate by assay, run match_taxa, write outputs --------------------

for (assay in ASSAYS) {
  yr_list <- results[[assay]]
  
  if (length(yr_list) == 0L) {
    message(sprintf("\nNo data found for assay: %s — skipping output.", assay))
    next
  }
  
  message(sprintf("\nAggregating assay: %s", assay))
  
  combined <- bind_rows(yr_list)
  
  # --- obistools::match_taxa() with local cache
  taxa_cache   <- load_taxa_cache()
  unique_names <- unique(na.omit(combined$scientificName))
  new_names    <- setdiff(unique_names, names(taxa_cache))
  
  if (length(new_names) > 0L) {
    message(sprintf(
      "  Running match_taxa() for %d new name(s) — please respond to prompts below.",
      length(new_names)
    ))
    worms_matches <- match_taxa(new_names)
    
    for (i in seq_along(new_names)) {
      taxa_cache[[new_names[[i]]]] <- worms_matches[i, , drop = FALSE]
    }
    save_taxa_cache(taxa_cache)
    
  } else {
    message("  All names already cached — skipping match_taxa().")
  }
  
  # Rebuild per-species lookup from cache (one row per name, first match kept)
  worms <- bind_rows(taxa_cache[unique_names])
  worms$scientificName <- unique_names
  worms <- worms %>%
    group_by(scientificName) %>%
    slice(1L) %>%
    ungroup()
  
  combined <- combined %>%
    select(-any_of("scientificNameID")) %>%
    left_join(worms, by = "scientificName")
  
  # Harmonise column naming after join (match_taxa may return scientificNameID.y)
  if ("scientificNameID.y" %in% colnames(combined)) {
    combined <- combined %>%
      rename(scientificNameID = scientificNameID.y) %>%
      select(-any_of("scientificNameID.x"))
  }
  
  # Final deduplication after match_taxa merge
  combined <- combined %>%
    group_by(occurrenceID) %>%
    slice_max(organismQuantity, n = 1L, with_ties = FALSE) %>%
    ungroup()
  
  # --- Uniqueness check (required for OBIS submission)
  n_total  <- nrow(combined)
  n_unique <- length(unique(combined$occurrenceID))
  
  if (n_total != n_unique) {
    warning(sprintf(
      "Duplicate occurrenceIDs detected in assay '%s': %d duplicate(s). Check join logic.",
      assay, n_total - n_unique
    ))
  } else {
    message(sprintf("  occurrenceID uniqueness check passed (%d records).", n_total))
  }
  
  # --- Select columns for each output file
  
  occurrence_out <- combined %>%
    select(
      occurrenceID, bibliographicCitation, materialSampleID, eventDate,
      decimalLatitude, decimalLongitude, scientificName, organismQuantity,
      organismQuantityType, sampleSizeUnit, eventID,
      basisOfRecord, locationID, recordedBy, country, datasetID,
      occurrenceStatus, minimumDepthInMeters, maximumDepthInMeters,
      language, month, year, taxonID, scientificNameID,
      kingdom, phylum, class, order, family, genus
    )
  
  dna_out <- combined %>%
    select(
      target_gene, target_subfragment, samp_name, samp_mat_process,
      size_frac, samp_size, project_name, env_broad_scale,
      env_local_scale, env_medium, seq_meth, otu_seq_comp_appr,
      otu_db, pcr_primer_forward, pcr_primer_reverse,
      pcr_primer_name_forward, pcr_primer_name_reverse,
      dna_sequence, occurrenceID, samp_collec_device
    )
  
  emof_out <- combined %>%
    select(
      seq_id, samp_category, filter_passive_active_0_1, checkls_ver,
      seq_kit, assay_name, assay_type, targetTaxonomicAssay,
      samp_collect_device, samp_size_unit, filter_material, filter_name,
      samp_store_sol, samp_store_temp, nucl_acid_ext_kit,
      nucl_acid_ext_modify, geo_loc_name, technical_rep_id,
      project_contact, protocol_ID, protocolVersion,
      occurrenceID, seq_run_id, project_id, pcr_0_1
    ) %>%
    distinct() %>%
    mutate(across(-occurrenceID, as.character)) %>%
    pivot_longer(
      cols      = -occurrenceID,
      names_to  = "measurementType",
      values_to = "measurementValue"
    ) %>%
    mutate(
      measurementTypeID  = dplyr::recode(
        measurementType,
        !!!EMOF_URI_MAP,
        .default = "https://github.com/FAIR-eDNA/FAIRe_checklist"
      ),
      measurementUnit    = NA_character_,
      measurementValueID = NA_character_,
      measurementUnitID  = NA_character_,
      measurementRemarks = NA_character_
    )
  
  # --- Write output CSVs
  slug <- gsub("-", "_", assay)
  
  write.csv(
    occurrence_out,
    here(paste0("OBIS_NSB_", slug, "_occurrence.csv")),
    row.names = FALSE
  )
  write.csv(
    dna_out,
    here(paste0("OBIS_NSB_", slug, "_dnaderiveddata.csv")),
    row.names = FALSE
  )
  write.csv(
    emof_out,
    here(paste0("OBIS_NSB_", slug, "_emof.csv")),
    row.names = FALSE
  )
  
  message(sprintf("  Outputs written for assay: %s", assay))
}

# --- 9. Session information ---------------------------------------------------
# Printed at the end for reproducibility reporting in supplementary materials.

message("\n=== Run complete ===\n")
message("Session information (record in supplementary materials):\n")
print(sessionInfo())
