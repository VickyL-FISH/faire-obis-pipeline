# install.packages("dplyr")
# install.packages("purrr")
# install.packages("worrms")
# install.packages("openxlsx")
# install.packages("lubridate")
# install.packages("tidyr")
# install.packages("obistools")

# PURPOSE: Standardize and aggregate DFO NSB eDNA data across years for OBIS submission
# Outputs: occurrence, dnaderiveddata, and eMoF CSVs per assay (2022-2024 combined)
# NOTE: Zero-read rows (organismQuantity == 0) are intentionally retained so that
#       GOTeDNA can calculate detection rates. These rows have occurrenceStatus = "absent".

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


##### Config — edit here as needed #####

# Set your working directory to the root of this repository, then set the paths below.
# All paths are relative to the repository root.

input_dir  <- "example_data/input"   # folder containing per-year FAIRe input files
output_dir <- "example_data/output"  # folder where OBIS CSVs will be written
cache_dir  <- "cache"                # folder for WoRMS and taxa match caches

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir,  showWarnings = FALSE, recursive = TRUE)

years  <- c(2022, 2023, 2024)
assays <- c("12s-mifish", "16s-mcinnes", "co1-geller")

assay_meta <- list(
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

url_map <- c(
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


##### Cleaning functions #####

clean_taxaOccu <- function(path) {
  df <- read.xlsx(path, startRow = 3)
  df <- df[, !colnames(df) %in% c("domain", "phylum", "class", "order", "family")]
  df$scientificName <- ifelse(is.na(df$scientificName), paste0(df$genus), df$scientificName)
  
  # pivot_longer keeps ALL sample x species combos, including zero-read rows.
  # These zeros are required for GOTeDNA detection rate calculations and will
  # appear as occurrenceStatus = "absent" in the final output.
  df_long <- pivot_longer(
    df,
    cols      = -c(seq_id, scientificName, genus),
    names_to  = "samp_name",
    values_to = "organismQuantity"
  )
  
  # Replace NA read counts with explicit 0 so absent rows are unambiguous
  df_long <- df_long %>%
    mutate(organismQuantity = replace_na(organismQuantity, 0))
  
  df_long %>% select(samp_name, everything())
}

clean_taxaRaw <- function(path) {
  df <- read.xlsx(path, startRow = 3)
  df[, c("seq_id", "dna_sequence")]
}

clean_metadata <- function(path) {
  df <- read.xlsx(path, startRow = 3, sheet = "sampleMetadata", detectDates = TRUE)
  df <- df[, c("samp_name", "region", "samp_category", "station_name", "eventDate",
               "decimalLongitude", "decimalLatitude", "filter_material", "filter_name",
               "prepped_samp_store_sol", "prepped_samp_store_temp", "nucl_acid_ext_kit",
               "nucl_acid_ext_modify", "env_broad_scale", "env_local_scale", "env_medium",
               "samp_collect_device", "filter_passive_active_0_1", "samp_size", "samp_size_unit",
               "minimumDepthInMeters", "maximumDepthInMeters", "samp_mat_process", "geo_loc_name")]
  
  # Coerce columns that are prone to type inconsistencies across years
  df <- df %>% mutate(
    prepped_samp_store_sol    = as.character(prepped_samp_store_sol),
    prepped_samp_store_temp   = as.numeric(prepped_samp_store_temp),
    filter_passive_active_0_1 = as.numeric(filter_passive_active_0_1),
    samp_size                 = as.numeric(samp_size),
    minimumDepthInMeters      = as.numeric(minimumDepthInMeters),
    maximumDepthInMeters      = as.numeric(maximumDepthInMeters)
  )
  
  return(df)
}

clean_run_meta <- function(path) {
  df <- read.xlsx(path, startRow = 3, sheet = "experimentRunMetadata")
  df[, c("samp_name", "assay_name", "seq_run_id")]
}


##### WoRMS cache #####

worms_cache_path <- file.path(cache_dir, "worms_cache.rds")

load_worms_cache <- function() {
  if (file.exists(worms_cache_path)) readRDS(worms_cache_path) else list()
}

save_worms_cache <- function(cache) {
  saveRDS(cache, worms_cache_path)
}

lookup_taxon_id <- function(name, cache) {
  if (name %in% names(cache)) return(cache[[name]])
  id <- tryCatch(wm_name2id(name, marine_only = FALSE), error = function(e) NA)
  if (is.null(id) || length(id) == 0) NA_integer_ else as.integer(id[1])
}

lookup_taxonomy <- function(id) {
  classif <- tryCatch(wm_classification(id), error = function(e) NULL)
  rec     <- tryCatch(wm_record(id),         error = function(e) NULL)
  .rank_val <- function(rank_name) {
    if (is.null(classif) || !all(c("rank", "scientificname") %in% names(classif))) return(NA_character_)
    v <- classif$scientificname[classif$rank == rank_name]
    if (length(v) == 0) NA_character_ else as.character(v[1])
  }
  tibble(
    taxonID          = id,
    scientificNameID = if (!is.null(rec) && "lsid" %in% names(rec)) rec$lsid else NA_character_,
    kingdom = .rank_val("Kingdom"), phylum = .rank_val("Phylum"),
    class   = .rank_val("Class"),   order  = .rank_val("Order"),
    family  = .rank_val("Family"),  genus  = .rank_val("Genus")
  )
}


##### match_taxa cache #####
# When running for different projects, delete the cache file at cache/taxa_match_cache.rds.
# Input y and then subsequently 1 for all species matching.

taxa_cache_path <- file.path(cache_dir, "taxa_match_cache.rds")

load_taxa_cache <- function() {
  if (file.exists(taxa_cache_path)) readRDS(taxa_cache_path) else list()
}

save_taxa_cache <- function(cache) {
  saveRDS(cache, taxa_cache_path)
}


##### Core processing function — one year + one assay #####

process_one <- function(year, assay) {
  project_id <- "dfo-msea-nsb"
  folder     <- file.path(input_dir, project_id, paste0(project_id, "-", year))
  prefix     <- paste0(project_id, "-", year)
  
  occ_path  <- file.path(folder, paste0("taxaOccurrence_", prefix, "_", assay, ".xlsx"))
  raw_path  <- file.path(folder, paste0("taxaRaw_",        prefix, "_", assay, ".xlsx"))
  meta_path <- file.path(folder, paste0(prefix, ".xlsx"))
  
  if (!file.exists(occ_path) || !file.exists(raw_path) || !file.exists(meta_path)) {
    message("Skipping ", year, " / ", assay, " — files not found")
    return(NULL)
  }
  
  message("Processing ", year, " / ", assay, " ...")
  
  OBIS_prep <- clean_taxaOccu(occ_path)
  OBIS_prep <- left_join(OBIS_prep, clean_taxaRaw(raw_path),   by = "seq_id")
  OBIS_prep <- left_join(OBIS_prep, clean_metadata(meta_path), by = "samp_name")
  OBIS_prep <- left_join(OBIS_prep, clean_run_meta(meta_path), by = "samp_name",
                         relationship = "many-to-many")
  
  # Drop rows missing coordinates — we cannot place these spatially
  OBIS_prep <- OBIS_prep %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude))
  
  # Build occurrenceID from eventID + assay_name + scientificName + seq_id so that
  # every species x sample combination (including zeros) gets a stable, unique ID.
  # seq_id scopes the ID to the exact OTU/ASV, preventing collisions when the same
  # species name appears under multiple seq_ids in the same sample.
  OBIS_prep <- OBIS_prep %>%
    mutate(
      eventID      = paste0("DFO_MSEA_NSB-", samp_name),
      occurrenceID = paste0(eventID, "-", assay_name, "-", scientificName, "-", seq_id)
    )
  
  am <- assay_meta[[assay]]
  
  occurrence <- OBIS_prep %>%
    mutate(
      materialSampleID         = as.character(samp_name),
      protocol_ID              = 1,
      protocolVersion          = 1,
      project_contact          = "Emily.Rubidge@dfo-mpo.gc.ca",
      bibliographicCitation    = "Not published",
      basisOfRecord            = "MaterialSample",
      locationID               = station_name,
      technical_rep_id         = "1",
      recordedBy               = "Emily Rubidge, Kristen Westfall",
      nameAccordingTo          = "WoRMS",
      country                  = "Canada",
      datasetID                = "DFO-MSEA-NSB",
      size_frac                = "0.22",
      samp_store_sol           = prepped_samp_store_sol,
      samp_store_temp          = prepped_samp_store_temp,
      organismQuantityType     = "DNA sequence reads",
      # Zero-read rows get occurrenceStatus = "absent" — required for GOTeDNA detection rates
      occurrenceStatus         = case_when(
        organismQuantity > 0  ~ "present",
        organismQuantity == 0 ~ "absent",
        TRUE                  ~ NA_character_
      ),
      language                 = "en",
      month                    = month(eventDate),
      year                     = year(eventDate),
      project_name             = "Northern Shelf Bioregion eDNA",
      targetTaxonomicAssay     = am$targetTaxonomicAssay,
      assay_type               = "metabarcoding",
      target_gene              = am$target_gene,
      target_subfragment       = assay_name,
      seq_meth                 = "ILLUMINA NovaSeq6000|NovaSeqX Plus",
      seq_kit                  = "SP (2 x 150 bp)|25B Reagent Kit (2 x 150 bp)",
      otu_seq_comp_appr        = "blastn v2.12",
      otu_db                   = am$otu_db,
      pcr_primer_forward       = am$pcr_primer_forward,
      pcr_primer_reverse       = am$pcr_primer_reverse,
      pcr_primer_name_forward  = am$pcr_primer_name_forward,
      pcr_primer_name_reverse  = am$pcr_primer_name_reverse,
      checkls_ver              = "1.0.2",
      sampleSizeUnit           = samp_size_unit,
      samp_collec_device       = samp_collect_device,
      project_id               = "DFO-MSEA-NSB",
      pcr_0_1                  = "1"
    ) %>%
    filter(
      !is.na(minimumDepthInMeters), !is.na(maximumDepthInMeters),
      !is.na(target_gene), !is.na(scientificName)
      # NOTE: organismQuantity == 0 rows are intentionally NOT filtered out here
    )
  
  # WoRMS taxon ID lookup with cache
  cache   <- load_worms_cache()
  species <- unique(na.omit(occurrence$scientificName))
  new_species <- setdiff(species, names(cache))
  
  if (length(new_species) > 0) {
    message("  Looking up ", length(new_species), " new species in WoRMS...")
    new_ids <- setNames(
      map_int(new_species, ~ lookup_taxon_id(.x, cache)),
      new_species
    )
    cache <- c(cache, as.list(new_ids))
    save_worms_cache(cache)
  }
  
  species_lookup <- tibble(
    scientificName = species,
    taxonID        = map_int(species, ~ { v <- cache[[.x]]; if (is.null(v)) NA_integer_ else as.integer(v) })
  )
  
  occurrence <- occurrence %>%
    select(-any_of(c("taxonID", "scientificNameID", "kingdom", "phylum",
                     "class", "order", "family", "genus", "aphiaID"))) %>%
    left_join(species_lookup, by = "scientificName") %>%
    mutate(worms_match = !is.na(taxonID), aphiaID = taxonID)
  
  taxonomy_lookup <- map_dfr(
    occurrence$taxonID %>% unique() %>% discard(is.na),
    lookup_taxonomy
  )
  
  occurrence <- occurrence %>%
    select(-any_of(c("kingdom", "phylum", "class", "order",
                     "family", "genus", "scientificNameID"))) %>%
    left_join(taxonomy_lookup, by = "taxonID")
  
  # Deduplicate only genuine many-to-many join duplicates (same occurrenceID from
  # run metadata fan-out). Among true duplicates keep the highest organismQuantity.
  # This does NOT drop zero-read rows that are unique combinations of species x sample.
  occurrence <- occurrence %>%
    group_by(occurrenceID) %>%
    slice_max(organismQuantity, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  list(occurrence = occurrence)
}


##### Run across all year + assay combinations #####

results <- list()

for (assay in assays) {
  results[[assay]] <- list()
  for (yr in years) {
    res <- process_one(yr, assay)
    if (!is.null(res)) results[[assay]][[as.character(yr)]] <- res$occurrence
  }
}


##### Aggregate by assay, run match_taxa, write outputs #####

for (assay in assays) {
  yr_list <- results[[assay]]
  if (length(yr_list) == 0) next
  
  combined <- bind_rows(yr_list)
  
  # --- match_taxa with caching ---
  taxa_cache   <- load_taxa_cache()
  unique_names <- unique(na.omit(combined$scientificName))
  new_names    <- setdiff(unique_names, names(taxa_cache))
  
  if (length(new_names) > 0) {
    message("Running match_taxa() for ", length(new_names), " new species in assay: ", assay)
    message("Please respond to the prompts below...\n")
    
    worms_matches <- match_taxa(new_names)
    
    for (i in seq_along(new_names)) {
      taxa_cache[[new_names[i]]] <- worms_matches[i, ]
    }
    save_taxa_cache(taxa_cache)
    
  } else {
    message("All species for ", assay, " already matched — using cache")
  }
  
  # Rebuild match table from cache, one row per species
  worms <- bind_rows(taxa_cache[unique_names])
  worms$scientificName <- unique_names
  
  # Deduplicate — keep first match per species
  worms <- worms %>%
    group_by(scientificName) %>%
    slice(1) %>%
    ungroup()
  
  # Drop any existing scientificNameID before merging
  combined <- combined %>% select(-any_of("scientificNameID"))
  combined <- left_join(combined, worms, by = "scientificName")
  
  if ("scientificNameID.y" %in% colnames(combined)) {
    colnames(combined)[colnames(combined) == "scientificNameID.y"] <- "scientificNameID"
    combined <- combined %>% select(-any_of("scientificNameID.x"))
  }
  
  # Final deduplication after match_taxa merge — keep highest organismQuantity per ID.
  # Zero-read rows that are unique (no duplicate occurrenceID) are unaffected.
  combined <- combined %>%
    group_by(occurrenceID) %>%
    slice_max(organismQuantity, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Uniqueness check
  n_total  <- length(combined$occurrenceID)
  n_unique <- length(unique(combined$occurrenceID))
  if (n_total != n_unique) {
    warning("Duplicate occurrenceIDs still detected in assay: ", assay,
            " (", n_total - n_unique, " duplicates). Check your merge.")
  } else {
    message("occurrenceID uniqueness check passed for: ", assay)
  }
  
  # Report absent vs present row counts for transparency
  n_absent  <- sum(combined$occurrenceStatus == "absent",  na.rm = TRUE)
  n_present <- sum(combined$occurrenceStatus == "present", na.rm = TRUE)
  message("  Rows — present: ", n_present, "  absent (zero reads): ", n_absent)
  
  # --- Select fields and write outputs ---
  
  cleaned_occ <- combined %>%
    select(
      occurrenceID, bibliographicCitation, materialSampleID, eventDate,
      decimalLatitude, decimalLongitude, scientificName, organismQuantity,
      organismQuantityType, sampleSizeUnit, eventID,
      basisOfRecord, locationID, recordedBy, country, datasetID,
      occurrenceStatus, minimumDepthInMeters, maximumDepthInMeters,
      language, month, year, taxonID, scientificNameID,
      kingdom, phylum, class, order, family, genus
    )
  
  cleaned_dna <- combined %>%
    select(
      target_gene, target_subfragment, samp_name, samp_mat_process, size_frac,
      samp_size, project_name, env_broad_scale, env_local_scale, env_medium,
      seq_meth, otu_seq_comp_appr, otu_db, pcr_primer_forward, pcr_primer_reverse,
      pcr_primer_name_forward, pcr_primer_name_reverse, dna_sequence,
      occurrenceID, samp_collec_device
    )
  
  emof <- combined %>%
    select(
      seq_id, samp_category, filter_passive_active_0_1, checkls_ver, seq_kit,
      assay_name, assay_type, targetTaxonomicAssay, samp_collect_device,
      samp_size_unit, filter_material, filter_name, samp_store_sol, samp_store_temp,
      nucl_acid_ext_kit, nucl_acid_ext_modify, geo_loc_name, technical_rep_id,
      project_contact, protocol_ID, protocolVersion,
      occurrenceID, seq_run_id, project_id, pcr_0_1
    ) %>%
    distinct() %>%
    mutate(across(-occurrenceID, as.character)) %>%
    pivot_longer(cols = -occurrenceID, names_to = "measurementType", values_to = "measurementValue") %>%
    mutate(
      measurementTypeID  = dplyr::recode(measurementType, !!!url_map,
                                         .default = "https://github.com/FAIR-eDNA/FAIRe_checklist"),
      measurementUnit    = NA_character_,
      measurementValueID = NA_character_,
      measurementUnitID  = NA_character_,
      measurementRemarks = NA_character_
    )
  
  slug <- gsub("-", "_", assay)
  write.csv(cleaned_occ, file.path(output_dir, paste0("OBIS_NSB_", slug, "_occurrence.csv")),     row.names = FALSE)
  write.csv(cleaned_dna, file.path(output_dir, paste0("OBIS_NSB_", slug, "_dnaderiveddata.csv")), row.names = FALSE)
  write.csv(emof,        file.path(output_dir, paste0("OBIS_NSB_", slug, "_emof.csv")),           row.names = FALSE)
  
  message("Written outputs for assay: ", assay, "\n")
}

message("All done.")