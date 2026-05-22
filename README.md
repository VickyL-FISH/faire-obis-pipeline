# faire-obis-pipeline

An R pipeline to standardise and aggregate DFO Northern Shelf Bioregion eDNA metabarcoding data across years into Darwin Core format for submission to OBIS, following FAIRe checklist standards.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Authors

**Code developed by:** Vicky Lee (Research Intern, Fisheries and Oceans Canada)  
**Project leads / data authors:** Emily Rubidge, Kristen Westfall (Fisheries and Oceans Canada)

---

## Overview

This pipeline processes FAIRe-formatted eDNA survey data from the DFO Northern Shelf Bioregion (NSB) across multiple years and assays, and produces three OBIS-ready output files per assay:

- `*_occurrence.csv` — Darwin Core occurrence records
- `*_dnaderiveddata.csv` — DNA-derived data extension
- `*_emof.csv` — Extended Measurement or Fact extension

The pipeline is designed to support any number of assays. The three assays used in the DFO NSB project are included as examples: 12s-mifish, 16s-mcinnes, co1-geller. Additional assays can be added by extending the assays vector and assay_meta list in the config section of the script, no other changes are needed.

> **Note:** Zero-read rows (`organismQuantity == 0`) are intentionally retained in the output so that GOTeDNA can calculate detection rates. These rows have `occurrenceStatus = "absent"`.

---

## Requirements

- R >= 4.4
- The following R packages:

```r
install.packages(c("dplyr", "tidyr", "purrr", "lubridate",
                   "openxlsx", "worrms", "stringr", "tibble", "Hmisc"))

# obistools is installed from GitHub:
# install.packages("remotes")
# remotes::install_github("iobis/obistools")
```

---

## Input data structure

Place your FAIRe input files under `example_data/input/` following this structure:

```
example_data/input/
└── dfo-msea-nsb/
    ├── dfo-msea-nsb-2022/
    │   ├── dfo-msea-nsb-2022.xlsx
    │   ├── taxaOccurrence_dfo-msea-nsb-2022_12s-mifish.xlsx
    │   ├── taxaRaw_dfo-msea-nsb-2022_12s-mifish.xlsx
    │   └── ... (one taxaOccurrence + taxaRaw per assay)
    ├── dfo-msea-nsb-2023/
    │   └── ...
    └── dfo-msea-nsb-2024/
        └── ...
```

An example dataset with a small synthetic subset of records is provided in `example_data/` to demonstrate the expected FAIRe input format and column structure. See [Example data](#example-data) below.

---

## Usage

1. Clone this repository:

```bash
git clone https://github.com/VickyL-FISH/faire-obis-pipeline.git
cd faire-obis-pipeline
```

2. Open `faire_obis_pipeline.R` and edit the config section at the top if needed:

```r
input_dir  <- "example_data/input"   # path to your FAIRe input files
output_dir <- "example_data/output"  # where OBIS CSVs will be written
cache_dir  <- "cache"                # where WoRMS and taxa match caches are stored
```

3. Run the script in R:

```r
source("faire_obis_pipeline.R")
```

4. Output CSVs will be written to `example_data/output/`.

> **WoRMS matching:** The first run will query the WoRMS API and prompt you to confirm taxon matches interactively. Results are cached in `cache/` so subsequent runs skip this step for already-matched species.

---

## Example data

The `example_data/input/` folder contains a small synthetic dataset that mirrors the real FAIRe input format. You can use it to test the pipeline end-to-end without access to the full DFO dataset.

Key columns in the FAIRe input format:

| Column | Description |
|---|---|
| `samp_name` | Unique sample identifier |
| `scientificName` | Species name (matched to WoRMS) |
| `organismQuantity` | DNA sequence read count (0 = absent) |
| `eventDate` | Date of sample collection |
| `decimalLatitude` / `decimalLongitude` | Sample coordinates |
| `minimumDepthInMeters` / `maximumDepthInMeters` | Depth range |
| `seq_id` | OTU/ASV identifier |
| `dna_sequence` | Representative DNA sequence |

---

## Outputs

| File | Contents |
|---|---|
| `OBIS_NSB_<assay>_occurrence.csv` | Core Darwin Core occurrence fields |
| `OBIS_NSB_<assay>_dnaderiveddata.csv` | DNA-derived data extension fields |
| `OBIS_NSB_<assay>_emof.csv` | Extended Measurement or Fact extension (long format) |

---

## Citation

This pipeline will be associated with a forthcoming publication. A Zenodo DOI will be added here upon release.

```
Lee, V., Rubidge, E., Westfall, K. (2026). faire-obis-pipeline [Software].
Fisheries and Oceans Canada. https://github.com/VickyL-FISH/faire-obis-pipeline
```

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
