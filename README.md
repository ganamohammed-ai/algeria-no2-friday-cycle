# Algerian NO₂ Friday cycle and short-event attribution

Publication code for **The Friday NO₂ decline across Algerian cities in weather-normalised TROPOMI observations: Replication and limits of short-event attribution**.

Authors: Mohamed GANA, Chaima TOUABA, Thierry SAINT-GERAND. Correspondence: gana.mohammed@umc.edu.dz. Mohamed GANA: https://orcid.org/0000-0003-0775-6954 . Affiliations are in CITATION.cff and .zenodo.json. Thierry's Caen affiliation was confirmed by the corresponding author during package preparation.

## Public data and software

- Dataset: https://doi.org/10.5281/zenodo.24445977
- Source repository: https://github.com/ganamohammed-ai/algeria-no2-friday-cycle
- Fixed software release and installable QGIS plugin: https://github.com/ganamohammed-ai/algeria-no2-friday-cycle/releases/tag/v1.0.0

The Zenodo record archives the processed data, selected upstream inputs and the code/plugin ZIPs. Its DOI identifies this combined supporting dataset; it is not the manuscript DOI or a separate software-release DOI.

## Scientific scope

Six study sites, six background sites and seven further cities; original record 7 August 2019–31 December 2025, with replication through 14 November 2026. Weekly-cycle discovery was exploratory. Subsequent predictions and deviations are documented in the original staged plans, preserved verbatim in scripts/analysis/. They are not represented as a public preregistration. The failed Annaba ranking prediction remains part of the record. The 2020 ranges and injection calibration remain exploratory.

## Installation

Python 3.12 was tested in a new isolated Windows environment. requirements-tested.txt records the exact package versions of this replay environment, **not** the undocumented original execution environment.

```sh
python -m venv .venv
# Windows
.venv\Scripts\python -m pip install -r requirements-tested.txt
# Linux/macOS: use .venv/bin/python instead; other OSes were not tested.
```

requirements.txt lists the direct dependencies for a fresh environment. The pinned tested file is preferred for the documented Windows replay. Direct raw-source downloaders additionally require the user's own Copernicus/CDS account and accepted product terms. Never commit credentials.

## Data layout and replay

The delivered desktop package has sibling GitHub/ and Zenodo/ directories. The default runner reads ../Zenodo/processed/. For an independent clone, download **01_processed_data_and_results.zip** from the dataset record, extract it into a directory, and pass that directory as --data-root. The public dataset and fixed software-release links are listed above.

```sh
python reproduce.py --mode quick --data-root ../Zenodo/processed --out reproduced_results
python reproduce.py --mode extended --data-root ../Zenodo/processed --out reproduced_results_extended
python reproduce.py --mode full --data-root ../Zenodo/processed --out reproduced_results_full
```

* quick: original Friday model, 2026 replication, seven cities, separate-period and pooled covariance checks, Table 8 export and four Friday subperiods.
* extended: quick plus containment OLS/reference-window/trend tables, six seasonal windows, wind/BLH decomposition, processor-step sensitivity, Level-2 elasticity comparison, ventilation models, dust/outlier/seasonal sensitivity and all five figures. Figures 4 and S1 use archived estimates/spatial summaries. Updated Figure 2 is generated explicitly; Figure 1 labels use computed slopes rather than fixed labels.
* full: extended plus the expensive original ML counterfactual, injections and calibration. These expensive stages were **not** fully rerun for this package. Their archived reference results are supplied, and their scientific limitations are unchanged. Added-control ML/counterfactual summaries require an explicit further call to lockdown_normalisation.py with the additional-control daily table as --controls; additional injection is included by full mode. Do not present a partial replay as full raw-data validation.

All replay outputs go outside the input tree. The runner checks reference CSVs and exits nonzero if a checked result differs beyond the recorded numeric tolerance. Machine-readable validation is in VALIDATION.json. The desktop Audit/ folder contains local logs and an additional test of the CLI/relocation behaviour; it is not included in the public code archive.

## Rebuilding spatial inputs

Extract upstream input archives 02–06 into one directory INPUTS. ZIP member names preserve the published data/... tree; the former local directory Nouveau dossier (4) is now data/07_replication/. Full spatial extraction was not rerun. One raw Algiers-2019 GeoTIFF extraction reproduced all 147 daily CSV records exactly.

```sh
python scripts/analysis/analyze_no2_ventilation.py --stage tables --tif-dir INPUTS/data/01_tropomi_L3_daily_geotiff --sites algiers --years 2019 --out raw_algiers_2019
python scripts/analysis/run_extra_sites.py C1 --input-root INPUTS --out-root rebuilt_extra
python scripts/analysis/run_extra_sites.py C2 --input-root INPUTS --out-root rebuilt_extra
python scripts/analysis/run_extra_sites.py D --input-root INPUTS --out-root rebuilt_extra
python scripts/analysis/l2_daily_table.py --clips INPUTS/data/02_tropomi_L2_orbit_clips INPUTS/data/02_tropomi_L2_orbit_clips_controls --table PROCESSED/outputs/no2_ventilation_daily_v2.csv PROCESSED/outputs/control/control_no2_ventilation_daily.csv --qa 0.75 --tag q75all --out rebuilt_l2
```

Run the other raw stages with their --help options and the supplied GEE scripts. Exact satellite exports are archived because later upstream reprocessing may change a re-download. One unique missing orbit (10372, 14 October 2019) is preserved in the download logs; repeated failed attempts are not separate missing orbits. The final single-version Level-2 2020 replication uses q75all, not q50 or the earlier study-only q75 files.

## Tables and figures

metadata/Manuscript_to_files.txt maps the article to input/result files. The new command for the descriptive range is:

```sh
python scripts/analysis/export_table8.py --results PROCESSED/outputs --out exported_tables
```

The runner's period-specific covariance checks write to separate subdirectories, avoiding overwrites. Full-period/city tables have different schemas from result tables; join daily data by site/date. See metadata/DATA_DICTIONARY.json and metadata/sites.csv.

## Release and citation

Code version: 1.0.0. MIT licence for code. CC BY 4.0 applies to author-owned derived data/results; original source data retain their terms (THIRD_PARTY_DATA.md). Source project files were not modified. Original plans remain unchanged; publication-only changes are listed in CHANGES_PUBLICATION.txt.

The v1.0.0 software release is fixed at commit 90ec32861ea90e0ad018bf80d1730a633b6eb621. All 146 published source/documentation files were checked byte-for-byte against the prepared package. Subsequent README and citation-link updates do not alter the scientific code in that release. The supporting Zenodo dataset also contains the original code ZIP and QGIS plugin ZIP. No separate software DOI has been minted.

## QGIS Processing tools

The QGIS/ directory contains a complete Processing plugin and Python equivalents of the supplied Earth Engine exports. See QGIS/README_QGIS.md. Source-product discrepancies discovered during conversion are documented in SOURCE_DISCREPANCIES.txt.
