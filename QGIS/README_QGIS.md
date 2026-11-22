# QGIS tools for the NO2 Algeria study

## Installation

The desktop package provides **QGIS_NO2_Tools.zip**, a standard QGIS plugin ZIP. In QGIS open **Plugins → Manage and Install Plugins → Install from ZIP**, select the ZIP and enable **NO2 Algeria Reproducibility**. Open **Processing → Toolbox**. No profile installation or remote jobs were performed automatically during preparation.

There are **48 Processing tools**: one isolated-environment preparer, one study replay, a native 19-site point layer, tools for every executable analysis/downloader, and Python equivalents of all 13 Earth Engine JavaScript export scripts. Shared sites.py is configuration, not an executable analysis, and is represented by the site-layer tool. Original BAT launchers are replaced by Processing tools; historical Markdown plans remain unchanged. SCRIPT_COVERAGE.json provides the complete mapping.

## Run the study

Use **01 Reproduce the study and verify reference results**. Set:

* Scientific Python executable: the package's local `Local_QGIS_Runtime/Scripts/python.exe`, or another environment installed from requirements-qgis-tested.txt. This is a separate Python process so QGIS's own numpy/scipy packages are not overwritten.
* Processed data root: `Zenodo/processed` in the desktop package, or a directory extracted from 01_processed_data_and_results.zip. It must contain outputs/ and data/.
* Replay scope: quick, extended or full. Extended was exercised from QGIS and passed 16 numeric reference comparisons, producing five figures. Full additionally launches expensive ML/injection/calibration experiments; these were not fully rerun for this release.
* Output folder: a new independent folder outside the processed inputs.

## Prepare a Python backend on another computer

The local runtime is for this Windows computer and is not included in the public ZIPs. Use **00 Prepare isolated scientific Python environment** with a standalone Python 3.12 (or an existing working 3.12 environment) with SSL. Choose an empty directory. The tool checks SSL/version before installing packages and does not modify QGIS Python. Some QGIS Python distributions lack standalone SSL support; use a normal Python distribution in that case. A successful setup records the backend executable in QGIS settings. The Python field may then be left blank. Authentication is never embedded.

## Individual tools

Each executable has a native Processing form and can be used in QGIS batch processing or the Model Designer. File inputs can be inferred from the processed-data root when left blank. For a raw raster stage, point the relevant input folder to the unpacked input archives. Semicolons separate multiple filenames, years or site names where a list is required. Parameters, filters and estimators are reused from the original backend; they are not new methods. The runner stages read-only copies for legacy CLIs which formerly read from their output directory. Logs and generated files remain in the chosen result folder; cancellation terminates the current backend process.

Use **02 Create the 19 study and reference sites as a QGIS layer** to create a point layer/GeoPackage in EPSG:4326. It contains study, background and replication-city groups, not validation observations.

## Earth Engine exports

Every original JS export has a Processing counterpart. The default writes an **offline task manifest**. Checking **Start remote Earth Engine export tasks** requires the Earth Engine API, your authenticated account and your authorized Cloud project. Authenticate in the prepared environment (for example `earthengine authenticate`) before submitting. Exports go to Google Drive, as the original scripts do; they are not downloaded directly into QGIS. No authentication or remote task submission was attempted here. Graph construction/serialization for all 13 definitions was tested against the official Earth Engine API test algorithm fixtures; remote dataset evaluation was not tested. Re-downloaded products may differ from the frozen archived exports.

The supplied static scripts use **WorldCover v200 (2021)** and a 2019–2025 monthly VIIRS mean. Those choices are preserved. Reconcile the manuscript's 2020/annual-product wording using SOURCE_DISCREPANCIES.txt. Do not silently change the input products.

## Validation and architecture

Native registration/execution, the 19-point layer, individual model/plot forms and all export previews were tested on installed QGIS 3.44.15 and 4.2.2. The extended QGIS replay passed 16 numerical comparisons; the new local backend also passed a QGIS-triggered quick replay. See VALIDATION_QGIS.json. Statistical kernels run in the bundled versioned backend using the selected Python executable. Thus these are native QGIS Processing entry points, not a claim that every numerical routine was rewritten using a different QGIS estimator.

Code is MIT; derived data are CC BY 4.0 subject to upstream source terms. Do not upload Local_QGIS_Runtime/ or Audit/ to GitHub or Zenodo.
