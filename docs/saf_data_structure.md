# SAF Q2 2025 Data Structure

This project expects the SRTR SAF folder to be available as `SAF Q2 2025` in Box.
The scripts auto-detect common locations on Windows and macOS, or use the
`SAF_Q2_2025_DIR` environment variable when set.

## Expected Root Folders

- Windows Box default: `C:/Users/<user>/Box/SAF Q2 2025`
- macOS Box default: `/Users/<user>/Library/CloudStorage/Box-Box/SAF Q2 2025`

Required subfolders:

- `pubsaf2506`: core public SAF SAS files.
- `SupplementalData2506`: supplemental geographic linkage files.

Reference files in the root folder:

- `dataDictionary.html`: SAF table and variable dictionary.
- `SAFsLinkingDiagram.pdf`: SAF table linking diagram.

## Tables Used By This Pipeline

Candidate registration tables in `pubsaf2506`:

- `cand_kipa.sas7bdat`: kidney/pancreas candidates.
- `cand_liin.sas7bdat`: liver/intestine candidates.
- `cand_thor.sas7bdat`: thoracic candidates.

The data dictionary describes the candidate tables as waitlist registration
records. The analysis uses `PERS_ID` and `PX_ID` as person and registration
identifiers, `WL_ORG` for organ, listing and activation dates, removal dates
and removal codes, candidate source, and demographic/model covariates such as
age at listing, sex, SRTR race, SRTR ethnicity, BMI, medical condition, and
listing center.

Candidate ZIP supplement in `SupplementalData2506`:

- `canzip2506.sas7bdat`: one candidate ZIP linkage file keyed by `PERS_ID`
  and `PX_ID`.

Fields used:

- `CAN_PERM_ZIP`: permanent candidate ZIP used for pollution exposure linkage.
- `CAN_WAIT_IN_PERM_ZIP`: whether the candidate waited in the permanent ZIP.

Status history tables in `pubsaf2506`:

- `stathist_kipa.sas7bdat`
- `stathist_liin.sas7bdat`
- `stathist_thor.sas7bdat`

The time-varying analysis joins status history by `PX_ID` and `WL_ORG`. It uses
`CANHX_BEGIN_DT`, `CANHX_END_DT`, and `CANHX_STAT_CD` to split follow-up into
active/inactive status intervals.

## Runtime Path Resolution

Shared path logic lives in `code/saf_paths.R`.

Override auto-detection by setting:

```r
Sys.setenv(SAF_Q2_2025_DIR = "C:/Users/<user>/Box/SAF Q2 2025")
```

or, from a shell before running `Rscript`:

```bash
export SAF_Q2_2025_DIR="/Users/<user>/Library/CloudStorage/Box-Box/SAF Q2 2025"
```
