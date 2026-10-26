## Preparation of ALSPAC NMR data for release (B4762)
## As described in the data note in Wellcome Open titled: "Life course NMR data in ALSPAC"

### Code to reformat and fix any issues with files provided by Nightingale, and summarise 

01_prep_metadata_files.R - code to reformat feature metadata file provided by Nightingale with F30 data so that it resembles a post 'ng_anno' merged data file (produced by metaboprep during processing). Essentially updating feature annotations. Input is worksheet 4 from the Nightingale F30 release. Output is a feature metadata file formatted to read into metaboprep.

02_fix_g2_ids.R - there was an issue with the sample IDs in the G2 dataset when it was returned from Nightingale originally. Some of the IDs had been replaced with a rounded value. They corrected this afterwards but only in the data file, not in the corresponding metadata files. Hence this script is to replace the IDs across all files in the set to match the updates.

03_extract_f30.R - script to deconstruct excel file provided by Nighitngale for F30. The format received was different for this one because it was a newly generated dataset - not an update like the other timepoints.

04_fix_tf4_ids.R - the TF4 dataset contains duplicates, i.e., samples that were analysed twice and therefore appear twuce in the dataset. This script identifies the duplicates and removes one of each pair (based on missingness).

05_omiprep.R - Code to run data through omiprep for summary and QC checks 


### Code to build files for return to ALSPAC

06_build_datasets.R - takes output from omiprep and converts to alspac compatible format (data linked to sample ID)

07_link_datasets.R - links sample IDs to internal ALSPAC participant IDs

08_check_dupes.R - looks for instances where multiple samples link to the same individual

09_validate.R - compares NMR derived to clinical assay derived measures where there is overlap. Select which sample data to keep when there are two or three samples linked ot the same individual.

10_build_alspac_datasets.R - build alspac compatible datasets linked to participant ID

11_build_feature_metadata.R - build feature metadata files for making openly available. Contain summary statistics from omiprep.


### Code to generate demographic summaries for the data note (age & sex)

alspac_demographics/g0_demog_data_extract.do - uses ALSPAC syntax from data extraction to extract age data from FOM1 and FOM2 clinics

alspac_demographics/g1_demog_data_extract.do - uses ALSPAC syntax from data extraction to extract age and sex data from G1 clinics

alspac_demographics/generate_sumstats.do - summarises age and sex by clinic for individuals with NMR data


