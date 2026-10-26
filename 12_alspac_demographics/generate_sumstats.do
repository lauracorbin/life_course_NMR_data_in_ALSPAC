*** Script to generate age and sex summaries for NMR data note
*** Author: Laura Corbin
*** Date created: 2026-10-11

*** set up ***

version 19.5
capture log close
clear all
set more off, permanently
set maxvar 120000

log using "../../../../results/alspac_demographics/generate_sumstats.log", replace text

*****************************************************************************************************************************************************************************************************************************.

** G0

** read in nmr data (built file from Gail White)

use "../../../../data/alspac_release/post_ALSPACprocessing/v2/G0_NMR_1a.dta", clear
count
sort aln
* keep flags and ids only
keep aln has_sample_*

** merge in demographics keeping only those in NMR file
merge 1:1 aln using "../../../../data/demogs/G0_demogs.dta", keep(match)

** summarise age at clinics
summ fm1a011 if fm1a011 <. & fm1a011 > 0 & has_sample_at_FOM1 == 1
summ fm2a011 if fm2a011 <. & fm2a011 > 0 & has_sample_at_FOM2 == 1

*****************************************************************************************************************************************************************************************************************************.

** G1

** read in nmr data (built file from Gail White)

use "../../../../data/alspac_release/post_ALSPACprocessing/v2/G1_NMR_1a.dta", clear
count
sort aln qlet
* remove duplicate caused by sample duplictae in tf4
*duplicates drop aln qlet, force

* keep flags and ids only
keep aln qlet has_sample_*

** merge in demographics keeping only those in NMR file
merge 1:1 aln qlet using "../../../../data/demogs/G1_demogs.dta", keep(match)

** convert to years where needed
gen f7_age = f7003c/12
gen tf3_age = fh0011a/12
gen bbs_age_yrs = age_BBS/12

** summarise age at clinics (years)
summ f7_age if f7_age <. & f7_age > 0 & has_sample_at_F7 == 1
summ tf3_age if tf3_age <. & tf3_age > 0 & has_sample_at_TF3 == 1
summ FJ003b if FJ003b <. & FJ003b > 0 & has_sample_at_TF4 == 1
summ FKAR0011 if FKAR0011 <. & FKAR0011 > 0 & has_sample_at_F24 == 1
summ FLAR0010 if FLAR0010 <. & FLAR0010 > 0 & has_sample_at_F30 == 1
summ bbs_age_yrs if bbs_age_yrs <. & bbs_age_yrs > 0 & (has_sample_at_BBS_Baseline == 1 | has_sample_at_BBS_PostGluc == 1)


** summarise sex distribution
tab kz021 if has_sample_at_F7 == 1, missing
tab kz021 if has_sample_at_TF3 == 1, missing
tab kz021 if has_sample_at_TF4 == 1, missing
tab kz021 if has_sample_at_F24 == 1, missing
tab kz021 if has_sample_at_F30 == 1, missing
tab kz021 if (has_sample_at_BBS_Baseline == 1 | has_sample_at_BBS_PostGluc == 1), missing


*****************************************************************************************************************************************************************************************************************************.
capture log close
     