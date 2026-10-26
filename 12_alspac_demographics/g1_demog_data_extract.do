*** Script to extract G1 demogs for NMR data note
*** Author: Ruby Tsang and Laura Corbin
*** Date created: 2026-10-11

*** set up ***

version 19.5
capture log close
clear all
set more off, permanently
set maxvar 120000

log using "../../../../results/alspac_demographics/g1_demog_data_extract.log", replace text

*****************************************************************************************************************************************************************************************************************************.
	
use "/Volumes/Data/current/other/cohort profile/G0/mother/mz_6a.dta", clear
sort aln
gen in_mz=1
merge 1:1 aln using "/Volumes/Data/current/quest/G0/mother/a_3e.dta", nogen
merge 1:1 aln using "/Volumes/Data/current/quest/G0/mother/b_4f.dta", nogen
merge 1:1 aln using "/Volumes/Data/current/quest/G0/mother/c_8a.dta", nogen
keep aln mz010a preg_in_alsp preg_in_core preg_enrol_status mz005l mz005m mz013 mz014 mz028b ///
a006 a525 ///
b032 b650 b663 - b667 ///
c645a c755 c765 ///
bestgest


*keep only those pregnancies enrolled in ALSPAC
keep if preg_enrol_status < 3

* Dealing with withdrawal of consent: For this to work additional variables required have to be inserted before bestgest, so replace the *** line above with additional variables. 
* If none are required remember to delete the *** line.
* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file. Note that mother based WoCs are set to .a

order aln mz010a, first
order bestgest, last

do "/Volumes/Data/Syntax/Withdrawal of consent/mother_WoC.do"

* Check withdrawal of consent frequencies=37 and baseline number is 15447
tab1 mz010a, mis

save "../../../../data/demogs/mother_temp.dta", replace


*****************************************************************************************************************************************************************************************************************************.
* G1 Child BASED files - in this section the following file types need to be placed:
* Mother completed Qs about YP
* Obstetrics file OA

* ALWAYS KEEP THIS SECTION EVEN IF ONLY CHILD COMPLETED REQUESTED, although you will need to remove the *****

use "/Volumes/Data/current/other/cohort profile/G1/cp_4a.dta", clear
sort aln qlet
gen in_kz=1

keep aln qlet kz021 kz030 kz011b ///
in_core in_alsp in_phase2 in_phase3 in_phase4 in_phase5 tripquad


* Dealing with withdrawal of consent: For this to work additional variables required have to be inserted before in_core, so replace the ***** line with additional variables.
* If none are required remember to delete the ***** line.
* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file. Note that child based WoCs are set to .b


order aln qlet kz021, first
order in_alsp tripquad, last

do "/Volumes/Data/Syntax/Withdrawal of consent/child_based_WoC.do"

* Check withdrawal of consent frequencies child based=39 (two mums of twins have withdrawn consent)
tab1 kz011b, mis

drop kz021
save "../../../../data/demogs/childB.dta", replace

*****************************************************************************************************************************************************************************************************************************.
* G1 Child COMPLETED files - in this section the following file types need to be placed:
* YP completed Qs
* Puberty Qs
* Child clinic data
* Child biosamples data
* School Qs
* Obstetrics file OC
* G1 IMD for years goes in this section e.g. jan1999imd2010_crimeq5_YP
* Child longitudinal data

* NOTE: Always keep this section even if no child completed variables requested. Remove the ***** line.
* NOTE: Keep kz021 tripquad just to make the withdrawal of consent work - these are dropped for this file as the ones in the child BASED file are the important ones and should take priority

use "/Volumes/Data/current/other/cohort profile/G1/cp_4a.dta", clear
sort aln qlet
gen in_kz=1
merge 1:1 aln qlet using "/Volumes/Data/current/clinic/G1/f07_6a.dta", nogen
merge 1:1 aln qlet using "/Volumes/Data/current/clinic/G1/tf3_5a.dta", nogen
merge 1:1 aln qlet using "/Volumes/Data/current/clinic/G1/tf4_7a.dta", nogen
merge 1:1 aln qlet using "/Volumes/Data/current/clinic/G1/F24_7a.dta", nogen
merge 1:1 aln qlet using "/Volumes/Data/current/clinic/G1/F30_G1_2a.dta", nogen
* file built by Kate Northstone
merge 1:1 aln qlet using "../../../../data/demogs/bbs_age_months.dta", nogen 

keep aln qlet kz021 kz030 ///
f7003c fh0011a FJ003b FKAR0011 FLAR0010 age_BBS /// Age at attendance
tripquad

* Dealing with withdrawal of consent: For this to work additional variables required have to be inserted before tripquad, so replace the ***** line with additional variables.
* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file.  Note that mother based WoCs are set to .b

order aln qlet kz021, first
order tripquad, last

do "/Volumes/Data/Syntax/Withdrawal of consent/child_completed_WoC.do"

* Check withdrawal of consent frequencies child completed=37
tab1 kz021, mis

drop tripquad
save "../../../../data/demogs/childC.dta", replace

*****************************************************************************************************************************************************************************************************************************.
** Matching all data together and saving out the final file*.
* NOTE: any linkage data should be added here*.

use "../../../../data/demogs/childB.dta", clear
merge 1:1 aln qlet using "../../../../data/demogs/childC.dta", nogen
merge m:1 aln using "../../../../data/demogs/mother_temp.dta", nogen

* Remove non-alspac children.
drop if in_alsp!=1.

* Remove trips and quads.
drop if tripquad==1

drop in_alsp tripquad
save "../../../../data/demogs/G1_demogs.dta", replace

*****************************************************************************************************************************************************************************************************************************.
* QC checks*
use "../../../../data/demogs/G1_demogs.dta", clear

* Check that there are 15,746 records.
count

*****************************************************************************************************************************************************************************************************************************.

** remove intermediate files
rm  "../../../../data/demogs/mother_temp.dta"
rm "../../../../data/demogs/childB.dta"
rm "../../../../data/demogs/childC.dta"


capture log close

