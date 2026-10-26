*** Script to extract G0 demogs for NMR data note
*** Author: Ruby Tsang and Laura Corbin
*** Date created: 2026-10-11

*** set up ***

version 19.5
capture log close
clear all
set more off, permanently
set maxvar 120000

log using "../../../../results/alspac_demographics/g0_demog_data_extract.log", replace text

*** extract data ***

* Mother based files - include here all files related to the mother
	
use "/Volumes/Data/current/other/cohort profile/G0/mother/mz_6a.dta", clear
sort aln
gen in_mz=1
merge 1:1 aln using "/Volumes/Data/current/quest/G0/mother/a_3e.dta", nogen
merge 1:1 aln using "/Volumes/Data/current/quest/G0/mother/b_4f.dta", nogen
merge 1:1 aln using "/Volumes/Data/current/quest/G0/mother/c_8a.dta", nogen
merge 1:1 aln using "/Volumes/Data/current/clinic/G0/mother/FOM1_3c.dta", nogen
merge 1:1 aln using "/Volumes/Data/current/clinic/G0/mother/FOM2_5b.dta", nogen


keep aln mz010a mum_in_alsp mum_in_core mum_enrol_status mz005l mz005m mz013 mz014 mz028b mum_and_preg_enrolled ///
a006 a525 ///
b032 b650 b663 - b667 ///
c645a c755 c765  ///
fm1a011 fm2a011 /// Age at attendance
bestgest

*keep in just those enrolled mothers and/or pregnancies
keep if mum_and_preg_enrolled <2

* Dealing with withdrawal of consent: For this to work additional variables required have to be inserted before bestgest, so please insert those as indicated above. 
* If none are required please delete the comments  line above bestgest.

* An additional do file is called in to set those withdrawing consent to missing so that this is always up to date whenever you run this do file. Note that mother based WoCs are set to .a

order aln mz010a, first
order bestgest, last

do "/Volumes/Data/Syntax/Withdrawal of consent/mother_WoC.do"

*Check withdrawal of consent frequencies=37
tab1 mum_enrol_status, mis

*Check there is a total of n=15,236 pregnancies/women   (14,833 individual women when we inlcude those who have withdrawn consent)
tab mum_enrol_status mz005l, mis


save "../../../../data/demogs/G0_demogs.dta", replace

* QC checks
use "../../../../data/demogs/G0_demogs.dta", clear

* Check that there are 15236 records.
count

capture log close


