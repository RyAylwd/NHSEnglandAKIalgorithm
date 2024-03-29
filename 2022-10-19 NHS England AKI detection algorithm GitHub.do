*NHS England AKI detection algorithm code 
*written by Ryan Aylward

*The following 3 variables are required in long format:
 *the variable 'creatinine' is C1 in the algorithm (umol/L)
 *the variable 'date' is the creatinine collection date of the blood sample. In this do file, date should be formated as %d. If date is in %tc format, then change all loops to: gen date_day = date/86400000 to convert milliseconds into days 
 *'pat' is the numeric patient ID: change this to your ID variable 
 
 *variables generated: rv1 ratio1 rv2 ratio2 rv ratio creat_diff dee warning_flag alert no_flag low_flag high_flag
 
********************************************************************************
*Data preparation

*read in original data
use "data.dta", clear // change this to your data file name 

*the long data format is first converted to a wide dataset and then re-combined with its own long data

*OPTIONAL: split dataset to reduce computational load if there are a large number of patients or creatinine values. If this is required, create a master do file that executes the algorithm in each split dataset separately, and then append all the split datasets once complete. 
cd "Split 1" // change working directory as desired 
keep if pat<= // for example if you have 100 000 patients, you might want to split into 2*50 000 each. it is faster to split and combine than to execute all patients at the same time
bysort pat (date) : gen n = _n // number of alerts per person

save "split1.dta", replace

 sum n // the maximum of n is the number of i=1/n that should be substituted for your specific data in the forvalues loop (below)
local number = r(max)

*Reshape wide, so you have date1, creatinine1, date2, creatinine2 etc
use "split1.dta", clear // Delete this line if splitting is not required
preserve 
reshape wide date creatinine, i(pat) j(n)
save "wide split1.dta", replace
restore

*Merge m:1 with the dataset in long format. So now you have each ID, date, creatinine, date1, creatinine1, date2, creatinine2 etc on each row
merge m:1 pat using "wide split1.dta"
drop _merge 

*clone dataset for each part below (will be re-combined at the end)
save "longwide split1.dta", replace // original
save "longwide 7d split1.dta", replace //0 - 7 days
save "longwide 365d split1.dta", replace //8 - 365 days

********************************************************************************
*1* Previous result within 365 days available?

*read in data
use "longwide split1.dta", clear

*determine if there is a creatinine value in the preceding 0-365days
gen date_day0 = date
forvalues i=1/`number' {  
gen date_day`i' = date`i'
} 

*loop over each creatinine to delete the creatinine value associated with itself
forvalues i=1/`number' {
replace creatinine`i'= . if creatinine`i'==creatinine & date_day`i'==date_day0 
}

*create variable of the difference in days
forvalues i=1/`number' { 
gen diff0365d`i' = date_day`i'-date_day0
}

*replace positive differences - become missing
forvalues i=1/`number' { 
replace diff0365d`i'=. if diff0365d`i'>0
}

*create a dummy variable that indicates difference is <365 days
forvalues i=1/`number' { 
gen diff365_yn`i' = 0
recode diff365_yn`i' (0=1) if diff0365d`i' >-365 & diff0365d`i' <0
}

*save data file
save "longwide split1.dta", replace

********************************************************************************
*2* Top right hand corner where there is no baseline creatinine available

*set lab reference interval - this depends on the assay used by the lab and patient sex, if sex is available as a variable, amend loops ie repeat for male and female
gen lab_upper_range = 
gen lab_lower_range = 

*if creatinine within RI and more than 365 days - no flag
gen no_flag = 0
forvalues i=1/`number' {
recode no_flag (0=1) if creatinine < lab_upper_range & creatinine >  lab_lower_range & diff365_yn`i'==0
} 

*if creatinine low - flag low
gen low_flag = 0
forvalues i=1/`number' { 
recode low_flag (0=1) if creatinine < lab_lower_range & diff365_yn`i'==0
}

*if creatinine high - flag high
gen high_flag = 0
forvalues i=1/`number' {  
recode high_flag (0=1) if creatinine > lab_upper_range & diff365_yn`i'==0 & creatinine!=.
}

*drop redundant variables
drop date_day0-diff365_yn`number'

*save data file
save "longwide topright split1.dta", replace

********************************************************************************
*3* 7 day criteria
*NB the variable 'rv' is the reference (preceding baseline) creatinine, either the lowest value within 7 days ('rv1') or median within 8 - 365 days ('rv2')

use "longwide 7d split1.dta", clear

gen date_day0 = date // if date is in milliseconds - divide /86400000
forvalues i=1/`number' { 
gen date_day`i' = date`i' // if date is in milliseconds - divide /86400000
} 

*loop over each creatinine to delete the creatinine value associated with itself
forvalues i=1/`number' {  
replace creatinine`i'= . if creatinine`i'==creatinine & date_day`i'==date_day0 
}

*create a new variable of the difference between dates in days
forvalues i=1/`number' { 
gen diff7d`i' = date_day`i'-date_day0
}

*replace 0 and positive differences because we are only interested in the 7days before the alert
forvalues i=1/`number' { 
replace diff7d`i'=. if diff7d`i'>=0
}

*create a dummy variable that indicates difference is 0 - 7 days
forvalues i=1/`number' {  
gen diff7_yn`i' = 0
recode diff7_yn`i' (0=1) if diff7d`i'>-8 & diff7d`i'<0
}

rename creatinine creat //so that creatinine is not included in creatinine* rowmin calculation 
*find minumum creatinine within 7 days
forvalues i=1/`number' { 
replace creatinine`i'=. if diff7_yn`i'==0 //only 0-7d creatinine values are used and rest are set to missing
}
egen rv1 = rowmin(creatinine*) //lowest creatinine value within 7d

*calculate the rv ratio
gen ratio1 = creat/rv1 

*drop redundant variables since these have been amended - will be restored later
drop creat creatinine1-diff7_yn`number'

*save dataset
save "longwide 7d split1.dta", replace

********************************************************************************
*4* 365 day criteria

use "longwide 365d split1.dta", clear

gen date_day0 = date // if date is in milliseconds - divide /86400000
forvalues i=1/`number' { 
gen date_day`i' = date`i' // if date is in milliseconds - divide /86400000
}  
 
*loop over each creatinine to delete the creatinine value associated with itself
forvalues i=1/`number' { 
replace creatinine`i'= . if creatinine`i'==creatinine & date_day`i'==date_day0 
}

*create a new variable of the difference between dates in days
forvalues i=1/`number' {  
gen diff8365d`i' = date_day`i'-date_day0
}

*replace 0 and positive differences
forvalues i=1/`number' {  
replace diff8365d`i'=. if diff8365d`i'>=0
}

*create a dummy variable that indicates difference is 8 - 365 days
forvalues i=1/`number' {  
gen diff8365_yn`i' = 0
recode diff8365_yn`i' (0=1) if diff8365d`i'>=-365 & diff8365d`i'<-8
}

rename creatinine creat // so that creatinine is not included in rowmedian creatinine* calculation 
*find median creat within 8 - 365 days
forvalues i=1/`number' { 
replace creatinine`i'=. if diff8365_yn`i'==0 //only 8-3665d creatinine values are used and any values outside of these dates are set to missing
}
egen rv2 = rowmedian(creatinine*) // median creatinine value within 8-365d

*calculate the rv ratio
gen ratio2 = creat/rv2 

*drop redundant variables - original values will be restored later
drop creat creatinine1-diff8365_yn`number'

*save dataset
save "longwide 365d split1.dta", replace

********************************************************************************
*combine datasets - 7d, 365d and original wide format (so all creatinine values and dates are restored)

merge m:1 pat n using "longwide 7d split1.dta" 
drop _merge
merge m:1 pat n using "longwide topright split1.dta"
drop _merge

***********
*5* compare ratio's  
gen ratio=ratio1 
replace ratio=ratio2 if ratio1==.
replace ratio=ratio2 if ratio2>ratio1 & ratio2!=. & ratio1!=. //use ratio2 if it is higher

*retain the creatinine rv1/2 as appropriate if ratio1/2 >1.5 and recode the value associated with lower ratio to missing
replace rv1=. if ratio2>ratio1 & ratio2!=. & ratio1!=.
replace rv2=. if ratio1>ratio2 & ratio1!=. & ratio2!=.

*create a single rv variable that takes on rv1/2, whichever ratio is highest
gen rv = rv1
replace rv = rv2 if ratio2>ratio1 & ratio2!=. & ratio1!=.
replace rv = rv2 if ratio==.

*calculate the difference between C1 and reference creatinine 
gen creat_diff = creatinine-rv 

*generate AKI alerts
gen alert=0 // no AKI

***********
*6* 48h criteria AKI 1 - called 'D' in algorithm, here called 'dee', 2 days is substituted for 48 hours. If %tc datetime available, convert milliseconds into hours by dividing date /3600000 and substituting days with hours. 

forvalues i=1/`number' {
gen date_day`i' = date`i'
}  

*create a new variable of the difference between dates
forvalues i=1/`number' {
gen diff2d`i' = date_day`i'-date
}

*replace 0 and positive differences with missings
forvalues i=1/`number' {
replace diff2d`i'=. if diff2d`i'>=0
}

*create a dummy variable that indicates difference is 0 - 2 days
forvalues i=1/`number' {
gen diff2d_yn`i' = 0
recode diff2d_yn`i' (0=1) if diff2d`i'>-2 & diff2d`i'<=0
}

rename creatinine creat //so that creatinine is not included in rowmin creatinine* calculation 
*find lowest creat within 48 h 
forvalues i=1/`number' {
replace creatinine`i'=. if diff2d_yn`i'==0 //only 0-2days creatinine values are used and any values outside of these dates are set to missing
}
egen low = rowmin(creatinine*) //lowest creatinine value within 2 days 
gen dee = creat - low // current creatinine minus the lowest creatinine in the last 2 days 

*generate alert flag AKI stage 1
forvalues i=1/`number' {
recode central_alert (0=1) if ratio<1.5 & dee>26 & dee!=. & diff2d`i'>=-2 & diff2d`i'<=0
} // AKI stage 1 

*create variable that indicates that increase within 7 days but not within 2 days 
gen warning_flag = 0
forvalues i=1/`number' {
recode warning_flag (0=1) if diff2d`i'>-7 & diff2d`i'<-2 & ratio<1.5 & creat_diff>26 
}

***********
*7* AKI 3
*if the difference between creatinine values C1 - rv1/2 is >354umol/l - flag as stage 3
recode alert (0=3) if ratio>=1.5 & creat_diff>354 & creat_diff!=. & ratio!=.
recode alert (0=3) if ratio>=3 & ratio!=.

***********
*8* AKI 2
recode alert (0=2) if ratio>= 2 & ratio<3 

***********
*9* AKI 1
recode alert (0=1) if ratio>= 1.5 & ratio<2

*drop redundant variables
keep pat date creat rv2 ratio2 rv1 ratio1 no_flag low_flag high_flag ratio rv central_alert warning_flag creat_diff dee 
rename creat creatinine 

**********
order pat date creatinine rv1 ratio1 rv2 ratio2 rv ratio creat_diff dee warning_flag alert no_flag low_flag high_flag

*save alert data file and delete files no longer required 

save "alert split1.dta", replace // append with other splits if necessary
erase "split1.dta"
erase "wide split1.dta"
erase "longwide split1.dta"
erase "longwide 7d split1.dta"
erase "longwide 365d split1.dta"
erase "longwide topright split1.dta"
exit, clear 
