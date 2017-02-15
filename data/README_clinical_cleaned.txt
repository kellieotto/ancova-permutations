This README_clinical_cleaned.txt file was created on 20161209 by Kellie Ottoboni

-------------------
GENERAL INFORMATION
-------------------


1. Title of Dataset 
	GERD clinical trial results

2. Author Information

    Kellie Ottoboni
	Dept. of Statistics, University of California, Berkeley
	367 Evans Hall
	Berkeley, CA 94720
	kellieotto@berkeley.edu
	
    Luigi Salmaso
	Department of Management and Engineering, University of Padova
	Stradella S. Nicola, 3
	36100 Vicenza, Italy
	luigi.salmaso@unipd.it
	
    Fraser Lewis
	Medical Affairs and Evidence Generation, Reckitt Benckiser
	103-105 Bath Road
	Slough, Berkshire SL1 3UH United Kingdom
	fraser.lewis@rb.com
	
3. Date of data collection 
	Unknown

4. Geographic location of data collection (where was data collected?): 
	8 anonymized sites in 2 countries

5. Information about funding sources that supported the collection of the data:
	Reckitt Benckiser collected this data for a clinical trial.


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 


1. Licenses/restrictions placed on the data:
TODO

2. Links to publications that cite or use the data:
TODO

3. Links to other publicly accessible locations of the data:
TODO

4. Links/relationships to ancillary data sets:
N/A

5. Was data derived from another source?
    Data was derived from a proprietary dataset belonging to RB.
	The R code used to derive this dataset from the original one is included in clinical_cleaned.R.

6. Recommended citation for the data:
TODO



---------------------
DATA & FILE OVERVIEW
---------------------


Filename: clinical_cleaned.csv
Description: 
	The data were collected in a clinical trial comparing the effectiveness of treatments A and B for treating symptoms of GERD.  The trial was conducted at 8 different sites in 2 countries.  Sites had varying numbers of patients.  Treatment was assigned at random to participants at each site.  Individuals were measured for 14 days: the first seven were baseline measures and the last seven were after treatment was administered.  Measurements were taken by a survey that asked questions about the severity and frequency of heartburn, regurgitation, and dyspepsia.



--------------------------
METHODOLOGICAL INFORMATION
--------------------------


1. Description of methods used for collection/generation of data: 
	Individuals were surveyed for 14 days. The survey is available at https://github.com/kellieotto/ancova-permutations/blob/master/data/questionnaire.pdf.


2. Methods for processing the data: 
	Preprocessing was done in the file clinical_cleaned.R.
	The original data included 14 observations per patient.
	The dataset we provide has 2 observations per patient, where each row is the average of 7 rows in the the original data.

3. Instrument- or software-specific information needed to interpret the data:
	Not applicable.
	
4. Standards and calibration information, if appropriate:
	Not applicable.
	
5. Environmental/experimental conditions:
	Two treatments, called A and B here.

6. Describe any quality-assurance procedures performed on the data:
	Not applicable.
	
7. People involved with sample collection, processing, analysis and/or submission:
	Not applicable.

---------------------------------------------------
DATA-SPECIFIC INFORMATION FOR clinical_cleaned.csv
---------------------------------------------------

1. Number of variables: 
	15, including column of subject IDs

2. Number of cases/rows: 
	273 rows including header; 136 subjects

3. Variable List
	SUBJID: Unique identifier for individuals. Integers 1, ..., 136.
	SITEID: Unique identifier for location. Integers 1, ..., 8.
	VISITNUM: Week of observation. 1 (pre-treatment) or 2 (post-treatment).
	tr: treatment. A or B.
	country: Unique identifier for country. 1 or 2.
	heart_sev: Average heartburn severity over the week. Range 0 to 3.
	regurg_sev: Average regurgitation severity over the week. Range 0 to 3.
	dysp_sev: Average dyspepsia severity over the week. Range 0 to 3.
	heart_freq: Average heartburn frequency over the week. Range 0 to 10.
	regurg_freq: Average regurgitation frequency over the week. Range 0 to 10.
	dysp_freq: Average dyspepsia frequency over the week. Range 0 to 10.
	daily_heart:
	daily_regurg:
	daily_hrdq: HRDQ score = severity * frequency/(frequency+3).
	daily_dysp:


4. Missing data codes:
	Not applicable. Missing data were removed in preprocessing.

5. Specialized formats of other abbreviations used
	Not applicable.