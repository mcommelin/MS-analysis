# MS-analysis

# Intro
This code provides an analysis template with specific functions, to streamline and homogenize the analysis of LC or GC MS analysis. All code is written with R. Because the analysis involves several expert decisions it is not fully automated. Each user will have to specify choices and prepare relevant meta-data files.

## Load data
As a first step the raw data from the analysis machine has to be loaded. Besides that a second file adding relevant variables is needed to perform the correct analysis.

### raw data files
As Raw data input a folder with raw data files is expected were the data per compound is all printed below each other. This is the default output from the MassLynx program.
To load this raw data the function 'load_raw_data' was developed. (see functions section)

### meta data file
The meta-data contains three relevant columns which add for the analysis. Because the content of this data is often user specific, this has to be prepared seperatly. The expected column names are (these names need to be exactly the same for follow up calculations)
  
  - tqs_code          The tqs_code is used as identifier for all samples from the raw data
  - an_type           The 'type' of the sample for analysis
  - correction_factor The dilution of the sample within during the extraction process and in the vial
  - dilution_factor   The extra dilutions done to dilute high concentrations to be within the measurement domain of the MS machine.
  
The analysis types (an_type) are divided in the following groups: sample, cal, standard, qc and blank.

  - sample    This are the samples in the measurement
  - cal       The calibration curve samples, the concentration of the specific calibration points should be part of the 'sample_text' variable e.g. (Cal_6.25_ng/mL)
  - standard  a single standard point to check stability of the machine, not included in the calibration curve
  - qc        Sample with spiked concentrations to check the recovery for specific matrixes
  - blank     Milliq water vials, to clean and/or check the bias of the machine.
  
Important: use these specific analysis type names and column names so follow up functions will work properly. The metadata file can be added to the raw data with the function 'meta_data_add' (see function section).

## Functions

### load_raw_data()

load_raw_data(dir, delim = "\t", tqs_code = "Name", sample_text = "Sample.Text", RT = "RT")

  - dir: the directory were the txt file(s) are stored. Each file is the result of 1 machine batch.
  - delim: the delimiter in the raw data files, default is "\t"
  - tqs_code, sample_text, RT - header names for these three columns, used to homogenize data for further analysis.

The function expects to find the word 'area' in every column header which stores peak area data, it will then rename these automatically.

### meta_data_add()

meta_data_add(meta_file, df_data, load = TRUE)

  - meta_file: either a file location e.g. ("data/meta_file.txt"), or a dataframe in R environment
  - df_data: the result of load_raw_data()
  - load: TRUE = load file from specified location, FALSE = the data is already in R.
