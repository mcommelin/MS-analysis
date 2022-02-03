# MS-analysis

# Intro
This code provides an analysis template with specific functions, to streamline and homogenize the analysis of LC or GC MS analysis. All code is written with R. Because the analysis involves several expert decisions it is not fully automated. Each user will have to specify choices and prepare relevant meta-data files.

## Load data
As a first step the raw data from the analysis machine has to be loaded. Besides that a second file adding relevant variables is needed to perform the correct analysis.

### raw data files
As Raw data input a .txt file is expected were the data per compound is all printed below each other. This is the default output from the MassLynx program.
To load this raw data the function 'load_raw_data' was developed.


## Functions

### load_raw_data()

load_raw_data(dir, delim = "\t", meta = TRUE, tqs_code = "Name", sample_text = "Sample.Text", RT = "RT")

  - dir: the directory were the txt file(s) are stored. Each file is the result of 1 machine batch.
  - delim: the delimiter in the raw data files, default is "\t"
  - meta: does the 'sample.text' column contain information on analysis type, default is TRUE
  - tqs_code, sample_text, RT - header names for these three columns, used to homogenize data for further analysis.

The function expects to find the word 'area' in every column header which stores peak area data, it will then rename these automatically.