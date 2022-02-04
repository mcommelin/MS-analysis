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

## Calculate concentrations
With the data loaded the conctrations for a specific compound can be calculated based on the area's measured by the machine. Several quality checks are built in to select the samples with a good response, these are: Ion ratio, Retention Time and Internal standard response. If a sample passes all three checks it is assumed to be valid.

The Ion ratio (IR) is defined as the ratio between the main and second response peaks, this is compound specific. The mean IR is calculated based on the IR values of the calibration curve standards. The IR check is done with: (1 - alpha) * IRmean < IR_sample < (1 + alpha) * IRmean; (default alpha = 0.3). The same type of evaluation formulas are used for retention time and internal standard response. Assuming that the standards have good response in the machine.

### lineair model fit
The calibration curve is used to fit a linear model to relate peak area to concentration level. For each compound the MS machine has a sensitive range in which the relation is linear, accurate analysis can only be done within the linear domain. As first step the linearity of the calibration measurements is tested against the response of the median value in the calibration range, sample that deviate more than the delta.linearity factor (range 0.2 - 0.3) will not be included in the linear model estimate. Besides that sample concentrations that fall outside the linear domain can not be analysed further.

A linear model is fit to the calibration points that fall within the linear domain. And the slope and intercept of this model are used to calculate the concentrations of the samples.

### correction for dilutions
Two different types of dilutions are used for the samples. The first (correction_factor) are the dilutions which are done during the standard laboratory preparations. For example: for soil matrix samples (S), 5 mL of H2O, 10 ml of ACN is added (A). From this extract only the ACN layer is pipetted, but a mixing with water of about 10% (B)is estimated. After that the extract is added 1:1 with milliq water (C) in a vial. This leads to a dilution factor of; 1/(S/(A*B)/C).

The second type of dilutions is done to let the expected concentration of the sample fall within the linear domain of the MS machine. This is done by additionally diluting the extract. Because the linear domain of the machine is at least a factor 10, dilutions can also be made in steps of a factor 10.


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
