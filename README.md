# Body Mass Index Analysis

## Introduction
The speed of heart rate recovery (HRR) is an important biomarker that has been recently studied to determine its efficacy as a prognostic tool for cardiovascular disease (CVD) and mortality. This MATLAB script tries correlate invidual BMI and HRR. The data that was used to perform this analysis is the ECG recodring that was collected thorughout the experiment. The subject in this experiment was at rest, lying on the floor, for 10 seconds. Next, the subject went to upright standing position. FInally, the subject stood still for 5 seconds. The ECG was collected thorughout the entire duration. The ECG was collected using the iWorx Kit Equipment. The signal analysis was completed in MATLAB

## Procedure
There are three main segments in this analysis:
1. Filtering
2. Calculation of Heart Rate
3. Plotting

### Filtering
- The Sampling Frequency of the DAQ is 200 samples/sec
- FFT analysis was performed to analyze the frequency components
  - After analyzing the frequency spectrum, the ECG data is a low frequencency data
- A low pass filter with a cutoff frequency of 60Hz was used to preseve the data and remove the noise

### Calculation of the Heart Rate
- The **R-waves** from the ECG trace were detected using peak dection algorithm
- Heart Rate was indirectly measured by counting the number of **R-waves** in 60 seconds

### Plotting
 - Data was visualized using graphs and figures for better understanding

## Usage
To run this program, download the *SravanaFinal.m* file and the *OrthoData_Sravana.mat* files. 

