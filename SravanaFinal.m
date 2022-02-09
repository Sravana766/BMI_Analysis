clear; clc; close all

%% Load the subject's data and define the variables
data = load('OrthoData.mat');
time = data.b1(:,1);
ecgraw = data.b1(:,2);
Fs = 200; %Sampling frequency of 200 samples/sec
height = input('Enter your height here [m]: '); %Create an input variable for subject to enter height
weight = input('Enter your weight here [kg]: '); %Create an input variable for subject to enter weight

%% Signal Analysis
figure(1)
plot(time,ecgraw); %Plot the original ECG signal
title('Raw ECG');
xlabel('Time [sec]');
ylabel('Amplitude [mV]');

%FFT for raw ECG
n_raw = length(ecgraw);
freq_raw = Fs*(0:n_raw-1)/n_raw; %Establish the frequency spectrum of the raw ECG signal
fftECG = fft(ecgraw,n_raw); %Compute the fft of the raw signal to determine frequencies of noise

figure(2)
plot(freq_raw,fftECG); %Plot fft vs. frequency
title('Frequency Spectrum of Raw ECG');
xlabel('Frequency [Hz]');
ylabel('Amplitude [mV]');

%Filtering
Wn_low = 60/(Fs/2); %Establish a cutoff frequency for the lowpass filter
Wn_high = 0.67/(Fs/2); %Establish a cutoff frequency for the highpass filter
lowPass = fir1(300, Wn_low,'low'); %Establish a lowpass FIR filter design
highPass = fir1(300,Wn_high,'high'); %Establish a highpass FIR filter design
low_pass_filter = filtfilt(lowPass,1,ecgraw); %Filter the raw ECG through the lowpass filter
filteredECG = filtfilt(highPass,1,low_pass_filter); %Filter the lowpass filtered ECG signal through the highpass filter

figure(3)
plot(time,filteredECG); %Plot the filtered ECG signal vs. time
title('Filtered ECG');
xlabel('Time [sec]');
ylabel('Amplitude [mV]');

%FFT for filtered ECG
n_filt = length(filteredECG);
freq_filt = Fs*(0:n_filt-1)/n_filt; %Establish the frequency spectrum of the raw ECG signal
fftECG_filtered = fft(filteredECG,n_filt); %Compute the fft of the filtered signal to check for noise

figure(4)
plot(freq_filt,fftECG_filtered); %Plot the filtered ECG signal's FFT vs. frequency
title('Frequency Spectrum of Filtered ECG');
xlabel('Frequency [Hz]');
ylabel('Amplitude [mV]');

%% Find peaks
%The following 3 sections of code are broken up into the three segments of
%data: resting, stand, and recovery. The 'MinPeakHeight,'
%'MinPeakDistance,' and indexes within the findpeaks function vary per
%subject and must be changed.

StandingTime = 119877;
EndofStanding = StandingTime+200*20;

%% Initial Steady State

[R_waves_resting,locs_resting] = findpeaks(filteredECG(1:StandingTime-1),'MinPeakHeight',0.09,'MinPeakDistance',125); %index varies from the first sample to the sample before initial rise
time_peaks_resting = time(locs_resting); %Change the location of each peak into seconds based on the given time scale
difference_resting = diff(time_peaks_resting); %Calculate the difference between peaks (= 1 heart beat); accounts for HRV

figure(5)
plot(time, filteredECG);
title('Filtered ECG with Peaks');
xlabel('Time [sec]');
ylabel('Amplitude [mV]');
hold on
plot(time_peaks_resting,R_waves_resting,'ro'); %Plot the location of peaks on top of the filtered ECG to check for accuracy

heart_rate_resting = 60./difference_resting; %Calculate heart rate in bpm

for  kk = 1:(length(heart_rate_resting))
    if kk == 1  || kk == 2 
        movingAvg_resting(kk) = heart_rate_resting(kk);
    elseif  kk == (length(heart_rate_resting) - 2) || kk == (length(heart_rate_resting) - 1) || kk == length(heart_rate_resting)
         movingAvg_resting(kk) = heart_rate_resting(kk);
    else
        movingAvg_resting(kk) = mean(heart_rate_resting(kk-2:kk+2));
    end   
end


Timeresting = 0:599.3750/length(movingAvg_resting):599.3750;

figure(6)
stairs(Timeresting(:,1:end-1),movingAvg_resting); %Plot heart rate on a stairs plot to better show variability
title('Heart Rate during Supine Position');
xlabel('Samples');
ylabel('Heart Rate [bpm]');


%% 10 s before stand + Time interval of stand + 10-20s Recovery = 600.868s - 620.868s

[R_waves_stand,locs_stand] = findpeaks(filteredECG(StandingTime:EndofStanding),'MinPeakHeight',0.1,'MinPeakDistance',73); %index varies from the first sample during rise to the 20th second post-rise
time_peaks_stand = time(StandingTime + locs_stand); %Add locations of peaks to the position of the first sample of this segment; correctly plots the new peaks
difference_stand = diff(time_peaks_stand); %Calculate the difference between peaks


figure(7)
plot(time, filteredECG);
title('Filtered ECG with Peaks');
xlabel('Time [sec]');
ylabel('Amplitude [mV]');
hold on
plot(time_peaks_stand,R_waves_stand,'ro');


heart_rate_stand = 60./difference_stand; %Calculate heart rate in bpm

for  kk = 1:(length(heart_rate_stand))
    if kk == 1  || kk == 2 
        movingAvg_stand(kk) = (heart_rate_stand(kk));
    elseif  kk == (length(heart_rate_stand) - 2) || kk == (length(heart_rate_stand) - 1) || kk == length(heart_rate_stand)
         movingAvg_stand(kk) = heart_rate_stand(kk);
    else
        movingAvg_stand(kk) = mean(heart_rate_stand(kk-2:kk+2));
    end   
end

Timestand = 599.38:(20/length(movingAvg_stand)):619.38; %changes for every subject


figure(8) 
stairs(Timestand(:,1:end-1),movingAvg_stand); %Plot heart rate on a stairs plot (this segment is strictly 20 seconds for each subject, but it is displayed in samples)
title('Heart Rate during Stand to 20 s Post-Stand');
xlabel('Samples');
ylabel('Heart Rate [bpm]');


%% Recovery Steady State
[R_waves_after,locs_after] = findpeaks(filteredECG(EndofStanding+1:end),'MinPeakHeight',0.13,'MinPeakDistance',107); %index varies from the first sample after the 20 s mark to the end of the signal
time_peaks_after = time(EndofStanding+1+locs_after); %Add locations of peaks to the position of the first sample of this segment; correctly plots the new peaks
difference_after = diff(time_peaks_after); %Calculate the difference between peaks

figure(9)
plot(time,filteredECG);
title('Filtered ECG with Peaks');
xlabel('Time [sec]');
ylabel('Amplitude [mV]');
hold on
plot(time_peaks_after, R_waves_after, 'ro');

heart_rate_after = 60./difference_after; %Calculate the heart rate

for  kk = 1:(length(heart_rate_after))
    if kk == 1  || kk == 2 
        movingAvg_after(kk) = (heart_rate_after(kk));
    elseif  kk == (length(heart_rate_after) - 2) || kk == (length(heart_rate_after) - 1) || kk == length(heart_rate_after)
         movingAvg_after(kk) = heart_rate_after(kk);
    else
        movingAvg_after(kk) = mean(heart_rate_after(kk-2:kk+2));
    end   
end

Timeafter = 619.385:(899.9450-619.385)/(length(movingAvg_after)):899.9450;

figure(10)
stairs(Timeafter(:,1:end-1),movingAvg_after);%Plot heart rate on a stairs plot 
title('Heart Rate after 20 s Post-Stand Mark');
xlabel('Samples');
ylabel('Heart Rate [bpm]');

xy = length(movingAvg_resting)+length(movingAvg_stand)+length(movingAvg_after);

WholeHeartRate = ones(1,xy); 
WholeHeartRate(1:length(movingAvg_resting)) = movingAvg_resting; 
WholeHeartRate(length(movingAvg_resting)+1:length(movingAvg_resting)+length(movingAvg_stand)) = movingAvg_stand(1:end); 
WholeHeartRate(length(movingAvg_resting)+length(movingAvg_stand)+1:end) = movingAvg_after(1:end); 

TimeWhole = ones(1,xy); 
TimeWhole(1:length(movingAvg_resting)) = Timeresting(1:end-1); 
TimeWhole(length(movingAvg_resting)+1:length(movingAvg_resting)+length(movingAvg_stand)) = Timestand(1:end-1); 
TimeWhole(length(movingAvg_resting)+length(movingAvg_stand)+1:end) = Timeafter(1:end-1); 

figure(11)
stairs(TimeWhole,WholeHeartRate);
title('Heart Rate Over 15 Minutes');
xlabel('Time [s]');
ylabel('Heart Rate [bpm]');

figure(14)
stairs(TimeWhole,WholeHeartRate,'k'); 
title('Heart Rate Over 15 Minutes');
xlabel('Time [s]');
ylabel('Heart Rate [bpm]');
axis([580 640 0 180]);
hold on
xline(599.38,'g'); 
xline(609.38,'r'); 
xline(619.38,'c');
hold off 
legend('Heart Rate','Stand Time', '10 s' , '20 s'); 


%% Determine heart rate values at specific points (only to plot)

% Average heart rate at resting steady state
HR_resting = mean(heart_rate_resting); %The mean is calculated for visualization purposes
fprintf('The Resting Heart Rate in Supine Position is: %f \n\n',HR_resting);
%Stand start time 
HR_standstart = movingAvg_stand(1); %Index corresponds to the intitial rise from orthostatic position; will be the first value in "heart_rate_stand"
fprintf('The Heart Rate at the Start of Stand is: %f \n\n',HR_standstart);
%Max Heart Rate
HR_max = max(movingAvg_stand); %Calculates maximum heart rate (should occur during stand)
fprintf('The Max Heart Rate is: %f \n\n',HR_max);
%Heart Rate at 10s
HR_10s = movingAvg_stand(15); %index changes for every subject; corresponds to the middle point in "heart_rate_stand" since the segment is only 20 s
fprintf('The Heart Rate at 10s after Stand is: %f \n\n',HR_10s);

%Heart Rate Recovery at 20s
HR_20s = movingAvg_stand(28); %index changes for every subject; corresponds to the last value in the 20 s segment
fprintf('The Heart Rate at 20s after Stand is: %f \n\n',HR_20s);

%Difference in heart rate recovery
dHR_10 = HR_10s - HR_resting; %Calculates the change in heart rate from baseline (resting) at 10 s
dHR_20 = HR_20s - HR_resting; %Calculates the change in heart rate from baseline (resting) at 20 s

HRR_speed_10_20 = dHR_20 - dHR_10; %Calculate the speed of heart rate recovery during 10 and 20 seconds post-stand
fprintf('The Heart Rate Recovery during 10 to 20s of standing is: %f \n\n',HRR_speed_10_20);

%Recovery Steady State
HR_recovery = mean(movingAvg_after); %The mean is calculated for visualization purposes
fprintf('The Average Heart Rate in a standing steady state is: %f \n\n',HR_recovery);

%% Plot HR values at certain times
HR_value = [HR_resting HR_resting HR_standstart HR_10s HR_20s HR_recovery HR_recovery]; %Establish a vector with heart rate values at the specified points
time_HR_value = [-300 -1 0 4 10 20 100]; %Establish a horizontal axis for time (in seconds) at the specified points

figure(12)
plot(time_HR_value,HR_value,'-*'); %Plot time vs HR to visualize heart recovery
title('HRR 10 to 20 seconds after the Orthostatic Challenge');
xlabel('Title [sec]');
ylabel('Heart Rate [bpm]');
hold on
xline(time_HR_value(4),'g'); %Create a vertical line to show where 10 seconds post-stand occurs
xline(time_HR_value(5),'r'); %Create a vertical line to show where 20 seconds post-stand occurs
legend('Heart Rate vs. Time', 'HR 10 s', 'HR 20 s');

%% BMI Calculation

BMI = weight/(height^2); %Calculate BMI

if (BMI >= 18.5) && (BMI <= 24.9) && (abs(HRR_speed_10_20) > 15.31)
    disp('Hypothesis accepted.');
else
    disp('Hypothesis rejected.');
end
