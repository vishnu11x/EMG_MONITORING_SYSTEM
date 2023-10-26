close all;
clear all;
clc;
Fs = 1000;  % sampling Freq
TS = 1/Fs;  % Time period
%% Import data

% signal1
data1 = importdata("D:\Project\EMG\signal\NORMAL\2Nsen.txt"); % import signal
signal1 = data1.data(:,3); % import channel 3 data
l1 = length(signal1); 
ts1 = 0:TS:(l1-1)/1000; % time 

% signal2
data2 = importdata("D:\Project\EMG\signal\NORMAL\3Nsen.txt");
signal2 = data2.data(:,3);
l2 = length(signal2);
ts2 = 0:TS:(l2-1)/1000;

% signal2
data3 = importdata("D:\Project\EMG\signal\NORMAL\5Nsen.txt");
signal3 = data3.data(:,3);
l3 = length(signal3);
ts3 = 0:TS:(l3-1)/1000;

% signal 4
data4 = importdata("D:\Project\EMG\signal\NORMAL\8Nsen.txt");
signal4 = data4.data(:,3);
l3 = length(signal4);
ts4 = 0:TS:(l3-1)/1000;

% signal 5
data5 = importdata("D:\Project\EMG\signal\NORMAL\9Nsen.txt");
signal5 = data5.data(:,3);
l3 = length(signal5);
ts5 = 0:TS:(l3-1)/1000;

% signal 6
data6 = importdata("D:\Project\EMG\signal\NORMAL\10Nsen.txt");
signal6 = data6.data(:,3);
l3 = length(signal6);
ts6 = 0:TS:(l3-1)/1000;

%% Butterworth filter

fnyq = Fs/2;  % Nyquist frequency
fcutlow = 400;
fcuthigh = 30;

% 4th order butterworth bandpass filter
[a,b] = butter(4,[fcuthigh,fcutlow]/fnyq,"bandpass");
signal1 = filtfilt(a,b,signal1);
signal2 = filtfilt(a,b,signal2);
signal3 = filtfilt(a,b,signal3);
signal4 = filtfilt(a,b,signal4);
signal5 = filtfilt(a,b,signal5);
signal6 = filtfilt(a,b,signal6);
 


%% Full wave Rectification

rct_signal1 = abs(signal1(:,1));
rct_signal2 = abs(signal2(:,1));
rct_signal3 = abs(signal3(:,1));
rct_signal4 = abs(signal4(:,1));
rct_signal5 = abs(signal5(:,1));
rct_signal6 = abs(signal6(:,1));

%% Root Mean Square
win_size = 50; % Window Size

rms_value1 = sqrt(movmean(rct_signal1.^2,win_size)); % Root Mean Square EMG
rms_value2 = sqrt(movmean(rct_signal2.^2,win_size));
rms_value3 = sqrt(movmean(rct_signal3.^2,win_size));
rms_value4 = sqrt(movmean(rct_signal4.^2,win_size));
rms_value5 = sqrt(movmean(rct_signal5.^2,win_size));
rms_value6 = sqrt(movmean(rct_signal6.^2,win_size));

%% RMS EMG SIGNALS
figure(1);
subplot(6,1,1),plot(ts1,rms_value1),xlabel("Time"),ylabel("Voltage");
subplot(6,1,2),plot(ts2,rms_value2),xlabel("Time"),ylabel("Voltage");
subplot(6,1,3),plot(ts3,rms_value3),xlabel("Time"),ylabel("Voltage");
subplot(6,1,4),plot(ts4,rms_value4),xlabel("Time"),ylabel("Voltage");
subplot(6,1,5),plot(ts5,rms_value5),xlabel("Time"),ylabel("Voltage");
subplot(6,1,6),plot(ts6,rms_value6),xlabel("Time"),ylabel("Voltage");

%% Average Amplitude
avgemg1 = mean(rms_value1);
avgemg2 = mean(rms_value2);
avgemg3 = mean(rms_value3);
avgemg4 = mean(rms_value4);
avgemg5 = mean(rms_value5);
avgemg6 = mean(rms_value5);

%% Reference mean value

 ref_val = (avgemg1 + avgemg2 + avgemg3 + avgemg4 + avgemg5 + avgemg6)/6;
 fprintf("Reference Amplitude: %.5g\n",ref_val);

 %% sEMG Comparision

 % signal1
ab_data1 = importdata("D:\Project\EMG\signal\ABNORMAL\1Asen.txt"); % import signal
ab_signal1 = ab_data1.data(:,3); % import channel 3 data
l = length(ab_signal1); 
ab_ts1 = 0:TS:(l-1)/1000; % time 

% signal2
ab_data2 = importdata("D:\Project\EMG\signal\ABNORMAL\2Asen.txt");
ab_signal2 = ab_data2.data(:,3);
l = length(ab_signal2);
ab_ts2 = 0:TS:(l-1)/1000;

% signal3
ab_data3 = importdata("D:\Project\EMG\signal\ABNORMAL\5Asen.txt");
ab_signal3 = ab_data3.data(:,3);
l = length(ab_signal3);
ab_ts3 = 0:TS:(l-1)/1000;

% signal 4
ab_data4 = importdata("D:\Project\EMG\signal\ABNORMAL\10Nsen.txt");
ab_signal4 = ab_data4.data(:,3);
l = length(ab_signal4);
ab_ts4 = 0:TS:(l-1)/1000;

% signal 5
ab_data5 = importdata("D:\Project\EMG\signal\ABNORMAL\8Nsen.txt");
ab_signal5 = ab_data5.data(:,3);
l = length(ab_signal5);
ab_ts5 = 0:TS:(l-1)/1000;

% signal 6
ab_data6 = importdata("D:\Project\EMG\signal\ABNORMAL\4Nsen.txt");
ab_signal6 = ab_data6.data(:,3);
l = length(ab_signal6);
ab_ts6 = 0:TS:(l-1)/1000;

%% Butterworth filter

fnyq = Fs/2;  % Nyquist frequency
fcutlow = 400;
fcuthigh = 30;

% 4th order butterworth bandpass filter
[a,b] = butter(4,[fcuthigh,fcutlow]/fnyq,"bandpass");
ab_signal1 = filtfilt(a,b,ab_signal1);
ab_signal2 = filtfilt(a,b,ab_signal2);
ab_signal3 = filtfilt(a,b,ab_signal3);
ab_signal4 = filtfilt(a,b,ab_signal4);
ab_signal5 = filtfilt(a,b,ab_signal5);
ab_signal6 = filtfilt(a,b,ab_signal6);

%% Full wave Rectification

ab_rct_signal1 = abs(ab_signal1(:,1));
ab_rct_signal2 = abs(ab_signal2(:,1));
ab_rct_signal3 = abs(ab_signal3(:,1));
ab_rct_signal4 = abs(ab_signal4(:,1));
ab_rct_signal5 = abs(ab_signal5(:,1));
ab_rct_signal6 = abs(ab_signal6(:,1));

%% Root Mean Square
win_size = 50; % Window Size

ab_rms_value1 = sqrt(movmean(ab_rct_signal1.^2,win_size)); % Root Mean Square EMG
ab_rms_value2 = sqrt(movmean(ab_rct_signal2.^2,win_size));
ab_rms_value3 = sqrt(movmean(ab_rct_signal3.^2,win_size));
ab_rms_value4 = sqrt(movmean(ab_rct_signal4.^2,win_size));
ab_rms_value5 = sqrt(movmean(ab_rct_signal5.^2,win_size));
ab_rms_value6 = sqrt(movmean(ab_rct_signal6.^2,win_size));

%% RMS EMG SIGNALS
figure(2);
subplot(6,1,1),plot(ab_ts1,ab_rms_value1),xlabel("Time"),ylabel("Voltage");
subplot(6,1,2),plot(ab_ts2,ab_rms_value2),xlabel("Time"),ylabel("Voltage");
subplot(6,1,3),plot(ab_ts3,ab_rms_value3),xlabel("Time"),ylabel("Voltage");
subplot(6,1,4),plot(ab_ts4,ab_rms_value4),xlabel("Time"),ylabel("Voltage");
subplot(6,1,5),plot(ab_ts5,ab_rms_value5),xlabel("Time"),ylabel("Voltage");
subplot(6,1,6),plot(ab_ts6,ab_rms_value6),xlabel("Time"),ylabel("Voltage");

%% average Amplitude
ab_avgemg1 = mean(ab_rms_value1);
ab_avgemg2 = mean(ab_rms_value2);
ab_avgemg3 = mean(ab_rms_value3);
ab_avgemg4 = mean(ab_rms_value4);
ab_avgemg5 = mean(ab_rms_value5);
ab_avgemg6 = mean(ab_rms_value5);

%% Comparison with ref_val

result1 = ( (abs (ref_val - ab_avgemg1) ) / ref_val) * 100;
result2 = ( (abs (ref_val - ab_avgemg2) ) / ref_val) * 100;
result3 = ( (abs (ref_val - ab_avgemg3) ) / ref_val) * 100;
result4 = ( (abs (ref_val - ab_avgemg4) ) / ref_val) * 100;
result5 = ( (abs (ref_val - ab_avgemg5) ) / ref_val) * 100;
result6 = ( (abs (ref_val - ab_avgemg6) ) / ref_val) * 100;

fprintf("Result of 1 sEMG: %.3g%% \n",result1);
fprintf("Result of 2 sEMG: %.3g%% \n",result2);
fprintf("Result of 3 sEMG: %.3g%%\n",result3);
fprintf("Result of 4 sEMG: %.3g%%\n",result4);
fprintf("Result of 5 sEMG: %.3g%%\n",result5);
fprintf("Result of 6 sEMG: %.3g%%\n",result6);










