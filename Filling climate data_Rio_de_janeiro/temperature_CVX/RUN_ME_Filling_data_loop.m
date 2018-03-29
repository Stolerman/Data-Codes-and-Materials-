clear all; close all; clc

%% Load data   
load temp_to_fill.mat
data_matrix=taver;

%save space for flled data
filled_data=zeros(size(data_matrix(:,1)));

%% Data completion procedure
for i=1
data = data_matrix(:,i)';
[g1,g2,g3]=Find_Gaps(data); % find gaps of missing data 
sig =Fill_Gaps(g3);         % fill gaps with CVX routine
filled_data(:,i)=sig;       
end

%% PLOT
figure, plot(filled_data(:,i),'r.-')
hold on
plot(g3,'b.-')


