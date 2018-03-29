
function [ INDEX,DATES ] = Index_and_Dates(ref_day,days_backward,days_forward,period_length,num_of_years)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%  

% Load Days vector
load Days.mat

% Number of initial dates in the loop.  
num_shift_interval_dates= days_forward + days_backward+1;

%calculate reference start index
reference_start_index= find(Days==ref_day);

     
%% 1. Setup vectors of  dates and indexes for the shift interval.


% Vector with the index of the shift interval dates
vec_start_index = zeros(num_shift_interval_dates,1);
vec_start_index(days_backward+1)= reference_start_index;
for j=1:days_backward
    vec_start_index(j)= reference_start_index - days_backward+(j-1);
end

for j=1:days_forward
   vec_start_index(j+days_backward+1)= reference_start_index+j;
end

% Vector with index of the end dates for the shift interval
vec_end_index = vec_start_index + period_length-1*ones(num_shift_interval_dates,1);

% Vector with index of both start and end dates of the shift interval.
index_shift_interval = [vec_start_index,vec_end_index];


% shift interval start & end days : DDMMAAAA
shift_interval_dates_years= Days(index_shift_interval);


% shift interval start & end days : DDMM
shift_interval_dates= floor(Days(index_shift_interval)/10000);

% auxiliar vector which will be used to obtain the dates across years.
year_to_sum=zeros(num_shift_interval_dates,1);
year_to_sum(index_shift_interval(:,1)<find(Days==01012003))=ones( length(find(index_shift_interval(:,1)<find(Days==01012003))),1);
year_to_sum(index_shift_interval(:,1)>=find(Days==01012003))=2*ones( length(find(index_shift_interval(:,1)>=find(Days==01012003))),1);

%% 2. Create matrix of dates across years for each simulation (row of initial_days)

INDEX = zeros(num_of_years,2,num_shift_interval_dates);
DATES = zeros(num_of_years,2,num_shift_interval_dates);

for j=1 :num_shift_interval_dates
    DATES(1,:,j) = shift_interval_dates_years(j,:);
    INDEX(1,:,j)= index_shift_interval(j,:);
 for d=2:num_of_years
 DATES(d,1,j)= shift_interval_dates(j,1)*10000 + 2000 + d + year_to_sum(j,1); 
 INDEX(d,1,j)= find(Days==DATES(d,1,j));
 INDEX(d,2,j) = INDEX(d,1,j) + period_length - 1;
 DATES(d,2,j)= Days(INDEX(d,2,j));
 end
 
end


 end
 
 

