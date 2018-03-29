function [ matrix_dengue, matrix_no_dengue] = fun_clim_statistics_rio(selected_data,num_of_years,num_of_vars,ind_dengue,ind_no_dengue, M)
% Explanation : This function compute physical quantities given a region of
% the ref_day X period_length plane. 


%% cleanning matrices of Physical variables for the statistics
matrix_mean_temp         = zeros(num_of_years,size(selected_data,1));
matrix_rate_precip       = zeros(num_of_years,size(selected_data,1));
matrix_mean_rain_peaks   = zeros(num_of_years,size(selected_data,1));



%% Begin Loop across the region in the Ref_day X period_length plane
for j1 =1:size(selected_data,1)

%% Set period length & ref_day for this step of the loop.     
period_length = selected_data(j1,1) ; 
ref_day       = selected_data(j1,2) ;


[INDEX,~]            = Index_and_Dates(ref_day,0,0,period_length,num_of_years);
index_dates          = INDEX(:,:,1);


%% Slicing the data and storing it in Blocks. 
%Cleaning Blocks
Blocks = zeros(period_length,num_of_vars,num_of_years);

% Loop for setup Blocks and evaluate statistics
for j2 = 1 : num_of_years

 %% Build block for ''year'' j2+1
 ind_start     =  index_dates(j2,1);
 ind_end       =  index_dates(j2,2);
 Blocks(:,:,j2) =  M(ind_start:ind_end,:);
          
 
%%  Compute Mean of average temperature in year j2+1
matrix_mean_temp(j2,j1) = mean(Blocks(:,1,j2));
  

%% Compute info about peaks of precipitation in year j2+1
thres_precip = 0;   % =0 means we accept rain peaks of any intensity.
% Find peaks and their localization.
[peaks,locs] = findpeaks(Blocks(:,2,j2),'MinPeakHeight',thres_precip);

% Mean of rain during precipitation events
matrix_mean_rain_peaks(j2,j1) = mean(peaks)/length(peaks);

%Rate of occurrence of precipitation peak events   
   if numel(locs)==0
       matrix_rate_precip (j2,j1)       = 0;
   else
       matrix_rate_precip (j2,j1)       = 1/mean([locs(1);diff(locs)]);
   end
   clear  locs peaks ind_start ind_end
  
                    
end
end


%% Put data in 'num_of_years' matrices of the form  N X 3 , where N = size(selected_data,1). 

clim_info = zeros(size(selected_data,1),3,num_of_years);
for s=1:num_of_years
clim_info(:,:,s) =  [matrix_mean_temp(s,:)',matrix_rate_precip(s,:)',matrix_mean_rain_peaks(s,:)'];
end


matrix_dengue = clim_info(:,:,ind_dengue(1));
for s2=ind_dengue(2:end)
matrix_dengue = [matrix_dengue;clim_info(:,:,s2)];
end


matrix_no_dengue = clim_info(:,:,ind_no_dengue(1));
for s3=ind_no_dengue(2:end)
matrix_no_dengue = [matrix_no_dengue;clim_info(:,:,s3)];
end




end

