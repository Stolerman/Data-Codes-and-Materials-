clear all;clc ;close all;clc

%% GENERAL DESCRIPTION
% This script produce the SVM heatmaps for Rio de Janeiro. In our paper we
% choose 2 different kernels: Linear and Rbf. Please read the information
% below before running the code

%% I. INFORMATION ABOUT HEATMAP KERNEL
% 1. For linear kernels
% - choose kernel variable as 'linear' at line 31
% - choose kernel_index variable as 1 at  line 32
% - choose option 'parula' in the colormap command at line 189

% 2. For Rbf kernels
% - choose kernel variable as 'rbf' at line 31
% - choose kernel_index variable as 1 at line 32 
% - choose option 'gray' in the colormap command at line 189


%% II.  INFORMATION ABOUT HEATMAP EDITION 
% The output figure in this script MUST be manually modified with
% Matlab's colormap editor, in order to correspond to those pictures of the
% paper.  
% Please procede as follows:
% 1. Open the figure file
% 2. Right click in the colorbar and choose the option 'open colormap editor'
% 3. Erase (Ctrl+X) all collor pins
% 4. Choose both darkest and lightest color pin and scroll both selected pins in order to match them in the desired threshold (in the paper: 0.8 for linear Kernel and 0.95 for Rbf Kernel) 

%% CHOOSE KERNEL OPTIONS
kernel       = 'linear';   % kernel for the SVM training step
kernel_index = 1;          %  kernel_ind  = 1 for linear and kernel_ind   = 2 for RBF.


%%  LOAD CLIM DATA FOR THE ANALYSIS AND SET FIXED PARAMETERS FOR THE SVM LOOP
load climdata_SVM_Rio_de_Janeiro.mat
M         = [taver precip]; % Define climate data Matrix  
clim_axis = [1 2];          % Fixed variable in our paper.
R1        = 5;              % rectangle dimension in the period length (p) axis
R2        = 6;              % rectangle dimension in the reference days (t0) (p) axis
num_iter  = 100;            % total number of cross-validations 
p_val     = 0.8;            % fraction of selected data for training procedure in the cross-validation steps. 



%% CODE STARTS HERE
SVM_matrix = zeros (length(vec_period_length),length(vec_ref_days));  % this will be our heatmap matrix.


num_ref_days = length(vec_ref_days); % define number of t0 values
num_p = length(vec_period_length);   % define number of p values

a_f = floor(num_p/R1);               % total number of rectangles in the p axis for the heatmap
b_f = floor(num_ref_days/R2);        % total number of rectangles in the t0 axis for the heatmap

%% START LOOP FOR SVM SCORES

%% 1. Loop on the perfectly divided area

for a = 1:a_f
  
    %index of bound period lengths
    index_bound_p  = (a-1)*R1+1:a*R1;
    
    for b=1:b_f

    %index of bound days
    index_bound_days = (b-1)*R2+1:b*R2;
    % transform index of a matrix in index of rows and columns.
    [r,c1] = ind2sub([length(index_bound_p),length(index_bound_days)],[1:length(index_bound_days)*length(index_bound_p)]'); 
    
    % Selected data in a matrix. It's ready for the analysis.
    selected_data  = [ vec_period_length(index_bound_p(r))', vec_ref_days(index_bound_days(c1))]; 
    
    % Loop for ref_days and period_lengths : Funcion
    % 'fun_clim_statistics'.Produce data for fishpots
    [ matrix_dengue, matrix_no_dengue] = fun_clim_statistics_rio(selected_data,num_of_years,num_vars,ind_dengue,ind_no_dengue, M);
    
    % LOOP for SVM score on the given rectangle
    mean_score  = SVM_loop(clim_axis,matrix_dengue,matrix_no_dengue,num_iter,p_val,kernel);

    % Fill SVM matrix with SVM scores for the rectangle     
    SVM_matrix(index_bound_p,index_bound_days) = mean_score.*ones(length(index_bound_p),length(index_bound_days));

    
    end
        
 end


clear a b index_bound_p index_bound_days

%% 2. Loop on upper region, fixed p-final-interval
index_bound_p = a_f*R1+1:num_p;
 
for b=1:b_f
     
%index of bound days
index_bound_days = (b-1)*R2+1:b*R2;
    
% transform index of a matrix in index of rows and columns.
[r,c1] = ind2sub([length(index_bound_p),length(index_bound_days)],[1:length(index_bound_days)*length(index_bound_p)]'); 
% Selected data in a matrix. It's ready for the analysis.
selected_data = [ vec_period_length(index_bound_p(r))', vec_ref_days(index_bound_days(c1))]; 
% Loop for ref_days and period_lengths : Funcion 'fun_clim_statistics'
[ matrix_dengue, matrix_no_dengue] = fun_clim_statistics_rio(selected_data,num_of_years,num_vars,ind_dengue,ind_no_dengue, M);
% LOOP for SVM score on the given rectangle
 mean_score  = SVM_loop(clim_axis,matrix_dengue,matrix_no_dengue,num_iter,p_val,kernel);
% Fill SVM matrix with SVM scores for the rectangle       
SVM_matrix(index_bound_p,index_bound_days) = mean_score.*ones(length(index_bound_p),length(index_bound_days));
end



 clear  b index_bound_p index_bound_days
 
 
 
 
 %% 3. Loop on right region, fixed ref_day-final-interval
 
 % index of bound days
index_bound_days = b_f*R2+1:num_ref_days;
 
for a=1:a_f

%index of bound period length
index_bound_p = (a-1)*R1+1:a*R1;
% transform index of a matrix in index of rows and columns.
[r,c1] = ind2sub([length(index_bound_p),length(index_bound_days)],[1:length(index_bound_days)*length(index_bound_p)]'); 
% Selected data in a matrix. It's ready for the analysis.
selected_data = [ vec_period_length(index_bound_p(r))', vec_ref_days(index_bound_days(c1))]; 
% Loop for ref_days and period_lengths : Funcion 'fun_clim_statistics'
[ matrix_dengue, matrix_no_dengue] = fun_clim_statistics_rio(selected_data,num_of_years,num_vars,ind_dengue,ind_no_dengue, M);
% LOOP for SVM score on the given rectangle
mean_score  = SVM_loop(clim_axis,matrix_dengue,matrix_no_dengue,num_iter,p_val,kernel);
% Fill SVM matrix with SVM scores for the rectangle        
SVM_matrix(index_bound_p,index_bound_days) = mean_score.*ones(length(index_bound_p),length(index_bound_days)); 

end




 clear  a index_bound_p index_bound_days
 
 %% 4. Remaining rectangule : up and right region
 
index_bound_p    = a_f*R1+1:num_p;
index_bound_days = b_f*R2+1:num_ref_days;
 
% transform index of a matrix in index of rows and columns.
[r,c1] = ind2sub([length(index_bound_p),length(index_bound_days)],[1:length(index_bound_days)*length(index_bound_p)]'); 
% Selected data in a matrix. It's ready for the analysis.
selected_data = [ vec_period_length(index_bound_p(r))', vec_ref_days(index_bound_days(c1))]; 
% Loop for ref_days and period_lengths : Funcion 'fun_clim_statistics'
[ matrix_dengue, matrix_no_dengue] = fun_clim_statistics_rio(selected_data,num_of_years,num_vars,ind_dengue,ind_no_dengue, M);
% LOOP for SVM score on the given rectangle
mean_score  = SVM_loop(clim_axis,matrix_dengue,matrix_no_dengue,num_iter,p_val,kernel);
% Fill SVM matrix with SVM scores for the rectangle     
SVM_matrix(index_bound_p,index_bound_days) = mean_score.*ones(length(index_bound_p),length(index_bound_days));

 

 %% PLOT Heatmap

fig=figure;
pcolor(SVM_matrix)
ax = gca;
xlabel('t_0','Fontsize',15)
ylabel('p','Fontsize',15)

%% Label axis 
%Label X axis
DATES=vec_ref_days;
num_initial_dates = length(DATES);
ax.XTick = 1:30:num_initial_dates;
ax.XTickLabel = {floor(DATES(1)/10000),floor(DATES(31)/10000),floor(DATES(61)/10000),...
                 floor(DATES(91)/10000),floor(DATES(121)/10000),floor(DATES(151)/10000),...
                 floor(DATES(181)/10000),floor(DATES(211)/10000),floor(DATES(241)/10000),...
                 };
ax.XTickLabelRotation = 45;
% Label Y axis
ax.YTick = 1:5:length(vec_period_length);
ax.YTickLabel = {vec_period_length(1):5:vec_period_length(end)};


%% choose colormap: parula for 'linear' Kernel and 'gray' for Rbf Kernel
colormap(fig,parula)
colorbar       % insert color bar
shading flat   % take off black lines around rectangles
caxis([0.6 1]) % choose axis
 
%% Title and Saving
suptitle('Rio de Janeiro')
title({['SVM score : ' kernel ' kernel']})
saveas(fig,sprintf('Result_Rio_SVM_kernel_%d_p_%f_R1_%d_R2_%d_numiter_%d.fig',...
                             kernel_index,p_val,R1,R2,num_iter))
