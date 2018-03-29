function [ mean_score ] = SVM_loop(clim_axis, matrix_dengue,matrix_no_dengue, num_iter,p,kernel )
%% LOOP for SVM

% initial sum of the scores
sum_score = 0;

% Start Loop
for count = 1:num_iter
 
% Shuffle indexes 
q1 = randperm(size(matrix_dengue,1));
q2 = randperm(size(matrix_no_dengue,1));

% define how many indexes  to be chosen
stop1=floor(p*size(matrix_dengue,1));
stop2=floor(p*size(matrix_no_dengue,1));

xtrain = [matrix_dengue(q1(1:stop1),clim_axis);matrix_no_dengue(q2(1:stop2),clim_axis)];
xtest  = [matrix_dengue(q1(stop1+1:end),clim_axis);matrix_no_dengue(q2(stop2+1:end),clim_axis)];
ctrain = [ones(stop1,1);2*ones(stop2,1)];


% SVM train and Classify
options.MaxIter = 100000;
svm = svmtrain(xtrain,ctrain,'kernel_function',kernel,'Options',options,'ShowPlot',false);
pre = svmclassify(svm,xtest);

% Auxiliar variables
dummy = [ones(size(matrix_dengue,1)-stop1,1);2*ones(size(matrix_no_dengue,1)-stop2,1)];
dummy2= find(abs(pre-dummy));

% Define score
score = 1- length(dummy2)/length(xtest);
sum_score = sum_score + score;

end

% Mean over iterations gives the final score of the rectangle
mean_score = sum_score ./ count;

end

