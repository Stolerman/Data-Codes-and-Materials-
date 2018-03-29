function [filled_data] = Fill_Gaps(data)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



a=find(~isnan(data));   % find index which are not NaN.
r1=randintrlv(a,1);     % Shuffle no-NaN index vector

ss=length(r1);          % Set subsample size 
perm=r1(1:ss);          % take the first 'ss' subsampled indexes

t2= data(perm);      
n=length(data);
D=dct(eye(n,n));
A=D(perm,:);

% Begin CVX routine
cvx_begin
variable x3(n);
minimize(norm(x3,1));
subject to
A*x3 == t2';
cvx_end
% end CVX routine

filled_data=dct(x3);    % set filled_data via discrete cosine trasform of the optimal solution
figure, plot(filled_data,'r.-')
hold on
plot(data,'b.-')
end




