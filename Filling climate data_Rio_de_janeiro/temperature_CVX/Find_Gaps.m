function [ length_gaps,gaps,ready_data] = Find_Gaps( data )
%% We find the gaps building the matrix gaps and then insert a middle point if the gap is big.
%% 1. build matrix gaps
a=find(isnan(data));
ready_data=data;
b=diff(a);

%
c=1;
i=1;
x=zeros(10,1);
x(i)=1;
for s1=1:length(b)
if b(s1)==1 
    c=c+1;x(i)=c;
else if b(s1)~=1
        i=i+1; c=1;x(i)=1;
    end
end
end

%
length_gaps=x(x~=0);
id=[1;cumsum(length_gaps)];
gaps=zeros(max(length_gaps),length(length_gaps));

gaps(1:length(a(id(1):id(2))),1)= a(id(1):id(2));

% matrix 'gaps' where each column is a list with the indexes of the holes
for s2=2:length(length_gaps)
gaps(1:length(a(id(s2)+1:id(s2+1))),s2)= a(id(s2)+1:id(s2+1));
end

%% 2. Routine to insert middle point in big gaps (>20  NaNs)
  if isempty(find(length_gaps>20, 1)) % If there is no big gap
     ready_data=data; % keep the data vector the same.

     
else if ~isempty(find(length_gaps>20, 1)) % if there are big gaps
        
% matrix with only big gaps          
big_gaps = gaps(:,length_gaps>20) ;

for s3=1:size(big_gaps,2)

% Before the gap
% pick the last 30 values before the gap    
aux_before = data(big_gaps(1,s3)-30:big_gaps(1,s3)-1); 
% evaluate the mean of the numbers in this period. Not include NaN's! 
b_hole     = mean(aux_before(~isnan(aux_before)));

% After the gap. 
gap_indexes = big_gaps(find(big_gaps(:,s3)),s3); % find gap indexes (get rid of zeros in the gap matrix). 
end_index =  gap_indexes(end); % find last index
%Choose the 30 values after and evaluate the mean in the same way we did before the gap.
aux_after = data(end_index+1:end_index+30);
a_hole =  mean(aux_after(~isnan(aux_after)));

% Evaluate the mean between before-the-gap and after-the-gap values.
m_hole = (b_hole+a_hole)./2;
ready_data(round(sum(gap_indexes([1;end]))./2))= m_hole ;

end



    end
   end

    
 end

