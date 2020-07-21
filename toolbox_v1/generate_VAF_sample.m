function VAF_sample=generate_VAF_sample(totalreads,vpr,nb_of_datapoints)
% This function generates a given number (nb_of_datapoints) of VAF values
% with a set variant probability V_pr (vpr) based on a vector of provided total reads (totalreads);
%
% Call function: VAF_sample=generate_VAF_sample(totalreads,vpr,nb_of_datapoints)
%
%      Input: 
%            totalreads - vector of total reads (integers)
%            vpr        - chosen value of variant probability
%            nb_of_datapoints - number of VAF values that will be gnerated 
%      Output:
%            VAF_sample - vector of VAF values
%        
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

nb_reads=numel(totalreads);
N1=ceil(nb_of_datapoints/2);
N2=floor(nb_of_datapoints/2);

p=vpr*ones(1,N1);
new_count_idx=randi(nb_reads,1,N1);
new_count=totalreads(new_count_idx);
y = binornd(new_count(:),p(:));
surr_reads1=y(:)./new_count(:);

p=1-vpr;
p=p*ones(1,N2);
new_count_idx=randi(nb_reads,1,N2);
new_count=totalreads(new_count_idx);
y = binornd(new_count(:),p(:));
surr_reads2=y(:)./new_count(:);

VAF_sample=[surr_reads1; surr_reads2];