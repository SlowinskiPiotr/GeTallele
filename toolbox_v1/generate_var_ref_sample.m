function [var_ref_reads,resampled_totalreads]=generate_var_ref_sample(totalreads,vpr,nb_of_datapoints)
% This function generates a given number (nb_of_datapoints) of varaiant/reference read values
% for a set variant probability V_pr (vpr) based on a vector of provided
% total reads (totalreads). It also returns resampled totalreads values
% used to generate varaiant/reference read values.
%
% Call function: VAF_sample=generate_VAF_sample(totalreads,vpr,nb_of_datapoints)
%
%      Input: 
%            totalreads - vector of total reads (integers)
%            vpr        - chosen value of variant probability
%            nb_of_datapoints - number of VAF values that will be gnerated 
%      Output:
%            var_ref_reads - vector of variant/reference values
%            resampled_totalreads - vector of resampled totalread values

%        
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

nb_reads=numel(totalreads); %number of data points - value necessery for resampling
% N1 of synthetic samples will be generated with vpr and N2 will be generated
% with 1-vpr so there is the same number of reference and variant samples.
% In case of odd values N1=N2+1;
N1=ceil(nb_of_datapoints/2);
N2=floor(nb_of_datapoints/2);

% generate variat samples using vpr
p=vpr*ones(1,N1);
% resample with replacment
new_count_idx=randi(nb_reads,1,N1);
new_count1=totalreads(new_count_idx);
% generate random numbers from binomial distribution
% new_count(:) varies between y values, p is fixed (it is a vector to take
% advantage of vectorisation)
y = binornd(new_count1(:),p(:));
surr_reads1=y(:);

% generate reference samples using 1-vpr
p=(1-vpr)*ones(1,N2);
new_count_idx=randi(nb_reads,1,N2);
new_count2=totalreads(new_count_idx);
y = binornd(new_count2(:),p(:));
surr_reads2=y(:);

% combine both sets of  values into a single sample
var_ref_reads=[surr_reads1; surr_reads2];
resampled_totalreads=[new_count1; new_count2];

% permute the values so the reference and variant reads are mixed together
% (matters practically only for later visualisation)
idx_rand_perm=randperm(nb_of_datapoints);

var_ref_reads=var_ref_reads(idx_rand_perm);
resampled_totalreads=resampled_totalreads(idx_rand_perm);

