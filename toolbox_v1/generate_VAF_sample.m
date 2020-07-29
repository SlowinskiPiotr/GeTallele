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
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------


nb_reads=numel(totalreads); %number of data points - value necessery for resampling
% N1 of synthetics sample will be generated with vpr and N2 will be generated
% with 1-vpr so there is the same number of VAF samples >0.5 and <0.5.
% In case of odd values N1=N2+1;
N1=ceil(nb_of_datapoints/2);
N2=floor(nb_of_datapoints/2); 

% generate VAF samples for vpr
p=vpr*ones(1,N1);
% resample with replacment`
new_count_idx=randi(nb_reads,1,N1); 
new_count=totalreads(new_count_idx);
% generate random numbers from binomial distribution
% new_count(:) varies between y values, p is fixed (it is a vector to take
% advantage of vectrisation)
y = binornd(new_count(:),p(:));
% convert into VAF values
surr_reads1=y(:)./new_count(:);

% generate VAF samples for 1-vpr
p=1-vpr;
p=p*ones(1,N2);
new_count_idx=randi(nb_reads,1,N2);
new_count=totalreads(new_count_idx);
y = binornd(new_count(:),p(:));
surr_reads2=y(:)./new_count(:);

% combine both sets of VAF values into a single sample
VAF_sample=[surr_reads1; surr_reads2];