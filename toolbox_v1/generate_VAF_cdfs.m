function [VAF_cdfs,VAF_samples]=generate_VAF_cdfs(totalreads,bin_edges,nb_of_datapoints,max_VAF)
% This function generates a set VAF cumulative denisty functions (CDFs) and
% samples for variant probabilites V_pr in range 0.5 to 1 with 0.01 step.
% Based on a provided vector of total reads (totalreads). All the CDFs have
% the same support - are centered at the same values (set as bin_edges) and are
% based on VAF samples that have a set number of datapoint
% (nb_of_datapoints). It is possibleb to exclude VAF values higher than
% max_VAF from the samples.
%
% Call function: [VAF_cdfs,VAF_samples]=generate_VAF_cdfs(totalreads,bin_edges,max_VAF,nb_of_datapoints)
%
%      Input: 
%            totalreads - vector of total reads (integers)
%            bin_edges  - bin_edges to generate CDF
%            nb_of_datapoints - number of VAF values that will be gnerated in each sample 
%            max_VAF - cut-off value for the  VAF values in the sample
%            (values bigger than max_VAF and values smaller than 1-max_VAF
%            are removed from the sample)
%      Output:
%            VAF_cdfs - a sparse matrix with VAF cdfs with size:
%            51x(bin_edges) (in histc The last bin will count any values of X that match bin_edges(end))
%            VAF_samples - a matlab structure with all the generated samples
%
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

% loop over vpr values from 1 to 0.5  
for k=100:-1:50 %loops goes in reverse direction to speed up performance (there is no need to pre-allocate variables)
    
    sample=generate_VAF_sample(totalreads,k/100,nb_of_datapoints); %VAF sample is being generated
    
    % if required remove VAF samples that are higher and lower the some
    % threshold (for example to keep only heterzygouse sites)
    if max_VAF<1
        idx=sample>max_VAF | sample<(1-max_VAF);
        sample(idx)=[];
    end     
    
    % insert samples into output structure
    VAF_samples(k-49).smpl=sample;
end

% preallocate output matrix
VAF_cdfs=zeros(51,numel(bin_edges));

% loop over all generated samples
for z=1:51
    % 'fold' VAF values to be in range from 0.5 to 1/ justified by symmetry of the VAF distribution
    sample=abs(VAF_samples(z).smpl-0.5)+0.5;
    % find bin counts to generate CDFs, implemented using histc because it is much faster than histcounts
    % implemented using sparse because for Faray sequences quite a lot of bins will be empty  
    Ni = sparse(histc(sample,bin_edges)); 
    % normalised cumulative sum of bin counts is the CDFs
    VAF_cdfs(z,:)=cumsum(Ni)/sum(Ni);
end

