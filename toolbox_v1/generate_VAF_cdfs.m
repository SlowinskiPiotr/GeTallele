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
%            nb_of_datapoints - number of VAF values that will be gnerated 
%            max_VAF - cut-off value for the  VAF values in the sample
%            (values bigger than max_VAF and values smaller than 1-max_VAF
%            are removed from the sample)
%      Output:
%            VAF_cdfs - a sparse matrix with VAF cdfs with size:
%            51x(bin_edges) (in histc The last bin will count any values of X that match bin_edges(end))
%            VAF_samples - a matlab structure with all the generated samples
%        
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

for k=100:-1:50
    
    sample=generate_VAF_sample(totalreads,k/100,nb_of_datapoints);
    
    if max_VAF<1
        idx=sample>max_VAF | sample<(1-max_VAF);
        sample(idx)=[];
    end     

    VAF_samples(k-49).smpl=sample;
end

VAF_cdfs=zeros(51,numel(bin_edges));

for z=51:-1:1
    sample=abs(VAF_samples(z).smpl-0.5)+0.5;
    Ni = sparse(histc(sample,bin_edges));
    li=sum(Ni);
    VAF_cdfs(z,:)=cumsum(Ni)/li;
end

