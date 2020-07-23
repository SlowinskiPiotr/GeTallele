function results=fit_vpr_in_segments(data,results_in,synth_VAF_vpr05,synth_VAF_cdf,bin_edges)
% This function takes a whole VAF dataset, and a segmentation of the the 
% region of interest generatd with find_segments_and_fit_vpr using a 
% different matached sequencing data. It estimates vprs values in these
% pre-specified segments and uses the Kolmogorov-Smirnov test to compare
% the observed VAF samples in the segments with a synthethic VAF sample
% with vpr=0.5.
%
% Call function: results=fit_vpr_in_segments(data,results_in,synth_VAF_vpr05,synth_VAF_cdf,bin_edges)
%
%      Input: 
%            data - 3xn matrix with n equal to the number of points in the 
%                   dataset, The columns are: 
%                   1. chromosome number, 
%                   2. base-pair location, 
%                   3. VAF value 
%            results_in  - structure with a segmentation of a region of
%                          interest; output of the find_segments_and_fit_vpr function 
%            synth_VAF_vpr05 - synthetic VAF sample with vpr=0.5; output of
%                              the generate_VAF_cdfs function
%            synth_VAF_cdf - matrix with CDFs of synthetic VAF sample used 
%                            to estimate vprs values in the segments;
%                            output of the generate_VAF_cdfs function
%            bin_edges - bin edges used to generate the CDFs 
%      Output:
%            results - a structure with fields:
%               results.data - a matrix with the analysed region of
%                              interest. The columns are: 
%                               1. chromosome number, 
%                               2. base-pair location, 
%                               3. VAF value 
%               results.fitted_vprs - vector of values of the vprs estimated in each
%                                     segment  
%               results.segment_edges - vector with indices of the first
%                                       element of the segment; contains
%                                       also the last index of the region
%                                       of intrest. The indices of the 
%                                       starting points of the segments are given as:
%                                       segment_edges(1:end-1), 
%                                       the end points of the segments are given as:
%                                       [segment_edges(2:end-1)-1 segment_edges(end)].
%               results.pv_vpr_neq_05 - p-values of the Kolmogorov-Smirnov
%                                       test that compares the VAF sample in the segment with
%                                       synthetic VAF sample with vpr=0.5
%
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------


% chromosome that contains the region of inteerest; value read from the output of the  provided segmentation
chr=results_in.data(1,1); 

%find indices that correspond to the chromosome that contains the region of interest
idx_chr=find(data(:,1)==chr); 

%quick test if a chromosome was in consecutive rows (if the data is ok)
if sum(diff(idx_chr)>1)>0
    disp(['warning: something is wrong with the chromosome: ' num2str(chr)])
end

% data on the chromosome of interest
data_chr=data(idx_chr,:); % part of the original data that contains only the chromosome of interest


% data in the range of interest
start_pos=results_in.data(1,2); % the base pair coordinate of the 1st data point in the region of interest
end_pos=results_in.data(end,2); % the base pair coordinate of the last data point in the region of interest

idx_data_start=find(data_chr(:,2)>=start_pos,1,'first'); %index of the first sample in the region of intrest
idx_data_end=find(data_chr(:,2)<=end_pos,1,'last'); %index of the last sample in the region of intrest

% in case there are no VAFs in the range of interest or data is corrupted assign NaN to all the
% outputs 
if isempty(idx_data_start) || isempty(idx_data_end) || data_chr(idx_data_end,1)<data_chr(idx_data_start,1)
    % the last condition can happen if data is corrupted
    results.data=NaN;
    results.fitted_vprs=NaN;
    results.pv_vpr_neq_05=NaN;
    results.segment_edges=NaN;
    return
end

% data in the range of interest
data_in_range=data_chr(idx_data_start:idx_data_end,:);

%%%%%%%%%% SEGMENTS FROM THE PRESPECIFIED SEGMENTATION %%%%%%%%%%
% eedges from the input
idx_segment_all_edges=results_in.segment_edges;
% indices of start and end points of the segments
idx_segment_start=idx_segment_all_edges(1:end-1);
idx_segment_end=[idx_segment_all_edges(2:end-1)-1 idx_segment_all_edges(end)];

nb_of_segments=numel(idx_segment_start);

%%%%%%%%%% FIT VPRS AND COMPUTE SEGMENT PROPERTIES %%%%%%%%%%
% we start by 'folding' the VAF values into 0.5 to 1 range because VAF values are generally symmetric. 
data_VAFs=abs(data_in_range(:,3)-0.5)+0.5;

for i=nb_of_segments:-1:1 % loop runs backward to automatically pre-allocate memory for the variables
    
    % indices of all the data points in the segment
    idcs_segment=idx_segment_start(i):idx_segment_end(i);
    
    % estimation of the vprs in the segments
    fitted_vprs(i)=fit_vpr(data_VAFs(idcs_segment),synth_VAF_cdf,bin_edges);
    
    %Kolmogorov-Smirnov test that compares the VAF sample in the segment with synthetic VAF sample with vpr=0.5
    [~,pv_vpr_neq_05(i)]=kstest2(data_VAFs(idcs_segment),abs(synth_VAF_vpr05-0.5)+0.5);
end

% assign output
results.data=data_in_range;
results.fitted_vprs=fitted_vprs;
results.pv_vpr_neq_05=pv_vpr_neq_05;
results.segment_edges=idx_segment_all_edges;
end

%%%%%%%%% FIT VPRS FUNCTIONS %%%%%%%%%
function vpr=fit_vpr(data_in_segment,synth_CDFs,bin_edges)

% find bin counts to generate CDFs, implemented using histc because it is much faster than histcounts
% implemented using sparse because for Faray sequences quite a lot of bins will be empty  
Nx=sparse(histc(data_in_segment,bin_edges));

% normalised cumulative sum of bin counts is the CDFs
csx=cumsum(Nx/sum(Nx))';

% replicate the CDF of the sample to take advantage of vectorisation (to
% avoid for loop)
csx_mat=repmat(csx,51,1);

% compute earth mover's distance (EMD) between the observed CDF and the 51 synthetic
% CDFs for vpr values from 0.5 to 1 (with 0.01 step)
emd_v=sum(abs(csx_mat-synth_CDFs),2);

% find the best fit, the smallest value of the EMD
[~,i_vpr]=min(emd_v);
% and convert it into vpr value
vpr=(i_vpr+49)/100;
end