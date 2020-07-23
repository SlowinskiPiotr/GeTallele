function results=find_segments_and_fit_vpr(data,chr,start_pos,end_pos,sensitivity,shortest_segment,synth_VAF_vpr05,synth_VAF_cdf,bin_edges)
% This function takes a whole VAF dataset, finds in it data points that are
% in the region of interest (specified with chromosome number and
% base-pairs coordinates), segments it into chromosomal segments 
% (continuous multi-SNV genomic regions), estimates vprs values in the
% segments and uses the Kolmogorov-Smirnov test to compare the observed
% VAF samples in the segments with a synthethic VAF sample with vpr=0.5.
%
% Call function: results=find_segments_and_fit_vpr(data,chr,start_pos,end_pos,sensitivity,shortest_segment,synth_VAF_vpr05,synth_VAF_cdf,bin_edges)
%
%      Input: 
%            data - 3xn matrix with n equal to the number of points in the 
%                   dataset, The columns are: 
%                   1. chromosome number, 
%                   2. base-pair location, 
%                   3. VAF value 
%            chr  - number of the chromosome that contains the region that 
%                   is meant to be segmented and on which the vprs are to be fitted
%            start_pos - location in base pairs of the begining of the
%                        region of interest. To include whole chromosome
%                        set it to 1.
%            end_pos - location in base pairs of the end of the
%                      region of interest. To include whole chromosome
%                      set it to Inf.
%            sensitivity - parameter 'MinThreshold' of the findchangepts function. Typically value can be set between 0.05-0.2 
%            shortest_segment - the lenght of the shortest segment that can
%                               be generated, parameter 'MinDistance' of
%                               the findchangepts function.
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

%find indices that correspond to the chromosome that contains the region of interest
idx_chr=find(data(:,1)==chr);

%quick test if a chromosome was in consecutive rows (if the data is ok)
if sum(diff(idx_chr)>1)>0
    disp(['warning: something is wrong with the chromosome: ' num2str(chr)])
end

% data on the chromosome of interest
data_chr=data(idx_chr,:); % part of the original data that contains only the chromosome of interest

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
nb_of_pts_data_in_range=size(data_in_range,1);

%%%%%%%%%% GENERATION OF SEGMENT %%%%%%%%%%
% we start by 'folding' the VAF values into 0.5 to 1 range because VAF values are generally symmetric. 
% also findchangepts uses mean to detect changes in the data
data_VAFs=abs(data_in_range(:,3)-0.5)+0.5; 

% findchangepts function:
% Killick, R. and Eckley I.A. changepoint: An R Package for Changepoint Analysis. Journal of Statistical Software, 2014;58(3):1–19. http://www.jstatsoft.org/v58/i03/.
% Killick, R., Fearnhead, P. and Eckley, I.A. Optimal Detection of Changepoints With a Linear Computational Cost. J Am Stat Assoc 2012;107(500):1590-1598.
% Lavielle, M. Using penalized contrasts for the change-point problem. Signal Process 2005;85(8):1501-1510.
iptsT=findchangepts(data_VAFs,'MinThreshold',sensitivity,'MinDistance',shortest_segment);

% edges from the input findchangepts
idx_segment_all_edges=[1 iptsT' nb_of_pts_data_in_range];
% indices of start and end points of the segments
idx_segment_start=idx_segment_all_edges(1:end-1); 
idx_segment_end=[idx_segment_all_edges(2:end-1)-1 idx_segment_all_edges(end)];

nb_of_segments=numel(idx_segment_start); %number of segments

%%%%%%%%%% FIT VPRS AND COMPUTE SEGMENT PROPERTIES %%%%%%%%%%
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