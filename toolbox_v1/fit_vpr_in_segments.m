function results=fit_vpr_in_segments(data,results_in,VAF_vpr05,generated_VAF_cdf,bin_edges)
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
%
% input:
%   data - matrix with chromosome, base pair positions and with VAF values 
%          for analysis from RNA2DNAlign
%
%   chr,start_pos,end_pos,layer - coordinates of the segment in the data:
%                                   chr=chosen chromosome,
%                                   start_pos=first base pair,
%                                   end_pos  =last base pair,
%                                   layer    =signal for windows generation                                  
%                                             (convention: Nex=1, Ntr=2, Tex=3, Ttr=4)
%
%
% output:
% results - structure with fields:

%
% If you have any questions please contact: p.m.slowinski@exeter.ac.uk

% START 
%here we find indices that correspond to the chromosome that contains the segment

chr=results_in.data(1,1);
idx_chr=find(data(:,1)==chr);

if sum(diff(idx_chr)>1)>0
    disp(['warning: something is wrong with the chromosome: ' num2str(chr)])
    %quick check if a chromosome was in consecutive rows
end

% data on the chromosome of interest
data_chr=data(idx_chr,:); % data_chr is a matrix with base-pairs and all layers

% data in the range of interest
start_pos=results_in.data(1,2);
end_pos=results_in.data(end,2);

idx_data_start=find(data_chr(:,2)>=start_pos,1,'first');
idx_data_end=find(data_chr(:,2)<=end_pos,1,'last');

% in case there are no VAFs in the range of interest
if isempty(idx_data_start) || isempty(idx_data_end) || data_chr(idx_data_end,1)<data_chr(idx_data_start,1)
    results.data=NaN;
    results.total_length=NaN;
    results.total_bp_length=NaN;
    results.fitted_vprs=NaN;
    results.seg_dp_length=NaN;
    results.seg_bp_length=NaN;
    results.pv_vpr_neq_05=NaN;
    results.segment_edges=NaN;
    return
end

data_in_range=data_chr(idx_data_start:idx_data_end,:);
nb_of_pts_data_in_range=size(data_in_range,1);

%%%%%%%%%% SEGMENTS FROM RESULTS %%%%%%%%%%
idx_segment_all_edges=results_in.segment_edges;
idx_segment_start=idx_segment_all_edges(1:end-1);
idx_segment_end=[idx_segment_all_edges(2:end-1)-1 idx_segment_all_edges(end)];

nb_of_segments=numel(idx_segment_start);

%%%%%%%%%% FIT VPRS AND COMPUTE SEGMENT PROPERTIES %%%%%%%%%%
data_VAFs=abs(data_in_range(:,3)-0.5)+0.5;

for i=nb_of_segments:-1:1
    idcs_segment=idx_segment_start(i):idx_segment_end(i);
    fitted_vprs(i)=fit_vpr(data_VAFs(idcs_segment),generated_VAF_cdf,bin_edges);
    seg_dp_length(i)=numel(idcs_segment);
    seg_bp_length(i)=data_in_range(idx_segment_end(i),2)-data_in_range(idx_segment_start(i),2);
    [~,pv_vpr_neq_05(i)]=kstest2(data_VAFs(idcs_segment),abs(VAF_vpr05-0.5)+0.5);
end

results.data=data_in_range;
results.total_length=nb_of_pts_data_in_range;
results.total_bp_length=end_pos-start_pos;
results.fitted_vprs=fitted_vprs;
results.seg_dp_length=seg_dp_length;
results.seg_bp_length=seg_bp_length;
results.pv_vpr_neq_05=pv_vpr_neq_05;
results.segment_edges=idx_segment_all_edges;
end

%%%%%%%%% FIT VPRS FUNCTIONS %%%%%%%%%
function vpr=fit_vpr(x1,models,bin_edges)

Nx=sparse(histc(x1,bin_edges));
lx=sum(Nx);
csx=cumsum(Nx/lx)';
csx_mat=repmat(csx,51,1);

emd_v=sum(abs(csx_mat-models),2);

[~,i_vpr]=min(emd_v);
vpr=(i_vpr+49)/100;
end