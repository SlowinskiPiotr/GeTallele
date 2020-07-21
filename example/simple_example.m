%% generate farey  bins
[bin_edges, bin_centres]=farey_bins(1000);
bin_edges=bin_edges(152097:end);

%%
example_data_Tex=importdata('example_data_Tex.tsv');
example_data_Nex=importdata('example_data_Nex.tsv');
load /Users/pms210/Dropbox/'GeTallele drafts'/data/chromosomes_lengths.mat
totalreads_Tex=example_data_Tex(:,3)+example_data_Tex(:,4);
totalreads_Nex=example_data_Nex(:,3)+example_data_Nex(:,4);

tic,
[generated_VAF_cdf_Tex,generated_VAF_samples_Tex]=generate_VAF_cdfs(totalreads_Tex,bin_edges,10000,1);
toc,
tic,
[generated_VAF_cdf_Nex,generated_VAF_samples_Nex]=generate_VAF_cdfs(totalreads_Nex,bin_edges,10000,0.9);
toc,

%%
data_VAF_Tex=convert_reads_to_VAF(example_data_Tex);
data_VAF_Nex=convert_reads_to_VAF(example_data_Tex);

chr_results=struct([]);

VAF_vpr05_Tex=generated_VAF_samples_Tex(1).smpl;
VAF_vpr05_Nex=generated_VAF_samples_Nex(1).smpl;

VAF_cdf_Tex=generated_VAF_cdf_Tex;
VAF_cdf_Nex=generated_VAF_cdf_Nex;

for chr_idx=1:22
    tic,
    resultsTex=find_segments_and_fit_vpr(data_VAF_Tex,chr_idx,1,chromosomes_lengths(chr_idx),0.2,10,VAF_vpr05_Tex,VAF_cdf_Tex,bin_edges);
    resultsNex=fit_vpr_in_segments(data_VAF_Nex,resultsTex,VAF_vpr05_Nex,VAF_cdf_Nex,bin_edges);
    toc,
    chr_results(chr_idx).Tex=resultsTex;
    chr_results(chr_idx).Nex=resultsNex;
end

%%
colr=lines(7);
colr=colr([1 6 2 3],:);

for chr_idx=1:22
    plot_VAF_and_vpr(chr_results(chr_idx).Tex,colr(3,:))
    hold one
    plot_VAF_and_vpr(chr_results(chr_idx).Nex,colr(1,:))
    hold off
    pause,
end

%%
% generate all possible mixtures
tic,
[all_vprs,populations,alleles]=gen_all_possible_vprs(4,4,0.01);
toc,
save -v7.3 vprs_p4e4 all_vprs populations alleles

VBP_est=NaN(1,72);
%
parfor dataset=1:72    
    [VBP_est(dataset),vbp(dataset).data,~,~]=VBP_estimator3(all_out(dataset).chr_results,all_vprs,populations);
end
%%
for i=1:72
    purity_estim(i,:)=all_data(i).purity;
end

