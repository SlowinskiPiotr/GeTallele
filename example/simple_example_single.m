%% generate farey  bins
[bin_edges, bin_centres]=farey_bins(1000);
bin_edges=bin_edges(152097:end);

%% load data and path to the toolbox
addpath(genpath('../toolbox_v1/')) % add path to the toolbox

% data structure:
% column 1 - chromosome
% column 2 - position in base pairs  on the chromosome
% column 3 - variant read count
% column 4 - reference read count
% for comparisons of vprs values matched sequencing sets are necessery. I.e. different signals have to have read counts at the same base pair
% postions. Matched signals can be obtained using RNA2DNAlign package
% (https://github.com/HorvathLab/RNA2DNAlign, https://pubmed.ncbi.nlm.nih.gov/27576531/)

example_data_Tex=importdata('example_data_Tex.tsv'); 

%% compute synthetic VAF distributions for vprs from 0.5 to 1 for Nex, Ntr, Tex and Ttr signals

totalreads_Tex=example_data_Tex(:,3)+example_data_Tex(:,4); 

tic,
[generated_VAF_cdf_Tex,generated_VAF_samples_Tex]=generate_VAF_cdfs(totalreads_Tex,bin_edges,10000,1);
toc,

%%

data_VAF_Tex=convert_reads_to_VAF(example_data_Tex);

VAF_vpr05_Tex=generated_VAF_samples_Tex(1).smpl;


VAF_cdf_Tex=generated_VAF_cdf_Tex;

chr_results=struct([]);

for chr_idx=1:22
    tic,
    resultsTex=find_segments_and_fit_vpr(data_VAF_Tex,chr_idx,1,Inf,0.1,10,VAF_vpr05_Tex,VAF_cdf_Tex,bin_edges);
    toc,
    
    chr_results(chr_idx).Tex=resultsTex;
end

%%
colr=lines(7);
colr=colr([1 6 2 3],:);
figure(1)
for chr_idx=1:22
    plot_VAF_and_vpr(chr_results(chr_idx).Tex,colr(3,:),0.025)
    pause,
end

% to make circos plots see script simple_example_circos.m