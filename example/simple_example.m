%% generate farey  bins
[bin_edges, bin_centres]=farey_bins(1000);
bin_edges=bin_edges(152097:end);

%% load data and path to the toolbox
addpath(genpath('../toolbox_v1/')) % add path to the toolbox
load('chromosomes_lengths.mat') % chromosomes_lengths in base pairs - values taken from https://en.wikipedia.org/wiki/Human_genome

% data structure:
% column 1 - chromosome
% column 2 - position in base pairs  on the chromosome
% column 3 - variant read count
% column 4 - reference read count
% for comparisons matched sequencing sets are necessery. I.e. different signals have to have reads at the same base pair
% postions. Matched signals can be obtained using RNA2DNAlign package
% (https://github.com/HorvathLab/RNA2DNAlign, https://pubmed.ncbi.nlm.nih.gov/27576531/)

example_data_Tex=importdata('example_data_Tex.tsv'); 
example_data_Nex=importdata('example_data_Nex.tsv');

%% compute synthetic VAF distributions for vprs from 0.5 to 1 for Nex and Tex signals
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
data_VAF_Nex=convert_reads_to_VAF(example_data_Nex);

chr_results=struct([]);

VAF_vpr05_Tex=generated_VAF_samples_Tex(1).smpl;
VAF_vpr05_Nex=generated_VAF_samples_Nex(1).smpl;

VAF_cdf_Tex=generated_VAF_cdf_Tex;
VAF_cdf_Nex=generated_VAF_cdf_Nex;

for chr_idx=1:22
    tic,
    resultsTex=find_segments_and_fit_vpr(data_VAF_Tex,chr_idx,1,chromosomes_lengths(chr_idx),0.1,10,VAF_vpr05_Tex,VAF_cdf_Tex,bin_edges);
    resultsNex=fit_vpr_in_segments(data_VAF_Nex,resultsTex,VAF_vpr05_Nex,VAF_cdf_Nex,bin_edges);
    toc,
    chr_results(chr_idx).Tex=resultsTex;
    chr_results(chr_idx).Nex=resultsNex;
end

%%
colr=lines(7);
colr=colr([1 6 2 3],:);

for chr_idx=1:22
    plot_VAF_and_vpr(chr_results(chr_idx).Tex,colr(3,:),0.025)
    hold on
    plot_VAF_and_vpr(chr_results(chr_idx).Nex,colr(1,:),-0.025)
    hold off
    pause,
end

%% make circos plots

%% vprs based purity estimation
% generate all possible mixtures; computation took >25 sec on machine wwith: 
% Intel Core i7 (2.2 GHz) processor and 32 GB memory
%
max_populations=3; %increasing the number of population to 4 will incereas time of computation to more than an 1.5h
max_events=5;
resolution=0.01;

tic,
[all_vprs,populations]=gen_all_possible_vprs(max_populations,max_events,resolution);
toc,

%% 
% estimate vprs based purity VBP; computation took >10 sec on machine with: 
% Intel Core i7 (2.2 GHz) processor and 32 GB memory
% for 4 populations and 5 events the computation time will icrease to > 8 minutes
[VBP_est,data]=VBP_estimator(chr_results,all_vprs,populations);


    