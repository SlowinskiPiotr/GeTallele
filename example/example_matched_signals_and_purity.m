% more detailed dscription of the steps can  be found in the
% simple_example.m script

%% add path to the toolbox and load example data 
addpath(genpath('../toolbox_v1/')) % add path to the toolbox

% data has the following structure:
% column 1 - chromosome
% column 2 - base pairs coordinates on the chromosomes
% column 3 - variant read count
% column 4 - reference read count

% for comparisons matched sequencing sets are necessery. I.e. different signals have to have reads at the same base pair
% postions. Matched signals can be obtained using RNA2DNAlign package
% (https://github.com/HorvathLab/RNA2DNAlign, https://pubmed.ncbi.nlm.nih.gov/27576531/)

example_data_Nex=importdata('example_data_Nex.tsv'); %Nex - normal exome
example_data_Ntr=importdata('example_data_Ntr.tsv'); %Ntr - normal transcriptome
example_data_Tex=importdata('example_data_Tex.tsv'); %Tex - tumor exome
example_data_Ttr=importdata('example_data_Ttr.tsv'); %Ttr - tumor transcriptome

%% generate farey  bins
[bin_edges, bin_centres]=farey_bins(1000);
bin_edges=bin_edges(152097:end);

%% compute synthetic VAF distributions for vprs from 0.5 to 1 for Nex, Ntr, Tex and Ttr signals
totalreads_Nex=example_data_Nex(:,3)+example_data_Nex(:,4);
totalreads_Ntr=example_data_Ntr(:,3)+example_data_Ntr(:,4);
totalreads_Tex=example_data_Tex(:,3)+example_data_Tex(:,4); 
totalreads_Ttr=example_data_Ttr(:,3)+example_data_Ttr(:,4); 

tic,
% Nex signal has been filter and only contain VAF values >0.1 and <0.9
% for this reason the max_VAF is set to 0.9
[synth_VAF_cdf_Nex,synth_VAF_samples_Nex]=generate_VAF_cdfs(totalreads_Nex,bin_edges,10000,0.9);
toc, % computation took ~1 sec on a machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory
tic,
[synth_VAF_cdf_Ntr,synth_VAF_samples_Ntr]=generate_VAF_cdfs(totalreads_Ntr,bin_edges,10000,1);
toc, % computation took ~10 sec on a machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory
tic,
[synth_VAF_cdf_Tex,synth_VAF_samples_Tex]=generate_VAF_cdfs(totalreads_Tex,bin_edges,10000,1);
toc, % computation took ~1 sec on a machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory
tic,
[synth_VAF_cdf_Ttr,synth_VAF_samples_Ttr]=generate_VAF_cdfs(totalreads_Ttr,bin_edges,10000,1);
toc, % computation took ~10 sec on a machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory

%%
data_VAF_Nex=convert_reads_to_VAF(example_data_Nex);
data_VAF_Ntr=convert_reads_to_VAF(example_data_Ntr);
data_VAF_Tex=convert_reads_to_VAF(example_data_Tex);
data_VAF_Ttr=convert_reads_to_VAF(example_data_Ttr);

VAF_vpr05_Nex=synth_VAF_samples_Nex(1).smpl;
VAF_vpr05_Ntr=synth_VAF_samples_Ntr(1).smpl;
VAF_vpr05_Tex=synth_VAF_samples_Tex(1).smpl;
VAF_vpr05_Ttr=synth_VAF_samples_Ttr(1).smpl;

VAF_cdf_Nex=synth_VAF_cdf_Nex;
VAF_cdf_Ntr=synth_VAF_cdf_Ntr;
VAF_cdf_Tex=synth_VAF_cdf_Tex;
VAF_cdf_Ttr=synth_VAF_cdf_Ttr;

for chr_idx=22:-1:1
    tic,
    % this function segments the chromosome and estimates vpr in each segment
    resultsTex=find_segments_and_fit_vpr(data_VAF_Tex,chr_idx,1,Inf,0.1,10,VAF_vpr05_Tex,VAF_cdf_Tex,bin_edges);
    % next function takes the segmentation from find_segments_and_fit_vpr
    % and uses it on a different signal to estimate vpr in each segment
    resultsNtr=fit_vpr_in_segments(data_VAF_Ntr,resultsTex,VAF_vpr05_Ntr,VAF_cdf_Ntr,bin_edges);
    resultsNex=fit_vpr_in_segments(data_VAF_Nex,resultsTex,VAF_vpr05_Nex,VAF_cdf_Nex,bin_edges);
    resultsTtr=fit_vpr_in_segments(data_VAF_Ttr,resultsTex,VAF_vpr05_Ttr,VAF_cdf_Ttr,bin_edges);
    toc, % computation took upto 4 secs on a machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory
    
    chr_results(chr_idx).Nex=resultsNex;
    chr_results(chr_idx).Ntr=resultsNtr;
    chr_results(chr_idx).Tex=resultsTex;
    chr_results(chr_idx).Ttr=resultsTtr;
end

%% example of a plotting function
colr=lines(7);
colr=colr([1 6 2 3],:);

for chr_idx=1:22
    subplot(2,1,1)
    plot_VAF_and_vpr(chr_results(chr_idx).Tex,colr(3,:),0.025)
    hold on
    plot_VAF_and_vpr(chr_results(chr_idx).Nex,colr(1,:),-0.025)
    hold off
    
    subplot(2,1,2)
    plot_VAF_and_vpr(chr_results(chr_idx).Tex,colr(3,:),-0.025)
    hold on
    plot_VAF_and_vpr(chr_results(chr_idx).Ttr,colr(4,:),0.025)
    hold off
    pause,
end

% to make circos plots see script simple_example_circos_4l

%% vprs based purity estimation
% generate all possible mixtures
%
max_populations=3; %increasing the number of population to 4 will incereas time of computation to more than an 1.5h
max_events=5;
resolution=0.01;

tic,
[all_vprs,populations,alleles]=gen_all_possible_vprs(max_populations,max_events,resolution);
toc, %computation took >25 sec on machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory

%% estimate vprs based purity VBP

% purity estimation requires the Nex (normal exome) and Tex (tumor exome)
% signals. This is because when we are generating all possible vprs
% (gen_all_possible_vprs.m) we assume that the normal population has
% vpr=0.5. We check this assumption and exclude any segments where vrp of
% Nex is >0.5. To this end we require (p-value of Kolmogorov-Smirnov 
% test that compares the VAF sample in the segment with synthetic VAF
% sample with vpr=0.5 to be <1e-6.

tic,
[VBP_est,lowest_complexity_mix,VBP_dist,all_prp]=VBP_estimator(chr_results,all_vprs,populations);
toc, %computation took >10 sec on machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory
% for 4 populations and 5 events the computation time will icrease to > 8 minutes

%% make a terenary plot for considered populations and events
plot_adm_mix(all_prp,3,5)