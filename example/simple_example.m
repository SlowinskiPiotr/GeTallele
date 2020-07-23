%% generate bins based on farey sequences
[bin_edges, bin_centres]=farey_bins(1000);
bin_edges=bin_edges(152097:end);

%% load data and path to the toolbox
addpath(genpath('../toolbox_v1/')) % add path to the toolbox

% data has the following structure:
% column 1 - chromosome
% column 2 - position in base pairs  on the chromosome
% column 3 - variant read count
% column 4 - reference read count

example_data_Tex=importdata('example_data_Tex.tsv'); 

%% compute synthetic VAF distributions for vprs from 0.5 to 1 for Nex, Ntr, Tex and Ttr signals

totalreads_Tex=example_data_Tex(:,3)+example_data_Tex(:,4); % total read count is a sum of variant and reference read counts
nb_of_syth_points=10000; % each synthetic sample will contain 10000 VAF values
max_VAF=1; 
% we will includ all VAF values (it migth be necessery to set max_VAF<1 if the orignal signal has been filtered 
% for example to keep only heterozygouse sites)
tic, 
[synth_VAF_cdf_Tex,synth_VAF_samples_Tex]=generate_VAF_cdfs(totalreads_Tex,bin_edges,nb_of_syth_points,max_VAF);
toc, % computation took ~1 sec on a machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory
% computation took ~10 sec for transcriptome signals

%% segmentation and  estimation of vprs values
data_VAF_Tex=convert_reads_to_VAF(example_data_Tex); %convert reads to VAFs
VAF_vpr05_Tex=synth_VAF_samples_Tex(1).smpl; %for this analysis we only need synthethic VAF sample with vpr=0.5

start_pos = 1;  % to segment whole chromosome it is set to 1.
end_pos = Inf;  % to segment whole chromosome it is set to Inf.
sensitivity = 0.1; % parameter 'MinThreshold' of the findchangepts function.
shortest_segment = 10; % lenght of the shortest segment that can be generated, parameter 'MinDistance' of the findchangepts function

for chr_idx=22:-1:1 % loop runs backward to automatically pre-allocate memory for the variables
    tic, 
    % this function segments the chromosome and estimates vpr in each segment
    resultsTex=find_segments_and_fit_vpr(data_VAF_Tex,chr_idx,start_pos,end_pos,sensitivity,shortest_segment,VAF_vpr05_Tex,synth_VAF_cdf_Tex,bin_edges);
    toc, % computation took upto 1 sec on a machine with: Intel Core i7 (2.2 GHz) processor and 32 GB memory
    
    chr_results(chr_idx).Tex=resultsTex;
end

%% Plot the VAF values with estimated vpr values
colr=lines(7);
colr=colr([1 6 2 3],:); %color scheme as in paper

figure(1)
% separatly for each chromosome
for chr_idx=1:22
    plot_VAF_and_vpr(chr_results(chr_idx).Tex,colr(3,:),0.025)
    pause, % pause after each plot
end

%% generate circos plot
% to install circos on MAC Os use http://circos.ca/tutorials/lessons/configuration/distribution_and_installation/
% and https://groups.google.com/forum/#!topic/circos-data-visualization/qOIjmy4_wnM
% worked in Feb 2020
% brew tap homebrew/science
% brew remove gd
% brew install gd
% brew install cpanminus
% sudo chown "$USER":admin /Library/Perl/5.18 # check your version on Mac OS
% sudo cpanm Config::General Font::TTF::Font Math::Bezier Math::VecStat Readonly Set::IntSpan Text::Format
% sudo cpanm --force GD::Polyline
% brew install circos

% move config files from toolbox/circos/circos_conf_files to
% your_circos_folder/etc e.g. /usr/local/Cellar/circos/0.69-9/libexec/etc
% folder with homebrew cirocs installation

% in the commands below change directory path to one where you want to keep data for generating circos plots
% edit the .conf file if you change names of the files with data

% insert results into a different structure
for chrm=22:-1:1    
chr_Tex(chrm).results=chr_results(chrm).Tex;
end

convert_vprs_to_circos_line(chr_Tex,0.57,0.28,0,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_1l_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(chr_Tex,0.57,0.28,1,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_1l_t.txt') %to plot vpr >0.5

convert_VAF_to_circos_data(data_VAF_Tex,'/usr/local/Cellar/circos/0.69-9/libexec/data/SAMPLE_1l.txt') %to plot VAF

cd ~
cd ../..
cd /usr/local/Cellar/circos/0.69-9/libexec/data % change path to directory with circos distribution

% to change color edit the vpr_VAF_1layer.conf file
system('../bin/circos -conf etc/vpr_VAF_1layer.conf -outputfile circos_example_1l.png'); 

% output files are saved in /usr/local/Cellar/circos/0.69-9/libexec/data/