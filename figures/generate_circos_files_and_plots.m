% to install circos on MAC Os use http://circos.ca/tutorials/lessons/configuration/distribution_and_installation/
% and https://groups.google.com/forum/#!topic/circos-data-visualization/qOIjmy4_wnM
% works Feb 2020
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

dataset=13;

mat_data=participant_BRCA(dataset).all_vprs_mat_Tex;

% change directory to one where you want to keep data for generating circos plots
% edit the .conf file if you change names of the files with data
convert_vprs_to_circos_line(mat_data,1,0.365,0.075,0,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Nex_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,2,0.535,0.075,0,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Ntr_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,3,0.705,0.075,0,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Tex_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,4,0.875,0.075,0,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Ttr_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,1,0.365,0.075,1,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Nex_t.txt') %to plot vpr >0.5
convert_vprs_to_circos_line(mat_data,2,0.535,0.075,1,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Ntr_t.txt') %to plot vpr >0.5
convert_vprs_to_circos_line(mat_data,3,0.705,0.075,1,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Tex_t.txt') %to plot vpr >0.5
convert_vprs_to_circos_line(mat_data,4,0.875,0.075,1,'/usr/local/Cellar/circos/0.69-9/libexec/data/vpr_Ttr_t.txt') %to plot vpr >0.5

dmCNA=participant_BRCA(dataset).CNA(:,[1:3 5]); % [chrm, bp start, bp end, VAF]
convert_CNA_to_circos_data(dmCNA,'/usr/local/Cellar/circos/0.69-9/libexec/data/SAMPLE_CNA.txt') %to plot CNA from theta

dmVAF=all_data(dataset).data(:,[1 2 3]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-9/libexec/data/SAMPLE_Nex.txt') %to plot VAF

dmVAF=all_data(dataset).data(:,[1 2 4]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-9/libexec/data/SAMPLE_Ntr.txt') %to plot VAF

dmVAF=all_data(dataset).data(:,[1 2 5]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-9/libexec/data/SAMPLE_Tex.txt') %to plot VAF

dmVAF=all_data(dataset).data(:,[1 2 6]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-9/libexec/data/SAMPLE_Ttr.txt') %to plot VAF

cd ~
cd ../..
cd /usr/local/Cellar/circos/0.69-9/libexec/data % change to directory with circos distribution
% command to call circos on OS X within Matlab
system('../bin/circos -conf etc/vpr_VAF_4layers.conf -outputfile circos_example_4l.png'); 

% output files are saved in /usr/local/Cellar/circos/0.69-9/libexec/data/