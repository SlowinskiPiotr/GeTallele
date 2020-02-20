dataset=1;
load dataset_BRCA
load all_data

mat_data=dataset_BRCA(dataset).all_modes_mat_Tex;

% change directory to one where you want to keep data for generating circos plots
% edit the .conf file if you change names of the files with data
convert_vprs_to_circos_line(mat_data,1,0.365,0.075,0,'/Users/pms210/circos-0.69-6/data/vpr_Nex_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,2,0.535,0.075,0,'/Users/pms210/circos-0.69-6/data/vpr_Ntr_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,3,0.705,0.075,0,'/Users/pms210/circos-0.69-6/data/vpr_Tex_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,4,0.875,0.075,0,'/Users/pms210/circos-0.69-6/data/vpr_Ttr_b.txt') %to plot vpr <0.5
convert_vprs_to_circos_line(mat_data,1,0.365,0.075,1,'/Users/pms210/circos-0.69-6/data/vpr_Nex_t.txt') %to plot vpr >0.5
convert_vprs_to_circos_line(mat_data,2,0.535,0.075,1,'/Users/pms210/circos-0.69-6/data/vpr_Ntr_t.txt') %to plot vpr >0.5
convert_vprs_to_circos_line(mat_data,3,0.705,0.075,1,'/Users/pms210/circos-0.69-6/data/vpr_Tex_t.txt') %to plot vpr >0.5
convert_vprs_to_circos_line(mat_data,4,0.875,0.075,1,'/Users/pms210/circos-0.69-6/data/vpr_Ttr_t.txt') %to plot vpr >0.5

dmCNA=dataset_BRCA(dataset).CNA(:,[1:3 5]); % [chrm, bp start, bp end, VAF]
convert_CNA_to_circos_data(dmCNA,'/Users/pms210/circos-0.69-6/data/SAMPLE_CNA.txt') %to plot CNA from theta

dmVAF=all_data(dataset).data(:,[1 2 3]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/Users/pms210/circos-0.69-6/data/SAMPLE_Nex.txt') %to plot VAF

dmVAF=all_data(dataset).data(:,[1 2 4]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/Users/pms210/circos-0.69-6/data/SAMPLE_Ntr.txt') %to plot VAF

dmVAF=all_data(dataset).data(:,[1 2 5]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/Users/pms210/circos-0.69-6/data/SAMPLE_Tex.txt') %to plot VAF

dmVAF=all_data(dataset).data(:,[1 2 6]); % [chrm bp VAF]
convert_VAF_to_circos_data(dmVAF,'/Users/pms210/circos-0.69-6/data/SAMPLE_Ttr.txt') %to plot VAF

cd ~
cd ../..
cd /Users/pms210/circos-0.69-6/data % change to directory with circos distribution
% command to call circos on OS X within Matlab
system(['../bin/circos -conf etc/vpr_VAF_4layers.conf -outputfile circos_example_4l.png']); 