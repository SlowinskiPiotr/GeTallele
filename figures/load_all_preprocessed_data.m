addpath(genpath('/Users/pms210/Dropbox/GeTallele/toolbox'));
addpath(genpath('/Users/pms210/Dropbox/GeTallele/data'));

chromosomes_lengths=importdata('chromosomes_lengths.mat');
gene_info=importdata('gene_coordinates.mat'); % .mat file generated using import_gene_names_and_positions_from_xlsx.m
load all_data_published.mat; %.mat file generated using load_preprocess_save_to_mat_data_script.m
load participant_BRCA.mat; %.mat file generated using load_preprocess_save_to_mat_data_script.m
load all_CNA_data.mat; %.mat file generated using load_preprocess_save_to_mat_data_script.m
load number_of_datapoints_in_genes.mat % .mat file generated using number_of_datapoints_in_genes.m

for i=1:72
    sampleID=participant_BRCA(i).sampleID;
    fun = @(s)any(strncmpi(s,sampleID,12));
    CNA_idx=cellfun(fun,all_CNA_data(1).id);
    participant_BRCA(i).CNA=all_CNA_data(CNA_idx).data;
    participant_BRCA(i).all_vprs_mat_Tex=all_vpr_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Tex');
    participant_BRCA(i).all_vprs_mat_Ttr=all_vpr_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Ttr');
end

% for vPR based purity figures
load vprs_p4e4; %saved output of function [dictionaries,populations,alleles]=gen_dict(4,4,0.01);
load VBP_all.mat; %.mat file generated...
load all_prp.mat; %.mat file generated...


