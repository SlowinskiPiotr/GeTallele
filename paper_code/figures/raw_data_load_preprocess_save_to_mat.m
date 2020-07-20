clear all
%% load purity and chromosome_lenghts
%0.8571 0.8700 0.9827 0.8500 0.9342
cd data
import_purity_from_xlsx;
save purity purity_BRCA
chromosomes_lengths=importdata('chromosomes_lengths.mat');
clc
%% load CNA
cd CNA_files

filelist=dir;
k=1;

for i=1:numel(filelist)
    tic,
    if filelist(i).isdir==0 && filelist(i).name(end)=='t'
        all_CNA_data(1).id{k}=filelist(i).name(1:12);
        all_CNA_data(k).data = import_CNA_files(filelist(i).name);
        k=k+1;
    end
    toc,
end

cd ..

save all_CNA_data all_CNA_data
%% load VAF
cd RC_known-rs

% columns in .matrix and .matrix_raw
% 1 chromosomes,
% 2 positions,
% 3:6 ratio(Nex, Ntr, Tex, Ttr),
% 7:10 SNVcount(Nex, Ntr, Tex, Ttr),
% 11:14 Refcount(Nex, Ntr, Tex, Ttr)
% .purity(ESTIMATE, ABSOLUTE, LUMP, IHC, CPE)

all_data=[];

for i=72:-1:1 % 16, 17, 46, 60 edited
    if i<10
        fname=['known_rs_00' num2str(i)  'readcounts.tsv'];
    else
        fname=['known_rs_0' num2str(i)  'readcounts.tsv'];
    end
    disp(fname)
    tic
    data=import_from_tsv(fname,purity_BRCA);
    all_data(i).data=data.matrix;
    all_data(i).data_raw=data.matrix_raw;
    all_data(i).sampleID=data.sampleID;
    all_data(i).purity=data.purity;
    all_data(i).our_label=data.our_label;
    toc
end
cd ..

save -v7.3 all_data_published all_data

cd ..
%% find modes in all datasets
%participant_BRCA=[];
bin_edges=farey_bins(1000);
bin_edges=bin_edges(152097:end);

for dataset=21:-1:1
    disp(dataset)
    participant_BRCA(dataset).sampleID=all_data(dataset).sampleID;
    participant_BRCA(dataset).purity=all_data(dataset).purity;

    tic,
    [participant_BRCA(dataset).model_05,...
        participant_BRCA(dataset).models_cdf]=gen_ideal_hist(all_data(dataset).data,bin_edges);
    toc,
    
    tic
    dataloop=all_data(dataset).data;
    participantloop=participant_BRCA(dataset);
    parfor chridx=1:22%22:-1:1
        chrm(chridx).data_and_stats_Nex=...
            find_win_and_vpr(dataloop,chridx,1,chromosomes_lengths(chridx),1,0.2,participantloop.model_05,participantloop.models_cdf,bin_edges);
        % in our convention layer 1 is Nex
        
        chrm(chridx).data_and_stats_Ntr=...
            find_win_and_vpr(dataloop,chridx,1,chromosomes_lengths(chridx),2,0.2,participantloop.model_05,participantloop.models_cdf,bin_edges);
        % in our convention layer 2 is Ntr
        
        chrm(chridx).data_and_stats_Tex=...
            find_win_and_vpr(dataloop,chridx,1,chromosomes_lengths(chridx),3,0.2,participantloop.model_05,participantloop.models_cdf,bin_edges);
        % in our convention layer 3 is Tex
        
        chrm(chridx).data_and_stats_Ttr=...
            find_win_and_vpr(dataloop,chridx,1,chromosomes_lengths(chridx),4,0.2,participantloop.model_05,participantloop.models_cdf,bin_edges);
        % in our convention layer 4 is Ttr
    end
    participant_BRCA(dataset).chrm=chrm;
    toc
end

save -v7.3 participant_BRCA participant_BRCA
%% generate all possible mixtures
tic,
[all_vprs,populations,alleles]=gen_all_possible_vprs(4,4,0.01);
toc,
save -v7.3 vprs_p4e4 all_vprs populations alleles