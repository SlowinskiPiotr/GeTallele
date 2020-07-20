%% generate all possible mixtures
% tic,
% [all_vprs,populations,alleles]=gen_all_possible_vprs(4,4,0.01);
% toc,
% save -v7.3 vprs_p4e4 all_vprs populations alleles
%%


%% generate farey  bins
[bin_edges, bin_centres]=farey_bins(1000);
bin_edges=bin_edges(152097:end);
%%

% % 1 chromosomes,
% % 2 positions,
% % 3:6 ratio(Nex, Ntr, Tex, Ttr),
% % 7:10 SNVcount(Nex, Ntr, Tex, Ttr),
% % 11:14 Refcount(Nex, Ntr, Tex, Ttr)
for dataset=72:-1:1
    data_from_tsv=all_data(dataset).data;
    
    idx_chr_1_22=data_from_tsv(:,1)<=22;
    totalreads_Tex=data_from_tsv(idx_chr_1_22,9)+data_from_tsv(idx_chr_1_22,13);
    totalreads_Nex=data_from_tsv(idx_chr_1_22,7)+data_from_tsv(idx_chr_1_22,11);
    
    
    tic,
    [all_out(dataset).generated_VAF_cdf_Tex,all_out(dataset).generated_VAF_samples_Tex]=generate_VAF(totalreads_Tex,bin_edges,1);
    toc,
    tic,
    [all_out(dataset).generated_VAF_cdf_Nex,all_out(dataset).generated_VAF_samples_Nex]=generate_VAF(totalreads_Nex,bin_edges,0.9);
    toc,
end
%%
for dataset=72:-1:1
    
    VAF_vpr05_Tex=all_out(dataset).generated_VAF_samples_Tex(1).smpl;
    VAF_vpr05_Nex=all_out(dataset).generated_VAF_samples_Nex(1).smpl;
    
    VAF_cdf_Tex=all_out(dataset).generated_VAF_cdf_Tex;
    VAF_cdf_Nex=all_out(dataset).generated_VAF_cdf_Nex;
    
    data_from_tsv=all_data(dataset).data;
    data_Tex=data_from_tsv(:,[1 2 5]);
    data_Nex=data_from_tsv(:,[1 2 3]);
    
    chr_results=struct([]);

    tic,
    for chr_idx=1:22
        resultsTex=find_segments_and_fit_vpr(data_Tex,chr_idx,1,chromosomes_lengths(chr_idx),0.2,10,VAF_vpr05_Tex,VAF_cdf_Tex,bin_edges);
        resultsNex=fit_vpr_in_segments(data_Nex,resultsTex,VAF_vpr05_Nex,VAF_cdf_Nex,bin_edges);
        chr_results(chr_idx).Tex=resultsTex;
        chr_results(chr_idx).Nex=resultsNex;
    end
    all_out(dataset).chr_results=chr_results;
    toc,
end    
%%
VBP_est=NaN(1,72);
%
parfor dataset=1:72    
    [VBP_est(dataset),vbp(dataset).data,~,~]=VBP_estimator3(all_out(dataset).chr_results,all_vprs,populations);
end
%%
for i=1:72
    purity_estim(i,:)=all_data(i).purity;
end
%%
% colr=lines(7);
% colr=colr([1 6 2 3],:);
%
% for chr_idx=1:22
%     plot_VAF_and_vpr_simple(chr_results(chr_idx).Tex,colr(3,:))
%     hold one
%     plot_VAF_and_vpr_simple(chr_results(chr_idx).Nex,colr(1,:))
%     hold off
%     pause,
% end
