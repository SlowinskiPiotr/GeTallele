%% Figure 1
%  uses subtightplot from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
sample=13;
bin_edges=farey_bins(1000);
bin_edges=bin_edges(152097:end);

% for visualisation we select gene with the highest number of VAFs
[pt,gene_idx]=nanmax(nb_pts,[],2);
index_gene=gene_idx(sample);
gene_name=gene_info.gene_names{index_gene};
chrm=gene_info.chr_pos(index_gene,1);

start_pos=gene_info.chr_pos(index_gene,2);
end_pos=gene_info.chr_pos(index_gene,3);

win_start=max([0,start_pos-round(chromosomes_lengths(chrm)*0.02)]);
win_end=min([chromosomes_lengths(chrm),end_pos+round(chromosomes_lengths(chrm)*0.02)]);

tic,
data_and_stats_chr=find_win_and_vpr(all_data(sample).data,chrm,1,chromosomes_lengths(chrm),3,0.2,participant_BRCA(sample).model_05,participant_BRCA(sample).models_cdf,bin_edges);
toc,
tic,
data_and_stats_part=find_win_and_vpr(all_data(sample).data,chrm,win_start,win_end,3,0.2,participant_BRCA(sample).model_05,participant_BRCA(sample).models_cdf,bin_edges);
toc,
tic,
data_and_stats_gene=find_win_and_vpr(all_data(sample).data,chrm,start_pos,end_pos,3,1,participant_BRCA(sample).model_05,participant_BRCA(sample).models_cdf,bin_edges);
toc,

% plotting
colr=lines(7);
colr=colr([1 6 2 3],:);
    
gap=[0.1 0.075];
m1=[0.075 0.075];
m2=[0.05 0.01];

subtightplot(2,2,1,gap,m1,m2)
cla
CNA=participant_BRCA(sample).CNA;
CNA=CNA(CNA(:,1)==1,:);
for k=1:size(CNA,1)
    
    start_bp=CNA(k,2);
    end_bp=CNA(k,3);
    val=CNA(k,5);
    
    patch([start_bp end_bp end_bp start_bp],[zeros(1,2) val*ones(1,2)],colr(1,:),'edgecolor',colr(1,:))
    hold on
end

idx_s=find(participant_BRCA(sample).all_vprs_mat_Tex(:,1)==chrm);
plot([min(participant_BRCA(sample).all_vprs_mat_Tex(idx_s,2)) max(participant_BRCA(sample).all_vprs_mat_Tex(idx_s,3))],[0 0],'k:')
axis([min(participant_BRCA(sample).all_vprs_mat_Tex(idx_s,2)) max(participant_BRCA(sample).all_vprs_mat_Tex(idx_s,3)) -0.5 0.5]);
ylabel('CNA') 
set(gca,'Ytick',[-0.5 0 0.5],'Yticklabel','','Xticklabel','')
box on
hold off


subtightplot(2,2,2,gap,m1,m2)
cla
rectangle('Position',[start_pos,0.005,end_pos-start_pos,0.99],...
          'FaceColor',[.75 .75 .75],'EdgeColor','none','LineWidth',1)
hold on
plot_VAF_and_vpr(data_and_stats_part,[3 4])
box on
title('Region ({\color[rgb]{0.85,0.325,0.098}Tex: 0.5}/{\color[rgb]{0.929,0.694,0.125} Ttr: 0.5})')
ylabel('VAF')
set(gca,'ytick',[0 0.5 1],'YTickLabel',{'0','','1'})
hold off

subtightplot(2,2,3,gap,m1,m2)
cla
rectangle('Position',[win_start,0.005,win_end-win_start,0.99],...
          'FaceColor',[.75 .75 .75],'EdgeColor','none','LineWidth',1)
hold on
plot_VAF_and_vpr(data_and_stats_chr,[3 4])
box on
title(['Chr ' num2str(chrm) ' ({\color[rgb]{0.85,0.325,0.098}Tex: 0.87, 0.5}/{\color[rgb]{0.929,0.694,0.125} Ttr: 0.9, 0.5})'],'interpreter','tex')
xlabel('position in base pairs ')
ylabel('VAF')
set(gca,'ytick',[0 0.5 1],'YTickLabel',{'0','','1'})
hold off

subtightplot(2,2,4,gap,m1,m2)
cla
plot_VAF_and_vpr(data_and_stats_gene,[3 4])
box on
title(['Gene - ' gene_name ' ({\color[rgb]{0.85,0.325,0.098}Tex: 0.5}/{\color[rgb]{0.929,0.694,0.125} Ttr: 0.5})'])
xlabel('position in base pairs')
ylabel('VAF')
set(gca,'ytick',[0 0.5 1],'YTickLabel',{'0','','1'})
hold off

