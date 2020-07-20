%% Figure 3
%  uses subtightplot from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
% data processing
sample=13;

[bin_edges1, bin_centres]=farey_bins(100);

tic,
[~,~,all_models]=gen_ideal_hist(all_data(sample).data,bin_edges);
toc,

layer=3; %Tex
n1=histcounts(all_models(layer).model(1).smpl,bin_edges1);
n3=histcounts(all_models(layer).model(38).smpl,bin_edges1);
n75=histcounts(all_models(layer).model(26).smpl,bin_edges1);
n66=histcounts(all_models(layer).model(17).smpl,bin_edges1);
n95=histcounts(all_models(layer).model(46).smpl,bin_edges1);

%% plotting
colr=lines(7);
colr=colr([1 6 2 3],:);

gap=[0.065 0.05];
m1=[0.05 0.05];
m2=[0.065 0.025];
subtightplot(5,2,10,gap,m1,m2)

n95(n95==0)=NaN;
bar(bin_centres,n95/5000,25,'facecolor','k','edgecolor','none')
hold on
title('Model VAF with v_{pr}= 0.95')

ylabel('Frequency')
set(gca,'xtick',[0 0.5 1],'ytick',[])
xlabel('VAF')
xlim([-0.01 1.01])
lims=axis;
ylim([0 lims(4)*1.05])
plot([0.5 0.5],[0 lims(4)*1.05],'k:')
box on
hold off

subtightplot(5,2,8,gap,m1,m2)
n3(n3==0)=NaN;
bar(bin_centres,n3/5000,25,'facecolor','k','edgecolor','none')
hold on
nTex=histcounts(abs(data_and_stats_chr.vpr_data(2).win_data(:,4)-0.5)+0.5,bin_edges1);
nTex(nTex==0)=NaN;
bar(bin_centres,nTex/nansum(nTex),50,'facecolor',colr(3,:),'edgecolor','none','linewidth',0.5,'facealpha',0.95)
title('Model VAF with v_{pr}= 0.87; {\color[rgb]{0.85,0.325,0.098}VAF_{Tex} with v_{pr}= 0.87}')

ylabel('Frequency')
set(gca,'xtick',[0 0.5 1],'ytick',[],'XTickLabel',[])
xlim([-0.01 1.01])
lims=axis;
ylim([0 lims(4)*1.05])
plot([0.5 0.5],[0 lims(4)*1.05],'k:')
box on
hold off

subtightplot(5,2,6,gap,m1,m2)
n75(n75==0)=NaN;
bar(bin_centres,n75/5000,25,'facecolor','k','edgecolor','none')
hold on
title('Model VAF with v_{pr}= 0.75')

ylabel('Frequency')
set(gca,'xtick',[0 0.5 1],'ytick',[],'XTickLabel',[])
xlim([-0.01 1.01])
lims=axis;
ylim([0 lims(4)*1.05])
plot([0.5 0.5],[0 lims(4)*1.05],'k:')
box on
hold off

subtightplot(5,2,4,gap,m1,m2)
n66(n66==0)=NaN;
bar(bin_centres,n66/5000,25,'facecolor','k','edgecolor','none')
hold on
title('Model VAF with v_{pr}= 0.66')

ylabel('Frequency')
set(gca,'xtick',[0 0.5 1],'ytick',[],'XTickLabel',[])
xlim([-0.01 1.01])
lims=axis;
ylim([0 lims(4)*1.05])
plot([0.5 0.5],[0 lims(4)*1.05],'k:')
box on
hold off

subtightplot(5,2,2,gap,m1,m2)
n1(n1==0)=NaN;
bar(bin_centres,n1/5000,25,'facecolor','k','edgecolor','none')
hold on
nTex=histcounts(abs(data_and_stats_chr.vpr_data(1).win_data(:,4)-0.5)+0.5,bin_edges1);
nTex(nTex==0)=NaN;
bar(bin_centres,nTex/nansum(nTex),50,'facecolor',colr(3,:),'edgecolor','none','linewidth',0.5,'facealpha',0.95)
title('Model VAF with v_{pr}= 0.5; {\color[rgb]{0.85,0.325,0.098}VAF_{Tex} with v_{pr}= 0.5}')

ylabel('Frequency')
set(gca,'xtick',[0 0.5 1],'ytick',[],'XTickLabel',[])
xlim([-0.01 1.01])
lims=axis;
ylim([0 lims(4)*1.05])
plot([0.5 0.5],[0 lims(4)*1.05],'k:')
box on
hold off

subtightplot(5,2,3,gap,m1,m2)
sample=13;
matrix=all_data(sample).data;
layer=3;
l1=layer+6;
l2=l1+4;
idx_c=find(matrix(:,1)<=22);
counts_vector=matrix(idx_c,l1)+matrix(idx_c,l2);

histogram(counts_vector,'facecolor',colr(3,:))
xlabel('total reads i.e., n_{REF}+n_{VAR}')
ylabel('Frequency')
set(gca,'YTickLabel','')
title('{\color[rgb]{0.85,0.325,0.098}Distribution of total reads in the Tex layer}')
xlim([-0.5 2000])

subtightplot(5,2,1,gap,m1,m2)
plot(counts_vector,'.','color',colr(3,:))
ylabel('n_{REF}+n_{VAR}')
xlabel('samples along chromosomes')
title('{\color[rgb]{0.85,0.325,0.098}Total reads in Tex layer (all chromosomes)}')
set(gca,'XTickLabel','','ytick',[1 2000],'YTickLabel',{'1','2e3'})
xlim([0 6780])
ylim([0 2000])

nb_reads=numel(counts_vector);
new_count_idx=randi(nb_reads,1,10000);
new_count=counts_vector(new_count_idx);

subtightplot(5,2,7,gap,m1,m2)
plot(new_count,'.','color','k')
ylabel('total reads')
xlabel('samples')
title('Resampled total reads for model')
set(gca,'XTickLabel','','ytick',[1 2000],'YTickLabel',{'1','2e3'})
axis([0 5000 0 2000])

%%
nb_reads=numel(counts_vector);
        
        p=75/100;
        p=p*ones(1,5000);
        new_count_idx=randi(nb_reads,1,5000);
        new_count=counts_vector(new_count_idx);
        y = binornd(new_count(:),p(:));
        surr1=y(:)./new_count(:);
        
        p=1-75/100;
        p=p*ones(1,5000);
        new_count_idx=randi(nb_reads,1,5000);
        new_count=counts_vector(new_count_idx);
        y = binornd(new_count(:),p(:));
        surr2=y(:)./new_count(:);