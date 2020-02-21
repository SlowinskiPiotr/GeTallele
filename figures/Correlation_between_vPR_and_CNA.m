%%
clear all
load participant_BRCA
load all_CNA_data
%%
for i=1:72
    sampleID=participant_BRCA(i).sampleID;
    fun = @(s)any(strncmpi(s,sampleID,12));
    CNA_idx=cellfun(fun,all_CNA_data(1).id);
    participant_BRCA(i).CNA=all_CNA_data(CNA_idx).data;
    participant_BRCA(i).all_modes_mat_Tex=all_vpr_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Tex');
    participant_BRCA(i).all_modes_mat_Ttr=all_vpr_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Ttr');
end

%% supplementary figure - correlations
% find significant correlations
thr=[];
ph=[];
pl=[];
ph7=[];
pl7=[];

for dataset=72:-1:1
    data=participant_BRCA(dataset).all_modes_mat_Tex;
    data(data(:,10)>1e-5 & data(:,6)>=0.58,:)=[];
    data(data(:,11)>1e-5 & data(:,7)>=0.58,:)=[];
%        data(data(:,10)>1e-5,6)=0.5;
%        data(data(:,11)>1e-5,7)=0.5;
%        data(data(:,6)<0.58,6)=0.5;
%        data(data(:,7)<0.58,7)=0.5;
    
    CNAs=data(:,18);
    
    idx_s_min=find(CNAs>=-0.3 & CNAs<=0.3);
    idx_th=find(data(idx_s_min,6)==min(data(idx_s_min,6)));
    
    thr(dataset)=median(data(idx_s_min(idx_th),18));
    thr_l=thr(dataset)-0.05;
    thr_h=thr(dataset)+0.05;
    
    idx_cr_low=(CNAs<thr_h);
    idx_cr_high=(CNAs>thr_l);
        
    [rh(dataset),ph(dataset)]=corr(data(idx_cr_high,18),data(idx_cr_high,6),'type','Pearson');
    [rl(dataset),pl(dataset)]=corr(data(idx_cr_low,18),data(idx_cr_low,6),'type','Pearson');
    [rh7(dataset),ph7(dataset)]=corr(data(idx_cr_high,18),data(idx_cr_high,7),'type','Pearson');
    [rl7(dataset),pl7(dataset)]=corr(data(idx_cr_low,18),data(idx_cr_low,7),'type','Pearson');
    
    idx_cmp_mid=(CNAs>=thr_l & CNAs<=thr_h);
    idx_cmp_low=(CNAs<thr_l);
    idx_cmp_high=(CNAs>thr_h);
    if any(idx_cmp_low)
        p_cmp_l(dataset)=ranksum(data(idx_cmp_low,6),data(idx_cmp_low,7));
        d_cmp_l(dataset)=median(data(idx_cmp_low,7)-data(idx_cmp_low,6));
    end
    if any(idx_cmp_high)
        p_cmp_h(dataset)=ranksum(data(idx_cmp_high,6),data(idx_cmp_high,7));
        d_cmp_h(dataset)=median(data(idx_cmp_high,7)-data(idx_cmp_high,6));
    end
    if any(idx_cmp_mid)
        p_cmp_m(dataset)=ranksum(data(idx_cmp_mid,6),data(idx_cmp_mid,7));
        d_cmp_m(dataset)=median(data(idx_cmp_mid,7)-data(idx_cmp_mid,6));
    end
end

fdrpl=mafdr(pl,'BHFDR',1);
fdrpl7=mafdr(pl7,'BHFDR',1);
fdrph=mafdr(ph,'BHFDR',1);
fdrph7=mafdr(ph7,'BHFDR',1);

sum(fdrph<0.05)
sum(fdrph7<0.05)
sum(fdrpl<0.05)
sum(fdrpl7<0.05)
%%
fdrpl([21 65 17 8])
fdrpl7([21 65 17 8])
fdrph([21 65 17 8])
fdrph7([21 65 17 8])
rl([21 65 17 8])
rl7([21 65 17 8])
rh([21 65 17 8])
rh7([21 65 17 8])
%%
fdrpl([17])
fdrpl7([17])
fdrph([17])
fdrph7([17])
%% plotting 4 correlations
figure(1)
colr=lines(3);
clf

thr=[];
gap=[0.1 0.05];
marg_h=[0.075 0.05];
marg_w=[0.075 0.05];

idx_dataset=[21 65 17 8];
for dataset=1:4
    
    subtightplot(2,2,dataset,gap,marg_h,marg_w)
    data=participant_BRCA(idx_dataset(dataset)).all_modes_mat_Tex;
    data(data(:,10)>1e-5 & data(:,6)>=0.58,:)=[];
    data(data(:,11)>1e-5 & data(:,7)>=0.58,:)=[];
    
    CNAs=data(:,18);  
    idx_s_min=find(CNAs>=-0.3 & CNAs<=0.3);
    idx_th=find(data(idx_s_min,6)==min(data(idx_s_min,6)));
    
    thr(dataset)=median(data(idx_s_min(idx_th),18));
    thr_l=thr(dataset)-0.05;
    thr_h=thr(dataset)+0.05;
    
    idx_s_low=(CNAs<=thr_h);
    idx_s_high=(CNAs>=thr_l);
    idx_s_mid=(CNAs>=thr_l & CNAs<=thr_h);
    
    plot(data(idx_s_low,18),data(idx_s_low,6),'s','color',colr(2,:),'markersize',3,'MarkerFaceColor',colr(2,:));
    hold on
    plot(data(idx_s_high,18),data(idx_s_high,6),'s','color',colr(2,:),'markersize',3,'MarkerFaceColor',colr(2,:));
    plot(data(idx_s_mid,18),data(idx_s_mid,6),'ks','markersize',3,'MarkerFaceColor','k');
    
    plot(data(idx_s_low,18),data(idx_s_low,7),'o','color',colr(3,:),'markersize',3);%,'linewidth',1);
    plot(data(idx_s_high,18),data(idx_s_high,7),'o','color',colr(3,:),'markersize',3);%,'linewidth',1);
    plot(data(idx_s_mid,18),data(idx_s_mid,7),'ko','markersize',3);%,'linewidth',1);
    
    plot([thr_l thr_l],[0.5 1],'k')
    plot([thr_h thr_h],[0.5 1],'k')
    
    plyh=polyfit(data(idx_s_high,18),data(idx_s_high,6),1);
    y=polyval(plyh,sort([data(idx_s_high,18); 3]));
    h(5)=plot(sort([data(idx_s_high,18); 3]),y,'-','color',colr(2,:),'linewidth',1);
    
    plyl=polyfit(data(idx_s_low,18),data(idx_s_low,6),1);
    y=polyval(plyl,sort([-1.5; data(idx_s_low,18)]));
    h(6)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-','color',colr(2,:),'linewidth',1);
    
    plyh7=polyfit(data(idx_s_high,18),data(idx_s_high,7),1);
    y=polyval(plyh7,sort([data(idx_s_high,18); 3]));
    h(2)=plot(sort([data(idx_s_high,18); 3]),y,'-','color',colr(3,:),'linewidth',1);
    
    plyl7=polyfit(data(idx_s_low,18),data(idx_s_low,7),1);
    y=polyval(plyl7,sort([-1.5; data(idx_s_low,18)]));
    if d_cmp_l(dataset)>=0.1 
        h(3)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-','color',colr(3,:),'linewidth',2);
    elseif d_cmp_l(dataset)>=0.05 && d_cmp_l(dataset)<0.1 && fdrpcl(dataset)<0.05
        h(3)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-.','color',colr(3,:),'linewidth',2);
    else
        h(3)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-','color',colr(3,:),'linewidth',1);
    end
    
    xlim([-1 2])
    nsig=0;
  
    if max(data(:,18))>0.3 || min(data(:,18))<-0.3
           
        if sum(idx_s_mid)==sum(idx_s_low) || fdrpl7(idx_dataset(dataset))>0.05 || plyl7(1)>0%rl(dataset)>0%0.05/72  pl_bts7(dataset)>(btsl_sz*0.05/72)
            delete(h(3));
            nsig=nsig-1;
        end

        if sum(idx_s_mid)==sum(idx_s_high) || fdrph7(idx_dataset(dataset))>0.05 || plyh7(1)<0%rh(dataset)<0%0.05/72 ph_bts7(dataset)>(btsh_sz*0.05/72)
            delete(h(2));
            nsig=nsig-1;
        end

        if sum(idx_s_mid)==sum(idx_s_low) || fdrpl(idx_dataset(dataset))>0.05 || plyl(1)>0%rl7(dataset)>0 pl_bts(dataset)>(btsl_sz*0.05/72)
            delete(h(6));
            nsig=nsig-1;
        end
        if sum(idx_s_mid)==sum(idx_s_high) || fdrph(idx_dataset(dataset))>0.05 || plyh(1)<0%rh7(dataset)<0 ph_bts(dataset)>(btsh_sz*0.05/72)
            delete(h(5));
            nsig=nsig-1;
        end
    else
        delete(h([2 3 5 6]))
    end
    
    ylim([0.5 1])
    set(gca,'YTick',0.5:0.1:1,'YTickLabel',[],'XTick',[-0.3 0 0.3],'XTickLabel',[])
    grid on
    %title([num2str(idx_dataset(dataset)) ': #sig ' num2str(nsig) ' ' num2str(participant_BRCA(idx_dataset(dataset)).purity(5),3)])
    if any([1 3]==dataset)
        ylabel('V_{pr}')
        set(gca,'YTick',0.5:0.1:1,'YTickLabel',0.5:0.1:1,'XTick',[-0.3 0 0.3],'XTickLabel',[])
    end
    if dataset==1
        text(-1.3,1.05,'A','fontsize',12,'FontWeight','bold')
    end
    if dataset==2
        text(-1.1,1.05,'B','fontsize',12,'FontWeight','bold')
    end
    if dataset==3
        text(-1.3,1.05,'C','fontsize',12,'FontWeight','bold')
        xlabel('CNA')
        set(gca,'YTick',0.5:0.1:1,'YTickLabel',0.5:0.1:1,'XTick',[-1 -0.3 0 0.3 2],'XTickLabel',{'-1','-0.3','0','0.3','2'})
    end
    if dataset==4
        text(-1.1,1.05,'D','fontsize',12,'FontWeight','bold')
        xlabel('CNA')
        set(gca,'YTick',0.5:0.1:1,'YTickLabel',[],'XTick',[-1 -0.3 0 0.3 2],'XTickLabel',{'-1','-0.3','0','0.3','2'})
    end
    drawnow
    hold off
end
%% plotting

colr=lines(3);
clf
p_bh=[];
p_s=[];
r2=[];
thr=[];
%colr=jet(72);
gap=[0.02 0.02];
marg_h=[0.035 0.02];
marg_w=[0.05 0.01];
% rhbts=[];
% phbts=[];
% rlbts=[];
% plbts=[];
% rhbts7=[];
% phbts7=[];
% rlbts7=[];
% plbts7=[];
btsh_sz=10000;
btsl_sz=10000;
for dataset=1:72
    subtightplot(12,6,dataset,gap,marg_h,marg_w)
    data=participant_BRCA(dataset).all_modes_mat_Tex;
    data(data(:,10)>1e-5 & data(:,6)>=0.58,:)=[];
    data(data(:,11)>1e-5 & data(:,7)>=0.58,:)=[];
    %    data(data(:,10)>1e-5,6)=0.5;
    %    data(data(:,11)>1e-5,7)=0.5;
    %    data(data(:,6)<0.58,6)=0.5;
    %    data(data(:,7)<0.58,7)=0.5;
    
    CNAs=data(:,18);
    
    idx_s_min=find(CNAs>=-0.3 & CNAs<=0.3);
    idx_th=find(data(idx_s_min,6)==min(data(idx_s_min,6)));
    
    thr(dataset)=median(data(idx_s_min(idx_th),18));
    thr_l=thr(dataset)-0.05;
    thr_h=thr(dataset)+0.05;
    
    idx_s_low=(CNAs<thr_h);
    idx_s_high=(CNAs>thr_l);
    idx_s_mid=(CNAs>=thr_l & CNAs<=thr_h);
    
    plot(data(idx_s_low,18),data(idx_s_low,6),'s','color',colr(2,:),'markersize',3,'MarkerFaceColor',colr(2,:));
    hold on
    plot(data(idx_s_high,18),data(idx_s_high,6),'s','color',colr(2,:),'markersize',3,'MarkerFaceColor',colr(2,:));
    plot(data(idx_s_mid,18),data(idx_s_mid,6),'ks','markersize',3,'MarkerFaceColor','k');
    
    if 1
        plot(data(idx_s_low,18),data(idx_s_low,7),'o','color',colr(3,:),'markersize',3);%,'linewidth',1);
        plot(data(idx_s_high,18),data(idx_s_high,7),'o','color',colr(3,:),'markersize',3);%,'linewidth',1);
        plot(data(idx_s_mid,18),data(idx_s_mid,7),'ko','markersize',3);%,'linewidth',1);
    end
    
    plot([thr_l thr_l],[0.5 1],'k')
    %plot([0 0],[0.5 1],'k-.')
    plot([thr_h thr_h],[0.5 1],'k')
    %xlim([-1.5 3])
    
    
    
    [rh(dataset),ph(dataset)]=corr(data(idx_s_high,18),data(idx_s_high,6),'type','spearman');
    [rl(dataset),pl(dataset)]=corr(data(idx_s_low,18),data(idx_s_low,6),'type','spearman');
    [rh7(dataset),ph7(dataset)]=corr(data(idx_s_high,18),data(idx_s_high,7),'type','spearman');
    [rl7(dataset),pl7(dataset)]=corr(data(idx_s_low,18),data(idx_s_low,7),'type','spearman');
    
    [plyh,s]=polyfit(data(idx_s_high,18),data(idx_s_high,6),1);
    y=polyval(plyh,sort([data(idx_s_high,18); 3]));
    h(5)=plot(sort([data(idx_s_high,18); 3]),y,'-','color',colr(2,:));%,'linewidth',1)
    nrh(dataset)=s.normr;
    [plyl,s]=polyfit(data(idx_s_low,18),data(idx_s_low,6),1);
    y=polyval(plyl,sort([-1.5; data(idx_s_low,18)]));
    h(6)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-','color',colr(2,:));%,'linewidth',1)
    nrl(dataset)=s.normr;
    [plyh7,s]=polyfit(data(idx_s_high,18),data(idx_s_high,7),1);
    y=polyval(plyh7,sort([data(idx_s_high,18); 3]));
    h(2)=plot(sort([data(idx_s_high,18); 3]),y,'-','color',colr(3,:));%,'linewidth',1)
    nrh7(dataset)=s.normr;
    [plyl7,s]=polyfit(data(idx_s_low,18),data(idx_s_low,7),1);
    nrl7(dataset)=s.normr;
    y=polyval(plyl7,sort([-1.5; data(idx_s_low,18)]));
    h(3)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-','color',colr(3,:));%,'linewidth',1)
    
%         dbts=data(idx_s_high,18);
%         nbtsh(dataset)=numel(dbts);
%     
%         tic,
%         parfor ibts=1:btsh_sz
%                 [rhbts(dataset,ibts),phbts(dataset,ibts)]=corr(dbts(randperm(nbtsh(dataset))),data(idx_s_high,6),'type','pearson');
%                 [rhbts7(dataset,ibts),phbts7(dataset,ibts)]=corr(dbts(randperm(nbtsh(dataset))),data(idx_s_high,7),'type','pearson');
%                 %[~,s]=polyfit(dbts(randperm(nbtsh(dataset))),data(idx_s_high,6),1);
%                 %nrhbts(dataset,ibts)=s.normr;
%                 %[~,s]=polyfit(dbts(randperm(nbtsh(dataset))),data(idx_s_high,7),1);
%                 %nrhbts7(dataset,ibts)=s.normr;
%         end%
%         toc,
%     
%     
%         dbts=data(idx_s_low,18);
%         nbtsl(dataset)=numel(dbts);
%         tic,
%         parfor ibts=1:btsl_sz
%                 [rlbts(dataset,ibts),plbts(dataset,ibts)]=corr(dbts(randperm(nbtsl(dataset))),data(idx_s_low,6),'type','pearson');
%                 [rlbts7(dataset,ibts),plbts7(dataset,ibts)]=corr(dbts(randperm(nbtsl(dataset))),data(idx_s_low,7),'type','pearson');
% %                 [~,s]=polyfit(dbts(randperm(nbtsl(dataset))),data(idx_s_low,6),1);
% %                 nrlbts(dataset,ibts)=s.normr;
% %                 [~,s]=polyfit(dbts(randperm(nbtsl(dataset))),data(idx_s_low,7),1);
% %                 nrlbts7(dataset,ibts)=s.normr;
%         end
%         toc,
    
%         pl_bts(dataset)=sum(plbts(dataset,:)<=pl(dataset));
%         ph_bts(dataset)=sum(phbts(dataset,:)<=ph(dataset));
%         pl_bts7(dataset)=sum(plbts7(dataset,:)<=pl7(dataset));
%         ph_bts7(dataset)=sum(phbts7(dataset,:)<=ph7(dataset));
    
%     pl_bts(dataset)=sum(nrlbts(dataset,:)<=nrl(dataset));
%     ph_bts(dataset)=sum(nrhbts(dataset,:)<=nrh(dataset));
%     pl_bts7(dataset)=sum(nrlbts7(dataset,:)<=nrl7(dataset));
%     ph_bts7(dataset)=sum(nrhbts7(dataset,:)<=nrh7(dataset));
            xlim([-1.5 3])
    nsig=0;
    if max(data(:,18))>0.3 || min(data(:,18))<-0.3
        %hls=lsline();
        
        if numel(h)>3
            nsig=4;
%             for i=1:3
%                 set(h(i),'linestyle','--');
%             end
%             for i=1:6
%                 set(h(i),'linewidth',1);
%             end
            delete(h(4));
            if sum(idx_s_mid)==sum(idx_s_low) || fdrpl7(dataset)>0.05 || plyl7(1)>0%rl(dataset)>0%0.05/72  pl_bts7(dataset)>(btsl_sz*0.05/72)
                delete(h(3));
                nsig=nsig-1;
            end
            if sum(idx_s_mid)==sum(idx_s_high) || fdrph7(dataset)>0.05 || plyh7(1)<0%rh(dataset)<0%0.05/72 ph_bts7(dataset)>(btsh_sz*0.05/72)
                delete(h(2));
                nsig=nsig-1;
            end
            
            if sum(idx_s_mid)==sum(idx_s_low) || fdrpl(dataset)>0.05 || plyl(1)>0%rl7(dataset)>0 pl_bts(dataset)>(btsl_sz*0.05/72)
                delete(h(6));
                nsig=nsig-1;
            end
            if sum(idx_s_mid)==sum(idx_s_high) || fdrph(dataset)>0.05 || plyh(1)<0%rh7(dataset)<0 ph_bts(dataset)>(btsh_sz*0.05/72)
                delete(h(5));
                nsig=nsig-1;
            end
        else
%             for i=1:3
%                 set(h(i),'linewidth',1);
%             end
            delete(h(1));
            if sum(idx_s_mid)==sum(idx_s_low) || fdrpl(dataset)>0.05 || plyl(1)>0%rl(dataset)>0%0.05/72
                delete(h(3));
            end
            if sum(idx_s_mid)==sum(idx_s_high) || fdrph(dataset)>0.05 || plyh(1)<0%rh(dataset)<0%0.05/72
                delete(h(2));
            end
        end
    else
        delete(h([2 3 5 6]))
    end
    
    ylim([0.5 1])
    set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[],'XTick',0,'XTickLabel',[])
    grid on
    title([num2str(dataset) ': #sig ' num2str(nsig)])
    if any(1:6:67==dataset)
        %ylabel('V_{pr}')
        set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[0.5 0.6 1],'XTick',0,'XTickLabel',[])
    end
    if dataset==67
        %ylabel('V_{pr}')
        %xlabel('CNA')
        set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[0.5 0.6 1],'XTick',[-1.5 0 3],'XTickLabel',{'-1.5','0','3'})
    end
    if any(68:72==dataset)
        %xlabel('CNA')
        set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[],'XTick',[-1.5 0 3],'XTickLabel',{'-1.5','0','3'})
    end
    if dataset==69
        %ylabel('V_{pr}')
        xlabel('CNA')
        %set(gca,'YTick',[0.5 0.58 1],'YTickLabel',[0.5 0.58 1],'XTick',[-1.5 0 3],'XTickLabel',[-1.5 0 3])
    end
    if dataset==31
        ylabel('V_{pr}')
        %xlabel('CNA')
        %set(gca,'YTick',[0.5 0.58 1],'YTickLabel',[0.5 0.58 1],'XTick',[-1.5 0 3],'XTickLabel',[-1.5 0 3])
    end
    %sum(phbts(dataset,:)<0.01))
    %title(participant_BRCA(dataset).purity(5))
    %title([num2str([rl(dataset),rh(dataset)],3) ', ' num2str([pl(dataset),ph(dataset)]*72,3)])
    drawnow
    hold off
end
%%
p_bh(:,1)=pl<0.05;%/72;%pl_bts<10000*0.05;%bonf_holm(pl,0.05)<0.05;%mafdr(p_s(:,1),'BHFDR',1);
p_bh(:,2)=ph<0.05;%/72;%ph_bts<10000*0.05;%ph<0.05/72;%bonf_holm(ph,0.05)<0.05;%bonf_holm(p_s(:,2),0.05)<0.05;%mafdr(p_s(:,2),'BHFDR',1)
sum(p_bh(:,1)==1 & rl'<0)
sum(p_bh(:,2)==1 & rh'>0)
sum(p_bh(:,2)==1 & p_bh(:,1)==1)
%%
for dataset=1:72
    data=participant_BRCA(dataset).all_modes_mat_Tex;
    data(data(:,10)>1e-5 & data(:,6)>=0.58,:)=[];
    data(data(:,11)>1e-5 & data(:,7)>=0.58,:)=[];
    %data(data(:,10)>1e-5,6)=0.5;
    %data(data(:,11)>1e-5,7)=0.5;
    %data(data(:,6)<0.58,6)=0.5;
    %data(data(:,7)<0.58,7)=0.5;
    %data(data(:,6)<0.58,:)=[];
    %data(data(:,7)<0.58,:)=[];
    
    CNAs=data(:,18);
    idx_s=(CNAs<0);
    idx_s_high=(CNAs>0);
    ex_or_tr=6;
    %[r2s,ps]=corr(data(:,15),data(:,6),'rows','pairwise','type','kendall');
    if ~all(idx_s==0)
        [r2i,pi]=corr(data(idx_s,18),data(idx_s,ex_or_tr),'rows','pairwise','type','spearman');
        if pi<0.05/72
            plot(data(idx_s,18),data(idx_s,ex_or_tr),'ro');%,'color',colr(dataset,:))
            %lsline()
            %pause,
            hold on
        else
            plot(data(idx_s,18),data(idx_s,ex_or_tr),'r.')
            hold on
        end
    else
        r2i=NaN;
        pi=NaN;
    end
    
    if ~all(idx_s_high==0)
        [r2ih,pih]=corr(data(idx_s_high,18),data(idx_s_high,ex_or_tr),'rows','pairwise','type','spearman');
        if pih<0.05/72
            plot(data(idx_s_high,18),data(idx_s_high,ex_or_tr),'bx');%,'color',colr(dataset,:))
            %lsline()
            hold on
        else
            plot(data(idx_s_high,18),data(idx_s_high,ex_or_tr),'b.')
            hold on
        end
    else
        r2ih=NaN;
        pih=NaN;
    end
    r2(dataset,1)=r2i;
    r2(dataset,2)=r2ih;
    %r2(dataset,3)=r2s;
    
    p_s(dataset,1)=pi;
    p_s(dataset,2)=pih;
    %p_s(dataset,3)=ps;
    
end
p_bh(:,1)=bonf_holm(p_s(:,1),0.05)<0.05;%p_s(:,1)<0.05/72;%mafdr(p_s(:,1),'BHFDR',1);
p_bh(:,2)=bonf_holm(p_s(:,2),0.05)<0.05;%p_s(:,2)<0.05/72;%bonf_holm(p_s(:,2),0.05)<0.05;%p_s(:,2)<0.05/72;%mafdr(p_s(:,2),'BHFDR',1);bonf_holm(p_s(:,2),0.05)>0.05;%
sum(p_bh(:,1)==1 & r2(:,1)<0)
sum(p_bh(:,2)==1 & r2(:,2)>0)
sum(p_bh(:,2)==1 & p_bh(:,1)==1)
%p_bh(:,3)=mafdr(p_s(:,3),'BHFDR',1);
%sum(p_bh(:,3)<0.05)

%h=lsline;
%for i=1:numel(h) set(h(i),'color',[0.5 0.5 0.5]); end
ylim([0.5 1])
xlabel('CNA')
ylabel('v_{PR,TEX}')
%plot([0.15 0.15],[0.5 1],'k--')
%plot([0 0],[0.5 1],'k--')
title('Significant correlations between CNA and v_{PR,TEX}')
grid on
%% subset of all the data for visualisation
dataset_BRCA=[];
for i=1:72
    dataset_BRCA(i).all_chrm_mat_Tex=participant_BRCA(i).all_modes_mat_Tex;
    dataset_BRCA(i).all_chrm_mat_Ttr=participant_BRCA(i).all_modes_mat_Ttr;
    
    dataset_BRCA(i).CNA=participant_BRCA(i).CNA;
    dataset_BRCA(i).purity=participant_BRCA(i).purity;
    %dataset_BRCA(i).fdr_corr=p_bh(i,:);
    %dataset_BRCA(i).r2_corr=r2(i,:);
end

%% visualisation
%figure,
colr=lines(3);

for dataset=1:72
    if dataset_BRCA(dataset).purity(5)>0
        for chrm=1:22
            if chrm>20
                subtightplot(8,3,chrm+1,[0.02 0.01])
                cla
            else
                subtightplot(8,3,chrm,[0.02 0.01])
                cla
            end
            
            %simple_layer(participant_BRCA(dataset).chrm(chrm).data_and_stats_Ttr,4,1)
            %hold on
            %simple_layer(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex,3,1)
            plot_mode_and_CNA_patch(dataset_BRCA(dataset).all_chrm_mat_Ttr,dataset_BRCA(dataset).CNA,chrm)
            %plot_mode_and_mode_patch(dataset_BRCA(dataset).all_chrm_mat_Tex,dataset_BRCA(dataset).all_chrm_mat_Ttr,chrm)
            
            %idx_chrm=(dataset_BRCA(dataset).all_chrm_mat(:,1)==chrm);
            %d_start=round(nanmean(abs(dataset_BRCA(dataset).all_chrm_mat(idx_chrm,10))),3)*100;
            %d_end=round(nanmean(abs(dataset_BRCA(dataset).all_chrm_mat(idx_chrm,11))),3)*100;
            %         if dataset_BRCA(dataset).chrm(chrm).crl(1,2)<0.05 && dataset_BRCA(dataset).chrm(chrm).crl(2,2)<0.05
            %         title(['Chr: ', num2str(chrm)...
            %                '; score: ' num2str(dataset_BRCA(dataset).chrm(chrm).scr)...
            %                '; corr 1: ' num2str(dataset_BRCA(dataset).chrm(chrm).crl(1,:))...
            %                '; corr 2: ' num2str(dataset_BRCA(dataset).chrm(chrm).crl(2,:))]);
            %         elseif dataset_BRCA(dataset).chrm(chrm).crl(1,2)<0.05
            %             title(['Chr: ', num2str(chrm)...
            %                '; score: ' num2str(dataset_BRCA(dataset).chrm(chrm).scr)...
            %                '; corr 1: ' num2str(dataset_BRCA(dataset).chrm(chrm).crl(1,:))...
            %                '; corr 2: ns' ]);
            %         elseif  dataset_BRCA(dataset).chrm(chrm).crl(2,2)<0.05
            %             title(['Chr: ', num2str(chrm)...
            %                '; score: ' num2str(dataset_BRCA(dataset).chrm(chrm).scr)...
            %                '; corr 1: ns'...
            %                '; corr 2: ' num2str(dataset_BRCA(dataset).chrm(chrm).crl(2,:))]);
            %         else
            %             title(['Chr: ', num2str(chrm)...
            %                '; score: ' num2str(dataset_BRCA(dataset).chrm(chrm).scr)]);
            %         end
            %'; ' num2str(dataset_BRCA(dataset).chrm(chrm).scores)...
            %, '; mean START error %chr: ' num2str(d_start) '; mean END error %chr: ' num2str(d_end)])
            
            if chrm==22
                xlabel(['Dataset: ' num2str(dataset) '; purity: ' num2str(dataset_BRCA(dataset).purity(5)) ])
            end
        end
        
        
        %     sigDNA05=dataset_BRCA(dataset).all_chrm_mat(:,6);
        %     sigRNA05=dataset_BRCA(dataset).all_chrm_mat(:,7);
        %     FDR=0.05/size(dataset_BRCA(dataset).all_chrm_mat,1);
        CNAs=dataset_BRCA(dataset).all_chrm_mat_Ttr(:,18);
        idx_s=(CNAs<0.2);
        idx_s_high=(CNAs>0);
        
        subtightplot(8,3,[21 24],[0.02 0.01])
        plot(dataset_BRCA(dataset).all_chrm_mat_Ttr(:,18),dataset_BRCA(dataset).all_chrm_mat_Ttr(:,7),'k.')
        hold on
        title_corr_str=[];
        title_corr_str1=[];
        title_corr_str2=[];
        
        if dataset_BRCA(dataset).fdr_corr(1)==1
            plot(dataset_BRCA(dataset).all_chrm_mat_Ttr(idx_s,18),dataset_BRCA(dataset).all_chrm_mat_Ttr(idx_s,7),'bo')
            title_corr_str1=['corr blue: ' num2str(dataset_BRCA(dataset).r2_corr(1)) '; ' num2str(dataset_BRCA(dataset).fdr_corr(1)) ';'];
        end
        if dataset_BRCA(dataset).fdr_corr(2)==1
            plot(dataset_BRCA(dataset).all_chrm_mat_Ttr(idx_s_high,18),dataset_BRCA(dataset).all_chrm_mat_Ttr(idx_s_high,7),'r+')
            title_corr_str2=['corr red: ' num2str(dataset_BRCA(dataset).r2_corr(2)) '; ' num2str(dataset_BRCA(dataset).fdr_corr(2)) ';'];
        end
        
        title_corr_str=[title_corr_str1 ' ' title_corr_str2];
        
        axis([-1 1 0.49 1.01])
        set(gca,'YAxisLocation','right')
        grid on
        hold off
        title(title_corr_str)
        pause,
        
    end
end
%% Tex
for dataset=1%:72
    tic,
    %convert_wval_to_circos_data(dataset_BRCA(dataset).all_chrm_mat,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_Ttr_1.txt')
    mat_data=participant_BRCA(dataset).all_vprs_mat_Tex;
    convert_vprs_to_circos_line(mat_data,1,0.365,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTex_b.txt')
    convert_vprs_to_circos_line(mat_data,2,0.535,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTtr_b.txt')
    convert_vprs_to_circos_line(mat_data,3,0.705,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPex_b.txt')
    convert_vprs_to_circos_line(mat_data,4,0.875,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPtr_b.txt')
    convert_vprs_to_circos_line(mat_data,1,0.365,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTex_t.txt')
    convert_vprs_to_circos_line(mat_data,2,0.535,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTtr_t.txt')
    convert_vprs_to_circos_line(mat_data,3,0.705,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPex_t.txt')
    convert_vprs_to_circos_line(mat_data,4,0.875,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPtr_t.txt')
    
    dmCNA=dataset_BRCA(dataset).CNA(:,[1:3 5]);
    convert_vprs_to_circos_data(dmCNA,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_CNA.txt')
    
    dmVAF=all_data(dataset).data(:,[1 2 3]);
    convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_NTex.txt')
    
    dmVAF=all_data(dataset).data(:,[1 2 4]);
    convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_NTtr.txt')
    
    dmVAF=all_data(dataset).data(:,[1 2 5]);
    convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_TPex.txt')
    
    dmVAF=all_data(dataset).data(:,[1 2 6]);
    convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_TPtr.txt')
    
    cd ~
    cd ../..
    cd /usr/local/Cellar/circos/0.69-6/mycircos
    system(['../bin/circos -conf etc/my_circos_all.conf -outputfile circos' num2str(dataset) 'Tex.png']);
    cd ~
    cd Dropbox/new_paper_Anelia/
    toc,
    disp(dataset)
end
% % Ttr
% for dataset=1%:72
%     tic,
%     %convert_wval_to_circos_data(dataset_BRCA(dataset).all_chrm_mat,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_Ttr_1.txt')
%     mat_data=dataset_BRCA(dataset).all_chrm_mat_Ttr;
%     convert_wval_to_circos_line(mat_data,1,0.365,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTex_b.txt')
%     convert_wval_to_circos_line(mat_data,2,0.535,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTtr_b.txt')
%     convert_wval_to_circos_line(mat_data,3,0.705,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPex_b.txt')
%     convert_wval_to_circos_line(mat_data,4,0.875,0.075,0,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPtr_b.txt')
%     convert_wval_to_circos_line(mat_data,1,0.365,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTex_t.txt')
%     convert_wval_to_circos_line(mat_data,2,0.535,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_NTtr_t.txt')
%     convert_wval_to_circos_line(mat_data,3,0.705,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPex_t.txt')
%     convert_wval_to_circos_line(mat_data,4,0.875,0.075,1,'/usr/local/Cellar/circos/0.69-6/mycircos/data/modes_TPtr_t.txt')
%     
%     dmCNA=dataset_BRCA(dataset).CNA(:,[1:3 5]);
%     convert_wval_to_circos_data(dmCNA,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_CNA.txt')
%     
%     dmVAF=all_data(dataset).data(:,[1 2 3]);
%     convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_NTex.txt')
%     
%     dmVAF=all_data(dataset).data(:,[1 2 4]);
%     convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_NTtr.txt')
%     
%     dmVAF=all_data(dataset).data(:,[1 2 5]);
%     convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_TPex.txt')
%     
%     dmVAF=all_data(dataset).data(:,[1 2 6]);
%     convert_VAF_to_circos_data(dmVAF,'/usr/local/Cellar/circos/0.69-6/mycircos/data/SAMPLE_TPtr.txt')
%     
%     cd ~
%     cd ../..
%     cd /usr/local/Cellar/circos/0.69-6/mycircos
%     system(['../bin/circos -conf etc/my_circos_all.conf -outputfile circos' num2str(dataset) 'Ttr.png']);
%     cd ~
%     cd Dropbox/new_paper_Anelia/
%     toc,
%     disp(dataset)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
cd ~
cd ../..
cd /usr/local/Cellar/circos/0.69-6/mycircos
system(['../bin/circos -conf etc/my_circos_all.conf -outputfile circos' num2str(dataset) '.png'])
cd ~
cd Documents/new_paper_Anelia/clean
%%
