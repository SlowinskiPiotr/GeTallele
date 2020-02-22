% first run Correlation_between_vPR_and_CNA.m
% uses subtightplot function from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

figure(1)
colr=lines(3);
clf

gap=[0.1 0.05];
marg_h=[0.075 0.05];
marg_w=[0.075 0.05];

idx_dataset=[21 65 17 8];
thr=[];

for dataset=1:4
    
    subtightplot(2,2,dataset,gap,marg_h,marg_w)
    data=participant_BRCA(idx_dataset(dataset)).all_vprs_mat_Tex;
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
    
    plot(data(idx_s_low,18),data(idx_s_low,7),'o','color',colr(3,:),'markersize',3);
    plot(data(idx_s_high,18),data(idx_s_high,7),'o','color',colr(3,:),'markersize',3);
    plot(data(idx_s_mid,18),data(idx_s_mid,7),'ko','markersize',3);
    
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

        if sum(idx_s_mid)==sum(idx_s_high) || fdrph7(idx_dataset(dataset))>0.05 || plyh7(1)<0
            delete(h(2));
            nsig=nsig-1;
        end

        if sum(idx_s_mid)==sum(idx_s_low) || fdrpl(idx_dataset(dataset))>0.05 || plyl(1)>0
            delete(h(6));
            nsig=nsig-1;
        end
        if sum(idx_s_mid)==sum(idx_s_high) || fdrph(idx_dataset(dataset))>0.05 || plyh(1)<0
            delete(h(5));
            nsig=nsig-1;
        end
    else
        delete(h([2 3 5 6]))
    end
    
    ylim([0.5 1])
    set(gca,'YTick',0.5:0.1:1,'YTickLabel',[],'XTick',[-0.3 0 0.3],'XTickLabel',[])
    grid on

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

disp('Panel B')
dts=65;
disp('vPR,Tex and CNA ampl: r, pFDR')
[rl(dts), fdrpl(dts)]
disp('vPR,Ttr and CNA ampl: r, pFDR')
[rl7(dts) , fdrpl7(dts)]
disp('vPR,Tex and CNA del: r, pFDR')
[rh(dts), fdrph(dts)]
disp('vPR,Ttr and CNA del: r, pFDR')    
[rh7(dts), fdrph7(dts)]
disp('Panel C')
dts=17;
disp('vPR,Tex and CNA ampl: r, pFDR')
[rl(dts), fdrpl(dts)]
disp('vPR,Ttr and CNA ampl: r, pFDR')
[rl7(dts) , fdrpl7(dts)]
disp('vPR,Tex and CNA del: r, pFDR')
[rh(dts), fdrph(dts)]
disp('vPR,Ttr and CNA del: r, pFDR')    
[rh7(dts), fdrph7(dts)]
disp('Panel D')
dts=8;
disp('vPR,Tex and CNA ampl: r, pFDR')
[rl(dts), fdrpl(dts)]
disp('vPR,Ttr and CNA ampl: r, pFDR')
[rl7(dts) , fdrpl7(dts)]
disp('vPR,Tex and CNA del: r, pFDR')
[rh(dts), fdrph(dts)]
disp('vPR,Ttr and CNA del: r, pFDR')    
[rh7(dts), fdrph7(dts)]