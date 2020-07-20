% uses subtightplot function from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
colr=lines(3);
clf

gap=[0.02 0.02];
marg_h=[0.035 0.02];
marg_w=[0.05 0.01];

for dataset=1:72
    subtightplot(12,6,dataset,gap,marg_h,marg_w)
    data=participant_BRCA(dataset).all_vprs_mat_Tex;
    data(data(:,10)>1e-5 & data(:,6)>=0.58,:)=[];
    data(data(:,11)>1e-5 & data(:,7)>=0.58,:)=[];
    CNAs=data(:,18);
    
    idx_s_min=find(CNAs>=-0.3 & CNAs<=0.3);
    idx_th=find(data(idx_s_min,6)==min(data(idx_s_min,6)));
    
    thr=median(data(idx_s_min(idx_th),18));
    thr_l=thr-0.05;
    thr_h=thr+0.05;
    
    idx_s_low=(CNAs<thr_h);
    idx_s_high=(CNAs>thr_l);
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
    h(5)=plot(sort([data(idx_s_high,18); 3]),y,'-','color',colr(2,:));
    
    plyl=polyfit(data(idx_s_low,18),data(idx_s_low,6),1);
    y=polyval(plyl,sort([-1.5; data(idx_s_low,18)]));
    h(6)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-','color',colr(2,:));
    
    plyh7=polyfit(data(idx_s_high,18),data(idx_s_high,7),1);
    y=polyval(plyh7,sort([data(idx_s_high,18); 3]));
    h(2)=plot(sort([data(idx_s_high,18); 3]),y,'-','color',colr(3,:));
    
    plyl7=polyfit(data(idx_s_low,18),data(idx_s_low,7),1);
    y=polyval(plyl7,sort([-1.5; data(idx_s_low,18)]));
    h(3)=plot(sort([-1.5; data(idx_s_low,18)]),y,'-','color',colr(3,:));
    
    nsig=0;
    if max(data(:,18))>0.3 || min(data(:,18))<-0.3
        if numel(h)>3
            nsig=4;
            delete(h(4));
            if sum(idx_s_mid)==sum(idx_s_low) || fdrpl7(dataset)>0.05 || plyl7(1)>0
                delete(h(3));
                nsig=nsig-1;
            end
            if sum(idx_s_mid)==sum(idx_s_high) || fdrph7(dataset)>0.05 || plyh7(1)<0
                delete(h(2));
                nsig=nsig-1;
            end            
            if sum(idx_s_mid)==sum(idx_s_low) || fdrpl(dataset)>0.05 || plyl(1)>0
                delete(h(6));
                nsig=nsig-1;
            end
            if sum(idx_s_mid)==sum(idx_s_high) || fdrph(dataset)>0.05 || plyh(1)<0
                delete(h(5));
                nsig=nsig-1;d
            end
        else
            delete(h(1));
            if sum(idx_s_mid)==sum(idx_s_low) || fdrpl(dataset)>0.05 || plyl(1)>0
                delete(h(3));
            end
            if sum(idx_s_mid)==sum(idx_s_high) || fdrph(dataset)>0.05 || plyh(1)<0
                delete(h(2));
            end
        end
    else
        delete(h([2 3 5 6]))
    end
    
    ylim([0.5 1])
    xlim([-1.5 3])

    set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[],'XTick',0,'XTickLabel',[])
    grid on
    
    title([num2str(dataset) ': #sig ' num2str(nsig)])
    if any(1:6:67==dataset)
        set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[0.5 0.6 1],'XTick',0,'XTickLabel',[])
    end
    if dataset==67

        set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[0.5 0.6 1],'XTick',[-1.5 0 3],'XTickLabel',{'-1.5','0','3'})
    end
    if any(68:72==dataset)
        set(gca,'YTick',[0.5 0.6 1],'YTickLabel',[],'XTick',[-1.5 0 3],'XTickLabel',{'-1.5','0','3'})
    end
    if dataset==69
        xlabel('CNA')
    end
    if dataset==31
        ylabel('V_{pr}')
    end
    drawnow
    hold off
end