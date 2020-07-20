% compute correlations between vPR and CNA

for dataset=72:-1:1
    data=participant_BRCA(dataset).all_vprs_mat_Tex;
    data(data(:,10)>1e-5 & data(:,6)>=0.58,:)=[];
    data(data(:,11)>1e-5 & data(:,7)>=0.58,:)=[];
    
    CNAs=data(:,18);
    
    idx_s_min=find(CNAs>=-0.3 & CNAs<=0.3);
    idx_th=find(data(idx_s_min,6)==min(data(idx_s_min,6)));
    
    thr=median(data(idx_s_min(idx_th),18));
    thr_l=thr-0.05;
    thr_h=thr+0.05;
    
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

disp('vPR,Tex and CNA ampl')
sum(fdrph<0.05)
disp('vPR,Ttr and CNA ampl')
sum(fdrph7<0.05)
disp('vPR,Tex and CNA del')
sum(fdrpl<0.05)
disp('vPR,Ttr and CNA del')
sum(fdrpl7<0.05)