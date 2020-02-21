clear all
load all_data_published.mat
load participant_BRCA
load all_CNA_data

%%
% for i=1:72
%     sampleID=participant_BRCA(i).sampleID;
%     fun = @(s)any(strncmpi(s,sampleID,12));
%     CNA_idx=cellfun(fun,all_CNA_data(1).id);
%     participant_BRCA(i).CNA=all_CNA_data(CNA_idx).data;
%     participant_BRCA(i).all_modes_mat_Tex=all_modes_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Tex');
%     participant_BRCA(i).all_modes_mat_Ttr=all_modes_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Ttr');
% end
for i=1:72
    sampleID=participant_BRCA(i).sampleID;
    fun = @(s)any(strncmpi(s,sampleID,12));
    CNA_idx=cellfun(fun,all_CNA_data(1).id);
    participant_BRCA(i).CNA=all_CNA_data(CNA_idx).data;
    participant_BRCA(i).all_modes_mat_Tex=all_vpr_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Tex');
    participant_BRCA(i).all_modes_mat_Ttr=all_vpr_to_matrix(participant_BRCA(i).chrm,participant_BRCA(i).CNA,'Ttr');
end

out_Tex=[];
for dataset=1:72
    out=participant_BRCA(dataset).all_modes_mat_Tex;
    out(:,21:25)=repmat(participant_BRCA(dataset).purity,size(out,1),1);
    out(:,26)=dataset;
    out_Tex=[out_Tex; out];
end
out_Tex=out_Tex(:,[26 1:25]);
out_Tex=sortrows(out_Tex,[1 2 3]);

out_Ttr=[];
for dataset=1:72
    out=participant_BRCA(dataset).all_modes_mat_Ttr;
    out(:,21:25)=repmat(participant_BRCA(dataset).purity,size(out,1),1);
    out(:,26)=dataset;
    out_Ttr=[out_Ttr; out];
end
out_Ttr=out_Ttr(:,[26 1:25]);
out_Ttr=sortrows(out_Ttr,[1 2 3]);
%%
clc
outD=out_Tex;
outR=out_Ttr;
size(outD,1)
% idx_nonD=outD(:,12)>1e-5 & outD(:,8)>=0.58 |...
%     outD(:,11)>1e-5 & outD(:,7)>=0.58;% |...
% %((outD(:,11)>1e-5)~=(outD(:,12)>1e-5)) & outD(:,18)>1e-5;
% outD(idx_nonD,:)=[];
% size(outD,1)

size(outR,1)
% idx_nonR=outR(:,12)>1e-5 & outR(:,8)>=0.58 |...
%     outR(:,11)>1e-5 & outR(:,7)>=0.58;% |...
% %((outR(:,11)>1e-5)~=(outR(:,12)>1e-5)) & outR(:,18)>1e-5;
% outR(idx_nonR,:)=[];
% size(outR,1)
%%
outD=out_Tex;
outR=out_Ttr;
k=1;
sDR=0;
sD=0;
sR=0;
sDj=0;
sRj=0;
mDR=0;
total_jmpD=0;
total_jmpR=0;
total_winD=0;
total_winR=0;
bts_prct_td=NaN(7,100000);
parfor bts=1:100000
    tic,
    k=1;
    td=[];
    for dataset=1:72
        %disp(dataset),
        all_data_pos=all_data(dataset).data(:,1:2);
        for chrm=1:22
            data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
            nb_pts=numel(data_pos);
            
            %disp(chrm)
            idxD=outD(:,1)==dataset & outD(:,2)==chrm;
            idxR=outR(:,1)==dataset & outR(:,2)==chrm;
            
            jmpD=outD(idxD,21);
            jmpR=outR(idxR,21);
            
            total_winD=total_winD+(numel(jmpD));
            
            total_winR=total_winR+(numel(jmpR));
            
            total_jmpD=total_jmpD+(numel(jmpD)-1);
            
            total_jmpR=total_jmpR+(numel(jmpR)-1);
            
            if numel(jmpD)>1 && numel(jmpR)>1
                jmpD(2:end)=randi([2 nb_pts],numel(jmpD)-1,1);
                jmpR(2:end)=randi([2 nb_pts],numel(jmpR)-1,1);
                
                mDR=mDR+1;
                if numel(jmpD)>numel(jmpR)
                    for i=1:numel(jmpR)
                        td(k)=min(abs(jmpD-jmpR(i)))/nb_pts;
                        k=k+1;
                    end
                elseif numel(jmpD)<numel(jmpR)
                    for i=1:numel(jmpD)
                        td(k)=min(abs(jmpR-jmpD(i)))/nb_pts;
                        k=k+1;
                    end
                else
                    for i=1:numel(jmpD)
                        td(k)=min(abs(jmpR-jmpD(i)))/nb_pts;
                        k=k+1;
                    end
                end
            elseif numel(jmpD)==1 && numel(jmpR)>1
                %disp('single in D; nb chg pts in R:')
                sD=sD+1;
                sDj=sDj+(numel(jmpR)-1);
            elseif numel(jmpD)>1 && numel(jmpR)==1
                %disp('single in R; nb chg pts in D:')
                sR=sR+1;
                sRj=sRj+(numel(jmpD)-1);
            elseif numel(jmpD)==1 && numel(jmpR)==1
                %disp('single in both')
                sDR=sDR+1;
            end
        end
        %pause,
    end
    toc,
    bts_prct_td(:,bts)=prctile(td,[2.5 25 50 75 80 90 97.5]);
end

%%
sDR=0;
sD=0;
sR=0;
sDj=0;
sRj=0;
mDR=0;
total_jmpD=0;
total_jmpR=0;
total_winD=0;
total_winR=0;
tic,
k=1;
td=[];
for dataset=1:72
    %disp(dataset),
    all_data_pos=all_data(dataset).data(:,1:2);
    for chrm=1:22
        data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
        nb_pts=numel(data_pos);
        
        %disp(chrm)
        idxD=outD(:,1)==dataset & outD(:,2)==chrm;
        idxR=outR(:,1)==dataset & outR(:,2)==chrm;
        
        jmpD=outD(idxD,21);
        jmpR=outR(idxR,21);
        
        total_winD=total_winD+(numel(jmpD));
        
        total_winR=total_winR+(numel(jmpR));
        
        total_jmpD=total_jmpD+(numel(jmpD)-1);
        
        total_jmpR=total_jmpR+(numel(jmpR)-1);
        
        if numel(jmpD)>1 && numel(jmpR)>1
            %jmpD(2:end)=randi([2 nb_pts],numel(jmpD)-1,1);
            %jmpR(2:end)=randi([2 nb_pts],numel(jmpR)-1,1);
            
            mDR=mDR+1;
            if numel(jmpD)>numel(jmpR)
                for i=1:numel(jmpR)
                    td(k)=min(abs(jmpD-jmpR(i)))/nb_pts;
                    k=k+1;
                end
            elseif numel(jmpD)<numel(jmpR)
                for i=1:numel(jmpD)
                    td(k)=min(abs(jmpR-jmpD(i)))/nb_pts;
                    k=k+1;
                end
            else
                for i=1:numel(jmpD)
                    td(k)=min(abs(jmpR-jmpD(i)))/nb_pts;
                    k=k+1;
                end
            end
        elseif numel(jmpD)==1 && numel(jmpR)>1
            %disp('single in D; nb chg pts in R:')
            sD=sD+1;
            sDj=sDj+(numel(jmpR)-1);
        elseif numel(jmpD)>1 && numel(jmpR)==1
            %disp('single in R; nb chg pts in D:')
            sR=sR+1;
            sRj=sRj+(numel(jmpD)-1);
        elseif numel(jmpD)==1 && numel(jmpR)==1
            %disp('single in both')
            sDR=sDR+1;
        end
    end
    %pause,
end
clc
total_jmpD
total_jmpR
mDR+sD+sR+sDR
mDR,
sD,
sR,
sDR,
mDR/1584
sD/1584
sR/1584
sDR/1584
mDR/1584+sD/1584+sR/1584
org_prct_td=prctile(td,[2.5 25 50 75 80 90 97.5])
median(bts_prct_td,2)'
sum(bts_prct_td'>org_prct_td)
%%
l2_Tex3Tex4=[];
l2_Ttr3Ttr4=[];
dl2_Tex3Tex4=[];
dl2_Ttr3Ttr4=[];

k=1;
k1=1;
for dataset=1:72%14%
    data=participant_BRCA(dataset);
    all_data_pos=all_data(dataset).data(:,1:2);
    for chrm=1:22%1%
        data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
        
        
        modes_based_on_l3=data.all_modes_mat_Tex;
        modes_based_on_l3(modes_based_on_l3(:,10)>1e-5 & modes_based_on_l3(:,6)>=0.58,6)=NaN;
        modes_based_on_l3(modes_based_on_l3(:,11)>1e-5 & modes_based_on_l3(:,7)>=0.58,7)=NaN;
        
        
        modes_based_on_l4=data.all_modes_mat_Ttr;
        modes_based_on_l4(modes_based_on_l4(:,10)>1e-5 & modes_based_on_l4(:,6)>=0.58,6)=NaN;
        modes_based_on_l4(modes_based_on_l4(:,11)>1e-5 & modes_based_on_l4(:,7)>=0.58,7)=NaN;
        
        modes_based_on_l3=modes_based_on_l3(modes_based_on_l3(:,1)==chrm,:);
        modes_based_on_l4=modes_based_on_l4(modes_based_on_l4(:,1)==chrm,:);
        
        if size(modes_based_on_l3,1)>1 || size(modes_based_on_l4,1)>1
            [modes_pos_bo_l3,ims_bo_l3]=unique([modes_based_on_l3(:,2); modes_based_on_l3(:,3)+eps]);
            all_mds_Ttr_bo_l3=[modes_based_on_l3(:,7); modes_based_on_l3(:,7)];
            all_mds_Ttr_bo_l3=all_mds_Ttr_bo_l3(ims_bo_l3);
            
            all_mds_Tex_bo_l3=[modes_based_on_l3(:,6); modes_based_on_l3(:,6)];
            all_mds_Tex_bo_l3=all_mds_Tex_bo_l3(ims_bo_l3);
            
            [modes_pos_bo_l4,ims_bo_l4]=unique([modes_based_on_l4(:,2); modes_based_on_l4(:,3)+eps]);
            all_mds_Ttr_bo_l4=[modes_based_on_l4(:,7); modes_based_on_l4(:,7)];
            all_mds_Ttr_bo_l4=all_mds_Ttr_bo_l4(ims_bo_l4);
            
            all_mds_Tex_bo_l4=[modes_based_on_l4(:,6); modes_based_on_l4(:,6)];
            all_mds_Tex_bo_l4=all_mds_Tex_bo_l4(ims_bo_l4);
            
            interp_mds_Tex_bo_l3=interp1(modes_pos_bo_l3,all_mds_Tex_bo_l3,data_pos,'nearest','extrap');
            interp_mds_Ttr_bo_l3=interp1(modes_pos_bo_l3,all_mds_Ttr_bo_l3,data_pos,'nearest','extrap');
            interp_mds_Tex_bo_l4=interp1(modes_pos_bo_l4,all_mds_Tex_bo_l4,data_pos,'nearest','extrap');
            interp_mds_Ttr_bo_l4=interp1(modes_pos_bo_l4,all_mds_Ttr_bo_l4,data_pos,'nearest','extrap');
            
            
            %         l2_Tex3Tex4(k)=nansum(abs(interp_mds_Tex_bo_l3-interp_mds_Tex_bo_l4))/numel(interp_mds_Tex_bo_l3);
            %         l2_Ttr3Ttr4(k)=nansum(abs(interp_mds_Ttr_bo_l3-interp_mds_Ttr_bo_l4))/numel(interp_mds_Ttr_bo_l4);
            %idx_nn=~isnan(interp_mds_Tex_bo_l3) & ~isnan(interp_mds_Tex_bo_l4);
            %l2_Tex3Tex4(k)=sqrt(trapz(data_pos(idx_nn),abs(interp_mds_Tex_bo_l3(idx_nn)-interp_mds_Tex_bo_l4(idx_nn)).^2))/(data_pos(end)-data_pos(1));%
            l2_Tex3Tex4(k)=nansum(abs(interp_mds_Tex_bo_l3-interp_mds_Tex_bo_l4))/numel(interp_mds_Tex_bo_l3);
            dl2_Tex3Tex4(k)=numel(interp_mds_Tex_bo_l3);

            %idx_nn=~isnan(interp_mds_Ttr_bo_l3) & ~isnan(interp_mds_Ttr_bo_l4);
            %l2_Ttr3Ttr4(k)=sqrt(trapz(data_pos(idx_nn),(interp_mds_Ttr_bo_l3(idx_nn)-interp_mds_Ttr_bo_l4(idx_nn)).^2))/(data_pos(end)-data_pos(1));%
            l2_Ttr3Ttr4(k)=nansum(abs(interp_mds_Ttr_bo_l3-interp_mds_Ttr_bo_l4))/numel(interp_mds_Ttr_bo_l4);
            dl2_Ttr3Ttr4(k)=numel(interp_mds_Ttr_bo_l4);

            
            if 0
                colr=lines(7);
                gap=[0.1 0.05];
                marg_h=[0.05 0.05];
                marg_w=[0.05 0.025];
                
                subtightplot(6,1,[1 2],gap,marg_h,marg_w)
                cla %dataset14 % hrm1
                plot(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,1),abs(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,5)-0.5)+0.5,...
                    '.','markersize',3,'color',colr(3,:))
                hold on,
                plot(data_pos,interp_mds_Ttr_bo_l3,'o-','markersize',3,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
                plot(data_pos,interp_mds_Ttr_bo_l4,'x-','markersize',3,'color',colr(3,:))
                set(gca,'xTickLabel','','yTick',[0.5 0.75 1])
                ylim([0.5 1])
                legend({'VAF_{TTR}','v_{PR,TTR} in windows based on VAF_{TEX}','v_{PR,TTR} in windows based on VAF_{TTR}'})
                
                subtightplot(6,1,3,gap,marg_h,marg_w)
                cla
                [~,ia] = unique(data_pos);
                bar(data_pos(ia),abs(interp_mds_Ttr_bo_l3(ia)-interp_mds_Ttr_bo_l4(ia)),100000,'k')
                ylim([0 0.3])
                set(gca,'xTickLabel','','yTick',[0 0.15 0.3])
                title('Absolute difference between v_{PR,TTR}''s in windows based on VAF_{TEX} and VAF_{TTR} signals.')
                
                
                subtightplot(6,1,[4 5],gap,marg_h,marg_w)
                cla
                plot(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,1),abs(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,4)-0.5)+0.5,...
                    '.','markersize',3,'color',colr(2,:),'MarkerFaceColor',colr(2,:))
                hold on,
                plot(data_pos,interp_mds_Tex_bo_l3,'o-','markersize',3,'color',colr(2,:),'MarkerFaceColor',colr(2,:))
                plot(data_pos,interp_mds_Tex_bo_l4,'x-','markersize',3,'color',[0.5 0.5 0.5])
                ylim([0.5 1])
                set(gca,'xTickLabel','','yTick',[0.5 0.75 1]) 
                legend({'VAF_{TEX}','v_{PR,TEX} in windows based on VAF_{TEX}','v_{PR,TEX} in windows based on VAF_{TTR}'})
                
                
                subtightplot(6,1,6,gap,marg_h,marg_w)
                cla
                bar(data_pos(ia),abs(interp_mds_Tex_bo_l3(ia)-interp_mds_Tex_bo_l4(ia)),100000,'k')
                ylim([0 0.3])
                set(gca,'xTickLabel','','yTick',[0 0.15 0.3])
                title('Absolute difference between v_{PR,TEX}''s in windows based on VAF_{TEX} and VAF_{TTR} signals.')
                xlabel('BP along a chromosome')
                hold off
                %title(num2str([modes_l3(:,6); varTexl3(k); c_Tex3Tex4(k)]))
                pause,
                
            end
            
            k=k+1,
        elseif size(modes_based_on_l3,1)==1 && size(modes_based_on_l4,1)==1
            k1=k1+1,
        end
    end
end
%%
btsl2_Tex3Tex4=NaN(1584,1000);
btsl2_Ttr3Ttr4=NaN(1584,1000);

for bts=1:1000
    bts,
    tic,
    sl2_Tex3Tex4=NaN(72,22);
    sl2_Ttr3Ttr4=NaN(72,22);
    
    for dataset=1:72
        data=participant_BRCA(dataset);
        all_data_pos=all_data(dataset).data(:,1:2);
        for chrm=1:22
            data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
            
            modes_based_on_l3=data.all_modes_mat_Tex;
            modes_based_on_l4=data.all_modes_mat_Ttr;
            
            modes_based_on_l3(:,[6 7 10 11])=modes_based_on_l3(randperm(size(modes_based_on_l3,1)),[6 7 10 11]);
            modes_based_on_l3(modes_based_on_l3(:,10)>1e-5 & modes_based_on_l3(:,6)>=0.58,6)=NaN;
            modes_based_on_l3(modes_based_on_l3(:,11)>1e-5 & modes_based_on_l3(:,7)>=0.58,7)=NaN;
            
            modes_based_on_l4(:,[6 7 10 11])=modes_based_on_l4(randperm(size(modes_based_on_l4,1)),[6 7 10 11]);
            modes_based_on_l4(modes_based_on_l4(:,10)>1e-5 & modes_based_on_l4(:,6)>=0.58,6)=NaN;
            modes_based_on_l4(modes_based_on_l4(:,11)>1e-5 & modes_based_on_l4(:,7)>=0.58,7)=NaN;
            
            modes_based_on_l3=modes_based_on_l3(modes_based_on_l3(:,1)==chrm,:);
            modes_based_on_l4=modes_based_on_l4(modes_based_on_l4(:,1)==chrm,:);
            
            if size(modes_based_on_l3,1)>1 || size(modes_based_on_l4,1)>1
                
                [modes_pos_bo_l3,ims_bo_l3]=unique([modes_based_on_l3(:,2); modes_based_on_l3(:,3)+eps]);
                all_mds_Ttr_bo_l3=[modes_based_on_l3(:,7); modes_based_on_l3(:,7)];
                all_mds_Ttr_bo_l3=all_mds_Ttr_bo_l3(ims_bo_l3);
                
                all_mds_Tex_bo_l3=[modes_based_on_l3(:,6); modes_based_on_l3(:,6)];
                all_mds_Tex_bo_l3=all_mds_Tex_bo_l3(ims_bo_l3);
                
                [modes_pos_bo_l4,ims_bo_l4]=unique([modes_based_on_l4(:,2); modes_based_on_l4(:,3)+eps]);
                all_mds_Ttr_bo_l4=[modes_based_on_l4(:,7); modes_based_on_l4(:,7)];
                all_mds_Ttr_bo_l4=all_mds_Ttr_bo_l4(ims_bo_l4);
                
                all_mds_Tex_bo_l4=[modes_based_on_l4(:,6); modes_based_on_l4(:,6)];
                all_mds_Tex_bo_l4=all_mds_Tex_bo_l4(ims_bo_l4);
                
                interp_mds_Tex_bo_l3=interp1(modes_pos_bo_l3,all_mds_Tex_bo_l3,data_pos,'nearest','extrap');
                interp_mds_Ttr_bo_l3=interp1(modes_pos_bo_l3,all_mds_Ttr_bo_l3,data_pos,'nearest','extrap');
                interp_mds_Tex_bo_l4=interp1(modes_pos_bo_l4,all_mds_Tex_bo_l4,data_pos,'nearest','extrap');
                interp_mds_Ttr_bo_l4=interp1(modes_pos_bo_l4,all_mds_Ttr_bo_l4,data_pos,'nearest','extrap');
                
                sl2_Tex3Tex4(dataset,chrm)=nansum(abs(interp_mds_Tex_bo_l3-interp_mds_Tex_bo_l4))/numel(interp_mds_Tex_bo_l3);
                sl2_Ttr3Ttr4(dataset,chrm)=nansum(abs(interp_mds_Ttr_bo_l3-interp_mds_Ttr_bo_l4))/numel(interp_mds_Ttr_bo_l4);
            end
        end
    end
    
    toc
    btsl2_Tex3Tex4(:,bts)=sl2_Tex3Tex4(:);
    btsl2_Ttr3Ttr4(:,bts)=sl2_Ttr3Ttr4(:);
end
%%
clc
disp('Tex3 vs Tex4')
nansum(l2_Tex3Tex4<0.03)/875%/1584
nansum(l2_Tex3Tex4==0)/875%/1584
max(sum(btsl2_Tex3Tex4<0.03)/875)%/1584)
max(sum(btsl2_Tex3Tex4==0)/875)%/1584)
prctile(l2_Tex3Tex4,[5,25,50,75,97.5])
sum(prctile(btsl2_Tex3Tex4,[5,25,50,75,97.5])'>prctile(l2_Tex3Tex4,[5,25,50,75,97.5]))
%
disp('Ttr3 vs Ttr4')
sum(l2_Ttr3Ttr4<0.03)/875%/1584
sum(l2_Ttr3Ttr4==0)/875%/1584
sum(l2_Ttr3Ttr4==0 & l2_Tex3Tex4==0)/875%/1584
max(sum(btsl2_Ttr3Ttr4<0.03)/875)%/1584)
max(sum(btsl2_Ttr3Ttr4==0)/875)%/1584)
max(sum(btsl2_Ttr3Ttr4==0 & btsl2_Tex3Tex4==0)/875)%/1584)
prctile(l2_Ttr3Ttr4,[5,25,50,75,97.5])
sum(prctile(btsl2_Ttr3Ttr4,[5,25,50,75,97.5])'>prctile(l2_Ttr3Ttr4,[5,25,50,75,97.5]))
%%
disp('Tex3 vs Ttr3')
sum(l2_Tex3Ttr3<0.03)/1584
sum(l2_Tex3Ttr3==0)/1584

disp('Tex4 vs Ttr4')
sum(l2_Tex4Ttr4<0.03)/1584
sum(l2_Tex4Ttr4==0)/1584
%%
clc
disp('Tex3 vs CNA')
s1=sum(nTexl3==1 & nCNA==1)
s2=sum(nTexl3==1 & nCNA>1 & ivarCNA<0.001)
s3=sum(nTexl3>1 & nCNA>1 & bonf_holm(p_TexCNAl3,0.05)<0.05 & abs(c_TexCNAl3)>0.65)
(s1+s2+s3)/1584
%sum(nTexl3>1 & nCNA==1)

disp('Ttr3 vs CNA')
s1=sum(nTtrl3==1 & nCNA==1)
s2=sum(nTtrl3==1 & nCNA>1 & ivarCNA<0.001)
s3=sum(nTtrl3>1 & nCNA>1 & bonf_holm(p_TtrCNAl3,0.05)<0.05 & abs(c_TtrCNAl3)>0.65)
(s1+s2+s3)/1584
%sum(nTtrl3>1 & nCNA==1)

disp('Tex4 vs CNV')
s1=sum(nTexl4==1 & nCNA==1)
s2=sum(nTexl4==1 & nCNA>1 & ivarCNA<0.001)
s3=sum(nTexl4>1 & nCNA>1 & bonf_holm(p_TexCNAl4,0.05)<0.05 & abs(c_TexCNAl4)>0.65)
(s1+s2+s3)/1584
%sum(nTexl4>1 & nCNA==1)

disp('Ttr4 vs CNV')
s1=sum(nTtrl4==1 & nCNA==1)
s2=sum(nTtrl4==1 & nCNA>1 & ivarCNA<0.001)
s3=sum(nTtrl4>1 & nCNA>1 & bonf_holm(p_TtrCNAl4,0.05)<0.05 & abs(c_TtrCNAl4)>0.65)
(s1+s2+s3)/1584
%sum(nTtrl4>1 & nCNA==1)
%%
subplot(2,2,1)
histogram(c_TexCNAl3(nTexl3>1 & nCNA>1 & p_TexCNAl3<(0.05/k)),-1:0.01:1)
subplot(2,2,2)
histogram(c_TexCNAl4(nTexl4>1 & nCNA>1 & p_TexCNAl4<(0.05/k)),-1:0.01:1)
subplot(2,2,3)
histogram(c_TtrCNAl3(nTtrl3>1 & nCNA>1 & p_TtrCNAl3<(0.05/k)),-1:0.01:1)
subplot(2,2,4)
histogram(c_TtrCNAl4(nTtrl4>1 & nCNA>1 & p_TtrCNAl4<(0.05/k)),-1:0.01:1)








