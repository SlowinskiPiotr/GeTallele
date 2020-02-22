out_Tex=[];
for dataset=1:72
    out=participant_BRCA(dataset).all_vprs_mat_Tex;
    out(:,21:25)=repmat(participant_BRCA(dataset).purity,size(out,1),1);
    out(:,26)=dataset;
    out_Tex=[out_Tex; out];
end
out_Tex=out_Tex(:,[26 1:25]);
out_Tex=sortrows(out_Tex,[1 2 3]);

out_Ttr=[];
for dataset=1:72
    out=participant_BRCA(dataset).all_vprs_mat_Ttr;
    out(:,21:25)=repmat(participant_BRCA(dataset).purity,size(out,1),1);
    out(:,26)=dataset;
    out_Ttr=[out_Ttr; out];
end
out_Ttr=out_Ttr(:,[26 1:25]);
out_Ttr=sortrows(out_Ttr,[1 2 3]);

clc
outTex=out_Tex;
outTtr=out_Ttr;
disp('windows based on Tex')
size(outTex,1)
disp('windows based on Ttr')
size(outTtr,1)

% chromosome-wise alignment of the start and end points of the windows
single_segment_Tex_Ttr=0;
single_segment_Tex=0;
single_segment_Ttr=0;
multiple_segments_Tex_Ttr=0;
total_segments_Tex=0;
total_segments_Ttr=0;

k=1;
difference_start_end=[];

for dataset=1:72
    all_data_pos=all_data(dataset).data(:,1:2);
    for chrm=1:22
        data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
        nb_pts=numel(data_pos);
        
        idxTex=outTex(:,1)==dataset & outTex(:,2)==chrm;
        idxTtr=outTtr(:,1)==dataset & outTtr(:,2)==chrm;
        
        chgpts_Tex=outTex(idxTex,21);
        chgpts_Ttr=outTtr(idxTtr,21);
        
        total_segments_Tex=total_segments_Tex+(numel(chgpts_Tex));
        
        total_segments_Ttr=total_segments_Ttr+(numel(chgpts_Ttr));
        
        if numel(chgpts_Tex)>1 && numel(chgpts_Ttr)>1
            
            multiple_segments_Tex_Ttr=multiple_segments_Tex_Ttr+1;
            if numel(chgpts_Tex)>numel(chgpts_Ttr)
                for i=1:numel(chgpts_Ttr)
                    difference_start_end(k)=min(abs(chgpts_Tex-chgpts_Ttr(i)))/nb_pts;
                    k=k+1;
                end
            elseif numel(chgpts_Tex)<numel(chgpts_Ttr)
                for i=1:numel(chgpts_Tex)
                    difference_start_end(k)=min(abs(chgpts_Ttr-chgpts_Tex(i)))/nb_pts;
                    k=k+1;
                end
            else
                for i=1:numel(chgpts_Tex)
                    difference_start_end(k)=min(abs(chgpts_Ttr-chgpts_Tex(i)))/nb_pts;
                    k=k+1;
                end
            end
        elseif numel(chgpts_Tex)==1 && numel(chgpts_Ttr)>1
            single_segment_Tex=single_segment_Tex+1;
        elseif numel(chgpts_Tex)>1 && numel(chgpts_Ttr)==1
            single_segment_Ttr=single_segment_Ttr+1;
        elseif numel(chgpts_Tex)==1 && numel(chgpts_Ttr)==1
            single_segment_Tex_Ttr=single_segment_Tex_Ttr+1;
        end
    end
end
%% bootstrap
bts_prct_difference_start_end=NaN(7,100000);
parfor bts=1:100000
    bts,
    k=1;
    bts_difference_start_end=[];
    for dataset=1:72
        all_data_pos=all_data(dataset).data(:,1:2);
        for chrm=1:22
            data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
            nb_pts=numel(data_pos);
            
            idxTex=outTex(:,1)==dataset & outTex(:,2)==chrm;
            idxTtr=outTtr(:,1)==dataset & outTtr(:,2)==chrm;
            
            chgpts_Tex=outTex(idxTex,21);
            chgpts_Ttr=outTtr(idxTtr,21);
            
            
            if numel(chgpts_Tex)>1 && numel(chgpts_Ttr)>1
                chgpts_Tex(2:end)=randi([2 nb_pts],numel(chgpts_Tex)-1,1);
                chgpts_Ttr(2:end)=randi([2 nb_pts],numel(chgpts_Ttr)-1,1);                     
                if numel(chgpts_Tex)>numel(chgpts_Ttr)
                    for i=1:numel(chgpts_Ttr)
                        bts_difference_start_end(k)=min(abs(chgpts_Tex-chgpts_Ttr(i)))/nb_pts;
                        k=k+1;
                    end
                elseif numel(chgpts_Tex)<numel(chgpts_Ttr)
                    for i=1:numel(chgpts_Tex)
                        bts_difference_start_end(k)=min(abs(chgpts_Ttr-chgpts_Tex(i)))/nb_pts;
                        k=k+1;
                    end
                else
                    for i=1:numel(chgpts_Tex)
                        bts_difference_start_end(k)=min(abs(chgpts_Ttr-chgpts_Tex(i)))/nb_pts;
                        k=k+1;
                    end
                end
            end
        end
    end
    bts_prct_difference_start_end(:,bts)=prctile(bts_difference_start_end,[2.5 25 50 75 80 90 97.5]);
end
%%
clc
disp('total number of windows/ segments based on Tex')
total_segments_Tex
disp('total number of windows/ segments based on Ttr')
total_segments_Ttr
disp('% of chromosomes with single windows/ segment in Tex and Ttr')
single_segment_Tex_Ttr/(72*22),
disp('% of chromosomes with many windows/ segment in Tex and many in Ttr')
multiple_segments_Tex_Ttr/(72*22),
disp('% of chromosomes with single windows/ segment in Tex and many in Ttr')
single_segment_Tex/(72*22),
disp('% of chromosomes with single windows/ segment in Ttr and many in Tex')
single_segment_Ttr/(72*22),
disp('% of chromosomes with many windows/ segment in Tex or in Ttr')
multiple_segments_Tex_Ttr/(72*22)+single_segment_Tex/(72*22)+single_segment_Ttr/(72*22)
disp('50 75 90 percentiles of % difference in terms of the number of data points in the chromosome')
org_prct_difference_start_end=prctile(difference_start_end,[50 75 90])
disp('median 50 75 90 percentiles of shuffeld data')
median(bts_prct_difference_start_end([3 4 6],:),2)'
disp('# of shuffeld percentiles > than original percentiles')
sum(bts_prct_difference_start_end([3 4 6],:)'>org_prct_difference_start_end)
%% comparison of the vPR values in the 55% of chromosomes where at least one signal produced more than one window
MAE_Tex=[];
MAE_Ttr=[];

k=1;
k1=1;
for dataset=1:72
    data=participant_BRCA(dataset);
    all_data_pos=all_data(dataset).data(:,1:2);
    for chrm=1:22
        data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
        
        vprs_based_on_Tex=data.all_vprs_mat_Tex;
        vprs_based_on_Tex(vprs_based_on_Tex(:,10)>1e-5 & vprs_based_on_Tex(:,6)>=0.58,6)=NaN;
        vprs_based_on_Tex(vprs_based_on_Tex(:,11)>1e-5 & vprs_based_on_Tex(:,7)>=0.58,7)=NaN;
        
        vprs_based_on_Ttr=data.all_vprs_mat_Ttr;
        vprs_based_on_Ttr(vprs_based_on_Ttr(:,10)>1e-5 & vprs_based_on_Ttr(:,6)>=0.58,6)=NaN;
        vprs_based_on_Ttr(vprs_based_on_Ttr(:,11)>1e-5 & vprs_based_on_Ttr(:,7)>=0.58,7)=NaN;
        
        vprs_based_on_Tex=vprs_based_on_Tex(vprs_based_on_Tex(:,1)==chrm,:);
        vprs_based_on_Ttr=vprs_based_on_Ttr(vprs_based_on_Ttr(:,1)==chrm,:);
        
        if size(vprs_based_on_Tex,1)>1 || size(vprs_based_on_Ttr,1)>1
            [vprs_pos_based_on_Tex,ims_based_on_Tex]=unique([vprs_based_on_Tex(:,2); vprs_based_on_Tex(:,3)+eps]);
            all_vprs_Ttr_based_on_Tex=[vprs_based_on_Tex(:,7); vprs_based_on_Tex(:,7)];
            all_vprs_Ttr_based_on_Tex=all_vprs_Ttr_based_on_Tex(ims_based_on_Tex);
            
            all_vprs_Tex_based_on_Tex=[vprs_based_on_Tex(:,6); vprs_based_on_Tex(:,6)];
            all_vprs_Tex_based_on_Tex=all_vprs_Tex_based_on_Tex(ims_based_on_Tex);
            
            [vprs_pos_based_on_Ttr,ims_based_on_Ttr]=unique([vprs_based_on_Ttr(:,2); vprs_based_on_Ttr(:,3)+eps]);
            all_vprs_Ttr_based_on_Ttr=[vprs_based_on_Ttr(:,7); vprs_based_on_Ttr(:,7)];
            all_vprs_Ttr_based_on_Ttr=all_vprs_Ttr_based_on_Ttr(ims_based_on_Ttr);
            
            all_vprs_Tex_based_on_Ttr=[vprs_based_on_Ttr(:,6); vprs_based_on_Ttr(:,6)];
            all_vprs_Tex_based_on_Ttr=all_vprs_Tex_based_on_Ttr(ims_based_on_Ttr);
            
            interp_vprs_Tex_based_on_Tex=interp1(vprs_pos_based_on_Tex,all_vprs_Tex_based_on_Tex,data_pos,'nearest','extrap');
            interp_vprs_Ttr_based_on_Tex=interp1(vprs_pos_based_on_Tex,all_vprs_Ttr_based_on_Tex,data_pos,'nearest','extrap');
            interp_vprs_Tex_based_on_Ttr=interp1(vprs_pos_based_on_Ttr,all_vprs_Tex_based_on_Ttr,data_pos,'nearest','extrap');
            interp_vprs_Ttr_based_on_Ttr=interp1(vprs_pos_based_on_Ttr,all_vprs_Ttr_based_on_Ttr,data_pos,'nearest','extrap');
            
            
            MAE_Tex(k)=nansum(abs(interp_vprs_Tex_based_on_Tex-interp_vprs_Tex_based_on_Ttr))/numel(interp_vprs_Tex_based_on_Tex);

            MAE_Ttr(k)=nansum(abs(interp_vprs_Ttr_based_on_Tex-interp_vprs_Ttr_based_on_Ttr))/numel(interp_vprs_Ttr_based_on_Ttr);
            
            k=k+1;
        end
    end
end
%% bootstrap
bts_MAE_Tex=NaN(1584,1000);
bts_MAE_Ttr=NaN(1584,1000);

for bts=1:1000
    bts,
    tic,
    bMAE_Tex=NaN(72,22);
    bMAE_Ttr=NaN(72,22);
    
    for dataset=1:72
        data=participant_BRCA(dataset);
        all_data_pos=all_data(dataset).data(:,1:2);
        for chrm=1:22
            data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);
            
            vprs_based_on_Tex=data.all_vprs_mat_Tex;
            vprs_based_on_Ttr=data.all_vprs_mat_Ttr;
            
            vprs_based_on_Tex(:,[6 7 10 11])=vprs_based_on_Tex(randperm(size(vprs_based_on_Tex,1)),[6 7 10 11]);
            vprs_based_on_Tex(vprs_based_on_Tex(:,10)>1e-5 & vprs_based_on_Tex(:,6)>=0.58,6)=NaN;
            vprs_based_on_Tex(vprs_based_on_Tex(:,11)>1e-5 & vprs_based_on_Tex(:,7)>=0.58,7)=NaN;
            
            vprs_based_on_Ttr(:,[6 7 10 11])=vprs_based_on_Ttr(randperm(size(vprs_based_on_Ttr,1)),[6 7 10 11]);
            vprs_based_on_Ttr(vprs_based_on_Ttr(:,10)>1e-5 & vprs_based_on_Ttr(:,6)>=0.58,6)=NaN;
            vprs_based_on_Ttr(vprs_based_on_Ttr(:,11)>1e-5 & vprs_based_on_Ttr(:,7)>=0.58,7)=NaN;
            
            vprs_based_on_Tex=vprs_based_on_Tex(vprs_based_on_Tex(:,1)==chrm,:);
            vprs_based_on_Ttr=vprs_based_on_Ttr(vprs_based_on_Ttr(:,1)==chrm,:);
            
            if size(vprs_based_on_Tex,1)>1 || size(vprs_based_on_Ttr,1)>1
                
                [vprs_pos_based_on_Tex,ims_based_on_Tex]=unique([vprs_based_on_Tex(:,2); vprs_based_on_Tex(:,3)+eps]);
                all_vprs_Ttr_based_on_Tex=[vprs_based_on_Tex(:,7); vprs_based_on_Tex(:,7)];
                all_vprs_Ttr_based_on_Tex=all_vprs_Ttr_based_on_Tex(ims_based_on_Tex);
                
                all_vprs_Tex_based_on_Tex=[vprs_based_on_Tex(:,6); vprs_based_on_Tex(:,6)];
                all_vprs_Tex_based_on_Tex=all_vprs_Tex_based_on_Tex(ims_based_on_Tex);
                
                [vprs_pos_based_on_Ttr,ims_based_on_Ttr]=unique([vprs_based_on_Ttr(:,2); vprs_based_on_Ttr(:,3)+eps]);
                all_vprs_Ttr_based_on_Ttr=[vprs_based_on_Ttr(:,7); vprs_based_on_Ttr(:,7)];
                all_vprs_Ttr_based_on_Ttr=all_vprs_Ttr_based_on_Ttr(ims_based_on_Ttr);
                
                all_vprs_Tex_based_on_Ttr=[vprs_based_on_Ttr(:,6); vprs_based_on_Ttr(:,6)];
                all_vprs_Tex_based_on_Ttr=all_vprs_Tex_based_on_Ttr(ims_based_on_Ttr);
                
                interp_vprs_Tex_based_on_Tex=interp1(vprs_pos_based_on_Tex,all_vprs_Tex_based_on_Tex,data_pos,'nearest','extrap');
                interp_vprs_Ttr_based_on_Tex=interp1(vprs_pos_based_on_Tex,all_vprs_Ttr_based_on_Tex,data_pos,'nearest','extrap');
                interp_vprs_Tex_based_on_Ttr=interp1(vprs_pos_based_on_Ttr,all_vprs_Tex_based_on_Ttr,data_pos,'nearest','extrap');
                interp_vprs_Ttr_based_on_Ttr=interp1(vprs_pos_based_on_Ttr,all_vprs_Ttr_based_on_Ttr,data_pos,'nearest','extrap');
                
                bMAE_Tex(dataset,chrm)=nansum(abs(interp_vprs_Tex_based_on_Tex-interp_vprs_Tex_based_on_Ttr))/numel(interp_vprs_Tex_based_on_Tex);
                bMAE_Ttr(dataset,chrm)=nansum(abs(interp_vprs_Ttr_based_on_Tex-interp_vprs_Ttr_based_on_Ttr))/numel(interp_vprs_Ttr_based_on_Ttr);
            end
        end
    end
    
    toc
    bts_MAE_Tex(:,bts)=bMAE_Tex(:);
    bts_MAE_Ttr(:,bts)=bMAE_Ttr(:);
end
%%
clc
disp('MAE between vPR,Tex based on Tex vs vPR,Tex based on Ttr')
disp('perfect agreement, MAE=0')
nansum(MAE_Tex==0)/(multiple_segments_Tex_Ttr+single_segment_Tex+single_segment_Ttr)
disp('50 75 97.5 percentiles of MAE')
prctile(MAE_Tex,[50,75,97.5])
disp('# of MAE of shuffeld > than original MAE')
sum(prctile(bts_MAE_Tex,[50,75,97.5])'>prctile(MAE_Tex,[50,75,97.5]))

disp('MAE between vPR,Ttr based on Tex vs vPR,Ttr based on Ttr')
disp('perfect agreement, MAE=0')
sum(MAE_Ttr==0)/(multiple_segments_Tex_Ttr+single_segment_Tex+single_segment_Ttr)
disp('perfect agreement in Tex and Ttr')
sum(MAE_Ttr==0 & MAE_Tex==0)/(multiple_segments_Tex_Ttr+single_segment_Tex+single_segment_Ttr)
disp('50 75 97.5 percentiles of MAE')
prctile(MAE_Ttr,[50,75,97.5])
disp('# of MAE of shuffeld > than original MAE')
sum(prctile(bts_MAE_Ttr,[50,75,97.5])'>prctile(MAE_Ttr,[50,75,97.5]))