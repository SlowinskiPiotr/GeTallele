% for each dataset collect vprs from all the chromosomes
vprs=[];
for dataset=numel(participant_BRCA):-1:1
    vprs(dataset).vpr1=[];
    vprs(dataset).vpr2=[];
    vprs(dataset).vpr3=[];
    vprs(dataset).vpr4=[];
    
    vprs(dataset).p1_05=[];
    vprs(dataset).p2_05=[];
    vprs(dataset).p3_05=[];
    vprs(dataset).p4_05=[];
    
    vprs(dataset).pNexTex=[];
    vprs(dataset).pNtrTtr=[];

    vprs(dataset).wl=[];

    for chr=1:22
        data_and_stats=participant_BRCA(dataset).chrm(chr).data_and_stats_Tex; % based on Tex
        vpr1=data_and_stats.fitted_vprs(1,:);
        vpr2=data_and_stats.fitted_vprs(2,:);
        vpr3=data_and_stats.fitted_vprs(3,:);
        vpr4=data_and_stats.fitted_vprs(4,:);
             
        p1_05=data_and_stats.pv_ks_05(1,:);
        p2_05=data_and_stats.pv_ks_05(2,:);
        p3_05=data_and_stats.pv_ks_05(3,:);
        p4_05=data_and_stats.pv_ks_05(4,:);
        
        pNexTex=data_and_stats.pv_ks_ll(2,:);
        pNtrTtr=data_and_stats.pv_ks_ll(5,:);

        wl=data_and_stats.win_dp_length;
        
        vprs(dataset).vpr1=[vprs(dataset).vpr1 vpr1];
        vprs(dataset).vpr2=[vprs(dataset).vpr2 vpr2];
        vprs(dataset).vpr3=[vprs(dataset).vpr3 vpr3];
        vprs(dataset).vpr4=[vprs(dataset).vpr4 vpr4];
        
        vprs(dataset).p1_05=[vprs(dataset).p1_05 p1_05];
        vprs(dataset).p2_05=[vprs(dataset).p2_05 p2_05];
        vprs(dataset).p3_05=[vprs(dataset).p3_05 p3_05];
        vprs(dataset).p4_05=[vprs(dataset).p4_05 p4_05];
        
        vprs(dataset).pNexTex=[vprs(dataset).pNexTex pNexTex];
        vprs(dataset).pNtrTtr=[vprs(dataset).pNtrTtr pNtrTtr];

        vprs(dataset).wl=[vprs(dataset).wl wl];
    end
end
%
m=[];
for dataset=numel(vprs):-1:1
    m=vprs(dataset).vpr3;
    
    m(vprs(dataset).vpr1>0.58)=NaN; 
    m(vprs(dataset).p3_05>1e-5)=NaN;
    m(vprs(dataset).wl<50)=NaN;

    m(isnan(m))=[];
    vprs(dataset).m=m;
    if isempty(vprs(dataset).m)
        disp(dataset)
    end
end
disp('done')

%%
VBP_all=[];
all_prp=[];

parfor dataset=1:72
    single_dataset_vprs=vprs(dataset).m;
    single_dataset_vprs(single_dataset_vprs<0.58)=[];
    
    in_out=[];
    if ~isempty(single_dataset_vprs)
        disp(dataset)
        tic,
        for i_pop=2:4 %loop over populations
            for i_al=1:4 %loop over events
                nb_prop=size(all_vprs(i_pop,i_al).vprs,3);
                for prop=1:nb_prop %loop over considered mixtures
                    vprs_p_e=all_vprs(i_pop,i_al).vprs(:,:,prop);

                    vprs_p_e=vprs_p_e(:); %all the vprs of a given mixture
                    vprs_p_e(vprs_p_e<0.5)=[]; %folded
                    vprs_p_e=unique(vprs_p_e(:)); %unique vprs
                    [~,dist]=knnsearch(vprs_p_e,single_dataset_vprs');
                    if any(dist>0.009) || isempty(dist) %search radius
                        in_out(i_pop,i_al).prp(prop)=NaN;
                    else
                        in_out(i_pop,i_al).prp(prop)=1;
                    end
                end
            end
        end
        toc
        
        prp_vec=zeros(4*(156849+4851+99),3);
        ev_vec=zeros(4*(156849+4851+99),1);
        pop_vec=zeros(4*(156849+4851+99),1);
        tic,
        k=1;
        for i_pop=2:4 %loop over populations
            for i_al=1:4 %loop over events
                nb_prop=size(all_vprs(i_pop,i_al).vprs,3);
                for prop=1:nb_prop %loop over considered mixtures
                    if isnan(in_out(i_pop,i_al).prp(prop))
                        prp_vec(k,:)=NaN;
                        ev_vec(k)=NaN;
                        pop_vec(k)=NaN;
                    else
                        prp_vec(k,1:i_pop)=populations(1,i_pop).p(prop,:);
                        ev_vec(k)=i_al;
                        pop_vec(k)=i_pop;
                    end
                    k=k+1;
                end
            end
        end
        toc
        
        all_prp(dataset).prp_vec_Tex=prp_vec;
        all_prp(dataset).ev_vec_Tex=ev_vec;
        all_prp(dataset).pop_vec_Tex=pop_vec;
        
        tic
        for p_vec=2:4
            for e_vec=1:4
                if sum(ev_vec==e_vec & pop_vec==p_vec)>0
                    idx=(ev_vec==e_vec & pop_vec==p_vec);
                    pr=round(1-prp_vec(idx,1),2);
                    VBP_all(dataset).pur_dist(p_vec,e_vec).pr=pr;
                end
            end
        end
        toc

    else
        disp(['Dataset: ' num2str(dataset) ' no vprs'])
        VBP_all(dataset).pur_dist=[];
        all_prp(dataset).prp_vec_Tex=[];
        all_prp(dataset).ev_vec_Tex=[];
        all_prp(dataset).pop_vec_Tex=[];
    end
end

save -v7.3 all_prp all_prp
save -v7.3 VBP_all VBP_all