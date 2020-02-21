
%%
% for each dataset collect modes from all the chromosomes

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
        
        % based on Tex
        data_and_stats=participant_BRCA(dataset).chrm(chr).data_and_stats_Tex;
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
tic,
[all_vprs,populations,alleles]=gen_all_possible_vprs(4,4,0.01);
toc,
save -v7.3 vprs_p4e4 all_vprs populations alleles
%%
VBP_all=[];
all_prp=[];

parfor dataset=1:72
    mds=vprs(dataset).m;%unique(round(mod(mod>0.5),2));
    mds(mds<0.58)=[];
    
    in_out=[];
    %clf
    if ~isempty(mds)
        disp(dataset)
        tic,
        for i_pop=2:4 %loop over populations
            for i_al=1:4 %loop over events
                nb_prop=size(all_vprs(i_pop,i_al).dict,3);
                for prop=1:nb_prop %loop over considered admixtures
                    dict_mat=all_vprs(i_pop,i_al).dict(:,:,prop);
                    if i_pop==3
                    if i_al==1
                        dict_mat(1:2,1:2)=0;
                    elseif i_al==2
                        dict_mat(1:3,1:3)=0;
                    elseif i_al==3
                        dict_mat(1:4,1:4)=0;
                    elseif i_al==4
                        dict_mat(1:5,1:5)=0;
                    end
                    elseif i_pop==4
                    if i_al==1
                        dict_mat(1:3,1:3)=0;
                    elseif i_al==2
                        dict_mat(1:6,1:6)=0;
                    elseif i_al==3
                        dict_mat(1:10,1:10)=0;
                    elseif i_al==4
                        dict_mat(1:15,1:15)=0;
                    end
                    end
                    modes_dict=dict_mat(:); %all the modes of a given admixture
                    modes_dict(modes_dict<0.5)=[]; %folded
                    modes_dict=unique(modes_dict(:)); %unique modes
                    [~,dist]=knnsearch(modes_dict,mds');
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
                nb_prop=size(all_vprs(i_pop,i_al).dict,3);
                for prop=1:nb_prop %loop over considered admixtures
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
        
        tic
        all_prp(dataset).prp_vec_Tex=prp_vec;
        all_prp(dataset).ev_vec_Tex=ev_vec;
        all_prp(dataset).pop_vec_Tex=pop_vec;
        toc
        
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
save VBP_all VBP_all