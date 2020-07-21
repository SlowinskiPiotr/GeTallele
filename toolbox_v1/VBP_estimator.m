function [VBP_est_keep,vbp_data_keep,VBP_dist,all_prp]=VBP_estimator(chr_results,all_vprs,populations)
% collate all vprs and disregard some

all_vpr_Nex=[];
all_vpr_Tex=[];
all_Tex_p_05=[];
all_seg_dp=[];

for chr_idx=1:22
    data_and_stats_Tex=chr_results(chr_idx).Tex; % based on Tex
    data_and_stats_Nex=chr_results(chr_idx).Nex; % fitted in Nex
    
    vpr_Nex=data_and_stats_Nex.fitted_vprs;
    
    vpr_Tex=data_and_stats_Tex.fitted_vprs;
    p_05=data_and_stats_Tex.pv_vpr_neq_05;
    seg_dp=data_and_stats_Tex.seg_dp_length;
    
    all_vpr_Nex=[all_vpr_Nex vpr_Nex];
    all_vpr_Tex=[all_vpr_Tex vpr_Tex];
    all_Tex_p_05=[all_Tex_p_05 p_05];
    
    all_seg_dp=[all_seg_dp seg_dp];
end

vprs_for_VBP=all_vpr_Tex;

vprs_for_VBP(all_vpr_Nex>0.58)=NaN;
vprs_for_VBP(all_Tex_p_05>1e-6)=NaN;
vprs_for_VBP(all_seg_dp<50)=NaN;

vprs_for_VBP(isnan(vprs_for_VBP))=[];
vprs_for_VBP(vprs_for_VBP<0.58)=[];

n_pop=size(all_vprs,1);
n_ev=size(all_vprs,2);

if ~isempty(vprs_for_VBP)
    tic,
    for i_pop=n_pop:-1:2 %loop over populations
        for i_ev=n_ev:-1:1 %loop over events
            nb_prop=size(all_vprs(i_pop,i_ev).vprs,3);
            for prop=1:nb_prop %loop over considered mixtures
                vprs_p_e=all_vprs(i_pop,i_ev).vprs(:,:,prop);
                
                vprs_p_e=vprs_p_e(:); %all the vprs of a given mixture
                vprs_p_e(vprs_p_e<0.5)=[]; %folded
                vprs_p_e=unique(vprs_p_e(:)); %unique vprs
                [~,dist]=knnsearch(vprs_p_e,vprs_for_VBP');
                if any(dist>0.009) || isempty(dist) %search radius
                    in_out(i_pop,i_ev).prp(prop)=NaN;
                else
                    in_out(i_pop,i_ev).prp(prop)=1;
                end
            end
        end
    end
    toc
    
    prp_times_ev_size=0;
    for i_pop=2:n_pop %loop over populations
        nb_prop=size(all_vprs(i_pop,2).vprs,3);
        prp_times_ev_size=prp_times_ev_size+n_ev*nb_prop;
    end
    
    prp_vec=zeros(prp_times_ev_size,n_pop);
    ev_vec=zeros(prp_times_ev_size,1);
    pop_vec=zeros(prp_times_ev_size,1);
    tic,
    k=1;
    for i_pop=2:n_pop %loop over populations
        for i_ev=1:n_ev %loop over events
            nb_prop=size(all_vprs(i_pop,i_ev).vprs,3);
            for prop=1:nb_prop %loop over considered mixtures
                if isnan(in_out(i_pop,i_ev).prp(prop))
                    prp_vec(k,:)=NaN;
                    ev_vec(k)=NaN;
                    pop_vec(k)=NaN;
                else
                    prp_vec(k,1:i_pop)=populations(1,i_pop).p(prop,:);
                    ev_vec(k)=i_ev;
                    pop_vec(k)=i_pop;
                end
                k=k+1;
            end
        end
    end
    toc
    
    all_prp.prp_vec=prp_vec;
    all_prp.ev_vec=ev_vec;
    all_prp.pop_vec=pop_vec;
    
    tic
    for p_vec=n_pop:-1:2
        for e_vec=n_ev:-1:1
            if sum(ev_vec==e_vec & pop_vec==p_vec)>0
                idx=(ev_vec==e_vec & pop_vec==p_vec);
                pr=round(1-prp_vec(idx,1),2);
                VBP_dist(p_vec,e_vec).pr=pr;
            else
                
                VBP_dist(p_vec,e_vec).pr=[];
                all_prp.prp_vec=[];
                all_prp.ev_vec=[];
                all_prp.pop_vec=[];
            end
        end
    end
    toc
    
else
    VBP_dist=[];
    all_prp.prp_vec=[];
    all_prp.ev_vec=[];
    all_prp.pop_vec=[];
end

if ~isempty(VBP_dist)
    to_keep=0;
    for p_i=2:n_pop
        for e_i=1:n_ev
            if ~isempty(VBP_dist(p_i,e_i).pr)
                vbp_data=VBP_dist(p_i,e_i).pr;
                VBP_est=min(vbp_data);
                disp('mixture found!')
                disp(['populations: ' num2str(p_i) '; events: ' num2str(e_i)]);
                disp(['Purity: ' num2str(VBP_est)]);
                if to_keep==0
                    disp('this mixture has the lowest complexity and this value will be returned');
                    vbp_data_keep=vbp_data;
                    VBP_est_keep=VBP_est;
                    to_keep=1; %this if statement will be entered only once
                end
            else
                disp(['no fit - populations: ' num2str(p_i) '; events: ' num2str(e_i)])
                vbp_data_keep=NaN;
                VBP_est_keep=NaN;
            end
        end
    end
else
    disp('no vprs')
    vbp_data_keep=NaN;
    VBP_est_keep=NaN;
end