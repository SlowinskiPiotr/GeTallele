function [VBP_est,vbp_data,VBP_dist,all_prp]=VBP_estimator(chr_results,all_vprs,populations)
% collate all vprs and remove some

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

in_out=[];
if ~isempty(vprs_for_VBP)
    tic,
    for i_pop=2:4 %loop over populations
        for i_al=1:4 %loop over events
            nb_prop=size(all_vprs(i_pop,i_al).vprs,3);
            for prop=1:nb_prop %loop over considered mixtures
                vprs_p_e=all_vprs(i_pop,i_al).vprs(:,:,prop);
                
                vprs_p_e=vprs_p_e(:); %all the vprs of a given mixture
                vprs_p_e(vprs_p_e<0.5)=[]; %folded
                vprs_p_e=unique(vprs_p_e(:)); %unique vprs
                [~,dist]=knnsearch(vprs_p_e,vprs_for_VBP');
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
    
    all_prp.prp_vec=prp_vec;
    all_prp.ev_vec=ev_vec;
    all_prp.pop_vec=pop_vec;
    
    tic
    for p_vec=2:4
        for e_vec=1:4
            if sum(ev_vec==e_vec & pop_vec==p_vec)>0
                idx=(ev_vec==e_vec & pop_vec==p_vec);
                pr=round(1-prp_vec(idx,1),2);
                VBP_dist(p_vec,e_vec).pr=pr;
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
    for p_i=2:4
        for e_i=1:4
            if ~isempty(VBP_dist(p_i,e_i).pr)
                vbp_data=VBP_dist(p_i,e_i).pr;
                VBP_est=min(vbp_data);
                break
            end
        end
    end
else
    disp('no vprs')
	vbp_data=NaN;
	VBP_est=NaN;  
end