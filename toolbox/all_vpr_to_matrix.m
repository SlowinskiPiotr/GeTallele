function out=all_vpr_to_matrix(chrm,CNA,swt)
%chr start_pos end_pos mode_l1 mode_l2 mode_l3 mode_l4 ...
k=1;
switch swt
    case 'Tex'
        for chr=1:22
            nb_modes=size(chrm(chr).data_and_stats_Tex.fitted_vprs,2);
            
            CNAchr=CNA(CNA(:,1)==chr,:);
            [pos_CNA,ims_CNA]=unique([CNAchr(:,2); CNAchr(:,3)+eps]);
            all_CNA=[CNAchr(:,5); CNAchr(:,5)];
            all_CNA=all_CNA(ims_CNA);
            
            for m=1:nb_modes
                win_idx=chrm(chr).data_and_stats_Tex.vpr_data(m).window;
                
                s_window=sort(win_idx,'ascend');
                idx_jump=find(diff(s_window)>1);
                
                nb_jumps=numel(idx_jump);
                win_edges_val=[];
                
                if nb_jumps>0
                    win_edges_idx=sort([s_window(1) s_window(idx_jump) s_window(idx_jump+1)  s_window(end)]);
                    win_edges_val=chrm(chr).data_and_stats_Tex.segment(win_edges_idx,1);
                else
                    win_edges_idx=[s_window(1) s_window(end)];
                    win_edges_val=chrm(chr).data_and_stats_Tex.segment(win_edges_idx,1);
                end
               
                
                for jmp=1:2:numel(win_edges_val)
                    out(k,1)=chr;
                    out(k,2)=win_edges_val(jmp);
                    out(k,3)=win_edges_val(jmp+1);
                    out(k,4:7)=chrm(chr).data_and_stats_Tex.fitted_vprs(:,m);
                    out(k,8:11)=chrm(chr).data_and_stats_Tex.pv_ks_05(:,m);
                    out(k,12:17)=chrm(chr).data_and_stats_Tex.pv_ks_ll(:,m);
                    
                    mds_pos=chrm(chr).data_and_stats_Tex.segment(win_edges_idx(jmp):win_edges_idx(jmp+1),1);
                    
                    interp_mds_CNA=interp1(pos_CNA,all_CNA,mds_pos,'nearest','extrap');
                    out(k,18)=median(interp_mds_CNA);
                    out(k,19)=numel(mds_pos);
                    out(k,20)=win_edges_idx(jmp);
                    k=k+1;
                end
            end
        end
    case 'Ttr'
        for chr=1:22
            nb_modes=size(chrm(chr).data_and_stats_Ttr.fitted_vprs,2);
            
            CNAchr=CNA(CNA(:,1)==chr,:);
            [pos_CNA,ims_CNA]=unique([CNAchr(:,2); CNAchr(:,3)+eps]);
            all_CNA=[CNAchr(:,5); CNAchr(:,5)];
            all_CNA=all_CNA(ims_CNA);
            
            for m=1:nb_modes
                win_idx=chrm(chr).data_and_stats_Ttr.vpr_data(m).window;
                
                s_window=sort(win_idx,'ascend');
                idx_jump=find(diff(s_window)>1);
                
                nb_jumps=numel(idx_jump);
                win_edges_val=[];
                
                if nb_jumps>0
                    win_edges_idx=sort([s_window(1) s_window(idx_jump) s_window(idx_jump+1)  s_window(end)]);
                    win_edges_val=chrm(chr).data_and_stats_Ttr.segment(win_edges_idx,1);
                else
                    win_edges_idx=[s_window(1) s_window(end)];
                    win_edges_val=chrm(chr).data_and_stats_Ttr.segment(win_edges_idx,1);
                end
                
                
                
                for jmp=1:2:numel(win_edges_val)
                    out(k,1)=chr;
                    out(k,2)=win_edges_val(jmp);
                    out(k,3)=win_edges_val(jmp+1);
                    out(k,4:7)=chrm(chr).data_and_stats_Ttr.fitted_vprs(:,m);
                    out(k,8:11)=chrm(chr).data_and_stats_Ttr.pv_ks_05(:,m);
                    out(k,12:17)=chrm(chr).data_and_stats_Ttr.pv_ks_ll(:,m);
                    
                    mds_pos=chrm(chr).data_and_stats_Ttr.segment(win_edges_idx(jmp):win_edges_idx(jmp+1),1);
                    
                    interp_mds_CNA=interp1(pos_CNA,all_CNA,mds_pos,'nearest','extrap');
                    out(k,18)=median(interp_mds_CNA);
                    out(k,19)=numel(mds_pos);
                    out(k,20)=win_edges_idx(jmp);

                    k=k+1;
                end
            end
        end
end
