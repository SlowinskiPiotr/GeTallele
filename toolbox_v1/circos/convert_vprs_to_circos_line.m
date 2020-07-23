function convert_vprs_to_circos_line(chr_results,loc_half,loc_width,top_or_bottom,fname)

fid=fopen(fname,'w');

for chrm=1:22
    
    chr_data=chr_results(chrm).results;
    seg_idx=chr_data.segment_edges;
    idx_start=seg_idx(1:end-1);
    idx_end=[seg_idx(2:end-1)+1 seg_idx(end)];

    for data_pts=1:numel(chr_data.fitted_vprs)
        
        start_pos=chr_data.data(idx_start(data_pts),2);
        end_pos=chr_data.data(idx_end(data_pts),2);
        value=chr_data.fitted_vprs(data_pts);
        p_ks_nq05=chr_data.pv_vpr_neq_05(data_pts);

        if value<=0.58
           value=0.5;
        end
        
        if p_ks_nq05>1e-5
           value=0.5;
        end
        
        if top_or_bottom==1
            value=loc_half+loc_width*(value-0.5)/0.5;
        else
            value=loc_half-loc_width*(value-0.5)/0.5;
        end
        
        datastr=['hs' num2str(chrm) ' ' num2str(start_pos) ' ' num2str(end_pos) ' ' num2str(value) ' ' ...
            'r0=' num2str(value,5) 'r-3p,' 'r1=' num2str(value,5) 'r+3p' '\n'];
        fprintf(fid,datastr);
    end
    
end
fclose(fid);
