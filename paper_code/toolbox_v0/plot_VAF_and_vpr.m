function plot_VAF_and_vpr(data,layer)
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
%
% input:
%   data - data structure from the function find_modes or compare modes
%   layer - single layer or two layers that should be ploted (convention for 4 layers: Nex=1, Ntr=2, Tex=3, Ttr=4
%                          convention for 2 layers: Ex=1, Tr=2)
%
% If you have any questions please contact: p.m.slowinski@exeter.ac.uk

nb_of_mwindows=numel(data.vpr_data);
nb_non_nan=sum(~isnan(data.fitted_vprs(:,1)));

% set some vprs to 0.5 by hand
data.fitted_vprs(data.fitted_vprs<0.58)=0.5;
data.fitted_vprs(data.pv_ks_05>1e-5)=0.5;

data.fitted_vprs,
if numel(layer)==1 %plot single layer (raw data + fitted modes)
    colr=lines(nb_non_nan); % up to 7 different colors
    k=1;
    for i=1:nb_of_mwindows
        s_data=data.vpr_data(i).win_data; % original data in the merged window in 5 x number of points matrix
        win_idx=data.vpr_data(i).window;
        if ~isempty(s_data)
            plot(s_data(:,1),s_data(:,layer+1),'o','color',colr(k,:),'markerfacecolor',colr(k,:))
            hold on
            s_window=sort(win_idx,'ascend');
            idx_jump=find(diff(s_window)>1);
            nb_jumps=numel(idx_jump);
            win_edges=[];
            if nb_jumps>0
                win_edges_idx=[s_window(1) s_window(idx_jump) s_window(idx_jump+1)  s_window(end)];
                win_edges_val=data.segment(sort(win_edges_idx),1);
                n_edges=3*(1+nb_jumps);
                win_edges(1:3:n_edges)=win_edges_val(1:2:end);
                win_edges(2:3:n_edges)=win_edges_val(2:2:end);
                win_edges(3:3:n_edges)=NaN;
            else
                win_edges=data.segment([s_window(1) s_window(end)],1);
            end
            plot(win_edges,ones(size(win_edges)).*data.fitted_vprs(layer,i),'-','color',colr(k,:),'linewidth',2)
            plot(win_edges,ones(size(win_edges)).*(1-data.fitted_vprs(layer,i)),'-','color',colr(k,:),'linewidth',2)
            k=k+1;
        end
    end
    plot([data.segment(1,1) data.segment(end,1)],[0.5 0.5],'k:')

    idx_fitted_vprs=find(~isnan(data.fitted_vprs(layer,:)));
    [s_fitted_modes,idx_sfm]=sort(data.fitted_vprs(layer,idx_fitted_vprs));
    if isfield(data,'pv_ks_05')
        s_sgnfc=data.pv_ks_05(1,idx_fitted_vprs);
        s_sgnfc=s_sgnfc(idx_sfm);
        s_sgnfc(s_sgnfc<1e-5)=0;
    else
        s_sgnfc=ones(1,numel(idx_fitted_vprs));
    end
    
    mode_str=[];
    for kk=1:numel(idx_fitted_vprs)
        mode_str=[mode_str num2str(s_fitted_modes(kk)) ', ' num2str(s_sgnfc(kk),1) '; '];
    end
    title(['Chrm ' num2str(data.chr) '; L: ' num2str(layer) '; Mds: ' mode_str])
    ylim([0 1])
    xlim([data.segment(1,1) data.segment(end,1)])
    set(gca,'XTickLabel',[])
    hold off
    
else %plot two layers
    colr=lines(7);
    colr=colr([1 6 2 3],:);
    colr=colr([layer(1) layer(2)],:);
    
    for i=1:nb_of_mwindows
        win_idx=data.vpr_data(i).window;
        if ~isempty(win_idx)
            s_window=sort(win_idx,'ascend');
            idx_jump=find(diff(s_window)>1);
            nb_jumps=numel(idx_jump);
            win_edges=[];
            if nb_jumps>0
                win_edges_idx=[s_window(1) s_window(idx_jump) s_window(idx_jump+1)  s_window(end)];
                win_edges_val=data.segment(sort(win_edges_idx),1);
                n_edges=3*(1+nb_jumps);
                win_edges(1:3:n_edges)=win_edges_val(1:2:end);
                win_edges(2:3:n_edges)=win_edges_val(2:2:end);
                win_edges(3:3:n_edges)=NaN;
            else
                win_edges=data.segment([s_window(1) s_window(end)],1);
            end
            
           
                plot(win_edges,ones(size(win_edges)).*data.fitted_vprs(layer(1),i),'-','color',colr(1,:),'linewidth',1.5)
                hold on
                plot(win_edges,ones(size(win_edges)).*(1-data.fitted_vprs(layer(1),i)),'-','color',colr(1,:),'linewidth',1.5)
                plot(win_edges,ones(size(win_edges)).*data.fitted_vprs(layer(2),i),'-','color',colr(2,:),'linewidth',1.5)
                plot(win_edges,ones(size(win_edges)).*(1-data.fitted_vprs(layer(2),i)),'-','color',colr(2,:),'linewidth',1.5)
        end
    end
    plot([data.segment(1,1) data.segment(end,1)],[0.5 0.5],'k:')  
    
    plot(data.segment(:,1),data.segment(:,layer(1)+1),'o','color',colr(1,:),'markerfacecolor',colr(1,:),'markersize',2)
    plot(data.segment(:,1),data.segment(:,layer(2)+1),'o','color',colr(2,:),'markerfacecolor',colr(2,:),'markersize',2)

    ylim([0 1])
    xlim([data.segment(1,1) data.segment(end,1)])
    set(gca,'XTickLabel',[])
    hold off
    %end
end

end