function plot_VAF_and_vpr(results,colr,text_offset)
% This function plot the VAF values together with fitted variant
% probabilities.
%
% Call function: plot_VAF_and_vpr(results,colr,text_offset)
%
%      Input: 
%            results - structure from the functions:
%            find_segments_and_fit_vpr or fit_vpr_in_segments
%            colr        - color for plotting points and lines
%            text_offset - optional control of vertical postion of the displayed Vpr values 
% 
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.       
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

if nargin==2
   text_offset=0;
end

VAFs=results.data(:,3);
bp_pos=results.data(:,2);
vprs=results.fitted_vprs;
no_different_from05=results.pv_vpr_neq_05;

segment_edges=results.segment_edges;
idx_start=segment_edges(1:end-1);
idx_end=[segment_edges(2:end-1)-1 segment_edges(end)];


plot(bp_pos,VAFs,'o','MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',2)
hold on

for i=1:numel(segment_edges)-1   
    plot([bp_pos(idx_start(i)) bp_pos(idx_end(i))],[vprs(i) vprs(i)],'-','Color',colr,'linewidth',3)
    plot([bp_pos(idx_start(i)) bp_pos(idx_end(i))],1-[vprs(i) vprs(i)],'-','Color',colr,'linewidth',3)
    if no_different_from05(i)>1e-5
        text(bp_pos(idx_start(i))+(bp_pos(idx_end(i))-bp_pos(idx_start(i)))/2,0.5+text_offset/10,'KS: 0.5, (p>1e-5)',...
            'Color',colr,'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle')
    else
        text(bp_pos(idx_start(i))+(bp_pos(idx_end(i))-bp_pos(idx_start(i)))/2,vprs(i)+text_offset,num2str(vprs(i),3),...
            'Color',colr,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold')
    end
end
hold off
end