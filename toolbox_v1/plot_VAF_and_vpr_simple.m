function plot_VAF_and_vpr_simple(results,colr)
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
%
% input:
%   results - data structure from the function find_segments_and_fit_vpr or compare modes
%
% If you have any questions please contact: p.m.slowinski@exeter.ac.uk

VAFs=results.data(:,3);
bp_pos=results.data(:,2);
vprs=results.fitted_vprs;
no_different_from05=results.pv_vpr_neq_05;

vprs(no_different_from05>1e-5)=0.5;

segment_edges=results.segment_edges;
idx_start=segment_edges(1:end-1);
idx_end=[segment_edges(2:end-1)-1 segment_edges(end)];


plot(bp_pos,VAFs,'o','MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',2)
hold on

for i=1:numel(segment_edges)-1   
    plot([bp_pos(idx_start(i)) bp_pos(idx_end(i))],[vprs(i) vprs(i)],'-','Color',colr,'linewidth',3)
    plot([bp_pos(idx_start(i)) bp_pos(idx_end(i))],1-[vprs(i) vprs(i)],'-','Color',colr,'linewidth',3)
end
hold off

end