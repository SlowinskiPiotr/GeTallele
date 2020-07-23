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
%            text_offset - optional control of vertical postion of the
%                          displayed estimated vpr values. If the estimated
%                          if the p-value of the Kolmogorov-Smirnov test
%                          >1e-5, when the observed sample cannot be
%                          differentiated from a synthetic VAF sample with
%                          vpr=0.5 the plot displays the following text
%                          instead of the estimated vpr value 
%                          'KS: 0.5, (p>1e-5)'.
%                        
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.       
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

% check if there is the text offset is set
if nargin==2
   text_offset=0;
end


% copy values from structure into vectors
VAFs=results.data(:,3); %observed VAF values
bp_pos=results.data(:,2); %base-pairs coordinates
vprs=results.fitted_vprs; %estimated vpr values
no_different_from05=results.pv_vpr_neq_05; %p-values of the Kolmogorov-Smirnov test
segment_edges=results.segment_edges; %indices of the change points

idx_start=segment_edges(1:end-1); %make indices of the start of the segments
idx_end=[segment_edges(2:end-1)-1 segment_edges(end)]; %make indices of the ends of the segments

% plot the VAF values as dots
plot(bp_pos,VAFs,'o','MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',2)
hold on

%for each segments
for i=1:numel(segment_edges)-1   
    % plot the vprs estimates as lines
    plot([bp_pos(idx_start(i)) bp_pos(idx_end(i))],[vprs(i) vprs(i)],'-','Color',colr,'linewidth',3)
    plot([bp_pos(idx_start(i)) bp_pos(idx_end(i))],1-[vprs(i) vprs(i)],'-','Color',colr,'linewidth',3)
    
    % add the text with the vpr value
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