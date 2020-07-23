function plot_adm_mix(all_prp,pop,ev)
% This function generates ternary plot showing all the admissible mixtures
% for a given number of population and with given maximal number of genetic
% events.
%
% Call function: plot_adm_mix(all_prp,pop,ev)
%
%      Input: 
%            all_prp - matrix with all the proportions (for any number of
%                      populations and events), returned by VBP_estimator function.
%            pop     - number of populations for which admissible mixtures
%                      are suppoose to be plotted
%            ev      - maximal number of geneti events for which admissible mixtures
%                      are suppoose to be plotted
% 
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
%
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

prp_vec=all_prp(:,3:4);
ev_vec=all_prp(:,2);
pop_vec=all_prp(:,1);

y_Tex = prp_vec(:,1)*sin(deg2rad(60));
x_Tex = prp_vec(:,2) + y_Tex*cot(deg2rad(60));

clor=lines(7);
clor=clor([3 6 5],:);

idx=(ev_vec==ev & pop_vec==pop);


ga_ternaxes(10,clor(1,:),clor(2,:),clor(3,:),8);
hold on,

if pop==4
    ga_ternlabel('T2+T3','T1','N', clor(2,:),clor(1,:),clor(3,:),10,0.05)
else
    ga_ternlabel('T2','T1','N', clor(2,:),clor(1,:),clor(3,:),10,0.05)
end

if sum(idx)<5 || pop==2
    xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
    scatter(xy(:,1),xy(:,2),60,'k','filled','markerfacealpha',0.25)
elseif sum(idx)<20 && pop>2
    xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
    scatter(xy(:,1),xy(:,2),40,'k','filled','markerfacealpha',0.25)
else
    xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
    scatter(xy(:,1),xy(:,2),10,'k','filled','markerfacealpha',0.25)
end


VBP=min(1-prp_vec(idx,1));

if ~isempty(VBP)
    text(0,0.8,{'Admissible mixtures for:',[num2str(pop) ' populations,'],[num2str(ev) ' events.'],['Purity:  ' num2str(100*VBP) '%']},...
    'Fontsize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle')
else
    text(0,0.8,{'NO admissible mixtures for:',[num2str(pop) ' populations,'],[num2str(ev) ' events.']},...
    'Fontsize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle','Color','Red')
end
hold off