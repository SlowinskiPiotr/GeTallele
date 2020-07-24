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
%                      columns are:  
%                       1. number of populations considered
%                       2. number of events considered
%                       3. proportion of the normal population
%                       4 to 5 (or 6). proportions of the tumor populations
%                       (depends on the pop parameter)
%            pop     - number of populations for which admissible mixtures
%                      are suppoose to be plotted
%            ev      - maximal number of geneti events for which admissible mixtures
%                      are suppoose to be plotted
% 
% Function uses modified functions from the ternplot toolbox: https://zenodo.org/badge/latestdoi/31416657
% which is include in the \external\ternplot folder.
%
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
%
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

% split input matrix into smaller varibles
prp_vec=all_prp(:,3:4); %
ev_vec=all_prp(:,2);
pop_vec=all_prp(:,1);

% transform proportions into the ternary coordinate system
y_Tex = prp_vec(:,1)*sin(deg2rad(60)); % proportio of the normal population N
x_Tex = prp_vec(:,2) + y_Tex*cot(deg2rad(60)); % proportion of the 1st tumor population
% proportion of the 2nd (or 2nd + 3rd) tumor population follows from the
% fact that proportions sum to 1.

% color convention from the paper
clor=lines(7);
clor=clor([3 6 5],:);


% select proportions that correspond to mixture of interest (number of populations and maximal number of events)
idx=(ev_vec==ev & pop_vec==pop);

% generate colored ternary axes (modified function ternaxes.m from the ternplot toolbox)
ga_ternaxes(10,clor(1,:),clor(2,:),clor(3,:),8);
hold on,

% set labels of the axes depending on thenumber of populations (modified function ternlabel.m from the ternplot toolbox)
if pop==4
    ga_ternlabel('T2+T3','T1','N', clor(2,:),clor(1,:),clor(3,:),10,0.05)
else
    ga_ternlabel('T2','T1','N', clor(2,:),clor(1,:),clor(3,:),10,0.05)
end

% plot the individual admissible mixtures as dots
% size of the dots depends on the number of admissible mixtures
if sum(idx)<5 || pop==2 % big dots if they are on the edges or there is <5 of them 
    xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
    scatter(xy(:,1),xy(:,2),60,'k','filled','markerfacealpha',0.25)
elseif sum(idx)<20 && pop>2 % smaller dots if there is <20 of them 
    xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
    scatter(xy(:,1),xy(:,2),40,'k','filled','markerfacealpha',0.25)
else % small dots if there is >=20 of them 
    xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
    scatter(xy(:,1),xy(:,2),10,'k','filled','markerfacealpha',0.25)
end

% conservative purity estimate: the smallest value of 1 - proportion of the normal population
VBP=min(1-prp_vec(idx,1));


% display plot title that inclueds purity estimate
if ~isempty(VBP)
    text(0,0.8,{'Admissible mixtures for:',[num2str(pop) ' populations,'],[num2str(ev) ' events.'],['Purity:  ' num2str(100*VBP) '%']},...
    'Fontsize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle')
else
    text(0,0.8,{'NO admissible mixtures for:',[num2str(pop) ' populations,'],[num2str(ev) ' events.']},...
    'Fontsize',12,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle','Color','Red')
end
hold off