% TERNAXES create ternary axis
%   HOLD_STATE = TERNAXES(MAJORS) creates a ternary axis system using the system
%   defaults and with MAJORS major tickmarks.

% Author: Carl Sandrock 20050211

% To Do

% Modifications
% 20160405 (SA) Added lines to change the order/direction of axes (i.e.
%               clockwise or counter-clockwise) cooresponding to user-specified 
%               option on terncoords

% Modifiers
% (CS) Carl Sandrock
% (SA) Shahab Afshari

function [hold_state, cax, next] = ga_ternaxes(majors,tc1,tc2,tc3,fs)
if nargin < 1
    majors = 10;
end


% TODO: Handle these as options
direction = 'clockwise';
percentage = false;

%TODO: Get a better way of offsetting the labels
xoffset = 0.25;
yoffset = 0.01;

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
axis(cax,'equal')

if nargin < 4
    tc1=tc;
    tc2=tc;
    tc3=tc;
end
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');

set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
         'DefaultTextFontName',   get(cax, 'FontName'), ...
         'DefaultTextFontSize',   get(cax, 'FontSize'), ...
         'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
         'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state
	%plot axis lines
	hold on;
	plot ([0 1 0.5 0],[0 0 sin(1/3*pi) 0], 'color', tc, 'linewidth',1,...
                   'handlevisibility','off');
	set(gca, 'visible', 'off');

    % plot background if necessary
    if ~ischar(get(cax,'color'))
       patch('xdata', [0 1 0.5 0], 'ydata', [0 0 sin(1/3*pi) 0], ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end

	% Generate labels
    majorticks = linspace(0, 1, majors + 1); 
    majorticks_lbl=majorticks(2:end);
    majorticks = majorticks(1:end-1);

    if percentage
        multiplier = 100;
    else
        multiplier = 1;
    end
    
    if ~strcmp(direction, 'clockwise')
        labels = num2str(majorticks'*multiplier);
    else
        kk=1;
        for kkk=numel(majorticks):-1:1
        %labels = num2str(majorticks(end:-1:1)'*multiplier);
        labels{kk} = num2str(majorticks(kkk)*multiplier);
        kk=kk+1;
        end
    end
    
    zerocomp = zeros(size(majorticks)); % represents zero composition
    
    if ~strcmp(direction, 'clockwise')
	% Plot right labels (no c - only b a)
    [lxc, lyc] = terncoords(1-majorticks, majorticks, zerocomp);
    text(lxc, lyc, [repmat('  ', length(labels), 1) labels], 'color', tc1,'Fontsize',fs);

	% Plot bottom labels (no b - only a c)
    [lxb, lyb] = terncoords(majorticks, zerocomp, 1-majorticks); % fB = 1-fA
    text(lxb, lyb, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'Middle', 'color', tc2,'Fontsize',fs);

	% Plot left labels (no a, only c b)
	[lxa, lya] = terncoords(zerocomp, 1-majorticks, majorticks);
    text(lxa, lya, labels,'HorizontalAlignment', 'right', 'color', tc3,'Fontsize',fs);
    else
    % Plot right labels (no c - only b a)
    [lxc, lyc] = terncoords(1-majorticks, majorticks, zerocomp);
    [lxcl, lycl] = terncoords(1-majorticks_lbl, majorticks_lbl, zerocomp);

	%text(lxc+0.05, lyc-0.025, [repmat('  ', length(labels), 1) labels], 'color', tc2);
    text(lxcl+0.02, lycl-0.02, labels,'HorizontalAlignment', 'left', 'VerticalAlignment', 'Middle', 'color', tc2, 'Fontsize',fs,'Rotation',-60);

	% Plot bottom labels (no b - only a c)
    [lxb, lyb] = terncoords(majorticks, zerocomp, 1-majorticks); % fB = 1-fA
    [lxbl, lybl] = terncoords(majorticks_lbl, zerocomp, 1-majorticks_lbl); % fB = 1-fA

	%text(lxb-0.115, lyb-0.065, labels, 'VerticalAlignment', 'Top', 'color', tc1);
    text(lxbl-0.03, lybl, labels,'HorizontalAlignment', 'right','VerticalAlignment', 'Middle', 'color',  tc3, 'Fontsize',fs);

	% Plot left labels (no a, only c b)
	[lxa, lya] = terncoords(zerocomp, 1-majorticks, majorticks);
    [lxal, lyal] = terncoords(zerocomp, 1-majorticks_lbl, majorticks_lbl);

	%text(lxa-0.035, lya+0.09, labels, 'color', tc3);
    text(lxal+0.02, lyal+0.02, labels,'HorizontalAlignment', 'left','VerticalAlignment', 'Middle', 'color',  tc1, 'Fontsize',fs,'Rotation',60);
    end
	
	nlabels = length(labels)-1;
	for i = 1:nlabels
        plot([lxa(i+1) lxb(nlabels - i + 2)], [lya(i+1) lyb(nlabels - i + 2)], ls, 'color', tc3, 'linewidth',0.25,...
           'handlevisibility','off');
        plot([lxb(i+1) lxc(nlabels - i + 2)], [lyb(i+1) lyc(nlabels - i + 2)], ls, 'color', tc2, 'linewidth',0.25,...
           'handlevisibility','off');
        plot([lxc(i+1) lxa(nlabels - i + 2)], [lyc(i+1) lya(nlabels - i + 2)], ls, 'color', tc1, 'linewidth',0.25,...
           'handlevisibility','off');
	end;
end;

% Reset defaults
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits', fUnits );
