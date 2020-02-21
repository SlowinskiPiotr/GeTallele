% TERNLABEL label ternary phase diagram
%   TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') labels a ternary phase diagram created using TERNPLOT
%   
%   H = TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') returns handles to the text objects created.
%   with the labels provided.  TeX escape codes are accepted.
%
%   See also TERNPLOT

% Author: Carl Sandrock 20020827

% To Do

% Modifications

% Modifiers

function h = ga_ternlabel(A, B, C, tc1, tc2, tc3, fs, offset)
if nargin < 6
    tc1='k';
    tc2='k';
    tc3='k';
end
r(1) = text(0.5, -0.05-offset/1.3, A, 'horizontalalignment', 'center','color',tc1,'FontSize',fs);
r(2) = text(1-0.45*sin(deg2rad(30))+offset, 0.5, B, 'rotation', -60, 'horizontalalignment', 'center','color',tc2,'FontSize',fs);
r(3) = text(0.45*sin(deg2rad(30))-offset, 0.5, C, 'rotation', 60, 'horizontalalignment', 'center','color',tc3,'FontSize',fs);

if nargout > 0
    h = r;
end;