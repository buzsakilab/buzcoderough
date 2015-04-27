function [ang,GoodRanges] = HeadDirection_Wrapper(fbasename)


% loads head-direction from a position file file (ending in .whl)
%
% USAGE
%     [ang,GoodRanges] = HeadDirectionWhl(fbasename)
%     
% INPUT:
%     fbasename: session file basename
%	
% OUTPUT
%     ang: a tsd object of angular values
%     GoodRanges: a intervalSet object where LEDs were successfully detected

% Adrien Peyrache 2011
  
[whl,t] = LoadPosition(fbasename);

t = t*10000;
dx = whl(:,1)-whl(:,3);
dy = whl(:,2)-whl(:,4);

ep = ~(isnan(dx) | isnan(dy));

tg = t;
tg(~ep) = -1;
tg(end) = -1;
GoodRanges = thresholdIntervals(tsd(t,tg),0,'Direction','Above');

ang = mod(atan2(dy(ep),dx(ep)),2*pi);
ang = tsd(t(ep),ang);


