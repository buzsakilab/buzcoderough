function [TSO,varargout] = oneSeries(S)

% Sorts all the points in each element in the tsdArray, and returns a single tsd with all the points 
%  
%  	USAGE:
%  	tso = oneSeries(tsa) 
%   [tso,spikeids] = oneSeries(tsa)... gives a second output spikeids which are like a clu file: cluster ID of origin of each spike 

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


  sz = size(S.C{1});
  sz = sz(2:end);
  if sz == 0
      sz = 1;
  end
  
  for i = 2:length(S.C)
    good = 0;
    try
      sz2 = size(S.C{i});
      sz2 = sz2(2:end);
      if sum(sz2) == 0 
          S.C{i} = ts(zeros([0 sz]));
          sz2 = size(S.C{i});
          sz2 = sz2(2:end);
      end
      
      if all(sz2 == sz) 
        good = 1;
      end
    catch
      ;
    end
    if ~good
      error(['incompatible data size  if nargout == 2 with tsdArray: ' num2str(i)]);
    end
  end
  
  
  d = Data(S.C{1});
  t = Range(S.C{1});
  cellids = 1+zeros(size(S.C{1}));
  
  for i = 2:length(S.C)
    t = [t; Range(S.C{i})];
    d  = cat(1, d, Data(S.C{i}));
    cellids = cat(1,cellids,i+zeros(size(S.C{i})));
  end
  
  [t, ix] = sort(t);
  
  d = SelectAlongFirstDimension(d, ix);
  
  TSO = tsd(t, d);

  if nargout == 2;
   varargout{1} =  cellids(ix);
  end
  
      