function A = scale(O, factor)

%  Scales intervalSet timepoints by a factor
%  	
%  	USAGE:
%  	A = shift(A, shift, TimeUnits) 
%  	
%  	INPUTS:
%  	O - an intervalSet object
%  	shift - the time to shift the intervalSet of
%  
%   OUTPUTS:
%   A - scaled intervalSet object

% copyright (c) 2014 Brendon O. Watson
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  A = intervalSet(O.start * factor, O.stop * factor);
  
  