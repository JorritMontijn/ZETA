function [C,I] = findmax(A,maxnum)
%findmax Finds the highest n values from a vector
%   syntax: [C,I] = findmax(A,maxnum)
%	input:
%	- A: data vector
%	- maxnum: number of maximum values to be returned
%	output:
%	- C: Vector containing highest values
%	- I: Vector containing indices of highest values
%
%	Note: findmax is just an easy way to use max() consecutively to
%	find the n highest values. If maxnum is larger than the amount
%	of valid numbers, it will have a NaN for every unvalid number.
%
%
%	Version history:
%	1.0 - April 21 2011
%	Created by Jorrit Montijn

if min(size(A)) > 1
	warning([mfilename ':InputMultiDim'],'Input Error: findmax() only works with vectors, but input is a matrix; transforming to matrix')
	A = A(:);
end
if nargin < 2
	maxnum = numel(A);
end
C = nan(1,maxnum);
I = nan(1,maxnum);
maxcount = 0;
boolStop = false;

while ~boolStop
	maxcount = maxcount + 1;
	
	[tC,tI] = max(A);
	
	A(tI) = nan(1,length(tI));
	if max(isnan(tC)) == 1
		boolStop = true;
	else
		C(maxcount) = tC;
		I(maxcount) = tI;
	end
	if maxcount >= maxnum
		boolStop = true;
	end
end