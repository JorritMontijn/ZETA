function [bcR, p, T, df] = bcdistcorr2(x, y)
	% BCDISTCORR computes the bias corrected distance correlation
	%
	%   [BCR,P,T,DF] = BCDISTCORR(X,Y), where X (size n-by-p) and Y (size
	%   n-by-q) are n random samples of observation.  The function returns the
	%   bias corrected distance correlation BCR and the corresponding p-value
	%   P, as well as the student t statistics T and its degree of freedom DF.
	%   Note that the The t-test of independence is unbiased for every n ? 4
	%   and any significance level.
	%
	%   This implementation is based on Székely, G. J., & Rizzo, M. L. (2013).
	%   The distance correlation t-test of independence in high dimension.
	%   Journal of Multivariate Analysis, 117, 193-213.
	%
	%   Date: 7.30.2016
	%   Author: Po-He Tseng (pohetsn@gmail.com)
	%   Date: 9 July 2019
	%   Rewrote A* calculation to be faster for large matrices (by Jorrit Montijn)
	
	assert(rows(x)==rows(y));
	n = rows(x);
	X = Astar(x);
	Y = Astar(y);
	XY = modified_distance_covariance(X, Y);
	XX = modified_distance_covariance(X, X);
	clear X;
	YY = modified_distance_covariance(Y, Y);
	clear Y;
	bcR = XY/sqrt(XX*YY);
	M = n*(n-3)/2;
	T = sqrt(M-1) * bcR / sqrt(1-bcR^2);
	df = M-1;
	p = 1 - tcdf(T, df);
	%fprintf('bias-corrected R = %.3f, p-value=%.3f, T(%d)=%.4f\n',...
	%	bcR, p, df, T);
end
function XY = modified_distance_covariance(X, Y)
	n = rows(X);
	XY = sum(sum(bsxfun(@times, X, Y)))...
		- (n/(n-2))*diag(X)'*diag(Y);
end
function matA = Astar(x)
	%use loop?
	boolNew = true;
	
	%get inputs
	matD = pdist2(x,x); %[n x n]
	n = rows(x);
	vecM = mean(matD); %[1 x n]
	vecTranspM = vecM'; %[n x 1]
	dblM = mean(matD(:)); %[1 x 1]
	
	%calculate A*
	if boolNew
		%loops are actually faster
		matA = ones(n,n);
		for i=1:n
			matA(:,i) = matD(:,i) - vecTranspM;
		end
		
		% A = A - ones(n,1)*m;
		for i=1:n
			matA(i,:) = matA(i,:) - vecM;
		end
		
		matA = matA + dblM;
	else
		%matrix operations require more memory
		matA = bsxfun(@minus, matD, bsxfun(@mtimes, vecTranspM, ones(1,n)));
		matA = bsxfun(@minus, matA, bsxfun(@mtimes, ones(n,1), vecM));
		matA = bsxfun(@plus, matA, dblM);
	end
	% A = A + M;
	
	matA = matA - matD/n;
	matA(1:n+1:end) = vecM - dblM;
	matA = (n/(n-1))*matA;
end
function r = rows(x)
	r = size(x,1);
end