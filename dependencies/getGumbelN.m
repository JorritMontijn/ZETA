function [dblP,dblZ] = getGumbelN(intN,dblX)
	%getGumbelN Calculate p-value and z-score for maximum value of N samples drawn from Gaussian
	%   [dblP,dblZ] = getGumbelN(intN,dblX)
	
	%syntax:  [dblP,dblZ] = getGumbel(intN,dblX)
	%	input:
	%	- intN: sample size
	%	- dblX: maximum value in sample
	%
	%	output:
	%	- dblP; p-value for dblX (chance that sample originates from Gaussian with mean 0 & variance 1
	%	- dblZ; z-score corresponding to P
	%
	%Note: this function assumes the underlying distribution is mean-zero
	%and variance 1; if this is not the case, you should rescale your values first
	%
	%Version history:
	%1.0 - March 11 2020
	%	Created by Jorrit Montijn
	%
	%Sources:
	%Baglivo (2005), ISBN: 9780898715668
	%Elfving (1947), https://doi.org/10.1093/biomet/34.1-2.111
	%Royston (1982), DOI: 10.2307/2347982
	%https://stats.stackexchange.com/questions/394960/variance-of-normal-order-statistics
	%https://stats.stackexchange.com/questions/9001/approximate-order-statistics-for-normal-random-variables
	%https://en.wikipedia.org/wiki/Extreme_value_theory
	%https://en.wikipedia.org/wiki/Gumbel_distribution
	
	%% calculate source statistics for N
	%calculate mean (E(X)) for sample size N; Elfving (1947)
	dblAlpha=pi/8;
	%fMaxE = @(N) norminv((1-dblAlpha)./((N*2)-2*dblAlpha+1));
	%dblE = -fMaxE(intN);
	fMaxE = @(N) norminv((N-dblAlpha)./(N-2*dblAlpha+1));
	dblE = fMaxE(intN);
	
	%calculate approximate variance for sample size N; Baglivo (2005)
	fMaxV = @(N) ((N) ./ (((N+1).^2).*(N+2))).*(1./((normpdf(norminv(N./(N+1)))).^2));
	dblV = fMaxV(intN);
	
	%% define Gumbel parameters from mean and variance
	%derive beta parameter from variance
	dblBeta = (sqrt(6).*sqrt(dblV))./(pi);
	
	%define Euler-Mascheroni constant
	dblEulerMascheroni = 0.5772156649015328606065120900824; %vpa(eulergamma)
	
	%derive mode from mean, beta and E-M constant
	dblMode = dblE - dblBeta.*dblEulerMascheroni;
	
	%define Gumbel cdf
	fGumbelCDF = @(x) exp(-exp(-((x(:)-dblMode)./dblBeta)));
	
	%% calculate output variables
	%calculate cum dens at X
	dblGumbelCDF = fGumbelCDF(dblX);
	%define p-value
	dblP = (1-dblGumbelCDF);
	%transform to output z-score
	dblZ = -norminv(dblP/2);
	
	% approximation for large X
	dblP(isinf(dblZ)) = exp( (dblMode-dblX(isinf(dblZ)))./dblBeta ) ;
	%transform to output z-score
	dblZ = -norminv(dblP/2);
end

