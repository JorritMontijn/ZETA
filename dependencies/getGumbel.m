function [dblP,dblZ,dblMode,dblBeta] = getGumbel(dblE,dblV,dblX)
	%getGumbel Calculate p-value and z-score for maximum value of N samples drawn from Gaussian
	%   [dblP,dblZ] = getGumbel(dblE,dblV,dblX)
	%
	%	input:
	%	- dblE: mean of distribution of maximum values
	%	- dblV: variance of distribution of maximum values
	%	- dblX: maximum value to express in quantiles of Gumbel
	%
	%	output:
	%	- dblP; p-value for dblX (chance that sample originates from distribution given by dblE/dblV)
	%	- dblZ; z-score corresponding to P
	%
	%Version history:
	%1.0 - March 11 2020
	%	Created by Jorrit Montijn
	%
	%Sources:
	%Baglivo (2005)
	%Elfving (1947), https://doi.org/10.1093/biomet/34.1-2.111
	%Royston (1982), DOI: 10.2307/2347982
	%https://stats.stackexchange.com/questions/394960/variance-of-normal-order-statistics
	%https://stats.stackexchange.com/questions/9001/approximate-order-statistics-for-normal-random-variables
	%https://en.wikipedia.org/wiki/Extreme_value_theory
	%https://en.wikipedia.org/wiki/Gumbel_distribution
	
	%% define constants
	%define Euler-Mascheroni constant
	dblEulerMascheroni = 0.5772156649015328606065120900824; %vpa(eulergamma)
	
	%define apery's constant
	%dblApery = 1.202056903159594285399738161511449990764986292;
	
	%% define Gumbel parameters from mean and variance
	%derive beta parameter from variance
	dblBeta = (sqrt(6).*sqrt(dblV))./(pi);
	%dblSkewness = (12*sqrt(6)*dblApery)/(pi.^3);
	
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

