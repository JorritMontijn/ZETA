function [vecIFR,sIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,intPlot)
	%getIFR Returns instaneous firing rate. Syntax:
	%   [vecIFR,sIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,intPlot)
	%Required input:
	%	- vecSpikeTimes [S x 1]: spike times (s)
	%	- vecEventStarts [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all spikes beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%
	%Optional inputs:
	%	- intSmoothSd: Gaussian SD of smoothing kernel (in # of bins) [default: 3]
	%	- dblMinScale: minimum derivative scale in seconds [default: 1/1000]
	%	- dblBase: critical value for locally dynamic derivative [default: 4]
	%	- intPlot: integer, plotting switch (0=none, 1=plot)
	%
	%Outputs:
	%	- vecIFR; Instantaneous firing rate
	%	- sIFR; structure with fields:
	%		- vecRate;
	%		- vecTime;
	%		- vecDiff;
	%		- vecScale; 
	%
	%Version history:
	%1.0 - 24 January 2019
	%	Created by Jorrit Montijn - split from getMultiScaleDeriv.m
	%1.1 - 24 June 2020
	%	Syntax cleanup [by JM]
	%1.2 - 3 July 2020
	%	Conform to ZETA naming [by JM]
	%1.3 - 14 October 2021
	%	Fixed bug [by JM]
	
	%% set default values
	if ~exist('intSmoothSd','var') || isempty(intSmoothSd)
		intSmoothSd = 5;
	end
	if ~exist('dblBase','var') || isempty(dblBase)
		dblBase = 1.5;
	end
	if ~exist('dblMinScale','var') || isempty(dblMinScale)
		dblMinScale = round(log(1/1000) / log(dblBase));
	end
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(vecEventStarts(:,1)));
	end
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	if size(vecEventStarts,2) > 2
		vecEventStarts = vecEventStarts';
	end
	
	%% prepare normalized spike times
	vecSpikeT = getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur);
	intSpikes = numel(vecSpikeT);
	
	%% get difference from uniform
	vecRealDiff = ...
		getTempOffset(vecSpikeT,vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);
	if numel(vecRealDiff) < 3
		return
	end
	
	%% get multi-scale derivative
	intMaxRep = size(vecEventStarts,1);
	dblMeanRate = (intSpikes/(dblUseMaxDur*intMaxRep));
	[vecIFR,sMSD] = getMultiScaleDeriv(vecSpikeT,vecRealDiff,intSmoothSd,dblMinScale,dblBase,intPlot,dblMeanRate,dblUseMaxDur);
	
	%% build output
	if nargout > 1
		sIFR = struct;
		sIFR.vecRate = vecIFR;
		sIFR.vecTime = vecSpikeT;
		sIFR.vecDiff = vecRealDiff;
		sIFR.vecScale = sMSD.vecScale;
	end
end

