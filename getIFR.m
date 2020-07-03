function [vecIFR,sIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,intPlot,boolVerbose)
	%getIFR Returns instaneous firing rate. Syntax:
	%   [vecIFR,sIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,intPlot,boolVerbose)
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
	%	- boolVerbose: boolean, switch to print messages
	%
	%Outputs:
	%	- vecIFR; Instantaneous firing rate
	%	- sIFR; structure with fields:
	%		- vecRate;
	%		- vecT;
	%		- vecD;
	%		- vecRate;
	%		- vecScale; 
	%
	%Version history:
	%1.0 - 24 January 2019
	%	Created by Jorrit Montijn - split from getMultiScaleDeriv.m
	%1.1 - 24 June 2020
	%	Syntax cleanup [by JM]
	%1.2 - 3 July 2020
	%	Conform to ZETA naming [by JM]
	
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
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = true;
	end
	
	%% prepare normalized spike times
	%pre-allocate
	intMaxRep = size(vecEventStarts,1);
	cellSpikeTimesPerTrial = cell(intMaxRep,1);
	
	%go through trials to build spike time vector
	for intEvent=1:intMaxRep
		%get times
		dblStartT = vecEventStarts(intEvent,1);
		dblStopT = dblStartT + dblUseMaxDur;
		
		% build trial assignment
		cellSpikeTimesPerTrial{intEvent} = vecSpikeTimes(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT) - dblStartT;
	end
	
	%get spikes in fold
	vecSpikeT = [0;sort(cell2vec(cellSpikeTimesPerTrial),'ascend');dblUseMaxDur];
	intSpikes = numel(vecSpikeT);
	
	%% get difference from uniform
	vecFracs = linspace(0,1,numel(vecSpikeT))';
	vecLinear = (vecSpikeT./max(vecSpikeT));
	vecDiff = vecFracs - vecLinear;
	vecDiff = vecDiff - mean(vecDiff);
	
	%% get multi-scale derivative
	dblMeanRate = (intSpikes/(dblUseMaxDur*intMaxRep));
	[vecIFR,sMSD] = getMultiScaleDeriv(vecSpikeT,vecDiff,intSmoothSd,dblMinScale,dblBase,intPlot,dblMeanRate,dblUseMaxDur);
	
	%% build output
	if nargout > 1
		sIFR = struct;
		sIFR.vecRate = vecIFR;
		sIFR.vecT = vecSpikeT;
		sIFR.vecD = vecDiff;
		sIFR.vecScale = sMSD.vecScale;
	end
end

