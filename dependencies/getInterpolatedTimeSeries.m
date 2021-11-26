function [vecRefT,matTracePerTrial] = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur,vecRefT)
	%getTraceInTrial Builds common timeframe
	%syntax: [vecRefT,matTracePerTrial] = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur,vecRefT)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%
	%Version history:
	%1.0 - June 26 2019
	%	Created by Jorrit Montijn
	
	%% assign data
	matTracePerTrial = nan(numel(vecEventStartT),numel(vecRefT));
	for intTrial=1:numel(vecEventStartT)
		%% get original times
		dblStartT = vecEventStartT(intTrial);
		intStartT = max([1 find(vecTimestamps > (dblStartT + vecRefT(1)),1) - 1]);
		intStopT = min([numel(vecTimestamps) find(vecTimestamps > (dblStartT + vecRefT(end)),1) + 1]);
		vecSelectSamples = intStartT:intStopT;
		
		%% get data
		vecUseTimes = vecTimestamps(vecSelectSamples);
		vecUseTrace = vecData(vecSelectSamples);
		
		%% interpolate
		vecUseInterpT = vecRefT+dblStartT;
		
		%get real fractions for training set
		matTracePerTrial(intTrial,:) = interp1(vecUseTimes,vecUseTrace,vecUseInterpT);
	end
end