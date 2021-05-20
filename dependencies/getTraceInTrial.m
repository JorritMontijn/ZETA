function [vecRefT,matTracePerTrial] = getTraceInTrial(vecTimestamps,vecTrace,vecEventStarts,dblSamplingInterval,dblUseMaxDur)
	%getTraceInTrial Builds common timeframe
	%syntax: [vecRefT,matTracePerTrial] = getTraceInTrial(vecTimestamps,vecTrace,vecEventStarts,dblSamplingInterval,dblUseMaxDur)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%
	%Version history:
	%1.0 - June 26 2019
	%	Created by Jorrit Montijn
	
	%% prepare
	%build common timeframe
	vecRefT = (dblSamplingInterval/2):dblSamplingInterval:dblUseMaxDur;
	
	%pre-allocate
	intTrialNum = numel(vecEventStarts);
	intTimeNum = numel(vecTimestamps);
	matTracePerTrial = nan(intTrialNum,numel(vecRefT));
	
	%% assign data
	for intTrial=1:intTrialNum
		%% get original times
		dblStartT = vecEventStarts(intTrial);
		dblStopT = dblStartT+dblUseMaxDur;
		intStartT = max([1 find(vecTimestamps > dblStartT,1) - 1]);
		intStopT = min([intTimeNum find(vecTimestamps > dblStopT,1) + 1]);
		vecSelectFrames = intStartT:intStopT;
		
		%% get data
		vecUseTimes = vecTimestamps(vecSelectFrames);
		vecUseTrace = vecTrace(vecSelectFrames);
		
		%% interpolate
		vecUseInterpT = vecRefT+dblStartT;
		
		%get real fractions for training set
		matTracePerTrial(intTrial,:) = interp1(vecUseTimes,vecUseTrace,vecUseInterpT);
	end
	
end

