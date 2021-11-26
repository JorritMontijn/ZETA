	function [vecThisDiff,vecThisFrac,vecThisFracLinear,vecRefT] = ...
			getTraceOffset3(vecTimestamps,vecData,vecEventStartT,vecRefT,dblUseMaxDur)
	%getTraceOffset Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisFrac,vecThisFracLinear] = ...
	%	getTraceOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,vecRefT,dblUseMaxDur)
	%
	%This is a subfunction for getZeta().
	
	%% prepare
	vecEventStartT = sort(vecEventStartT);
	if isempty(vecRefT)
		%pre-allocate
		intTrialNum = numel(vecEventStartT);
		intTimeNum = numel(vecTimestamps);
		%build common timeframe
		cellRefT = cell(1,intTrialNum);
		for intTrial=1:intTrialNum
			% get original times
			dblStartT = vecEventStartT(intTrial);
			dblStopT = dblStartT+dblUseMaxDur;
			intStartT = max([1 find(vecTimestamps > dblStartT,1) - 1]);
			intStopT = min([intTimeNum find(vecTimestamps > dblStopT,1)]);
			vecSelectSamples = intStartT:intStopT;
			
			%% get data
			cellRefT{intTrial} = vecTimestamps(vecSelectSamples)-dblStartT;
		end
		
		%set tol
		dblSampInterval = median(diff(vecTimestamps));
		dblTol = dblSampInterval/100;
		vecRefT = uniquetol(sort(cell2vec(cellRefT)),dblTol);
	end
	
	%get data
	%build interpolated data
	[vecRefT,matTracePerTrial] = getInterpolatedTimeSeries(vecTimestamps,vecData,vecEventStartT,dblUseMaxDur,vecRefT);
	vecMeanTrace = nanmean(matTracePerTrial,1)';
	vecThisFrac = cumsum(vecMeanTrace) / sum(vecMeanTrace);
	
	%get linear fractions
	vecThisFracLinear = linspace(mean(vecMeanTrace),sum(vecMeanTrace),numel(vecMeanTrace))' / sum(vecMeanTrace);
	
	%assign data
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

