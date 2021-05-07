function [vecThisDiff,vecThisFrac,vecThisFracLinear] = ...
			getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur)
	%getTempOffset Calculate temporal offset vectors. Syntax:
	%[vecThisDiff,vecThisFrac,vecThisFracLinear] = ...
	%	getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur)
	%
	%This is a subfunction for getZeta().
	
	%% get temp diff vector
	%pre-allocate
	vecThisSpikeTimes = unique(getSpikeT(vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur));
	vecThisSpikeFracs = linspace(1/numel(vecThisSpikeTimes),1,numel(vecThisSpikeTimes))';
	vecThisFrac = interp1(vecThisSpikeTimes,vecThisSpikeFracs,vecSpikeT);
	
	%get linear fractions
	vecThisFracLinear = (vecSpikeT./dblUseMaxDur);
	
	%calc difference
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

