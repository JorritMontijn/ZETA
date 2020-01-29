function [vecThisDiff,vecThisFrac,vecThisFracLinear] = ...
			getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur)
	%getTempOffset Calculate temporal offset vectors across folds and offsets. Syntax:
	%[vecThisDiff,vecThisFrac,vecThisFracLinear] = ...
	%	getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur)
	%
	%This is a subfunction for getZeta().
	
	%% get inputs
	intMaxRep = numel(vecStimUseOnTime);
	
	%% get temp diff vector
	%pre-allocate
	cellSpikeTimesPerTrial = cell(intMaxRep,1);
	
	%go through trials to build spike time vector
	for intEvent=1:intMaxRep
		%get times
		dblStartT = vecStimUseOnTime(intEvent);
		dblStopT = dblStartT + dblUseMaxDur;
		
		% build trial assignment
		cellSpikeTimesPerTrial{intEvent} = vecSpikeTimes(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT) - dblStartT;
	end
	
	%get spikes in fold
	vecThisSpikeT = unique(cell2vec(cellSpikeTimesPerTrial));
	
	%get real fractions for training set
	vecThisSpikeTimes = sort([0;vecThisSpikeT(:);dblUseMaxDur],'ascend');
	vecThisSpikeFracs = linspace(0,1,numel(vecThisSpikeTimes))';
	vecThisFrac = interp1(vecThisSpikeTimes,vecThisSpikeFracs,vecSpikeT);
	
	%get linear fractions
	vecThisFracLinear = (vecSpikeT./dblUseMaxDur);
	
	%calc difference
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	vecThisDiff = vecThisDiff - mean(vecThisDiff);
end

