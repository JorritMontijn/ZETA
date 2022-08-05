function vecSpikeT = getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur)
	%getSpikeT Summary of this function goes here
	%   vecSpikeT = getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur)
	
	%pre-allocate
	vecSpikeT = nan(numel(vecSpikeTimes),1);
	intIdx = 1;
	intMaxRep = size(vecEventStarts,1);
	
	%go through trials to build spike time vector
	for intEvent=1:intMaxRep
		%get times
		dblStartT = vecEventStarts(intEvent,1);
		dblStopT = dblStartT + dblUseMaxDur;
		
		% build trial assignment
		vecTempSpikes = vecSpikeTimes(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT) - dblStartT;
		intTempSpikeNr = numel(vecTempSpikes);
		vecSpikeT(intIdx:(intIdx+intTempSpikeNr-1)) = vecTempSpikes;
		intIdx = intIdx + intTempSpikeNr;
	end
	vecSpikeT(intIdx:end) = [];
	
	%get spikes in fold
	vecSpikeT = [0;sort(vecSpikeT(:),'ascend');dblUseMaxDur];
end

