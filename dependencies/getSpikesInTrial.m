function [vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts)
	%getSpikesInTrial Retrieves spiking times per trial
	%syntax: [vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%
	%Version history:
	%1.0 - June 26 2019
	%	Created by Jorrit Montijn
	
	
	%sort spikes
	vecTrialPerSpike = nan(size(vecSpikes));
	vecTimePerSpike = nan(size(vecSpikes));
	for intSpike=1:numel(vecSpikes)
		%% build trial assignment
		vecTrialPerSpike(intSpike) = sum(vecTrialStarts < vecSpikes(intSpike));
		if vecTrialPerSpike(intSpike) > 0
			dblRemTime = vecTrialStarts(vecTrialPerSpike(intSpike));
		else
			dblRemTime = 0;
		end
		vecTimePerSpike(intSpike) = vecSpikes(intSpike) - dblRemTime;
	end
	
end

