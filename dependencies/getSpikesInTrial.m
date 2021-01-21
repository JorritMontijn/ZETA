function [vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts,dblMaxDur)
	%getSpikesInTrial Retrieves spiking times per trial
	%syntax: [vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts,dblMaxDur)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%	- dblTrialDur: (optional) if supplied, uses trial duration for
	%			computations instead of assigning spikes directly to trials 
	%
	%Version history:
	%1.0 - June 26 2019
	%	Created by Jorrit Montijn
	%2.0 - February 7 2020
	%	Added overlap by using dblMaxDur [by JM]
	
	%check inputs
	if nargin < 3
		dblMaxDur = [];
	end
	
	if ~isempty(dblMaxDur)
		%sort spikes
		intTrials = numel(vecTrialStarts);
		cellTrialPerSpike = cell(size(vecTrialStarts));
		cellTimePerSpike = cell(size(vecTrialStarts));
		for intTrial=1:intTrials
			%get spikes
			vecTheseSpikes = vecSpikes((vecSpikes >= vecTrialStarts(intTrial)) & vecSpikes < (vecTrialStarts(intTrial)+ dblMaxDur));
			vecTheseSpikes = vecTheseSpikes - vecTrialStarts(intTrial);
			
			% assign
			cellTrialPerSpike{intTrial} = intTrial*ones(size(vecTheseSpikes));
			cellTimePerSpike{intTrial} = vecTheseSpikes;
		end
		%transform to vectors & reorder
		vecTimePerSpike = cell2vec(cellTimePerSpike);
		vecTrialPerSpike = cell2vec(cellTrialPerSpike);
	else
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
end

