function plotRaster(vecSpikes,vecTrialStarts)
	%plotRaster Makes raster plot
	%syntax: plotRaster(vecSpikes,vecTrialStarts)
	%	input:
	%	- vecSpikes; spike times (s)
	%	- vecTrialStarts: trial start times (s)
	%
	%Version history:
	%1.0 - June 18 2019
	%	Created by Jorrit Montijn
	
	%get spike times in trials
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts);
	dblTrialDur = median(diff(vecTrialStarts));
	
	%plot per trial
	hold all;
	for intTrial=1:numel(vecTrialStarts)
		vecTimes = vecTimePerSpike(vecTrialPerSpike==intTrial);
		vecTimes(vecTimes>dblTrialDur)=[];
		line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color','k','LineWidth',1.5);
	end
	hold off
	
	%set fig props
	ylim([0 numel(vecTrialStarts)]);
	xlim([0 dblTrialDur]);
	xlabel('Time from trial start (s)');
	ylabel('Trial #');
	fixfig;
end

