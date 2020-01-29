function [vecIdx,vecDur,vecEnergy] = getSustainedPeaks(vecT,vecData,vecSlope,dblCutOffZ)
	%% get durations
	intN = numel(vecData);
	vecPeakDuration = nan(intN-1,1);
	for i=1:(intN-1)
		if vecSlope(i) > 0
			intPeakEnd = find(vecData(i) > vecData((i+1):end),1);
		else
			intPeakEnd = find(vecData(i) < vecData((i+1):end),1);
		end
		if isempty(intPeakEnd)
			intPeakEnd=i;
		else
			intPeakEnd = intPeakEnd + i;
		end

		vecPeakDuration(i) = abs(vecT(intPeakEnd) - vecT(i));
	end

	%% left slope * duration
	vecLeftProminence = vecPeakDuration .* vecSlope;
	vecLPZ = zscore(vecLeftProminence);
	
	%% energy
	vecLP_Rect = vecLeftProminence;
	vecLP_Rect(abs(vecLPZ) < dblCutOffZ) = 0;
	[dummy,vecPeakLocs] = findpeaks(vecSlope);
	
	vecSignChange = [0;diff(vecLP_Rect > 0);0];
	vecStartPeaks = find(vecSignChange == 1);
	vecStopPeaks = find(vecSignChange == -1);
	%ensure results start with a start, end with a stop, and contain at least one peak
	if isempty(vecPeakLocs),vecPeakLocs=1;end
	if isempty(vecStopPeaks),vecStopPeaks=intN;end
	if isempty(vecStartPeaks),vecStartPeaks=1;end
	if vecStopPeaks(1) < vecStartPeaks(1),vecStopPeaks(1) = [];end
	if numel(vecStopPeaks) > numel(vecStartPeaks),vecStartPeaks = [1;vecStartPeaks];end
	if numel(vecStartPeaks) > numel(vecStopPeaks),vecStopPeaks = [vecStopPeaks;intN];end
	%analyze peaks
	intPeaks = numel(vecStartPeaks);
	vecLocs = nan(intPeaks,1);
	vecEnergy = nan(intPeaks,1);
	vecDur = nan(intPeaks,1);
	for intPeakIdx = 1:intPeaks
		vecIncluded = vecStartPeaks(intPeakIdx):vecStopPeaks(intPeakIdx);
		intUsePeak = find(ismember(vecPeakLocs,vecIncluded),1);
		if isempty(intUsePeak)
			[dummy,intUsePeak] = min(abs(vecIncluded(1)-vecPeakLocs));
		end
		vecLocs(intPeakIdx) = vecPeakLocs(intUsePeak);
		vecEnergy(intPeakIdx) = sum(vecData(vecIncluded));
		vecDur(intPeakIdx) = abs(vecT(vecStopPeaks(intPeakIdx)) - vecT(vecStartPeaks(intPeakIdx)));
	end
	% sort peaks
	[vecLocs,vecKeep] = unique(vecLocs);
	vecEnergy = vecEnergy(vecKeep);
	vecDur = vecDur(vecKeep);
	
	[vecEnergy,vecReorder] = sort(vecEnergy,'descend');
	vecIdx = vecLocs(vecReorder);
	vecDur = vecDur(vecReorder);
end