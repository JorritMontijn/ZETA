function [vecPseudoTime,vecPseudoData,vecPseudoStartT] = getPseudoTimeSeries(vecTimestamps,vecData,vecEventT,dblWindowDur)
	
	%ensure sorting and alignment
	[vecTimestamps,vecReorder] = sort(vecTimestamps(:));
	vecData = vecData(:);
	vecData = vecData(vecReorder);
	vecEventT = sort(vecEventT(:));
	
	%% pre-allocate
	dblMedianDur = median(diff(vecTimestamps));
	intSamples = numel(vecTimestamps);
	intTrials = numel(vecEventT);
	cellPseudoTime = cell(1,intTrials);
	cellPseudoData = cell(1,intTrials);
	vecPseudoStartT = nan(1,intTrials);
	dblStartNextAtT = 0;
	intLastUsedSample = 0;
	%run
	for intTrial=1:intTrials
		dblEventT = vecEventT(intTrial);
		intStartSample = (find(vecTimestamps > dblEventT,1)-1);
		intEndSample = (find(vecTimestamps > (dblEventT+dblWindowDur),1));
		if isempty(intEndSample)
			intEndSample = intStartSample;
		end
		vecEligibleSamples = intStartSample:intEndSample;
		indRemSamples = (vecEligibleSamples <= 0) | (vecEligibleSamples > intSamples);
		vecUseSamples = vecEligibleSamples(~indRemSamples);
		
		%check if beginning or end
		if intTrial==1
			vecUseSamples = 1:vecUseSamples(end);
		end
		if intTrial==intTrials
			vecUseSamples = vecUseSamples(1):intSamples;
		end
		
		vecUseT = vecTimestamps(vecUseSamples);
		indOverlap = (vecUseSamples <= intLastUsedSample);
		if any(indOverlap)
			vecUseSamples = vecUseSamples(~indOverlap);
			vecUseT = vecTimestamps(vecUseSamples);
		end
		
		
		if isempty(vecUseSamples)
			vecLocalPseudoT = [];
			vecLocalPseudoV = [];
			dblPseudoEventT = dblEventT - vecTimestamps(intLastUsedSample) + dblStartNextAtT;
		else
			intLastUsedSample = vecUseSamples(end);
			vecLocalPseudoV = vecData(vecUseSamples);
			vecLocalPseudoT = vecUseT - vecUseT(1) + dblStartNextAtT;
			dblPseudoEventT = dblEventT - vecUseT(1) + dblStartNextAtT;
			
			if numel(vecTimestamps) > intLastUsedSample
				dblStepEnd = vecTimestamps(intLastUsedSample+1) - vecTimestamps(intLastUsedSample);
			else
				dblStepEnd = dblMedianDur;
			end
			dblStartNextAtT = vecLocalPseudoT(end) + dblStepEnd;
		end
		if intTrial==1
			intFirstSample = vecUseSamples(1);
			dblPseudoT0 = vecLocalPseudoT(1);
		end
		
		
		cellPseudoTime{intTrial} = vecLocalPseudoT;
		cellPseudoData{intTrial} = vecLocalPseudoV;
		vecPseudoStartT(intTrial) = dblPseudoEventT;
		
	end
	
	%add beginning
	dblT1 = vecTimestamps(intFirstSample);
	intT0 = find(vecTimestamps > (dblT1 - dblWindowDur),1);
	if ~isempty(intT0) && intFirstSample > 1
		dblStepBegin = vecTimestamps(intFirstSample) - vecTimestamps(intFirstSample-1);
		vecSampAddBeginning = intT0:(intFirstSample-1);
		cellPseudoTime = cat(2,{vecTimestamps(vecSampAddBeginning) - vecTimestamps(vecSampAddBeginning(1)) + dblPseudoT0 - dblStepBegin - range(vecTimestamps(vecSampAddBeginning))},cellPseudoTime);
		cellPseudoData = cat(2,{vecData(vecSampAddBeginning)},cellPseudoData);
	end
	
	%add end
	dblTn = vecTimestamps(intLastUsedSample);
	intTn1 = find(vecTimestamps > (dblTn + dblWindowDur),1);
	if ~isempty(intTn1)
		vecSampAddEnd = (intLastUsedSample+1):intTn1;
		cellPseudoTime = cat(2,cellPseudoTime,{vecTimestamps(vecSampAddEnd) - vecTimestamps(vecSampAddEnd(1)) + dblStartNextAtT});
		cellPseudoData = cat(2,cellPseudoData,{vecData(vecSampAddEnd)});
	end
	
	%recombine into vector
	vecPseudoTime = cell2vec(cellPseudoTime);
	vecPseudoData = cell2vec(cellPseudoData);
end

