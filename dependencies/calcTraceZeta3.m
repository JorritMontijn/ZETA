function [vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcTraceZeta3(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingInterval,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	%calcTraceZeta Calculates neuronal responsiveness index zeta
	%[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcTraceZeta(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingInterval,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	
	%% check inputs and pre-allocate error output
	vecRefT = [];
	vecRealDiff = [];
	vecRealFrac = [];
	vecRealFracLinear = [];
	matRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	
	%% reduce data
	if size(vecEventStarts,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	%discard leading/lagging data
	vecEventStarts = vecEventStarts(:,1);
	dblPreUse = -dblUseMaxDur*((dblJitterSize-1)/2);
	dblPostUse = dblUseMaxDur*((dblJitterSize+1)/2);
	
	dblStartT = min(vecEventStarts) + dblPreUse*2 - dblJitterSize*dblUseMaxDur;
	dblStopT = max(vecEventStarts) + dblPostUse*2 + dblJitterSize*dblUseMaxDur;
	indRemoveEntries = (vecTraceT < dblStartT) | (vecTraceT > dblStopT);
	vecTraceT(indRemoveEntries) = [];
	vecTraceAct(indRemoveEntries) = [];
	%calculate reference time
	intTrials = numel(vecEventStarts);
	
	%stitch trials
	vecTimestamps = vecTraceT;
	vecData = vecTraceAct;
	vecEventT = vecEventStarts+dblPreUse;
	dblWindowDur = dblPostUse-dblPreUse;
	[vecPseudoTime,vecPseudoData,vecPseudoStartT] = getPseudoTimeSeries(vecTimestamps,vecData,vecEventT,dblWindowDur);
	vecPseudoData = vecPseudoData - min(vecPseudoData);
			
	%old
	if numel(vecPseudoTime) < 3
		return;
	end
	
	%% get trial responses
	[vecRealDiff,vecRealFrac,vecRealFracLinear,vecRefT] = ...
		getTraceOffset3(vecPseudoTime,vecPseudoData,vecPseudoStartT',[],dblUseMaxDur);
	intSamples = numel(vecRealDiff);
	
	%% run bootstraps; try parallel, otherwise run normal loop
	matRandDiff = nan(intSamples,intResampNum);
	vecStartOnly = vecPseudoStartT(:);
	vecJitterPerTrial = dblJitterSize*linspace(dblUseMaxDur/intTrials,dblUseMaxDur,intTrials)';
	matJitterPerTrial = nan(intTrials,intResampNum);
	for intResampling=1:intResampNum
		matJitterPerTrial(:,intResampling) = vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
	end
	try
		parfor intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTraceOffset3(vecPseudoTime,vecPseudoData,vecStimUseOnTime,vecRefT,dblUseMaxDur);
			
			%assign data
			matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
			
		end
	catch
		for intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTraceOffset3(vecPseudoTime,vecPseudoData,vecStimUseOnTime,vecRefT,dblUseMaxDur);
			if any(isnan(vecRandDiff))
				pause
			end
			%assign data
			matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
		end
	end
	
	%% calculate measure of effect size (for equal n, d' equals Cohen's d)
	%find highest peak and retrieve value
	vecMaxRandD = max(abs(matRandDiff),[],1);
	dblRandMu = mean(vecMaxRandD);
	dblRandVar = var(vecMaxRandD);
	[dblMaxD,intZETALoc]= max(abs(vecRealDiff));
	
	if boolDirectQuantile
		%calculate statistical significance using empirical quantiles
		%define p-value
		dblZetaP = 1 - (sum(dblMaxD>vecMaxRandD)/(1+numel(vecMaxRandD)));
		
		%transform to output z-score
		dblZETA = -norminv(dblZetaP/2);
	else
		%calculate statistical significance using Gumbel distribution
		[dblZetaP,dblZETA] = getGumbel(dblRandMu,dblRandVar,dblMaxD);
		%fprintf('Pre-correction d=%.3f,post-correction z=%.3f (p=%.3f)\n',dblD,dblZETA,dblP);
	end
	
end

