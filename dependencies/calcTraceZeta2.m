function [vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcTraceZeta2(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingInterval,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
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
	dblPreUse = -dblUseMaxDur*((dblJitterSize-1)/2);
	dblPostUse = dblUseMaxDur*((dblJitterSize+1)/2);
	
	dblStartT = min(vecEventStarts(:,1)) + dblPreUse*2;
	dblStopT = max(vecEventStarts(:,1)) + dblPostUse*2;
	indRemoveEntries = (vecTraceT < dblStartT) | (vecTraceT > dblStopT);
	vecTraceT(indRemoveEntries) = [];
	vecTraceAct(indRemoveEntries) = [];
	%calculate reference time
	vecWideT = (dblPreUse+dblSamplingInterval/2):dblSamplingInterval:dblPostUse;
	intT0 = find(vecWideT>=0,1);
	intSamples = numel(vecWideT);
	intTrials = numel(vecEventStarts(:,1));
	
	%stitch trials
	[vecRefT2,matWideTrace] = getTraceInTrial(vecTraceT,vecTraceAct,vecEventStarts+dblPreUse,dblSamplingInterval,dblPostUse-dblPreUse);
	vecPseudoStartT = vecRefT2(1):vecRefT2(end):(vecRefT2(end)*(intTrials-1)+vecRefT2(1)+eps);
	matPseudoT = bsxfun(@plus,vecPseudoStartT',vecRefT2);
	vecPseudoTrace = matWideTrace(:)';
	[vecPseudoT,vecReorder] = sort(matPseudoT(:)');
	vecPseudoTrace = vecPseudoTrace(vecReorder);
	if numel(vecPseudoT) < 3
		return;
	end
	vecPseudoTrace = vecPseudoTrace - min(vecPseudoTrace(:));
	
	%% get trial responses
	[vecRealDiff,vecRealFrac,vecRealFracLinear,vecRefT] = ...
		getTraceOffset(vecPseudoT,vecPseudoTrace,vecPseudoStartT',dblSamplingInterval,dblUseMaxDur);
	
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
			vecRandDiff = getTraceOffset(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,dblSamplingInterval,dblUseMaxDur);
			
			%assign data
			matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
			
		end
	catch
		for intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTraceOffset(vecPseudoT,vecPseudoTrace,vecStimUseOnTime,dblSamplingInterval,dblUseMaxDur);
			
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

