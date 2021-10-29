function [vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcTraceZeta(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingFreq,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	%calcTraceZeta Calculates neuronal responsiveness index zeta
	%[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcTraceZeta(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingFreq,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	
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
	dblStartT = max([vecTraceT(1) min(vecEventStarts(:,1))-dblUseMaxDur*5*dblJitterSize]);
	dblStopT = max(vecEventStarts(:,1))+dblUseMaxDur*5*dblJitterSize;
	indRemoveEntries = vecTraceT < dblStartT | vecTraceT > dblStopT;
	vecTraceT(indRemoveEntries) = [];
	if numel(vecTraceT) < 3
		return;
	end
	
	%% get trial responses
	[vecRealDiff,vecRealFrac,vecRealFracLinear,vecRefT] = ...
		getTraceOffset(vecTraceT,vecTraceAct,vecEventStarts(:,1),dblSamplingFreq,dblUseMaxDur);
	
	%% run bootstraps; try parallel, otherwise run normal loop
	intSpikes = numel(vecRealDiff);
	matRandDiff = nan(intSpikes,intResampNum);
	vecStartOnly = vecEventStarts(:,1);
	try
		parfor intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2);
			
			%get temp offset
			vecRandDiff = getTraceOffset(vecTraceT,vecTraceAct,vecStimUseOnTime,dblSamplingFreq,dblUseMaxDur);
			
			%assign data
			matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
		end
	catch
		for intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2);
			
			%get temp offset
			vecRandDiff = getTraceOffset(vecTraceT,vecTraceAct,vecStimUseOnTime,dblSamplingFreq,dblUseMaxDur);
			
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

