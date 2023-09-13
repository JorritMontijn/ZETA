function [vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize)
	%calcZeta Calculates neuronal responsiveness index zeta
	%[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum,boolDirectQuantile)

	%% check inputs and pre-allocate error output
	vecSpikeT = [];
	vecRealDiff = [];
	vecRealFrac = [];
	vecRealFracLinear = [];
	matRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end
	
	%% reduce spikes
	if size(vecEventStarts,2)>2,error([mfilename ':IncorrectMatrixForm'],'Incorrect input form for vecEventStarts; size must be [m x 1] or [m x 2]');end
	dblStartT = max([vecSpikeTimes(1) min(vecEventStarts(:,1))-dblUseMaxDur*5*dblJitterSize]);
	dblStopT = max(vecEventStarts(:,1))+dblUseMaxDur*5*dblJitterSize;
	vecSpikeTimes(vecSpikeTimes < dblStartT | vecSpikeTimes > dblStopT) = [];
	if numel(vecSpikeTimes) < 3
		return;
	end
	
	%% prepare interpolation points
	vecSpikeT = getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur);
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFrac,vecRealFracLinear] = ...
		getTempOffset(vecSpikeT,vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur);
	if numel(vecRealDiff) < 3
		return
	end
	vecRealDiff = vecRealDiff - mean(vecRealDiff);
	
	%% run bootstraps; try parallel, otherwise run normal loop
	intSpikes = numel(vecRealDiff);
	matRandDiff = nan(intSpikes,intResampNum);
	vecStartOnly = vecEventStarts(:,1);
	intTrials = numel(vecStartOnly);
	vecJitterPerTrial = dblJitterSize*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2); %original zeta
	%vecJitterPerTrial = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials)'; %new
	matJitterPerTrial = nan(intTrials,intResampNum);
	for intResampling=1:intResampNum
		matJitterPerTrial(:,intResampling) = vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
	end
	try
		parfor intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
			
			%assign data
			matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
		end
	catch
		for intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + matJitterPerTrial(:,intResampling);
			
			%get temp offset
			vecRandDiff = getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
			
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

