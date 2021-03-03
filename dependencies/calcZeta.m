function [vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum)
	%calcZeta Calculates neuronal responsiveness index zeta
	%[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	%	calcZeta(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intResampNum)

	%% check inputs and pre-allocate error output
	vecSpikeT = [];
	vecRealDiff = [];
	vecRealFrac = [];
	vecRealFracLinear = [];
	matRandDiff = [];
	dblZetaP = 1;
	dblZETA = 0;
	intZETALoc = nan;
	
	%% reduce spikes
	dblStartT = max([vecSpikeTimes(1) min(vecEventStarts(:,1))-dblUseMaxDur*5]);
	dblStopT = max(vecEventStarts(:,1))+dblUseMaxDur*5;
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
	
	%% run bootstraps; try parallel, otherwise run normal loop
	intSpikes = numel(vecRealDiff);
	matRandDiff = nan(intSpikes,intResampNum);
	vecStartOnly = vecEventStarts(:,1);
	try
		parfor intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + 2*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2);
			
			%get temp offset
			vecRandDiff = getTempOffset(vecSpikeT,vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur);
			
			%assign data
			matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
		end
	catch
		for intResampling=1:intResampNum
			%% get random subsample
			vecStimUseOnTime = vecStartOnly + 2*dblUseMaxDur*((rand(size(vecStartOnly))-0.5)*2);
			
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
	[dblPosD,intZETALoc]= max(abs(vecRealDiff));
	
	%calculate statistical significance using Gumbel distribution
	[dblZetaP,dblZETA] = getGumbel(dblRandMu,dblRandVar,dblPosD);
	%fprintf('Pre-correction d=%.3f,post-correction z=%.3f (p=%.3f)\n',dblD,dblZETA,dblP);
	
end

