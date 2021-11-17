function [dblZetaP,sZETA] = getTraceZeta(vecTraceT,vecTraceAct,matEventTimes,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize)
	%getTraceZeta Calculates neuronal responsiveness index zeta for traces
	%syntax: [dblZetaP,sZETA] = getTraceZeta(vecTraceT,vecTraceAct,vecEventTimes,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize)
	%	input:
	%	- vecTraceT [N x 1]: time (s) corresponding to entries in vecTraceAct
	%	- vecTraceAct [N x 1]: activation trace (e.g., calcium imaging dF/F0)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all values beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 100)
	%	- intPlot: integer, plotting switch (0=none, 1=traces only, 2=activity heat map as well) (default: 0)
	%	- boolDirectQuantile; boolean, switch to use the empirical
	%							null-distribution rather than the Gumbel approximation (default: false)
	%	- dblJitterSize; scalar, sets the temporal jitter window relative to dblUseMaxDur (default: 2)
	%
	%	output:
	%	- dblZetaP; Zenith of Event-based Time-locked Anomalies: responsiveness z-score (i.e., >2 is significant)
	%	- sZETA; structure with fields:
	%		- dblZETA; FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%		- dblD; temporal deviation value underlying ZETA
	%		- dblP; p-value corresponding to ZETA
	%		- dblPeakT; time corresponding to ZETA
	%		- intPeakIdx; entry corresponding to ZETA
	%		- dblMeanD; Cohen's D based on mean-rate stim/base difference
	%		- dblMeanP; p-value based on mean-rate stim/base difference
	%		- vecTraceT: timestamps of trace entries (corresponding to vecZ)
	%		- vecD; temporal deviation vector of data
	%		- matRandD; baseline temporal deviation matrix of jittered data
	%
	%Version history:
	%1.0 - 2021 October 29
	%	Created by Jorrit Montijn
	
	%% prep data
	%ensure orientation
	vecTraceT = vecTraceT(:);
	vecTraceAct = vecTraceAct(:);
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanD = nan;
	if size(matEventTimes,2) > 2
		matEventTimes = matEventTimes';
	end
	if size(matEventTimes,2) == 2
		boolStopSupplied = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(matEventTimes(:,1)));
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 250;
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolVerbose
	if ~exist('boolDirectQuantile','var') || isempty(boolDirectQuantile)
		boolDirectQuantile = false;
	end
	
	%get dblJitterSize
	if ~exist('dblJitterSize','var') || isempty(dblJitterSize)
		dblJitterSize = 2;
	end
	
	%sampling frequency
	dblSamplingFreq = median(diff(vecTraceT));
	
	%% build onset/offset vectors
	vecEventStarts = matEventTimes(:,1);
	vecEventStops = matEventTimes(:,2);
	
	%% gettacezeta
	[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,matRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
		calcTraceZeta(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingFreq,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize);
	
	%get location
	dblMaxDTime = vecRefT(intZETALoc);
	dblD = vecRealDiff(intZETALoc);
	
	%% calculate mean-rate difference
	if boolStopSupplied
		%pre-allocate
		intMaxRep = size(vecEventStarts,1);
		vecStimAct = zeros(intMaxRep,1);
		vecBaseAct = zeros(intMaxRep,1);
		dblMedianBaseDur = median(vecEventStarts(2:end) - vecEventStops(1:(end-1)));
		intTimeNum = numel(vecTraceT);
		
		%go through trials to build spike time vector
		for intEvent=1:intMaxRep
			%% get original times
			dblStartT = vecEventStarts(intEvent,1);
			dblStopT = dblStartT+dblUseMaxDur;
			dblPreT = dblStartT - dblMedianBaseDur;
			
			intStartT = max([1 find(vecTraceT > dblStartT,1) - 1]);
			intStopT = min([intTimeNum find(vecTraceT > dblStopT,1) + 1]);
			intPreT = max([1 find(vecTraceT > dblPreT,1) - 1]);
			vecSelectFramesBase = intPreT:(intStartT-1);
			vecSelectFramesStim = intStartT:intStopT;
			
			%% get data
			vecUseBaseTrace = vecTraceAct(vecSelectFramesBase);
			vecUseStimTrace = vecTraceAct(vecSelectFramesStim);
			
			%% get activity
			vecBaseAct(intEvent) = mean(vecUseBaseTrace);
			vecStimAct(intEvent) = mean(vecUseStimTrace);
		end
		
		%get metrics
		dblMeanD = abs(nanmean(vecStimAct - vecBaseAct)) / ( (nanstd(vecStimAct) + std(vecBaseAct))/2);
		indUseTrials = ~isnan(vecStimAct) & ~isnan(vecBaseAct);
		[h,dblMeanP]=ttest(vecStimAct(indUseTrials),vecBaseAct(indUseTrials));
	end
	
	%% plot
	if intPlot
		%plot maximally 100 traces
		intPlotIters = min([size(matRandDiff,2) 100]);
		
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		if intPlot > 1
			[vecRefT,matTracePerTrial] = getTraceInTrial(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingFreq,dblUseMaxDur);
 			subplot(2,3,1)
			imagesc(vecRefT,1:size(matTracePerTrial,1),matTracePerTrial);
			colormap(hot);
			xlabel('Time from event (s)');
			ylabel('Trial #');
			title('Z-scored activation');
			fixfig;
			grid off;
		end
		
		%plot
		subplot(2,3,2)
		sOpt = struct;
		sOpt.handleFig =-1;
		sOpt.vecWindow = [0 dblUseMaxDur];
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecTraceT,vecTraceAct,vecEventStarts(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecMean,vecSEM);
		%ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean value over trials'));
		xlabel('Time from event (s)');
		ylabel('Trace value');
		fixfig
		
		subplot(2,3,3)
		plot(vecRefT,vecRealFrac)
		hold on
		plot(vecRefT,vecRealFracLinear,'color',[0.5 0.5 0.5]);
		hold off
		title(sprintf('Real data'));
		xlabel('Time from event (s)');
		ylabel('Fractional position of value in trial');
		fixfig
		
		subplot(2,3,4)
		cla;
		hold all
		for intOffset=1:intPlotIters
			plot(vecRefT,matRandDiff(:,intOffset),'Color',[0.5 0.5 0.5]);
		end
		plot(vecRefT,vecRealDiff,'Color',lines(1));
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), z(mean)=%.3f (p=%.3f)',dblZETA,dblZetaP,dblMeanD,dblMeanP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblZetaP));
		end
		fixfig
		
		
	end
	
	%% build optional output structure
	if nargin > 1
		sZETA = struct;
		sZETA.dblZETA = dblZETA;
		sZETA.dblD = dblD;
		sZETA.dblP = dblZetaP;
		sZETA.dblPeakT = dblMaxDTime;
		sZETA.intPeakIdx = intZETALoc;
		if boolStopSupplied
			sZETA.dblMeanD = dblMeanD;
			sZETA.dblMeanP = dblMeanP;
		end
		sZETA.vecTraceT = vecTraceT;
		sZETA.vecD = vecRealDiff;
		sZETA.matRandD = matRandDiff;
	end
end
