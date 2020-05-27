function [dblZetaP,vecLatencies,sZETA,sRate] = getTraceZeta(vecTraceT,vecTraceAct,varEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolVerbose)
	%getTraceZeta Calculates neuronal responsiveness index zeta for traces
	%syntax: [dblZetaP,vecLatencies,sZETA,sRate] = getTraceZeta(vecTraceT,vecTraceAct,varEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolVerbose)
	%	input:
	%	- vecTraceT [N x 1]: time (in seconds) corresponding to entries in vecTraceAct
	%	- vecTraceAct [N x 1]: activation trace (e.g., calcium imaging dF/F0)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times to calculate difference of on/off averages
	%	- dblUseMaxDur: float (s), window length for calculating ZETA: ignore all values beyond this duration after event onset
	%								[default: median of event onset to event onset]
	%	- intResampNum: integer, number of resamplings (default: 50)
	%	- intPlot: integer, plotting switch (0=none, 1=traces only, 2=activity heat map as well) (default: 0)
	%	- intLatencyPeaks: integer, maximum number of latency peaks to return (default: 3)
	%	- boolVerbose: boolean, switch to print messages
	%
	%	output:
	%	- dblZetaP; p-value based on Zenith of Event-based Time-locked Anomalies
	%	- vecLatencies; different latency estimates, number determined by intLatencyPeaks. If no peaks are detected, it returns NaNs:
	%		1) Latency of ZETA
	%		2) Latency of largest z-score with inverse sign to ZETA
	%		3) Peak time of instantaneous firing rate
	%		4) Onset time of above peak, defined as the first crossing of peak half-height
	%	- sZETA; structure with fields:
	%		- dblZETA; FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%		- dblD; temporal deviation value underlying ZETA
	%		- dblP; p-value corresponding to ZETA
	%		- dblPeakT; time corresponding to ZETA
	%		- intPeakIdx; entry corresponding to ZETA
	%		- dblMeanD; Cohen's D based on averaged stim/base difference
	%		- dblMeanP; p-value based on averaged stim/base difference
	%		- vecRefT: timestamps of data points (corresponding to vecD)
	%		- vecD; temporal deviation vector of data
	%		- matRandD; baseline temporal deviation matrix of jittered data
	%		- dblD_InvSign; largest peak of inverse sign to ZETA (i.e., -ZETA)
	%		- dblPeakT_InvSign; time corresponding to -ZETA
	%		- intPeakIdx_InvSign; entry corresponding to -ZETA
	%		- dblUseMaxDur; window length used to calculate ZETA
	%	- sRate; structure with fields: (only if intLatencyPeaks > 0)
	%		- vecRate; instantaneous spiking rates (like a PSTH)
	%		- vecT; time-points corresponding to vecRate (same as sZETA.vecRefT)
	%		- vecM; Mean of multi-scale derivatives
	%		- vecScale; timescales used to calculate derivatives
	%		- matMSD; multi-scale derivatives matrix
	%		- vecV; values on which vecRate is calculated (same as sZETA.vecD)
	%		Data on the peak:
	%		- dblPeakTime; time of peak (in seconds)
	%		- dblPeakWidth; duration of peak (in seconds)
	%		- vecPeakStartStop; start and stop time of peak (in seconds)
	%		- intPeakLoc; spike index of peak (corresponding to sZETA.vecSpikeT)
	%		- vecPeakStartStopIdx; spike indices of peak start/stop (corresponding to sZETA.vecSpikeT)
	%		Additionally, it will return peak onset latency (first crossing of peak half-height) using getOnset.m:
	%		- dblOnset: latency for peak onset
	%
	%Version history:
	%0.9 - 27 June 2019
	%	Created by Jorrit Montijn
	%1.0 - 24 September 2019
	%	New procedure to determine statistical significance [by JM]
	%2.0 - 27 January 2020
	%	New peak detection procedure using multi-scale derivatives [by JM]
	%2.1 - 5 February 2020
	%	Minor changes and bug fixes [by JM]
	%2.2 - 11 February 2020
	%	Peak width, analytical ZETA correction [by JM]
	%2.3 - 13 March 2020
	%	MSD-based instantaneous spiking rates, onset latency [by JM]
	%2.4 - 13 March 2020
	%	Closed-form statistical significance using Gumbel distribution [by	JM]
	%2.5 - 27 May 2020
	%	Standardized syntax and variable names [by JM]

	
	%% prep data
	%ensure orientation
	vecTraceT = vecTraceT(:);
	vecTraceAct = vecTraceAct(:);
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblMeanD = nan;
	if size(varEventTimes,2) > 2
		varEventTimes = varEventTimes';
	end
	if size(varEventTimes,2) == 2
		boolStopSupplied = true;
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = median(diff(varEventTimes(:,1)));
	end
	
	%get resampling num
	if ~exist('intResampNum','var') || isempty(intResampNum)
		intResampNum = 50;
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get boolPlot
	if ~exist('intLatencyPeaks','var') || isempty(intLatencyPeaks)
		intLatencyPeaks = 4;
	end
	
	%get boolPlot
	if ~exist('vecRestrictRange','var') || isempty(vecRestrictRange)
		vecRestrictRange = [-inf inf];
	end

	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	
	%sampling frequency
	dblSamplingFreq = median(diff(vecTraceT));
	
	%% build onset/offset vectors
	vecEventStarts = varEventTimes(:,1);
	
	%% get trial responses
	intMaxRep = size(vecEventStarts,1);
	[vecRefT,matTracePerTrial] = getTraceInTrial(vecTraceT,vecTraceAct,vecEventStarts,dblSamplingFreq,dblUseMaxDur);
	
	%% run normal
	%get data
	vecMeanTrace = nanmean(matTracePerTrial,1)';
	vecRealFrac = cumsum(vecMeanTrace) / sum(vecMeanTrace);
	
	%get linear fractions
	vecRealFracLinear = linspace(median(vecMeanTrace),sum(vecMeanTrace),numel(vecMeanTrace))' / sum(vecMeanTrace);
	
	%assign data
	vecRealDiff = vecRealFrac - vecRealFracLinear;
	vecRealDiff = vecRealDiff - mean(vecRealDiff);
	
	%% run bootstraps
	matRandDiff = nan(numel(vecMeanTrace),intResampNum);
	hTic = tic;
	for intResampling=1:intResampNum
		%% msg
		if boolVerbose && toc(hTic) > 5
			fprintf('Now at resampling %d/%d\n',intResampling,intResampNum);
			hTic = tic;
		end
		
		%% get random subsample
		vecStimUseOnTime = vecEventStarts(:,1) + 2*median(diff(vecEventStarts(:,1)))*rand(size(vecEventStarts(:,1)));
		
		%get trial responses
		[vecRefT,matRandTracePerTrial] = getTraceInTrial(vecTraceT,vecTraceAct,vecStimUseOnTime,dblSamplingFreq,dblUseMaxDur);
		
		%get data
		vecMeanRandTrace = nanmean(matRandTracePerTrial,1)';
		vecRandFrac = cumsum(vecMeanRandTrace) / sum(vecMeanRandTrace);
		
		%get linear fractions
		vecRandFracLinear = linspace(median(vecMeanRandTrace),sum(vecMeanRandTrace),numel(vecMeanRandTrace))' / sum(vecMeanRandTrace);
		
		%assign data
		vecRandDiff = vecRandFrac - vecRandFracLinear;
		vecRandDiff = vecRandDiff - mean(vecRandDiff);
		
		%assign data
		matRandDiff(:,intResampling) = vecRandDiff;
	end
	
	
	%% calculate measure of effect size (for equal n, d' equals Cohen's d)
	if numel(vecRealDiff) < 3
		dblZETA = 0;
		sZETA = [];
		vecLatencies = [];
		sRate = [];
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		
		%build placeholder outputs
		if numel(vecLatencies) < intLatencyPeaks
			vecLatencies(end+1:intLatencyPeaks) = nan;
		end
		sZETA = struct;
		sZETA.dblZETA = 0;
		sZETA.dblD = 0;
		sZETA.dblP = 1;
		sZETA.dblPeakT = nan;
		sZETA.intPeakIdx = [];
		if boolStopSupplied
			sZETA.dblMeanD = 0;
			sZETA.dblMeanP = 1;
		end
		sZETA.vecSpikeT = [];
		sZETA.vecD = [];
		sZETA.matRandD = [];
		
		sZETA.dblD_InvSign = 0;
		sZETA.dblPeakT_InvSign = nan;
		sZETA.intPeakIdx_InvSign = [];
		sZETA.dblUseMaxDur = nan;
		sZETA.vecLatencyVals = [];
		return
	end
	
	%find highest peak and retrieve value
	vecMaxRandD = max(abs(matRandDiff),[],1);
	dblRandMu = mean(vecMaxRandD);
	dblRandVar = var(vecMaxRandD);
	[dblPosD,intZETALoc]= max(abs(vecRealDiff));
	
	%get location
	dblMaxDTime = vecTraceT(intZETALoc);
	dblD = vecRealDiff(intZETALoc);
	
	%calculate statistical significance using Gumbel distribution
	[dblZetaP,dblZETA] = getGumbel(dblRandMu,dblRandVar,dblPosD);
	%fprintf('Pre-corr d=%.3f,post-corr z=%.3f (p=%.3f)\n',dblD,dblZETA,dblP);
	
	%find peak of inverse sign
	[dummy,intPeakLocInvSign] = max(-sign(dblD)*vecRealDiff);
	dblMaxDTimeInvSign = vecTraceT(intPeakLocInvSign);
	dblD_InvSign = vecRealDiff(intPeakLocInvSign);

	%build common timeframe
	vecRefT = (dblSamplingFreq/2):dblSamplingFreq:dblUseMaxDur;
	
	if boolStopSupplied
		%% calculate mean-rate difference
		%pre-allocate
		vecEventStops = varEventTimes(:,2);
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
		%plot maximally 50 traces
		intPlotIters = min([size(matRandDiff,2) 50]);
		
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		if intPlot > 2
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
		scatter(dblMaxDTime,vecRealDiff(intZETALoc),'bx');
		scatter(dblMaxDTimeInvSign,vecRealDiff(intPeakLocInvSign),'b*');
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), d(Hz)=%.3f (p=%.3f)',dblZETA,dblZetaP,dblMeanD,dblMeanP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblZetaP));
		end
		fixfig
		
		
	end
	
	%% calculate MSD
	if intLatencyPeaks > 0
		dblMeanRate = mean(vecTraceAct);
		[vecRate,sRate] = getMultiScaleDeriv(vecRefT,vecRealDiff,0,[],[],intPlot,dblMeanRate,dblUseMaxDur);
	else
		sRate = [];
	end
	
	%% calculate MSD statistics
	if ~isempty(sRate) && intLatencyPeaks > 0
		%get MSD peak
		[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = getPeak(vecRate,vecRefT,vecRestrictRange);
		
		sRate.dblPeakRate = dblPeakRate;
		sRate.dblPeakTime = dblPeakTime;
		sRate.dblPeakWidth = dblPeakWidth;
		sRate.vecPeakStartStop = vecPeakStartStop;
		sRate.intPeakLoc = intPeakLoc;
		sRate.vecPeakStartStopIdx = vecPeakStartStopIdx;
		
		
		if ~isnan(dblPeakTime)
			%assign array data
			if intLatencyPeaks > 3
				%get onset
				[dblOnset,dblOnsetVal] = getOnset(vecRate,vecRefT,dblPeakTime,vecRestrictRange);
				sRate.dblOnset = dblOnset;
				vecLatencies = [dblMaxDTime dblMaxDTimeInvSign dblPeakTime dblOnset];
				vecLatencyVals = [vecRate(intZETALoc) vecRate(intPeakLocInvSign) vecRate(intPeakLoc) dblOnsetVal];
			else
				sRate.dblOnset = [nan];
				vecLatencies = [dblMaxDTime dblMaxDTimeInvSign dblPeakTime];
				vecLatencyVals = [vecRate(intZETALoc) vecRate(intPeakLocInvSign) vecRate(intPeakLoc)];
			end
			vecLatencies = vecLatencies(1:intLatencyPeaks);
			vecLatencyVals = vecLatencyVals(1:intLatencyPeaks);
			if intPlot > 0
				hold on
				scatter(dblPeakTime,vecRate(intPeakLoc),'gx');
				scatter(dblMaxDTime,vecRate(intZETALoc),'bx');
				scatter(dblMaxDTimeInvSign,vecRate(intPeakLocInvSign),'b*');
				if intLatencyPeaks > 3
					scatter(dblOnset,dblOnsetVal,'rx');
					title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms,On=%.2fms',dblMaxDTime*1000,dblMaxDTimeInvSign*1000,dblPeakTime*1000,dblOnset*1000));
				else
					title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms',dblMaxDTime*1000,dblMaxDTimeInvSign*1000,dblPeakTime*1000));
				end
				hold off
				fixfig;
				
				if intPlot > 3
					vecHandles = get(gcf,'children');
					ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
					axes(ptrFirstSubplot);
					vecY = get(gca,'ylim');
					hold on;
					if intLatencyPeaks > 3,plot(dblOnset*[1 1],vecY,'r--');end
					plot(dblPeakTime*[1 1],vecY,'g--');
					plot(dblMaxDTime*[1 1],vecY,'b--');
					plot(dblMaxDTimeInvSign*[1 1],vecY,'b-.');
					hold off
				end
			end
		else
			%placeholder peak data
			sRate.dblOnset = [nan];
			vecLatencies = [nan nan nan nan];
			vecLatencyVals = [nan nan nan nan];
		end
	else
		vecLatencies = [];
		vecLatencyVals = [];
	end
	
	%check number of latencies
	if numel(vecLatencies) < intLatencyPeaks
		vecLatencies(end+1:intLatencyPeaks) = nan;
		vecLatencyVals(end+1:intLatencyPeaks) = nan;
	end
	
	%% build optional output structure
	if nargout > 2
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
		sZETA.vecRefT = vecRefT;
		sZETA.vecD = vecRealDiff;
		sZETA.matRandD = matRandDiff;
		
		sZETA.dblD_InvSign = dblD_InvSign;
		sZETA.dblPeakT_InvSign = dblMaxDTimeInvSign;
		sZETA.intPeakIdx_InvSign = intPeakLocInvSign;
		sZETA.dblUseMaxDur = dblUseMaxDur;
		sZETA.vecLatencyVals = vecLatencyVals;
	end
end

