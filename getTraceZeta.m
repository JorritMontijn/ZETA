function [dblZETA,vecLatencies,sZETA,sMSD] = getTraceZeta(vecTraceT,vecTraceAct,varEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,boolVerbose)
	%getTraceZeta Calculates neuronal responsiveness index zeta for traces
	%syntax: [dblZETA,vecLatencies,sZETA,sMSD] = getTraceZeta(vecTraceT,vecTraceAct,varEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,boolVerbose)
	%	input:
	%	- vecTraceT [N x 1]: time (s) corresponding to entries in vecTraceAct
	%	- vecTraceAct [N x 1]: activation trace (e.g., calcium imaging dF/F0)
	%	- vecEventTimes [T x 1]: event on times (s), or [T x 2] including event off times
	%	- dblUseMaxDur: float (s), ignore all values beyond this duration after stimulus onset
	%								[default: median of trial start to trial start]
	%	- intResampNum: integer, number of resamplings (default: 50)
	%	- intPlot: integer, plotting switch (0=none, 1=traces only, 2=activity heat map as well) (default: 0)
	%	- intLatencyPeaks: integer, maximum number of latency peaks to return (default: 4)
	%	- boolVerbose: boolean, switch to print messages
	%
	%	output:
	%	- dblZETA; Zenith of Event-based Time-locked Anomalies: FDR-corrected responsiveness z-score (i.e., >2 is significant)
	%	- vecLatencies; different latency estimates, number if determined by intLatencyPeaks:
	%		1) Latency of ZETA
	%		2) Latency of largest z-score with inverse sign to ZETA
	%		3) Peak time of multi-scale derivatives
	%		4...N) Onsets of largest sustained peaks
	%	- sZETA; structure with fields:
	%		- dblZ; uncorrected peak z-score 
	%		- dblP; p-value corresponding to zeta
	%		- dblPeakT; time corresponding to ZETA
	%		- dblHzD; Cohen's D based on mean-rate stim/base difference
	%		- dblHzP; p-value based on mean-rate stim/base difference
	%		- vecTraceT: timestamps of trace entries (corresponding to vecZ)
	%		- vecZ; z-score for all time points corresponding to vecTraceT
	%		- vecRealDiff: real offset of spikes relative to uniform rate
	%		- matRandDiff; matrix of shuffled runs with offset to uniform
	%	- sMSD; structure with fields: (only if intLatencyPeaks > 0)
	%		- vecMSD; Multi-scale derivative
	%		- vecScale; timescales used to calculate derivatives
	%		- matSmoothMSD; smoothed multi-scale derivatives matrix
	%		- matMSD; raw multi-scale derivatives matrix
	%		Additionally, it will return peak latencies using findsustainedpeaks.m:
	%		- vecTime: latencies for activation
	%		- vecIdx:  indices corresponding to peak locations in vecSpikeT/vecZ entries
	%		- vecDuration: peak durations
	%		- vecEnergy: peak energies
	%
	%Version history:
	%0.9 - June 27 2019
	%	Created by Jorrit Montijn
	%1.0 - September 24 2019
	%	New procedure to determine statistical significance [by JM]
	%2.0 - January 27 2020
	%	New peak detection procedure using multi-scale derivatives [by JM]
	
	%% prep data
	%ensure orientation
	vecTraceT = vecTraceT(:);
	vecTraceAct = vecTraceAct(:);
	
	%calculate stim/base difference?
	boolStopSupplied = false;
	dblHzD = nan;
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
		intLatencyPeaks = 2;
	end
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	
	%sampling frequency
	dblSamplingFreq = median(diff(vecTraceT));
	
	%% build onset/offset vectors
	vecEventStarts = varEventTimes(:,1);
	vecEventStops = varEventTimes(:,2);
	
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
	%define plots
	vecRandMean = nanmean(matRandDiff,2);
	vecRandSd = nanstd(matRandDiff,[],2);
	vecZ = ((vecRealDiff-mean(vecRandMean))./mean(vecRandSd));
	if numel(vecZ) < 3
		dblZETA = 0;
		sZETA = struct;
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate zeta');
		return
	end
	
	%find highest peak and retrieve value
	[dummy,intPeakLoc]= max(abs(vecZ));
	dblMaxZTime = vecTraceT(intPeakLoc);
	dblZ = vecZ(intPeakLoc);
	dblCorrectionFactor = 2/3.5;
	dblZETA = dblZ*dblCorrectionFactor;
	dblP=1-(normcdf(abs(dblZETA))-normcdf(-abs(dblZETA)));
	%find peak of inverse sign
	[dummy,intPeakLocInvSign] = max(-sign(dblZ)*vecZ);
	dblMaxZTimeInvSign = vecTraceT(intPeakLocInvSign);
	dblZ_InvSign = vecZ(intPeakLocInvSign);
	
	%build common timeframe
	vecRefT = (dblSamplingFreq/2):dblSamplingFreq:dblUseMaxDur;
		
	if boolStopSupplied
		%% calculate mean-rate difference
		%pre-allocate
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
		dblHzD = abs(nanmean(vecStimAct - vecBaseAct)) / ( (nanstd(vecStimAct) + std(vecBaseAct))/2);
		indUseTrials = ~isnan(vecStimAct) & ~isnan(vecBaseAct);
		[h,dblHzP]=ttest(vecStimAct(indUseTrials),vecBaseAct(indUseTrials));
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
		
		if intPlot == 2
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
		scatter(dblMaxZTime,vecRealDiff(intPeakLoc),'bx');
		scatter(dblMaxZTimeInvSign,vecRealDiff(intPeakLocInvSign),'b*');
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		if boolStopSupplied
			title(sprintf('ZETA=%.3f (p=%.3f), d(Hz)=%.3f (p=%.3f)',dblZETA,dblP,dblHzD,dblHzP));
		else
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,dblP));
		end
		fixfig
		
		
	end
	
	%% calculate MSD
	if intLatencyPeaks > 0
		[vecMSD,sMSD] = getMultiScaleDeriv(vecRefT,vecRealDiff,[],[],[],intPlot);
	else
		sMSD = [];
	end
	
	%% calculate MSD statistics
	if ~isempty(sMSD) && intLatencyPeaks > 0
		%get sustained peak onset
		[vecPeakTime,vecPeakIdx,vecPeakDuration,vecPeakEnergy,vecSlope] = findsustainedpeaks(vecMSD,vecRefT,intLatencyPeaks);
		%get MSD most prominent peak time
		[vecVals,vecLocs,vecsWidth,vecProms]=findpeaks(vecMSD);
		[dummy,intIdxMSD] = max(vecVals);
		intPeakLocMSD = vecLocs(intIdxMSD);
		dblPeakTimeMSD = vecTraceT(intPeakLocMSD);
		
		sMSD.intPeakLocMSD = intPeakLocMSD;
		sMSD.dblPeakTimeMSD = dblPeakTimeMSD;
		sMSD.vecTime = vecPeakTime;
		sMSD.vecIdx = vecPeakIdx;
		sMSD.vecDuration = vecPeakDuration;
		sMSD.vecEnergy = vecPeakEnergy;
		vecLatencies = [dblMaxZTime dblMaxZTimeInvSign dblPeakTimeMSD vecPeakTime(:)'];
		vecLatencies = vecLatencies(1:intLatencyPeaks);
		if intPlot > 0
			hold on
			scatter(vecPeakTime(1),vecMSD(vecPeakIdx(1)),'rx');
			scatter(dblPeakTimeMSD,vecMSD(intPeakLocMSD),'gx');
			scatter(dblMaxZTime,vecMSD(intPeakLoc),'bx');
			scatter(dblMaxZTimeInvSign,vecMSD(intPeakLocInvSign),'b*');
			hold off
			title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms,Sh=%.0fms',dblMaxZTime*1000,dblMaxZTimeInvSign*1000,dblPeakTimeMSD*1000,vecPeakTime(1)*1000));
			fixfig;
		

			vecHandles = get(gcf,'children');
			ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
			axes(ptrFirstSubplot);
			vecY = get(gca,'ylim');
			hold on;
			plot(vecPeakTime(1)*[1 1],vecY,'r--');
			plot(dblPeakTimeMSD*[1 1],vecY,'g--');
			plot(dblMaxZTime*[1 1],vecY,'b--');
			plot(dblMaxZTimeInvSign*[1 1],vecY,'b-.');
			
			hold off
		end
	else
		vecLatencies = [];
	end
		
	%% build optional output structure
	if nargin > 1
		sZETA = struct;
		sZETA.dblZ = dblZ;
		sZETA.dblP = dblP;
		sZETA.dblPeakT = dblMaxZTime;
		sZETA.intPeakIdx = intPeakLoc;
		if boolStopSupplied
			sZETA.dblHzD = dblHzD;
			sZETA.dblHzP = dblHzP;
		end
		sZETA.vecTraceT = vecTraceT;
		sZETA.vecZ = vecZ;
		sZETA.vecRealDiff = vecRealDiff;
		sZETA.matRandDiff = matRandDiff;
		
		sZETA.dblZ_InvSign = dblZ_InvSign;
		sZETA.intPeakIdx_InvSign = intPeakLocInvSign;
		sZETA.dblPeakT_InvSign = dblMaxZTimeInvSign;
	end
end

