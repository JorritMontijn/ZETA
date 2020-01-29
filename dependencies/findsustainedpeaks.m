function [vecTime,vecIdx,vecDuration,vecEnergy,vecSlope] = findsustainedpeaks(vecData,vecT,intMaxPeaks,dblCutOffZ)
	%findsustainedpeaks Returns the left-hand onset time of sustained peaks. Syntax:
	%    [vecTime,vecIdx,vecDurations,vecEnergy,vecSlope] = findsustainedpeaks(vecData,vecT,intMaxPeaks,dblCutOffZ)
	%
	%Required input:
	%	- vecData [N x 1]: values
	%
	%Optional inputs:
	%	- vecT [N x 1]: timestamps corresponding to vecData
	%	- intMaxPeaks: maximum number of peaks to return
	%	- dblCutOffZ: removes all peaks under this z-score (default: 1);
	%					higher values improve performance, but lower values may be more	accurate
	%
	%Outputs:
	%	- vecTime: positive peak times
	%	- vecIdx: positive peak indices
	%	- vecDuration: positive peak durations
	%	- vecEnergy: positive peak energies
	%	- vecSlope: MSD slopes of vecData
	%
	%Note: requires getMultiScaleDeriv.m
	%
	%Version history:
	%1.0 - January 27 2019
	%	Created by Jorrit Montijn
	
	%% check inputs
	vecData=vecData(:);
	intN = numel(vecData);
	if ~exist('vecT','var')
		vecT = 1:intN;
	end
	if ~exist('intMaxPeaks','var')
		intMaxPeaks = 1;
	end
	if ~exist('dblCutOffZ','var')
		dblCutOffZ = 0;
	end
	
	%% calculate slope using MSD
	vecSlope = getMultiScaleDeriv(vecT,vecData);
	vecSlope = vecSlope(2:end);
	
	%% positive peaks
	[vecIdx,vecDuration,vecEnergy] = getSustainedPeaks(vecT,vecData,vecSlope,dblCutOffZ);
	
	%% remove minor peaks
	intKeepPosPeaks = min([intMaxPeaks numel(vecEnergy)]);
	vecEnergy = vecEnergy(1:intKeepPosPeaks);
	vecIdx = vecIdx(1:intKeepPosPeaks);
	vecDuration = vecDuration(1:intKeepPosPeaks);
	vecTime = vecT(vecIdx);
	
	% plot
	%plot(vecT,vecData);
	%hold on
	%scatter(vecT(vecPosIdx),vecData(vecPosIdx),'gx');
	%scatter(vecT(vecNegIdx),vecData(vecNegIdx),'rx');
	%hold off
end
