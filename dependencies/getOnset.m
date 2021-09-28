function [dblOnset,dblValue,dblBaseVal,dblPeakT,dblPeakVal] = getOnset(vecData,vecT,dblPeakT,vecRestrictRange,intSwitchZ)
	%getOnset Returns peak onset. Syntax:
	%    [dblOnset,dblValue] = getOnset(vecData,vecT,dblPeakT,vecRestrictRange)
	%
	%Required input:
	%	- vecData [N x 1]: values
	%
	%Optional inputs:
	%	- vecT [N x 1]: timestamps corresponding to vecData (default: [1:N])
	%	- dblPeakT (float): timestamp corresponding to peak
	%	- vecRestrictRange [2 x 1]: restrict peak to lie within vecRestrictRange(1) and vecRestrictRange(end)
	%
	%Outputs:
	%	- dblOnset: time of peak onset (first crossing half-height of peak)
	%	- dblValue: value at peak onset
	%
	%Version history:
	%1.0 - February 26 2020
	%	Created by Jorrit Montijn
%%
	%check inputs
	if ~exist('vecT','var') || isempty(vecT)
		vecT = 1:numel(vecData);
	end
	if ~exist('vecRestrictRange','var') || isempty(vecRestrictRange)
		vecRestrictRange = [min(vecT) min(vecT)+range(vecT)];
	end
	if ~exist('intSwitchZ','var') || isempty(intSwitchZ)
		intSwitchZ = 1;
	end
	
	%z-score
	if intSwitchZ == 1
		vecDataZ = zscore(vecData);
	elseif intSwitchZ == 2
		dblMu = mean(vecData(vecT<0.02));
		vecDataZ = (vecData - dblMu)/std(vecData);
	else
		vecDataZ = vecData;
	end
	if ~exist('dblPeakT','var') || isempty(dblPeakT)
		[dummy,dblPeakT] = getPeak(vecDataZ,vecT,vecRestrictRange,0);
	end
	
	%remove time points outside restricted range
	indRemove = vecT < vecRestrictRange(1) | vecT > vecRestrictRange(end);
	vecCropT = vecT(~indRemove);
	vecDataZ = vecDataZ(~indRemove);
	
	%calculate first timepoint crossing half-height of peak 
	[dummy,intPeakIdx]=min(abs(vecCropT-dblPeakT));
	dblPeakVal = vecDataZ(intPeakIdx);
	dblBaseVal = vecDataZ(1);
	dblThresh = (dblPeakVal - dblBaseVal)/2 + dblBaseVal;
	if dblThresh > 0
		intOnsetIdx = find(vecDataZ >= dblThresh,1,'first');
	else
		intOnsetIdx = find(vecDataZ <= dblThresh,1,'first');
	end
	dblOnset = vecCropT(intOnsetIdx);
	dblValue = vecData(find(vecT > dblOnset,1));
	
	%check if empty
	if isempty(intOnsetIdx)
		dblOnset = nan;
		dblValue = nan;
	end
end
