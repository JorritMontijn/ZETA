function [cellOut] = vec2cell(vecIn,cellSize)
	%UNTITLED8 Summary of this function goes here
	%   Detailed explanation goes here
	if ~exist('cellSize','var') || isempty(cellSize),cellSize=cellfill([1 1],size(vecIn));end
		
	cellOut = [];
	intElCounter = 1;
	for intCell=1:numel(cellSize)
		vecSize = cellSize{intCell};
		intElements = prod(vecSize);
		vecElements = intElCounter:(intElCounter+intElements-1);
		intElCounter = intElCounter + intElements;
		
		matThisCell = reshape(vecIn(vecElements),vecSize);
		cellOut{intCell} = matThisCell;
	end
end

