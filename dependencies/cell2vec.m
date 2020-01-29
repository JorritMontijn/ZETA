function [vecOut,cellSize] = cell2vec(cellIn)
	%UNTITLED7 Summary of this function goes here
	%   Detailed explanation goes here
	
	vecOut = [];
	for intCell=1:numel(cellIn)
		cellSize{intCell} = size(cellIn{intCell});
		vecOut = [vecOut;cellIn{intCell}(:)];
	end
end

