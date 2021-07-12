function vecOut = cell2vec(cellIn)
	%cell2vec Transforms cell array to vector. Syntax:
	%   vecOut = cell2vec(cellIn)
	
    try
        intElements = sum(cellfun(@numel,cellIn),'all');
    catch
        intElements = sum(cellfun(@numel,cellIn(:)));
    end
	vecOut = nan(intElements,1);
	intElCounter = 1;
	for intCell=1:numel(cellIn)
		vecData = cellIn{intCell}(:);
		intCellSize = numel(vecData);
		vecOut(intElCounter:(intElCounter+intCellSize-1)) = vecData;
		intElCounter = intElCounter + intCellSize;
	end
end

