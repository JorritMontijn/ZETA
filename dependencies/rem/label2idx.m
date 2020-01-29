function [varDataOut,vecUnique,vecCounts,cellSelect,vecRepetition] = label2idx(varData)
	%label2idx Transforms label-entries to consecutive integer-based data
	%Syntax: [varDataOut,vecUnique,vecCounts,cellSelect,vecRepetition] = label2idx(varData)
	
	varDataOut = nan(size(varData));
	vecUnique = unique(varData);
	vecCounts = zeros(size(vecUnique));
	vecIdx = 1:length(vecUnique);
	cellSelect = cell(1,length(vecUnique));
	vecRepetition = zeros(1,length(vecUnique));
	for intIdx=vecIdx
		indEntries = varData==vecUnique(intIdx);
		cellSelect{intIdx} = indEntries;
		varDataOut(indEntries) = intIdx;
		vecCounts(intIdx) = sum(indEntries);
		vecRepetition(indEntries) = 1:sum(indEntries);
	end
end

