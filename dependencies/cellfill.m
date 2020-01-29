function cellArray = cellfill(varContent,vecSize)
	%cellfill Fills cell array of size [vecSize] with content in
	%[varContent]; [varContent] can be any kind of data that can be
	%replicated using the function repmat
	%   Syntax: cellArray = cellfill(varContent,vecSize)
	
	
	%create data
	matPreCell = repmat(varContent,vecSize);
	
	%create syntax
	strSyntaxLoop = '';
	for intDim=1:length(vecSize)
		strSyntaxLoop = [strSyntaxLoop ',size(varContent,' num2str(intDim) ')*ones(1,vecSize(' num2str(intDim) '))'];
	end
	strLine = ['mat2cell(matPreCell' strSyntaxLoop ');'];
	
	%transform to cell array
	cellArray= eval(strLine);
end

