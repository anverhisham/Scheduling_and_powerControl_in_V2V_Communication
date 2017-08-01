%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% This file contains GENERIC-FUNCTIONS(gf), which are too small to place in a seperate file %%%%%%%%%%%%%%% 
%
% -Function to add:
%     insert/deletion function for both cell & double array
%%%%% templocal

classdef gf

    methods(Static)

        
        %% -Addin paths for Generic functions, in case this program run in Glenn Cluster
        function addpaths()
            addpath(genpath('~/Programming/Matlab/MyMatlabFunctions/yalmip'));  %% -Adding Gurobi
            addpath('~/Programming/Matlab/MyMatlabFunctions/');
            addpath('~/Programming/Matlab/MyMatlabFunctions/gurobi_mex_v1.61');
            addpath('~/Programming/Matlab/ALLtoALL_Communication/Library');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -MATLAB INBUILT Functions tweaks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%- Variant of MATLAB inbuilt function 'error', the difference is that this function prints all inputs unlike MATLAB error function!
        function error(varargin) 
            gf.Disp('Error: '); gf.Disp(varargin); error(' ');
        end
        %% %%%%%%%%%%%- Function to display everything Passed -%%%%%%%%%%%%%%%%%%%%
        function Disp(varargin) 
            for iArg=1:nargin disp(varargin{iArg}); end;
        end
        
        %% Similar to MATLAB round, but here one can specify the precision
        %   gf.round([2,3,3.4,4.1,8.79],0.2)  ->  [2  3  3.4  4  8.79]
        function in = round(in,precision)
            inRounded = round(in);
            residues = in-inRounded;
            isCandidate = abs(residues) <= precision;
            in(isCandidate) = inRounded(isCandidate);
        end
        
        %% MATLAB 'pow2db' throws error upon NaN values. Here we return NaN itself
        function out = pow2db(in)
            out = nan(size(in));
            out(~isnan(in)) = pow2db(in(~isnan(in)));
        end
        
        %% MATLAB 'db2pow' throws error upon NaN values. Here we return NaN itself
        function out = db2pow(in)
            out = nan(size(in));
            out(~isnan(in)) = db2pow(in(~isnan(in)));
        end
        
        
        %% %%%%% Convert input to a Column vector %%%%%%%%%
        function out = vec(in)
           out = in(:);
        end
        
        %% Extension of gf.vec, but here we reduce the matrics to nDims, by appending along 'appendingDim'
        %   Note: gf.slack(in,1,1) is same as gf.vec(in)
        function out = slack(in,nDims,appendingDim)
           assert(appendingDim<=nDims,'Error in inputs');
           if(nDims==1) out = in(:); return; end
           if(nDims>=ndims(in)) out = in; return; end
           permuteDimOrder = 1:ndims(in);
           permuteDimOrder(nDims) = appendingDim;  permuteDimOrder(appendingDim) = nDims;
           tempin = permute(in,permuteDimOrder);
           tempinSizePerDim = size(tempin);
           tempoutSizePerDim = tempinSizePerDim(1:nDims);
           tempoutSizePerDim(nDims) = tempoutSizePerDim(nDims)*prod(tempinSizePerDim(nDims+1:end));
           tempout = reshape(tempin,tempoutSizePerDim);
           outPermuteDimOrder = 1:ndims(tempout);
           outPermuteDimOrder(nDims) = appendingDim;  outPermuteDimOrder(appendingDim) = nDims;
           out = permute(tempout,outPermuteDimOrder);           
        end
        
        %% Same as MATLAB size, but here we are computing 'Standard Size'.   Eg.[1,2,3]'  -> [3] (Not [3,1] as in MATLAB size)
        %   Definition:
        %       'Standard Size': ignore all trailing singleton dimensions, but return size as a row (instead of a column)
        function sizeV = size(in,dims)
            if(gf.isempty('dims'))
                sizeV = gf.sizeConvert(size(in));
            else
                sizeV = nan(1,numel(dims));
                for iDim = 1:numel(dims)
                    sizeV(iDim) = size(in,dims(iDim));
                end
            end
        end
        
        %% Function to convert matlab size vector ([3,1]) into 'Standard Size' vector ([3]), by removing trailing '1's
        function outSize = sizeConvert(inSize)
            isNotOne = inSize~=1;
            n = find(isNotOne,1,'last');
            if(isempty(n))
                n=1;
            end
            outSize = inSize(1:n);
        end
        
        
        %% Same as MATLAB ndims, but here it returns '1' for column vector 
        function nDim = ndims(in)
            nDim = ndims(in);
            if(nDim == 2 && size(in,2)==1)
               nDim = 1; 
            end
        end
        
        %% -Same as MATLAB num2str/str2num, except that this function can also accept array of cells as inputs!  
        function out = num2str(in)
            out = in;
            if(gf.isempty('in'))
                out = '';
                return;
            end
            if(isa(in,'cell'))
                if(numel(in)>1)
                    for i=1:numel(in)
                        out{i} = gf.num2str(in{i});
                    end
                else
                    out{1} = gf.num2str(in{1});
                end
            else        %% Array of doubles
                out = num2str(in);
            end
        end
        function [out,st] = str2num(in)     %% Advantage of str2num over str2double, is the presence of 2nd output 'st' :-) 
            out = in; st = in;
            if(gf.isempty('in'))
                out = '';
                st = 0;
                return;
            end
            if(isa(in,'cell'))
                if(numel(in)>1)
                    for i=1:numel(in)
                        [out{i},st{i}] = gf.str2num(in{i});
                    end
                    st = min(cell2mat(gf.vec(st)));
                else
                    [out{1},st] = gf.str2num(in{1});
                end
            else        %% Array of doubles
                [out,st] = str2num(in);
            end
        end
                
        
        %% Same as MATLAB strjoin, but this one concatenate equal sized in1 & in2
        function out = strjoin(in1,in2)
            if(numel(in1)==numel(in2))
                in1{end+1} = '';
            end
            out = strjoin(in1,in2);
        end
        
        
        %% Translate a literal string to latex format (so that it appears as such while using interpreter=latex)
        %   Eg:  'alg_Heuristic'   ->   'alg\_Heuristic'
        function out = translateToLatexString(in)
            out = regexprep(in,'_','\\_');
        end
        
        
        %% SAME AS MATLAB 'isempty', but input can be anything (including structure fields and non-existing variables) %%%%%%%
        %   NOTE: input must be name of a variable (string format)!
        function out = isempty(inputString)
            if(~isempty(strfind(inputString,':')) && ~isempty(regexp(inputString,':\d*}(?=\s*$)','once')))      %% For 'a{:}', we need to ensure that, this function is not tricked with a{1} alone.  NOTE: 'contains' not supported in old(Glenn) matlab versions :-(
                to = regexp(inputString,'}(?=\s*$)','once');
                from = gf.findMatchingParentheses(inputString,to);
                if(gf.findComplicatedParentheses(inputString(from:to),1))
                    inputString = ['{',inputString,'}'];
                end
            end
            try
                out = evalin('caller',['isempty(',inputString,')']);
            catch
               out = true;                
           end
        end
        
        %% -Same as MATLAB built in numel function except that this one supports multiple inputs
        function out = numel(varargin)
            if(nargin>1)
               out = nan(nargin,1);
               for i=1:nargin
                out(i) = gf.numel(varargin{i});
               end
            else
               out = numel(varargin{1});
            end
        end
        
        %% -Obsolete: Use isa(in,'integer')
        function out = isInteger(in)
            out = (rem(in,1)==0);
        end
        
        %% -Matlab variant of 'all'. Here we find out if all elements in multi-dimensional matrix are '1's
        function outBool = all(in)
            outBool = all(in(:));
        end
        
        %% -Matlab variant of 'all'. Here we find out if all elements in multi-dimensional matrix are '1's
        function outBool = any(in)
            outBool = any(in(:));
        end
        
        %% -Same as MATLAB 'sum', except that here you can specify list of dimensions, instead of single dimension :-)  
        function in = sum(in,dims)
            if(gf.isempty('dims')) dims = 1:numel(size(in)); end
            if(isa(in,'sym')) fh = @gf.sumSymbolic; else fh = @sum; end
            for dim = gf.convertColumnVectorsToRowVectors(dims)
                in1=fh(in,dim);
                in = in1;
            end
        end
        
        %% -Same as MATLAB 'nansum', except that here you can specify list of dimensions, instead of single dimension :-)  
        function in = nansum(in,dims)
            if(gf.isempty('dims')) dims = 1:numel(size(in)); end
            if(isa(in,'sym')) fh = @gf.sumSymbolic; else fh = @nansum; end
            for dim = gf.convertColumnVectorsToRowVectors(dims)
                in1=fh(in,dim);
                in = in1;
            end
        end
        
        
        %% -Same as MATLAB 'mean', except that here you can specify list of dimensions, instead of single dimension :-)  
        function in = mean(in,dims)
            if(gf.isempty('dims')) dims = 1:ndims(in); end
            for dim = gf.convertColumnVectorsToRowVectors(dims)
                in = mean(in,dim);
            end
        end
        
        %% MATLAB's inbuilt 'sum' function doesn't support summing along 3rd or 4th ... dimensions. Hence this function.
        function M = sumSymbolic(M,dim)
            s=size(M);
            if(numel(s)<dim) s(dim) = 1; end
            M=permute(M,[setdiff(1:ndims(M),dim),dim]);
            M=reshape(M,[],s(dim));
            M=sum(M,2);
            s(dim)=1;
            M=reshape(M,s);
        end
        
        %% Find max of all values in a matrix
        function in = max(in,dims)
            if(gf.isempty('dims')) dims = 1:ndims(in); end
            for dim = gf.convertColumnVectorsToRowVectors(dims)
                in = max(in,[],dim);
            end
        end
        
        %% Find max of all values in a matrix
        function in = nanmax(in,dims)
            if(gf.isempty('dims')) dims = 1:ndims(in); end
            for dim = gf.convertColumnVectorsToRowVectors(dims)
                in = nanmax(in,[],dim);
            end
        end
        
        %% Find max of all values in a matrix
        function in = min(in,dims)
            if(gf.isempty('dims')) dims = 1:ndims(in); end
            for dim = gf.convertColumnVectorsToRowVectors(dims)
                in = min(in,[],dim);
            end
        end

        %% Find max of all values in a matrix
        function in = nanmin(in,dims)
            if(gf.isempty('dims')) dims = 1:ndims(in); end
            for dim = gf.convertColumnVectorsToRowVectors(dims)
                in = nanmin(in,[],dim);
            end
        end


        %% Same as MATLAB inbuilt 'squeeze' except that this will squeeze cell array too
        function in =  squeeze(in)
            if(isa(in,'cell'))
               indicesTodelete = [];
               for i = 1:numel(in)
                   if(isempty(in{i}))
                       indicesTodelete = [indicesTodelete,i];
                   end
               end
               in(indicesTodelete) = [];
            else
               in = squeeze(in); 
            end
        end
        
        %% Same as MATLAB inbuilt 'squeeze' except that this function also squeezes row vector to column vector 
        function out =  squeezeAndTransposeToRowVector(in)
            out = squeeze(in);
            if(numel(size(in))==2 && size(in,1)==1)
                out = in.';
            end
        end
        
        %% Similar to matlab ind2sub, but here 'indicesPerDim' is an array of cells, with each cell containing a set of indices 
        %   ind2sub([2,3],5)  <=>  gf.ind2sub({1:2,1:3},5)
        function varargout = ind2sub(indicesPerDim,index)
            lengthPerDimension = cell2mat(cellfun(@numel,indicesPerDim,'uni',false));
            out = cell(1,numel(lengthPerDimension));
            [out{:}] = ind2sub(lengthPerDimension,index);
            for iDimension = 1:numel(indicesPerDim)
                varargout{iDimension} = indicesPerDim{iDimension}(out{iDimension});
            end
        end
        
        %% Same as MATLAB inbuilt 'cell2mat' except that this function compensate for incompatible sizes of c{i} by placing 'Nan'  
        %   'in' is a matrix of cells, with each cell containing a matrix along high dimensions.
        % -Algorithm:
        %   1. Find the lowest dimension 'currentDim' with length >1
        %   2. Split the 'in' along the dimension 'currentDim', and call cell2mat recursively for each sub 'in' 
        %   3. Concatenate all the results got in the above step. 
        function out = cell2mat(in)
            nDim = numel(size(in));
            currentDim = find(size(in)>1,1);
            if(~isempty(currentDim))
                exString = repmat(':,',1,nDim);     %% Create string as ':,:,12,:,:'
                exString(end) = [];         %% Removing trailing ','
                outSub = cell(size(in,currentDim),1);
                for i = 1:size(in,currentDim)
                    exStringCurrent = [exString(1:2*(currentDim-1)),num2str(i),exString(2*currentDim:end)];
                    outSub{i} = gf.cell2mat(eval(['in(',exStringCurrent,')']));
                end
                out = gf.cat(currentDim,outSub{:});
            else
                out = in{1};
            end
        end
        
        %% Same as MATLAB inbuilt 'cell2mat' except that this function concatenate incompatible sized by appending 'Nan's 
        %   gf.cat(1,A,B,'-appendValue',123) will append '123' for missing values 
        function out = cat(dim,varargin)
            optionIndex = gf.contains(varargin,'-appendValue');
            if(~isempty(optionIndex))
                optionIndex = optionIndex{1}(1);
                appendValue = varargin{optionIndex+1};
                varargin(optionIndex:optionIndex+1) = [];
            else
                appendValue = nan;
            end
            sizePerInput = [];
            for i = 1:numel(varargin)
                if(~isempty(sizePerInput) && size(sizePerInput,1) ~= numel(size(varargin{i})))
                    sizePerInput = gf.cat(2,sizePerInput,size(varargin{i})'); %% Used gf.cat instead of cat, since varargin{n} can have higher dimensions than varargin{1},varargin{2},... varargin{n-1}
                else
                    sizePerInput = cat(2,sizePerInput,size(varargin{i})');    %% To avoid inifinite recursion
                end
            end
            maxSizePerDim = nanmax(sizePerInput,[],2);  
            in = cell(numel(varargin),1);
            for i = 1:numel(varargin)
                desiredSizePerDim = maxSizePerDim;
                desiredSizePerDim(dim) = size(varargin{i},dim);
                if(isnan(appendValue))
                    in{i} = gf.augmentMatrix(varargin{i},desiredSizePerDim);                %% For speeding up
                else
                    in{i} = gf.augmentMatrix(varargin{i},desiredSizePerDim,appendValue);
                end
            end
            out = cat(dim,in{:});
        end
        
        %% Enlarge a numeric matrix by adding 'NaN's 
        %  NOTES:
        %       1. NaN values 'outSize' indicate, don't touch those dimensions.
        %       1. 'isIgnoreSmallerOutSize' = true  =>  Program ignore ith dimension for which outSize(i) < size(in,i)
        %   Inputs: 
        %       1. 'in' must be a matrix, shouldn't contain any cells 
        %       2. 'outSize' must be a column vector, and in 'Matlab size' format (Not 'Standard Size', ie. [3,1] instead of [3])
        %               NaN in 'outSize' indicates "do not touch" that dimension
        %       3. 'outIndicesPerDim' must be an array of cells. outIndicesPerDim{i} indicates the indices in output to consider while filling along ith dimension.
        %               A cell containing ':' indicates, no reordering along corresponding dimension.
        %   Warnings:
        %       1. Changing default inputs affect gf.cat()
        function out = augmentMatrix(in,outSize,appendValue,isIgnoreSmallerOutSize,outIndicesPerDim)
            if(nargin<3 || gf.isempty('appendValue')) appendValue = nan; end
            if(nargin<4 || gf.isempty('isIgnoreSmallerOutSize')) isIgnoreSmallerOutSize = false; end
            % -Rectify inputted sizes
            outSize = gf.convertRowVectorsToColumnVectors(outSize);
            inSize = size(in)';
            if(numel(inSize)<numel(outSize))
                inSize = gf.augmentMatrix(inSize,size(outSize),1);
            elseif(numel(inSize)>numel(outSize))
                outSize = gf.augmentMatrix(outSize,size(inSize),1);
            end
            outSize(isnan(outSize)) = inSize(isnan(outSize));
            if(isIgnoreSmallerOutSize)
                outSize(outSize<inSize) = inSize(outSize<inSize);
            end
            if(nargin<5 || gf.isempty('outIndicesPerDim')) outIndicesPerDim = num2cell(repmat(':',numel(inSize),1)); end
            for iDim = 1:numel(inSize)
                if(strcmp(outIndicesPerDim{iDim},':'))
                    outIndicesPerDim{iDim} = [1:inSize(iDim)];
                end
            end
            % -Verification
            assert(all(inSize==gf.numel(outIndicesPerDim{:})),'Error: outSize must match to the corresponding number of element in ''outIndicesPerDim'' !!')
            assert(numel(outSize)>=2,'Error: ''outSize'' must be in MATLAB Format!!');
            assert(all(inSize<=outSize),'Error: Desired out size is less than ''in'' size!!');
            % -Assigning to output
            out = appendValue*ones(outSize');
            out(outIndicesPerDim{:}) = in;
        end
        
        %% -Returns true, if input is a column vector
        function out = isColumnVector(in)
            if(numel(size(in))==2 && size(in,2)==1)   %% If in is a column vector
                out = true;
            else
                out = false;
            end
        end
        
        %% -Returns true, if input is a column vector
        function out = isRowVector(in)
            if(numel(size(in))==2 && size(in,1)==1)   %% If in is a row vector
                out = true;
            else
                out = false;
            end
        end
        
        %% -Function to quantize a value to {offset, offset+stepSize, offset+2*stepSize, ...etc}
        function out = quantize(in,stepSize,offset)
            if(gf.isempty('stepSize')) stepSize = 1; end    %% By default, quantize to integers.
            if(gf.isempty('offset')) offset = 0; end
            in = in - offset;
            inSign = sign(in);
            in = abs(in);
            remainder = rem(in,stepSize);
            count = floor(in/stepSize);
            out = (count+(remainder>=stepSize/2))*stepSize;
            out = inSign.*out;
            out = out + offset;
        end
        
        %% Create a unique positive real number 'c' from two natural number 'a' and 'b'
        % Ref: Cantor pairing function
        function c = pairNaturalNumbers(varargin)
            if(nargin>2)
                a = gf.pairNaturalNumbers(varargin{1:end-1});
                b = varargin{end};
                c = gf.pairNaturalNumbers(a,b);
            else
                a = varargin{1};
                b = varargin{2};
                c = 0.5*(a+b)*(a+b+1)+b;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Containers- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% -Container Map
        %   Example: algorithms = gf.Map('alg1','alg2'); Now algorithms('alg2')=2     
        %   To check the presence of a key =>  'isKey(algorithms,'alg3') = 0'
        function out = Map(varargin)
            valueSet = 1:numel(varargin);
            out = containers.Map(varargin,valueSet);
        end
        
        %%
        % -Get map values for multiple keys
        %
        function varargout = getMapValues(inMap,varargin)
            if(nargin>2)
                for i=1:nargin
                    varargout{i} = gf.getMapValues(varargin{i});
                end
            else
                if(isa(varargin{1},'cell'))
                    keys = varargin{1};
                    out = cell(size(keys));
                    for i=1:numel(keys)
                        out{i} = gf.getMapValues(inMap,keys{i});
                    end
                else
                    out = inMap(varargin{1});
                end
                varargout{1} = out;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Vector Functions- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% -Matlab computes norm of a matrix, whereas this function computes norm of each column, and returns a column vector (after transposing)
        function out = norm(in)
            assert(numel(size(in))==2,'Input Matrix must be 2 Dimensional')
            out = nan(size(in,2),1);
            for iCol = 1:size(in,2)
               out(iCol) = norm(in(:,iCol)); 
            end
        end
        
        
        %% -Sort w.r.t specified weights
        %   2016/07/15: Renamed function from 'weightedSort' to 'sortWeighted'  
        function outVector = sortWeighted(inVector,weights)
            [~,I] = sort(weights);
            outVector = inVector(I);
        end
        
        %% Same as MATLAB's reshape, but here we trim/expand input matrix, if total number of elements in desired size is not matching with number of elements in input
        function out = reshape(in,varargin)
            % -Get '-appendValue'
            if(nargin>=4)       %% For runtime optimization.
                optionIndex = gf.contains(varargin,'-appendValue');
            else
                optionIndex = [];
            end
            if(~isempty(optionIndex))
                optionIndex = optionIndex{1}(1);
                appendValue = varargin{optionIndex+1};
                varargin(optionIndex:optionIndex+1) = [];
            else
                appendValue = nan;
            end        
            % -Reshape
            outNumel = prod([varargin{:}]);
            inNumel = numel(in);
            if(inNumel<outNumel) 
                in = [in(:);appendValue*ones(outNumel-inNumel,1)]; 
            elseif(inNumel>outNumel) 
                in = in(1:outNumel); 
            end
            out = reshape(in,varargin{:});
        end
        
        
        %% %%%%% Make input size same as the size specified in 'outSize', by circular affixing
        %         Note: If 'outSize' is less than input size, then crop input 
        function out = reshapeByCircularAffix(in,outSize)
            inSize = [size(in),ones(1,numel(outSize)-numel(size(in)))];
            repPerDim = ceil(outSize./inSize(1:numel(outSize)));    %% - Repetition factor per dimension
            inExpanded = repmat(in,repPerDim);
            %% -Create cropExpression
            cropExpression = 'inExpanded(';
            for iDim = 1:numel(outSize)
                cropExpression = [cropExpression,'1:',num2str(outSize(iDim)),','];
            end
            cropExpression(end)=[];                                 %% -Removing last ','
            cropExpression = [cropExpression,')'];
            %% -Evaluate cropExpression
            out = eval(cropExpression);
        end
        
        
        
        %% -Bound input values within the limit specified,
        %   Note: If 'maxAmplitude' is a vector/matrix, then it's applied column/element individually :-) 
        function out = boundAmplitude(in,maxAmplitude)
            if(~gf.isempty('maxAmplitude'))
                outAbs = abs(in);
                if(numel(maxAmplitude)==1)
                    outAbs(outAbs>maxAmplitude) = maxAmplitude;                    
                else
                    maxAmplitude = gf.reshapeByCircularAffix(maxAmplitude,size(in));    %% -Making the size equal
                    diffMatrix = outAbs - maxAmplitude;
                    outAbs(diffMatrix>0) = outAbs(diffMatrix>0) - diffMatrix(diffMatrix>0);
                end
                out = outAbs.*exp(1i*angle(in));
            end          
        end
        
        %% - Equivalent of MATLAB min function, but here we return random index, when multiple min values are same!
        function [maxValue,index] = maxRandIndex(varargin)
            [maxValue,~] = max(varargin{:});
            isMaxValue = bsxfun(@ge,varargin{1},maxValue);      %% Placed '>=' instead of '=~' to avoid any rounding off error.
            isMaxValue = isMaxValue + 0.9*rand(size(isMaxValue));
            [~,index] = max(isMaxValue,varargin{2:end});
        end
        
        %% - Equivalent of MATLAB min function, but here we return random index, when multiple min values are same!
        function [minValue,index] = minRandIndex(varargin)
            [minValue,~] = min(varargin{:});
            isNotMinValue = bsxfun(@gt,varargin{1},minValue);      %% Placed '>' instead of '=~' to avoid any rounding off error.
            isNotMinValue = isNotMinValue + 0.9*rand(size(isNotMinValue));
            [~,index] = min(isNotMinValue,varargin{2:end});
        end
        
% %         ** FOLLOWING FUNCTIONS GOT OBSOLETED by above functions **      
% %         %% %%%%% -If multiple max/min values exists, then output any of the indices of max value - %%%%%%%%%%
% %         function [indices,maxValues] = maxRandIndex(input,~,dim)
% %             if(nargin<3 || strcmp(dim,''))
% %                 dim = 1;
% %             end
% %             inputSize = size(input);
% %             outputSize = inputSize; outputSize(dim) = 1;
% %             indices = nan(outputSize);
% %             permutedDimensions = [dim, setdiff(1:numel(inputSize),dim)];
% %             %% Converting input to 2D-matrix with 1st dimension as 'dim' dimension of input
% %             modifiedInput = reshape(permute(input,permutedDimensions),inputSize(dim),[]);
% %             maxValues = max(modifiedInput);            
% %             for iCol=1:size(modifiedInput,2)
% %                 indices(iCol) = gf.randsample(find(modifiedInput(:,iCol)==maxValues(iCol)),1);
% %             end
% %         end
% %         
% %         function [indices,minValues] = minRandIndex(input,~,dim)
% %             if(nargin<3 || strcmp(dim,''))
% %                 dim = 1;
% %             end
% %             inputSize = size(input);
% %             outputSize = inputSize; outputSize(dim) = 1;
% %             indices = nan(outputSize);
% %             permutedDimensions = [dim, setdiff(1:numel(inputSize),dim)];
% %             %% Converting input to 2D-matrix with 1st dimension as 'dim' dimension of input
% %             modifiedInput = reshape(permute(input,permutedDimensions),inputSize(dim),[]);
% %             minValues = min(modifiedInput);            
% %             for iCol=1:size(modifiedInput,2)
% %                 indices(iCol) = gf.randsample(find(modifiedInput(:,iCol)==minValues(iCol)),1);
% %             end
% %         end
        
        %% -Permute the elements in inVector
        function out = randPermute(in)
            out = in(randperm(numel(in)));
            out = reshape(out,size(in));
        end
        
        %% -Matlab randsample function output a number <= sampleInputs, if sampleInputs is a scalar.
        %    Example matlab 'randsample(7,1)' outputs a number from [1:7], whereas we need output '7'
        function outputSample = randsample(sampleInputs,numSamples,varargin)
            if(~isempty(varargin))
               replacement = varargin{1};
            else
                replacement = false;
            end
            if(numel(sampleInputs)==0)
                outputSample = nan(1,numSamples);
            elseif(numel(sampleInputs)==1)
                outputSample = randsample([sampleInputs,sampleInputs],numSamples,replacement);
            else
                outputSample = randsample(sampleInputs,numSamples,replacement);
            end
        end
        
        %% Let's a vector 'x' is upsampled to 'xu', this function gives the newIndices in the upsampled output xu
        %   Assumption: 'oldIndices' is a column vector
        function newIndices = upsampleIndex(oldIndices,upsamplingFactor)
            assert(size(oldIndices,2)==1,'Error(gf.upsampleIndex): oldIndices must be column vector!')
            tempM = ones(numel(oldIndices),1)*[upsamplingFactor-1:-1:0];
            newIndices = upsamplingFactor*oldIndices;
            newIndices = bsxfun(@minus,newIndices,tempM)';
            newIndices = newIndices(:);
        end

        %% Like MATLAB inbuilt find, but here it returns indices (and only indices) for multiple dimensions input
        %   i.e., varargout{i} corresponds to indices along i^th dimension
        function varargout = findIndices(in,dims)
            if(nargin<2 || gf.isempty('dims'))
                if(nargout==1 && gf.ndims(in)==1)   %% If 'in' column vector, then single output must be supported
                    dims = [1];
                else
                    dims = 1:nargout;
                end
            end
            assert(nargout==numel(dims),'Error: Number of outputs must be equal or greater than ndims(in) !!');
            indicesPerDim = cell(ndims(in),1);
            [indicesPerDim{:}] = ind2sub(size(in),find(in(:)));
            if(max(dims)>numel(indicesPerDim))      %% Appending for trailing dimensions
                [indicesPerDim{numel(indicesPerDim)+1:max(dims)}] = deal(ones(numel(indicesPerDim{1}),1));
            end
            varargout = indicesPerDim(dims);
        end
        
        %% 'outIndices' is having same size as 'valuesToSearch'
        %  outIndices(i) is the first index of valuesToSearch(i)
        function outIndices = findIndicesOfSearchValues(referenceVector,valuesToSearch)
            outIndices = nan(size(valuesToSearch));
            for iValueToSearch=1:numel(valuesToSearch)
                outIndex = find(referenceVector==valuesToSearch(iValueToSearch),1,'first');
                if(~isempty(outIndex))
                    outIndices(iValueToSearch) = outIndex;
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -MATRIX CREATION FUNCTIONS- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %% - out=getRowSumMatrix(prod(size(X)), =>  out*X(:) = sum(X,1)
        function out = getRowSumMatrix(matrixSize)
            nRows = matrixSize(1);
            nCols = matrixSize(2);
            out = kron(ones(1,nCols),eye(nRows));
        end
        
        %% - out=getColSumMatrix(prod(size(X)), =>  out*X(:) = sum(X,2)
        function out = getColSumMatrix(matrixSize)
            nRows = matrixSize(1);
            nCols = matrixSize(2);
            out = kron(eye(nRows),ones(1,nCols));
        end
        
        %% -If input is ([1,2]',[2,3]') then => out = [1,2;1,3;2,2;2,3]'
        function out = createAllCombinationOfIntegerVectors(startVector,endVector)
            assert(~any(size(startVector)-size(endVector)));         %% -Both inputs must be of same size
            assert(size(startVector,2)==1 && size(endVector,2)==1); %% -Both inputs must be column vectors
            if(numel(startVector)>1)            %% -RECURSIVE CALL
                subOut = gf.createAllCombinationOfIntegerVectors(startVector(2:end),endVector(2:end));
                topValueOriginal = startVector(1):endVector(1);
                topValues = ones(size(subOut,2),1)*topValueOriginal;
                out = [topValues(:)';repmat(subOut,[1,numel(topValueOriginal)])];
            else                                %% -TERMINATION
                out = startVector:endVector;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -MATRIX CHECKING FUNCTIONS- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        function out = isPermutationMatrix(in)
            out = gf.all(size(in,1)==size(in,2)) & gf.all(in>=0) & gf.all(sum(in,1)==1) & gf.all(sum(in,2)==1) & gf.all(sum(in.^2,1)==1);
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -MATRIX MANIPULATION FUNCTIONS- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %%%%%%%% Covert [2,3] -> [1,1,2,2,2] %%%%%%%%
        function out = exapandWeightsToObjects(in)
            out = cell2mat(arrayfun(@(x) x*ones(1,in(x)), 1:numel(in), 'UniformOutput', false));
        end
        
        %%%%%%% Delete the all nan values & return it %%%%%%%%%
        function out = noNaNs(in)
           out = in(~isnan(in)); 
        end 
        
        %%%%%%%%%% Upper bound all values in input %%%%%%%%%%%%
        %   Warning:
        %       1. There's no check for overshooting size(in) > size(upperBound)
        function out = ceil(in,upperBound)
            out = in;
            if(numel(upperBound)==1)
                out(in>upperBound) = upperBound;
            else
                if(~isequal(size(in),size(upperBound)))
                   upperBound = gf.reshapeByCircularAffix(upperBound,size(in));
                end
                out(in>upperBound) = upperBound(in>upperBound);
            end
        end
        
        %%%%%%%%%% Lower bound all values in input %%%%%%%%%%%%
        function out = floor(in,lowerBound)
            out = in;
            out(out<lowerBound) = lowerBound;
        end
        
        %%%%%%% Delete the specified elements in input & return it %%%%%%%%%
        function out = deleteElements(in,elementsToDelete)
            out = in;
            for element=elementsToDelete
                out(out==element) = [];
            end
        end 
        
        %%%%%%% Make all non-NaN values '0' & return it %%%%%%%%%
        function out = clearNonNaNs(in)
            out = nan(size(in));
            out(~isnan(in)) = 0; 
        end 
        
        %% -Matlab cat throws error while concatenating incompatible matrices, whereas following cat function concatenate after appending 'Nan's
        function out = myCat(dim,M1,M2)
            M1s = size(M1);
            M2s = size(M2);
            %% -Equalize the size vectors by appending '1's
            if(numel(M1s)>numel(M2s)) 
                M2s = [M1s,ones(1,(numel(M1s)-numel(M2s)))]; 
            elseif(numel(M1s)<numel(M2s)) 
                M1s = [M2s,ones(1,(numel(M2s)-numel(M1s)))]; 
            end
            %% -Equalizing all dimensions except concatenating dimension, of M1 & M2 by appending NaN
            for iDim = setdiff(1:numel(M1s),dim)
                if(M1s(iDim)>M2s(iDim))
                    S = M2s;
                    S(iDim) = M1s(iDim) - M2s(iDim);
                    M2 = cat(iDim,M2,nan(S));
                elseif(M1s(iDim)<M2s(iDim))
                    S = M1s;
                    S(iDim) = M2s(iDim) - M1s(iDim);
                    M1 = cat(iDim,M1,nan(S));
                end
            end
            %% -Do the final concatenation now,,
            out = cat(dim,M1,M2);
        end
        
        %% -Same as MATLAB intersect, except that here we chack whether an element is present in 'thresholdPercentage'
        %  -Warning: Each of varargin must be column vectors, (Also expect different behaviour for repetitive values!)
        function out = intersectPercentage(thresholdPercentage,varargin)
            numInputs = numel(varargin);
            input = cat(1,varargin{:});
            inputUnique = unique(input);
            countPerInUnique = histcounts(input,[inputUnique;max(inputUnique)+1]);
            out = inputUnique(countPerInUnique>=numInputs*(thresholdPercentage/100));
        end
        
        %% -Output only the values with some frequency(count) inside input vector.
        function out = getHighFrequentValues(frequencyThreshold,input)
            if(isempty(input))
                out = [];
            else
                inputUnique = unique(input);
                countPerInUnique = histcounts(input,[inputUnique;max(inputUnique)+1]);
                out = inputUnique(countPerInUnique>=frequencyThreshold);    
            end
        end
        
        %% - Same as MATLAB built-in 'intersect' function, except that multiple inputs are supported!
        function out = intersect(varargin)
            out = intersect(varargin{1},varargin{2});
            for iInput=3:nargin
                out = intersect(out,varargin{iInput});
            end
        end
        
        %% -Same as MATLAB built in setdiff, except that it preserves repeated values in first input!
        function out = setdiffNonUnique(in1,in2)
            out = in1(~ismember(in1,in2));
        end
        
        %%%%%%% -Convert input Matrix to a binary matrix- %%%%%%%%%%%%%
        %%%%%%% Note: nan vallues are preserved as such.. %%%%%%%%%%%%%
        function output = convertToBinaryFormat(input,binaryFormat)
            %%%%%% -Verifying/Correct Inputs- %%%%%%
            if(nargin<2) binaryFormat = [0,1]; end
            if(~isa(binaryFormat,'double') || numel(binaryFormat)~=2) gf.error('binaryFormat must be a vector of double!! binaryFormat = ',binaryFormat); end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            output = input;
            output(input==0) = binaryFormat(1);
            output(input~=0 & ~isnan(input)) = binaryFormat(2);
        end
        
        %%%%%%% Returns the first DimensionIndices, for which ALL elements in submatrix of inputMatrix is non-NaN values %%%%%%%%%
        function output = findNoNanIndicesForFirstDimension(inputMatrix)
            for iDimension=fliplr(2:numel(size(inputMatrix)))
               inputMatrix = sum(inputMatrix,iDimension);
            end
            output = find(~isnan(inputMatrix));
        end
        
        %%%%%% Make both input matrices with same size, by cyclic suffixing the smaller matrix
        function [out1,out2] = equalizeMatrices(in1,in2)
            out1=in1; out2=in2;
            if(numel(in1)<numel(in2))
                out1=nan(size(in2));
                out1(:) = in1(rem(1:numel(out1),numel(in1))+1);
            elseif(numel(in1)>numel(in2))
                out2=nan(size(in1));
                out2(:) = in2(rem(1:numel(out2),numel(in2))+1);                
            end
        end
        
        %%%%%% -Extract values from inputMatrix, & place it in outputMatrix in the specifiedDiemnsios
        %%%%%% -NOTE: Currently supported upto 3D matrix only ...
        function outputMatrix = extractFromMatrix(inputMatrix,varargin)
            % -Ensuring all input indices are vectors (not matrices)
            for iInput = 1:numel(varargin)
                if(sum(size(varargin{iInput})>1)>1) error('Error(extractFromMatrix): All index inputs must be vectors !!! '); end
            end
            
            inputDimensions = 1:numel(size(inputMatrix));
            outputDimensions = inputDimensions;
            for iInput = 1:numel(varargin)
               [~,outputDimensions(iInput)] = max(size(varargin{iInput}));  %% -Dimension of output is the direction(dimension) of each input index argument 
            end
            % -Trim input matrix as per indices
            trailingDimensions = num2cell(repmat(':',1,numel(size(inputMatrix))-numel(varargin)));
            inputMatrix = inputMatrix(varargin{:},trailingDimensions{:});  %% -Trimming ...
            %% -Do compression, looking at duplicate dimensions in input, (Assumption: Max 3D matrix)
            [~,compressionDimensions] = gf.getDuplicateValues(outputDimensions);
            compressionDimensions = find(compressionDimensions);
            %%%% TODO: Nowe rotate the matrix shape &  pick the diagonal
            %%%% elements for each 2D matrix ...
            
            %%%% -TODO: Rotate & scale compress the dimensios...
            
        end
        
        %% -Input is a 2D matrix,
        function in = removeSymmetricColumns(in)
            iCol = 1;
            while iCol <= size(in,2)
                colToRemove = find(sum(abs(bsxfun(@minus,in,in(:,iCol))))==0);
                if(numel(colToRemove)>1)            %% -Ofcourse the first element = iCol
                   in(:,colToRemove(2:end))= [];
                end
                iCol = iCol+1;
            end
        end
        
        %% -Divide each element of input matrix, by row & column sum
        %   Example: [1,2;3,4]  ->  [1/6,2/7;3/8,4/9]
        function out = normalizeWithColumnRowSummation(in)
            W = bsxfun(@plus,ones(1,size(in,1))*in,in*ones(size(in,2),1)); %% - Column sum + Row sum 
            out = in./(W-in);
        end
        
        %% -[1,2,3,1,4,2]  ->  [1,1,1, 2,1,2]
        function out = findOrderOfOccurrence(in)
            out = nan(size(in));
            inMembers = unique(in);
            for i=1:numel(inMembers)
               logicalIndices = ismember(in,inMembers(i));
               out(logicalIndices) = cumsum(logicalIndices(logicalIndices==true));
            end
        end
        
        %% Squeeze 'in' matrix along the singleton dimension of 'preservingIndices', using indices kept in 'preservingIndices'
        %
        function out = squeezeWithIndices(in,preservingIndices)
            % -Permute inputs
            assert(sum(size(in)>size(preservingIndices))<=1,'Error: preservingIndices must have same shape as in, except along one of the dimensions!!');
            dim = find(size(in)>size(preservingIndices));
            if(size(in,dim)==1) return; end
            assert(size(preservingIndices,dim)==1,'Error: preservingIndices must have unit lenght along squeezinf dimension');
            permOrder = 1:ndims(in); permOrder([1,dim]) = [dim,1];
            in = permute(in,permOrder);
            preservingIndices = permute(preservingIndices,permOrder);
            % -Find the gloabl preserving indices
            out = nan(size(preservingIndices));
            sizeI = size(in);
            modifiedIndices = reshape(0:prod(sizeI(2:end))-1,[1,sizeI(2:end)]);
            modifiedIndices = repmat(size(in,1)*modifiedIndices,size(preservingIndices,1),1) + preservingIndices;
            modifiedIndices = gf.noNaNs(modifiedIndices);
            % -Assign values in relevant indices to 'out'
            out = in(modifiedIndices);
            out = ipermute(out,permOrder);    
        end        
        
        %% NOTE:
        %   1. 0/NaN in 'preservingIndices' indicates ignoring indices
        %  DESCRIPTION:
        %   1. Look at 'in' and 'preservingIndices' column by column (or as per 'dim'). Preserve the values of 'in', whose indices are included in 'preservingIndices'.
        %       Clear remaining locations in 'in' using 'clearingValue'
        function out = preserveDataOnSpecifiedIndices(in,preservingIndices,clearingValue,dim)
            % -Permute inputs
            if(gf.isempty('dim')) dim = 1; end
            preservingIndices(preservingIndices==0) = nan;
            permOrder = 1:ndims(in); permOrder([1,dim]) = [dim,1];
            in = permute(in,permOrder);
            preservingIndices = permute(preservingIndices,permOrder);
            % -Find the gloabl preserving indices
            sizeI = size(in);
            modifiedIndices = reshape(0:prod(sizeI(2:end))-1,[1,sizeI(2:end)]);
            modifiedIndices = repmat(size(in,1)*modifiedIndices,size(preservingIndices,1),1) + preservingIndices;
            modifiedIndices = gf.noNaNs(modifiedIndices);
            % -Assign values in relevant indices to 'out'
            out = clearingValue*ones(size(in));
            out(modifiedIndices) = in(modifiedIndices);
            out = ipermute(out,permOrder);    
        end
        
        % Same as gf.preserveDataOnSpecifiedIndices(), but here we clear values instead of preserving
        function out = clearDataOnSpecifiedIndices(in,clearingIndices,clearingValue,dim)
            % -Permute inputs
            if(gf.isempty('dim')) dim = 1; end
            clearingIndices(clearingIndices==0) = nan;
            permOrder = 1:ndims(in); permOrder([1,dim]) = [dim,1];
            in = permute(in,permOrder);
            clearingIndices = permute(clearingIndices,permOrder);
            % -Find the gloabl clearing indices
            sizeI = size(in);
            modifiedIndices = reshape(0:prod(sizeI(2:end))-1,[1,sizeI(2:end)]);
            modifiedIndices = repmat(size(in,1)*modifiedIndices,size(clearingIndices,1),1) + clearingIndices;
            modifiedIndices = gf.noNaNs(modifiedIndices);
            % -Clear values of out in relevant indices
            out = in;
            out(modifiedIndices) = clearingValue;
            out = ipermute(out,permOrder);    
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -MATRIX STRUCTURE FUNCTIONS- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %% Populate all the fields in the structer into caller's workspace
        function struct2workspace(inStruct)
            if(isa(inStruct,'char'))
                inputname1 = inStruct;
            elseif(isa(inStruct,'struct'))
                inputname1 = inputname(1);
            else
                error('Error(struct2workspace): Input is not a structure!!');
            end
            fields = evalin('caller',['fieldnames(',inputname1,')']);
            for field = fields'
                field_ = field{1};
                evalin('caller',[field_ ' = ',inputname1,'.' field_ ';'])
            end
        end
        
        %% Convert a structure to string
        function str = struct2string(in,isTransposeColumnVector)
            if(gf.isempty('isTransposeColumnVector')) isTransposeColumnVector = false; end
            fieldNames = fieldnames(in);
            str = [inputname(1),'::',10];
            for fieldName = fieldNames'
                sizeV = gf.size(in.(fieldName{1}));
                if(isTransposeColumnVector)
                    if(numel(sizeV)==1)                                 %% Column vector
                        fieldContent = evalc(['disp(in.',fieldName{1},''' );']);
                    elseif(numel(sizeV)==2 && sizeV(1)==1)              %% Row vector
                        fieldContent = evalc(['disp(in.',fieldName{1},');']);
                    elseif(prod(sizeV)==0)                              %% Empty vector
                        fieldContent = ['',10];
                    end
                else
                    fieldContent = [10,evalc(['disp(in.',fieldName{1},');'])];
                end
                str = [str,'   ',fieldName{1},':  ',fieldContent];
            end
        end
        %% -Saves fields of a structure, instead of saving the structure itself!
        % EXAMPLE:
        %  s.a=12; s.b=34; save('file1.mat','s') => Save structure 's'
        %   "   " gf.saveStruct('file1.mat','s') => Save all fields in the structure 's' 
        function saveStruct(filename,in,appendOption)
            if(gf.isempty('appendOption') || ~strcmpi(appendOption,'-append'))
                appendOption = '';
            end
            fields = fieldnames(in);
            gf.struct2workspace(in);
            save(filename,fields{:},appendOption);            
        end
        
        %% -Add notes to an existing file, as an extra cell element
        %   'notes' must be a string
        function outStrings = addNotes(filename,notesIn)
            in = load(filename,'notes');	%% loading all elements in a .mat file will be slow for big .mat files
            if(~isa(notesIn,'cell')) notesIn = {notesIn}; end
            if(gf.isempty('in.notes')) 
                in.notes = {}; 
            elseif(~isa(in.notes,'cell'))
                in.notes = {in.notes};
            end
            notes = [in.notes;notesIn];
            %fields = fieldnames(in);
            %gf.struct2workspace(in);
            save(filename,'notes','-append');
            outStrings = in.notes;
        end
        

        %% 'in' must be either a filename or a structure
        %  'outStrings' is a column vector of cells, with each cell containing a string. First string will always be variable's name 
        %   Assumption: in.notes must be eitheir a string, or a column vector of cells of strings 
        function outStrings = grepNotes(in,keyword,inputname1)
            if(gf.isempty('keyword'))  keyword = 'notes'; end
            if(gf.isempty('inputname1'))  inputname1 = inputname(1); end
            if(isa(in,'char'))
                in = load(in);
                outStrings = gf.grepNotes(in,keyword);
            else
                outStrings = {};
                % -If in has got notes, then store it to outStrings 
                if(any(ismember(fieldnames(in)',keyword)))
                    outStrings = [outStrings;{[inputname1,'.',keyword]}];
                    if(isa(in.(keyword),'cell'))
                        outStrings = [outStrings;in.(keyword)(:)];
                    elseif(isa(in.(keyword),'char'))
                        outStrings = [outStrings;{in.(keyword)}];  
                    else
                       error(['Error: ',inputname1,'.',keyword,' is neither a cell nor a string!!!']); 
                    end
                end
                % -Look inside all other fields, and see if 'notes' exists or not, 
                for field = setdiff(fieldnames(in),keyword)'
                    if(~isempty(field))             %% Since even (0x1) cell cause entrance to this 'for' loop
                        field_ = field{1};
                        if(isa(in.(field_),'struct'))
                            subOut = gf.grepNotes(in.(field_),keyword,[inputname1,'.',field_]);
                            outStrings = [outStrings;subOut]; 
                        end
                    end
                end
            end
        end

        %% -Function to convert all row vectors to column vectors..
        % -NOTES:
        %      1. #Outputs = #Inputs
        %      2. If input is not a vector, then don't touch it.
        %      3. Inputting multi-dimensional matrix would be resulted in throwing error
        % -LOG:
        %      1. (2017/02/27): Splitted the existing function gf.convertRowVectorsToColumnVectors() into 2 function:
        %           gf.convertRowVectorsToColumnVectorsRecursively() and gf.convertRowVectorsToColumnVectors()
        function varargout = convertRowVectorsToColumnVectors(varargin)
            %% -For multiple inputs, call the function recursively for each input
            if(nargin>1)
                for iInput=1:nargin
                    varargout{iInput} = gf.convertRowVectorsToColumnVectors(varargin{iInput});
                end
            else
                varargout{1} = varargin{1};
                if(isa(varargout{1},'cell') || isa(varargout{1},'double') || isa(varargout{1},'logical'))
                    if(size(varargout{1},1)==1 && size(varargout{1},2)>1) 
                        varargout{1} = varargout{1}';
                    end
                end
            end
        end
        
        function varargout = convertColumnVectorsToRowVectors(varargin)
            %% -For multiple inputs, call the function recursively for each input
            if(nargin>1)
                for iInput=1:nargin
                    varargout{iInput} = gf.convertColumnVectorsToRowVectors(varargin{iInput});
                end
            else
                varargout{1} = varargin{1};
                if(isa(varargout{1},'cell') || isa(varargout{1},'double') || isa(varargout{1},'logical'))
                    if(size(varargout{1},2)==1 && size(varargout{1},1)>1) 
                        varargout{1} = varargout{1}';
                    end
                end
            end
        end
        
        
        %%%%%% -Function to convert all row vectors to column vectors..
        %%%%%% -Features
        %%%%%%      1. #Outputs = #Inputs
        %%%%%%      2. If input is cell/structure, go into the individual element, then do the job.
        %%%%%% -TODO
        %%%%%%      1. Has to come up with handling inputted class        
        function varargout = convertRowVectorsToColumnVectorsRecursively(varargin)
            %% -For multiple inputs, call the function recursively for each input
            if(nargin>1)
                for iInput=1:nargin
                    varargout{iInput} = gf.convertRowVectorsToColumnVectorsRecursively(varargin{iInput});
                end
            else
                if(isa(varargin{1},'cell'))
                    varargout{1} = varargin{1};
                    for iCell=1:numel(varargout{1})
                        varargout{1}{iCell} = gf.convertRowVectorsToColumnVectorsRecursively(varargin{1}{iCell});
                    end
                elseif(isa(varargin{1},'struct'))
                    fieldsInInputStructure = fieldnames(varargin{1});
                    varargout{1} = varargin{1};
                    for iField=1:length(fieldsInInputStructure)     
                        varargout{1}.(fieldsInInputStructure{iField}) = gf.convertRowVectorsToColumnVectorsRecursively(varargin{1}.(fieldsInInputStructure{iField}));
                    end
                elseif(isa(varargin{1},'double') || isa(varargin{1},'logical'))
                    if(size(varargin{1},1)==1) 
                        varargout{1} = varargin{1}';
                    else
                        varargout{1} = varargin{1};
                    end
                elseif(isa(varargin{1},'char'))
                    varargout{1} = varargin{1};
                else
                    varargout{1} = varargin{1};
                    % Error(['Unknown input type varargin = ',varargin]);   %% Commented this on 2017/02/24
                end
            end
        end
        
        function varargout = convertColumnVectorsToRowVectorsRecursively(varargin)
            %% -For multiple inputs, call the function recursively for each input
            if(nargin>1)
                for iInput=1:nargin
                    varargout{iInput} = gf.convertColumnVectorsToRowVectorsRecursively(varargin{iInput});
                end
            else
                if(isa(varargin{1},'cell'))
                    varargout{1} = varargin{1};
                    for iCell=1:numel(varargout{1})
                        varargout{1}{iCell} = gf.convertColumnVectorsToRowVectorsRecursively(varargin{1}{iCell});
                    end
                elseif(isa(varargin{1},'struct'))
                    fieldsInInputStructure = fieldnames(varargin{1});
                    varargout{1} = varargin{1};
                    for iField=1:length(fieldsInInputStructure)     
                        varargout{1}.(fieldsInInputStructure{iField}) = gf.convertColumnVectorsToRowVectorsRecursively(varargin{1}.(fieldsInInputStructure{iField}));
                    end
                elseif(isa(varargin{1},'double') || isa(varargin{1},'logical'))
                    if(size(varargin{1},2)==1) 
                        varargout{1} = varargin{1}';
                    else
                        varargout{1} = varargin{1};
                    end
                else
                    varargout{1} = varargin{1};
                    % Error(['Unknown input type varargin = ',varargin]);   %% Commented this on 2017/02/24
                end
            end
        end
        
        %% -Function to transpose all member of inputted structure
        function input = transposeStructure(input)
            if(~isa(input,'struct'))    %% -If input is not a structure, then transpose & terminate
                if(isa(input,'cell') || isa(input,'double') || isa(input,'logical'))
                    input = input';
                end
                return; 
            else                        %% -If input is a structure, then recursively call this function for each member of the structure,
                fieldsInInput = fieldnames(input);
                for iField=1:numel(fieldsInInput)
                    input.(fieldsInInput{iField}) = gf.transposeStructure(input.(fieldsInInput{iField}));
                end
            end
        end

        %% -Generic function to merge 2 structures, with priority given to inputStructure
        %- Detail: Fill all the empty structure in inputStructure by correspodig values from defaultStructure...
        %- If no check is needed for certain fields, then make corresponding element in possibleValuesPerField empty cell (ie {})
        % if 'isAllowNewMembers'=true, then 'inputStructure' can contain fields not present in 'defaultStructure'  
        function outputStructure = mergeTwoStructures(defaultStructure,inputStructure,possibleValuesPerField,isAllowNewMembers)
            % If 'defaultStructure'/'inputStructure' is a string, then bring the corresponding variable from caller's scope  
            if(isa(defaultStructure,'char')) 
                if(evalin('caller',['~gf.isempty(''',defaultStructure,''')'])) 
                    defaultStructure = evalin('caller',defaultStructure); 
                else
                    defaultStructure = struct(); 
                end;
            end;
            if(isa(inputStructure,'char')) 
                if(evalin('caller',['~gf.isempty(''',inputStructure,''')'])) 
                    inputStructure = evalin('caller',inputStructure); 
                else
                    inputStructure = struct(); 
                end;
            end;
            %%%%%%%%%%%% -If input structure is empty, then return defaultStructure- %%%%%%%% 
            if(gf.isempty('inputStructure')) outputStructure = defaultStructure; return; end
            if(gf.isempty('isAllowNewMembers')) isAllowNewMembers = false; end
            %%%%%%%%%%%% -Input Verification (Check if inputStructure contains fields which is not present in defaultStructure)- %%%%%%%%%%%%%%% 
            fieldsInDefaultStructure = fieldnames(defaultStructure);
            fieldsInInputStructure = fieldnames(inputStructure);
            if(~isAllowNewMembers)
                violatedFields = setdiff(fieldsInInputStructure,fieldsInDefaultStructure);
                if(~isempty(violatedFields)) 
                    gf.error('mergeTwoStructures: inputStructure contains fields, which is not present in defaultStructure as follows, ',violatedFields{:}); 
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%% -Fill the empty fields of inputStructure with that of defaultStructure- %%%%%%%%%%%%%% 
%             fieldsInDefaultStructure = fieldnames(defaultStructure);
            for iField=1:length(fieldsInDefaultStructure)
                if(~any(ismember(fieldsInInputStructure,fieldsInDefaultStructure{iField})) || isempty(inputStructure.(fieldsInDefaultStructure{iField})))
                    inputStructure.(fieldsInDefaultStructure{iField}) = defaultStructure.(fieldsInDefaultStructure{iField});
                elseif(isa(defaultStructure.(fieldsInDefaultStructure{iField}),'struct'))                                  %% -If field is present & it's a structure
                    inputStructure.(fieldsInDefaultStructure{iField}) = gf.mergeTwoStructures(defaultStructure.(fieldsInDefaultStructure{iField}),inputStructure.(fieldsInDefaultStructure{iField}),'',isAllowNewMembers);                    
                elseif(~gf.isempty('possibleValuesPerField'))    %% -Validate Non-empty fields  (TODO_Future: Extend this to include substructures as well)
                    if(~isempty(possibleValuesPerField{iField}) && ~any(ismember(possibleValuesPerField{iField},inputStructure.(fieldsInDefaultStructure{iField})))) %% -TODO_Future: Flip the inputs to ismember, and remove 'any' 
                        error(['Error(gf.mergeTwoStructures): invalid field in input structure,, ',fieldsInDefaultStructure{iField},' = ',inputStructure.(fieldsInDefaultStructure{iField})]);
                    end
                end
            end            
            outputStructure = inputStructure;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        %%   ****THIS FUNCTION IS OBOLETE***, Use function gf.mergeTwoStructures() instead. (Not deleting this, since many old matlab files are still using it)         
        %%%%%%%%%%%%%- Function to create Default Structure, from input Structure & defaultStructure %%%%%
        %%%%%%%%%%%%%- Detail: Fill all the emoty structure in inputStructure by correspodig values from defaultStructure...
        %%%%%%%%%%%%%- If no check is needed for certain fields, then make corresponding element in possibleValuesPerField empty cell (ie {})
        function outputStructure = createDefaultStructure(inputStructure,defaultStructure,possibleValuesPerField)
            %%%%%%%%%%%% -If input structure is empty, then return defaultStructure- %%%%%%%% 
            if(gf.isNullString(inputStructure)) outputStructure = defaultStructure; return; end
            
            %%%%%%%%%%%% -Input Verification- %%%%%%%%%%%%%%% 
            if(nargin>3) gf.error('createDefaultStructure: Number of Inputs must be 2. Currently nargin = ',nargin); end
            fieldsInInputStructure = fieldnames(inputStructure);
            violatedFields = {};
            for iField=1:length(fieldsInInputStructure)
                if(~isfield(defaultStructure,fieldsInInputStructure{iField}))
                    violatedFields = {violatedFields,fieldsInInputStructure{iField}};
                end
            end
            if(length(violatedFields)) gf.error('createDefaultStructure: inputStructure contains fields, which is not present in defaultStructure as follows, ',violatedFields{:}); end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%% -Fill the empty fields of inputStructure with that of defaultStructure- %%%%%%%%%%%%%% 
            fieldsInDefaultStructure = fieldnames(defaultStructure);
            for iField=1:length(fieldsInDefaultStructure)
                if(~isfield(inputStructure,fieldsInDefaultStructure{iField}) || gf.isNullString(inputStructure.(fieldsInDefaultStructure{iField}))) %% -If field is not present
                    inputStructure.(fieldsInDefaultStructure{iField}) = defaultStructure.(fieldsInDefaultStructure{iField});
                elseif(isa(defaultStructure.(fieldsInDefaultStructure{iField}),'struct'))                                  %% -If field is present & it's a structure
                    inputStructure.(fieldsInDefaultStructure{iField}) = gf.createDefaultStructure(defaultStructure.(fieldsInDefaultStructure{iField}),inputStructure.(fieldsInDefaultStructure{iField}));
                elseif(nargin>2)    %% -Validate Non-empty fields
                    if(~isempty(possibleValuesPerField{iField}) && ~any(ismember(possibleValuesPerField{iField},inputStructure.(fieldsInDefaultStructure{iField}))))
                        error(['Error(gf.createDefaultStructure): invalid field in input structure,, ',fieldsInDefaultStructure{iField},' = ',inputStructure.(fieldsInDefaultStructure{iField})]);
                    end
                end
            end            
            outputStructure = inputStructure;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        %% - Convert {struct{'a',12,'b',34},struct{'a',567,'b',890}}  ->  struct('a',{12,567},'b',{34,890})
        function out = convertCellsOfStructToStructOfCells(in)
            out = struct();
            for iCell = 1:numel(in)
                for field = fieldnames(in{iCell})'
                    field = field{1};       %% -Removing the 'cell' outer layer
                    if(isfield(out,field))
                        out.(field) = {out.(field){:},in{iCell}.(field)};   %% -NOTE: Horizontal concatenation chosen, since vertical concatenation is not supported in MATLAB!
                    else
                        out.(field) = {in{iCell}.(field)};
                    end
                end
            end
        end
        
        %% - Combine structures by giving last inputs as highest priority (ie 2nd input will overwrite common fields in 1st input)
        %   NOTE: Difference w.r.t function 'gf.mergeTwoStructures' is that, this function don't require all the fields to be present in 1st input!
        function out = mergeStructures(struct1,struct2,varargin)
           if(numel(varargin)>0)        %% -Recursion for multiple inputs
               tempOut = gf.mergeStructures(struct1,struct2);
               out = gf.mergeStructures(tempOut,varargin{:});
           else                         %% -For 2 inputs case
               out = struct1;
               for field = fieldnames(struct2)'
                  field = field{1};
                  if(~isfield(out,field) || isempty(out.(field)))
                     out.(field) = struct2.(field); 
                  end
               end
           end
        end
        
        %% -Function to vertcat 2 structures
        %   Same as MATLAB vertcat, except that this vercat function can concatenate 2 dissimilar structures! %% TODO_Refactoring (2017 Feb 27)
        function out = vertcat(in1,in2)
            assert(isa(in1,'struct') && isa(in2,'struct'),'Error(vertcat): Inputs must be structures');
            in1 = in1(:);
            in2 = in2(:);
            out = in1;
            if(isempty(fieldnames(in1))) in1Length = 0; else in1Length = numel(in1); end
            for i = 1:numel(in2)
                names = fieldnames(in2);
                for inames = 1:numel(names)
                   eval(['out(in1Length+i).',names{inames},' = in2(i).',names{inames},';']);
                end
            end
        end
        
        %% 'a{:}.b{:}.c{:}'  ->   cell(#a,#b,#c)
        % NOTE: 
        %   1. Here number of dimensions = number of ':' appearing in structExprStr. ie. each output-dimension indicate the corresponding ':' 
        %   2. Single ':' in structExprStr implies output is a column vector of cells  
        % EXAMPLE:
        %   gf.fetchFromStruct('out1.simOut(:).rowOut(:).simulationTime')
        %   gf.fetchFromStruct('arc{:}{:}')
        function out = fetchFromStruct(structExprStr)
            colonLocations = strfind(structExprStr,':');
            structExprStrFirstElement = structExprStr;
            structExprStrFirstElement(colonLocations) = '1';
            numelPerDim=[];
            for iDim = 1:numel(colonLocations)
                numelPerDim(iDim) = evalin('caller', ['numel(',structExprStrFirstElement(1:colonLocations(iDim)-2),')']);
            end
            out = cell([numelPerDim,1]);      %% '1' is appended to avoid making 2D matrix, while 'numelPerDim' is single element vector
            indexPerDimension = cell(size(numelPerDim));
            for iElement = 1:prod(numelPerDim)
                [indexPerDimension{:}] = ind2sub(numelPerDim,iElement);
                strArray = strsplit(structExprStr,':');
                strArrayCat = [strArray(1:end-1)',gf.num2str(indexPerDimension')];  %% Not using matlab num2str(), since given a column vector, matlab-num2str() prepends empty spaces before every number to have it proper rectangular shape.
                strArrayCat = strArrayCat';
                out{iElement} = evalin('caller',[strArrayCat{:},strArray{end}]);               
            end
        end
        
        %%
        % NOTE:
        %   1. One way to check if this function works or not, is to execute your input by replacing all variable indices with '1's 
        %       i.e.  s{1}.s1(1).s2{1}.a
        %   2. 'isFillWithNan' = true (Default) means, absent fields would be regarded as NaN 
        %   3. Here number of dimensions = number of brackets appearing in 'inString' + dimension of last field  
        %   4. Single ':' in structExprStr implies output is a column vector of cells  
        %   5. s{[1,2;3,4],[1,2]}.a will result in 1D output
        % PROGRAMMING NOTES:
        %   1. There are some design errors in MATLAB software, as follows
        %       a). MATLAB ignores, tralining singleton dimensions. 
        %       b). 'a(2:5)' and 'a{:}' are outputted as row-vectors, instead of column vecotors
        %       c). ndims of column vector is 2 (instead of 1)!!. 
        %     Solution:
        %       The problem (a) is rectified using 'nDimPerScope'.  To overcome (c), I implemented gf.ndims().
        % TODO:
        %   1. Include accepting 'end' keyword.
        %   2. Nonworking Input: gf.congregate('in{:,2,3}{1}.F')
        function out = congregate(inString,isFillWithNan,squeezeLevel)
           if(gf.isempty('isFillWithNan')) isFillWithNan = true; end
           if(gf.isempty('squeezeLevel')) squeezeLevel = 0; end
           %% -USER calling.
           %    i.e., Just bring the USER's indexing variables to this scope, rename it (also rename inString), and pass it to next scope  
               inString = strtrim(inString);
               %% Fetch variables from USER scope
               [variableNames,remainder] = regexp(inString,'(?<=^|\(|\[|\{|,|:)[a-zA-Z]\w*','match','split'); %% variableNames{1} is base object
               remainder(1)=[];
               for i=1:numel(variableNames)
                  userVariables{i} = evalin('caller',variableNames{i}); 
                  userVariablesNewname{i} = ['userVariables{',num2str(i),'}'];
               end
               inStringNew = gf.strjoin(userVariablesNewname,remainder);
               [out,nDimPerScope] = gf.congregate_Recursive(inStringNew,isFillWithNan);
               out = gf.cellOfCell2mat(out,squeezeLevel,nDimPerScope);
        end
        
        %% NOTE: Calling of this function is restricted to only gf.congregate()
        %   'nDimsCurrent' : Local variable indicating whether current scope deserve how many dimensions
        %            NaN   => Discard this value in gf.cellOfCell2mat   
        %            0     => 0 dimension in gf.cellOfCell2mat   
        function [out,nDimPerScope] = congregate_Recursive(inString,isFillWithNan,recursiveCallingScope) 
               if(nargin<3) recursiveCallingScope = 1; end                      %% Avoided using gf.isempty() for speed
               userVariables = evalin('caller','userVariables');
               % Get 's(1).s1(p)' from 's(1).s1(p).s2(q).a'
               [from, to] = gf.findComplicatedParentheses(inString,1);
               %  TODO_Verify(2016/11/3):   Removed && to~=numel(inString) from below expression
               if(~isempty(to) && ~gf.isempty(inString(1:to)))  %% TODO_Verify(2016/11/3): Changed from gf.isempty(['[',inString(1:to),']']) to gf.isempty(['inString(1:to)'])
                   inStringnew = inString;
                   inStringnew(from) = '(';  inStringnew(to) = ')';
                   tmpIn = eval(inStringnew(1:to));                             %% Changing the base object
                   [userVariables{1},nDimsCurrent] = congregate_reshape(tmpIn); %% Changing vector 'tempIn' to intended matrix 'userVariables{1}'
                   out = cell(size(userVariables{1}));
                   for i=1:numel(out)
                       if(isa(userVariables{1},'cell'))
                        inStringNew = ['userVariables{1}{',num2str(i),'}',inString(to+1:end)];
                       else
                        inStringNew = ['userVariables{1}(',num2str(i),')',inString(to+1:end)];
                       end
                       [out{i},nDimsSub{i}] = gf.congregate_Recursive(inStringNew,isFillWithNan,recursiveCallingScope+1);
                   end
                   nDimsSub = nanmax(gf.cellOfCell2mat(nDimsSub)',[],2);
                   nDimPerScope = [nDimsCurrent;nDimsSub];      %% Since sometimes 'inStringNew' would be empty cells
               else
                  if(recursiveCallingScope==1 && ~isempty(to) && gf.isempty(['[',inString(1:to),']']))  %% TODO_Verify(2016/11/3): Changed from gf.isempty(['[',inString(1:to),']']) to gf.isempty(inString(1:to))
                     error(['Error: Non-existent parent object  ']);
                  end
                  if(eval(['~gf.isempty(''',inString,''')']))
                     out = eval(inString);
                  else                  %% Missing or Empty field
                     if(isFillWithNan) 
                         out = NaN;
                     else
                         error(['Non-existant field: ',inString]);
                     end
                  end
                  nDimPerScope = gf.ndims(out);
               end
            %% This function is written in order to make sure that 'a{:,:,:}' is considered as a 3D matrix instead of doing vectorization! 
            %   'nDimsCurrent': 
            %       0   => tmpIn is a single object without any ':'
            %       NaN => Ignore this value by any clients
            %           for each ':', and '1:1' is taken as extra dimension
            % PROGRAMMING NOTE:
            %   1. MATLAB create 'a(:)' and 'a{:}' as a row-vector, but ideally it should be column vectors. So here we are always sticking with the ideal.  
            % TODO_Future: 
            %   Program to handle every complicated expressions too 
            function [tmpIn,nDimsCurrent] = congregate_reshape(tmpIn)
                indicesInputStr = regexp(strtrim(inString(from+1:to-1)),',','split');   %% Get entries inside last bracket 
                nDimsCurrent = 0;
%                 if(numel(indicesInputStr)>1)      %% TODO: Clarify: Why did I placed this constraint in the first place? 
                    clear finalSizePerDim;
                    for i=1:numel(indicesInputStr)
                        if(regexp(indicesInputStr{i},'^\s*:\s*$'))
                          if(numel(indicesInputStr)==1)     %% If 'a(:)'
                            finalSizePerDim(i) = numel(eval(inString(1:from-1)));
                          else
                            finalSizePerDim(i) = size(eval(inString(1:from-1)),i);
                          end
                          nDimsCurrent = nDimsCurrent + 1;
                        else
                          finalSizePerDim(i) = numel(eval(indicesInputStr{i}));
                          if(finalSizePerDim(i)>1 || ~isempty(regexp(indicesInputStr{i},':')))    
                            nDimsCurrent = nDimsCurrent + 1;
                          end
                        end
                    end
                    tmp1 = reshape(tmpIn,[finalSizePerDim,1]);     %% '1' is appended to avoid the error when numel(finalSizePerDim)==1   
                    tmpIn = tmp1;
%                 end
%                 nDimsCurrent = numel(finalSizePerDim);
            end            
        end
           

        %% Function to find '([1,2,3])' in 's(1).s1([1,2,3])' (Note from & to corresponds to bracket location)
        %       '(1,2)', '('abcd')' are not complicated parentheses.
        %   'count' determines how many outputs needed. If count>1, then an array of cells is outputted   
        function [from,to] = findComplicatedParentheses(str,count)
            if(nargin<2) count = Inf; end                                               %% Avoided using gf.isempty() for speed
            if(count>1)
                % Count the number of brackets in str
                i = 1; to_last = 0;
                while i<=count && to_last < numel(str)
                   [from{i},to{i}] = gf.findComplicatedParentheses(str(to_last+1:end),1);
                   if(isempty(from{i})) 
                       break;           
                   else
                       from{i} = from{i} + to_last;
                       to{i} = to{i} + to_last;
                       to_last = to{i};
                       i = i+1;
                   end
                end
            else
                indices = regexp(str,'{|(');
                for i=1:numel(indices)
                   from = indices(i);
                   to = gf.findMatchingParentheses(str,from);
                   if(isempty(regexp(str(from+1:to-1),'^((\s*\d+\s*)|(''.*''))$')))    %% If number or a string
                       return;
                   end
                end
                from = []; to = [];
            end
        end
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -CELL FUNCTIONS- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %% {{[1,2],[3,4,5]},{[11,12]'}}  ->  []  
        %   Brief Description:
        %       Each {} introduces an extra dimension (eventhough it's singleton dimension). To avoid this feature, input 'nDimPerScope' carefully. 
        %   -'squeezeLevel'
        %       0: No Squeezing. Expect redundant dimensions for row-vectors   
        %       1: 'in' (which is a matrix of cells) and all cells within is (in all scope) are squeezed 
        %       2: Both 'in' and all the numerical matrix content inside any cells are squeezed
        %   - 'nDimPerScope'
        %       This would make sure that, ith scope has got nDimPerScope(i) dimensions, by inserting singleton dimensions. (Highly useful in gf.congregate_Recursive) 
        %       nDimPerScope=nan  =>  Ignore this parameter
        function out = cellOfCell2mat(in,squeezeLevel,nDimPerScope)
            if(nargin < 2 || gf.isempty('squeezeLevel')) squeezeLevel = 2; end
            if(nargin < 3 || gf.isempty('nDimPerScope')) nDimPerScope = [NaN]; end
            if(numel(nDimPerScope)==1) nDimPerScope = [nDimPerScope;NaN]; end
            %
            if(isa(in,'cell'))
                outTemp = cell(size(in));
                if(squeezeLevel>0)
                    outTemp = gf.squeezeAndTransposeToRowVector(outTemp);
                end
                offsetDim = gf.ndims(outTemp);
                for iCell = 1:numel(in)
                    outTemp{iCell} = gf.cellOfCell2mat(in{iCell},squeezeLevel,nDimPerScope(2:end));
                    outTemp{iCell} = shiftdim(outTemp{iCell},-offsetDim);
                end
                out = gf.cell2mat(outTemp);         %% TODO_Next: Write another program called gf.cell2matRaw() which wouldn't loose any dimensions. Also clearly define 'squeezeLevel'   
                if(~isnan(nDimPerScope(1)) && gf.ndims(outTemp)<gf.ndims(out))   %% Insert singleton dimensions based on 'nDimPerScope' 
                    if(nDimPerScope(1)>0 && gf.ndims(outTemp)>nDimPerScope(1))
                        error('Error(): Inputted nDimPerScope is not matching with');
                    end
                    if(gf.ndims(outTemp)<nDimPerScope(1))   %% Adding extra dimensions to the current scope
                        newDimIndices = [1:gf.ndims(outTemp), gf.ndims(out)+1:gf.ndims(out)+1+(nDimPerScope(1)-gf.ndims(outTemp)-1), ...
                             gf.ndims(outTemp)+1:gf.ndims(out)];
                        out = permute(out,newDimIndices);
                    elseif(nDimPerScope(1)==0)              %% Deleting current dimension, by appending it the subscope dimensions 
                        outSize = size(out);
                        outSize(gf.ndims(outTemp)+1) = prod(outSize(1:gf.ndims(outTemp)+1));
                        outSize(1:gf.ndims(outTemp)) = [];
                        out = reshape(out,[outSize,1]);
                    end
                end
            else
               if(squeezeLevel>1)
                out = gf.squeezeAndTransposeToRowVector(in);
               else
                out = in; 
               end
            end
        end
        
        %% Rename fields of a structure using regular expressions
        function out = renameFields(in,repString,subString)
            data = struct2cell(in);
            fields = fieldnames(in);
            fields = regexprep(fields,repString,subString);
            out = cell2struct(data, fields);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -STRING MANIPULATION- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        function out = convertToMatlabVarName(in,ignoreRules)
            if(gf.isempty('ignoreRules')) ignoreRules = false; end
            if(isa(in,'cell'))
                out = cell(size(in));
                for i=1:numel(in)
                   out{i} = gf.convertToMatlabVarName(in{i});
                end
            elseif(isa(in,'char'))
                out = regexprep(in,'\s+|\.|&|\?','_');
                if(~ignoreRules)
                    out = regexprep(out,'^(?=[0-9])','a');   %% Make variable name to start with a letter.
                end
            else
                error('Error: input must be either a string, or an array of strings');
            end
        end
        
        function out = isNullString(inputString)
            if(isa(inputString,'char'))
                out = ~isempty(regexp(inputString,'^\s*$','emptymatch'));
            else
                out=0;
            end
        end
        
        %%- Example getField('1.2.3',2) --> output=2;
        function outNumber = getField(inputString,fieldIndices)
            outNumber = nan(size(fieldIndices));
            for i=1:numel(fieldIndices)
                [st,fieldValue]=system(['echo ',inputString,' | perl -pe ''s/(\d+\.){',num2str(fieldIndices(i)-1),'}//'' | grep -Po ''\d+'' |head -1']);
                if(st) error(['Error(): Can''t grep integer fields']); end
                outNumber(i) = str2double(fieldValue);
            end
        end
        
        %% Function to compare different types of objects
        %   'indexPerin2' is an array of cells, with indexPerin2{i}=indices => in1{indices}==in2{i}  
        function indexPerin2 = contains(in1,in2,isCaseInsensitive,isRegex)
            if(nargin < 3 || gf.isempty('isCaseInsensitive')) isCaseInsensitive = true; end
            if(nargin < 4) isRegex = false; end
            if(~isa(in1,'cell')) in1 = {in1}; end
            if(~isa(in2,'cell')) in2 = {in2}; end
            indexPerin2 = cell(numel(in2),1);
            for i=1:numel(in2)
                if(isa(in2{i},'char'))
                    if(isRegex)
                        fh = @regexp;
                    elseif(isCaseInsensitive) 
                        fh=@strcmpi; 
                    else
                        fh=@strcmp; 
                    end
                else
                    fh=@isequal;
                end
                for j=1:numel(in1)
                    if(strcmp(class(in1{j}),class(in2{i})))
                        if(fh(in1{j},in2{i}))
                            indexPerin2{i} = [indexPerin2{i};j];
                        end
                    end
                end
            end
            if(numel(indexPerin2)==1 && isempty(indexPerin2{1}))    %% Make it fully empty output.
                indexPerin2 = cell(numel(in2),0);
            end
        end
        
        
        %% 'in' must be array of cells, with each cell containing a string
        % [an,ind]=gf.findregexp({'abcd','efab','jkl','ab12'},'ab.*')  ->  an = {'abcd','ab12'}, ind = [1,4]
        function [out,ind] = findregexp(in,expr)
            expr = ['^',expr,'$'];
            out = regexp(in,expr,'match','once');
            if(nargout>1)
                isEmpty = cellfun(@isempty,out);
                isEmpty = gf.convertRowVectorsToColumnVectors(isEmpty);
                ind = find(~isEmpty);
            end
            out = gf.squeeze(out);
        end
        
        function [out,ind] = findregexpNot(in,expr)
            [~,ind] = gf.findregexp(in,expr);
            ind = setdiff([1:numel(in)]',ind);
            out = in(ind);
        end
        
        
        %% -Same as MATLAB inbuilt function ismember, except that,
        %   1. All string inputs are considered case insensitive
        %   2. Strings and numbers also compared instead of throwing an error! (TODO: Implement this!)
        function varargout = ismemberi(A,B)
            A = gf.lower(A);
            B = gf.lower(B);
            if(nargout>1)
                [varargout{1},varargout{2}] = ismember(A,B);
            else
                varargout{1} = ismember(A,B);                
            end
        end
        
        %% -Same as MATLAB built in 'lower' function, except that it leaves everything except 'string', as such!!
        function out = lower(in)
            out = in;
            if(isa(in,'cell'))
               for i=1:numel(in)
                  out{i} = gf.lower(in{i}); 
               end
            elseif(isa(in,'char'))
                out = lower(in);
            end
        end

        %% Find matching parentheses
        %   findMatchingParentheses('a(bcd)',2)  ->  6  
        %   findMatchingParentheses('a(bcd)',6)  ->  1  
        % TODO: Publish it
        % https://www.mathworks.com/matlabcentral/cody/problems/111-find-matching-parenthesis
        function index = findMatchingParentheses(str,inIndex)
            parentheses = str(inIndex);
            switch(parentheses)
                case '{', antiParentheses = '}';
                case '(', antiParentheses = ')';
                case '[', antiParentheses = ']';
                case '}', antiParentheses = '{';
                case ')', antiParentheses = '(';
                case ']', antiParentheses = '[';
                otherwise, error('Index in the inputted str doesn''t corresponds to a bracket!!');
            end
            isParentheses = zeros(numel(str),1);
            ind1 = regexp(str,parentheses);
            ind2 = regexp(str,antiParentheses);
            isParentheses(ind1) = 1;
            isParentheses(ind2) = -1;
            if(regexp(parentheses,'{|(|['))
                index = find(cumsum(isParentheses(inIndex:end))==0,1,'first');
                index = inIndex-1 + index;
            else
                index = find(cumsum(isParentheses(1:inIndex),1,'reverse')==0,1,'last');
            end
        end
        
        %  gf.repairFunctionHandleString('alg_edgeRemoval('' '',[],3)','a','b')  ->    'alg_edgeRemoval(a,[],3)'  
        function outString = repairFunctionHandleString(fh,varargin)
            %[inVariables,remainders] = regexp(fh,'((?<=\().*?(?=\,))|((?<=,).*?(?=\,))|((?<=,).*?(?=\)))','match','split');
            if(isempty(regexp(fh,'\)\s*$')))        %% In case input is function name alone, without any arguments.
                fh = [fh,'(',repmat(' '''',',[1,numel(varargin)])];
                fh = [fh(1:end-1),')'];
            end
            [inVariables,remainders] = regexp(fh,'((?<=\().*?(?=\,))|((?<=,).*?(?=\,))|((?<=,).*?(?=\)))','match','split');
            for i = 1:numel(inVariables)
                if(regexp(strtrim(inVariables{i}),'^\''\s*\''$|^\[\s*\]$'))
                   inVariables{i} = varargin{i};
                end
            end
            outString = strjoin(remainders,inVariables);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -FILE HANDLING- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
        %% Same as MATLAB inbuilt 'dir' function but this function remove all hidden files (ie. files starting with '.') 
        function files = dir(foldername,filePattern)
           if(gf.isempty('filePattern')) filePattern = '.*'; end            
            tmp = dir(foldername);
            files = {tmp.name};
            ind1 = gf.contains(files,'^[^.].*',false,true);     %% Remove all filenames starting with '.'
            ind2 = gf.contains(files,filePattern,false,true);     %% Remove all filenames starting with '.'
            if(~isempty(ind1) && ~isempty(ind2))
                ind3 = intersect(ind1{1},ind2{1});
                files = files(ind3);
            end
        end


        %% Assigns data from '.csv'/'.xlsx' files
        % Wrapper function for 'gf.assignData_' function  
        function assignData(fileName,LHS,RHS,rowIndices,delimiter)
           if(gf.isempty('delimiter')) delimiter = ';'; end
           if(gf.isempty('RHS')) RHS = ''; end           
           LHS = strtrim(LHS);
           RHS = strtrim(RHS); 
           %% Get base-object name, and assign it an empty cell if it doesn't exist in caller's scope 
           baseObjName = regexp(LHS,'^.*?(?=\(|\{|\.|$)','match');
           baseObjName = baseObjName{1};
           if(evalin('caller',['gf.isempty(''',baseObjName,''')']))
            assignin('caller',baseObjName,{});
           end
           %% Bring variables used in LHS & RHS to userVariablesLHS{} & userVariablesRHS{} respectively  
           userVariablesLHS = {}; userVariablesRHS = {};
           [LHS,cmdLHS] = gf.bringVariablesFromCallersScope(LHS,'userVariablesLHS');    %% Assumption userVariablesLHS{1} would contain base object 
           [RHS,cmdRHS] = gf.bringVariablesFromCallersScope(RHS,'userVariablesRHS');
           eval(cmdLHS); eval(cmdRHS);
           %% -Get data from file
           if(exist(fileName,'file')) 
               extension = regexp(fileName,'\..*?(?=\s*)$','match');
               switch(extension{1})
                   case '.csv'
                    data = gf.csvread(fileName,delimiter);
                   case '.xlsx'
                    [~,~,data] = xlsread(fileName);
                   otherwise
                    error('Invalid fileName extension');
               end
           elseif(isa(fileName,'cell'))
               data = fileName;
           else
               error('First input must be either an existing fileName, or matrix of cells');
           end
           if(~gf.isempty('rowIndices'))    %% If 'rowIndices' specified, then trim off the remaining rows 
               if(isa(rowIndices,'double'))
                data = data(rowIndices,:);
               else
                data = eval(['data(',rowIndices,',:)']); 
               end
           end  
           % Fetch values from fileName
           if(~gf.isempty('RHS'))
            gf.assignData_(data,LHS,RHS);
           else
            userVariablesLHS{1} = data;
           end
           % Assign object in the caller's scope
           assignin('caller',baseObjName,userVariablesLHS{1});            
        end
        
        % Both LHS and RHS are strings
        %   'data' must be matrix of cells
        % Handle the cases of ('s(<A>).s1(<B>).<field1*>','<field1*>'),         -> Done 
        %                     ('s(<A+10>).s1(<B>).<field1*>','<field1*>'),
        %                     ('s(<A>).s1(<B>).<1:10>','<1:10>'),               -> Done 
        %                     ('s(<A>).s1(<B>).s2(<1:10>)','<1:10>'),           -> Done 
        % Details:
        %   LHSVarIndexNames: Indicate the header of columns, which should be used for corresponding index {<.*>} 
        %   RHS:              Indicate how many columns must be fetched
        %   LHSFieldIndexNames: Indicate how many columns must be fetched as d
        % Warning:
        %   1. XL-sheet header names must always start with a letter (not a number).
        %   2. Remove all unwanted space in the input strings. (Replace with comma for vectors ) 
        % NOTE:
        %   1. This function can't be called independently, it can be called only from gf.assignData
        % Assumption:
        %   1. 'LHS' must always start with 'userVariablesLHS{1}'
        % TODO:
        %   1. Handle the variable indexing (ie p1 & p2 in 'patient{<User>}.session{<Session>}.affectedArm.pointTest.<p1:p2>' ) 
        function assignData_(data,LHS,RHS)
            userVariablesLHS = evalin('caller','userVariablesLHS');
            userVariablesRHS = evalin('caller','userVariablesRHS');
            header = data(1,:);
%             % Bring baseObj into current scope from caller's scope 
%             baseObjName = regexp(LHS,'^.*?(?=\(|\{|\.|$)','match');
%             baseObjName = baseObjName{1};
%             tmp = evalin('caller',baseObjName);
%             eval([baseObjName,'= tmp;']);
            % -Split LHS
            [LHSVarIndexNames,from,to] = regexp(LHS,'(?<=(\(|\{)<\s*)[a-zA-Z].*?(?=\s*>(\}|\)))','match');
            %% Handle LHSVarColumnNames here,  {i.e., handle '<User>' part in 's.s1(<User>).s2(<1:10>)'} 
            if(~isempty(LHSVarIndexNames))    
                columnIndexPerLHSVarIndexNames = gf.contains(header,LHSVarIndexNames,false,true);
                if(any(cellfun(@isempty,columnIndexPerLHSVarIndexNames)))
                    missingInd =  find(cellfun(@isempty,columnIndexPerLHSVarIndexNames));
                    error(['Following titles are missing from the column headers!!  ',strjoin(LHSVarIndexNames(missingInd),',  ')]);
                end
                columnIndexPerLHSVarIndexNames = cell2mat(columnIndexPerLHSVarIndexNames);
                LHSSegmentPrime = {};
                for i=1:numel(LHSVarIndexNames)
                    if(i==1)
                       LHSSegmentPrime = {LHSSegmentPrime{:},LHS(1:from(i)-2)};
                    else
                       LHSSegmentPrime = {LHSSegmentPrime{:},LHS(to(i-1)+2:from(i)-2)};                
                    end
                end
                LHSSegmentPrime = {LHSSegmentPrime{:},LHS(to(end)+2:end)};
                while size(data,1)>1
                    valuesPerLHSVarIndexNames_string = gf.num2str(data(2,columnIndexPerLHSVarIndexNames));      %% NOTE: for csv file, data{i,j} is always a string; whereas for xlsx files, data{i,j} can be a number too 
                    valuesPerLHSVarIndexNames = cell2mat(gf.str2num(valuesPerLHSVarIndexNames_string));
                    rowIndicesToFetch = find(ismember(cell2mat(data(2:end,columnIndexPerLHSVarIndexNames)),valuesPerLHSVarIndexNames,'rows')) +1; %% Find row indices which have got same 'valuesPerLHSVarIndexNames' 
                    LHSnew = gf.strjoin(LHSSegmentPrime,valuesPerLHSVarIndexNames_string);
                    gf.assignData_(data([1;rowIndicesToFetch],:),LHSnew,RHS);
                    data(rowIndicesToFetch,:) = [];
                end
            %% Handle LHSVarFieldNames here,  {i.e., handle '<1:10>' part in 's.s1(1).s2(<1:10>)}     
            else                                                
                % Dereferencing in RHS
                RHSVarFieldNames = regexp(RHS,'(?<=<).*?(?=>)','match');
                if(isempty(RHSVarFieldNames)) RHSVarFieldNames = RHS; else RHSVarFieldNames = RHSVarFieldNames{1}; end
                if(regexp(RHSVarFieldNames,'[a-zA-Z]'))
                    RHSColIndices = cell2mat(gf.contains(header,RHSVarFieldNames,false,true));
                elseif(regexp(RHSVarFieldNames,'[0-9]'))        %% For the case LHSVarFieldNames='1:10'
                    RHSColIndices = eval(RHSVarFieldNames);
                else
                    error('Specify proper RHS Filds or columns indices');
                end
                %% Dereferencing LHS
                [LHSVarFieldNames,from,to] = regexp(LHS,'(?<=<).*?(?=>)','match');
                if(numel(LHSVarFieldNames)>1)
                    error('Error: Multiple number of variable fieldnames in LHS');
                elseif(numel(LHSVarFieldNames)==1)
                    % Dereference all the fieldnames/numerical indices in LHS 
                    if(regexp(LHSVarFieldNames{1},'[a-zA-Z]'))  %% For the case s.s1(1).<.*Target.*> (Fetch all columns with header with pattern '.*Target.*')
                        tmpIndices = cell2mat(gf.contains(header,LHSVarFieldNames{1},false,true));
                        LHSFieldNames = header(1,tmpIndices);
                    elseif(regexp(LHS(from-2),'\(|\{'))         %% For the case s.s1(1).s2(<1:10>)  (i.e., when LHSVarFieldNames='1:10') 
                        LHSFieldNames = gf.num2str(num2cell(eval(LHSVarFieldNames{1})));
                    else                                        %% For the case s.s1(1).<1:10> (i.e., when LHSVarFieldNames='1:10')
                        tmpIndices = eval(LHSVarFieldNames{1});
                        LHSFieldNames = header(1,tmpIndices);
                    end
                    LHSFieldNames = gf.convertToMatlabVarName(LHSFieldNames,true);
                end
                %% Assign LHS with RHS
                if(gf.isempty('LHSFieldNames'))
                   try
                    eval([LHS,' = cell2mat(data(2:end,RHSColIndices));']); 
                   catch
                    eval([LHS,' = data(2:end,RHSColIndices);']); 
                   end
                elseif(numel(LHSFieldNames)==numel(RHSColIndices))
                   LHS(end+1) = ' ';
                   for iCol=1:numel(RHSColIndices)
                       if(to(1)+2 <= numel(LHS))
                        LHSnew = [LHS(1:from(1)-2),LHSFieldNames{iCol},LHS(to(1)+2:end)];
                       else
                        LHSnew = [LHS(1:from(1)-2),LHSFieldNames{iCol}];                   
                       end
                       currentData = data(2:end,RHSColIndices(iCol));
                       if(isa(currentData,'cell'))  
                           if(isa(currentData{1},'double'))
                               eval([LHSnew,'= cell2mat(currentData);']);
                           else
                               eval([LHSnew,'= currentData;']);                               
                           end
                       else
                           error('Input data must be a cell array!!');
                       end                     
                   end
                else
                    error('Error: Number of field names in LHS is not matching with number of field names in RHS');
                end
            end
%             assignin('caller',baseObjName,eval(baseObjName));
            assignin('caller','tmp_for_assignData_function',userVariablesLHS{1});   
            evalin('caller','userVariablesLHS{1} = tmp_for_assignData_function; clear(''tmp_for_assignData_function'');');   
        end
        
        %% When a function tries to execute a commandString, it has to bring all the necessary variables from caller's scope apriori.
        %      This function helps to do that.
        % Inputs:
        %   'callerVariableNewName': is a string specifying, to which array of cells we need to bring the values from caller's scope  
        function [outString,cmdToExecute] = bringVariablesFromCallersScope(inString,callerVariableNewName)
            if(gf.isempty('callerVariableNewName')) callerVariableNewName = 'userVariables'; end
            %
            [variableNames,remainder] = regexp(inString,'(?<=^|\(|\[|\{|,|:)[a-zA-Z]\w*','match','split');
            cmdToExecute = '';
            %
            userVariablesNewname = {};
            for i=1:numel(variableNames)
              cmdToExecute = [cmdToExecute,callerVariableNewName,'{',num2str(i),'} = evalin(''caller'',''',variableNames{i},''');']; 
              userVariablesNewname{i} = [callerVariableNewName,'{',num2str(i),'}'];
            end
            %
            outString = gf.strjoin(remainder,userVariablesNewname);
        end     
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Indexing Functions- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Currently only 2D-Matrix is supported.. Has to extend it in future.. %%%%%%%%%%%%  
        function out = convertLogicalIndexToPhysicalIndex(logicalIndexForIndices,options)
            %%%%% Create Proper/Default Options %%%%%%%
            if(nargin>1 && isfield(options,'preserveShape') && strcmp(options.preserveShape,'false'))
                options.preserveShape = false;
            else
                options.preserveShape = true;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            inputMatrixSize = size(logicalIndexForIndices);
            inputIndices = reshape(1:prod(inputMatrixSize),inputMatrixSize);
            if(options.preserveShape)   %%% output must be same shape as input...
                out = nan(inputMatrixSize);    
                for iRow = 1:inputMatrixSize(1)
                    temp = inputIndices(iRow,logicalIndexForIndices(iRow,:));
                    out(iRow,1:length(temp)) = temp;
                end
            else                %%% output is just a vector...
                out = inputIndices(logicalIndexForIndices);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Same as MATLAB's sub2ind, except that, 
        %   1. Here we ignore the '0' indices. (Used when we use matrixGeneralized)
        %   2. If an index is singleton, then we extend it by repeatation.    
        %  Example:
        %   1. a(gf.sub2ind(size(a),:,:)) is diagonal vector
        %   2. a(gf.sub2ind(size(a),:,2)) is 2nd column vector of 'a'
        function outIndices = sub2ind(inSize,varargin)
            % -First, get rid of ':' in varargin
            for iInput = 1:numel(varargin)
               if(strcmp(varargin{iInput},':'))
                   varargin{iInput} = [1:inSize(iInput)]';
               else
                   varargin{iInput} = varargin{iInput}(:); %% -Converting to column matrix.                   
               end
            end  
            % -Second, make singleton elements repeating, thus making all input indices to same size
            numelPerVarargin = gf.numel(varargin{:});
            if(all(ismember(numelPerVarargin,[0,1]')))      %% -Exception, in case all non-fixed indices are of zero length
                outIndices = [];
                return; 
            end
            idealNumIndices = max(numelPerVarargin);
            assert(all(numelPerVarargin==1 | numelPerVarargin==idealNumIndices),'Number of indices in subIndex must be same (or singleton)');
            for iInput = find(numelPerVarargin==1)'
                varargin{iInput} = repmat(varargin{iInput},idealNumIndices,1);
            end            
            % -Third, remove 0 indices
            for i=1:numel(varargin)
                for j=setdiff(1:numel(varargin),i)
                    if(numel(varargin{i}) == numel(varargin{i}))
                        varargin{j}(varargin{i}==0) = [];                        
                    end
                end
                varargin{i}(varargin{i}==0) = [];
            end
            outIndices = sub2ind(inSize,varargin{:});
        end
        
        %% -Following is a modified version of 'sub2ind' function, in which 
        %    1. ':' can also be inputted, but unlike gf.sub2ind(), here it means to take all the possible combination of values from that dimension
        %       Example: a=magic(4); a(gf.mysub2ind(size(a),:,[2,3])) outputs 8x1 vector
        %       -Note: (Handling ':'), First all dimensions corresponding to ':' are considered as of length 1.
        %       Absolute indices are computed, then extended to the dimensions of ':'...
        %   TODO:
        %    2. TODO_Future: Create a new advanced function, which can specify which all dimensions are filtering dimensions (as in matlab sub2ind), and which all are the re-ordering dimensions
        function outIndices = mysub2ind(sizeOfInputMatrix,varargin) 
            outIndices = [1];
            if(numel(varargin)>numel(sizeOfInputMatrix))        %% Append trailing dimensions
                sizeOfInputMatrix(end+1:numel(varargin)) = 1;
            end
            sizeOfInputMatrix_cumProd = [1,cumprod(sizeOfInputMatrix)];
            for iInput=1:numel(varargin)
               if(~strcmp(varargin{iInput},':'))
                   varargin{iInput} = varargin{iInput}(:);                              %% -Converting to column matrix.
                   if(numel(outIndices)>1 && numel(outIndices)~=numel(varargin{iInput})) 
                       error(['Error(mysub2ind): number of indices in subIndex not matching to the previous size']);
                   end
                   outIndices = outIndices+(varargin{iInput}-1)*sizeOfInputMatrix_cumProd(iInput);
               end
                outIndices=outIndices(:);
            end
            for iInput=1:numel(varargin)
               if(strcmp(varargin{iInput},':'))
                   outIndices = outIndices*ones(1,sizeOfInputMatrix(iInput)) + ones(size(outIndices))*[0:sizeOfInputMatrix(iInput)-1]*sizeOfInputMatrix_cumProd(iInput);
               end
                outIndices=outIndices(:);
            end
            outIndices = sort(outIndices);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% -Function to Swap dimensions- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Following function has become obsolete, since matlab built-in function permute does the same job!!!
        function outputMatrix = swapDimension(inputMatrix,dimension1,dimension2)
            if(~isa(inputMatrix,'double'))
                error(['Error: input must be a Matrix!!! ']);
            end
            Dmin=min(dimension1,dimension2); Dmax=max(dimension1,dimension2);
            inputSize = size(inputMatrix); outputSize = inputSize;
            outputSize(Dmin)=inputSize(Dmax); outputSize(Dmax)=inputSize(Dmin);
            outputMatrix = nan(outputSize);
            
            leftColons=num2cell(repmat(':',1,Dmin-1)); midColons=num2cell(repmat(':',1,Dmax-Dmin-1));
            rightColons=num2cell(repmat(':',1,length(inputSize)-Dmax)); %%if(length(rightColons)) rightColons(end)=':'; end
            for iDmin = 1:inputSize(Dmin)
                for iDmax = 1:inputSize(Dmax)
                    outputMatrix(leftColons{:},iDmax,midColons{:},iDmin,rightColons{:}) = inputMatrix(leftColons{:},iDmin,midColons{:},iDmax,rightColons{:});
                end
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%% -Function to Swap Indices & value- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function outputArray = swapIndexAndValue(inputVector)
            if(~(size(inputVector,1)==1 || size(inputVector,2)==1))
                error(['Error: input must be a column/row vector!!! Currently size(inputVector) = ',num2str(size(inputVector))]);
            end
            outputArray = cell(max(inputVector),1);
            for iInputVector=1:length(inputVector)
                outputArray{inputVector(iInputVector)} = [outputArray{inputVector(iInputVector)},iInputVector];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%% Function does 2 checkings,
        %%%%% 1. Check if each row starts from failSequence..
        %%%%% 2. Check the occurance of [preludeSquence,failSequence], in each row.. 
        %%%%% NOTE: 
        %%%%%   1. NaN values will be preserved.
        %%%%% Warnings: 
        %%%%%   1. All inputs must be bvec/bmats (To solve this, 1. correlation has to take care  of '0's. 2. Non-boolean correlations have same value for different input sequence..)  
        %%%%% TODO_Future: Currently input supported is matrix of upto 2D-dimension
        %%%%% TODO_Future: Scope for extension for multiple searchSequence.. 
        function startIndices = findStartIndexOfSequence(input,preludeSquence,searchSequence,options)
            %%%%%%%%% -Defaulting Inputs/Options- %%%%%%%%
            if(nargin<4) options=''; end
            if(gf.isNullString(preludeSquence)) preludeSquence=[]; end
            if(size(searchSequence,1)~=1) searchSequence = searchSequence'; end
            options = gf.createDefaultStructure(options,struct('isOverlappingSupported',true,'preserveShape',false)); %%% TODO: Has to implement these options..
            %%% -Converting all inputs to binary vectors...
            input = gf.convertToBinaryFormat(input,[-1,1]);
            preludeSquence = gf.convertToBinaryFormat(preludeSquence,[-1,1]);
            searchSequence = gf.convertToBinaryFormat(searchSequence,[-1,1]);
            %%%%%%%%% -Input Verification- %%%%%%%% 
            if(~isa(input,'double')) gf.error('input must be a matrix of double !!!'); end
            if(~isa(preludeSquence,'double')) gf.error('preludeSquence must be a matrix of double !!!'); end
            if(~isa(searchSequence,'double')) gf.error('searchSequence must be a matrix of double !!!'); end
            if(min(size(searchSequence))~=1) gf.error('searchSequence must be a 1D matrix !!!'); end
            %%%%%%%%%% -Create Proper Correlational matrix- %%%%%%%%%%
            searchSequence_proper = [preludeSquence,searchSequence];
            searchSequence_proper(searchSequence_proper==0) = -1;
            correlationMatrix = conv2(input,fliplr(searchSequence_proper));
            correlationMatrix = correlationMatrix(:,length(searchSequence):end-length(searchSequence_proper)+1);
            %%%%%%%%%% -Now search for proper indices- %%%%%%%%%
            comparisonMatrix = [ones(size(correlationMatrix,1),1)*[length(searchSequence):length(searchSequence_proper)],length(searchSequence_proper)*ones(size(correlationMatrix,1),size(correlationMatrix,2)-(length(searchSequence_proper)-length(searchSequence)+1))];
            startIndices = find(correlationMatrix==comparisonMatrix);
        end
        
        function mostCorrelatedColumnIndex = getMostCorrelatedColumnIndex(baseColumns,vectorToCorrelate,correlationMethod)
            vectorToCorrelate = gf.convertRowVectorsToColumnVectors(vectorToCorrelate);
            if(strcmpi(correlationMethod,'difference'))
                correlationValuePerBaseColumn = sum(abs(baseColumns-repmat(vectorToCorrelate,[1,size(baseColumns,2)])));
                [~,mostCorrelatedColumnIndex] = min(fliplr(correlationValuePerBaseColumn));
                mostCorrelatedColumnIndex = numel(correlationValuePerBaseColumn)+1-mostCorrelatedColumnIndex;   %% -To make sure tht highest index is outputted, when multiple min values exists...
            else
                error(['Specify Proper correlationMethod!! Currently it''s ',correlationMethod]);
            end
        end
        
        %%%% -Function to get a submatrix with row & column indices as specified in input,,
        function output = getIndicedSubMatrix(inputMatrix,rowIndices,columnIndices)
            output = inputMatrix(sub2ind(size(inputMatrix),rowIndices,columnIndices));
        end

        %%%%%%%%%%%%% -Function to Convert singleDimension-Address to it's original Address- %%%%%%%% 
        %%%% output: Each row corresponds to each value of indicesOfElement... 
        function output = convertSingleIndexToOriginalAddress(matrixSize,indicesOfElement)
            %%% -Input Correction- %%%
            if(sum(size(indicesOfElement)~=1)>1) gf.error('indicesOfElement must be column/row vector!!! indicesOfElement = ',indicesOfElement); end
            %%% -Convert indicesOfElement to Column Vetor...
            if(size(indicesOfElement,2)>1) indicesOfElement = indicesOfElement'; end
            if(max(indicesOfElement)>prod(matrixSize)) gf.error('Out of Index Error for matrixSize = ',matrixSize,',, indicesOfElement = ',indicesOfElement); end
            cumprodFromLeft = cumprod(matrixSize);
            output = rem(floor(((indicesOfElement-1)*ones(1,length(matrixSize)))./(ones(length(indicesOfElement),1)*[1,cumprodFromLeft(1:end-1)])),ones(length(indicesOfElement),1)*matrixSize)+1;
        end
        
        %%%%%%%%%%%%% -Function to Convert singleDimension-Address to it's original Address- %%%%%%%% 
        function output = convertOriginalAddressToSingleIndex(matrixSize,multiDimensionalAddressOfElement)
            cumprodFromLeft = cumprod(matrixSize);
            output = sum((multiDimensionalAddressOfElement-1).*[1,cumprodFromLeft(1:end-1)]) +1;
        end
        
        %%%%%%%%%%%%% -Function to Compute number of elements, per Input- %%%%%%%% 
        function output = getNumElementsPerInput(dimensionsPerInput)
            for iInput=1:length(dimensionsPerInput)
                if(isa(dimensionsPerInput,'cell'))
                    output(iInput) = prod(dimensionsPerInput{iInput});
                else
                    output(iInput) = prod(dimensionsPerInput(iInput));
                end
            end
        end
        
        %%%%%%%%%%%%% -Function to ....- %%%%%%%% 
        function cellsOfCells = convertCellsOfDoubleToCellsOfCells(cellsOfDouble);
            for iRow=1:size(cellsOfDouble,1)
                for iCol=1:size(cellsOfDouble,2)
                    cellsOfCells{iRow,iCol} = num2cell(cellsOfDouble{iRow,iCol});
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Drawing Functions- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% -Create curved arrow,
        %   NOTE: 1. 'from' & 'to' are Nx2 matrix for 2D plot, Nx3 matrix for 3D plot.
        %         2. if 'curvatureRatio = 0', then a flat arrow is drawn!
        function curvedArrow(from,to,spec,curvatureRatio)
            if(gf.isempty('spec')) spec = cell(0,0); end
            if(gf.isempty('curvatureRatio')) curvatureRatio = 0.5; end
%             spec = {'color','b',spec{:}};
            %% -Find distance, direction & midpoint,
            distance = sqrt(sum((to-from).^2,2));
            direction = bsxfun(@rdivide,(to-from),distance);     %% -Unitary vector representing direction
            midpoint = (from+to)/2;
            midpoint(:,1:2) = midpoint(:,1:2) + bsxfun(@times,(curvatureRatio*distance),[-direction(:,2),direction(:,1)]);
            %% -Calculate 'x','y' & 'z' values,
            x = [from(:,1),to(:,1)];
            y = [from(:,2),to(:,2)];
            if(size(from,2)>=3)
                z = [from(:,3),to(:,3)];
            end
            %% -Draw curves,
            hasCallerHoldOnPlot = ishold;
            for iArrow = 1:size(from,1)
                xx = spline([0,0.5,1],[x(iArrow,1),midpoint(iArrow,1),x(iArrow,2)],[0:0.1:1]);
                yy = spline([0,0.5,1],[y(iArrow,1),midpoint(iArrow,2),y(iArrow,2)],[0:0.1:1]);
                if(exist('z','var'))     %% For 3D plot
                    zz = spline([0,0.5,1],[z(iArrow,1),midpoint(iArrow,3),z(iArrow,2)],[0:0.1:1]);
                    plot3(xx,yy,zz,spec{:}); hold on;
                    quiver3(xx(end-1),yy(end-1),zz(end-1),xx(end)-xx(end-1),yy(end)-yy(end-1),zz(end)-zz(end-1),'MaxHeadSize',2,spec{:});
                else                    %% For 2D plot
                    line(xx,yy,spec{:}); hold on;
                    h = annotation('arrow');
                    set(h,'parent', gca,'position', [xx(end-1),yy(end-1),xx(end)-xx(end-1),yy(end)-yy(end-1)], ...      %% -NOTE: 'spec' gets priority over all other parameters
                        'HeadLength', 10, 'HeadWidth', 10, 'HeadStyle', 'cback1',spec{:});
                end
            end
            if(~hasCallerHoldOnPlot) hold off; end
        end

        
        %%%%%%%%%%%%% - Function to plot the diagonal grids - %%%%%%%%%%%%%%%%
        %%%% -Note: Diagonal => line from bottom-left (south-west) to top-right (north-east)
        function gridDiagonal(xlimits,ylimits,xStep,yStep,isCrossDiagonal,color1,color2)
            %% -Input Correction ...
            if(nargin<1)
                error('Specity minimum 1 inputs!!');
            end
            if(nargin<2) ylimits = xlimits; end
            if(nargin<3) xStep = 1; end
            if(nargin<4) yStep = (ylimits(2)-ylimits(1))/ceil((xlimits(2)-xlimits(1))/xStep); end  %% -Making number of Y-points same as X-points ...
            if(nargin<5) isCrossDiagonal = false; end
            if(nargin<6) color1 = [160,150,150]/255; end  %% -[145,255,110]/255;
            if(nargin<7) color2 = [155,160,220]/255; end
            %% -Preserve Caller's 'Hold Status' ON/OFF, even after this function,
            hasCallerHoldOnPlot = ishold();
            hold on;
            %% -Aspect ratio is Width/Height of the graph
            width = xlimits(2)-xlimits(1);
            height = ylimits(2)-ylimits(1);
            aspectRatio = width/height;
            %% -Plot Diagonal grids (bottom-left to top-right), {NOTE: 'X' & 'Y' are computed after rooting to bottom-left point to origin}
            xPoints = 0:xStep:width;
            yPoints = height:-yStep:0;
            X = [xPoints;width*ones(1,numel(xPoints))];
            Y = [zeros(1,numel(yPoints));yPoints];
            plot(xlimits(1)+X,ylimits(1)+Y,'-.','color',color1);                           %% -Potting for right-bottom corner of the graph
            plot(xlimits(1)+Y*aspectRatio,ylimits(1)+X/aspectRatio,'-.','color',color1);   %% -Potting for left-top corner of the graph
            %% -Plot Cross-Diagonal grids (top-left to bottom-right)
            if(isCrossDiagonal)
                xPoints = width:-xStep:0;
                yPoints = height:-yStep:0;
                X = [zeros(1,numel(xPoints));xPoints];
                Y = [yPoints;zeros(1,numel(yPoints))];
                plot(xlimits(1)+X,ylimits(1)+Y,'-.','color',color2);                                             %% -Potting for left-bottom corner of the graph
                plot(xlimits(1)+(height-Y)*aspectRatio,ylimits(1)+(height-X)/aspectRatio,'-.','color',color2);   %% -Potting for right-top corner of the graph
            end
            %% -Bring back caller's hold ON/OFF status ...
            if(~hasCallerHoldOnPlot) hold off; end
        end
        
        %% -Function to connect points in a matrix by it's value in ascending order...
        %   NOTE: '0', 'Nan' and '-ve' values are ignored ...
        function connectingPlot(M,varargin)
            xLimits = xlim; yLimits = ylim;
            M = flipud(M);                  %% -To make plot looks oriented with input matrix,,,
            M(isnan(M) & M<0) = 0;
            [Y,X] = ind2sub(size(M),find(M));
            values = M(M>0);
            [~,order] = sort(values);
            X = X(order); Y = Y(order);
            U = X(2:end); V = Y(2:end);
            X = X(1:end-1); Y = Y(1:end-1);
            U = U-X; V = V-Y;
            if(~isempty(X))
                quiver(X,Y,U,V,0,varargin{:}); %% - '0' is used for no scaling
            end
            %% -Trimming off the arrow heads across the axes boundaries, (Warning: 'xlim manual' & 'ylim manual' are displaying blank figure, whenever ther's no breakpoint!)
            xlim(xLimits); ylim(yLimits);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Binning/Plotting Functions- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% -
        function [nValuesPerBinCenters,binIndexPerValue] = binValues(values,binCenters)
            if(isempty(binCenters))
                error('Error(binValues): 2nd input argument must not be empty');
            elseif(isempty(values))
                nValuesPerBinCenters = zeros(size(binCenters));
                binIndexPerValue = [];
            else
                values1 = gf.convertColumnVectorsToRowVectors(values);
                binCenters1 = gf.convertColumnVectorsToRowVectors(binCenters);
                [binCenters1,sortIndices] = sort(binCenters1);
                [nValuesPerBinCenters,binIndexPerValue] = histc(values1, [-inf binCenters1(1:end-1) + diff(binCenters1)/2 Inf]);
                nValuesPerBinCenters = nValuesPerBinCenters(1:end-1);
                binIndexPerValue = reshape(sortIndices(binIndexPerValue),size(values));
                nValuesPerBinCenters(sortIndices) = nValuesPerBinCenters;
            end
        end
        
        %%
        %   'color' and 'lineStyle' must be a column vector of cells
        %   dims = []
        function plot(xMatrix,yMatrix,dims,colorOrLineStyle1,colorOrLineStyle2)
            if(isempty(yMatrix))  return;  end
            % -Re-arrange 'xMatrix' and 'yMatrix'
            permutationOrder = [dims, setdiff([1:max(ndims(yMatrix),max(dims))],dims)];
            yMatrix = permute(yMatrix,permutationOrder);
            if(~gf.isempty('xMatrix'))
                xMatrix = gf.convertRowVectorsToColumnVectors(xMatrix);
                if(gf.ndims(xMatrix)>1)
                    xMatrix = gf.reshapeByCircularAffix(xMatrix,size(yMatrix));
                    xMatrix = permute(xMatrix,permutationOrder);
                else
                    xMatrix = repmat(xMatrix,1,size(yMatrix,2),size(yMatrix,3));
                end
            else
                xMatrix = [1:size(yMatrix,1)]';
                xMatrix = repmat(xMatrix,1,size(yMatrix,2),size(yMatrix,3));
            end
            assert(ndims(yMatrix)<=3,'Error: We can plot maximum a 3D matrix!!');
            
            % -Create 'lineStyle' and 'colors' for each plot
            isColorOrderFirst = true;
            colors = get(0, 'defaultAxesColorOrder');
            colors = mat2cell(colors,ones(1,size(colors,1)));
            lineStyles = {'-','--','--+','-^','-*','-v'}';
            if(~gf.isempty('colorOrLineStyle1'))
                if(isa(colorOrLineStyle1{1},'double'))
                    colors = colorOrLineStyle1;
                    isColorOrderFirst = true;
                else
                    lineStyles = colorOrLineStyle1;
                    isColorOrderFirst = false;
                end
            end
            if(~gf.isempty('colorOrLineStyle2'))
                if(isa(colorOrLineStyle2{1},'double'))
                    colors = colorOrLineStyle2;
                    isColorOrderFirst = false;
                else
                    lineStyles = colorOrLineStyle2;
                    isColorOrderFirst = true;
                    
                end
            end
            if(size(colors,2)>1) colors = colors'; end
            if(size(lineStyles,2)>1) lineStyles = lineStyles'; end
            if(isColorOrderFirst)
                colors = repmat(colors',[1,1,size(xMatrix,3)]);
                lineStyles = repmat(shiftdim(lineStyles,-2),[1,size(xMatrix,2),1]);
            else
                colors = repmat(shiftdim(colors,-2),[1,size(xMatrix,2)]);                
                lineStyles = repmat(lineStyles',[1,1,size(xMatrix,3)]);                
            end
            
            assert(all(gf.size(xMatrix,[2,3])<=gf.size(lineStyles,[2,3])),'Error: Not enough lineStyles');
            assert(all(gf.size(xMatrix,[2,3])<=gf.size(colors,[2,3])),'Error: Not enough colors');
            % -PLOT
            hasCallerHoldOnPlot = ishold();
            for i3rdDimension = 1:size(xMatrix,3)
                for i2ndDimension = 1:size(xMatrix,2)
                    plot(xMatrix(:,i2ndDimension,i3rdDimension),yMatrix(:,i2ndDimension,i3rdDimension),lineStyles{1,i2ndDimension,i3rdDimension},'color',colors{1,i2ndDimension,i3rdDimension});
                    hold on;
                end
            end
            if(hasCallerHoldOnPlot) hold on; else hold off; end
        end
        
        
        
        function pdfplot(values,binCenters,varargin)
            %%%%%%%%%%%% -Input Check- %%%%%%%%%%%
            if(sum(diff(binCenters)<0)) Error('binCenters must be a increasing 1-D vector of values!!!, binCenters = ',binCenters); end
            % %  -Obsolete:   probablityPerBinCenters = histc(values, [-inf binCenters(1:end-1) + diff(binCenters)/2])/length(values);
            probablityPerBinCenters = gf.binValues(values,binCenters)/length(values);
            plot(binCenters,probablityPerBinCenters,varargin{:});
        end
        
        %%%%%%% -If number of values is less, then we can use following function for cdfplot- %%%%%%%%
        function cdfplotAfterInterpolation(inputValues,varargin)
            interpolationFactor = 10;
            x= sort(inputValues);
            y=1/length(x):1/length(x):1;
            xnew=linspace(min(x),max(x),interpolationFactor*length(x));
            if(isempty(varargin))
                plot(xnew,interp1(x,y,xnew));
            else
                plot(xnew,interp1(x,y,xnew),varargin);
            end
        end
        
        %% -Merge multiple plots-
        %%%% -Note1: First value of input is considered as Sink figure & remaining as Source figures.
        %%%% -Note2: All Source figures are closed.
        %%%% -Note3: Everything in Sink figure get preserved except legend.
        function mergeFigures(figureIndices)
            sinkFigureIndex = figureIndices(1);
            sourceFigureIndices = figureIndices(2:end);
            ax1 = get(figure(sinkFigureIndex), 'Children');
            for sourceFigureIndex=sourceFigureIndices
                ax2 = get(figure(sourceFigureIndex), 'Children');
                for i = 1 : numel(ax2) 
                   ax2Children = get(ax2(i),'Children');
                   copyobj(ax2Children, ax1(i));
                end
                close(figure(sourceFigureIndex));
            end
        end

        %% -Function to plot matrix on a 2D plane
        %   mapColor can be 'gray', 'bone', 'pink' ...etc
        function plotMatrixIn2D(mat,mapColor)
            if(gf.isempty('mapColor'))
                mapColor = 'gray';
            end
            [r,c] = size(mat);                           %# Get the matrix size
            imagesc((1:c)+0.5,(1:r)+0.5,mat);            %# Plot the image
            colormap(flipud(eval(mapColor)));                              %# Use a gray colormap
            axis equal                                   %# Make axes grid sizes equal
            set(gca,'XTick',1:(c+1),'YTick',1:(r+1),...  %# Change some axes properties
                    'XLim',[1 c+1],'YLim',[1 r+1],...
                    'GridLineStyle','-','XGrid','on','YGrid','on');
        end
        
        %% Just like gf.plotMatrixIn2D, but in this function we can input multi-dimensional matrix
        % -NOTES: 
        %   1. Here address is written in reverse order, i.e.,  India, Kerala, Thrissur.
        %   2. addressPerParticle{i} is i^th particle address
        %      addressPerParticle{i} = {11,12,[],14,15}  =>  11th partition in 1st dimension, 12th partition in 2nd dimension
        % -EXAMPLE
        %   1. M=cat(3,[1,0,0;0,0,1],[0,1,0;1,0,1]); gf.plotMatrixLayout(M)
        function plotMatrixLayout(inMatrix,finalDimPerInputDim)
            if(gf.isempty('finalDimPerInputDim'))
                finalDimPerInputDim = rem(0:ndims(inMatrix)-1,2)'+1;
            end
            assert(numel(finalDimPerInputDim) == ndims(inMatrix),'Error: Improper finalDimPerInputDim');
            permDims = [flipud(find(finalDimPerInputDim==1));flipud(find(finalDimPerInputDim==2))];
            outMatrix = reshape(permute(inMatrix,permDims),prod(gf.size(inMatrix,find(finalDimPerInputDim==1))),[]);
            gf.plotMatrixIn2D(outMatrix);
        end
        
        %%%%%%% -Compute the probablity of success, given the values for succeeded & failed case...- %%%%%%%%%%
        function successProbablityPerBinCenters = getSuccessProbablityPerBinCentres(succeededValues,failedValues,binCenters)
            %%%%%%%%%%%% -Input Check- %%%%%%%%%%%
            if(sum(diff(binCenters)<0)) Error('binCenters must be a increasing 1-D vector of values!!!, binCenters = ',binCenters); end
            nSuccessValuesPerBinCenters = histc(succeededValues, [-inf binCenters(1:end-1) + diff(binCenters)/2]);
            nFailedValuesPerBinCenters = histc(failedValues, [-inf binCenters(1:end-1) + diff(binCenters)/2]);
            nFailedValuesPerBinCenters(find((nSuccessValuesPerBinCenters+nFailedValuesPerBinCenters)==0)) = 1;  %%% Avoiding division by '0' error by replacing 0's with 1's..
            successProbablityPerBinCenters = nSuccessValuesPerBinCenters./(nSuccessValuesPerBinCenters+nFailedValuesPerBinCenters);
            successProbablityPerBinCenters = successProbablityPerBinCenters(1:end-1);
        end
        
        %%%% -NOTE: Input must be a vector
        function varargout = getDuplicateValues(input)
            [~,inputDimension] = max(size(input));
            %% -Converting input to a Column vector & get duplicate values
            input = input(:);    
            input = sort(input);
            duplicateElementsIndices = [diff(input)==0;false] | [true;diff(input)==0];
            output = unique(input(duplicateElementsIndices));
            % -Reshaping output to similar to input shape,
            permuteIndices = 1:numel(size(input));
            permuteIndices([1,inputDimension]) = permuteIndices([inputDimension,1]);
            output = permute(output,permuteIndices);
            duplicateElementsIndices = permute(duplicateElementsIndices,permuteIndices);
            % -Assigning to output
            varargout{1} = output;
            if(nargout>1) varargout{2} = duplicateElementsIndices; end
        end
        
        function varargout = getNonRepeatingValues(input)
           [varargout{1:nargout}] = gf.getDuplicateValues(input);
           varargout{1} = setdiff(input,varargout{1});
        end
        
        
        %% -Function to delete rows,columns,3rd dimensions .... etc, from inputMatrix
        %%% -varargin{1} indicate rowIndices to delete, similarly varargin{n} indicate n'th dimension's indices to delete.
        function outputMatrix = trimMatrix(inputMatrix,varargin)
            outputMatrix = inputMatrix;
            numOfDimensions = numel(size(outputMatrix));
            cellsOfColons = num2cell(repmat(':',1,numOfDimensions));
            for iDimension=1:length(varargin)
                outputMatrix(cellsOfColons{1:iDimension-1},varargin{iDimension},cellsOfColons{iDimension+1:end}) = [];
            end
        end        

        %% -For Multiplying multii-dimensional matrices...
        function C = mmat(A,B,dim)
            % Simple matrix multiplication of multidimensional arrays.
            %
            % Input:
            % A, B      Multidimensional input arrays.
            % dim       Contains two number, that selects two dimensions.
            %
            % The multiplication is a standard matrix multiplication. The two matrices
            % are selected by dim:
            %   AB = A(... dim(1) ... dim(2) ...) * B(... dim(1) ... dim(2) ...)
            % The necessary condition that the multiplication can be performed:
            %   size(A,dim(2)) = size(B,dim(1))
            %
            % Singleton dimensions in both A and B matrices are supported.
            %
            % The default value for dim is dim = [1 2].
            %
            % Examples:
            % For 2D matrices mmat is identical to the Matlab built-in multiplication:
            % A = [1 2 3 4];
            % B = [1;2;3;4];
            % C = mmat(A,B)
            %
            % C will be 30.
            %
            % For multidimensional arrays:
            % A = repmat([1 2 3 4],[1 1 5]);
            % B = [1 2 3 4]';
            % C = mmat(A,B)
            % C will be an array with dimensions of 1x1x5 and every element is 30.
            %

            if nargin == 0
                help mmat;
                return;
            end

            if (nargin < 3)
                dim = [1 2];
            end

            if numel(dim)~=2
                error('sw:sw_matmult:WrongInput','dim has to be a two element array!');
            end

            if size(A,dim(2)) ~= size(B,dim(1))
                error('sw:sw_matmult:WrongInput','Wrong input matrix sizes!');
            end

            nDA = ndims(A);
            nDB = ndims(B);
            nD = max(nDA,nDB);

            nA = [size(A),ones(1,nD-nDA)]; nA = nA(dim); 
            nB = [size(B),ones(1,nD-nDB)]; nB = nB(dim);

            % form A matrix
            % (nA1) x (nA2) x nB2
            A = repmat(A,[ones(1,nD) nB(2)]);
            % form B matrix
            % nA1 x (nB1) x (nB2)
            idx = 1:nD+1; idx([dim end]) = idx([end dim]);
            repv = ones(1,nD+1); repv(dim(1)) = nA(1);

            B = repmat(permute(B,idx),repv);

            % multiply with expanding along singleton dimensions
            C = sum(bsxfun(@times,A,B),dim(2));

            idx2 = 1:nD+1; idx2([dim end]) = idx2([dim(1) end dim(2)]);

            % permute back the final result to the right size
            C = permute(C,idx2);

        end
        
        %% -Similar to above function mmat,  (TODO_Future: Extend for all dimensions)
        function output = kron(input1,input2)
            if(size(input1,2)==1 && size(input2,1)==1 && size(input2,3)==1)             % When input1 is along 1st & 3rd dimension & input2 is present only in 2nd dimension
                output = gf.mmat(input1,input2);
            elseif(size(input1,3)==1 && size(input2,1)==1 && size(input2,2)==1)         % When input1 is 2D & input2 is present only in 3rd dimension
                output = gf.mmat(permute(input1,[1,3,2]),permute(input2,[1,3,2]));
                output = permute(output,[1,3,2]);
            end
        end
        
        
        %% -Function to generate all combinations of +ve integers whose sum is a Fixed number
        %%%% output: Each column sum would be finalSum
        function output = getAllNumberSetsWithFixedSum(finalSum,numOfIntegers)
            output = [];
            %% -TERMINATION
            if(numOfIntegers<=1) output = [finalSum]; return; end
            if(finalSum==0) output = zeros(numOfIntegers,1); return; end
            %% -Recursive Loop
            for currentValue=0:finalSum
               childOutput =  gf.getAllNumberSetsWithFixedSum(finalSum-currentValue,numOfIntegers-1);
               output = [output, [currentValue*ones(1,size(childOutput,2));childOutput]];
            end
        end
        
        %% -Problem Statement: U have 'N' liquids with different concentration. U want to create a final liquid with a desired concentration(within tolerance limit) by mixing 'N' liquids.
        %%%% .. The restriction here is that, each liquid can be taken (for mixing) with integer count only (which includes '0').
        %%%% -Get all the combination of liquids for which desiredConcentration is achieved within tolerance limit
        %%%% -Each column in output contains volumes per liquid-type for mixing.
        function output = getAllCombinationOfMixing(concentrationPerLiquidType,totalVolume,desiredConcentration,tolerance)
            output = [];
            if(nargin<4 || strcmp(tolerance,'')) tolerance = 0; end
            volumePerLiquidTypePerCombination = gf.getAllNumberSetsWithFixedSum(totalVolume,numel(concentrationPerLiquidType));
            concentrationPerCombination = (gf.convertColumnVectorsToRowVectors(concentrationPerLiquidType)*volumePerLiquidTypePerCombination)./sum(volumePerLiquidTypePerCombination);
            desiredCombinationIndices = find(concentrationPerCombination>=desiredConcentration*(1-tolerance) & concentrationPerCombination<=desiredConcentration*(1+tolerance));
            %% -sort desiredCombinationIndices as best combinationIndex first.
            [~,tempIndices] = sort(abs(concentrationPerCombination(desiredCombinationIndices)-desiredConcentration));
            desiredCombinationIndices = desiredCombinationIndices(tempIndices);
            %% -Create output Matrix, with each column corresponds to a combination
            for desiredCombinationIndex = desiredCombinationIndices
                currentOutput=[];
                for iLiquid=1:size(volumePerLiquidTypePerCombination,1)
                    currentOutput = [currentOutput;iLiquid*ones(volumePerLiquidTypePerCombination(iLiquid,desiredCombinationIndex),1)];
                end
                output = [output currentOutput];
            end
        end

        %% -Same as above fucntions, but returns only the first combination
        function output = getFirstCombinationOfMixing(concentrationPerLiquidType,totalVolume,desiredConcentration,tolerance)
            if(nargin<4 || strcmp(tolerance,'')) tolerance = 0; end
            [concSorted,sortedIndex] = sort(concentrationPerLiquidType);
            if(desiredConcentration<concSorted(1) || desiredConcentration>concSorted(end))
               error(['error(getFirstCombinationOfMixing): desiredConcentration is too high/low to achieve!!!']); 
            end
            averageVolumeForLiquid = totalVolume/numel(concentrationPerLiquidType);
            numLiquidTypesWithCeiledVolume = rem(totalVolume,numel(concentrationPerLiquidType));
            volumePerLiquidType = [ceil(averageVolumeForLiquid)*ones(numLiquidTypesWithCeiledVolume,1);floor(averageVolumeForLiquid)*ones(numel(concentrationPerLiquidType)-numLiquidTypesWithCeiledVolume,1)];
            closestLiquidTypeIndex = max(find(concSorted<=desiredConcentration));    %% -Closest & less than desiredConcentration
            if(closestLiquidTypeIndex==numel(concentrationPerLiquidType)) %% -Make sure upper partition has got some elements
                closestLiquidTypeIndex = closestLiquidTypeIndex-1;
            end
            itr=1;
            while(itr<10e4)
                itr = itr+1;
                achievedConcentations = sum(volumePerLiquidType.*concSorted)/sum(volumePerLiquidType);
                if(achievedConcentations>=desiredConcentration*(1-tolerance) && achievedConcentations<=desiredConcentration*(1+tolerance))
                    output(sortedIndex)=volumePerLiquidType;
                    return;
                else
                    lowerLiquidIndex = randi([1,closestLiquidTypeIndex],1);
                    upperLiquidIndex = randi([closestLiquidTypeIndex+1,numel(concentrationPerLiquidType)],1);
                    if(achievedConcentations<desiredConcentration)
                        if(volumePerLiquidType(lowerLiquidIndex)<=0) continue; end
                        volumePerLiquidType(lowerLiquidIndex) = volumePerLiquidType(lowerLiquidIndex)-1;
                        volumePerLiquidType(upperLiquidIndex) = volumePerLiquidType(upperLiquidIndex)+1;
                    else
                        if(volumePerLiquidType(upperLiquidIndex)<=0) continue; end
                        volumePerLiquidType(lowerLiquidIndex) = volumePerLiquidType(lowerLiquidIndex)+1;
                        volumePerLiquidType(upperLiquidIndex) = volumePerLiquidType(upperLiquidIndex)-1;                        
                    end
                end
            end
            disp(['Warning(getFirstCombinationOfMixing): All iterations over!!,, Couldn''t find a combination!!']);
        end
                    
        
        %% -Pick each column of input matrix, replace that column with all permutations of that column.
        %%%% -Note: Outputted permutations would be unique
        function output = expandEachColumnToallPermutations(input)
           output = [];
           for iColumn = 1:size(input,2)
               output = [output,perms(input(:,iColumn))'];
           end            
           output = unique(output','rows')';
        end
        
        %% - Input is a cell array of column vectors (ie {[1,2,3]',[4,5]',[6]'})
        %%%% Output is a matrix of size (#Parameters,#TotalCombination) ie [[1,4,6]',[1,5,6]',[2,4,6]',[2,5,6]',[3,4,6]',[3,5,6]']
        function valuePerParameterPerCombination = getAllParameterCombinations(valuesPerParameter)
            % - First generate the index matrix ...
            numParametersPerParameter = cellfun(@(X) size(X,1),valuesPerParameter)';
            numParametersPerParameter=gf.convertRowVectorsToColumnVectors(numParametersPerParameter);
            numParameters = numel(valuesPerParameter);
            numCombinations = prod(numParametersPerParameter);
            [parameterIndexPerParameterPerCombination{1:numParameters}] = ind2sub(fliplr(numParametersPerParameter'),[1:numCombinations]');
            parameterIndexPerParameterPerCombination = fliplr(cell2mat(parameterIndexPerParameterPerCombination))';  %% -Same as [1,1,1; 1,1,2; ...; 3,3,3;];  %% Assumption nUE doesn't get changed ...
            % -Create output matrix
            valuePerParameterPerCombination = nan(numParameters,numCombinations);
            for iParam = 1:numParameters
                valuePerParameterPerCombination(iParam,:) = valuesPerParameter{iParam}(parameterIndexPerParameterPerCombination(iParam,:));
            end
        end
        
        %% -TODO: Make inverse function of the above function, ie given a column vector, find out the column index
        %%%% NOTE: valuesPerParameter is an cell array of column vectors, and valuePerParameter is a column vectors
        function outputIndex = getIndexFromAllParameterCombinations(valuesPerParameter,valuePerParameter)
            valuePerParameterPerCombination = gf.getAllParameterCombinations(valuesPerParameter);
            outputIndex = find(sum(bsxfun(@eq,valuePerParameterPerCombination,valuePerParameter)) == numel(valuePerParameter));
        end
        
        %% If you want to rearrange rows or columns or along any dimension, this function is for you.
        function out = permuteAlongAdimension(in,permSequence,dimIndex)
            if(gf.isempty('dimIndex'))  dimIndex = 1; end
            colonPrefix = repmat(':,',[1,dimIndex-1]);
            colonSuffix = repmat(',:',[1,ndims(in)-dimIndex]);
            eval(['out = in(',colonPrefix,'permSequence',colonSuffix,');']);
        end
            
        %% -Same as Sub-block interleaver in LTE. 
        %       Input are filled into a matrix with nCols=blockLength in each row-by-row.
        %           then it's read each column by column.
        function out = interleaveBlock_old(in,blockLengths)
            out = nan(numel(in),numel(blockLengths));
            for iblockLengths = 1:numel(blockLengths)
                nCols = blockLengths(iblockLengths);
                nRows = ceil(numel(in)/nCols);
                buffer = nan(nCols,nRows);
                buffer(1:numel(in)) = in(:);
                out(:,iblockLengths) = gf.noNaNs(buffer');
            end
        end
        
        %% -Same as Sub-block interleaver in LTE. 
        %       Input are filled into a matrix with nCols=blockLength in each row-by-row.
        %           then it's read each column by column.
        % -Inputs:
        %   'in' is a matrix (or column vector), and we are permuting columns of this matrix.
        % -Mandatory Inputs:
        %   blockLengths = 0  =>  Error exit
        function out = interleaveBlock(in,blockLength)
            m = size(in,1);
            mDesired = blockLength*ceil(m/blockLength);
            permSequence = gf.vec(reshape(1:mDesired,blockLength,[])');
            permSequence(permSequence>m) = [];
            out = gf.permuteAlongAdimension(in,permSequence);
        end
            
        %% Input 'in' must be a matrix (or a single column), and the function interleaves each column
        %       [1 2 3 4 5 6 7 8] -> [2,1,4,3  6,5,8,7]  ->  [2,4,1,3  6,8,5,7]  ->  [2,6,4,8, 1,5,7,3]
        function out = interleaveMaxSpread(in)
            m = size(in,1);  n = size(in,2);
            % -Concatenate NaN rows to make up #Rows = 2^x where x is an integer
            mDesired = 2^(ceil(log2(m)));
            permSequence = [1:mDesired];
            permSequence = gf.vec([2:2:numel(permSequence); 1:2:numel(permSequence)]);    %% First level of permutation
            for i = 2:log2(mDesired)                                %% Each stage of permutation
                permSequence = reshape(permSequence,2^i,[]);
                permSequence = gf.interleaveBlock(permSequence,2^(i-1));
            end
            permSequence(permSequence>m) = [];
            out = gf.permuteAlongAdimension(in,permSequence);
        end
        
        %% Given an input vector, we want to choose 'outLength' elements from it, which are maximally spreaded out within the input
        function out = chooseMaximallySpreadElements(in,outLength)
            inLength = size(in,1);
            if(outLength<inLength)
                inSize = size(in);
                in = reshape(in,inLength,[]);       %% Make 'in' a 2D matrix
                w = (inLength-1)/(outLength-1);
                out = in(round(1:w:inLength),:);
                out = reshape(out,[outLength,inSize(2:end)]);
            else
                out = in;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - ACCESSORY FUNCTIONS - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% - Output of following function is a string, used while saving files with time-stamp 
        %   Update: Following is an OBSOLETE function. Use 'clock' and 'etime' instead.  
        function timeStamp = time()
            timeStamp = datestr(now,'yyyy-mm-dd-HHMMSS');
        end
        
        %% - Output a string indicating current time, 
        %   NOTE: Difference from gf.time() is that, this one has got high precision in 'milliseconds'.
        %   Warning: if date-format is changed, then do corresponding changes in function 'computeTimeDifference'
        %   Update: Following is an OBSOLETE function. Use 'clock' and 'etime' instead.  
        function currentTime = getCurrentTime()
            currentTime = datestr(now,'yyyy-mm-dd HH:MM:SS.FFF');
        end

        %% - Output a string indicating time difference,
        %   Update: Following is an OBSOLETE function. Use 'clock' and 'etime' instead.  
        function timeDifference = computeTimeDifference(dateStr1,dateStr2)
            t1=datetime(dateStr1,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            t2=datetime(dateStr2,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            timeDifference = t1 - t2;
            timeDifference = datestr(timeDifference,'HH:MM:SS.FFF');    %% -Convert 'duration' object to 'string'
        end
        
        %% -Save all variables upto 'numUpScope' parent scopes, into the filename specified. Very usefull while debugging :-)
        %   NOTE: The first input argument fileName must contain '.mat' extension
        % -TODO: Untill now, I couldn't find out a way for saving variables in all scope, to implement 'numUpScope' feature!
        function saveAllVariables(fileName,numUpScope)
           if(gf.isempty('fileName')) 
                fileName = [gf.time,'_',gf.functionNameCaller,'.mat'];    %% TODO: Change to caller function name!
            end
            if(gf.isempty('numUpScope')) numUpScope = 1; end
            evalin('caller',['save(''',fileName,''')']);
        end
        
        %% -Complementary of gf.saveAllVariables.  Very usefull while debugging :-)
        %   Load all variables from recent saved .mat file
        function loadAllVariables(fileName)
           if(gf.isempty('fileName')) 
               [~,filesAll]=system('ls -1 -t *.mat');
               filenamePattern = '[0-9a-zA-Z](\w|\.|\-)*';
               filesAll = regexp(filesAll,filenamePattern,'match');
               fileName = filesAll{1};
           end
            evalin('caller',['load(''',fileName,''')']);
        end
        
        %% - Save all variables in caller's workspace into a structure ...
        function out = workspaceToStruct()
            variables = evalin('caller','who');
            for ivariable = 1:numel(variables)
               out.(variables{ivariable}) = evalin('caller',variables{ivariable});
            end
        end
        
        %% -Same as MATLAB inbuilt error function, except that it will also save all variables in caller's function to a file named '<DATE>_<FUNCTION-NAME>_error.mat'
        %   NOTE: Inputs must be strings to print.
        function Error(varargin)
            fileName = [gf.time,'_',gf.functionNameCaller,'-Error.mat'];
            evalin('caller',['save(''',fileName,''')']);
            save(fileName,'fileName','-append');    %% -Appending the '.mat' file name also into the saved workspace
            gf.Disp(['Error(',fileName,'): '],varargin{:});
            error(' ');     %% NOTE: blank space is mandatory for error function to work ...
        end
        
        %% Details
        %   1. This function saves all the local variables upon error catching.
        %       if 'ME' object is rethrown (after calling this function), then this function brings all the variables in the .mat file created during previous call,
        %           and create a new .mat file after deleting the existing one.
        %   2. Has to be called in catch block in exceptions as,
        %       ME = gf.saveWorkspaceUponException(ME);
        function out = saveWorkspaceUponException(in)
            workspaceOfCaller = evalin('caller','gf.workspaceToStruct()');
            workspaceOfCaller.ME = in;
            if(exist(in.message,'file')) 
                workspacePerScope = load(in.message,'workspacePerScope');
                workspacePerScope = workspacePerScope.workspacePerScope;
                delete(in.message);
            else
                workspacePerScope = {};
            end
            workspacePerScope = {workspaceOfCaller, workspacePerScope{:}};
            workspacePerScope = workspacePerScope';
            % Save to file
            uniqueID = datestr(now,'yyyymmddHHMMssFFF');
            dbstack_ = dbstack(); functionName = dbstack_(2).name;
            filenameToSave = [pwd,'/',uniqueID,'-error-',functionName,'.mat'];
            save(filenameToSave,'workspacePerScope');
            out = MException('MATLAB:function_call:gf:saveWorkspaceUponException',filenameToSave);
        end
        
        %% Function to read csv file as matrix of cells.
        %   Got this code from http://stackoverflow.com/questions/4747834/import-csv-file-with-mixed-data-types 
        function lineArray = csvread(fileName,delimiter)
          if(gf.isempty('delimiter')) delimiter = ','; end
          fid = fopen(fileName,'r');   %# Open the file
          lineArray = cell(100,1);     %# Preallocate a cell array (ideally slightly
                                       %#   larger than is needed)
          lineIndex = 1;               %# Index of cell to place the next line in
          nextLine = fgetl(fid);       %# Read the first line from the file
          while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
            lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
            lineIndex = lineIndex+1;          %# Increment the line index
            nextLine = fgetl(fid);            %# Read the next line from the file
          end
          fclose(fid);                 %# Close the file
          lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
          for iLine = 1:lineIndex-1              %# Loop over lines
            lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                                'Delimiter',delimiter);
            lineData = lineData{1};              %# Remove cell encapsulation
            if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
              lineData{end+1} = '';                     %#   ends with a delimiter
            end
            lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
          end
          %% Following is added by Anver, in order to conver str to number, wherever possible 
          isNumber = ~isnan(str2double(lineArray));
          lineArrayDouble = num2cell(str2double(lineArray));
          lineArray(isNumber) = lineArrayDouble(isNumber);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - FUNCTIONS HANDLING FUNCTIONS - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% -Output the function name. Useful for printing while error-exiting!
        function functionNameStr = functionName()
            a = dbstack();
            functionNameStr = a(2).name;
        end
        
        %% -Output the function name. Useful for printing while error-exiting!
        function functionNameCallerStr = functionNameCaller()
            a = dbstack();
            if(numel(a)<3)      %% -Calling from workspace
                functionNameCallerStr = 'MATLAB-Console';
            else
                functionNameCallerStr = a(3).name;
            end
            functionNameCallerStr = regexprep(functionNameCallerStr,'/','_');   %% If caller function is a nested function, then it would be of format "function/nestedFunction" 
        end
        
        %%  Warning: Presently our assumption is functionHandle returns only a single output or array of cells.  
        %   To call functionHandle(a,b) and functionHandle(c,d),  make inArgumentCellMatrix = {{a,b},{c,d}} 
        %   'ordering' = {'-forward','-reverse'} tells which order the iterations have to be done!
        function outArgumentCellMatrix = executeIteratively(functionHandle,inArgumentCellMatrix,nThreads,ordering)
            % Create order of simulation
            if(gf.isempty('ordering')) ordering = '-forward'; end
            n = numel(inArgumentCellMatrix);
            order = 1:n;
            if(strcmpi(ordering,'-reverse')) order = fliplr(order); end
            [~,orderReverse] = sort(order);
            inArgumentCellMatrix(orderReverse) = inArgumentCellMatrix;      %% Used orderReverse in LHS in order to preserve the shape of inArgumentCellMatrix
            % Run each threads
            if(gf.isempty('nThreads') || nThreads<=1) nThreads = 0; end
            outArgumentCellMatrix = cell(size(inArgumentCellMatrix));
            if(nThreads<=1)
                for i = 1:n         %% To capture error in local machine
                    outArgumentCellMatrix{i} = functionHandle(inArgumentCellMatrix{i}{:});
                end                
            else
                parfor(i = 1:n,nThreads)
                    outArgumentCellMatrix{i} = functionHandle(inArgumentCellMatrix{i}{:});
                end
            end
            outArgumentCellMatrix(order) = outArgumentCellMatrix;            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - SYMBOLIC FUNCTIONS - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Function useful while creating Linear Programming Models  
        % Parameters:
        % Inputs:
        %   'symExpressions': Array of symblic expressions
        %   'variables'     : Array of symbolic variables
        % Output:
        %   'coefficientsPerVariable': Numerical matrix with ith row corresponds to ith expression in 'symExpressions', jth column corresponds to jth variable 
        function coefficientsPerVariable = findCoefficientsPerVariable(symExpressions,variablesALL)
            coefficientsPerVariable = zeros(numel(symExpressions),numel(variablesALL));
            for i = 1:numel(symExpressions)     %% For each expression
                [coefficients,variables] = coeffs(symExpressions(i));
                [tmp,indices] = ismember(variablesALL,variables);                   %% Observation: 10seconds for simulation when numel(variablesALL)=9600 and numel(variables)=4
                coefficientsPerVariable(i,find(tmp)) = coefficients(indices(indices>0));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - MATLAB HIGH LEVEL FUNCTIONS TWEAKS - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% - Same as linprog, except that if the problem is infeasible, then find a solution satisfying max equations
        %   Details:
        %       1. Problem 'Ax <= b' is converted into '[A,-K*I][x,y]' <= b' ,, whereas 'y' is the infeasiblity indicator (y=1 => that equation is infeasible)
        %   -NOTE:
        %       1. Both inputs and outputs are same as MATLAB linprog
        %   -Warnings: 
        %       1. Objective function scaling (f) must be choosen very less than '1' to have the objectivity of maximizing number of equations
        %       2. When using Gurobi, with params params.FeasibilityTol = params.IntFeasTol = 1e-9,  I observed around 5% constraints are not satisfied,
        %           and observed infeasibility tolerance for those constraints lies in [10^-12, 10^-10]
        function [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,varargin)
            
            [m,n] = size(A);
            
            % -Create proper 'K' value
            K = ceil(max(sum(abs(A),2) + abs(b)));
            if(~isempty(lb) && ~isempty(ub)) K = ceil(K+max(abs(A)*abs(mean([lb,ub],2)))); end
            K = 10^ceil(log10(K));      %% ceiling to the closest 10^x
            
            % See if '-gurobi' is specified in inputs or not. Assign inputted params 
            useGurubi = false;
            tempIndex = gf.contains(varargin,'-gurobi',true,true);
            if(~isempty(tempIndex))
                useGurubi = true;
                varargin(tempIndex{1}) = [];
            end
            params = struct();
            tempIndex = gf.contains(varargin,'-params',true,true);
            if(~isempty(tempIndex))
                params = varargin{tempIndex{1}+1};
                varargin(tempIndex{1}:tempIndex{1}+1) = [];
            end
            
            % -Manipulate 'b'  (Compensating for 'IntFeasTol' and 'FeasibilityTol')     %% TODO: Do similar thing for matlab intlinprog()
            bOriginal = b;
            b = b - K*params.IntFeasTol - params.FeasibilityTol;      
            
            % -Manipulate 'A' matrix,
            A = [A,-K*eye(size(A,1))];

            % -Manipulate f
            if(gf.isempty('f')) f = zeros(1,n); end
            f = [f,ones(1,size(A,1))];
                        
            % -Create lb & ub
            if(gf.isempty('Aeq')) Aeq = []; end
            if(gf.isempty('beq')) beq = []; end
            if(gf.isempty('lb')) lb = -Inf(n,1); end
            if(gf.isempty('ub')) ub = Inf(n,1); end
            lb = [lb; zeros(m,1)];
            ub = [ub; ones(m,1)];
            
            % -Create 'intcon'
            intcon = n+1:m+n;
            % -Solve Mixed Integer Linear Program
            if(useGurubi)
                    model = struct('obj',f,'A',sparse(A),'rhs',b,'lb',lb,'ub',ub);
                    model.vtype(1:m+n) = 'C';
                    model.vtype(intcon) = 'B';
                    [x,exitflag,output] = gurobi_(model,params);
                    % -Debug
                    %assert(all(x(intcon)==0 | x(intcon)==1),'Error: Non integer variables !!');
                    %assert(all(A*x<=bOriginal),'Error: Constrains violation');
                    %output = [];
            else
                [x,~,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,varargin{:});
            end
            
%             %%%% -Debuging (ie. Looking for the observed constraints infeasibility tolerance)
%             excessPowerPerLink=b-A*x;
%             excessPowerPerSuccessClaimedLinks = excessPowerPerLink(x(21:end)==0);
%             excessPowerPerSuccessClaimedAndFailedLinks = excessPowerPerSuccessClaimedLinks(excessPowerPerSuccessClaimedLinks<0);
%             figure; cdfplot(log10(-excessPowerPerSuccessClaimedAndFailedLinks));
            
            % -Assigning output
            if(exitflag~=1)
                x = zeros(numel(f),1);           %% -Dummy output ...
            end
            
            x = x(1:n); 
            fval = f(1:n)*x;    %% NOTE: 'f(n+1:end)*x(n+1:end)' is the number of satisfied constraints.
        
        end
        
        
        %% input/output are in the same format as MATLAB inbuilt function 'intlinprog'
        %   Here we do exhaustive enumeration of all integer combinations of 'x', and solve the resulting linear-programs for all the enumerations. This is very time taking procedure!! 
        %   Split A = [B,Asub] where columns in 'B' corresponds to integer-constraints, and columns in 'C' corresponds to non-integer constrains
        %   Warning(TODO): Currently inputted 'Aeq', 'Beq', 'options' are discarded!
        function varargout = intlinprogExhaustive(f,intcon,A,b,Aeq,beq,lb,ub,options)
            isIntCons = zeros(size(A,2),1);
            isIntCons(intcon) = 1;
            isIntCons = logical(isIntCons);
            if(isempty(f)), f = zeros(size(A,2),1); end
            fSub = f(~isIntCons);
            B = A(:,isIntCons);
            Asub = A(:,~isIntCons);
            lbSub = lb(~isIntCons);
            ubSub = ub(~isIntCons);
            xBCombinations = gf.createAllCombinationOfIntegerVectors(lb(isIntCons),ub(isIntCons));     %% (#Integer-Constraints x #Combinations)
            if(~gf.isempty('options.isRemoveSymmetry') && options.isRemoveSymmetry)
                xBCombinations = gf.removeSymmetricColumns(xBCombinations);
            end
            optimumValue = Inf;
            x = nan(size(A,2),1);
            for iCombination = 1:size(xBCombinations,2)
                xB = xBCombinations(:,iCombination);
                bSub = b-B*xB;
                [xSub,~,exitflag,output] = linprog(fSub,Asub,bSub,[],[],lbSub,ubSub,'',struct('Display','none'));
                x(isIntCons) = xB;   x(~isIntCons) = xSub;
                val = f'*x;
                if(exitflag==1 && val<optimumValue) 
                    optimumX = x;
                    optimumExitFlag = exitflag;
                    optimumOutput = output;
                    optimumValue = val;
                    if(~any(f))       %% -If 'f' is all-zero vector, then it's only feasiblity check problem!
                        break;
                    end
                end
            end
            x = optimumX;
            varargout{1} = x;
            if(nargout>=2) varargout{2} = optimumValue; end
            if(nargout>=3) varargout{3} = optimumExitFlag; end
            if(nargout>=4) varargout{4} = optimumOutput; end
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -TRIAL FUNCTIONS - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        
    end 
end
