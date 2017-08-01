%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -BRIEF DESCRIPTION:
%   One can control what should be the output of zero/nan indexing by using 
%       the properties 'zeroIndexedElement' and 'nanIndexedElement' of this class.   
%   NOTE: If both 'zeroIndexedElement' and 'nanIndexedElement' are empty, then this class would act as a regualar 'double' class!
%
%%%% -DETAILED DESCRIPTION:
%       1. 'unAssignedValue': MATLAB assigns 0s for un-assigned values. Here we can set that value manuallly
%       2. 'isPreferColumnVector': MATLAB creates column vector upon assigning new indices for a scalar. Setting isPreferColumnVector=true, changes the behaviour
%
%%%% -DEFINITIONS
%       1. Proper Matrix: A matrix, not a vector along any dimension
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1. a = matrixGeneralized([11,12,13],121,122); 
%           Place '121' on zero indices, and 122 on 'NaN' indices 
%       2. b = matrixGeneralized(a,0,nan); 
%           Same as 'a' but with different outputs upon 0/NaN indices 
%       3. b = matrixGeneralized([],[],[],-1); 
%           b(3)=123  =>  b = [-1,-1,123]
%       4. b=matrixGeneralized([],a)
%           Creates empty matrix 'b' with properties of 'a'
%          b=matrixGeneralized(a)
%           Creates matrix 'b' with same data of 'a', but with default properties
%          b=matrixGeneralized(a,a)
%           Creates matrix 'b' with same data of 'a', also with properties of 'a'. (Same as "b=a")
%
%%%% -NOTES:
% 
%
%%%% -NOTES (Programming):
%
%%%% -TODO:
%       1. Handle single element input
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -
classdef matrixGeneralized < double
    properties
        zeroIndexedElement;
        nanIndexedElement;
        unAssignedValue;
        isPreferColumnVector = false;
    end
    
    methods(Access=private)
        
        function this = setPropertiesUsingPropertyHolder(this,propertyHolder)
            this.zeroIndexedElement = propertyHolder.zeroIndexedElement;
            this.nanIndexedElement = propertyHolder.nanIndexedElement;
            this.unAssignedValue = propertyHolder.unAssignedValue;
            this.isPreferColumnVector = propertyHolder.isPreferColumnVector;
        end
    
    end
    
    methods
        
        % varargin = {zeroIndexedElement,nanIndexedElement,unAssignedValue,isPreferColumnVector}
        function obj = matrixGeneralized(data,varargin)
            if(gf.isempty('data')) data = []; end
            obj = obj@double(data);     %% -NOTE: This kind of construction prevents any properties defined inside this class!!
            obj = setProperties(obj,varargin{:});
        end
       
        %% -Call this as a = setProperties(a)
        %   NOTE: 
        %       1. propertyHolder is an object of type 'matrixGeneralized', from which we extract the properties.
        %   Warnings:
        %       1. '[]' is not a empty input, but instead used to clear that argument.
        function this = setProperties(this,zeroIndexedElementOrPropertyHolder,nanIndexedElement,unAssignedValue,isPreferColumnVector)
            if(nargin>=2)                                                       %% Used nargin instead of gf.isempty(), since we need to accept []
                if(isa(zeroIndexedElementOrPropertyHolder,'matrixGeneralized'))
                    propertyHolder = zeroIndexedElementOrPropertyHolder;
                    this = setPropertiesUsingPropertyHolder(this,propertyHolder);
                    assert(nargin<=2,'Error: Extra inputs would be disgarded, while setting properties using propertyHolder!!');
                else
                    zeroIndexedElement = zeroIndexedElementOrPropertyHolder;
                    if(nargin>=2 && ~gf.isNullString(zeroIndexedElement))  this.zeroIndexedElement = zeroIndexedElement; end
                    if(nargin>=3 && ~gf.isNullString(nanIndexedElement))  this.nanIndexedElement = nanIndexedElement; end
                    if(nargin>=4 && ~gf.isNullString(unAssignedValue))  this.unAssignedValue = unAssignedValue; end
                    if(nargin>=5 && ~gf.isNullString(isPreferColumnVector))  obj.isPreferColumnVector = isPreferColumnVector; end
                end
            end
        end
        
        %% -Call this as a = clearAllProperties(a)
        function this = clearAllProperties(this)
            propertyHolder = matrixGeneralized();
            this = setProperties(this,propertyHolder);
        end
        
        function objDouble = subsref(obj,s)
            switch s(1).type
                case '()'
                    objDouble = double(obj);
                    currentSize = size(objDouble);                
                    if(numel(s(1).subs)==1)     %% Global Indexing
                        if(~strcmp([s(1).subs{1}],':'))
                            assert(max(s(1).subs{1})<=prod(currentSize),'Error: Index is out of range!!');
                            objDouble = gf.cat(numel(currentSize),objDouble,0,0,'-appendValue',0);
                            if(~isempty(obj.zeroIndexedElement))
                                zeroIndexTo = numel(objDouble)-1;
                                objDouble(zeroIndexTo) = obj.zeroIndexedElement;
                                s(1).subs{1}(s(1).subs{1}==0) = zeroIndexTo;
                            end
                            if(~isempty(obj.nanIndexedElement))
                                nanIndexTo = numel(objDouble);
                                objDouble(nanIndexTo) = obj.nanIndexedElement;
                                s(1).subs{1}(isnan(s(1).subs{1})) = nanIndexTo;
                            end
                        end
                    else                        %% Indexing per dimensions
                        for iDim = 1:numel(currentSize)
                            if(~strcmp([s(1).subs{iDim}],':'))
                                if(~isempty(obj.zeroIndexedElement))
                                    objDouble = gf.cat(iDim,objDouble,obj.zeroIndexedElement,'-appendValue',obj.zeroIndexedElement);
                                    zeroIndexTo = size(objDouble,iDim);
                                    s(1).subs{iDim}(s(1).subs{iDim}==0) = zeroIndexTo;
                                end
                                if(~isempty(obj.nanIndexedElement))
                                    objDouble = gf.cat(iDim,objDouble,obj.nanIndexedElement,'-appendValue',obj.nanIndexedElement);
                                    nanIndexTo = size(objDouble,iDim);
                                    s(1).subs{iDim}(isnan(s(1).subs{iDim})) = nanIndexTo;
                                end
                            end
                        end
                    end
                    objDouble = objDouble(s(1).subs{:});
                    
                case '.'                                    %% Accessing properties like 'zeroIndexedElement', 'nanIndexedElement', ...etc
                    objDouble = eval(['obj.',s(1).subs]);
                    
                otherwise
                    error(['Error(matrixGeneralized): Use bracket always... Currently used ',s(1).type,' !!!']);
            end
        end
       
        
      %%
      % -PROCEDURE:
      %     Crop both 0/NaN index both in LHS indices, and in RHS 'input' matrix    
      function obj = subsasgn(obj,s,input)
           isInputScalar = numel(input)==1;
           switch s(1).type
               case '()'
                objDouble = double(obj);
                currentSize = size(objDouble)';                
                if(ndims(objDouble)<ndims(input))
                    currentSize = gf.augmentMatrix(currentSize,[ndims(input),1],1);
                end
                if(numel(s(1).subs)==1)     %% Global Indexing
                    if(~strcmp([s(1).subs{1}],':'))
                        if(~isempty(obj.zeroIndexedElement))
                            if(~isInputScalar)
                                input(s(1).subs{1}==0) = [];
                            end
                            s(1).subs{1}(s(1).subs{1}==0) = [];
                        end
                        if(~isempty(obj.nanIndexedElement))
                            if(~isInputScalar)
                                input(isnan(s(1).subs{1})) = [];
                            end
                            s(1).subs{1}(isnan(s(1).subs{1})) = [];
                        end
                    end
                else                        %% Indexing per dimensions
                    colStr = repmat(',:',1,numel(currentSize)-1);
                    for iDim = 1:numel(currentSize)
                        if(~strcmp([s(1).subs{iDim}],':'))
                            dimOrder = 1:numel(currentSize);
                            dimOrder([1,iDim]) = [iDim,1];
                            if(~isInputScalar)
                                input = permute(input,dimOrder);
                            end
                            if(~isempty(obj.zeroIndexedElement))            %% Remove '0' indexed elements from RHS
                                if(~isInputScalar)
                                    eval(['input(s(1).subs{iDim}==0',colStr,') = [];']);
                                end
                                s(1).subs{iDim}(s(1).subs{iDim}==0) = [];
                            end
                            if(~isempty(obj.nanIndexedElement))             %% Remove 'NaN' indexed elements from RHS
                                if(~isInputScalar)
                                    eval(['input(isnan(s(1).subs{iDim})',colStr,') = [];']);
                                end
                                s(1).subs{iDim}(isnan(s(1).subs{iDim})) = [];
                            end
                            if(~isInputScalar)
                                input = ipermute(input,dimOrder); 
                            end
                        end
                    end
                end   
                
                if(~isempty(obj.unAssignedValue) && ~(numel(s(1).subs)==1 && isa(s(1).subs{1},'char'))) %% Avoiding the case "a(:) = b"
                    s(1).subs = gf.convertColumnVectorsToRowVectorsRecursively(s(1).subs);
                    for iInputDim = 1:numel(s(1).subs)
                       if(isa(s(1).subs{iInputDim},'char'))
                           s(1).subs{iInputDim} = 1:max([size(objDouble,iInputDim),size(input,iInputDim),1]);
                       end
                    end
                    desiredSize = nanmax(gf.cat(1,s(1).subs{:}),[],2);
                    if(numel(desiredSize)==1)      
                        assert(sum(currentSize>1)<=1,'Error: We can''t augment a proper matrix by using global indexing!! Confusion is along which dimension we need to augment');
                        [~,mainDimension] = max(currentSize);
                        if(mainDimension>1)
                            desiredSize = [ones(1,mainDimension-1),desiredSize];
                        else
                            desiredSize = [1;desiredSize];
                        end
                    else
                        desiredSize = nanmax(gf.cat(2,currentSize,desiredSize),[],2);
                    end
                    objDouble = gf.augmentMatrix(objDouble,desiredSize,obj.unAssignedValue);                        
                end
                
                objDouble(s(1).subs{:}) = input;
                if(obj.isPreferColumnVector && numel(obj)<=1 && numel(input)==1)
                    objDouble = gf.convertRowVectorsToColumnVectors(objDouble);
                end
                obj = matrixGeneralized(objDouble,obj.zeroIndexedElement,obj.nanIndexedElement,obj.unAssignedValue);
                
               case '.'                                    %% Mutating properties like 'zeroIndexedElement', 'nanIndexedElement', ...etc
                   eval(['obj.',s(1).subs,'=input;']);
                
               otherwise
                   error(['Error(matrixGeneralized): Use bracket always... Currently used ',s(1).type,' !!!']);
           end
      end
       
   end
end
