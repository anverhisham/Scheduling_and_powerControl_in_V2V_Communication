%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% This file contains SUPPORTING-FUNCTIONS(sf) %%%%%%%%%%%%%%% 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef sf

    methods(Static)
         
        
    %%    
    %   -INPUTS:
    %       'isSuccessPerRBPerTimePerTx' :  NaN indicate there was no attempted reception.   
    function isSuccessPerTxPerRx = getTxRxSuccessStatusFromReceptionStatus(txVUEPerRBPerTime,isSuccessPerRBPerTimePerRx)
        [nRB,nTimeslots,nRx] = size(isSuccessPerRBPerTimePerRx);
        nTx = nanmax(nRx,nanmax(txVUEPerRBPerTime(:)));
        isSuccessPerTxPerRx = nan(nTx,nRx);
        txVUEPerRBPerTimePerRx = repmat(txVUEPerRBPerTime,[1,1,nRx]);
        rxVUEPerRBPerTimePerRx = repmat(shiftdim(1:nRx,-1),[nRB,nTimeslots,1]);
        % Assign failed links
        isSuccessPerTxPerRx(sub2ind(size(isSuccessPerTxPerRx),txVUEPerRBPerTimePerRx(isSuccessPerRBPerTimePerRx==0),rxVUEPerRBPerTimePerRx(isSuccessPerRBPerTimePerRx==0))) = 0;
        % Assign successful links (NOTE: this should be done after assigning failed links! )
        isSuccessPerTxPerRx(sub2ind(size(isSuccessPerTxPerRx),txVUEPerRBPerTimePerRx(isSuccessPerRBPerTimePerRx==1),rxVUEPerRBPerTimePerRx(isSuccessPerRBPerTimePerRx==1))) = 1;
    end
        
    %%
    %  F = #RBs, acirPerRB = [0.001,0.001,0.001,0.001,0.0005, ...] for LTE-ACIR
    %  Warning: if 'acirPerRB' is inputted, then 'fh' is ignored
    % -Execution:
    %   sf.findSchedulingRBSequence(10, @(x) sum(sqrt(x)))
    function schOrderPerRB = findSchedulingRBSequence(F,fh)     %,acirPerRB)
        if(gf.isempty('fh'))
            fh = @(x) prod(x);
        end
        schOrderPerRB = [1;nan(F-1,1)];
        for iScheduling = 2:F
            scheduledRBs = find(~isnan(schOrderPerRB));
            unScheduledRBs = setdiff([1:F]',scheduledRBs);
            distPerScheduledRBPerUnscheduledRB = abs(bsxfun(@minus,scheduledRBs,unScheduledRBs'));
            utilityPerUnscheduledRB = nan(numel(unScheduledRBs),1);
            for j=1:numel(unScheduledRBs)
%                 if(~gf.isempty('acirPerRB'))
%                     utilityPerUnscheduledRB(j) = -sum(acirPerRB(distPerScheduledRBPerUnscheduledRB(:,j)));
%                 else
                    utilityPerUnscheduledRB(j) = fh(distPerScheduledRBPerUnscheduledRB(:,j));
%                 end
            end
            winningRB = unScheduledRBs(utilityPerUnscheduledRB==max(utilityPerUnscheduledRB));
            %[~,winningRB] = max(utilityPerUnscheduledRB);
            %% If there are multiple RBs with same "max utility", then resovle the tie
            if(numel(winningRB)>1)  
                tieBreakingAlgIndex = 2;                %% Observation (2017/06/21): tieBreakingAlgIndex=2 is outperforming compared to 3, by 0.3%/0.34%/-0.03% for equalPower/KangWang/OptimalPower
                switch(tieBreakingAlgIndex)
                    case 1
                        winningRB = max(winningRB);   %% Pick the highest RB
                    case 2
                        [~,iwinningRB] = max(mean(abs(bsxfun(@minus,winningRB,scheduledRBs')),2));     %% Maximizing the mean distance to all scheduled RBs
                        winningRB = winningRB(iwinningRB);
                    case 3
                        [~,iwinningRB] = max(min(abs(bsxfun(@minus,winningRB,scheduledRBs')),[],2));      %% Maximizing the min distance to all scheduled RBs
                        winningRB = winningRB(iwinningRB);
                end
            end
            schOrderPerRB(winningRB) = iScheduling;
        end
    end
    
    
    function aci = findACIInAnRB(distPerScheduledRB,acirType)
        acirPerRBGap = acir.acirPerRBGap(acirType);
        aci = sum(acirPerRBGap(distPerScheduledRB));
    end
        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MORE GENERIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Find best permutation matrix 'W' s.t. W*LHS > RHS is maximized!
    %   Examples:
    %       1. If f=@gt and isfromRight=true, then LHS*W > RHS
    %       2.  "    "       "      "   false,  "  W*LHS > RHS      (Default)
    %   Rotate LHS from XY to ZY (ie X->Z) plane, and RHS from XY to ZX (ie X->Z & Y->X) plane
    function W = findBestPermutation(f,LHS,RHS,isfromRight)
        global I;
        if(~gf.isempty('isfromRight') && isfromRight)    %% Transforming to 'isfromRight=true' format
            LHS = transpose(LHS);
            RHS = transpose(RHS);
        end
        LHS = permute(LHS,[1,3,2]);
        RHS = permute(RHS,[3,1,2]);
        nSuccesfulLinksPerRBPerVUE = sum(bsxfun(f,LHS,RHS),3);
        [VUEIndexPerRB,~] = hungarianAlgorithm(-nSuccesfulLinksPerRBPerVUE);
        W = I(:,VUEIndexPerRB);
        W = transpose(W);
        if(~gf.isempty('isfromRight') && isfromRight)
            W = transpose(W);
        end
    end
    
    %% - LHS and RHS are vectors. 
    %    Find best lambda, which would maximizes the number of elements in lambda*LHS >= RHS  
    function [lambda,nSuccess] = findBestScalingFactorForInequality(LHS,RHS,bound)
        % -Input patching
        if(~gf.isempty('bound(1)') && ~isnan(bound(1))) lamda_min = bound(1); else lamda_min = -Inf; end
        if(~gf.isempty('bound(2)') && ~isnan(bound(2))) lamda_max = bound(2); else lamda_max = Inf; end
        [LHS,RHS] = gf.convertRowVectorsToColumnVectorsRecursively(LHS,RHS);
        % -Threshold value computations
        [thresholdValues,sortIndices] = sort(RHS./LHS);
        % -Pick values within {lamda_min,lamda_max} only.
        trimIndices = find(thresholdValues<lamda_min | thresholdValues>lamda_max);
        thresholdValues(trimIndices) = [];
        sortIndices(trimIndices) = [];
        thresholdValuesDirection = sign(LHS(sortIndices));
        if(isempty(thresholdValues))
            lambda = (lamda_min + lamda_max)/2;     %% Tweaking here. (Changed 'lamda_min')
            nSuccess = sum(lambda*LHS >= RHS);
            return;
        end
        nSuccessPerThresholdValues = nan(numel(thresholdValues),1);
        for i=1:numel(thresholdValues)
            nSuccessPerThresholdValues(i) = sum(thresholdValuesDirection(1:i)>=0) + sum(thresholdValuesDirection(i:end)<0);
        end
        [nSuccess,ind] = max(nSuccessPerThresholdValues);
        lambda = thresholdValues(ind);
        %% Tweaking below, to get reasonable(middle) lambda values
        if(thresholdValuesDirection(ind)>=0)
            if(ind==numel(thresholdValues(ind)))
                lambda = mean([lambda,lamda_max]);
            else
                lambda = mean([lambda,thresholdValues(ind+1)]);
            end
        else
            if(ind==1)
                lambda = mean([lambda,lamda_min]);
            else
                lambda = mean([lambda,thresholdValues(ind-1)]);
            end
        end
    end
    
    %% -Each column of the output 'vertices' corresponds to a vertex
    function vertices=getVerticesOfPolyhedron(A,b,lb,ub)
        %% -Incorporate 'lb' & 'ub' in A & b
        nDimensions = size(A,2);
        if(~gf.isempty('lb'))
            if(numel(lb)==1) 
                lb = lb*ones(nDimensions,1); 
            end
            A = [A;-eye(nDimensions)];
            b = [b;-lb];
        end
        if(~gf.isempty('ub'))
            if(numel(ub)==1) 
                ub = ub*ones(nDimensions,1); 
            end
            A = [A;eye(nDimensions)];
            b = [b;ub];
        end
        %% -Now solve the problem
        vertices = (lcon2vert(A,b))';
    end
    
    %% -Wrapper function for finding the point in polyhedra Ax<=b which is closest to any permutation matrix
    function [closestPermutationMatrix,closestPoint] = findClosestPermutationMatrixToConvexHull(A,b,isReturnPointAsMatrix,algorithm)
        if(gf.isempty('isReturnPointAsMatrix')) isReturnPointAsMatrix = false; end
        if(gf.isempty('algorithm')) algorithm = 'squeezeDirecting1'; end
        if(strcmpi(algorithm,'squeezeDirecting1'))
            [closestPermutationMatrix,closestPoint] = sf.findClosestPermutationMatrixToConvexHull_squeezeDirecting1(A,b,isReturnPointAsMatrix);
        elseif(strcmpi(algorithm,'vertexEnumeration'))
            [closestPermutationMatrix,closestPoint] = sf.findClosestPermutationMatrixToConvexHull_vertexEnumeration(A,b,isReturnPointAsMatrix);
        else
           error(['Error(findClosestPermutationMatrixToConvexHull): Specify proper algorithm! Presently  algorithm = ',algorithm]); 
        end
        
    end
    
    %% -Algorithm
    %       min -tr(XX')
    %        A*X(:) <= b
    function [closestPermutationMatrix,closestPoint] = findClosestPermutationMatrixToConvexHull_squeezeDirecting1(A,b,isReturnPointAsMatrix)
        N = sqrt(size(A,2));
        assert(gf.isInteger(N),'Error: N is not an integer');
        I = eye(N);
%         H = -2*kron(ones(N,N),eye(N));                        %% sum(X'AY) = xT*kron(ones(m,n),A)*y, x=X(:), y=Y(:) and m=#Cols(X), n=#Cols(Y) ["Matrix-Cook book" -Eq.524]  
        H = -2*eye(N^2);                                        %% trace(X'X) = xT*I*x, x=X(:),
%         hessianFcn = @(x,lambda) H;
%         closestPoint = quadprog(H,'',A,b);
%         closestPoint = fmincon(@(x) x.'*H*x/2,I(:),A,b);      %% -Division by 2, is done to compensate the earlier scaling by '2' of H  
        closestPoint = fmincon(@(x) x.'*H*x/2,zeros(N^2,1),A,b,[],[],[],[],[], ...      %% -Division by 2, is done to compensate the earlier scaling by '2' of H  %%TODO: How to find a better starting point.
            optimoptions('fmincon','Display','none','Algorithm','interior-point' ));
%         , 'GradObj','on','GradConstr','on','Hessian','user-supplied','HessFcn',@(x,lambda) sf.quadhess(x,lambda,H)));     %% Using this function, throws error "Error using  / Too many output arguments."  !!!     
        closestPermutationMatrix = sf.findClosestPermutationMatrixToPoint(closestPoint);
        if(isReturnPointAsMatrix)
            closestPoint = reshape(closestPoint,N,N);
        end
    end
    
    function hess = quadhess(x,lambda,H)
        hess = H;
    end
    
    %% -Each column of 'vertices' is a vertex of the convex hull
    % Algorithm: 
    %   Step1: Find the closest permutation matrix to the set of vertices 
    %   Step2: Find the closest point in the convex hull to the chosen permutation matrix in step1 
    % Conclusion: This program won't work for dimensions > 25 !!!
    function [closestPermutationMatrix,closestPoint] = findClosestPermutationMatrixToConvexHull_vertexEnumeration(A,b,isReturnPointAsMatrix)
        vertices = sf.getVerticesOfPolyhedron(A,b);
        if(isempty(vertices))       %% -When convex hull is empty set
            closestPermutationMatrix = [];
            closestPoint = [];
            return;
        end
        N = sqrt(size(vertices,1)); assert(gf.isInteger(N),'Error: N is not a integer!!!');
        distance = Inf;
        %% -Find closestPermutationMatrix to the vertices of convex hull 
        for iVertex = 1:size(vertices,2)
            vertex = vertices(:,iVertex);
            X = sf.findClosestPermutationMatrixToPoint(vertex);            
            [distance,index] = gf.minRandIndex([distance,norm(vertex-X(:))]);
            if(index==2) closestPermutationMatrix = X; end
        end
        %% -Now find the closest point in convexHull to closestPermutationMatrix
        C = eye(numel(closestPermutationMatrix));
        closestPoint = lsqlin(C,closestPermutationMatrix(:),A,b);
        if(isReturnPointAsMatrix) closestPoint = reshape(closestPoint,N,N); end
    end    
    
    %% -Find the closest permutation matrix to the point 
    function closestPermutationMatrix = findClosestPermutationMatrixToPoint(point)
        N = sqrt(numel(point));
        I = eye(N);
        closestPointMatrix = reshape(point,N,N);
        [colIndexPerRow,~] = hungarianAlgorithm(-closestPointMatrix);
        closestPermutationMatrix = I(:,colIndexPerRow);
    end
    
    %% Since Ax <= b not feasible, Solve min sum(z) s.t. Ax-z <= b, z >= 0.
    %   then return the indices of blocking constraints after sorting (z>0) in descending order (the most blocking as first indexed).
    %  NOTE: If 'isBoundsSpecified=true', then we ignore (while computing blocking constraints) last rows of 'A' which corresponds to bounding constraints. 
   function blockingConstraints = findBlockingConstraints(A,b,isBoundsSpecified)
        nDim = size(A,2);
        if(isBoundsSpecified)
            nBounds = 2*nDim;
            nConstr = size(A,1) - nBounds;
            A_boundPart = A(nConstr+1:end,:);
            A = A(1:nConstr,:);
        else
            nBounds = 0;
            nConstr = size(A,1);
            A_boundPart = zeros(0,nDim);
        end
        A = [A,eye(nConstr);A_boundPart,zeros(nBounds,nConstr)];
        c = [zeros(1,nDim),ones(1,nConstr)];
        xAndz = linprog(c,A,b);
        z = xAndz(nDim+1:end);
        [~,blockingConstraints] = sort(z,'descend');
        blockingConstraints = blockingConstraints(1:sum(z>0));
    end
    
    
    %% Since Ax <= b not feasible, Solve min t, s.t. Ax-t <= b,
    %   then return the indices of blocking constraints after sorting (z>0) in descending order (the most blocking as first indexed).
    %  NOTE: If 'isBoundsSpecified=true', then we ignore (while computing blocking constraints) last rows of 'A' which corresponds to bounding constraints. 
    function [blockingConstraints,x,t] = findBlockingConstraints_t(A,b,isBoundsSpecified,lb,ub)
        nDim = size(A,2); tol = 1e-10;
        [A,b] = sf.normaliseConstraints(A,b);
        if(gf.isempty('lb')) lb = []; elseif(numel(lb)==1) lb = [lb*ones(nDim,1);-Inf]; end
        if(gf.isempty('ub')) ub = []; elseif(numel(ub)==1) ub = [ub*ones(nDim,1);Inf]; end        
        if(isBoundsSpecified)
            nBounds = 2*nDim;
            nConstr = size(A,1) - nBounds;
            A_boundPart = A(nConstr+1:end,:);
            A = A(1:nConstr,:);
        else
            nBounds = 0;
            nConstr = size(A,1);
            A_boundPart = zeros(0,nDim);
        end
        A = [A,-ones(nConstr,1);A_boundPart,zeros(nBounds,1)];
        c = [zeros(1,nDim),1];
        xAndt = sf.linprogGurobi(c,A,b,[],[],lb,ub,optimoptions('linprog','Algorithm','dual-simplex','Display','none'));
        x = xAndt(1:end-1); t = xAndt(end);
        [~,blockingConstraints] = sort(b-A*xAndt);
        blockingConstraints = intersect(blockingConstraints,find(A*xAndt>=b-tol));
%         blockingConstraints = gf.randPermute(blockingConstraints);
    end
    
    
    %% 
    % NOTE: 'weakConstraints' and 'ignoreConstraints' must be logical vectors 
    %       but the output 'constraints' must be an integer vector containing indices of the constraints 
    function [constraints,xStar,exitStatus] = sortConstraintsInLP(A,b,weakConstraints,mandatoryConstraints,method)
        exitStatus = 0;
        if(gf.isempty('ignoreConstraints')) ignoreConstraints = ~(weakConstraints | mandatoryConstraints); end
        c = [zeros(size(A,2),1);1];
        if(strcmpi(method,'BestConstraintSlack'))       %% Finding Strongest link
            A = [-A,-weakConstraints];
            b = -b;
        elseif(strcmpi(method,'WorstConstraintSlack'))  %% Finding Weakest link
            A = [A,-weakConstraints];           
        else
            error(['Error(sortConstraintsInLP); Unkown method = ',method]);
        end
        x0 = [zeros(size(A,2)-1,1);100];    %% Without this, linprog is throwing infeasibility error sometimes!
        xStar = linprog(c,A(~ignoreConstraints,:),b(~ignoreConstraints),[],[],[],[],x0,optimoptions('linprog','Display','none','Algorithm','dual-simplex','MaxIterations',1e5));
        if(~isempty(xStar))
            [~,constraints] = sort(b-A*xStar);
        else
           constraints = randperm(numel(b))';       %% TODO: May be this is what leads to removing 
           exitStatus = 2^9;     %% Setting 10th bit
        end
        constraints = intersect(constraints,find(weakConstraints),'stable'); %% 'stable' is placed to maintain the order
    end

    %% TODO_Urgent: Change AXB >= gamma*Y  into  AXB >= gamma*Y
    function [isFeasible,closestPermutationMatrix,Xcap] = isYfeasible(A,B,gamma,Y)
        assert(numel(size(A))==2 && size(A,1)==size(A,2),'Error: A must be a square matrix!!');
        % -Convert AXB problem into Dx<=b problem
        [~,D,b] = AXBproblem_to_Dxb_translate(A,B,gamma,Y,2);   
        [closestPermutationMatrix,Xcap]= findClosestPermutationMatrixToPolytope(D,b,true);  %%sf.findClosestPermutationMatrixToConvexHull(A,b,true,'squeezeDirecting1'); %% 
        if(~isempty(Xcap) && norm(closestPermutationMatrix(:)-Xcap(:)) < 10^-10)
            isFeasible = true;
        else
            isFeasible = false;
        end
    end    
    
    
% %     %% -STEPS:
% %     %       1. Translate B=(I-gammaT*aci)'  ->   A=(I-(gammaT/1+kappa*gammaT)(aci+kappa*I))  
% %     %                    gamma              ->   (gammaT/1+kappa*gammaT)*noiseVariance
% %     %       2. Translate problem AXB >=gamma*Y   ->   Dx<=b
% %     %       3. Chop off unwanted constraints: D -> D(Y),  b -> b(Y)
% %     %       4. Compute any feasible point (closestPointAsMatrix) in Dx <= b, and map that point to closestPermutationMatrix
% %     %       5. If closestPermutationMatrix is inside polytope Dx <= b, then closestPointAsMatrix = closestPermutationMatrix  
% %     %  -ASSUMPTIONS
% %     %       1. Inputted gamma = gammaT*noiseVariance
% %     function [isFeasible,closestPermutationMatrix,closestPointAsMatrix] = isYfeasible(A,B,gamma,Y)
% %         assert(numel(size(A))==2 && size(A,1)==size(A,2),'Error: A must be a square matrix!!');
% %         global noiseVariance;
% %         gammaT = gamma/noiseVariance;
% %         N = size(A,2);
% %         I = eye(N);
% %         aci = (I-B')/gammaT;
% %         isFeasible = false;
% %         tol1 = 10^-6; tol2 = 10^-5;
% %         kappa_lower = 10^-4;    %%10^-2; 
% %         kappa_upper = 10^8;    %% Assumption: kappa_lower>=1, Otherwise remove the rounding while taking geometric mean!
% %         itr=1; maxItr = 100;
% %         while(itr<=maxItr)
% %             if(itr==1)
% %                kappa = 0;
% %             elseif(itr==2)
% %                kappa = kappa_lower;
% %             else
% %                kappa = sqrt(kappa_upper*kappa_lower);  %% Increasing cappa helps to reach more close to PermutationMatrix, but also increases the chance of infeasibility!                
% %             end
% %             [isFeasible,closestPermutationMatrixNew,closestPointAsMatrix,exitType] = isYfeasibleWithkappa(A,B,gamma,Y,kappa);
% %             if(isFeasible || kappa_upper-kappa_lower<10^-9 || (itr<=2 && exitType==1))
% %                closestPermutationMatrix = closestPermutationMatrixNew;
% %                break;
% %             elseif(itr>2 && exitType==1)     %% Empty Polytope
% %                kappa_upper = kappa;
% %             elseif(itr>2 && exitType==2)     %% Polytope is too large to pinpoint permutation-matrix in it. It's also possible that the Polytope doesn't contain any permutation matrix.
% %                kappa_lower = kappa;
% %                closestPermutationMatrix = closestPermutationMatrixNew;
% %             elseif(itr>2)
% %                error('Error: Unknown exitType!!!');
% %             end
% %             itr = itr+1;
% %         end
% %         
% %         function [isFeasible,closestPermutationMatrix,closestPointAsMatrix,exitType] = isYfeasibleWithkappa(A,B,gamma,Y,kappa)
% %             isFeasible = false;
% %             %% -Step 1
% %             B = (I-(gammaT/(1+kappa*gammaT))*(aci+kappa*I))';
% %             gamma = noiseVariance*gammaT/(1+kappa*gammaT);
% %             %% -Step 2
% %             [~,D,b] = AXBproblem_to_Dxb_translate(A,B,gamma,true);   
% %             %% -Step 3
% %             D = [D(find(Y),:);D(N^2+1:end,:)];      %% Removing non-candidates, and Re-appending Birkhoff Polytope constraints
% %             b = [b(find(Y));b(N^2+1:end)];
% %             %% -Step 4
% % %             [closestPointAsMatrix,~,exitFlag] = linprog([],D,b,[],[],[],[],optimoptions('linprog','Display','none','algorithm','dual-simplex','MaxIterations',1e6));    %% Algorithm is specified to avoid any info-warnings
% %             clear model params; model.A=sparse(D); model.rhs=b; model.obj=zeros(N^2,1); model.sense='<'; params.outputflag = 0;
% %             result = gurobi(model,params);                          %% Using Gurobi, since matlab linprog is saying infeasible in many cases!
% %             exitFlag = sf.mapGurobiStatusToexitFlag(result.status);
% %             if(exitFlag==1)  closestPointAsMatrix = result.x;  end
% %             
% %             if(exitFlag==1)
% %                 closestPointAsMatrix = reshape(closestPointAsMatrix,N,N);
% %                 closestPermutationMatrix = sf.findClosestPermutationMatrixToPoint(closestPointAsMatrix);
% %                 %% -Step 5
% %                 if(all(D*closestPermutationMatrix(:)<=b+tol2))
% %                     isFeasible = true;
% %                     closestPointAsMatrix = closestPermutationMatrix;
% %                     exitType = 0;
% %                 else
% %                     exitType = 2;           %% How come   gf.all(A*B<=gamma)=1 ?? for itr 101/25              
% %                 end
% %             else                                    %% For empty Polytope
% %                 closestPermutationMatrix = eye(N);
% %                 closestPointAsMatrix = [];
% %                 exitType = 1;
% %             end  
% %             if(exitFlag==0)
% %                disp('Warning(isYfeasibleWithkappa): MaxIteration Exhausted !!!') 
% %             end
% %         end
% %     end
    
    %% Convert Gurobi status string to  Matlab linprog exitFlag 
    %   Gurobi Status messages: http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html#sec:StatusCodes  
    %   Matlab exitFlags      : https://se.mathworks.com/help/optim/ug/linprog.html#outputarg_exitflag  
    function exitFlag = mapGurobiStatusToexitFlag(gurobiStatus)
        switch(gurobiStatus)
            case 'OPTIMAL'
                exitFlag = 1;
            case 'INFEASIBLE'
                exitFlag = -2;
            case 'ITERATION_LIMIT'
                exitFlag = 0;
            otherwise
                exitFlag = 100;
        end
    end
    
    %% Inputs must be of same format as MATLAB 'linprogs' 
    function x = linprogGurobi(c,A,b,~,~,lb,ub,~)
        model.A = sparse(A);
        model.obj = c;
        model.rhs = b;
        model.sense = '<';
        model.vtype = 'C';
        model.modelsense = 'min';
        params.TimeLimit = 300;
        params.IntFeasTol = 1e-9;
        params.outputflag = 0;
        if(~gf.isempty('lb')) model.lb = lb; end
        if(~gf.isempty('ub')) model.ub = ub; end
        result = gurobi(model, params);
        x = result.x;    
    end
    

    function [D,b] = normaliseConstraints(D,b)
        normalisingFactor = nanmean(abs(D),2);
        D = bsxfun(@rdivide,D,normalisingFactor);
        b = b./normalisingFactor;
    end
    
    
    
    
    end
end
