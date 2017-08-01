%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Find the max norm point in the polytope Ax<=b
%       ie.  max |x|^2  subject to: Ax<=b 
%
% -Assumption
%   1. The polytope Ax<=b fully lies within Birkhoff Polytope. 
%
% -EXAMPLES
%   1. To test accuracy of this program, Comment replacing closestPoint with closestPermutationMatrix paert, and do following, 
%        n=4^2; A=[eye(n);-eye(n)]; b=[ones(n,1);zeros(n,1)]; [perMatrix,point]=findClosestPermutationMatrixToPolytope(A,b); norm(point)
%
% -OBSERVATION & CONCLUSION
%   1. For 100 dimensional box constraints (2016/10/12):
%       6.7 seconds with maxItr=10, solution is 6.16 (~=10)
%       75 seconds with maxItr=100, solution is 9.36 (~=10)
%      Conclusion: CAN'T USE THIS CODE, since 
%           1. it's taking too high time, 
%           2. Throwing error 'Error in findMaxVolumeEllipsoidInPolytope>mve_presolve (line 180)  R = chol(B);'  
%
% -Debug
%   1. Use plotRegion() while debugging
%
% -TODO
%   1. 
%
% -REFERENCE
%   1. http://web.stanford.edu/class/cme334/docs/2011-11-17-Rubira_normmax.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [closestPermutationMatrix,closestPoint] = findClosestPermutationMatrixToPolytope(A,b,isReturnPointAsMatrix)

warningState = warning;
warning('off','all');

    [m,n] = size(A);
    [xStar,~,exitFlag] = linprog([zeros(n,1);1],[A,-ones(m,1)],b,[],[],[],[],[],optimoptions('linprog','Display','none'));  %% min t, s.t Ax-t <= b
    xStar(end) = [];

    Anew = A; bnew = b;
    maxItr = 100; itr = 1;
    while(itr<maxItr)
        [xStarNew,Adelta,bdelta,exitStatus] = findMaxNormPointInPolytope(Anew,bnew);
        Anew = [Anew;Adelta];  bnew = [bnew;bdelta];
        if(exitStatus || (norm(xStarNew)-norm(xStar))/norm(xStarNew) < 1e-3)         %% 'ellipseVolume' is not a good stopping criteria, since ellipsoid can have zero volume when it lies on a plane (ie. subspace)       
            %Anew(end,:) = []; bnew(end,:) = [];   %% This is to handle the case, when ellipseVolume<10^-16, function 'findMaxVolumeEllipsoidInPolytope' may return errored xStarNew
            break;
        end
        xStar = xStarNew;
        itr = itr +1;
    end
    closestPoint = xStar;
    closestPermutationMatrix = sf.findClosestPermutationMatrixToPoint(closestPoint);
    %%  -If closestPermutationMatrix lies in Ax<=b 
    if(all(A*closestPermutationMatrix(:)<=b+eps))
        closestPoint = closestPermutationMatrix(:);
    end
    
    if(~gf.isempty('isReturnPointAsMatrix') && isReturnPointAsMatrix)
        closestPoint = reshape(closestPoint,sqrt(numel(xStar)),[]);
    end

warning(warningState);    

end





function [xStar,Adelta,bdelta,exitStatus] = findMaxNormPointInPolytope(A,b)

    exitStatus = 0; Q = nan(size(A,2));
    xStar = []; Adelta = []; bdelta = [];

    %% 1st stage: Compute max-volume Ellipsoid Qu+r  s.t. |u|?1
    try
        [r,Q,exitStatus] = findMaxVolumeEllipsoidInPolytope(A,b);   %% TODO: Why Q is not symmetric?? (hence not positive definite?)
    catch
        exitStatus = exitStatus + 2^14;
        return;
    end
    n = size(Q,1);  %% n = N^2 (ie #VUE^2)
    ellipseVolume = (abs(det(Q)))^(1/n);
    if(ellipseVolume<10^-9) exitStatus = exitStatus + 2; end
    if(~all(A*r<=b)) exitStatus = exitStatus + 4; end           %% r not lies in the feasible region! 
    if(exitStatus)
        return;
    end

    %% 2nd stage find max-norm point on the computed ellipse
    Qinv = Q^-1;
%     Qinvs = Qinv^2;
    Qinvs = (Q*Q')^-1;
    lambda = sdpvar(1);
    gamma = sdpvar(1);
    I = eye(size(Q));
    
%     constraints = [lambda >= 0,[-I+lambda*Qinvs, -lambda*Qinvs'*r; -lambda*r'*Qinvs, lambda*r'*Qinvs*r-lambda-gamma] >= 0];
    %% Refer Stephan-Boyd: 4.6.2
    F{1} = zeros(n+1); F{1}(1:n,1:n) = -eye(n);
    F{2} = [Qinvs,-Qinvs'*r;-r'*Qinvs,r'*Qinvs*r-1];
    F{3} = zeros(size(F{1})); F{3}(end,end) = -1;
%     constraints = [lambda >= 0,F{1} + lambda*F{2} + gamma*F{3} >=0];              
%     diagonostic = optimize(constraints,-gamma,sdpsettings('verbose',0));    %% TODO: Why this command is failing, even for sufficiently large ellipsoid?? (ie 16 dimensional box constraints problem)
    [temp, info] = solveSDPstandard(F,[0,-1]',[-1,0],0);
    lambda = temp(1);
    xStar = pinv(-I+lambda*Qinvs)*lambda*Qinvs*r;    
    if(~all(A*xStar<=b))        %% xStar not lies in the feasible region!
        exitStatus = exitStatus + 8;
    end
%     if(diagonostic.problem==0)
%         lambda = value(lambda);
%         xStar = pinv(-I+lambda*Qinvs)*lambda*Qinvs*r;
%     else                                        %% If yalmip couldn't solve the problem, then go along the best eigen vector direction
%         [eigVectors,eigValues] = eig(Q);
%         Xcandidates = bsxfun(@plus,r,[eigVectors*eigValues,-eigVectors*eigValues]);
%         [~,XcandidatesIndex] = max(gf.norm(Xcandidates));
%         xStar = Xcandidates(:,XcandidatesIndex);
%     end
%     disp(['norm(xStar) = ',num2str(norm(xStar))]);
    %% 3rd stage: Create Excluding halfspace constraint (Complement of supporting halfspace)
    Adelta = -(2*xStar'*Qinvs-2*r'*Qinvs);    
%     Adelta = -xStar';         %% Observation: Adelta = -xStar' improves solution??  
    bdelta = Adelta*xStar;
end

