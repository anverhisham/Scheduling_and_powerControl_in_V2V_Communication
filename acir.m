%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Supporting function for linkRemovalAlgorithms.m
%
%
% -WARNINGS
%   1. Make sure that diagonal elements of 'acirPerDesiredTxPerInterferer' is zero, otherwise none of the scheduling algorithms will work
%
%
% -LOG:
%   1. 2016/11/24: Changed ' acirPerRBGap = [1;acirPerRBGap]; '     to      '   acirPerRBGap = [0;acirPerRBGap];  '  
%   2. 2016/11/24: Changed dummy ACIR   'acirPerRBGap = db2pow(-30:-10:-10000)';'   to    'acirPerRBGap = db2pow(-30:-5:-10000)';' 
%
% -\Author: Anver <anver@chalmers.se>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef acir

    properties (Constant)
        acirIndexPerType = gf.Map('LTEmask','scFDMA','OFDMA','measuredByJessica');

    end    
    
    
    methods(Static)    
        
    %% - NOTE/WARNING: 1 Gap => Same RB, 2 Gap => Adjacent RB
    function acirPerRBGap = acirPerRBGap(acirType)
        % -Find ACIR values
        switch(acirType)
            case 'LTEmask'
                acirPerRBGap = [db2pow(-30)*ones(4,1); db2pow(-43)*ones(1000,1);]; %% - STANDARD ACIR TABLE (ACIR Values Ref:36.942-5.1.1.3) (Warning: Self ACIR (ie a RB to itself) is kept as '1', Changing it create multiple problems!)
            case 'scFDMA'
                acirPerRBGap = db2pow([-43.0019, -50.0159, -59.4531, -64.6828, -68.7968, -72.1095, -74.5685, -77.1738, -78.6474, -80.7189, -81.4768, -82.9241, -82.9202*ones(1,50)])';     %% (2017/04/11): aci.simulateScenario1(); config.isSCFDMA = true; FFTSize=8192; nScheduledSC=300; percentageClipping=1, 
            case 'OFDMA'
                acirPerRBGap = db2pow([-38.4885, -50.1101, -58.7904, -64.1905, -68.2060, -71.4465, -74.1372, -76.3748, -78.2752, -79.8940, -81.0874, -82.0181*ones(1,50)])';     %% (2017/03/11): aci.simulateScenario1(); FFTSize=8192; nScheduledSC=300; percentageClipping=1,
            case 'measuredByJessica'
                acirPerRBGap = db2pow([-24,-62,-78,-96,-110,-120,-128,-134,-140,-145,-149,-152,   -155,-158:-2:-200])';     %% (2017/03/11): ACIR until -152dB is found out from Jessica's files.
            otherwise
                error(['Specify proper acirType = ',acirType]);
        end
        
        % -Appending with dummy zeros
        if(numel(acirPerRBGap)<1000)
            acirPerRBGap = [acirPerRBGap;zeros(1000-numel(acirPerRBGap),1)];
        end
    end
      
    
    
    function acirPerDesiredTxPerInterferer = getACIRfromInterfererToDesiredLink(txIDPerRB,acirType)

        if(gf.isempty('acirType'))  acirType = 'LTEmask'; end
        acirPerRBGap = acir.acirPerRBGap(acirType);
        acirPerRBGap = [0;acirPerRBGap];        % -NOTE/WARNING: 1 Gap => Same RB, 2 Gap => Adjacent RB (Adopted this way, since 0 index not supported in MATLAB)!! %% TODO_Check: What would happen if prepended by '1' instead of '0'

        rxRBPerRxRBPerTxRB = [1:100]'*ones(1,100);
        txRBPerRxRBPerTxRB = ones(100,1)*[1:100];
        acirPerRxRBPerTxRB = acirPerRBGap(abs(txRBPerRxRBPerTxRB-rxRBPerRxRBPerTxRB)+1);

        %% - INPUT PARAMETER
        N = numel(txIDPerRB);
        [~,rbPerTxID] = sort(txIDPerRB);

        %% - OUTPUT COMPUTATION ...
        acirPerDesiredTxPerInterferer = acirPerRxRBPerTxRB(rbPerTxID,rbPerTxID);
    end
        
    end
end
