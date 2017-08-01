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
%   
%
%%%% -DETAILED DESCRIPTION:
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1. 
%
%%%% -NOTES:
%
%%%% -NOTES (Programming):
%
%
%%%% -EXECUTION:
%       1. fetchResult.majorAlgorithms:  Takes 20min for fetching six ".mat" files
%
%%%% -TODO:
%       1. 
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




classdef fetchResult

methods(Static)

    
    function data = majorAlgorithms(fetchConfig,foldernameOrFiles)
        fetchConfig = fetchResult.patchFetchConfig('fetchConfig');
        files = fetchResult.getFiles(foldernameOrFiles);

        data = fetchResult.majorAlgorithms_(fetchConfig,files);
    end
    
    

    function data = majorAlgorithms_(fetchConfig,files)

    % -Make sure un-assigned values in 'data' are NaNs
    data.dataPerDuplexTypePerACIRtype = cell(2,2);
    for ii = 1:numel(data.dataPerDuplexTypePerACIRtype)
        data.dataPerDuplexTypePerACIRtype{ii}.nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = matrixGeneralized([],[],[],nan);
        data.dataPerDuplexTypePerACIRtype{ii}.txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = matrixGeneralized([],[],[],nan);
    end

    %% -Process each files
    for ifile = 1:numel(files)
        out = load(files{ifile});
        opt = out.opt;
        if(~gf.isempty('out.config.configChannel.isEnableFullDuplex'))      %% Legacy support!
            out.config.isEnableFullDuplex = out.config.configChannel.isEnableFullDuplex;
        end    
        if(numel(opt.powerControlAlgorithmsToRun)==1)
            opt.powerControlAlgorithmsToRun = repmat(opt.powerControlAlgorithmsToRun,1,numel(opt.schAlgorithmsToRun));
        end
        % -Compute nSuccess
        nSuccessPerFsPerTsPerNsPerAlgorithm = permute(nanmean(gf.nansum(gf.congregate('out.out{:,:,:}(:).isSuccessPerTxPerRxPerAlgorithm')==1,[5,6]),4),[1,2,3,7,4,5,6]);   %% 'gf.nansum' & 'nanmean' are placed to accomodate the case when scheduling algorithm("sch_gurobi.m") failed to schedule any VUEs
        nSuccessPerFsPerTsPerNsPerAlgorithm = bsxfun(@rdivide,nSuccessPerFsPerTsPerNsPerAlgorithm,shiftdim(opt.Ns,-2));     %% Normalising to #successfulLinksPerVUE
        % -Compute txPower (total power of a VUE across FxT frequency-time grid) averaged across all Tx-VUEs
        % (Obsolete Code): Correct 'powerPerTxPerTimePerAlgorithm' using 'txVUEPerRBPerTimePerAlgorithm'
        % %     temp1 = gf.congregate('out.out{:,:,:}(:).powerPerTxPerTimePerAlgorithm');
        % %     temp2 = gf.congregate('out.out{:,:,:}(:).txVUEPerRBPerTimePerAlgorithm');
        % %     temp1 = gf.preserveDataOnSpecifiedIndices(temp1,temp2,nan,5);
        % %     txPowerPerFsPerTsPerNsPerAlgorithm = permute(nanmean(gf.nansum(temp1,[5,6]),4),[1,2,3,7,4,5,6]);   %% 'gf.nansum' & 'nanmean' are places to accomodate the "Scheduling failure" case of sch_gurobi.m
        txPowerPerFsPerTsPerNsPerAlgorithm = permute(nanmean(gf.nansum(gf.congregate('out.out{:,:,:}(:).powerPerTxPerTimePerAlgorithm'),[5,6]),4),[1,2,3,7,4,5,6]);   %% 'gf.nansum' & 'nanmean' are places to accomodate the "Scheduling failure" case of sch_gurobi.m
        txPowerPerFsPerTsPerNsPerAlgorithm = bsxfun(@rdivide,txPowerPerFsPerTsPerNsPerAlgorithm,shiftdim(opt.Ns,-2));     %% Normalising to #successfulLinksPerVUE
        % -Assigning copmuted values to 'dataPerDuplexTypePerACIRtype'
        isEnableFullDuplex = out.config.isEnableFullDuplex;
        iDuplex = 1 + isEnableFullDuplex;               %% 1 -> halfDuplex,  2 -> fullDuplex
        acirTypeIndex = acir.acirIndexPerType(out.config.acirType);
        for iAlgorithm = 1:numel(opt.schAlgorithmsToRun)
            schIndex = fetchConfig.indexPerSchAlgorithm(opt.schAlgorithmsToRun{iAlgorithm});
            powerControlIndex = fetchConfig.indexPerPowerControlAlgorithm(opt.powerControlAlgorithmsToRun{iAlgorithm});
            data.dataPerDuplexTypePerACIRtype{iDuplex,acirTypeIndex}.nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(opt.Fs,opt.Ts,opt.Ns,schIndex,powerControlIndex) ...
                = nSuccessPerFsPerTsPerNsPerAlgorithm(:,:,:,iAlgorithm);
            data.dataPerDuplexTypePerACIRtype{iDuplex,acirTypeIndex}.txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(opt.Fs,opt.Ts,opt.Ns,schIndex,powerControlIndex) ...
                = txPowerPerFsPerTsPerNsPerAlgorithm(:,:,:,iAlgorithm);
        end
    end
    
    for ii = 1:numel(data.dataPerDuplexTypePerACIRtype)
        data.dataPerDuplexTypePerACIRtype{ii}.nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = double(data.dataPerDuplexTypePerACIRtype{ii}.nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl);
        data.dataPerDuplexTypePerACIRtype{ii}.txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = double(data.dataPerDuplexTypePerACIRtype{ii}.txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl);
    end

    end
    

    %% -
    function files = getFiles(foldernameOrFiles)
        if(isa(foldernameOrFiles,'cell'))
            files = foldernameOrFiles;
        elseif(exist(foldernameOrFiles,'dir')==7)
            foldername = foldernameOrFiles;
            temp1=dir([foldername,'*_architecture2_wrapper.mat']);
            files = {temp1.name}';
            files = gf.findregexpNot(files,'.*error.*');
        elseif(isa(foldernameOrFiles,'char') && exist(foldernameOrFiles,'file'))
            files{1} = foldernameOrFiles;
        else
            error('Error: Please input nonempty ''foldernameOrFiles'' ...');
        end
    end
        
    
    %% Patch 'fetchConfig'
    % Warning: 'fetchConfig.indexPerSchAlgorithm' must be matching to the 'plotConfig.schAlgInTitle' order in plotResult.m file    
    function fetchConfig = patchFetchConfig(fetchConfigString)
        if(evalin('caller',['gf.isempty(''',fetchConfigString,''')']))
            fetchConfig = struct();
        else
            fetchConfig = evalin('caller',fetchConfigString);
        end;
%         fetchConfigDefault.indexPerSchAlgorithm = containers.Map({'sch_naturalScheduling','sch_blockInterleaver1','sch_blockInterleaver2','sch_blockInterleaver3','sch_Heuristic1','sch_Heuristic1([],[],[],[],true)', ...
%             'sch_gurobi','sch_gurobi('' '','' '',txVUEPerRBPerTime_0)','sch_gurobi(config,topology,txVUEPerRBPerTime_0)','sch_optimizedBlockInterleaver1'},{1,2,3,4,5,6,7,7,7,8});
%         fetchConfigDefault.indexPerPowerControlAlgorithm = containers.Map({'powerControl_equalPower','powerControl_kangWang','powerControl_gurobi'},{1,2,3});
        fetchConfigDefault.indexPerSchAlgorithm = Co.indexPerSchAlgorithm;
        fetchConfigDefault.indexPerPowerControlAlgorithm = Co.indexPerPowerControlAlgorithm;        
        fetchConfig = gf.mergeTwoStructures(fetchConfigDefault,fetchConfig,'',true);
    end



end
end


