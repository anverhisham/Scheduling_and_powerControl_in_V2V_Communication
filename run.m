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
%%%% -TODO:
%       1. 
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef run

methods(Static)

    function out = sch_interleavers_widths()
        opt.powerControlAlgorithmsToRun =  {'powerControl_equalPower'};
        opt.Fs = [6,12,20]';  opt.Ts = [1,2,3,4,5]';  opt.Ns = [3:20]';
        opt.nTrials = 1;

        opt.schAlgorithmsToRun = {'sch_blockInterleaver1([],[],[],[],1)','sch_blockInterleaver1([],[],[],[],2)','sch_blockInterleaver1([],[],[],[],3)','sch_blockInterleaver1([],[],[],[],4)', ...
            'sch_blockInterleaver1([],[],[],[],5)','sch_blockInterleaver1([],[],[],[],6)'};    
        tic; out{1} = architecture2_wrapper([],opt); toc

        opt.schAlgorithmsToRun = {'sch_blockInterleaver2([],[],[],[],1)','sch_blockInterleaver2([],[],[],[],2)','sch_blockInterleaver2([],[],[],[],3)','sch_blockInterleaver2([],[],[],[],4)', ...
            'sch_blockInterleaver2([],[],[],[],5)','sch_blockInterleaver2([],[],[],[],6)'};    
        tic; out{2} = architecture2_wrapper([],opt); toc

        opt.schAlgorithmsToRun = {'sch_blockInterleaver3([],[],[],[],1)','sch_blockInterleaver3([],[],[],[],2)','sch_blockInterleaver3([],[],[],[],3)','sch_blockInterleaver3([],[],[],[],4)', ...
            'sch_blockInterleaver3([],[],[],[],5)','sch_blockInterleaver3([],[],[],[],6)'};    
        tic; out{3} = architecture2_wrapper([],opt); toc
    end

    %% TODO (2017/04/13): why 'powerControl_kangWang' outperforms 'powerControl_gurobi' ??
    function data = test1()
    
        figureHeadingPerDuplexTypePerACIRtype = {'halfDuplexLTEacir','halfDuplexDummyACIR'; 'fullDuplexLTEacir','fullDuplexDummyACIR'};
        opt = struct('nTrials',100,'Ns',[20]','Fs',[6,12,20]','Ts',[1,3,5]', 'schAlgorithmsToRun', {{'sch_naturalScheduling','sch_naturalScheduling','sch_naturalScheduling'}}, ...
                'powerControlAlgorithmsToRun',{{'powerControl_equalPower','powerControl_kangWang','powerControl_gurobi'}});
        for iDuplex = 1:2
            for acirTypeIndex = cell2mat(acir.acirIndexPerType.values({'LTEmask','scFDMA'}))
                isEnableFullDuplex = iDuplex-1;
                acirTypesALL = acir.acirIndexPerType.keys;  acirType = acirTypesALL{acirTypeIndex};
                config = df.v2struct(isEnableFullDuplex,acirType);
                tic, out = architecture2_wrapper(config,opt); toc
                nSuccessPerFsPerTsPerNsPerAlgorithm = permute(nanmean(gf.nansum(gf.congregate('out{:,:,:}(:).isSuccessPerTxPerRxPerAlgorithm')==1,[5,6]),4),[1,2,3,7,4,5,6]);   %% 'gf.nansum' & 'nanmean' are placed to accomodate the case when scheduling algorithm("sch_gurobi.m") failed to schedule any VUEs
                nSuccessPerFsPerTsPerNsPerAlgorithm = bsxfun(@rdivide,nSuccessPerFsPerTsPerNsPerAlgorithm,shiftdim(opt.Ns,-2));
                plotConfig.figureHeading = figureHeadingPerDuplexTypePerACIRtype{iDuplex,acirTypeIndex};
                plotResult.plotnSuccessfulLinksAndAvgTxPower(plotConfig,opt,nSuccessPerFsPerTsPerNsPerAlgorithm);
                data.dataPerDuplexTypePerACIRtype{iDuplex,acirTypeIndex}.nSuccessPerFsPerTsPerNsPerAlgorithm = nSuccessPerFsPerTsPerNsPerAlgorithm;
            end
        end
        
    end
    
    
    


    %% 2017/04/06: Obsolete function (Use architecture2_optimizedBlockInterleaver.m instead), I'm not able to understand this...
    function findBestBlockInterleaverWidth1(isCompute,isPlot)

        if(gf.isempty('isCompute')) isCompute = false; end
        if(gf.isempty('isPlot')) isPlot = false; end
        foldername = 'Datafiles/';

        % -Parameter specification
        opt.schAlgorithmsToRun = {'sch_blockInterleaver1([],[],[],[],1)','sch_blockInterleaver1([],[],[],[],2)','sch_blockInterleaver1([],[],[],[],3)','sch_blockInterleaver1([],[],[],[],4)', ...
            'sch_blockInterleaver1([],[],[],[],5)','sch_blockInterleaver1([],[],[],[],6)','sch_blockInterleaver1([],[],[],[],7)','sch_blockInterleaver1([],[],[],[],8)','sch_blockInterleaver1([],[],[],[],9)'};
        opt.powerControlAlgorithmsToRun =  {'powerControl_equalPower'};
        opt.Fs = [6,12,20]';  opt.Ts = [1:20]';  opt.Ns = [3:20]';
        opt.nTrials = 50;

       if(isCompute) 
            for isEnableFullDuplex = [false,true]
                for acirType = acir.acirIndexPerType.keys
                    acirType = acirType{1};
                    config = df.v2struct(isEnableFullDuplex,acirType);
                    tic; out = architecture2_wrapper(config,opt); toc
                    nSuccessPerFsPerTsPerNsPerAlgorithm = permute(nanmean(gf.nansum(gf.congregate('out{:,:,:}(:).isSuccessPerTxPerRxPerAlgorithm')==1,[5,6]),4),[1,2,3,7,4,5,6]); 
                    nSuccessPerFsPerTsPerNsPerAlgorithm = bsxfun(@rdivide,nSuccessPerFsPerTsPerNsPerAlgorithm,shiftdim(opt.Ns,-2));

                    %% Plot Best Algorithm index
                    [nSuccessfulLinksForBestAlgPerFsPerTsPerNs,bestInterleaverWidthPerFsPerTsPerNs_] = max(nSuccessPerFsPerTsPerNsPerAlgorithm,[],4);
                    percentageImprovementForBestAlgPerFsPerTsPerNs = (nSuccessfulLinksForBestAlgPerFsPerTsPerNs./nSuccessPerFsPerTsPerNsPerAlgorithm(:,:,:,3) -1)*100;
                    bestInterleaverWidthPerFsPerTsPerNsPerDuplexTypePerACIRtype(opt.Fs,opt.Ts,opt.Ns,isEnableFullDuplex+1,acir.acirIndexPerType(acirType)) = bestInterleaverWidthPerFsPerTsPerNs_;
                    
                    % - Result Verifications 
                    if(isEnableFullDuplex==false) % When number of VUEs is less, higher interleaver width shouldn't matter
                        iN = find(opt.Ns==8); iT = find(opt.Ts==3);
                        assert(gf.all(diff(nSuccessPerFsPerTsPerNsPerAlgorithm(:,iT,iN,[3,4,5,6]),1,4)==0),'Error: For lower number of VUEs per timeslots, higher interleaver widths shouldn''t make any difference!!');
                    end
                end
            end
            save([foldername,'sch_optimizedBlockInterleaver1.mat'],'bestInterleaverWidthPerFsPerTsPerNsPerDuplexTypePerACIRtype');
       end
        
        %% Warning: this one, plot only for the last instance of "bestInterleaverWidthPerFsPerTsPerNs_"
        if(isPlot)
            load([foldername,'sch_optimizedBlockInterleaver1.mat'],'bestInterleaverWidthPerFsPerTsPerNsPerDuplexTypePerACIRtype');
            figure; gf.plot(opt.Ns,bestInterleaverWidthPerFsPerTsPerNs_,[3]); h =gca;
            titleString = ['T=1 ''-'' ,   ','T=2 ''- -'' ,  ','T=3 ''--+'' ,   ','T=4 ''-$\Delta$'' ,   ','T=5 ''-*''    '];       %   ,'T=6 ''-$\nabla$''    '
            title(titleString,'Interpreter','latex'); xlabel('Number VUEs (N)'); ylabel('Best interleaver width'); h.XTick = sort(opt.Ns); h.YTick = 1:gf.nanmax(bestInterleaverWidthPerFsPerTsPerNs_); grid on; 
            h = legend(num2str(opt.Fs),'Location','southeast'); v = get(h,'title'); set(v,'string','Number of RBs');

            figure; gf.plot(opt.Ns,percentageImprovementForBestAlgPerFsPerTsPerNs,[3]); h =gca;
            titleString = ['T=1 ''-'' ,   ','T=2 ''- -'' ,  ','T=3 ''--+'' ,   ','T=4 ''-$\Delta$'' ,   ','T=5 ''-*''    '];       %   ,'T=6 ''-$\nabla$''    '
            title(titleString,'Interpreter','latex'); xlabel('Number VUEs (N)'); ylabel('Percentage improvement of best interleaver w.r.t width = 3'); h.XTick = sort(opt.Ns); grid on; 
            h = legend(num2str(opt.Fs),'Location','southeast'); v = get(h,'title'); set(v,'string','Number of RBs');

            %% Standard plotting with architecture2_plot()
            nSuccessPerFsPerTsPerNsPerAlgorithm_ = matrixGeneralized([],'','',nan);
            nSuccessPerFsPerTsPerNsPerAlgorithm_(opt.Fs,opt.Ts,opt.Ns,:) = nSuccessPerFsPerTsPerNsPerAlgorithm;
            plotConfig.figureFolder = foldername;
            plotConfig.schAlgInTitle = {'Block interleaver1, width = 1','width = 2','width = 3','width = 4','width = 5','width = 6'};
            architecture2_plot.plotnSuccessfulLinksAndAvgTxPower(plotConfig,'',nSuccessPerFsPerTsPerNsPerAlgorithm_);  
        end
    
    end 
    
    
    
    

end
end


