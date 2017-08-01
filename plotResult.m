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
%       RESULT ANALSIS OVER MULTIPLE ALGORITHMS
%
%%%% -DETAILED DESCRIPTION:
%       1. Output 'isSuccessPerRBPerTimePerRxVUEPerAlgorithm' contains NaN, which indicates that no VUE transmitted on that frequency-time slot 
%
%%%% -OPTIONS:
%
%
%%%% -ACRONYMS
%       1. 'pc'  =>  Power Control
%       2. 'Alg' =>  Algorithm
%
%%%% -EXAMPLES:
%       1. plotResult.plot('/Volumes/Glenn/Programming/Matlab/ALLtoALL_Communication/Architecture2/2016121827220_architecture2_wrapper.mat')
%          plotResult.plot('~/Programming/Matlab/ALLtoALL_Communication/Architecture2/Results/2016Dec24/2016122052796_architecture2_wrapper.mat')
%%%% -NOTES:
% 
%
%%% -Execution
%       1. Fetching time: 'data = plotResult.fetchData(fetchConfig,files(1));'  takes 24*1.5 = 36minutes
%
%
%
%%%% -NOTES (Programming):
%       1. 'txVUE' = 0  =>  No VUE is transmitting
%       2. 'successStatus' = 1      =>  Success
%                          = 0      =>  Failure
%                          = NaN    =>  No attempted transmission or receptions
%
%%%% -WARNINGS (Programming):
%       1. Matlab won't distinguish between a structure and, "an array of structers with numel '1' " 
%       2. In this program, never try to change Default fetchConfig/plotConfig/optConfig
%
%%%% -TODO:
%       1. TODO_Urgent: Fix bug in matrixGeneralized
%           a= matrixGeneralized(magic(3)); a(:,:)=123
%           b=magic(3); b(:,:)=123
%       2. TODO_Urgent: Correct 'powerPerTxPerTimePerAlgorithm' using 'txVUEPerRBPerTimePerAlgorithm'
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef plotResult

methods(Static)

%%
% NOTE: One can also input 'nSuccessPerFsPerTsPerNsPerSchAlgorithm' instead of 'nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl', 
%           then opt.powerControlAlgorithmsToRun = 'powerControl_equalPower' is assumed
function plotnSuccessfulLinksAndAvgTxPower(plotConfig,opt,nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl,txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl,figNumbersToPlot)

    if(gf.isempty('figNumbersToPlot'))  figNumbersToPlot = [1:1000]; end

    plotConfig = plotResult.patchPlotConfig('plotConfig');

    opt.Fs = sort(opt.Fs,'descend');
    opt.Ts = sort(opt.Ts,'descend');
    opt.Ns = sort(opt.Ns,'descend');
    opt.PAs = sort(opt.PAs,'descend');       %% Power Control algorithm indices
    plotConfig.pcAlgInLegend = flipud(plotConfig.pcAlgInLegend);
    % nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(opt.Fs,2:end,end,end,:) = nan;      %% Setting NaN for all T>0 for 'alg_Gurobi'   %% Update(2017/02/24): Fixed bug in matrixGeneralized.m
    nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(:,:,:,plotConfig.schAlgIndicesToPlot,:);
    txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(:,:,:,plotConfig.schAlgIndicesToPlot,:);

    
    if(~gf.isempty('nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl'))
    %% Plot nSuccessfulLinks for "Equal-Powercontrol" 
        nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = gf.augmentMatrix(nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl,[max(opt.Fs),max(opt.Ts),max(opt.Ns),nan,max(opt.PAs)],nan,true);
        nSuccessPerFsPerTsPerNsPerSchAlgorithm = nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(:,:,:,:,1);
        %% 1. Fix N=20. Plot #Successful links/VUE  Vs F, each curve for each T
        figNumber = 1; 
        if(ismember(figNumber,figNumbersToPlot))
            figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
            gf.plot(opt.Fs,nSuccessPerFsPerTsPerNsPerSchAlgorithm(opt.Fs,opt.Ts,end,:),[1,2,4],[],plotConfig.lineStylePerSchAlgIndex);
            title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of RBs (F)'); ylabel('Number of successful Links per VUE'); h.XTick = sort(opt.Fs); grid on; 
            h = legend(num2str(opt.Ts),'Location','southeast'); v = get(h,'title'); set(v,'string','Number of time slots');
            plotResult.saveFigure(plotConfig,figNumber);
        end
        
        %% 2. Fix N=20. Plot #Successful links/VUE  Vs T, each curve for each F
        figNumber = 2; 
        if(ismember(figNumber,figNumbersToPlot))
            figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
            gf.plot(opt.Ts,nSuccessPerFsPerTsPerNsPerSchAlgorithm(opt.Fs,opt.Ts,end,:),[2,1,4],[],plotConfig.lineStylePerSchAlgIndex);
            title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of timeslots (T)'); ylabel('Number of successful Links per VUE'); h.XTick = sort(opt.Ts); grid on; 
            h = legend(num2str(opt.Fs),'Location','southeast'); v = get(h,'title'); set(v,'string','Number of RBs');
            plotResult.saveFigure(plotConfig,figNumber);
        end

        %% 3. Fix F=20. Plot #Successful links/VUE  Vs N, each curve for each T
        figNumber = 3;
        if(ismember(figNumber,figNumbersToPlot))
        figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
            gf.plot(opt.Ns,nSuccessPerFsPerTsPerNsPerSchAlgorithm(end,opt.Ts,opt.Ns,:),[3,2,4],[],plotConfig.lineStylePerSchAlgIndex);
            title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of VUEs (N)'); ylabel('Number of successful Links per VUE'); h.XTick = sort(opt.Ns); grid on; 
            h = legend(num2str(opt.Ts),'Location','southeast'); v = get(h,'title'); set(v,'string','Number of time slots');
            plotResult.saveFigure(plotConfig,figNumber);
        end

        %% 4. Fix T=1. Plot #Successful links/VUE  Vs N, each curve for each F
        figNumber = 4;
        if(ismember(figNumber,figNumbersToPlot))
            figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
            gf.plot(opt.Ns,nSuccessPerFsPerTsPerNsPerSchAlgorithm(opt.Fs,1,opt.Ns,:),[3,1,4],[],plotConfig.lineStylePerSchAlgIndex);
            title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of VUEs (N)'); ylabel('Number of successful Links per VUE'); h.XTick = sort(opt.Ns); grid on; 
            h = legend(num2str(opt.Fs),'Location','northwest'); v = get(h,'title'); set(v,'string','Number of RBs');
            plotResult.saveFigure(plotConfig,figNumber);
        end

    %% Plot nSuccessfulLinks w.r.t "Power-Control"
        %% 5. Fix T=1, N=20. Plot #Successful links/VUE  Vs F, each curve for each "Power-control Algorithms"
        if(size(nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl,5)>1 && gf.any(nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(:,:,:,:,2:end)>0))
            figNumber = 5;
            if(ismember(figNumber,figNumbersToPlot))
                figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
                gf.plot(opt.Fs,nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(opt.Fs,1,end,:,opt.PAs),[1,5,4],[],plotConfig.lineStylePerSchAlgIndex);
                title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of RBs (F)'); ylabel('Number of successful Links per VUE'); h.XTick = sort(opt.Fs); grid on; 
                h = legend(plotConfig.pcAlgInLegend,'Location','southeast'); v = get(h,'title'); set(v,'string','Power control algorithms');
                plotResult.saveFigure(plotConfig,figNumber);
            end

            %% 6. Fix F=20, N=20. Plot #Successful links/VUE  Vs F, each curve for each "Power-control Algorithms"
            figNumber = 6;
            if(ismember(figNumber,figNumbersToPlot))
                figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
                gf.plot(opt.Ts,nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(end,opt.Ts,end,:,opt.PAs),[2,5,4],[],plotConfig.lineStylePerSchAlgIndex);
                title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of timeslots (T)'); ylabel('Number of successful Links per VUE'); h.XTick = sort(opt.Ts); grid on; 
                h = legend(plotConfig.pcAlgInLegend,'Location','southeast'); v = get(h,'title'); set(v,'string','Power control algorithms');
                plotResult.saveFigure(plotConfig,figNumber);
            end

            %% 7. Fix F=20, T=1,  Plot #Successful links/VUE  Vs F, each curve for each "Power-control Algorithms"
            figNumber = 7;
            if(ismember(figNumber,figNumbersToPlot))
                figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
                gf.plot(opt.Ns,nSuccessPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(end,1,opt.Ns,:,opt.PAs),[3,5,4],[],plotConfig.lineStylePerSchAlgIndex);
                title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of VUEs (N)'); ylabel('Number of successful Links per VUE'); h.XTick = sort(opt.Ns); grid on; 
                h = legend(plotConfig.pcAlgInLegend,'Location','southeast'); v = get(h,'title'); set(v,'string','Power control algorithms');
                plotResult.saveFigure(plotConfig,figNumber);
            end
        end
    end   
    
    
    %% Plot "Average Tx-Power" w.r.t "Power-Control"
    if(~gf.isempty('txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl'))
        
        txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = gf.augmentMatrix(txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl,[max(opt.Fs),max(opt.Ts),max(opt.Ns),nan,max(opt.PAs)],nan,true);
        
        % -Coverting power to dBm
        txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl = gf.pow2db(txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl)+30;
        
        %% 8. Fix N. Plot #Successful links/VUE  Vs F, each curve for each T
        figNumber = 8;
        if(ismember(figNumber,figNumbersToPlot))
            figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
            gf.plot(opt.Fs,txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(opt.Fs,1,end,:,opt.PAs),[1,5,4],[],plotConfig.lineStylePerSchAlgIndex);
            title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of RBs (F)'); ylabel('Average tx power (dBm)'); h.XTick = sort(opt.Fs); grid on;
            h = legend(plotConfig.pcAlgInLegend,'Location','southeast'); v = get(h,'title'); set(v,'string','Power control algorithms');
            plotResult.saveFigure(plotConfig,figNumber);
        end

        %% 9. Fix N. Plot #Successful links/VUE  Vs F, each curve for each T
        figNumber = 9;
        if(ismember(figNumber,figNumbersToPlot))
            figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
            gf.plot(opt.Ts,txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(end,opt.Ts,end,:,opt.PAs),[2,5,4],[],plotConfig.lineStylePerSchAlgIndex);
            title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of timeslots (T)'); ylabel('Average tx power (dBm)'); h.XTick = sort(opt.Ts); grid on;
            h = legend(plotConfig.pcAlgInLegend,'Location','northwest'); v = get(h,'title'); set(v,'string','Power control algorithms');
            plotResult.saveFigure(plotConfig,figNumber);
        end

        %% 10. Fix N. Plot #Successful links/VUE  Vs F, each curve for each T
        figNumber = 10;
        if(ismember(figNumber,figNumbersToPlot))
            figure('Name',[num2str(figNumber),'. ',plotConfig.figureHeading]); h = gca;
            gf.plot(opt.Ns,txPowerPerFsPerTsPerNsPerSchAlgorithmPerPowerControl(end,1,opt.Ns,:,opt.PAs),[3,5,4],[],plotConfig.lineStylePerSchAlgIndex);
            title(plotConfig.titleString,'Interpreter','latex'); xlabel('Number of VUEs (N)'); ylabel('Average tx power (dBm)'); h.XTick = sort(opt.Ns); grid on;
            h = legend(plotConfig.pcAlgInLegend,'Location','southeast'); v = get(h,'title'); set(v,'string','Power control algorithms');
            plotResult.saveFigure(plotConfig,figNumber);
        end
    end    
    

    
end

%%
% NOTE: indices of 'plotConfigDefault.schAlgInTitle' must be of the same order as 'fetchConfig.indexPerSchAlgorithm' in analyzeResults.majorAlgorithms()
function plotConfig = patchPlotConfig(plotConfigString)
    if(evalin('caller',['gf.isempty(''',plotConfigString,''')']))
        plotConfig = struct();
    else
        plotConfig = evalin('caller',plotConfigString);
    end;

    plotConfigDefault.foldername = './';
    plotConfigDefault.schAlgInTitle = {'sch_naturalScheduling','sch_interleaver1','sch_interleaver2','sch_interleaver3','sch_Heuristic1','sch_Heuristic2','sch_gurobi','sch_optimizedBlockInterleaver'};
    plotConfigDefault.pcAlgInLegend = {'Equal-Power','Heuristic','Optimal'}';
    plotConfigDefault.schAlgIndicesToPlot = [1:numel(plotConfigDefault.schAlgInTitle)]';
    plotConfigDefault.lineStylePerSchAlgIndex = {'-','--','--+','-^','-*','-x','-v','-^'};
    plotConfigDefault.lineStyleInTitlePerSchAlgIndex = {'-','- -','--+','-$\Delta$','-*','-x','-$\nabla$','-$\Delta$'};   %% NOTE: we actually ran out of symbols, after 6 plots!
    plotConfigDefault.figureHeading = '';
    plotConfigDefault.isSaveFile = false;
    
    plotConfig = gf.mergeTwoStructures(plotConfigDefault,plotConfig,'',true);
    plotConfig.schAlgInTitle = plotConfig.schAlgInTitle(plotConfig.schAlgIndicesToPlot);
    plotConfig.lineStylePerSchAlgIndex = plotConfig.lineStylePerSchAlgIndex(plotConfig.schAlgIndicesToPlot);
    plotConfig.lineStyleInTitlePerSchAlgIndex = plotConfig.lineStyleInTitlePerSchAlgIndex(plotConfig.schAlgIndicesToPlot);
    
    titleString = '';
    for i = 1:numel(plotConfig.schAlgInTitle)
        titleString = [titleString,gf.translateToLatexString(plotConfig.schAlgInTitle{i}),': ''',plotConfig.lineStyleInTitlePerSchAlgIndex{i},''',\quad   '];
        if(rem(i,3)==0) titleString = [titleString,10]; end
    end
    plotConfig.titleString = titleString;
    
    if(exist(plotConfig.foldername,'dir')~=7)
        mkdir(plotConfig.foldername);
    end
    
end


function nSuccessfulLinksForBestAlgPerFsPerTsPerNs = sch_interleavers_widths(plotConfig,opt,nSuccessPerFsPerTsPerNsPerAlgorithm)

    opt.Fs = sort(opt.Fs,'descend'); opt.Ts = sort(opt.Ts,'descend'); opt.Ns = sort(opt.Ns,'descend');
    
    %% Plot Best Algorithm index
    [nSuccessfulLinksForBestAlgPerFsPerTsPerNs,bestAlgorithmIndexPerFsPerTsPerNs] = max(nSuccessPerFsPerTsPerNsPerAlgorithm,[],4);
%     percentageImprovementForBestAlgPerFsPerTsPerNs = (nSuccessfulLinksForBestAlgPerFsPerTsPerNs./nSuccessPerFsPerTsPerNsPerAlgorithm(:,:,:,3) -1)*100;

%     figure; gf.plot(opt.Ns,percentageImprovementForBestAlgPerFsPerTsPerNs(opt.Fs,opt.Ts,opt.Ns),[3]); h =gca;
%     titleString = ['T=1 ''-'' ,   ','T=2 ''- -'' ,  ','T=3 ''--+'' ,   ','T=4 ''-$\Delta$'' ,   ','T=5 ''-*''    '];       %   ,'T=6 ''-$\nabla$''    '
%     title(titleString,'Interpreter','latex'); xlabel('Number VUEs (N)'); ylabel('Percentage improvement of best interleaver w.r.t width = 3'); h.XTick = sort(opt.Ns); grid on; 
%     h = legend(num2str(opt.Fs),'Location','southeast'); v = get(h,'title'); set(v,'string','Number of RBs');

    figure; gf.plot(opt.Ns,bestAlgorithmIndexPerFsPerTsPerNs(opt.Fs,opt.Ts,opt.Ns),[3]); h =gca;
    titleString = ['T=1 ''-'' ,   ','T=2 ''- -'' ,  ','T=3 ''--+'' ,   ','T=4 ''-$\Delta$'' ,   ','T=5 ''-*''    '];       %   ,'T=6 ''-$\nabla$''    '
    title(titleString,'Interpreter','latex'); xlabel('Number VUEs (N)'); ylabel('Best interleaver width'); h.XTick = sort(opt.Ns); h.YTick = 1:gf.nanmax(bestAlgorithmIndexPerFsPerTsPerNs); grid on; 
    h = legend(num2str(opt.Fs),'Location','southeast'); v = get(h,'title'); set(v,'string','Number of RBs');
    plotResult.saveFigure(setfield(plotConfig,'figureHeading','bestInterleaverWidth'));

    %% Standard plotting with architecture2_plot()
    plotConfig.schAlgInTitle = {'Block interleaver1, width = 1','width = 2','width = 3','width = 4','width = 5','width = 6'};
    plotResult.plotnSuccessfulLinksAndAvgTxPower(plotConfig,opt,nSuccessPerFsPerTsPerNsPerAlgorithm);

end

    
function saveFigure(plotConfig,figNumber)
    if(gf.isempty('figNumber')) figNumber = 0; end
    if(plotConfig.isSaveFile)
        figFileName = [plotConfig.foldername,plotConfig.figureHeading,'_',num2str(figNumber)]; 
        savefig(figFileName); 
        saveas(gcf,figFileName,'epsc');
    end
end

        
end 
end
