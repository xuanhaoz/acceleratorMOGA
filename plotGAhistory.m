function out = plotGAhistory(in,varargin)
    plotSol = getoption(varargin,'plotSol',0);
    smoothPlot = getoption(varargin,'smoothPlot',1);

    nGens = length(in.gaHistory);
    [popsize,nObjs] = size(in.scores);
    figNum = 1234;

    scores = zeros(nGens,popsize,nObjs);
    for n = 1:nGens
        scores(n,:,:) = in.gaHistory{n}.score;
    end

    plotVars = {};
    plotVars.mean = mean(scores,2);
    % plotVars.rms = rms(scores,2);
    plotVars.min = max(min(scores,[],2),1e-10);

    plotNames = fieldnames(plotVars);
    figure(figNum)
    clf;
    for n = 1:length(plotNames)
        subplot(length(plotNames),1,n);
        plotVar(plotVars.(plotNames{n}),plotNames{n},figNum+n);
    end
    fig=gcf;
    fig.Position(3:4) = [600,600];
    set(fig,'PaperPosition',[0 0 6 6]);
    set(fig,'PaperSize',[6 6]);
    movegui('center');

    % figure(figNum+1)
    % clf;
    % [popsize,nObjs] = size(in.fval);
    % for n = 1:3
    %     subplot(3,1,n);
    %     plot(1:popsize,in.fval(:,n),'LineWidth',2);
    %     ylabel(sprintf('Obj. %d',n));
    %     grid on
    % end

    % figure(figNum+2)
    % clf;
    % for n = 1:4
    %     subplot(4,1,n);
    %     plot(1:popsize,in.fval(:,n+3),'LineWidth',2);
    %     ylabel(sprintf('Obj. %d',n+3));
    %     grid on
    % end

    figure(figNum+3)
    clf;
    if smoothPlot
        xq1 = 1:0.01:nObjs;
        for i = 1:length(in.fval(:,1))
            s(i,:) = pchip(1:nObjs,in.fval(i,:),xq1);
        end
        plot(xq1,s,'Color','#6b6b6b');
    else
        plot(1:nObjs,in.fval,'Color','#6b6b6b');
    end

    lines = [];
    if plotSol
        cmap = parula(5);
        colororder('default');
        hold on
        for n = 1:length(plotSol)
            sol = plotSol(n);
            if smoothPlot
                line = plot(xq1,s(sol,:),'LineWidth',2,'DisplayName',sprintf('Sol. %d',sol));
            else
                line = plot(1:nObjs,in.fval(sol,:),'LineWidth',2,'DisplayName',sprintf('Sol. %d',sol));
            end
            lines(n) = line;
        end
        legend(lines)
    end
    grid on

    xticklabels({'DA0','DA+','DA-','ADTSx','ADTSy','MDTS+','MDTS-','DAnegx'});
    ylabel('Obj. value')

end

function plotVar(varargin)
    fval = varargin{1};
    fname = varargin{2};
    
    [nGens,~,nObjs] = size(fval);

    for n = 1:nObjs
        plot(1:nGens,fval(:,1,n),'LineWidth',2,'DisplayName',sprintf('Obj. %d',n));
        hold on
    end
    grid on
    legend();
    yscale('log');
    xlabel('Generations')
    ylabel(fname)
end
