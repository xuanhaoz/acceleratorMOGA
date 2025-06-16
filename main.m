function out = main(varargin)
    verbose = getoption(varargin,'verbose',1);
    debug = getoption(varargin,'debug',0);
    getFuncs = getoption(varargin,'getFuncs',0);

    if getFuncs
        out.getObjective = @getObjective;
        out.getParam = @getParam;
        out.applyChange = @applyChange;
        out.outputFunction = @outputFunction;
        out.objectiveFunction = @objectiveFunction;
        return
    end

    global gaHistory
    
    % ring = AS2v333;
    % version = 331;
    ring = AS2v625_sol58;
    version = 621;
    % ring = AS2v622;
    % version = 621;
    % ring = AS2v335;
    % version = 331;
    
    % load custom initial ring from previous run
    %
    % ring = load('sol_20250213.mat').ring;

    % paramList = {'OC1';'OC2';'OC3';'OC4'};
    paramList = {'OC1';'OC2';'OC3';'OC4';
        'SF2';'SF3';'SF4';'SF5';'SD2';'SD3';'SD4';'SD5'};

    lb1 = ones(1,length(paramList));
    ub1 = ones(1,length(paramList));

    % octupole constraints
    %
    lb1(1:4) = -1.6e5;
    ub1(1:4) =  1.6e5;

    % sextupole limits
    %
    lb1(5:8)  = 0;
    ub1(5:8)  = 1300;
    lb1(9:12) = -1300;
    ub1(9:12) = 0;

    % custom theta distribution for more even coverage
    %
    thetas = [0 0.05 0.1 0.3 0.5 0.7 0.9 0.95 1]*pi;

    x0 = [];
    for n = 1:length(paramList)
        x0(n) = getParam(ring,paramList{n});
    end

    chromRM = getProjectedSextupoleStrength(ring,x0(5:12)','getRM',1);

    DAlinear = geometricAcceptance(ring);
    binarySearchDA_funcs = binarySearchDA(ring,'getFuncs',1);
    R0 = zeros(1,length(thetas));
    for n = 1:length(thetas)
        R0(n) = binarySearchDA_funcs.findStartingRadius(thetas(n),DAlinear);
    end

    objList = {
        'DA0', 0, 1;
        'DAdp+', 0, 1;
        'DAdp-', 0, 1;
        'unstableAx', 0, 1;
        'unstableAy', 0, 1;
        'unstableDppos', 0, 1;
        'unstableDpneg', 0, 1;
        'DAnegx', 0, 1;
    };

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    TolX = 1E-6;
    TolFun = 1E-9;
    MaxFunEvals = 10000;
    MaxIter = 1000;
    walltime = 18;      % hours
    MaxTime = walltime*60*60*0.98;  % 95% of requested wall time
    PopulationSize = 200;
    FunctionTolerance = 1e-5;
    nTurns = 128;

    if debug
        MaxIter = 3;
        PopulationSize = 100;
        nTurns = 1;
    end

    InitialPopulationMatrix = [];
    InitialPopulationMatrix(1,:) = x0;
    for i = 2:fix(PopulationSize/2)
        InitialPopulationMatrix(i,:) = (1 + rand * 0.1)*x0;
    end

    fobj = @(x) objectiveFunction(ring,paramList,objList,x,...
                                    'thetas',thetas,...
                                    'nTurns',nTurns,...
                                    'chromRM',chromRM,...
                                    'version',version,...
                                    'R0',R0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple fminsearch
    %
    % options = optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,...
    %                     'Display','iter','PlotFcns',{@optimplotx @optimplotfval});
    % [sol,fval] = fminsearch(fun,x0,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options = optimoptions('gamultiobj','UseParallel',true,'Display','diagnose',...
                            'MaxGenerations',MaxIter,'PopulationSize',PopulationSize,'MaxTime',MaxTime,...
                            'FunctionTolerance',FunctionTolerance,'UseParallel',true,...
                            'InitialPopulationMatrix',InitialPopulationMatrix,...
                            'OutputFcn',@outputFunction);
                            % 'PlotFcns',{@gaplotbestf @gaplotbestindiv @gaplotrange},...

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single objective genetic algorithm
    %
    % [sol,fval,exitflag,output,population,scores] = ga(fobj,4,[],[],[],[],lb1,ub1,@getConstraint,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiobj GA
    %
    [sol,fval,exitflag,output,population,scores] = gamultiobj(fobj,length(paramList),[],[],[],[],lb1,ub1,options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = [];
    for i = 1:size(sol,1)
        m = sol(i,5:12)'; % free varying sextupole strengths
        out = getProjectedSextupoleStrength(ring,m,'in',chromRM);
        t(i,:) = out.t;
    end
    
    tFams = {'SF1','SD1'};
    paramList = {paramList{:}, tFams{:}}';

    nVars = size(sol,2);
    sol(:,nVars+1:nVars+length(tFams)) = t;

    ring = applyChange(ring,paramList,sol);

    out.ring = ring;
    out.sol = sol;
    out.fval = fval;
    out.output = output;
    out.population = population;
    out.scores = scores;
    out.gaHistory = gaHistory;

    if verbose
        fprintf('-----------------------------------------\n')
        fprintf('Final fval: %.4e\n',fval);

        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        % need to change <getObjective> function
        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % sprintf('Final objectives:\n');
        % T = table;
        % DA = 0;
        % for n = 1:nObj
        %     [objName,target,weight] = objList{n,:};
        %     res = getObjective(ring,objName,target,weight,'DA',DA,'thetas',thetas);
        %     if length(res) > 1
        %         objVal = res{1};
        %         DA = res{2};
        %     else 
        %         objVal = res;
        %     end
        %     T = [T;{objName,objVal}];
        % end
        % T.Properties.VariableNames = {'Objective','Loss value'};
        % disp(T);

        sprintf('Final setpoints:\n');
        T = table;
        for n = 1:length(paramList)
            T = [T;{paramList{n},getParam(ring,paramList{n})}];
            % fprintf('%s: %.4f\n',paramList{n},getParam(dsms,uc,paramList{n}));
        end
        T.Properties.VariableNames = {'Param','Value'};
        disp(T);
    end
end

function [state,options,optchanged] = outputFunction(options,state,flag)
    persistent history
    optchanged = false;
    global gaHistory

    switch flag
        case 'init'
            res.population = state.Population;
            res.score = state.Score;
            history{1} = res;
        case 'iter'
            res.population = state.Population;
            res.score = state.Score;
            history{length(history)+1} = res;
        case 'done'
            res.population = state.Population;
            res.score = state.Score;
            history{length(history)+1} = res;
    end
    gaHistory = history;

end

function val = getParam(varargin)
    ring    = varargin{1};
    param   = varargin{2};

    switch param(1:2)
        case {'OC'}
            ord = atgetcells(ring,'FamName',param);
            val = atgetfieldvalues(ring,ord,'PolynomB',{4});
            val = val(1);
        case {'SF','SD'}
            ord = atgetcells(ring,'FamName',param);
            val = atgetfieldvalues(ring,ord,'PolynomB',{3});
            val = val(1);
        otherwise
            error('Unsupported param')
    end
end


function val = getObjective(varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deprecated
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ring    = varargin{1};
    objName = varargin{2};
    target  = varargin{3};
    weight  = varargin{4};
    npoints = getoption(varargin,'npoints',20);
    dp      = getoption(varargin,'dp',0.04);

    DA = getoption(varargin,'DA',0);

    switch objName
        case {'DA0','DApos','DAneg'}
            switch objName
                case 'DA0'
                    DA = binarySearchDA(ring,'epsilon',1e-5,'verbose',0,varargin{:});
                case 'DApos'
                    DA = binarySearchDA(ring,'dp',dp,'epsilon',1e-5,'verbose',0,varargin{:});
                case 'DAneg'
                    DA = binarySearchDA(ring,'dp',-dp,'epsilon',1e-5,'verbose',0,varargin{:});
            end
            % val{1} = weight * (1/DA.area - target);
            f = DA.linearDA - DA.RMAXs;
            f = max(zeros(length(f),1), f);
            f = f./DA.linearDA;

            % adding extra weights to DA horizontal DA
            %
            weights = ones(length(DA.thetas),1);
            idx = find(DA.thetas <= 0.05 | DA.thetas >= 0.95*pi);
            weights(idx) = weights(idx)*10;
            f = f.*weights;

            val{1} = sum(f)/sum(weights);
            val{2} = DA;
        case 'ADTS'
            if ~isstruct(DA)
                DA = binarySearchDA(ring,'epsilon',1e-5,'verbose',0,'thetas',[0 pi/2],varargin{:});
            end
            % find lines to calculate tune footprint
            %
            idx = find(DA.thetas == 0 | DA.thetas == pi/2);
            ADTS = [];
            N_outside = [];
            for i = idx
                res = computeADTS(ring,DA.thetas(i),DA.RMAXs(i),'npoints',npoints);
                ADTS(end+1) = res.distance;
                N_outside(end+1) = npoints - res.N_viable;
            end
            val = sum(ADTS) + sum(N_outside);
        otherwise
            error('Unsupported objName')
    end
end

function ring = applyChange(varargin)
    ring    = varargin{1};
    paramList   = varargin{2};
    valList     = varargin{3};

    for n = 1:length(paramList)
        param = paramList{n};
        val = valList(n);
        switch param(1:2)
            case {'OC'}
                ord = atgetcells(ring,'FamName',param);
                ring = atsetfieldvalues(ring,ord,'PolynomB',{4},val);
            case {'SF','SD'}
                ord = atgetcells(ring,'FamName',param);
                ring = atsetfieldvalues(ring,ord,'PolynomB',{3},val);
            otherwise
                error('Unsupported param')
        end
    end
end

function F = objectiveFunction(varargin)
    ring        = varargin{1};
    paramList   = varargin{2};
    objList     = varargin{3};
    sol         = varargin{4};
    chromRM     = getoption(varargin,'chromRM',0);
    npoints     = getoption(varargin,'npoints',20);
    dp          = getoption(varargin,'dp',0.03);
    maxdp       = getoption(varargin,'dp',0.05);
    version     = getoption(varargin,'version',331);

    maxk2 = 1300;   % sextupole field strength limit [m-2]

    % find chromcaticity correction sextupole strengths 
    %
    m = sol(5:12)'; % free varying sextupole strengths
    out = getProjectedSextupoleStrength(ring,m,'in',chromRM);
    t = out.t;

    tFams = {'SF1','SD1'};
    for n = 1:length(tFams)
        paramList{end+1} = tFams{n};
        sol(end+1) = t(n);
    end

    ring = applyChange(ring,paramList,sol);

    % objective function
    %
    F = zeros(length(objList),1);

    % check if chromaticity sextupoles strength constraint
    %
    k2_penalty = max(0,abs(t)-maxk2);
    if any(k2_penalty > 0)
        F = F + 2;
        F(1:length(k2_penalty)) = F(1:length(k2_penalty)) + k2_penalty;

        return
    end

    % compute DAs
    %
    DA0   = binarySearchDA(ring,'dp',0, 'epsilon',1e-5,'verbose',0,varargin{:});
    DApos = binarySearchDA(ring,'dp',dp, 'epsilon',1e-5,'verbose',0,varargin{:});
    DAneg = binarySearchDA(ring,'dp',-dp, 'epsilon',1e-5,'verbose',0,varargin{:});

    % get DA size along theta = pi
    %
    idx = find(DA0.thetas == pi);
    DAnegx = DA.linearDA(idx) - DA0.RMAXs(idx);
    DAnegx = DAnegx/DA.linearDA(idx);

    % compute ADTS
    %
    idx = find(DA0.thetas == pi | DA0.thetas == pi/2);
    ADTS = [];
    unstableADTS = [];
    for i = 1:length(idx)
        res = computeADTS(ring,DA0.thetas(idx(i)),DA0.RMAXs(idx(i)),'npoints',npoints,'version',version);
        ADTS(i) = res.distance;
        unstableADTS(i) = (npoints - res.N_viable)/npoints;
    end

    % compute MDTS    
    %
    res = computeMDTS(ring,'version',version,'maxdp',maxdp,'npoints',npoints);
    MDTS(1) = res.distance;
    unstableMDTS(1) = (npoints - res.N_viable)/npoints;

    res = computeMDTS(ring,'version',version,'maxdp',-maxdp,'npoints',npoints);
    MDTS(2) = res.distance;
    unstableMDTS(2) = (npoints - res.N_viable)/npoints;

    % construct objective function
    %
    function objectiveValue = calculateDAobjectiveValue(DA)
        f = DA.linearDA - DA.RMAXs;
        f = max(zeros(length(f),1), f);
        f = f./DA.linearDA;

        weights = ones(length(DA.thetas),1);
        idx = find(DA.thetas <= 0.05 | DA.thetas >= 0.95*pi);
        weights(idx) = weights(idx)*10;
        f = f.*weights;

        objectiveValue = sum(f)/sum(weights);
    end

    F(1) = calculateDAobjectiveValue(DA0);
    F(2) = calculateDAobjectiveValue(DApos);
    F(3) = calculateDAobjectiveValue(DAneg);

    F(4) = unstableADTS(1);
    F(5) = unstableADTS(2);

    F(6) = unstableMDTS(1);
    F(7) = unstableMDTS(2);

    F(8) = DAnegx;

    F = F + 0.01*sum(MDTS) + 1*sum(ADTS);

end