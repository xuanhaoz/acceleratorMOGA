function out = matchDSMS(varargin)
    verbose = getoption(varargin,'verbose',1);
    debug = getoption(varargin,'debug',0);
    getFuncs = getoption(varargin,'getFuncs',0);

    if getFuncs
        out.getParam = @getParam;
        out.applyChange = @applyChange;
        out.outputFunction = @outputFunction;
        out.objectiveFunction = @objectiveFunction;
        return
    end
    
    % master input option list
    %
    opts = {};

    uc0 = AS2v336('line','UC0');
    dsms = AS2v336('line','DSMS');
    opts.uc0 = uc0;
    opts.dsms = dsms;

    [ed,rd] = atlinopt(uc0,0,1:length(uc0)+1);
    opts.twiss_in = ed(end);

    paramList = {'QMS1','QMS2','QMS3','QMS4'};
    opts.paramList = paramList;

    lb1 = ones(1,length(paramList));
    ub1 = ones(1,length(paramList));

    % quadrupole constraints
    %
    lb1(1:4) = -14;
    ub1(1:4) =  14;

    x0 = [];
    for n = 1:length(paramList)
        x0(n) = getParam(paramList{n},opts);
    end

    objList = {
        'betax/y', 1.8, 1;
        'betay/x', 1.8, 1;
        'betax=y', 0.05, 1;
    };
    opts.objList = objList;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    TolX = 1E-6;
    TolFun = 1E-9;
    MaxFunEvals = 10000;
    MaxIter = 1000;
    walltime = 15;      % hours
    MaxTime = walltime*60*60*0.98;  % 95% of requested wall time
    PopulationSize = 200;
    FunctionTolerance = 1e-5;
    nTurns = 128;

    if debug
        MaxIter = 3;
        PopulationSize = 4;
        nTurns = 1;
    end

    InitialPopulationMatrix = [];
    % InitialPopulationMatrix(1,:) = x0;
    % for i = 2:fix(PopulationSize/2)
    %     InitialPopulationMatrix(i,:) = (1 + rand * 0.1)*x0;
    % end

    fobj = @(x) objectiveFunction(x,opts);

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

    ring = applyChange(sol,opts);

    out.ring = ring;
    out.sol = sol;
    out.fval = fval;
    out.output = output;
    out.population = population;
    out.scores = scores;

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
            T = [T;{paramList{n},getParam(paramList{n},opts)}];
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
    param   = varargin{1};
    opts    = varargin{2};
    ring    = opts.dsms;

    switch param(1:3)
        case {'QMS'}
            ord = atgetcells(ring,'FamName',param);
            val = atgetfieldvalues(ring,ord,'PolynomB',{2});
            val = val(1);
        otherwise
            error('Unsupported param')
    end
end

function ring = applyChange(varargin)
    valList     = varargin{1};
    opts        = varargin{2};
    ring        = opts.dsms;
    paramList   = opts.paramList;

    for n = 1:length(paramList)
        param = paramList{n};
        val = valList(n);
        switch param(1:3)
            case {'QMS'}
                ord = atgetcells(ring,'FamName',param);
                ring = atsetfieldvalues(ring,ord,'PolynomB',{2},val);
            otherwise
                error('Unsupported param')
        end
    end
end

function F = objectiveFunction(varargin)
    sol         = varargin{1};
    opts        = varargin{2};
    dsms        = opts.dsms;
    paramList   = opts.paramList;
    objList     = opts.objList;

    dsms = applyChange(sol,opts);

    % objective function
    %
    F = zeros(length(objList),1);

    oc1 = find(atgetcells(dsms,'FamName','OC1'));
    oc2 = find(atgetcells(dsms,'FamName','OC2'));
    oc3 = find(atgetcells(dsms,'FamName','OC3'));

    ords = [oc1 oc2 oc3];

    [ed,rd] = atlinopt(dsms,0,1:length(dsms)+1,'twiss_in',opts.twiss_in);

    % calcualte final matching alpha
    %
    alpha = ed(end).alpha;
    if any(abs(alpha) > 1e-6)
        F(1) = 2 + abs(alpha(1));
        F(2) = 2 + abs(alpha(2));
        F(3) = 2;
        return
    end

    % calculate beta fraction at octupole locations
    %
    beta = cat(1,ed.beta);

    beta_xdivy = beta(ords,1)./beta(ords,2);
    beta_ydivx = beta(ords,2)./beta(ords,1);

    [~,id_xdivy] = max(beta_xdivy);
    [~,id_ydivx] = min(beta_xdivy);

    id_xeqy = [1 2 3];
    id_xeqy([id_xdivy id_ydivx]) = [];

    F(1) = beta_ydivx(id_xdivy);
    F(2) = beta_xdivy(id_ydivx);
    F(3) = abs(1 - beta_xdivy(id_xeqy));

    F = F + 1000 * sum(abs(alpha));

end