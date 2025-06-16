function out = plotSol(varargin)
    in = varargin{1};
    sol = varargin{2};
    ring = varargin{3};
    version = getoption(varargin,'version',621);

    paramList = {'OC1';'OC2';'OC3';'OC4';
        'SF2';'SF3';'SF4';'SF5';'SD2';'SD3';'SD4';'SD5';
        'SF1';'SD1'};

    f = main('getFuncs',1);

    ring = f.applyChange(ring,paramList,in.sol(sol,:));
    out.ring = ring;

    thetas = [0 0.05 0.1 0.3 0.5 0.7 0.9 0.95 1]*pi;

    DA0 = binarySearchDA(ring,'thetas',thetas,'plot',1,'label',sprintf('Sol. %d',sol),'epsilon',1e-5);
    out.DA0 = DA0;
    
    idx = find(thetas == pi);
    out.ADTSx = computeADTS(ring,DA0.thetas(idx),DA0.RMAXs(idx),'plot',1,'version',version,'parallel',1);
    idx = find(thetas == 0.5*pi);
    out.ADTSy = computeADTS(ring,DA0.thetas(idx),DA0.RMAXs(idx),'plot',1,'version',version,'parallel',1);

    out.MDTSpos = computeMDTS(ring,'plot',1,'maxdp',0.05,'version',version);
    out.MDTSneg = computeMDTS(ring,'plot',1,'maxdp',-0.05,'version',version);
end