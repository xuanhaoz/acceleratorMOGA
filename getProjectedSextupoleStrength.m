function out = getProjectedSextupoleStrength(varargin)
    % calculate the chromaticity correction sextupole strengths
    % from free sextupole parameters
    % using constant tuning matrices
    % X1 = M*m + T*t + X0, where:
    %   X1: target chromaticity
    %   X0: natural chromaticity
    %   m: free varying sextupole strengths
    %   t: chromaticity correction sextupole strengths
    %   M: response matrix for m
    %   T: response matrix for t
    %
    ring = varargin{1};
    m = varargin{2};
    getRM = getoption(varargin,'getRM',0);
    X1 = getoption(varargin,'X1',[1;1]);
    in = getoption(varargin,'in',0);

    X0 = 0;
    M = 0;
    T = 0;
    if isstruct(in)
        X0 = in.X0;
        M  = in.M;
        T  = in.T;
    end

    tFams = {'SF1','SD1'};
    mFams = {'SF2','SF3','SF4','SF5',...
             'SD2','SD3','SD4','SD5'};
            
    if ~(length(X0) == 2)
        sx = find(atgetcells(ring,'Class','Sextupole'));
        ring0 = atsetfieldvalues(ring,sx,'PolynomB',{3},0);

        [rd,~] = atlinopt6(ring0,'get_chrom');
        X0 = rd.chromaticity(1:2)';
        out.X0 = X0;
    end

    if getRM 
        sx = find(atgetcells(ring,'Class','Sextupole'));
        ring0 = atsetfieldvalues(ring,sx,'PolynomB',{3},0);

        M = getChromRM(ring0,'fams',mFams,'dK2',10);
        T = getChromRM(ring0,'fams',tFams,'dK2',10);
        out.M = M;
        out.T = T;
    else
        assert(size(T,2) == length(tFams),'T matrix does not match tFam size');
        assert(size(M,2) == length(mFams),'M matrix does not match tFam size');
    end

    A = X1 - M*m - X0;
    t = pinv(T)*A;
    out.t = t;
    
end