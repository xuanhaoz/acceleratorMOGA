function chromRM = getChromRM(ring,varargin)
    % calculate the chromaticity response matrix
    %
    SXfams = getoption(varargin,'fams',0);
    deltak2 = getoption(varargin,'dk2',0.1);

    if ~iscell(SXfams)
        SXfams = {
            'SF1','SF2','SF3','SF4','SF5',...
            'SD1','SD2','SD3','SD4','SD5'};
    end

    nFams = length(SXfams);
    chromRM = zeros(2,nFams);

    [rd,~] = atlinopt6(ring,'get_chrom');
    oldChrom = rd.chromaticity(1:2);

    for n = 1:nFams
        sx = atgetcells(ring,'FamName',SXfams{n});
        oldk2 = atgetfieldvalues(ring,sx,'PolynomB',{3});
        ring1 = atsetfieldvalues(ring,sx,'PolynomB',{3},oldk2(1)+deltak2);

        [rd,~] = atlinopt6(ring1,'get_chrom');
        newChrom = rd.chromaticity(1:2);
        dChrom = newChrom - oldChrom;
        chromRM(:,n) = dChrom/deltak2;
    end
end