function tunes = fixTuneTurnAround(rawTunes,varargin)
    % script to identify and fix any artifical turn around lines in tune space
    % cased by the tracked tunes crossing resonance lines
    %
    makeplot = getoption(varargin,'plot',0);
    WP = getoption(varargin,'WP',[0 0]);

    % remove nans
    %
    idx = isnan(rawTunes);
    idx = any(idx,2);
    rawTunes = rawTunes(~idx,:);

    % identify any turning around points
    %

    % fix fake turning points
    %
    correctedTunes = correctFakeTurningPoint(rawTunes);

    % scale tunes to WP
    %
    tunes = scaleTunes(correctedTunes,WP);

    if makeplot
        figure(7021)
        plot(tunes(:,1),tunes(:,2),'Marker','o','LineWidth',2);
        grid on
    end
end

function correctedTunes = correctFakeTurningPoint(varargin)
    rawTunes    = varargin{1};
    epsilon     = getoption(varargin,'epsilon',0.05);

    idx = diff(sign(diff(rawTunes)));

    correctedTunes = rawTunes;
    if ~any(reshape(idx,1,[]))
        % correctedTunes = rawTunes;
        return
    end

    % get fractional part of tunes
    %
    fracPart = mod(abs(rawTunes),1);
    intPart = fix(rawTunes);

    turningPoints = find(any(idx,2))+1; % +1 to account for diff function index shift

    % flip turning points in reverse order
    %
    for tp = flip(turningPoints')
        axes = find(idx(tp-1,:));
        for axis = axes
            tpTune = fracPart(tp,axis);
            if abs(tpTune) <= 0.5   % tp between 0 ~ 0.5
                if abs(tpTune) < epsilon    % near 0 line
                    correctedTunes(tp+1:end,axis) = intPart(tp+1:end,axis) - fracPart(tp+1:end,axis);
                elseif abs(tpTune - 0.5) < epsilon  % near 1/2 line
                    correctedTunes(tp+1:end,axis) = intPart(tp+1:end,axis) - fracPart(tp+1:end,axis) + 1;
                else    % not classed as fake turning point
                    continue
                end
            else    % tp between 0.5 ~ 1
                if abs(tpTune - 0.5) < epsilon  % near 1/5 line
                    correctedTunes(tp+1:end,axis) = intPart(tp+1:end,axis) - fracPart(tp+1:end,axis) + 1;
                elseif abs(tpTune - 1) < epsilon    % near 1 line
                    correctedTunes(tp+1:end,axis) = intPart(tp+1:end,axis) - fracPart(tp+1:end,axis) + 2;
                else    % not classed as fake turning point
                    continue
                end
            end
        end
    end
end

function tunes = scaleTunes(tunes,WP)
    integerPart = floor(WP);
    integerPart = integerPart .* ones(size(tunes));

    tunes = tunes + integerPart;
end