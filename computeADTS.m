function out = computeADTS(varargin)
    % function to compute the amplitude dependent tune shift along search vectors
    %
    ring     = varargin{1};
    theta    = varargin{2};
    Rmax     = varargin{3};
    npoints  = getoption(varargin,'npoints',20);
    makeplot = getoption(varargin,'plot',0);
    dp       = getoption(varargin,'dp',0);
    nturns   = getoption(varargin,'nturns',128);
    version  = getoption(varargin,'version',621);
    parallel = getoption(varargin,'parallel',0);

    if parallel
        parforArg = Inf;
    else
        parforArg = 0;
    end

    % define half integer triangle
    %

    % fractional tune only
    %
    % a = [0.05 0.05];
    % b = [0.45 0.45];
    % c = [0.45 0.05];

    switch version
        case 621
            a = [61.05 19.05];
            b = [61.45 19.05];
            c = [61.45 19.45];
        case 331
            a = [62.05 20.05];
            b = [62.45 20.05];
            c = [62.45 20.45];
        case 333
            a = [62.05 20.05];
            b = [62.45 20.05];
            c = [62.45 20.45];
        case 225
            a = [62.55 19.55];
            b = [62.95 19.55];
            c = [62.95 19.95];
        otherwise
            error('Lattice version not supported');
    end

    T = [a;b;c];

    % get cartesian coordiantes for tracking points
    %
    r = linspace(3e-5,Rmax,npoints);
    thetas = ones(1,length(r))*theta;
    [x,y] = pol2cart(thetas,r);

    parfor (i = 1:length(x), parforArg)
        tunes(i,:) = trackTune(ring,x(i),y(i),varargin{:});
    end

    % ----------------------------
    % find correct fractional tune according to <atnuampl> function
    %
    [ld,frac0] = atlinopt(ring,dp,1:length(ring)+1);
    tune0 = ld(end).mu/2/pi;
    offs = [1 -1];

    [~,k] = min([frac0 - tunes(1,:); 1 - frac0 - tunes(1,:)]);
    np = offs(k);
    offset = round(tune0 - np.* tunes(1,:));
    tunes = np.*tunes + offset;

    out.rawTunes = tunes;

    tunes = fixTuneTurnAround(tunes);

    % check points in allowed triangle
    %
    [in,on] = inpolygon(tunes(:,1),tunes(:,2),T(:,1),T(:,2));
    N_viable = length(find(in));

    distance = 0;
    if any(~in)
        idx = find(~in);
        for i = 1:length(idx)
            distance(i) = findDistanceToTriangle(tunes(i,:),T);
        end
        distance = sum(distance.^2);
    end

    if makeplot
        figure(6021)
        T = [a;b;c;a];
        ax_handle = gca;
        lines = get(ax_handle,'Children');
        tags = [];
        if length(lines) > 0
            tags = [lines.Tag];
        end

        if ~strcmp(tags,'triangle')
            plot(T(:,1),T(:,2),'LineWidth',2,'Tag','triangle');
        end
        hold on
        plot(tunes(in,1),tunes(in,2),'b+');
        plot(tunes(~in,1),tunes(~in,2),'ro');
        plot(tunes(:,1),tunes(:,2),'DisplayName',sprintf('theta = %.1f deg',rad2deg(theta)),'LineWidth',1.5);
        xlim([min(T(:,1))-0.05 max(T(:,1))+0.05]);
        ylim([min(T(:,2))-0.05 max(T(:,2))+0.05]);
        xlabel('Qx'); ylabel('Qy')
        grid on
    end

    out.tunes = tunes;
    out.distance = distance;
    out.N_viable = N_viable;
end


function tune = trackTune(ring,x,y,varargin)
    % wrapper for calling AT <findtune> function
    %
    dp = getoption(varargin,'dp',0);
    nturns = getoption(varargin,'nturns',128);

    epsilon = 1e-4;
    x = fsign(x)*max(abs(x),epsilon);
    y = fsign(y)*max(abs(y),epsilon);

    p0 = [x;0;y;0;dp;0];
    p1 = atpass(ring,p0,1,nturns);

    % clunky way to suppress fprintf from findtune
    %
    x = p1(1,:)';
    y = p1(3,:)';

    [~,vx] = evalc('findtune(x)');
    [~,vy] = evalc('findtune(y)');

    % tune = [findtune(p1(1,:)') findtune(p1(3,:)')];
    tune = [vx vy];
end

function out = fsign(in)
    if in == 0
        out = 1;
    else
        out = sign(in);
    end
end

function dist = findDistanceToTriangle(tune,T)
    % calculate nearest distance between point and triangle outter boundary
    %
    a = T(1,:);
    b = T(2,:);
    c = T(3,:);

    dist(1) = pointSegmentDistance(tune,a,b);
    dist(2) = pointSegmentDistance(tune,a,c);
    dist(3) = pointSegmentDistance(tune,b,c);

    dist = min(dist);
end

function d = pointSegmentDistance(P,A,B)
    % calculate smallest distance between point P and segment AB
    %
    AB = B-A;
    AP = P-A;
    t = dot(AP,AB)/dot(AB,AB);
    t = max(0,min(1,t));
    nearestPoint = A + t*AB;
    d = norm(P-nearestPoint);
end
