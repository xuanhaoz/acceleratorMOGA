function out = binarySearchDA(ring,varargin)
    % binary search algorithm for DA
    % from M. Kranjcevic, et. al., 2020
    % "Multi-objective optimization of the the dynamic aperture 
    % for the Swiss Light Source Upgrade"
    %
    nLines      = getoption(varargin,'nLines',15);   % number of lines to evaluate DA
    DA0         = getoption(varargin,'DA0',[]);         % linear DA (no non-linear elements)
    epsilon     = getoption(varargin,'epsilon',1e-6);
    b           = getoption(varargin,'b',0.25);
    nTurns      = getoption(varargin,'nTurns',128);
    dp          = getoption(varargin,'dp',0);
    makePlot    = getoption(varargin,'plot',0);
    verbose     = getoption(varargin,'verbose',1);
    parallel    = getoption(varargin,'parallel',1);
    thetas      = getoption(varargin,'thetas',0);
    getFuncs    = getoption(varargin,'getFuncs',0);
    R0          = getoption(varargin,'R0',false);    
    label       = getoption(varargin,'label','');
    weights     = [];

    if getFuncs
        out.findStartingRadius = @findStartingRadius;
        out.findIntersects = @findIntersects;
        out.isStable = @isStable;
        return
    end

    if parallel
        parforArg = Inf;
    else
        parforArg = 0;
    end

    % workflow
    % 1. determine search lines from angle and search limits from geometric acceptance
    % 2. along each search line follow algorithm from SLS-2 paper
    % 3. calculate area and output

    if ~isstruct(DA0)
        try
            DA0 = geometricAcceptance(ring,varargin{:});
        catch ME
            disp(ME);
            fprintf('Error while calculating geometric acceptance\n');
            return
        end
    end

    if ~thetas
        thetas = [0:nLines-1]*pi/(nLines-1);
    end

    R = zeros(1,length(thetas));
    parfor (i = 1:length(thetas), parforArg)
        % track 1 turn to initiate ringpass
        %
        p0 = [3e-5;0;3e-5;0;0;0];
        p1 = fatpass(ring,p0,1,1);

        theta = thetas(i);
        if ~R0
            L0 = findStartingRadius(theta,DA0);
        else
            L0 = R0(i);
        end
        linearDA(i) = L0;
        r_pos = L0;
        r_neg = 0;

        while (r_pos - r_neg) >= epsilon
            r = b*r_neg + (1-b)*r_pos;
            if isStable(theta,r,ring,nTurns,dp)
                r_neg = r;
            else
                r_pos = r;
            end
        end

        R(i) = (r_neg + r_pos)/2;
        if verbose
            fprintf('theta: %.1f deg, rho: %.1f mm\n',rad2deg(theta),R(i)*1e3);
        end
    end

    if makePlot
        figure(1311)
        [x,y] = pol2cart(thetas,R);
        plot(x*1e3,y*1e3,'LineWidth',2,'Marker','x','DisplayName',label);
        xlabel('x [mm]'); ylabel('y [mm]')
        grid on
    end

	dthetas = diff(thetas);
	r0 = R(1:(end-1));
	r1 = R(2:end);
	DA = sum(sin(dthetas) .* r0 .* r1 / 2.);

    out.area = DA;
    out.thetas = thetas;
    out.RMAXs = R';
    out.linearDA = linearDA';

end

function out = findStartingRadius(theta,DA0)
    % function to return the intersection point from DA search ray defined by theta
    % and the linear acceptance given as a set of x,y bounding coordinates

    x = DA0.x;
    y = DA0.y;

    if theta <= pi/2
        pos = find(x >= 0);
    else
        pos = find(x < 0);
    end

    x = x(pos);
    y = y(pos);
    A = [x' y'];

    % create line segment that's at least longer than max(abs(x))
    % 
    r = 0:max(abs(x))*0.1:max(abs(x))*1.1;
    t = theta*ones(1,length(r));
    [xx,yy] = pol2cart(t,r);
    B = [xx' yy'];

    try
        [xi,yi] = findIntersects(A,B);
        assert(length(xi) == length(yi) && length(xi) == 1, 'Wrong number of intersects found.');
        [~,out] = cart2pol(xi,yi);
    catch ME
        fprintf('find intersect failed for theta %.3f\n',theta);
        fprintf('Reverting to 0.1 [m] as starting radius\n');
        out = 0.1;
    end

end

function [xi,yi] = findIntersects(A,B)
    % find the cartesion coordinate of intersection points of two line segements
    % defined by matrix of points A = [x1;y1]; B = [x2;y2];

    xv = linspace(min(B(:,1)), max(B(:,1)));
    Ai = interp1(A(:,1), A(:,2), xv,'spline');
    Bi = interp1(B(:,1), B(:,2), xv);

    idx = find(diff(sign(Ai - Bi)));
    for i = 1:numel(idx)
        idxreg = max(1, idx(i)-2) : min(numel(xv), idx(i)+2);
        xi(i) = interp1(Ai(idxreg) - Bi(idxreg), xv(idxreg), 0);
        yi(i) = interp1(xv(idxreg), Bi(idxreg), xi(i));
    end
    [xi,i,~] = unique(xi);
    yi = yi(i);
end

function out = isStable(varargin)
    % tracking particle through ring and return binary depending on stability definition
    %
    theta = varargin{1};
    r = varargin{2};
    ring = varargin{3};
    nTurns = varargin{4};
    dp = varargin{5};

    minOffset = 3e-5;

    [x,y] = pol2cart(theta,r);
    if abs(x) < minOffset
        x = sign(x)*minOffset;
    end

    if abs(y) < minOffset
        y = sign(y)*minOffset;
    end

    p0 = [x;0;y;0;dp;0];

    % p1 = ringpass(ring,p0,nTurns,'KeepLattice');
    p1 = atpass(ring,p0,0,nTurns);
    p1 = p1(:,end); % get final turn particle coordinate

    % <<<<<<<<<<<<<<< check if need to define physical apertures and check if the particle
    % ever moves beyond the apertures
    %
    if any(isnan(p1))
        out = 0;
    else
        out = 1;
    end
end

function r = fatpass(varargin)
	% fatpass - fake atpass to circumvent some madlab parfor shizzle. I don't care anymore.
	r = atpass(varargin{:});
end

