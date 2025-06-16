function out = geometricAcceptance(ring,varargin)
    % calculate the linear geometric acceptance of the lattice
    % referenced in "insideOPA"
    %
    s = getoption(varargin,'s',0);
    ax = getoption(varargin,'ax',10e-3);
    ay = getoption(varargin,'ay',10e-3);
    % k = getoption(varargin,'k',0:0.1:1);      % h/v coupling ratio
    k = getoption(varargin,'k',linspace(0,1,10));      % h/v coupling ratio
    dp = getoption(varargin,'dp',0);
    makeplot = getoption(varargin,'plot',0);      % h/v coupling ratio
    verbose = getoption(varargin,'v',0);

    [~,ed] = atlinopt4(ring,1:length(ring)+1,'dp',dp);
    beta = cat(1,ed.beta);
    betax = beta(:,1);
    betay = beta(:,2);
    x0 = [ed.ClosedOrbit];
    x0 = x0(1,:)';
    y0 = 0;     % assume flat lattice

    n = @(k) (ay^2 * betax * (1-k) + ax^2 * betay * k);
    p = @(k) ((ay^2 * abs(x0) * sqrt(1-k) .* sqrt(betax))./n(k));
    q = @(k) ((ay^2 * (x0.^2 - ax^2))./n(k));

    sqrtA = @(k) (-p(k) + sqrt(p(k).^2 - q(k)));

    A = sqrtA(k); 
    [A,Amin] = min(A);
    A = A.^2;

    % get optics at position s to evaluate the geometric acceptance
    %
    spos = findspos(ring,1:length(ring)+1);
    [~,index] = min(abs(s - spos));
    [~,ed] = atlinopt4(ring,index,'dp',dp);
    betax_s = ed.beta(1);
    betay_s = ed.beta(2);
    x0_s = ed.ClosedOrbit(1);

    xA = [x0_s - sqrt((1-k).*A*betax_s); x0_s + sqrt((1-k).*A*betax_s)];
    yA = sqrt(k.*A*betay_s);

    if verbose
        fprintf('max x acceptance: %.2e, %.2e mm\n',min(xA(1))*1e3,max(xA(2))*1e3);
        fprintf('max y acceptance: %.2e mm\n',max(yA)*1e3);
    end

    out.xA = xA;
    out.yA = yA;

    x = [xA(1,:) flip(xA(2,:))];
    y = [yA flip(yA)];

    [x,ix,~] = unique(x);
    y = y(ix);

    out.x = x;
    out.y = y;

    if makeplot
        figure(1311)
        hold on
        ax_handle = gca;
        lines = get(ax_handle,'Children');
        tags = [];
        if length(lines) > 0
            tags = [lines.Tag];
        end

        if ~strcmp(tags,'linearDA')
            plot(out.x*1e3, out.y*1e3,'LineWidth',2,'DisplayName','Geometric acceptance','Tag','linearDA');
            grid on
            xlabel('x [mm]'); ylabel('y [mm]');
            legend();
        end
    end
end