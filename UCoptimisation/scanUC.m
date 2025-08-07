function res = scanUC(uc,varargin)
    % DESCRIPTION
    % -----------
    % function to scan the main dipole focusing strength of the unit cell 
    % for a MBA type lattice
    %
    % INPUTS
    % ------
    % 'uc'::  AT lattice for the unit cell
    %
    % OPTIONS
    % -------
    % 'plot'    (0)::
    %   plot the lattice functions for the whole scan range
    % 'dipoleFamName'   ('B1')::
    %   fam name of the main dipole to scan
    % 'bk1'   ([1:-0.05:-2.5])::
    %   scan range for the quadrupole k1 [m-1] component
    % 'targetTune'  ([0.4,0.1])::
    %   target betatron tune value to match
    % 'Q1'  ('CF')::
    %   fam name of the first quadrupole family used to match tune
    % 'Q2'  ('CD')::
    %   fam name of the second quadrupole family used to match tune
    % 'SF'  ('SF1')::
    %   fam name of the focusing sextupole in the unit cell
    % 'SD'  ('SD1')::
    %   fam name of the defocusing sextupole in the unit cell
    % 'k1Limit'  (10)::
    %   upper limit for allowed quadrupole gradient [m-1]
    %
    % workflow
    % 1. Scan main dipole focusing angle
    % 2. For each angle, use CF and CD to match cell tune to 0.4/0.1
    % 3. Get lattice parameters for each focusing angle
    %
    % Author: F. Zhang - ANSTO, 2025
    % fzhang@ansto.gov.au
    %

    warning('off','AT:InconsistentK');
    plotFlag = getoption(varargin,'plot',0);
    dipoleFamName = getoption(varargin,'dipoleFamName','B1');
    bk1 = getoption(varargin,'bk1',[1:-0.05:-2.5]);
    targetTune = getoption(varargin,'targetTune',[0.4,0.1]);
    Q1 = getoption(varargin,'Q1','CF');
    Q2 = getoption(varargin,'Q2','CD');
    SF = getoption(varargin,'SF','SF1');
    SD = getoption(varargin,'SD','SD1');
    k1Limit = getoption(varargin,'k1Limit',10);

    colors = {};
    colors.blue   = "#0072BD";
    colors.green  = "#77AC30";
    colors.orange = "#D95319";
    colors.yellow = "#EDB120";

    colors.lightBlue = [0.47, 0.89, 1.0];
    colors.darkBlue  = [0.0, 0.24, 0.60];
    colors.lightGreen = [0.76, 1.0, 0.46];
    colors.darkGreen  = [0.14, 0.53, 0.00];
    colors.lightOrange = [1.0, 0.73, 0.43];
    colors.darkOrange  = [0.74, 0.29, 0.0];

    % reduce nominal k1 value to zero slowly to ensure closed orbit is preserved
    %
    b1_ord = atgetcells(uc,'FamName',dipoleFamName);
    b1_ord = find(b1_ord);
    bk1_init = atgetfieldvalues(uc,b1_ord(1),'PolynomB',{2});
    for i = [bk1_init:sign(max(bk1)-bk1_init)*0.05:max(bk1)]
        uc = atsetfieldvalues(uc,b1_ord,'PolynomB',{2},i);
        uc = atfittune(uc,[0.4,0.1],Q1,Q2);
    end

    res = {};
    for i = bk1
        uc = atsetfieldvalues(uc,b1_ord,'PolynomB',{2},i);
        uc = atfittune(uc,targetTune,Q1,Q2);

        % uc = atfitchrom(uc,[0,0],SF,SD);
        % use custom function is faster, based on linear chromaticity response matrix
        %
        uc = fitChromaticity(uc,{SF,SD});

        res{end+1} = uc;
    end

    for i = 1:length(res)
        ords = find(atgetcells(uc,'FamName',SF));
        k2 = atgetfieldvalues(res{i},ords,'PolynomB',{3});
        m1ls(i) = k2(1);

        ords = find(atgetcells(uc,'FamName',SD));
        k2 = atgetfieldvalues(res{i},ords,'PolynomB',{3});
        m2ls(i) = k2(1);
    end

    % ANALYTICAL EXPRESSION FOR SEXTUPOLE STRENGTHS
    %
    % m1ls = (-4*pi./etax1).*((chromx.*betay2) - (chromy.*betax2))./((betax1.*betay2) - (betax2.*betay1));
    % m2ls = (-4*pi./etax2).*((chromx.*betay1) - (chromy.*betax1))./((betax1.*betay2) - (betax2.*betay1));

    for i = 1:length(res)
        rp = atsummary(res{i});
        emx(i) = rp.naturalEmittance;
        Jx(i) = rp.damping(1);
        alphac(i) = rp.compactionFactor;
        Eloss(i) = rp.radiation*1e6;    % convert from GeV to keV

        ord = atgetcells(res{i},'FamName',Q1);
        val = atgetfieldvalues(res{i},ord,'PolynomB',{2});
        CFk1(i) = val(1);

        ord = atgetcells(res{i},'FamName',Q2);
        val = atgetfieldvalues(res{i},ord,'PolynomB',{2});
        CDk1(i) = val(1);

        cell = res{i};
        ord = atgetcells(cell,'Class','Sextupole');
        cell = atsetfieldvalues(cell,ord,'PolynomB',{3},0);
        [~,chrom] = tunechrom(cell,'get_chrom');
        chromaticity(i,:) = chrom;
    end

    % create patch for forbidden region based on CF/CD quad strength
    %
    CFlimit = find(abs(CFk1)>k1Limit);

    diff = CFlimit(2:end) - CFlimit(1:end-1);
    diff = find(diff>1);
    if isempty(diff)
        CFlimit = CFlimit(end);
    else
        CFlimit = CFlimit(diff);
    end

    CDlimit = find(abs(CDk1)>k1Limit);

    diff = CDlimit(2:end) - CDlimit(1:end-1);
    diff = find(diff>1);
    if isempty(diff)
        CDlimit = CDlimit(end);
    else
        CDlimit = CDlimit(diff);
    end

    bk1Limit = max(CFlimit,CDlimit);
    disp(bk1Limit);
    bk1Limit = bk1(bk1Limit);
    fprintf('bk1 limit: %.2f 1/m\n',bk1Limit);
    patchX = [bk1Limit max(bk1) max(bk1) bk1Limit];

    % ----------
    %
    figure(8071)
    clf
    t = tiledlayout('vertical','TileSpacing','compact','Padding','compact');
    nexttile

    yyaxis left
    plot(bk1,CFk1,'LineWidth',2,'DisplayName','Q1 k_1'); hold on
    plot(bk1,CDk1,'LineWidth',2,'DisplayName','Q2 k_1'); hold on
    ylabel('k_1 [m^{-1}]')
    ylim_val = [min([CFk1 CDk1]) max([CFk1 CDk1])];
    ylim(ylim_val);
    patchY = [ylim_val(1) ylim_val(1) ylim_val(2) ylim_val(2)];
    patch(patchX,patchY,'red','EdgeColor','none','FaceAlpha',0.2,'DisplayName','k_1 limit');

    yyaxis right
    plot(bk1,1e12*emx,'LineWidth',2,'DisplayName','\epsilon_x'); hold on
    % ylim([80 120])
    ylabel('\epsilon_x [pm]')
    grid on
    legend('Orientation','horizontal');

    nexttile

    plot(bk1,abs(m1ls),'LineWidth',2,'DisplayName','SF','LineStyle','-'); hold on
    plot(bk1,abs(m2ls),'LineWidth',2,'DisplayName','SD'); hold on
    rmsk2 = sqrt(m1ls.^2+m2ls.^2);
    plot(bk1,rmsk2,'LineWidth',2,'DisplayName','Sum'); hold on
    ylabel('|k_2| [m^{-2}]');
    ylim_val = [min([abs(m1ls) abs(m2ls) rmsk2]) max([abs(m1ls) abs(m2ls) rmsk2])];
    ylim(ylim_val);
    patchY = [ylim_val(1) ylim_val(1) ylim_val(2) ylim_val(2)];
    patch(patchX,patchY,'red','EdgeColor','none','FaceAlpha',0.2,'DisplayName','k_1 limit');
    grid on
    xlabel('dipole k1 [m^{-1}]');
    legend('Orientation','horizontal');

    % ----------
    %
    figure(8072)
    t = tiledlayout('vertical','TileSpacing','compact','Padding','compact');
    nexttile

    yyaxis left
    plot(bk1,1e12*emx,'LineWidth',2,'DisplayName','\epsilon_x'); hold on
    ylabel('\epsilon_x [pm]')
    ylim_val = [1e12*min(emx) 1e12*max(emx)];
    ylim(ylim_val);
    patchY = [ylim_val(1) ylim_val(1) ylim_val(2) ylim_val(2)];
    patch(patchX,patchY,'red','EdgeColor','none','FaceAlpha',0.2,'DisplayName','k_1 limit');

    yyaxis right
    plot(bk1,Eloss,'LineWidth',2,'DisplayName','E_{loss}/turn'); hold on
    ylabel('E_{loss}/turn [keV]')
    grid on
    xlabel('dipole k1 [m^{-1}]');
    legend('Orientation','horizontal')

    nexttile
    
    yyaxis left
    plot(bk1,Jx,'LineWidth',2,'DisplayName','J_x'); hold on
    ylabel('J_x')
    ylim_val = [min(Jx) max(Jx)];
    ylim(ylim_val);
    patchY = [ylim_val(1) ylim_val(1) ylim_val(2) ylim_val(2)];
    patch(patchX,patchY,'red','EdgeColor','none','FaceAlpha',0.2,'DisplayName','k_1 limit');

    yyaxis right
    plot(bk1,alphac,'LineWidth',2,'DisplayName','\alpha_c'); hold on
    ylabel('\alpha_c')
    grid on
    xlabel('dipole k1 [m^{-1}]');
    legend('Orientation','horizontal')

    % ----------
    %
    if plotFlag
        nGrad = length(res);
        colors.gradBlue = [...
            linspace(colors.lightBlue(1),colors.darkBlue(1),nGrad)',...
            linspace(colors.lightBlue(2),colors.darkBlue(2),nGrad)',...
            linspace(colors.lightBlue(3),colors.darkBlue(3),nGrad)',...
            ];
        colors.gradGreen = [...
            linspace(colors.lightGreen(1),colors.darkGreen(1),nGrad)',...
            linspace(colors.lightGreen(2),colors.darkGreen(2),nGrad)',...
            linspace(colors.lightGreen(3),colors.darkGreen(3),nGrad)',...
            ];
        colors.gradOrange = [...
            linspace(colors.lightOrange(1),colors.darkOrange(1),nGrad)',...
            linspace(colors.lightOrange(2),colors.darkOrange(2),nGrad)',...
            linspace(colors.lightOrange(3),colors.darkOrange(3),nGrad)',...
            ];
        for n = 1:length(res)
            uc = res{n};

            figure(4579)
            hold off
            curve = atbaseplot(uc);

            figure(4580)
            hold on
            yyaxis left
            plot(curve.left(1).XData,curve.left(1).YData,'LineWidth',1,'Marker','None','LineStyle','-','color',colors.gradBlue(n,:)); hold on
            plot(curve.left(2).XData,curve.left(2).YData,'LineWidth',1,'Marker','None','LineStyle','-','color',colors.gradGreen(n,:)); hold on
            ylabel('\beta [m]')
            yyaxis right
            plot(curve.right(1).XData,curve.right(1).YData,'LineWidth',1,'Marker','None','LineStyle','-','color',colors.gradOrange(n,:)); hold on
            ylabel('\eta_x [m]')
        end
        grid on
        xlabel('s [m]'); 
        yyaxis left
        ax = gca;
        la = atplotsyn(ax,ucIn);
        legend({'','','\beta_x','\beta_y'})
    end

end

function chromRM = getChromRM(uc,SXfams)
    deltak2 = 0.1;
    chromRM = zeros(2,2);

    [rd,~] = atlinopt4(uc,'get_chrom');
    oldChrom = rd.chromaticity;

    for n = 1:length(SXfams)
        sx = atgetcells(uc,'FamName',SXfams{n});
        oldk2 = atgetfieldvalues(uc,sx,'PolynomB',{3});
        uc1 = atsetfieldvalues(uc,sx,'PolynomB',{3},oldk2(1)+deltak2);

        [rd,~] = atlinopt4(uc1,'get_chrom');
        newChrom = rd.chromaticity;
        dChrom = newChrom - oldChrom;
        chromRM(:,n) = dChrom/deltak2;
    end
end

function uc1 = fitChromaticity(uc,SXfams)
    targetX = [1,1];

    chromRM = getChromRM(uc,SXfams);
    chromRMinv = pinv(chromRM);

    [rd,~] = atlinopt4(uc,'get_chrom');
    oldChrom = rd.chromaticity;

    dX = targetX - oldChrom;
    dK = chromRMinv * dX';

    uc1 = uc;
    for n = 1:length(SXfams)
        sx = atgetcells(uc,'FamName',SXfams{n});
        oldk2 = atgetfieldvalues(uc,sx,'PolynomB',{3});
        uc1 = atsetfieldvalues(uc1,sx,'PolynomB',{3},oldk2(1)+dK(n));
    end
end
