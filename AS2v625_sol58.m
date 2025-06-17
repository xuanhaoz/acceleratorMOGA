function ring = AS2v625
    % created 2025/04/09
    % modified from AS2v625
    % changed SX and OC fields based on MOGA results for sol 58
    %
    disp('Loading AS2 version 625 lattice');

    bdsl = 0.024;
    bdst = 0;   % variation to dispersion suppressor bending angle
    rb = -0.04;     % reverse bending angle
    cdb = -0.04;
    bdstratio = 1;

    CF   = atrbend('CF', 'Length',0.12, 'BendingAngle',deg2rad(rb),...
                    'K',10.194625,  'NumIntSteps',100,'MaxOrder',1);
    CD   = atrbend('CD', 'Length',0.10, 'BendingAngle',deg2rad(cdb-rb),...
                    'K',-1.503839, 'NumIntSteps',100,'MaxOrder',1);
    B1   = atrbend('B1', 'Length',0.703, 'BendingAngle',deg2rad(3-2*cdb),...
                    'K',-1.9,      'NumIntSteps',100,'MaxOrder',1);
    CFDS = atrbend('CFDS', 'Length',0.12, 'BendingAngle',deg2rad(rb-bdst/bdstratio),...
                    'K',10.809547,  'NumIntSteps',100,'MaxOrder',1);
    CDDS = atrbend('CDDS', 'Length',0.12, 'BendingAngle',deg2rad(cdb-rb-bdst*(bdstratio-1)/bdstratio),...
                    'K',-3.748819, 'NumIntSteps',100,'MaxOrder',1);
    BDS  = atrbend('BDS','Length',0.663/2+bdsl,'BendingAngle',deg2rad(3/2-cdb+bdst),...
                    'K',-2.560427, 'NumIntSteps',100,'MaxOrder',1);

    % 'PassMethod','BndMPoleSymplectic4E2Pass',...
    CF.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    CD.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    B1.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    CFDS.('PassMethod')  = 'BndMPoleSymplectic4Pass';
    CDDS.('PassMethod')  = 'BndMPoleSymplectic4Pass';
    BDS.('PassMethod')   = 'BndMPoleSymplectic4Pass';

    QMS1   = atquadrupole('QMS1','Length',0.12, 'K', 11.000000,'NumIntSteps',10,'MaxOrder',1);
    QMS2   = atquadrupole('QMS2','Length',0.12, 'K', -2.825421,'NumIntSteps',10,'MaxOrder',1);
    QMS3   = atquadrupole('QMS3','Length',0.12, 'K', -9.689775,'NumIntSteps',10,'MaxOrder',1);
    QMS4   = atquadrupole('QMS4','Length',0.12, 'K',  9.391910,'NumIntSteps',10,'MaxOrder',1);

    % 'PassMethod','StrMPoleSymplectic4RadPass',...
    QMS1.('PassMethod') = 'StrMPoleSymplectic4Pass';
    QMS2.('PassMethod') = 'StrMPoleSymplectic4Pass';
    QMS3.('PassMethod') = 'StrMPoleSymplectic4Pass';
    QMS4.('PassMethod') = 'StrMPoleSymplectic4Pass';

    SF0 = atsextupole('SF1','Length',0.065/2,'S',  864.556548,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF1 = atsextupole('SF1','Length',0.065,'S', 1.290273315171157e+03,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF2 = atsextupole('SF2','Length',0.065,'S', 7.926408497275762e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF3 = atsextupole('SF3','Length',0.065,'S', 9.552922156909447e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF4 = atsextupole('SF4','Length',0.065,'S', 7.347556021189174e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF5 = atsextupole('SF5','Length',0.065,'S', 1.079179012976061e+03,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD1 = atsextupole('SD1','Length',0.065,'S',-7.294464109757021e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD2 = atsextupole('SD2','Length',0.065,'S',-8.654593820710700e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD3 = atsextupole('SD3','Length',0.065,'S',-7.802059753795896e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD4 = atsextupole('SD4','Length',0.065,'S',-9.404580882851570e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD5 = atsextupole('SD5','Length',0.065,'S',-7.621566621766000e+02,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');

    SF1.('PolynomB')(3)  = SF1.('S');
    SF2.('PolynomB')(3)  = SF2.('S');
    SF3.('PolynomB')(3)  = SF3.('S');
    SF4.('PolynomB')(3)  = SF4.('S');
    SF5.('PolynomB')(3)  = SF5.('S');
    SD1.('PolynomB')(3)  = SD1.('S');
    SD2.('PolynomB')(3)  = SD2.('S');
    SD3.('PolynomB')(3)  = SD3.('S');
    SD4.('PolynomB')(3)  = SD4.('S');
    SD5.('PolynomB')(3)  = SD5.('S');

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % octupole strength conversion from OPA need to change from integrated strength to normalised
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OC1 = atmultipole('OC1','Length',0.055,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0  1.165246318823110e+04],'PassMethod','StrMPoleSymplectic4Pass');
    OC2 = atmultipole('OC2','Length',0.055,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0  2.047949528948905e+04],'PassMethod','StrMPoleSymplectic4Pass');
    OC3 = atmultipole('OC3','Length',0.055,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0 -3.120386650214563e+04],'PassMethod','StrMPoleSymplectic4Pass');
    OC4 = atmultipole('OC4','Length',0.055,'PolynomA',[0 0 0 0],'PolynomB',[0 0 0  6.209939804813499e+03],'PassMethod','StrMPoleSymplectic4Pass');
     
    D1     = atdrift('Drift',0.119,'PassMethod','DriftPass');
    D1h    = atdrift('Drift',0.119/2,'PassMethod','DriftPass');
    D2     = atdrift('Drift',0.119,'PassMethod','DriftPass');
    D2h    = atdrift('Drift',0.119/2,'PassMethod','DriftPass');
    D3     = atdrift('Drift',0.119,'PassMethod','DriftPass');
    D3h    = atdrift('Drift',0.119/2,'PassMethod','DriftPass');
    D4     = atdrift('Drift',0.050,'PassMethod','DriftPass');
    D4h    = atdrift('Drift',0.050/2,'PassMethod','DriftPass');

    D25    = atdrift('Drift',0.0250,'PassMethod','DriftPass');
    D37    = atdrift('Drift',0.0375,'PassMethod','DriftPass');
    D75    = atdrift('Drift',0.0750,'PassMethod','DriftPass');
    D55    = atdrift('Drift',0.0550,'PassMethod','DriftPass');
    D50    = atdrift('Drift',0.0500,'PassMethod','DriftPass');
    D57    = atdrift('Drift',0.0570,'PassMethod','DriftPass');
    D022   = atdrift('Drift',0.055/2,'PassMethod','DriftPass');
    D150   = atdrift('Drift',0.1500,'PassMethod','DriftPass');
    D100   = atdrift('Drift',0.1000,'PassMethod','DriftPass');
    D130   = atdrift('Drift',0.1300,'PassMethod','DriftPass');
    D250   = atdrift('Drift',0.2500,'PassMethod','DriftPass');
    DDS3   = atdrift('Drift',0.25967+0.00833-bdsl,'PassMethod','DriftPass');
    D2500  = atdrift('Drift',2.5000,'PassMethod','DriftPass');
    DMS    = atdrift('Drift',0.119,'PassMethod','DriftPass');
    DMSh   = atdrift('Drift',0.119/2,'PassMethod','DriftPass');

    BPM = @(name) atmonitor(name, 'IdentityPass');
    MARK = @(name) atmarker(name, 'IdentityPass');

    COR = @(name) atcorrector(name);

    RFC  = atrfcavity('RFCav',...
        'Energy',3E9);
        % 'Voltage',2e6,...
        % 'HarmNumber',758,...
        % 'Frequency',500e6,...
    	% 'PassMethod','RFCavityPass');

    AP  =   ataperture('AP',[-10 10 -10 10]*1e-3,'AperturePass');


    % Unit Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    UC1 = {D1 SF1 D1 CF D2 CD D3h BPM('BPM5') D3h SD2 D4 B1 D4 SD2 D3h BPM('BPM6')  D3h CD D2h COR('COR4') D2h CF};
    UC2 = {D1 SF2 D1 CF D2 CD D3h MARK('GirderEnd') MARK('GirderStart') BPM('BPM7')  D3h SD3 D4 B1 D4 SD3 D3h BPM('BPM8')  D3h CD D2h COR('COR5') D2h CF};
    UC3 = {D1 SF3 D1 CF D2 CD D3h BPM('BPM9')  D3h SD4 D4 B1 D4 MARK('GirderEnd') MARK('GirderStart') SD4 D3h BPM('BPM10') D3h CD D2h COR('COR6') D2h CF};
    UC4 = {D1 SF4 D1 CF D2 CD D3h BPM('BPM11') D3h SD5 D4 B1 D4 SD5 D3h BPM('BPM12') D3h CD D2h COR('COR7') D2h CF};
    UC5 = {D1 SF5 D1};

    arc = [UC1 UC2 UC3 UC4 UC5];

    % Dispersion Suppressor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    dsPartA = {CFDS D2h COR('COR3') D2h CDDS D3h BPM('BPM4') D3h SD1 MARK('GirderStart') D4 MARK('GirderEnd') DDS3 BPM('BPM3') BDS};
    dsPartB = {CFDS D2 CDDS D3h MARK('GirderEnd') MARK('GirderStart') D3h SD1 D4 DDS3 BPM('BPM13') BDS};

    dsA = [flip(dsPartA)];
    dsB = [dsPartB];

    % Matching Straight %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    msPartA = {DMS QMS1 DMS OC1 DMS QMS2 DMS OC2 DMSh COR('COR2') BPM('BPM2')  DMSh QMS3 DMS OC3 DMS QMS4 DMS OC4 COR('COR1') BPM('BPM1')  MARK('GirderStart') D2500};
    msPartB = {DMS QMS1 DMS OC1 DMS QMS2 DMS OC2 DMSh COR('COR8') BPM('BPM14') DMSh QMS3 DMS OC3 DMS QMS4 DMS OC4 COR('COR9') BPM('BPM15') MARK('GirderEnd')   D2500};

    msA = [flip(msPartA)];
    msB = [msPartB];

    msPartAnoOct = {DMS QMS1 DMS D55 DMS QMS2 DMS D55 DMSh COR('COR2') BPM('BPM2')  DMSh QMS3 DMS D55 DMS QMS4 DMS D55 COR('COR1') BPM('BPM1')  MARK('GirderStart') D2500};
    msPartBnoOct = {DMS QMS1 DMS D55 DMS QMS2 DMS D55 DMSh COR('COR8') BPM('BPM14') DMSh QMS3 DMS D55 DMS QMS4 DMS D55 COR('COR9') BPM('BPM15') MARK('GirderEnd')   D2500};

    msAnoOct = [flip(msPartAnoOct)];
    msBnoOct = [msPartBnoOct];

    % Ring %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %

    sector = [msA dsA arc dsB msB];
    ring = [msAnoOct dsA arc dsB msB repmat(sector,1,22) msA dsA arc dsB msBnoOct];

    RP = atringparam('AS2',3.0e9);
    ring = [{RP}; {AP}; {RFC}; ring'];
    % --------------------------------------
    % additional cavity at midway 29/04/2024
    %
    % ring = [{RP}; {RFC}; repmat(sector,12,1); {RFC}; repmat(sector,12,1)];
    % --------------------------------------
    ring = atsetcavity(ring,3e6,1,784);
    % ring = [{RP}; {RFC}; ds0'; msB' ];

    % ring = atsetcavity(ring,'frequency','nominal',...
    %         'harmnumber',758,'voltage',2e6);

    % ring = atsetringproperties(ring);
    % ring = atsetrfcavity(ring,2e6,1,758,0);

end