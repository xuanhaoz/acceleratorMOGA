function ring = UC_AS2v625
    % created 2025/04/03
    % modified from UC_AS2v625
    % extended main dipole closer to SD, gap size is now 50 mm
    %
    disp('Loading AS2 version 625 lattice');

    bdsl = 0.11;
    bdst = 0;   % variation to dispersion suppressor bending angle
    rb = -0.08;     % reverse bending angle
    cdb = -0.08;
    bdstratio = 1;

    CF   = atrbend('CF', 'Length',0.12, 'BendingAngle',deg2rad(rb),...
                    'K',9.74920,  'NumIntSteps',100,'MaxOrder',1);
    CD   = atrbend('CD', 'Length',0.10, 'BendingAngle',deg2rad(cdb-rb),...
                    'K',-0.50340, 'NumIntSteps',100,'MaxOrder',1);
    B1   = atrbend('B1', 'Length',0.703, 'BendingAngle',deg2rad(3-2*cdb),...
                    'K',-1.970943,      'NumIntSteps',100,'MaxOrder',1);

    % 'PassMethod','BndMPoleSymplectic4E2Pass',...
    CF.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    CD.('PassMethod')    = 'BndMPoleSymplectic4Pass';
    B1.('PassMethod')    = 'BndMPoleSymplectic4Pass';

    SF1h= atsextupole('SF1','Length',0.065/2,'S',  548.371265,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD1h= atsextupole('SD1','Length',0.065/2,'S', -588.555764,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF1 = atsextupole('SF1','Length',0.065,'S',  548.371265,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF2 = atsextupole('SF2','Length',0.065,'S',  801.265453,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SF3 = atsextupole('SF3','Length',0.065,'S', 1107.546315,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD1 = atsextupole('SD1','Length',0.065,'S', -588.555764,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD2 = atsextupole('SD2','Length',0.065,'S', -532.404607,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD3 = atsextupole('SD3','Length',0.065,'S',-1007.911263,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD4 = atsextupole('SD4','Length',0.065,'S',-1007.900266,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');
    SD5 = atsextupole('SD5','Length',0.065,'S', -875.698257,'NumIntSteps',10,'MaxOrder',2,'PassMethod','StrMPoleSymplectic4Pass');

    SF1.('PolynomB')(3)  = SF1.('S');
    SF2.('PolynomB')(3)  = SF2.('S');
    SF3.('PolynomB')(3)  = SF3.('S');
    SD1.('PolynomB')(3)  = SD1.('S');
    SD2.('PolynomB')(3)  = SD2.('S');
    SD3.('PolynomB')(3)  = SD3.('S');
    SD4.('PolynomB')(3)  = SD4.('S');
    SD5.('PolynomB')(3)  = SD5.('S');

    D1     = atdrift('Drift',0.119,'PassMethod','DriftPass');
    D1h    = atdrift('Drift',0.119/2,'PassMethod','DriftPass');
    D2     = atdrift('Drift',0.119,'PassMethod','DriftPass');
    D2h    = atdrift('Drift',0.119/2,'PassMethod','DriftPass');
    D3     = atdrift('Drift',0.119,'PassMethod','DriftPass');
    D3h    = atdrift('Drift',0.119/2,'PassMethod','DriftPass');
    D4     = atdrift('Drift',0.05,'PassMethod','DriftPass');
    D4h    = atdrift('Drift',0.05/2,'PassMethod','DriftPass');

    D25    = atdrift('Drift',0.0250,'PassMethod','DriftPass');
    D37    = atdrift('Drift',0.0375,'PassMethod','DriftPass');
    D75    = atdrift('Drift',0.0750,'PassMethod','DriftPass');
    D50    = atdrift('Drift',0.0500,'PassMethod','DriftPass');
    D57    = atdrift('Drift',0.0570,'PassMethod','DriftPass');
    D022   = atdrift('Drift',0.055/2,'PassMethod','DriftPass');
    D150   = atdrift('Drift',0.1500,'PassMethod','DriftPass');
    D100   = atdrift('Drift',0.1000,'PassMethod','DriftPass');
    D130   = atdrift('Drift',0.1300,'PassMethod','DriftPass');
    D114   = atdrift('Drift',0.1140,'PassMethod','DriftPass');
    D250   = atdrift('Drift',0.2500,'PassMethod','DriftPass');
    DDS3   = atdrift('Drift',0.25967-bdsl,'PassMethod','DriftPass');
    D2500  = atdrift('Drift',2.5000,'PassMethod','DriftPass');

    BPM = @(name) atmonitor(name, 'IdentityPass');
    MARK = @(name) atmarker(name, 'IdentityPass');

    COR = @(name) atcorrector(name);

    RFC  = atrfcavity('RFCav',...
        'Energy',3E9);
        % 'Voltage',2e6,...
        % 'HarmNumber',758,...
        % 'Frequency',500e6,...
    	% 'PassMethod','RFCavityPass');

    M1  = atmarker('M1','IdentityPass');
    M2  = atmarker('M2','IdentityPass');
    M3  = atmarker('M2','IdentityPass');
    M4  = atmarker('M2','IdentityPass');
    MSD  = atmarker('MSD','IdentityPass');


    RP = atringparam('AS2',3.0e9);
    % ring = [{RP}; {RFC}; repmat(sector,22,1)];
    % --------------------------------------
    % additional cavity at midway 29/04/2024
    %
    % ring = [{RP}; {RFC}; repmat(sector,12,1); {RFC}; repmat(sector,12,1)];
    % --------------------------------------
    % ring = atsetcavity(ring,3e6,1,750);

    UC0 = {M1 SF1h D1 CF D2 CD D3h BPM('BPM5') D3h SD1h M2 SD1h D4 B1 D4 SD1h M3 SD1h D3h BPM('BPM6')  D3h CD D2h COR('COR4') D2h CF D1 SF1h M4};
    ring = [{RP}; UC0'];

    % ring = atsetcavity(ring,'frequency','nominal',...
    %         'harmnumber',758,'voltage',2e6);

    % ring = atsetringproperties(ring);
    % ring = atsetrfcavity(ring,2e6,1,758,0);

end