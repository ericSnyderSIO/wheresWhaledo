% Set up data and whatnot

encounterStart = datenum([18 06 11 10 30 00]);
encounterEnd = datenum([18 06 11 12 30 00]);

% % load in DOA data:
% DOA{1} = load('D:\SOCAL_E_63\tracking\forRyan\SOCAL_E_63_EE\SOCAL_E_63_EE_C4_disk07_TDOA_DOA.mat');
% DOA{2} = load('D:\SOCAL_E_63\tracking\forRyan\SOCAL_E_63_EW\SOCAL_E_63_EW_C4_disk07_TDOA_DOA.mat');
% 
% for idoa = 1:numel(DOA)
%     Idoa{idoa} = find(DOA{idoa}.TDet>=encounterStart & DOA{idoa}.TDet<=encounterEnd);
% end

c = 1488.4;

%% Detect clicks on all instruments

% load in xwav tables
XH{1} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EE_C4_xwavLookupTable.mat');
XH{2} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EW_C4_xwavLookupTable.mat');
XH{3} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EN_xwavLookupTable.mat');
XH{4} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_ES_xwavLookupTable.mat');

%% load in H matrices

hyd1 = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EE_Hmatrix_fromHydLocInversion_210702.mat');
hyd2 = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EW_Hmatrix_fromHydLocInversion_210702.mat');

% HEW = H;

% Reorder hydrophones to fit new TDOA order
H{1} = [hyd1.hydPos(2,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(3,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(3,:)-hyd1.hydPos(2,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(2,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(3,:)];
     
H{2} = [hyd2.hydPos(2,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(3,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(3,:)-hyd2.hydPos(2,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(2,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(3,:)];


%% Run Detector
% load('detections.mat')

% [DET{1}] = detectClicks_4ch_SOCAL_E_63(encounterStart, encounterEnd, XH{1}.xwavTable, H{1}, c, 'detClicks_4ch.params');
% [DET{2}] = detectClicks_4ch_SOCAL_E_63(encounterStart, encounterEnd, XH{2}.xwavTable, H{2}, c, 'detClicks_4ch.params');
 
% [DET{3}] = detectClicks_1ch(encounterStart, encounterEnd, XH{3}.xwavTable, 'detClicks_1ch.params');
% [DET{4}] = detectClicks_1ch(encounterStart, encounterEnd, XH{4}.xwavTable, 'detClicks_1ch.params');
 
% save('detections', 'DET')


%% run brushDOA
% load('detections.mat')
% load('detections_brushDOA.mat')
% [det1, det2] = brushDOA(DET{1}, DET{2});
% 
% DET{1} = det1;
% DET{2} = det2;
 
% save('detections_brushDOA.mat', 'DET')

%% step through 30 seconds at a time, brushdet, validate detections, click train corr, coarse localization
load('detections_brushDOA.mat')
M200 = load('D:\SOCAL_E_63\tracking\experiments\largeApertureTDOA\TDOAmodel_200m.mat'); 

arrno = 2; % array with best detections
wn = 1;
spd = 24*60*60;
fsct = 100e3;
N = 1024;
Wk = hann(N/4);


I = find(DET{arrno}.color==wn+2);
tstart = DET{arrno}.TDet(I(1));
tend = DET{arrno}.TDet(I(end));

t1 = tstart;
t2 = tstart + 30/spd;
while t2<tend
    
    tct = t1:1/(spd*fsct):t2;

    XCT = zeros([length(tct), 4]);

    for ih = 1:4 % iterate through each instrument
        if ih == arrno % if 
            I30 = find(DET{ih}.TDet>=t1 & DET{ih}.TDet<=t2 & DET{ih}.color==wn+2);
            tdet = DET{ih}.TDet(I30);
        else
            
            I30 = find(DET{ih}.TDet>=t1 & DET{ih}.TDet<=t2);
            tdet = DET{ih}.TDet(I30);
        end
        
        IndDet = round((tdet - tct(1))*spd*fsct) + 1;
        
%         XCT(:,ih) = zeros(size(tct));
        XCT(IndDet,ih) = 1;

        XCT2(:, ih) = conv(XCT(:,ih), Wk);
    end

    ok = 1;
end



%% Step through each detection, validate, select good arrays, coarse and fine localization

%% Localize each detection

%% brush localizations