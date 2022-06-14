% step through 30 seconds at a time, brushdet, validate detections, click train corr, coarse localization
load('detections_brushDOA180611_1030')

arrno = 2; % array with best detections
wn = 1;
spd = 24*60*60;
fsct = 100e3;
maxLag = round(fsct*2000/1500);
twin = 15; % window length for click train
N = 1024;
Wk = hann(N/4);
xcovInd = [2,3,4,7,8,12];

xct = zeros(twin*fsct, 4);

encStart = min([DET{1}.TDet(1), DET{2}.TDet(1)]);
encEnd = min([DET{1}.TDet(end), DET{2}.TDet(end)]);


TDet = [];
TDOA = [];
% tstart = (encStart+encEnd)/2;
tstart = encStart;
tend = tstart + twin/spd;
ndet = 0;
while tend <= encEnd

    uniqueLabels = unique([DET{1}.color; DET{2}.color]);

    for ih = 1:4 % iterate through each instrument
        tct = tstart:1/(spd*fsct):tend;
        I{ih} = find(DET{ih}.TDet>=tstart & DET{ih}.TDet<=tend);
        tdet = DET{ih}.TDet(I{ih});


        for idet = 1:length(tdet)
            [~, ind] = min(abs(tct-tdet(idet)));
            xct(ind, ih) = 1;
        end
        xctHann(:, ih) = conv(Wk, xct(:, ih));

    end
    if 0 % plot or not
        figure(1)
        for ih = 1:4
            subplot(4,1,ih)
            plot(xctHann(:, ih))
        end
    end

    [xctCorr, lags] = xcorr(xctHann, maxLag);

    tdoa = -99.*ones(length(uniqueLabels), length(xcovInd)); % initialize with -99 so I can easily remove faulty TDOAs
    for np = 1:length(xcovInd)
        numWhales = length(uniqueLabels);

        [pks, locs] = findpeaks(xctCorr(:, xcovInd(np)), 'SortStr','descend', 'NPeaks', numWhales);

        if max(pks)>150
            tdoa(:, np) = lags(locs)./fsct;
        end

        if 0 % plot or not
            figure(2)
            subplot(6,1,np)
            plot(xctCorr(:, xcovInd(np)))
        end
    end
    
    tdet = tstart.*ones(length(uniqueLabels), 1);
    TDet = [TDet, tdet.'];
    TDOA = [TDOA; tdoa];

    tstart = tend;
    tend = tstart + twin/spd;
    
end

%%
fig = figure(3);
for np = 1:6
    subplot(6, 1, np)
    Iplt = find(TDOA(:, np)~=-99);
    plot(TDet(Iplt), TDOA(Iplt, np), '.')
    title(['Pair ', num2str(np), ', window length = ', num2str(twin), ' s'])
    ylabel('TDOA')
    datetick
end


saveas(fig, ['windowLength_', num2str(twin)], 'fig')