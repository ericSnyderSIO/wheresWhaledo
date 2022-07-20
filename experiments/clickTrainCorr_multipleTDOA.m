
trackName = 'track43_180327_084016'
% trackName = '180611_1030';

arrno = 2; % which array is primary array

saveFileName = [trackName, '_CTC_', 'Array', num2str(arrno)]

encounterStart = min(DET{arrno}.TDet);
encounterEnd = max(DET{arrno}.TDet);

tdir = dir(['*det*', trackName, '*.mat']);
load(tdir.name)

% encounterStart = min([min(DET{1}.TDet), min(DET{2}.TDet)]);
% encounterEnd = max([max(DET{1}.TDet), max(DET{2}.TDet)]);

%% Determine if detector needs to be run for single channels
if numel(DET)==2
    % load in xwav tables
    XH{1} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EE_C4_xwavLookupTable.mat');
    XH{2} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EW_C4_xwavLookupTable.mat');
    XH{3} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EN_xwavLookupTable.mat');
    XH{4} = load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_ES_xwavLookupTable.mat');

    DET{1} = outDet1;
    DET{2} = outDet2;


    % SOCAL_E_63_EN
    [DET{3}] = detectClicks_1ch(encounterStart, encounterEnd, XH{3}.xwavTable, 'detClicks_1ch.params');

    % SOCAL_E_63_ES
    [DET{4}] = detectClicks_1ch(encounterStart, encounterEnd, XH{4}.xwavTable, 'detClicks_1ch.params');

    DET = fixAngle(DET);

    save(tdir.name, 'DET')

end

%% Set parameters for CTC

global brushing
loadParams('brushing.params')
spd = 24*60*60;
fsct = 10e3;
c = 1488.4;
maxLag = round(fsct*(2000/1500 + .4));
twin = 30; % window length for click train
whaleSpeed = 3; % m/s, overestimate of whale speed
Nhann = ceil((twin*whaleSpeed/c)*fsct/2);
if mod(Nhann, 2)==1 % if Nhann is odd, add one to make it even
    Nhann = Nhann + 1;
end
Wk = hann(Nhann);
maxNumTDOA = 4; % maximum possible TDOAs per pair per detection

xct = zeros(twin*fsct, 4);

encStart = min([DET{1}.TDet(1), DET{2}.TDet(1)]);
encEnd = min([DET{1}.TDet(end), DET{2}.TDet(end)]);


% assign indices of other arrays
if arrno==1
    otherArrays = [2,3,4];
    xcovInd = [2,3,4];
    hpair{1} = '1-2';
    hpair{2} = '1-3';
    hpair{3} = '1-4';
elseif arrno==2
    otherArrays = [1,3,4];
    xcovInd = [2,7,8];
    hpair{1} = '1-2';
    hpair{2} = '2-3';
    hpair{3} = '2-4';
end
%% remove detections that are too close together
mindt = Nhann/fsct;
for ih = 1:numel(DET)
    Irem = find(diff(DET{ih}.TDet)<mindt/spd);
    DET{ih}(Irem, :) = [];
end

%%


for wn = unique(DET{arrno}.color).'
    Ilab = find(DET{arrno}.color==wn); % detections labeled wn

    whale{wn}.TDOA = -99.*ones(length(Ilab), length(xcovInd), maxNumTDOA);
    whale{wn}.XCTpk = whale{wn}.TDOA;
    whale{wn}.SNR = whale{wn}.TDOA;
    whale{wn}.label = ones(length(Ilab), 1);

    for ndet = 1:length(Ilab) % iterate through each detection
        %         tic
        tdet = DET{arrno}.TDet(Ilab(ndet));

        tstart = tdet - twin/(spd*2);
        tend = tdet + twin/(spd*2);

        tct = tstart:1/(spd*fsct):tend;
        xct = zeros(length(tct), 4);

        % indices of detections within window & labeled as wn
        Iwn = find(DET{arrno}.TDet>=tstart & DET{arrno}.TDet<=tend & DET{arrno}.color==wn);

        for i = 1:length(Iwn)
            [~, ind] = min(abs(tct-DET{arrno}.TDet(Iwn(i))));
            xct(ind, arrno) = 1;
        end

        xctHann(:, arrno) = conv(Wk, xct(:, arrno));

        for ia = 1:3
            Iwn = find(DET{otherArrays(ia)}.TDet>=tstart & DET{otherArrays(ia)}.TDet<=tend);
            for i = 1:length(Iwn)

                [~, ind] = min(abs(tct-DET{otherArrays(ia)}.TDet(Iwn(i))));

                xct(ind, otherArrays(ia)) = 1;

            end

            xctHann(:, otherArrays(ia)) = conv(Wk, xct(:, otherArrays(ia)));
        end
        %         if 0 % plot or not
        %             figure(1)
        %             for ih = 1:4
        %                 subplot(4,1,ih)
        %                 plot(xctHann(:, ih))
        %             end
        %         end


        [tdoa, xcpk, snr] = calcTDOA_CTC(xctHann, maxLag, xcovInd, maxNumTDOA, Nhann, fsct);

        %         detCount = detCount + 1;
        whale{wn}.TDet(ndet) = tdet;
        whale{wn}.TDOA(ndet, :, :) = tdoa;
        whale{wn}.XCTpk(ndet, :, :) = xcpk;
        whale{wn}.SNR(ndet, :, :) = snr;
        whale{wn}.label(ndet) = wn;

        %         elapsedTime(ndet) = toc;
    end

end

%%

fig = figure(3);
for np = 1:3
    subplot(3, 1, np)
    for wn = 2:numel(whale)

        Iplt = find(whale{wn}.TDOA(:, np, 1)~=-99);
        if ~isempty(Iplt)

            for nt = 1:maxNumTDOA
                Iplt = find(whale{wn}.TDOA(:, np, nt)~=-99);
                scatter(whale{wn}.TDet(Iplt), whale{wn}.TDOA(Iplt, np, nt), 80.*whale{wn}.XCTpk(Iplt, np, nt)./max(max(max(whale{wn}.XCTpk))), brushing.params.colorMat(wn, :).*ones(length(Iplt), 3), 'filled')
                hold on
            end
        end
    end
    hold off
    title(['Pair ', hpair{np}, ', window length = ', num2str(twin), ' s'])
    ylabel('TDOA')
    datetick
end
% test whether label color legend exists
figCol = findall(0, 'Type', 'figure', 'name', 'Legend of Label Colors');
if isempty(figCol)
    generateColorSchemeLegend(brushing) % if legend doesn't exist, generate it
end



saveas(fig, saveFileName, 'fig')
% save(['clickByClick_roughTDOA_', trackName], 'TDet', 'TDOA', 'XCTpk', 'label')
save(saveFileName, 'whale')



%%
function [tdoa, xcpk, snr] = calcTDOA_CTC(xct, maxLag, xcovInd, maxNumTDOA, Nhann, fsct)
% [tdoa, xcpk, snr] = calcTDOA_CTC(xctCorr, xcovInd, maxNumTDOA, fsct)
% calculates the tdoa, peak of xcorr (xcpk), and SNR of a click train
% inputs:
% xct is the click train
% maxLag is the maximum lag of the xcorr
% xcovInd is the index of the cross-correlation needed for this hydrophone pair
% maxNumTDOA is the maximum number of TDOAs to save for this detection
% Nhann is the length of the window used in place of the clicks
% fsct is the sampling rate of the click train

% initialize outputs:
tdoa = -99.*ones(1, length(xcovInd), maxNumTDOA); % initialize with -99 so I can easily remove faulty TDOAs
xcpk = zeros(1, length(xcovInd), maxNumTDOA);
snr = xcpk;

[xctCorr, lags] = xcorr(xct, maxLag);

for ipair = 1:length(xcovInd)
    if max(xctCorr(:, xcovInd(ipair)))>200

        [pks, locs] = findpeaks(xctCorr(:, xcovInd(ipair)), 'SortStr', 'descend', 'minPeakHeight', 200, 'NPeaks', maxNumTDOA);
        if length(pks)>1
            if pks(1)>.8*pks(2)
                npksBig = find(pks>.8*pks(1)); % find indices of peaks bigger than .8 of max peak
                tdoa(1, ipair, 1:max(npksBig)) = lags(locs(npksBig))/fsct;
                xcpk(1, ipair, 1:max(npksBig)) = pks(npksBig);

                % calculate SNR:
                snr = zeros(size(xcpk));
                for ipk = 1:length(npksBig)
                    noiseInd = 1:length(xctCorr); % indices of 'noise'
                    detInd = locs(npksBig(ipk));
                    sigInd = max([1, detInd-Nhann/2]):min([length(xctCorr), detInd+Nhann/2]); % indices of signal
                    noiseInd(sigInd) = []; % remove indices of signal from noise

                    sigPow = mean(xctCorr(sigInd, xcovInd(ipair)).^2);
                    noisePow =  mean(xctCorr(noiseInd, xcovInd(ipair)).^2);

                    snr(1, ipair, ipk) = sigPow/noisePow;
                end
            else
                tdoa(1, ipair, 1) = lags(locs(1))./fsct;
                xcpk(1, ipair, 1) = pks(1);

                % calculate SNR:
                noiseInd = 1:length(xctCorr); % indices of 'noise'
                detInd = locs(1);
                sigInd = max([1, detInd-Nhann/2]):min([length(xctCorr), detInd+Nhann/2]); % indices of signal
                noiseInd(sigInd) = []; % remove indices of signal from noise

                sigPow = mean(xctCorr(sigInd).^2);
                noisePow =  mean(xctCorr(noiseInd).^2);

                snr(1, ipair, 1) = sigPow/noisePow;
            end
        elseif length(pks)==1
            tdoa(1, ipair, 1) = lags(locs)./fsct;
            xcpk(1, ipair, 1) = pks;

            % calculate SNR:
            noiseInd = 1:length(xctCorr); % indices of 'noise'
            detInd = locs;
            sigInd = max([1, detInd-Nhann/2]):min([length(xctCorr), detInd+Nhann/2]); % indices of signal
            noiseInd(sigInd) = []; % remove indices of signal from noise

            sigPow = mean(xctCorr(sigInd).^2);
            noisePow =  mean(xctCorr(noiseInd).^2);

            snr(1, ipair, 1) = sigPow/noisePow;
        end

    end
end
end
