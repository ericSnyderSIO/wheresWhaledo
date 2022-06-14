% step through 30 seconds at a time, brushdet, validate detections, click train corr, coarse localization
load('detections_brushDOA180611_1030')
global brushing
loadParams('brushing.params')
spd = 24*60*60;
fsct = 100e3;
maxLag = round(fsct*2000/1500);
twin = 30; % window length for click train
N = 1024;
Wk = hann(N/4);
xcovInd = [2,3,4,7,8,12];

xct = zeros(twin*fsct, 4);

encStart = min([DET{1}.TDet(1), DET{2}.TDet(1)]);
encEnd = min([DET{1}.TDet(end), DET{2}.TDet(end)]);

arrno = 2; % which array's labels are being used
TDet = [];
TDOA = [];
detCount = 0;
for wn = unique(DET{arrno}.color).'
    Ilab = find(DET{arrno}.color==wn);

    for ndet = 1:length(Ilab); % iterate through each detection
        tdet = DET{arrno}.TDet(Ilab(ndet));

        tstart = tdet - twin/(spd*2);
        tend = tdet + twin/(spd*2);

        tct = tstart:1/(spd*fsct):tend;
        xct = zeros(length(tct), 4);

        I = find(DET{arrno}.TDet>=tstart & DET{arrno}.TDet<=tend & DET{arrno}.color==wn);
        for i = 1:length(I)
            [~, ind] = min(abs(tct-DET{arrno}.TDet(I(i))));
            xct(ind, arrno) = 1;
        end

        xctHann(:, arrno) = conv(Wk, xct(:, arrno));

        if arrno==1
            otherArrays = [2,3,4];
        elseif arrno==2
            otherArrays = [1,3,4];
        end

        for ia = 1:3
            I = find(DET{otherArrays(ia)}.TDet>=tstart & DET{otherArrays(ia)}.TDet<=tend);
            for i = 1:length(I)
                [~, ind] = min(abs(tct-DET{otherArrays(ia)}.TDet(I(i))));
                xct(ind, otherArrays(ia)) = 1;
            end
            xctHann(:, otherArrays(ia)) = conv(Wk, xct(:, otherArrays(ia)));
        end
        if 0 % plot or not
            figure(1)
            for ih = 1:4
                subplot(4,1,ih)
                plot(xctHann(:, ih))
            end
        end
        [xctCorr, lags] = xcorr(xctHann, maxLag);

        tdoa = -99.*ones(1,length(xcovInd)); % initialize with -99 so I can easily remove faulty TDOAs

        for np = 1:length(xcovInd)
            [pk, imax] = max(xctCorr(:, xcovInd(np)));
            if pk>150
                tdoa(:, np) = lags(imax)./fsct;
            end
        end
        detCount = detCount + 1;
        TDet(detCount) = tdet;
        TDOA(detCount, :) = tdoa;
        label(detCount) = wn;

    end

end
%%
fig = figure(3);
for np = 1:6
    subplot(6, 1, np)
    Iplt = find(TDOA(:, np)~=-99);
    scatter(TDet(Iplt), TDOA(Iplt, np), 20, brushing.params.colorMat(label(Iplt), :), 'filled')
    title(['Pair ', num2str(np), ', window length = ', num2str(twin), ' s'])
    ylabel('TDOA')
    datetick
end


saveas(fig, ['clickByClick_windowLength_', num2str(twin)], 'fig')