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

%%
% If an instrument has -99 in the TDOA then don't use that instrument in the final localization

% Localize on 100 m grid
load('clickByClick_roughTDOA_180611_1030.mat')
M200 = load('D:\SOCAL_E_63\tracking\experiments\largeApertureTDOA\TDOAmodel_100m')

% set up parameters for LMS (CHECK WITH JOHN ABOUT THESE)
sigmaH_sml = .1e-3; % uncertainty in small ap hydrophone locations
sigmaX_sml = .05e-3; % uncertainty in small ap TDOA
sig_sml = sqrt(sigmaH_sml^2 + sigmaX_sml^2);

sigmaH_lrg = 50e-3; % uncertainty in large ap hydrophone locations
sigmaX_lrg = 30e-3; % uncertainty in large ap TDOA
sig_lrg = sqrt(sigmaH_lrg^2 + sigmaX_lrg^2);

parr = 2; % principal array, aka which one was used as basis of ctc above

iloc = 0;
for idet = 1:length(TDet)

    tdoaToUse = 1:6;

    % determine which large ap pairs to use:
    Irem = find(TDOA(idet, :)==-99);
    if length(Irem)==6 % insufficient good TDOAs
        continue
    else
        iloc = iloc+1;
        if ~isempty(Irem)
            M = 6 - length(Irem); % number of large ap TDOAs to use
            tdoaToUse(Irem) = []; % remove large ap tdoas
            Llrg = (2*pi*sig_lrg^2)^(-M/2).*exp((-1/(2*sig_lrg^2)).*sum((M200.TDOA(:, tdoaToUse + 12) - TDOA(idet, tdoaToUse)).^2, 2));
        else
            M = 6; % number of large ap TDOAs to use
            Llrg = (2*pi*sig_lrg^2)^(-M/2).*exp((-1/(2*sig_lrg^2)).*sum((M200.TDOA(:, 13:18) - TDOA(idet, :)).^2, 2));
        end

        % determine which small ap pairs to use:
        if TDOA(idet, 1)==-99 % if pair 1-2 is a faulty TDOA
            N = 6; % number of small aperture TDOAs used
            if parr==1
                [~, ind] = min(abs(DET{1}.TDet-TDet(idet)));
                tdoasml = -DET{1}.TDOA(ind, :);
                Lsml = (2*pi*sig_sml^2)^(-N/2).*exp((-1/(2*sig_sml^2)).*sum((M200.TDOA(:, 1:6) - tdoasml).^2, 2));
            elseif parr==2
                [~, ind] = min(abs(DET{2}.TDet-TDet(idet)));
                tdoasml = -DET{2}.TDOA(ind, :);
                Lsml = (2*pi*sig_sml^2)^(-N/2).*exp((-1/(2*sig_sml^2)).*sum((M200.TDOA(:, 7:12) - tdoasml).^2, 2));
            end
        else % use both small ap
            N = 12; % number of small aperture TDOAs used
            if parr==1
                [~, ind1] = min(abs(DET{1}.TDet-TDet(idet)));
                [~, ind2] = min(abs(DET{2}.TDet-(TDet(idet) - TDOA(1)/spd)));
                tdoasml = -[DET{1}.TDOA(ind1, :), DET{2}.TDOA(ind2, :)];
                Lsml = (2*pi*sig_sml^2)^(-N/2).*exp((-1/(2*sig_sml^2)).*sum((M200.TDOA(:, 1:12) - tdoasml).^2, 2));
            elseif parr==2
                [~, ind1] = min(abs(DET{1}.TDet-(TDet(idet) + TDOA(1)/spd)));
                [~, ind2] = min(abs(DET{2}.TDet-TDet(idet)));
                tdoasml = -[DET{1}.TDOA(ind1, :), DET{2}.TDOA(ind2, :)];
                Lsml = (2*pi*sig_sml^2)^(-N/2).*exp((-1/(2*sig_sml^2)).*sum((M200.TDOA(:, 1:12) - tdoasml).^2, 2));
            end
        end
        L = Lsml.*Llrg;

        [~, im] = max(L);

        wloc_course(iloc, :) = M200.wloc(im, :);
        locLabel(iloc) = label(idet);
    end
end

figure(101)
scatter3(wloc_course(:, 1), wloc_course(:, 2), wloc_course(:, 3), 33, brushing.params.colorMat(locLabel, :), 'filled')

%% WELP That didn't work. Let's try using expected TDOA and aligning clicks?

% Steps: 
% 1) iterate through each detection on the primary array. 
% 2) Take all detections on other arrays within the window of the click train
% 3) Determine if a click in primary array aligns with each click in click
% train of other array
% 4) if a primary array click is close enough to a click in the examined
% array, assign it the TDOA value determined by the click train
% 5) after a click on the examined array has been assigned all possible
% TDOA values from all the click trains in which it appeared, do some sort
% calculation to determine most likely true TDOA
% 6) Extract time frame around click from most likely TDOA and xcorr for
% more precise TDOA

Ndet = length(TDet); % number of detections on primary array
NTDOA = 6; % number of TDOAs (no. of hydrophone pairs)

for idet = 1:Ndet
    
end