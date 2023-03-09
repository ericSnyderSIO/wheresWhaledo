function whaleOut = localize(whaleIn, hloc, H1, H2, driftPoly, varargin)
% whaleOut = localize(whaleIn, hloc, H1, H2, driftPoly)
% whaleOut = localize(whaleIn, hloc, H1, H2, driftPoly, paramFile)
% takes the TDOAs from a "whale" struct (output of calcTDOAfromCTC) and
% performs iterative model generation to determine the most likely whale location.
% INPUTS:
% -whaleIn = the "whale" struct
% -hloc = hydrophone locations (m).
% -H1 and H2 = H matrices for 4ch instruments (m).
% -driftPoly{i} contains the polynomials for the relative drift between
%  instruments. If drift is linear, driftPoly{i} will be two values
%  (y=mx+b).

global LOC

if nargin==6
    loadParams(varargin{1})
elseif nargin==5
    loadParams('localize.params')
else
    fprintf('Incorrect number of inputs')
end

% Overall process:
% 1) Iterate through each detection
%    a) Use LMSE to determine most likely whale location based on TDOA in
%    courseModel 

whaleOut = whaleIn;
Nlrg = numel(driftPoly); % number of large ap TDOAs

sig2sml = LOC.sig_sml^2; % variance of small ap
sig2lrg = LOC.sig_lrg^2; % variance of large ap

for wn = 1:numel(whaleOut)
    if ~isempty(whaleOut{wn})
        % initialize output loc files
        whaleOut{wn}.wlocCoarse = nan(length(whaleOut{wn}.TDet), 3);
        whaleOut{wn}.wloc = nan(length(whaleOut{wn}.TDet), 3);  
        whaleOut{wn}.CIx = nan(length(whaleOut{wn}.TDet), 2);  
        whaleOut{wn}.CIy = whaleOut{wn}.CIx;
        whaleOut{wn}.CIz = whaleOut{wn}.CIx;

        % want either both four-channels or one four-channel and two or
        % more large ap TDOAs
        detUse = find(sum(~isnan(whaleOut{wn}.TDOA(:, 13:end)),2)>1 | sum(~isnan(whaleOut{wn}.TDOA(:, 1:12)),2)==12);

        for ndet = 1:length(detUse)
            detInd = detUse(ndet); % index of detection currently being used

            % coarse localization with model:
            
            TDOA = whaleOut{wn}.TDOA(detInd, :);
            drift = zeros(1, Nlrg);
            for ntdoa = 1:Nlrg
                drift(ntdoa) = polyval(driftPoly{ntdoa}, whaleOut{wn}.TDet(detInd));
            end
            TDOA(13:end) = TDOA(13:end) + LOC.driftSign.*drift;

            Isml = find(~isnan(TDOA(1:12))); % indices of small ap used
            Ilrg = find(~isnan(TDOA(13:end)))+12; % indices of large ap used
            
            Asml = (2*pi*sig2sml)^(-length(Isml)/2); % coefficient of small ap
            Alrg = (2*pi*sig2lrg)^(-length(Ilrg)/2); % coefficient of large ap
            
            xv = linspace(LOC.M{1}.x(1), LOC.M{1}.x(2), LOC.M{1}.Nx);
            yv = linspace(LOC.M{1}.y(1), LOC.M{1}.y(2), LOC.M{1}.Ny);
            zv = linspace(LOC.M{1}.z(1), LOC.M{1}.z(2), LOC.M{1}.Nz);
             
            for niter = 1:LOC.niter-1 
                % iteratively generate model and perform max likelihood
                [mTDOA, mwloc] = makeModel(xv, yv, zv, hloc, H1, H2, LOC.c);

                Lsml = Asml*exp(-1./(2.*sig2sml).*sum((mTDOA(:,1:12)-TDOA(1:12)).^2, 2, 'omitnan'));
                Llrg = Alrg*exp(-1./(2.*sig2lrg).*sum((mTDOA(:,13:end)-TDOA(13:end)).^2, 2, 'omitnan'));

                L = Lsml.*Llrg;
                [~, I] = max(L);
                wlocEst = mwloc(I, :);
                xv = linspace(LOC.M{niter+1}.x(1) + wlocEst(1), LOC.M{niter+1}.x(2) + wlocEst(1), LOC.M{niter+1}.Nx);
                yv = linspace(LOC.M{niter+1}.y(1) + wlocEst(2), LOC.M{niter+1}.y(2) + wlocEst(2), LOC.M{niter+1}.Ny);
                zv = linspace(LOC.M{niter+1}.z(1) + wlocEst(3), LOC.M{niter+1}.z(2) + wlocEst(3), LOC.M{niter+1}.Nz);

            end

            % final iteration (outside for-loop to avoid unecessary
            % calculations)

            [mTDOA, mwloc] = makeModel(xv, yv, zv, hloc, H1, H2, LOC.c);

            Lsml = Asml*exp(-1./(2.*sig2sml).*sum((mTDOA(:,1:12)-TDOA(1:12)).^2, 2, 'omitnan'));
            Llrg = Alrg*exp(-1./(2.*sig2lrg).*sum((mTDOA(:,13:end)-TDOA(13:end)).^2, 2, 'omitnan'));

            L = Lsml.*Llrg;
            [~, Ibest] = max(L);
            whaleOut{wn}.wloc(detInd, :) = mwloc(Ibest, :);

            % calculate confidence intervals:
            [CIx, CIy, CIz] = calcCI(TDOA, whaleOut{wn}.wloc(detInd, :), hloc, H1, H2, Asml, Alrg, LOC);
            whaleOut{wn}.CIx(detInd, :) = CIx;
            whaleOut{wn}.CIy(detInd, :) = CIy;
            whaleOut{wn}.CIz(detInd, :) = CIz;

        end
    end
end




end

%%
function [TDOA, wloc] = makeModel(xv, yv, zv, h, H1, H2, c)

% make wloc (matrix of whale positions)
[cx, cy, cz] = ndgrid(xv, yv, zv);
wloc = [cx(:), cy(:), cz(:)];

s1 = wloc-h(1, :);
r1 = sqrt(sum(s1.^2, 2)); % range to instrument 1
s1 = s1./r1; % direction vector to instrument 1

s2 = wloc-h(2, :);
r2 = sqrt(sum(s2.^2, 2)); % range to instrument 2
s2 = s2./r2; % direction vector to instrument 2

s3 = wloc-h(3, :);
r3 = sqrt(sum(s3.^2, 2)); % range to instrument 3

s4 = wloc-h(4, :);
r4 = sqrt(sum(s4.^2, 2)); % range to instrument 4

% small aperture TDOAs
TDOA(:, 1:6) = (s1*H1.')./c;
TDOA(:, 7:12) = (s2*H2.')./c;

% large aperture TDOAs
TDOA(:, 13) = (r1-r2)./c;
TDOA(:, 14) = (r1-r3)./c;
TDOA(:, 15) = (r1-r4)./c;
TDOA(:, 16) = (r2-r3)./c;
TDOA(:, 17) = (r2-r4)./c;
TDOA(:, 18) = (r3-r4)./c;

end

%%
function [CIx, CIy, CIz] = calcCI(TDOA, wloc, h, H1, H2, Asml, Alrg, LOC)
CIx = nan(1, 2);
CIy = nan(1, 2);
CIz = nan(1, 2);
% vectors used for CI calculations:
xv = LOC.M{1}.x(1):LOC.M{1}.x(2);
yv = LOC.M{1}.y(1):LOC.M{1}.y(2);
zv = LOC.M{1}.z(1):LOC.M{1}.z(2);

% calculate CIx:
[mTDOA, mwloc] = makeModel(xv, wloc(2), wloc(3), h, H1, H2, LOC.c);
Lsml = Asml*exp(-1./(2.*LOC.sig_sml^2).*sum((mTDOA(:,1:12)-TDOA(1:12)).^2, 2, 'omitnan'));
Llrg = Alrg*exp(-1./(2.*LOC.sig_lrg^2).*sum((mTDOA(:,13:end)-TDOA(13:end)).^2, 2, 'omitnan'));
Lx = Lsml.*Llrg;
% Cx = cumsum(Lx)./sum(Lx);
Cx = cumsum(Lx);
% Cx = Cx-min(Cx);
Cx = Cx./max(Cx);
I = find(Cx<=.025, 1, 'last');
if ~isempty(I)
    CIx(1) = xv(I);
else
    CIx(1) = LOC.wxlim(1);
end
I = find(Cx>=.975, 1, 'first');
if ~isempty(I)
    CIx(2) = xv(I);
else
    CIx(2) = LOC.wxlim(2);
end

% calculate CIy:
[mTDOA, mwloc] = makeModel(wloc(1), yv, wloc(3), h, H1, H2, LOC.c);
Lsml = Asml*exp(-1./(2.*LOC.sig_sml^2).*sum((mTDOA(:,1:12)-TDOA(1:12)).^2, 2, 'omitnan'));
Llrg = Alrg*exp(-1./(2.*LOC.sig_lrg^2).*sum((mTDOA(:,13:end)-TDOA(13:end)).^2, 2, 'omitnan'));
Ly = Lsml.*Llrg;
Cy = cumsum(Ly);
% Cy = Cy - min(Cy);
Cy = Cy./max(Cy);
I = find(Cy<=.025, 1, 'last');
if ~isempty(I)
    CIy(1) = yv(I);
else
    CIy(1) = LOC.wylim(1);
end
I = find(Cy>=.975, 1, 'first');
if ~isempty(I)
    CIy(2) = yv(I);
else
    CIy(2) = LOC.wylim(2);
end

% calculate CIz:
[mTDOA, mwloc] = makeModel(wloc(1), wloc(2), zv, h, H1, H2, LOC.c);
Lsml = Asml*exp(-1./(2.*LOC.sig_sml^2).*sum((mTDOA(:,1:12)-TDOA(1:12)).^2, 2, 'omitnan'));
Llrg = Alrg*exp(-1./(2.*LOC.sig_lrg^2).*sum((mTDOA(:,13:end)-TDOA(13:end)).^2, 2, 'omitnan'));
Lz = Lsml.*Llrg;
Cz = cumsum(Lz);
% Cz = Cz-min(Cz);
Cz = Cz./max(Cz);
I = find(Cz<=.025, 1, 'last');
if ~isempty(I)
    CIz(1) = zv(I);
else
    CIz(1) = LOC.wzlim(1);
end
I = find(Cz>=.975, 1, 'first');
if ~isempty(I)
    CIz(2) = zv(I);
else
    CIz(2) = LOC.wzlim(2);
end

end