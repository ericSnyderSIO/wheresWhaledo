% build TDOA model

% hewLL = [32.65646  -119.48815 -1330.1631];
% heeLL = [32.65879  -119.47705 -1319.6305];
% henLL = [32.66221  -119.48424 -1321.2775];
% hesLL = [32.65352  -119.48446 -1331.3959];
% 
% hoLL = mean([hewLL; heeLL; henLL; hesLL]);
% 
% [hew(1), hew(2)] = latlon2xy(hewLL(1), hewLL(2), hoLL(1), hoLL(2));
% hew(3) = hewLL(3)-hoLL(3);
% [hee(1), hee(2)] = latlon2xy(heeLL(1), heeLL(2), hoLL(1), hoLL(2));
% hee(3) = heeLL(3)-hoLL(3);
% [hen(1), hen(2)] = latlon2xy(henLL(1), henLL(2), hoLL(1), hoLL(2));
% hen(3) = henLL(3)-hoLL(3);
% [hes(1), hes(2)] = latlon2xy(hesLL(1), hesLL(2), hoLL(1), hoLL(2));
% hes(3) = hesLL(3)-hoLL(3);

load('D:\SOCAL_E_63\xwavTables\instrumentLocs.mat')  % calculated in D:\MATLAB_addons\gitHub\wheresWhaledo\experiments\calcSigma.m

hee = [0,0,0];
hew = h(1,:);
hen = h(2,:);
hes = h(3,:);


hyd1 = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EE_Hmatrix_fromHydLocInversion_210702.mat')
hyd2 = load('D:\SOCAL_E_63\tracking\experiments\inverseProblem\matfiles\SOCAL_E_63_EW_Hmatrix_fromHydLocInversion_210702.mat')

% HEW = H;

% Reorder hydrophones to fit new TDOA order
HEE = [hyd1.hydPos(2,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(3,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(1,:);
         hyd1.hydPos(3,:)-hyd1.hydPos(2,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(2,:);
         hyd1.hydPos(4,:)-hyd1.hydPos(3,:)];
     
HEW = [hyd2.hydPos(2,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(3,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(1,:);
         hyd2.hydPos(3,:)-hyd2.hydPos(2,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(2,:);
         hyd2.hydPos(4,:)-hyd2.hydPos(3,:)];

c = 1488.4;

%% coarse resolution

coarseRes = 200; % grid size in coarse grid

xvecCoarse = -4500:coarseRes:3500;
yvecCoarse = -4000:coarseRes:4000;
zvecCoarse = -coarseRes:coarseRes:1400; % height above array

wloc = zeros(length(xvecCoarse)*length(yvecCoarse)*length(zvecCoarse), 3);
TDOA = zeros(length(xvecCoarse)*length(yvecCoarse)*length(zvecCoarse), 18);
n = 0;
for nx = 1:length(xvecCoarse)
    for ny = 1:length(yvecCoarse)
        for nz = 1:length(zvecCoarse)
            n = n+1;
            wloc(n, :) = [xvecCoarse(nx), yvecCoarse(ny), zvecCoarse(nz)];
            
            % find received time if source time is 0:
            see = wloc(n,:) - hee;
            ree = sqrt(sum(see.^2));
            see = see./ree;
            tee = ree/c;
            
            sew = wloc(n,:) - hew;
            rew = sqrt(sum(sew.^2));
            sew = sew./rew;
            tew = rew/c;
            
            sen = wloc(n,:) - hen;
            ren = sqrt(sum(sen.^2));
            ten = ren/c;
            
            ses = wloc(n,:) - hes;
            res = sqrt(sum(ses.^2));
            tes = res/c;
            
            TDOA(n, 1:6) = (HEE*see.')./c; % small aperture TDOA for EE
            TDOA(n, 7:12) = (HEW*sew.')./c; % small aperture TDOA for EW
            % large aperture TDOAs:
            TDOA(n, 13) = tee - tew;
            TDOA(n, 14) = tee - ten;
            TDOA(n, 15) = tee - tes;
            TDOA(n, 16) = tew - ten;
            TDOA(n, 17) = tew - tes;
            TDOA(n, 18) = ten - tes;
            
        end
    end
    
end


save('B:\TDOAmodel_200m', 'TDOA', 'wloc')

%% 10 m resolution

wlocCoarse = wloc;
fineRes = 10;

% create vectors for all points in grid
% ranges are 200 more than previous so that each file contains +/2 from
% 100m resolution whale loc


for nCoarse = 1:length(wlocCoarse)
        
    xvecFine = (wlocCoarse(nCoarse, 1)-2*coarseRes):fineRes:(wlocCoarse(nCoarse, 1)+2*coarseRes);
    yvecFine = (wlocCoarse(nCoarse, 2)-2*coarseRes):fineRes:(wlocCoarse(nCoarse, 2)+2*coarseRes);
    zvecFine = (wlocCoarse(nCoarse, 3)-2*coarseRes):fineRes:(wlocCoarse(nCoarse, 3)+2*coarseRes);
    
    wloc = zeros(length(xvecFine)*length(yvecFine)*length(zvecFine), 3);
    TDOA = zeros(length(xvecFine)*length(yvecFine)*length(zvecFine), 18);
    n = 0;
    
    for nx = 1:length(xvecFine)
        for ny = 1:length(yvecFine)
            for nz = 1:length(zvecFine)
                n = n+1;
                wloc(n, :) = [xvecFine(nx), yvecFine(ny), zvecFine(nz)];
                
                % find received time assuming source time is 0:
                see = wloc(n,:) - hee;
                ree = sqrt(sum(see.^2));
                see = see./ree;
                tee = ree/c; % received time EE
                
                sew = wloc(n,:) - hew;
                rew = sqrt(sum(sew.^2));
                sew = sew./rew;
                tew = rew/c;
                
                sen = wloc(n,:) - hen;
                ren = sqrt(sum(sen.^2));
                ten = ren/c;
                
                ses = wloc(n,:) - hes;
                res = sqrt(sum(ses.^2));
                tes = res/c;
                
                TDOA(n, 1:6) = (HEE*see.')./c; % small aperture TDOA for EE
                TDOA(n, 7:12) = (HEW*sew.')./c; % small aperture TDOA for EE
                % large aperture TDOAs:
                TDOA(n, 13) = tee - tew;
                TDOA(n, 14) = tee - ten;
                TDOA(n, 15) = tee - tes;
                TDOA(n, 16) = tew - ten;
                TDOA(n, 17) = tew - tes;
                TDOA(n, 18) = ten - tes;
                
            end
        end
        
    end
    
    
    save(['B:\modelFiles_10mFrom200m\TDOAmodel_10m_n=', sprintf('%05d', nCoarse)], 'TDOA', 'wloc')
end



