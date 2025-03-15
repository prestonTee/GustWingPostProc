function [tPhase, uPhase, pPhase, uMed] = processVartData(vDat, tOverTau, iShift, nPhases, nPerPhase)
% processVartData Pulls u and p data from raw vart data
%   This fucntion also reshaped the data for the purposes of the unsteady
%   OTS simulations. The period and nubmer of cycles must be passed to the
%   function for it to behave correctly

% vDat [2D array]   ->  Contains the raw vart data from extractVartinfo()
% tOverTau          ->  Time vector that measures location in the cycle from the reference point use to set the velocity signal phase
% iShift            ->  Indexes to shift the data by to get to the a phase start
% nPhases           ->  Number of phases to keep
% nPerPhase         ->  Number of timesteps per phase

%% Code
%Arrays for ease of use
uInp = vDat(:, 3);
pInp = vDat(:, 2);
uMed = (max(uInp)+min(uInp))/2;
%     uNorm = uInp/((max(uInp)+min(uInp))/2);
% ======= Changing the denominator here is the easiest way to achive a vertical shift ======
%     uNorm = uInp/U_mean; %Dasha had more data going into the mean so use her's as we cannot reproduce
% modTTau = mod(tOverTau, 1);

% uArr(:, i) = uInp;

uDat = reshape(uInp(iShift:iShift+nPhases*nPerPhase-1), [nPerPhase, nPhases])';
pDat = reshape(pInp(iShift:iShift+nPhases*nPerPhase-1), [nPerPhase, nPhases])';

%Pick the good phases of data (uncommoent figure 5 for raw data)
%     tPhase = mod(tOverTau(1801:3600), 1);
% tPhase = mod(tOverTau(i_shift:i_shift+nPerPhase-1), 1);
tPhase = mod(tOverTau(iShift:iShift+nPerPhase-1), 1);
uPhase = mean(uDat)';
uPhaseSTD = std(uDat)';
pPhase = mean(pDat)';
pPhaseSTD = std(pDat)';

% uPhaseArr(:, i) = uPhase;
% dat = table(tPhase, uPhase, pPhase);
% d.(['Pos', num2str(i)]) = dat;

%     save(['CompCenterlineData_', nameList{i}, '.mat'], 'dat', 'ts')

% Rearrange data to match time
[tPhase, I] = sort(tPhase);

uPhase = uPhase(I);
pPhase = pPhase(I);

uDat = uDat(:,I);
pDat = pDat(:,I);

uNorm = uDat/uMed;

%% Mean filtering
% This is for pressure only
% I think this was filtering added by Ken
pPhaseMod = [pPhase(end-49:end); pPhase; pPhase(1:50)];
windowSize = 10; 
 b = (1/windowSize)*ones(1,windowSize);
 a = 1;
 y = filter(b,a,pPhaseMod);
 pPhaseMod = filter(b,a,y);
 y = filter(b,a,pPhaseMod);
 pPhaseMod = filter(b,a,y);
 y = filter(b,a,pPhaseMod);
 pPhaseMod = filter(b,a,y);
 y = filter(b,a,pPhaseMod);
pPhase = y(51:end-50);

%% Signal Filtering
pDatF = pInp;
windowSize = 10; 
 b = (1/windowSize)*ones(1,windowSize);
 a = 1;
 y = filter(b,a,pDatF);
 pDatF = filter(b,a,y);
 y = filter(b,a,pDatF);
 pDatF = filter(b,a,y);
 y = filter(b,a,pDatF);
 pDatF = filter(b,a,y);
 y = filter(b,a,pDatF);
pDatF = reshape(y(iShift:iShift+nPhases*nPerPhase-1), [nPerPhase, nPhases])';
end

