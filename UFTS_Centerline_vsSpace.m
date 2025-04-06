clear; clc; close all;

%% Used for u vs x plots - currently CFD only

%% Wind Tunnel Geometry
D = 0.762; % [m]

%% Setup - Comp
% Non-comparitive selectior
% Choose 3 or 2 Hz first, next selector only matters if 3Hz
f = 3; % Hz
simple = true;

if f == 2
    vFileLoc = 'Unsteady2Hz/varts/';
    P = 1/6.279*pi;
    ts = 1.85313511962e-4;
    caseName = '2Hz';
elseif f == 3 && ~simple
    vFileLoc = 'Unsteady3Hz/varts/';
    P = 1/9.415*pi;
    ts = 1.85377509506e-4;
    caseName = 'Full 3Hz';
elseif f == 3 && simple
    vFileLoc = 'Unsteady3Hz_S/varts/';
    P = 1/9.415*pi;
    ts = 1.85377509506e-4;
    caseName = 'Simple 3Hz';
end
delta_ts = 1;

vPre = 'varts.';
vSuffix = '.dat';
locFileName = 'xyzts_noheader.dat';

dat = load([vFileLoc, locFileName]); % Positional data
datNormByD = dat./D;

vList = [1, 85:107]; % along centerline
vLoc =  dat(vList, 1);


%% ////////////////////////////////////////////////////////////////////////
% Computational Plotting
%  ////////////////////////////////////////////////////////////////////////
%% My old Plotting

vPrefix = ['varts.', num2str(vList(1)), '.'];
vDat = extractVartInfo(vFileLoc, vPrefix, vSuffix, delta_ts);
vDat = vDat(1:end, :);
uInp = vDat(:, 3);

%Time data
tsDat = vDat(:, 1);
t = tsDat*ts;
nPerPhase = ceil(P/ts);

%Match Dasha's shift
% ======= Changing the denominator here is the easiest way to achive a phase shift ======
if simple 
    i_shift = find(uInp==max(uInp))-floor(nPerPhase/4);
else
    i_shift = find(uInp==max(uInp))-floor(nPerPhase/3.65);
end
% iShift = mod(i_shift, 1800);
if f == 3 && ~simple
    iShift = mod(i_shift, nPerPhase) + nPerPhase; % Extra 180 throws away the first cycle
else
    iShift = mod(i_shift, nPerPhase) + nPerPhase*2; % Extra 180 throws away the first cycle
end

% Relocate t=0 to match dasha's data
t = t - t(i_shift);
tOverTau = t/P;

uPhaseArr = nan(nPerPhase, length(vList));
uArr = nan(length(uInp), length(vList));
for i = 1:length(vList)
    % Read the data
    vPrefix = ['varts.', num2str(vList(i)), '.'];
    vDat = extractVartInfo(vFileLoc, vPrefix, vSuffix, delta_ts);
    vDat = vDat(1:end, :);
    nPhases = 40;

[tPhase, uPhase, pPhase, uMed] = processVartData(vDat, tOverTau, iShift, nPhases, nPerPhase);

uPhaseArr(:, i) = uPhase;
pPhaseArr(:, i) = pPhase;
    
%     color = [0.49 0.18 0.55] - ([0.49 0.18 0.55]*dList(i));
end

%% 
dVLoc = diff(vLoc);
dudxPhaseArr = diff(uPhaseArr, 1, 2)/dVLoc(1);
dpdxPhaseArr = diff(pPhaseArr, 1, 2)/dVLoc(1);

vLocDiff = vLoc(1:end-1) + dVLoc/2;

%%
% save(['PhaseArrData', num2str(f), 'Hz.mat'], 'uArr', 'uPhaseArr', 'vLoc')
locPhaseAcc = [.7, .75, .8, .85, .9, .95, 0, .05, .1, .15, .2, .25];
iPhaseAcc = locPhaseAcc*nPerPhase +1;
locPhaseDec = [.3, .35, .4, .45, .5, .55, .6, .65];
iPhaseDec = floor(locPhaseDec*nPerPhase +1);

figure(10)
hold on
grid on 
for i = 1:length(iPhaseAcc)
    name = ['$t/\tau = $', num2str(locPhaseAcc(i))];
    % color = [0 0.4470 0.7410] - ([0 0.4470 0.7410]*((i-1)/(length(iPhaseAcc)-1))); % Blue to black
    color = [0.1914 0.6367 0.3281] - ([0.1914 0.6367 0.3281]*((i-1)/(length(iPhaseAcc)-1))); % Green to black
    % color = [0.4660 0.6740 0.1880] - ([0.4660 0.6740 0.1880]*((i-1)/(length(iPhaseAcc)-1))); % Green (MATLAB) to black
    plot(vLoc/D, uPhaseArr(iPhaseAcc(i), :), 'DisplayName', name, 'color', color, 'linewidth', 2.5)
end
xlabel('Normalized Downstream Distance, x/D', 'Interpreter', 'latex')
ylabel('Phase Averaged Velocity, $\langle u \rangle$ [m/s]', 'Interpreter', 'latex')
set(gca, "FontSize", 12)
% title('Velocity Profiles: Acceleration', 'FontSize', 10)
legend('location','eastoutside', 'Interpreter', 'latex')
% set(gca, "FontSize", 12)
f = gcf;
% exportgraphics(f,'SavedPlots/UvsX_Acel.png','Resolution',400)


figure(11)
hold on
grid on
for i = 1:length(iPhaseDec)
    name = ['$t/\tau = $', num2str(locPhaseDec(i))];
    % color = [0 0.4470 0.7410] - ([0 0.4470 0.7410]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Blue
    % color = [0.8984 0.3320 0.0508] - ([0.8984 0.3320 0.0508]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Orange
    color = [0.8500 0.3250 0.0980] - ([0.8500 0.3250 0.0980]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Orange (MATLAB)
    plot(vLoc/D, uPhaseArr(iPhaseDec(i), :), 'DisplayName', name, 'color', color, 'linewidth', 2.5)
end
xlabel('Normalized Downstream Distance, x/D', 'Interpreter', 'latex')
ylabel('Phase Averaged Velocity, $\langle u \rangle$ [m/s]', 'Interpreter', 'latex')
set(gca, "FontSize", 12)
% title('Velocity Profiles: Deceleration', 'FontSize', 10)
legend('location','eastoutside', 'Interpreter', 'latex')
% set(gca, "FontSize", 12)
f = gcf;
% exportgraphics(f,'SavedPlots/UvsX_Decc.png','Resolution',400)

%%
figure(20)
hold on
grid on 
for i = 1:length(iPhaseAcc)
    name = ['$t/\tau = $', num2str(locPhaseAcc(i))];
    % color = [0 0.4470 0.7410] - ([0 0.4470 0.7410]*((i-1)/(length(iPhaseAcc)-1))); % Blue to black
    color = [0.1914 0.6367 0.3281] - ([0.1914 0.6367 0.3281]*((i-1)/(length(iPhaseAcc)-1))); % Green to black
    % color = [0.4660 0.6740 0.1880] - ([0.4660 0.6740 0.1880]*((i-1)/(length(iPhaseAcc)-1))); % Green (MATLAB) to black
    plot(vLocDiff/D, dudxPhaseArr(iPhaseAcc(i), :), '--', 'DisplayName', name, 'color', color, 'linewidth', 2.5)
end
xlabel('Normalized Downstream Distance, x/D', 'Interpreter', 'latex')
ylabel('$\langle \frac{\partial u}{\partial x} \rangle$ [m/s]', 'Interpreter', 'latex')
set(gca, "FontSize", 12)
% title('Velocity Profiles: Acceleration', 'FontSize', 10)
legend('location','eastoutside', 'Interpreter', 'latex')
ylim([-10, 10])
% set(gca, "FontSize", 12)
f = gcf;
% exportgraphics(f,'SavedPlots/UvsX_Acel.png','Resolution',400)


figure(21)
hold on
grid on
for i = 1:length(iPhaseDec)
    name = ['$t/\tau = $', num2str(locPhaseDec(i))];
    % color = [0 0.4470 0.7410] - ([0 0.4470 0.7410]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Blue
    % color = [0.8984 0.3320 0.0508] - ([0.8984 0.3320 0.0508]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Orange
    color = [0.8500 0.3250 0.0980] - ([0.8500 0.3250 0.0980]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Orange (MATLAB)
    plot(vLocDiff/D, dudxPhaseArr(iPhaseDec(i), :), '--', 'DisplayName', name, 'color', color, 'linewidth', 2.5)
end
xlabel('Normalized Downstream Distance, x/D', 'Interpreter', 'latex')
ylabel('$\langle \frac{\partial u}{\partial x} \rangle$ [m/s]', 'Interpreter', 'latex')
set(gca, "FontSize", 12)
% title('Velocity Profiles: Deceleration', 'FontSize', 10)
legend('location','eastoutside', 'Interpreter', 'latex')
ylim([-10,10])
% set(gca, "FontSize", 12)
f = gcf;
% exportgraphics(f,'SavedPlots/UvsX_Decc.png','Resolution',400)

figure(30)
hold on
grid on 
for i = 1:length(iPhaseAcc)
    name = ['$t/\tau = $', num2str(locPhaseAcc(i))];
    % color = [0 0.4470 0.7410] - ([0 0.4470 0.7410]*((i-1)/(length(iPhaseAcc)-1))); % Blue to black
    color = [0.1914 0.6367 0.3281] - ([0.1914 0.6367 0.3281]*((i-1)/(length(iPhaseAcc)-1))); % Green to black
    % color = [0.4660 0.6740 0.1880] - ([0.4660 0.6740 0.1880]*((i-1)/(length(iPhaseAcc)-1))); % Green (MATLAB) to black
    plot(vLocDiff/D, dpdxPhaseArr(iPhaseAcc(i), :), '-.', 'DisplayName', name, 'color', color, 'linewidth', 2.5)
end
xlabel('Normalized Downstream Distance, x/D', 'Interpreter', 'latex')
ylabel('$\langle \frac{\partial p}{\partial x} \rangle$ [m/s]', 'Interpreter', 'latex')
set(gca, "FontSize", 12)
% title('Velocity Profiles: Acceleration', 'FontSize', 10)
legend('location','eastoutside', 'Interpreter', 'latex')
ylim([-50, 40])
% set(gca, "FontSize", 12)
f = gcf;
% exportgraphics(f,'SavedPlots/UvsX_Acel.png','Resolution',400)


figure(31)
hold on
grid on
for i = 1:length(iPhaseDec)
    name = ['$t/\tau = $', num2str(locPhaseDec(i))];
    % color = [0 0.4470 0.7410] - ([0 0.4470 0.7410]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Blue
    % color = [0.8984 0.3320 0.0508] - ([0.8984 0.3320 0.0508]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Orange
    color = [0.8500 0.3250 0.0980] - ([0.8500 0.3250 0.0980]*((length(iPhaseDec)-(i-1)-1)/(length(iPhaseDec)-1))); % Black to Orange (MATLAB)
    plot(vLocDiff/D, dpdxPhaseArr(iPhaseDec(i), :), '-.', 'DisplayName', name, 'color', color, 'linewidth', 2.5)
end
xlabel('Normalized Downstream Distance, x/D', 'Interpreter', 'latex')
ylabel('$\langle \frac{\partial p}{\partial x} \rangle$ [m/s]', 'Interpreter', 'latex')
set(gca, "FontSize", 12)
% title('Velocity Profiles: Deceleration', 'FontSize', 10)
legend('location','eastoutside', 'Interpreter', 'latex')
ylim([-50, 40])
% set(gca, "FontSize", 12)
f = gcf;
% exportgraphics(f,'SavedPlots/UvsX_Decc.png','Resolution',400)


%%
% Create a custom colormap (Blue -> White -> Red)
n = 256; % Number of colors
half_n = round(n / 2);

% Define the blue to white transition
blue = [0.1, 0.15, 0.75]; % Blue
white = [1, 1, 1]; % White
red = [0.7, 0, 0];   % Red

% Interpolate colors
cmap1 = [linspace(blue(1), white(1), half_n)', ...
         linspace(blue(2), white(2), half_n)', ...
         linspace(blue(3), white(3), half_n)'];

cmap2 = [linspace(white(1), red(1), half_n)', ...
         linspace(white(2), red(2), half_n)', ...
         linspace(white(3), red(3), half_n)'];

% Combine into one colormap
custom_cmap = [cmap1; cmap2];


%%
meanU = mean(uPhaseArr, 1);
uAdj = uPhaseArr - uMed;

[X, T] = meshgrid(vLoc/D, (1:nPerPhase)/nPerPhase);
figure(100)
hold on 
grid on
colormap(custom_cmap);
levels = linspace(9, 16, 200);
contourf(X, T, uPhaseArr, levels, 'LineColor','none','HandleVisibility','off')
% set(gca, 'YScale', 'log');
xlabel('x/D')
ylabel('t/$\tau$', 'interpreter', 'latex')
c = colorbar;
c.Label.String = 'Velocity [m/s]';
% contour(X, T, uAdj, [0,0], 'k', 'DisplayName', 'Mean Velocity')
% legend()
caxis([min(uPhase),max(uPhase)]);
title(caseName)
hold off
f = gcf;
% exportgraphics(f,'SavedPlots/XvsT.png','Resolution',400)

meanU = mean(uPhaseArr, 1);
uAdj = uPhaseArr - uMed;