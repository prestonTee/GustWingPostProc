clear; clc; close all;

%% Used for making compaion plots of centerline velocity
% Note that this is all done automatically and switching happens with
% simple toggles of boolians at the top of this script

% This was also used with vertical lines to make animations at one point

%% Wind Tunnel Geometry
D = 0.762; % [m]

%% Setup - Comp
% Input selections of what to plot here:
i2 = false;
i3 = true;
i3s = false;
iExp3 = true;

% Choose linestyles
l2 = ':';
l3 = '-';
l3s = '--';
lExp3 = '--';

% --- Below this should not need to be touched ---
% vFileLoc = '/nobackup/uncompressed/Models/GustWing/OTS/FullRoomAgain/PostProc/v1_NoWing/Unsteady/varts/';
vFileLoc2 = 'Unsteady2Hz/varts/';
vFileLoc3 = 'Unsteady3Hz/varts/';
vFileLoc3s = 'Unsteady3Hz_S/varts/';
dFileLocExp3 = 'Experimental/';
% vFileLoc = '/nobackup/uncompressed/Models/GustWing/OTS/Coflow/PostProc/NoWing/Unsteady/Varts/';
vPre = 'varts.';
vSuffix = '.dat';
locFileName = 'xyzts_noheader.dat';

% Setup stuff
delta_ts = 1; % Increase to skip timesteps in CFD data
lineVec = {l2, l3, l3s, lExp3};
nameVec = {'i2', 'i3', 'i3s', 'iExp3'};
titleVec = {'2Hz', 'CFD 3Hz', 'Simple 3Hz', 'Experimental 3Hz'};
plotVec = [i2, i3, i3s, iExp3];
iTrue = find(plotVec); % Selects only those set to true above

figure(1);
hold on
% % clf;
set(gcf,'units','in','position',[1 1 12 6]);



% Our global color vector
ColorVecG = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0 0 0]};

% For all those set to true above
for k = iTrue
    % Pull out some data about the current selction 
    config = nameVec{k};
    n = titleVec{k};
    line = lineVec{k};
    if strcmp('i2', config) % Check if this is the 2Hz case and set things 
        P = 1/6.279*pi;
        ts = 1.85313511962e-4;
        
        numTsPerP = 2700;
        
        nPhases = 40;
        vFileLoc = vFileLoc2;
    elseif strcmp('iExp3', config) %Check if we are plotting experimental data
        % Nothing to do at the moment

    else % Otherwise this is one of the 3Hz cases
        P = 1/9.415*pi;
        ts = 1.85377509506e-4;
        
        numTsPerP = 1800;
        if strcmp('i3', config)
            nPhases = 40;
            vFileLoc = vFileLoc3;
            ColorVec = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0 0 0]};

        elseif strcmp('i3s', config)
            nPhases = 40;
            vFileLoc = vFileLoc3s;
            ColorVec = {[0.8500 0.3250 0.0980], [0.6350 0.0780 0.1840], [0 0 0]};
        end
    end
    
    figure(10+k)
    hold on
    set(gcf,'units','in','position',[1 1 12 6]);
    

    % New loop here for experimental data pulling and plotting
    if strcmp('iExp3', config)
        load([dFileLocExp3, 'OTS_3Hz_12ms.mat']);

        % I believe this came from Dasha's code
        ind = [11, 30];% 0, D/2, D (zero implied)
        U2_mean = mean(U2);
        U_mean = mean(U2_mean);
        U2_mncyc = mean(reshape(U2_mean(1:Ncyc*Nwt),Nwt,Ncyc),2); % Seems to be phase averaging
        ind_shift = find(U2_mncyc==max(U2_mncyc))-floor(Nwt/4); % Shifting, looks like defineing the peak as pi/2
        U2_mncyc = circshift(U2_mncyc,-ind_shift);

        figure(1);
        plot((1:length(U2_mncyc))/length(U2_mncyc),U2_mncyc,'--','color',ColorVec{1},'linewidth',2.5);
        col_num = length(ind)+1;
        for i = 1:(col_num-1)
            U3_mncyc(i,:) = mean(reshape(U3(ind(i),1:Ncyc*Nwt),Nwt,Ncyc),2);
            U3_mncyc(i,:) = circshift(U3_mncyc(i,:),-ind_shift);
            t_tau = (1:length(U2_mncyc))/length(U2_mncyc);
            % For normalized:
            % plot(t_tau,U3_mncyc(i,:)/U_mean,'color',colors_WT(i+1,:),'linewidth',1.5);
            % For absolute:
            plot(t_tau,U3_mncyc(i,:),'--','color',ColorVecG{i+1},'linewidth',2.5);
        end

        figure(10+k)
        plot((1:length(U2_mncyc))/length(U2_mncyc),U2_mncyc,'--','color',ColorVec{1},'linewidth',2.5);
        for i = 1:(col_num-1)
            U3_mncyc(i,:) = mean(reshape(U3(ind(i),1:Ncyc*Nwt),Nwt,Ncyc),2);
            U3_mncyc(i,:) = circshift(U3_mncyc(i,:),-ind_shift);
            t_tau = (1:length(U2_mncyc))/length(U2_mncyc);
            % For normalized:
            % plot(t_tau,U3_mncyc(i,:)/U_mean,'color',colors_WT(i+1,:),'linewidth',1.5);
            % For absolute:
            plot(t_tau,U3_mncyc(i,:),'--','color',ColorVecG{i+1},'linewidth',2.5);
        end
        title(n)

    else
        dat = load([vFileLoc, locFileName]); % Positional data
        datNormByD = dat./D;
        
        % vList = [1, 22, 43, 64]; % = [~0, D/4, D/2, D] along centerline
        vList = [1, 43, 64]; % = [~0, D/2, D] along centerline
        vLoc =  dat(vList, 1);
        
        figure(10+k)
        hold on
        % % clf;
        set(gcf,'units','in','position',[10 5 12 6]);
        
        
        ColorVecG = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0 0 0]};
            
        %% ////////////////////////////////////////////////////////////////////////
        % Computational Plotting
        %  ////////////////////////////////////////////////////////////////////////
        %% My old Plotting
        
        nameList = {'zero', 'dOver4', 'dOver2', 'd'};      
        
        vPrefix = ['varts.', num2str(vList(1)), '.'];
        vDat = extractVartInfo(vFileLoc, vPrefix, vSuffix, delta_ts);
        vDat = vDat(1:end, :);
        uInp = vDat(:, 3);
        
        %Time data
        tsDat = vDat(:, 1);
        t = tsDat*ts;
        nPerPhase = ceil(P/ts);
        
        % Match Dasha's shift
        % ======= Changing the denominator (s) here is the easiest way to achive a phase shift ======
        s = 3.65;
        if strcmp('i3s', config)
            s = 4;
        end
        i_shift = find(uInp==max(uInp))-floor(nPerPhase/s);
        % iShift = mod(i_shift, numTsPerP);
        if strcmp('i3', config)
            iShift = mod(i_shift, numTsPerP) + numTsPerP; % Last term determines how many cycles to throw away
        else
            iShift = mod(i_shift, numTsPerP) + numTsPerP*2; 
        end
        
        % Relocate t=0 to match dasha's data
        t = t - t(i_shift);
        tOverTau = t/P;
        tsStart = tsDat(i_shift)-numTsPerP*23;
        
        
        for i = 1:length(vList)
            % Read the data
            vPrefix = ['varts.', num2str(vList(i)), '.'];
            vDat = extractVartInfo(vFileLoc, vPrefix, vSuffix, delta_ts);
            vDat = vDat(1:end, :);

            [tPhase, uPhase, pPhase, uMean, uPhaseSTD] = processVartData(vDat, tOverTau, iShift, nPhases, nPerPhase);
            % uNorm = uDat/uMed;
            
            %% Plot it
            %     color = [0.49 0.18 0.55] - ([0.49 0.18 0.55]*dList(i));
            figure(1)
            %     plot(tPhase, uPhase, '--', 'Color', color, 'LineWidth',1.5)
            plot(tPhase, uPhase, line, 'Color', ColorVecG{i}, 'LineWidth',k+1)
            fill([tPhase; flip(tPhase)], [uPhase+uPhaseSTD; flip(uPhase-uPhaseSTD)], ColorVecG{i}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')

            figure(10+k)
            %     plot(tPhase, uPhase, '--', 'Color', color, 'LineWidth',1.5)
            plot(tPhase, uPhase, line, 'Color', ColorVecG{i}, 'LineWidth',4)
            fill([tPhase; flip(tPhase)], [uPhase+uPhaseSTD; flip(uPhase-uPhaseSTD)], ColorVecG{i}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
            title(n)
        end        
    end
end