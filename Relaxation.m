clearvars -except drive_path; close all; clc;

%% This plots relaxation plots (and a lot of other stuff)
% This is originally a Dasha script that I just shoehorned in CFD data (via
% phaseArrData.mat). Updates to the original script did occour to change
% the way the cross-correlation was being applied (it is now localized and
% the mean is removed before computation)

%% Methodology - UWT: OTS: sigma, lambda vs x
sz_fnt = 12;
U = '12ms';
figure(1);
hold on
h = gca; colors = h.ColorOrder;
col_num = 6;
base_col = colors(1,:)-0.8*colors(1,:);
colors_WT = ([linspace(colors(1,1),base_col(1),col_num)', linspace(colors(1,2),base_col(2),col_num)', linspace(colors(1,3),base_col(3),col_num)']);
colors_WT(4,:) = colors_WT(3,:);
colors_WT(2,:) = colors_WT(1,:);

lin_type = {'-','--',':'};
symb = 'ods';
D = 0.762;

figure(2);
hold on
plot([0 1.5],[1 1],'r:','linewidth',2,'handlevisibility','off');
for ii = 2:2:6
    switch ii
        case 1; f = 0.5; A = 'A';
        case 2; f = 1.0; A = 'B';
        case 3; f = 1.5; A = 'C';
        case 4; f = 2.0; A = 'D';
        case 5; f = 2.5; A = 'E';
        case 6; f = 3.0; A = 'F';
    end
    ind = [1:17 27:38]; % 17 33 38
    tau = 1/f;
    % cd 'C:\Users\dglou\OneDrive - UCB-O365\Research\Projects\WT_Characterization\OJWT_Characterization\Hotwire_Measurements'
    file1 = ['Experimental/Hotwire_Measurements\OJWT_Unsteady_' U '_X1_' A '.tdms'];
    file2 = ['Experimental/Hotwire_Measurements\OJWT_Unsteady_' U '_X2_' A '.tdms'];
    file3 = ['Experimental/Hotwire_Measurements\OJWT_Unsteady_' U '_X3_' A '.tdms'];
    Ncyc = 33;
    f_WT = 2000;
    Nwt = floor(f_WT*tau);

    [U2_X1,U3_X1,~,~,x_X1,z_X1,D,~] = fcn_OJWT_Unsteady_Analysis(file1,'1');
    [U2_X2,U3_X2,~,~,x_X2,z_X2,~,~] = fcn_OJWT_Unsteady_Analysis(file2,'2');
    [U2_X3,U3_X3,~,~,x_X3,z_X3,~,~] = fcn_OJWT_Unsteady_Analysis(file3,'3');
    U2 = [U2_X1; U2_X2; U2_X3]; x = [x_X1, x_X2, x_X3]; clear U2_X1 U2_X2 U2_X3 x_X1 x_X2 x_X3
    U3 = [U3_X1; U3_X2; U3_X3]; z = [z_X1, z_X2, z_X3]; clear U3_X1 U3_X2 U3_X3 z_X1 z_X2 z_X3
    clear file1 file2 file3

    U2_mean = mean(U2);
    t_WT = (1:length(U2_mean))/f_WT;
    U2_ideal = mean(U2_mean) * (1 + (max(U2_mean)-min(U2_mean))/2/mean(U2_mean)*sin(2*pi*f*t_WT));
    [C,lag] = xcorr(U2_ideal,U2_mean);
    ind_shift1 = lag(C==max(C));
    U2_mean = circshift(U2_mean,ind_shift1);
    U2_cyc = reshape(U2_mean(1:Ncyc*Nwt),Nwt,Ncyc);
    U2_mncyc = mean(U2_cyc(:,7:end-3),2);
    U_mean2 = (max(U2_mncyc)+min(U2_mncyc))/2;
    sig2 = (max(U2_mncyc)-min(U2_mncyc))/U_mean2/2;

    clear U3_m sig_x ind_shiftp k
    ind_shift1 = 0;
    t.tau = (1:length(U2_mncyc))/length(U2_mncyc);
    k = 0;
    theta_o = 0;
    for i = 1:length(ind)
        %         U3(ind(i),:) = circshift(U3(ind(i),:),ind_shift1);
        ind_U = floor(tau*f_WT*2):length(U3(ind(i),:));

        figure(10+ii+1);
        if i == 2
            subplot(2,1,1); plot(t_WT/tau,U3(ind(i),:),'k'); hold on;
        elseif i == length(ind)
            subplot(2,1,2); plot(t_WT/tau,U3(ind(i),:),'k'); hold on;
        end

        U3_m(i) = (max(U3(ind(i),ind_U))+min(U3(ind(i),ind_U)))/2;
        sig_x(i) = (max(U3(ind(i),ind_U))-min(U3(ind(i),ind_U)))/(2*U3_m(i));
        %         u_xt = U3_m(i) * (1 + sig_x(i) * sin( -2*pi*f*t_WT + theta_o)); dx = x(ind(i));
        u_xt = U3(ind(1),:); dx = x(ind(i));
        %         if i < 2; u_xt = U3(ind(1),:); else; u_xt = U3(ind(i-2),:); end; dx = x(ind(i))-x(ind(i-1));
        if i == 1
            % [C,lag] = xcorr(U3(ind(i),ind_U),u_xt(ind_U));
            [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),u_xt(ind_U)-mean(u_xt(ind_U)));
            ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
            dx = x(ind(i));
        elseif i == length(ind)
            % [C,lag] = xcorr(U3(ind(i),ind_U),U3(ind(i-1),ind_U));
            [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),U3(ind(i-1),ind_U)-mean(U3(ind(i-1),ind_U)));  
            ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
            dx = x(ind(i))-x(ind(i-1));
        elseif i < 3 || i >  length(ind)-2
            % [C,lag] = xcorr(U3(ind(i),ind_U),U3(ind(i-1),ind_U));
            [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),U3(ind(i-1),ind_U)-mean(U3(ind(i-1),ind_U)));
            ind_shiftp1 = lag(C==max(C));  % shift 2nd signal to phase of first
            % [C,lag] = xcorr(U3(ind(i+1),ind_U),U3(ind(i),ind_U));
            [C,lag] = xcorr(U3(ind(i+1),ind_U)-mean(U3(ind(i+1),ind_U)),U3(ind(i),ind_U)-mean(U3(ind(i),ind_U))); 
            ind_shiftp2 = lag(C==max(C));  % shift 2nd signal to phase of first
            ind_shiftp(i) = ind_shiftp2+ind_shiftp1;  % shift 2nd signal to phase of first
            dx = x(ind(i+1))-x(ind(i-1));
        else
            % [C,lag] = xcorr(U3(ind(i),ind_U),U3(ind(i-2),ind_U));
            [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),U3(ind(i-2),ind_U)-mean(U3(ind(i-2),ind_U)));
            ind_shiftp1 = lag(C==max(C));  % shift 2nd signal to phase of first
            % [C,lag] = xcorr(U3(ind(i+2),ind_U),U3(ind(i),ind_U)); 
            [C,lag] = xcorr(U3(ind(i+2),ind_U)-mean(U3(ind(i+2),ind_U)),U3(ind(i),ind_U)-mean(U3(ind(i),ind_U))); 
            ind_shiftp2 = lag(C==max(C));  % shift 2nd signal to phase of first
            ind_shiftp(i) = ind_shiftp2+ind_shiftp1;  % shift 2nd signal to phase of first
            dx = x(ind(i+2))-x(ind(i-2));
        % else
        %     [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),u_xt(ind_U)-mean(u_xt(ind_U)));
        %     ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
        end
        % [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),u_xt(ind_U)-mean(u_xt(ind_U))); % 1st, 2nd
        % ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
        if i == 1
            theta_o = ind_shiftp(i)/(f_WT*tau)*2*pi;
        else
            k(i) = ind_shiftp(i)/(f_WT*tau)*2*pi/dx;
        end
        %         u_xt = U3_m(i) * (1 + sig_x(i) * sin( k(i)*x(ind(i)) - 2*pi*f*t_WT + theta_o));
        u_xt = circshift(u_xt,ind_shiftp(i));
        if (i == 2) || (i == length(ind))
            plot(t_WT/tau,u_xt,'k--');
        end
        axis([8 10 8 18]); xlabel('$t / \tau$');
        ylabel('$u$ [m/s]'); grid on; %title(['$x/D = ' num2str(x(ind(i))/D) '$']);
    end
    figure(1);
    hold on
    plot(x(ind)/D,(sig_x), '--', 'color',colors_WT(ii,:), 'linewidth',2);
    set(gca,'fontsize',sz_fnt);

    lambda = 2*pi./k;
    figure(2);
    hold on
    plot(x(ind)/D,(lambda)*f./U3_m, '--', 'color',colors_WT(ii,:), 'linewidth',2);
    set(gca,'fontsize',sz_fnt);

    figure(3)
    plot(nan,nan, '--', 'color',colors_WT(ii,:), 'linewidth',2); hold on;

    figure(4)
    hold on 
    plot(x(ind)/D,U3_m, '--', 'color',colors_WT(ii,:), 'linewidth',2);

end


% cd(['C:\Users\dglou\OneDrive - UCB-O365\Research\Dissertation Defense\figs']);
% print(figure(1),['sig_lambda_OTS.png'],'-dpng','-r300');

%% Repeat for computational data:

% Load data
% Probably smart to write out my phase arr and spatial data and just
% load them here to save time
% Make some dummy data to play with
% n = 30;
% close all
% phaseArr = ones(1800, 3);
% uBar = 12;
f = 3;
% omega = 2*pi*f;
% sigma = linspace(.2, .3, n);
% kappa = logspace(-1, 5, n);
% lambda_ref = 2*pi./kappa;
%
%
% D = 0.762;


% x = linspace(0, 0.5, n);
%
% figure(100)
% hold on
% for i = 1:n
%     for j = 1:1800
%        phaseArr(j, i) = uBar*(1+sigma(i)*sin(kappa(i)*x(i) - omega*t_WT(j)));
%     end
%     plot(linspace(0, 1, 1800), phaseArr(:, i))
% end
%
% figure(200)
% hold on
% for i = 1:100:900
%     plot(x, phaseArr(i, :))
% end
clear ind_U U3_m sig_x u_xt ind_shiftp k u_xt
load("PhaseArrData.mat")
ts = 0.000185377509506;
f_WT = 1/ts;

tau = 1/f;
U3 = uArr';
x = vLoc;
% t_WT = 0:ts:(ts*(size(U3, 2)-1));
t_WT = (1:size(U3, 2))/f_WT;

ind = 1:size(U3, 1);
ind_shift1 = 0;
k = 0;
theta_o = 0;
for i = 1:length(ind)
    %         U3(ind(i),:) = circshift(U3(ind(i),:),ind_shift1);
    ind_U = 1:length(U3(ind(i),:));

    figure(600+1);
    if i == 2
        subplot(2,1,1); plot(t_WT/tau,U3(ind(i),:),'k'); hold on;
    elseif i == length(ind)
        subplot(2,1,2); plot(t_WT/tau,U3(ind(i),:),'k'); hold on;
    end

    U3_m(i) = (max(U3(ind(i),ind_U))+min(U3(ind(i),ind_U)))/2;
    sig_x(i) = (max(U3(ind(i),ind_U))-min(U3(ind(i),ind_U)))/(2*U3_m(i));
    %         u_xt = U3_m(i) * (1 + sig_x(i) * sin( -2*pi*f*t_WT + theta_o)); dx = x(ind(i));
    u_xt = U3(ind(1),:); dx = x(ind(i));
    %         if i < 2; u_xt = U3(ind(1),:); else; u_xt = U3(ind(i-2),:); end; dx = x(ind(i))-x(ind(i-1));
        if i == 1
            % [C,lag] = xcorr(U3(ind(i),ind_U),u_xt(ind_U));
            [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),u_xt(ind_U)-mean(u_xt(ind_U)));
            ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
            dx = x(ind(i));
        elseif i == length(ind)
            % [C,lag] = xcorr(U3(ind(i),ind_U),U3(ind(i-1),ind_U));
            [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),U3(ind(i-1),ind_U)-mean(U3(ind(i-1),ind_U)));  
            ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
            dx = x(ind(i))-x(ind(i-1));
        else
            % [C,lag] = xcorr(U3(ind(i),ind_U),U3(ind(i-1),ind_U));
            [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),U3(ind(i-1),ind_U)-mean(U3(ind(i-1),ind_U)));
            ind_shiftp1 = lag(C==max(C));  % shift 2nd signal to phase of first
            % [C,lag] = xcorr(U3(ind(i),ind_U),U3(ind(i-1),ind_U));
            [C,lag] = xcorr(U3(ind(i+1),ind_U)-mean(U3(ind(i+1),ind_U)),U3(ind(i),ind_U)-mean(U3(ind(i),ind_U))); 
            ind_shiftp2 = lag(C==max(C));  % shift 2nd signal to phase of first
            ind_shiftp(i) = ind_shiftp2+ind_shiftp1;  % shift 2nd signal to phase of first
            dx = x(ind(i+1))-x(ind(i-1));
        % else
        %     [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),u_xt(ind_U)-mean(u_xt(ind_U)));
        %     ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
        end
        % [C,lag] = xcorr(U3(ind(i),ind_U)-mean(U3(ind(i),ind_U)),u_xt(ind_U)-mean(u_xt(ind_U))); % 1st, 2nd
        % ind_shiftp(i) = lag(C==max(C));  % shift 2nd signal to phase of first
         if i == 1
        theta_o = ind_shiftp(i)/(f_WT*tau)*2*pi;
    else
        k(i) = ind_shiftp(i)/(f_WT*tau)*2*pi/dx;
    end
    %         u_xt = U3_m(i) * (1 + sig_x(i) * sin( k(i)*x(ind(i)) - 2*pi*f*t_WT + theta_o));
    u_xt = circshift(u_xt,ind_shiftp(i));
    if (i == 2) || (i == length(ind))
        plot(t_WT/tau,u_xt,'k--');
    end
    axis([8 10 8 18]); xlabel('$t / \tau$');
    ylabel('$u$ [m/s]'); grid on; %title(['$x/D = ' num2str(x(ind(i))/D) '$']);
end
uMean = mean(uArr, 1);
uMean = (max(uArr, [], 1)+min(uArr, [], 1))./2;
figure(1)
plot(x/D,(sig_x), 'k', 'linewidth',2);
% set(gca,'fontsize',sz_fnt);

lambda = 2*pi./(k);
lambda(1) = nan;
figure(2)
plot(x/D,(lambda)*f./uMean, 'k', 'linewidth',2);
% set(gca,'fontsize',sz_fnt);

figure(4)
plot(x/D, uMean, 'k', 'linewidth',2);

figure(3)
plot(nan,nan, 'k', 'linewidth',2);

figure(1)
set(gca, "FontSize", 12)
xlabel('$x/D$','interpreter','latex');
ylabel('Amplitude, $\sigma (x)$','interpreter','latex');
% title('Amplitude', FontSize=10)
axis([0 1.5 0 0.5]);
xlim([0, 1.15])
grid on;

figure(2)
set(gca, "FontSize", 12)
xlabel('$x/D$','interpreter','latex');
ylabel('Normalized Convective Speed, $\lambda(x) f / \overline{u_1} $','interpreter','latex');
% title('Normalized Convective Speed', FontSize=10)
% yticks(0:1:20); yticklabels({'0','','2','','4','','6','','8','','10'});
axis([0 1.15 0 4]);
% xlim([0, 1.15])
grid on;

% ltr = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
% for jj = 1:2
%     subplot(1,2,jj);
%     fg1 = gca;
%     if jj == 1
%         fg1.Position(1) = fg1.Position(1) - 0.0250; pause(0.1);
%     else
%         fg1.Position(1) = fg1.Position(1) + 0.0250; pause(0.1);
%     end
%     fg1.Position(4) = fg1.Position(4) - 0.075; pause(0.1);
%     % annotation('textbox',[fg1.Position(3)/2 + fg1.Position(1)-0.02   0 0.1 0.1],'String',ltr{jj},'EdgeColor','none','FontSize', sz_fnt);
% end

figure(1)
l = legend('$f = 1Hz$ Exp.','$f = 2Hz$ Exp.','$f = 3Hz$ Exp.', '$f = 3Hz$ CFD', 'interpreter','latex','location','northwest');
% l.Box = 'off';
% l.Position
% plt = gcf;
% l.Position([1 2]) = [0.3 0.92];
% l.ItemTokenSize = [15 18];
legend;
set(gca, "FontSize", 12)
fi = gcf;
exportgraphics(fi,'../NoWing_Production/SavedPlots/AmplvsSpace.png','Resolution',400)

figure(2)
l = legend('$f = 1Hz$ Exp.','$f = 2Hz$ Exp.','$f = 3Hz$ Exp.', '$f = 3Hz$ CFD', 'interpreter','latex','location','northeast');
title('Mean removed, local xcorr')
% l.Box = 'off';
% l.Position
% plt = gcf;
% l.Position([1 2]) = [0.3 0.92];
% l.ItemTokenSize = [15 18];
legend;
set(gca, "FontSize", 12)
% set(gcf,'Position',(get(l,'Position')...
%     .*[0, 0, 1, 1].*get(gcf,'Position')));
% set(l,'Position',[0,0,1,1]);
% set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);
fi = gcf;
exportgraphics(fi,'../NoWing_Production/SavedPlots/NormConvSpeed.png','Resolution',400)


figure(4)
l = legend('$f = 1Hz$ Exp.','$f = 2Hz$ Exp.','$f = 3Hz$ Exp.', '$f = 3Hz$ CFD', 'interpreter','latex','location','northeast');
% l.Box = 'off';
% l.Position
% plt = gcf;
% l.Position([1 2]) = [0.3 0.92];
% l.ItemTokenSize = [15 18];
legend(Location="south");
set(gca, "FontSize", 12)


%%
function [U_f52, U_f53, U_C2, U_C3, x, z, D, vane] = fcn_OJWT_Unsteady_Analysis(file,X)
% cd 'C:\Users\dglou\OneDrive - UCB-O365\Research\Projects\WT_Characterization\OJWT_Characterization\Hotwire_Measurements'
data_WT = TDMS_getStruct(file);
D = 0.762;
N = (length(fields(data_WT))-1);
clear rho U RPM dP P T Vel
for i = 1:N
    subfile = ['DataPt_' num2str(i)];
    RPM(i)  = data_WT.(subfile).MotorSpeed.Props.Mean;
    T(i)    = data_WT.(subfile).TotalTemperature.Props.Mean;
    rho(i)  = data_WT.(subfile).Density.Props.Mean;
    U(i,:)    = data_WT.(subfile).Velocity.data;
    vane(i,:) = data_WT.(subfile).VaneAngle.data;
    dP(i,:)   = data_WT.(subfile).DifferentialPressure.data;
    P(i,:)    = data_WT.(subfile).StaticPressure.data;
    RH(i,:)   = data_WT.(subfile).RelativeHumidity.data;
    HW1(i,:)  = data_WT.(subfile).Hotwire1.data;
    HW2(i,:)  = data_WT.(subfile).Hotwire2.data;
    HW3(i,:)  = data_WT.(subfile).Hotwire3.data;
    x(i)    = data_WT.(subfile).X.Props.Mean/1000;
    z(i)    = data_WT.(subfile).Z.Props.Mean/1000;
end
% Calculate viscosity using Sutherland's law
for j = 1:length(T)
    [mu(j)] = Sutherland(T(j));     % dynamic viscosity [kg/m-s]
    % !!!! should be taken using T_film = 1/2(Ta+Tw) !!!!
end
nu = mu./rho;       % kinematic viscosity [m2/s]



if X==2
    file_cal2 = ['Experimental/Hotwire_Calibrations\Hotwire_Calibration_2021_09_14_HW2_X2.tdms'];
    % [date_data] = fcn_date_cal(data_WT);
    date_data = '2021_09_14'
    file_cal3 = ['Experimental/Hotwire_Calibrations\Hotwire_Calibration_' date_data '.tdms'];
else
    % [date_data] = fcn_date_cal(data_WT);
    date_data = '2021_09_14'
    file_cal2 = ['Experimental/Hotwire_Calibrations\Hotwire_Calibration_' date_data '.tdms'];
    file_cal3 = ['Experimental/Hotwire_Calibrations\Hotwire_Calibration_' date_data '.tdms'];
end

%     cd 'Q:\ZOC_Test\Hotwire_Calibrations'
% cd 'C:\Users\dglou\OneDrive - UCB-O365\Research\Projects\AFOSR_Gust_Wing_Project\Pressure_Wing_Data\Hotwire_Calibrations'
data_CL2 = TDMS_getStruct(file_cal2);
data_CL3 = TDMS_getStruct(file_cal3);
ind_CL = 2:(length(fields(data_CL2))-1);
for i_CL = 1:length(ind_CL)
    subfile = ['DataPt_' num2str(ind_CL(i_CL))];
    rho2(i_CL)   = data_CL2.(subfile).Density.Props.Mean;
    U_cal2(i_CL) = data_CL2.(subfile).Velocity.Props.Mean;
    es2(i_CL)   = data_CL2.(subfile).Hotwire2.Props.Mean;
    T2(i_CL)     = data_CL2.(subfile).TotalTemperature.Props.Mean;
end
% [f52,B2,nu2,k_a2] = fcn_HWA_Fit(es2,U_cal2,rho2,T2,'5th'); %'Cim'
f52 = polyfit(es2,U_cal2,5);
ind_CL = 2:(length(fields(data_CL3))-1);
for i_CL = 1:length(ind_CL)
    subfile = ['DataPt_' num2str(ind_CL(i_CL))];
    rho3(i_CL)   = data_CL3.(subfile).Density.Props.Mean;
    U_cal3(i_CL) = data_CL3.(subfile).Velocity.Props.Mean;
    es3(i_CL)   = data_CL3.(subfile).Hotwire3.Props.Mean;
    T3(i_CL)     = data_CL3.(subfile).TotalTemperature.Props.Mean;
end
% [f53,B3,nu3,k_a3] = fcn_HWA_Fit(es3,U_cal3,rho3,T3,'5th'); %'Cim'
f53 = polyfit(es3,U_cal3,5);

es_ex = -10:0.1:10;
%     figure();
%     subplot(1,2,1); plot(es_ex, polyval(f52,es_ex),'b'); hold on; plot(es2,U_cal2,'kx')
% %         plot(es2,nu2.*((B2(1).*es2+B2(3)).^2./(k_a2*(B2(5)-T3))+B2(2)).^B2(4),'r');
%         title('HW2'); grid on; xlabel('E [V]'); ylabel('U [m/s]'); %ylim([0 35]);
%     subplot(1,2,2); plot(-10:0.1:10, polyval(f53,-10:0.1:10),'b'); hold on; plot(es3,U_cal3,'kx');
%         title('HW3'); grid on; xlabel('E [V]'); ylabel('U [m/s]'); %ylim([0 35]);


U_f52 = polyval(f52,HW2);
U_f53 = polyval(f53,HW3);

for i = 1:N
    % U_C2(i,:)  = nu(i).*((B2(1).*HW2(i,:)+B2(3)).^2./(k_a2*(B2(5)-T(i)))+B2(2)).^B2(4);
    % U_C3(i,:)  = nu(i).*((B3(1).*HW3(i,:)+B3(3)).^2./(k_a3*(B3(5)-T(i)))+B3(2)).^B3(4);
end
U_C2 = 0;
U_C3 = 0;
end

function [mu] = Sutherland(T)
mu0 = 1.716e-5;
T0 = 273;
S = 111;

mu = (T/T0)^(3/2)*(T0+S)/(T+S)*mu0;
end