clearvars; clc; close all;
do_gifs = 1;

heavi = @(x,a,b) (a-b)*(x > 0)+b;
% ----------- simulate love waves. Need a cylindrical domain... -----------------
par5.speeds = [0.5, 1]; % top is face 1, bottom is face 2
par5.k = 1;
par5.v = 0.63429190138857311;
% par5.v = 0.50;
par5.hboom = 1;
par5.depth1 = 1; % depth of the interface
par5.omega1 = par5.k*sqrt((par5.v/par5.speeds(1))^2-1);
par5.omega2 = par5.k*sqrt(1-(par5.v/par5.speeds(2))^2);
par5.depthf = @(z) (z > 0).*cos(par5.omega1*(z-par5.depth1)) + ...
    (z <= 0).*cos(par5.omega1*par5.depth1).*exp(par5.omega2*z);
par5.initf  = @(x,y) cos(par5.k*x).*par5.depthf(y);
par5.initfp = @(x,y) -par5.k*par5.v*sin(par5.k*x).*par5.depthf(y);
par5.speedf = @(y) heavi(y, par5.speeds(1)^2, par5.speeds(2)^2);
par5.tfinal = 4.5;
par5.depth2 = 5; % size of the second domain
par5.width  = 2*pi;
par5.reltol = 1.0e-3; par5.abstol = 1.0e-6; % time tolerances
% par5.step = 0.0075; % refinement in space
par5.step = 0.01;
par5.nonlinearity = @(x) 0;


par = par5;
[t,x,y,u] = periodic_lovewave( ...
    par.speedf, ...
    par.tfinal, ...
    par.depth1, par.depth2, par.width, ...
    par.initf, par.initfp, ...
    par.reltol, par.abstol, ...
    par.step, ...
    par.nonlinearity);

%% visualizations
zmax = par.hboom;
Ninter = floor(par.depth1/par.step);
mindelay = 0.015;
utop   = u(:,1,:);
uinter = u(:,Ninter,:);
min_prom = 0.05;
% res = [3,3];


if do_gifs
figure
set(gcf,'Position',[100, 100, 1500, 1000])
delete surf.gif
gif('surf.gif','nodither','DelayTime',mindelay);
delta_t = 0.1; tlast = -100;
nrm = 0;
for k = 1:length(t)
    if t(k) - tlast > delta_t
        tlast = t(k);
        surf(x,y(1:Ninter),squeeze(u(k,1:Ninter,:)),EdgeColor='none');
        hold on
        surf(x,y(Ninter:end-nrm),squeeze(u(k,Ninter:end-nrm,:)),EdgeColor='none',FaceAlpha=0.8);
        hold off
        % pcolor(x,y,squeeze(u(i,:,:)));
        colorbar
        % view(0,90)
        %c = max(abs(min(u(k,:,:),[],"all")), abs(max(u(k,:,:),[],"all")));
        %clim([-c,c]);
        %zlim([-zmax/2,zmax]);
        gif
    end
end

% plots of top of domain and interface
figure
delete interface_vs_top.gif
gif('interface_vs_top.gif','nodither','DelayTime',mindelay);
delta_t = 0.1; tlast = -100;
finterp = @(x,y,xs) spline(x,y,xs);
xref = linspace(x(1),x(end),750*length(x)); % refined xinputs
peak_top = 0;
peak_inter = 0;
tvals = 0;
ii = 0;
for k = 1:length(t)
    if t(k) - tlast > delta_t
        ii = ii + 1;
        tlast = t(k);
        tvals(ii) = tlast;

        TF = islocalmax(finterp(x,utop(k,:),xref),'MinProminence',min_prom);
        if sum(TF) ~= 0
            indices = 1:length(TF); I = indices(TF);
            peak_top(ii) = xref(I(1));
        else
            peak_top(ii) = 0;
        end
        TF = islocalmax(finterp(x,uinter(k,:),xref),'MinProminence',min_prom);
        if sum(TF) ~= 0
            indices = 1:length(TF); I = indices(TF);
            peak_inter(ii) = xref(I(1));
        else
            peak_inter(ii) = 0;
        end

        subplot(2,1,1);
        % plot(x,utop(k,:));
        plot(x,utop(k,:),peak_top(ii),finterp(x,utop(k,:),peak_top(ii)),"x");
        title("u(x,y,t) at the top of the domain");

        subplot(2,1,2);
        % plot(x,uinter(k,:));
        plot(x,uinter(k,:),peak_inter(ii),finterp(x,uinter(k,:),peak_inter(ii)),"x");
        title("u(x,y,t) at the interface of the domain");
        gif
    end
end

% tee = [tvals(2),tvals(end)];
% c = par.speeds;
% 
% figure
% subplot(2,1,1);
% plot(tvals(2:end), abs(diff(peak_top)./diff(tvals)), ...
%     tee, c(1)*[1,1], tee, c(2)*[1,1]);
% ylim([0,max(c)+0.3])
% title("Wave speed at top over time")
% legend("Wave speed","speed in top domain","speed in bottom domain")
% subplot(2,1,2);
% plot(tvals(2:end), abs(diff(peak_inter)./diff(tvals)),...
%     tee, c(1)*[1,1], tee, c(2)*[1,1]);
% ylim([0,max(c)+0.3])
% title("Wave speed at interface over time")
end

% Track the maximum of a piecewise spline of the data
finterp = @(x,y,xs) spline(x,y,xs);
xref = linspace(x(1),x(end),750*length(x)); % refined xinputs

peak_top   = zeros(1,length(t));
peak_inter = zeros(1,length(t));
for ii = 1:length(t)
    TF = islocalmax(finterp(x,utop(ii,:),xref),'MinProminence',min_prom);
    if sum(TF) ~= 0
        indices = 1:length(TF); I = indices(TF);
        peak_top(ii) = xref(I(1));
    end
    TF = islocalmax(finterp(x,uinter(ii,:),xref),'MinProminence',min_prom);
    if sum(TF) ~= 0
        indices = 1:length(TF); I = indices(TF);
        peak_inter(ii) = xref(I(1));
    end
end

tee = [t(2),t(end)];
c = par.speeds;

figure
subplot(2,1,1);
plot(t(2:end), abs(diff(peak_top)./diff(t')), tee, c(1)*[1,1], tee, c(2)*[1,1]);
ylim([0,max(c)+0.3])
title("Wave speed at top over time")
subplot(2,1,2);
plot(t(2:end), abs(diff(peak_inter)./diff(t')), tee, c(1)*[1,1], tee, c(2)*[1,1]);
ylim([0,max(c)+0.3])
title("Wave speed at interface over time")
legend("Wave speed","speed in top domain","speed in bottom domain",'Location','northwest')

delete wavespeed_time.png
exportgraphics(gcf,'wavespeed_time.png', 'Resolution', 1000)

