clearvars; clc; close all;
do_gifs = 0;

heavi = @(x,a,b) (a-b)*(x > 0)+b;
% ----------- simulate love waves. Need a cylindrical domain... -----------------
par5.speeds = [0.5, 1]; % top is face 1, bottom is face 2
par5.k = 2;
par5.v = 0.56115873910026276;
par5.amp = 1;
par5.depth1 = 1; % depth of the interface
par5.omega1 = par5.k*sqrt((par5.v/par5.speeds(1))^2-1);
par5.omega2 = par5.k*sqrt(1-(par5.v/par5.speeds(2))^2);
par5.depthf = @(z) (z > 0).*cos(par5.omega1*(z-par5.depth1)) + ...
    (z <= 0).*cos(par5.omega1*par5.depth1).*exp(par5.omega2*z);
par5.initf  = @(x,y) par5.amp*cos(par5.k*(x)).*par5.depthf(y);
par5.initfp = @(x,y) -par5.amp*par5.k*par5.v*sin(par5.k*(x)).*par5.depthf(y);
par5.speedf = @(y) heavi(y, par5.speeds(1)^2, par5.speeds(2)^2);
par5.tfinal = pi/par5.v;
par5.depth2 = 5.5; % size of the second domain
par5.width  = 2*pi;
par5.reltol = 1.0e-4; par5.abstol = 1.0e-7; % time tolerances
par5.step = 0.015; % refinement in space
% par5.step = 0.01;
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
zmax = par.amp;
Ninter = floor(par.depth1/par.step);
mindelay = 0.015;
utop   = u(:,1,:);
uinter = u(:,Ninter,:);
min_prom = 0.005;
% res = [3,3];


% Track the maximum of a spline of u on the interface and top of the domain
finterp = @(x,y,xs) spline(x,y,xs);
xref = linspace(x(1),x(end),1000*length(x)); % refined xinputs
peak_top   = zeros(1,length(t));
peak_inter = zeros(1,length(t));
for ii = 1:length(t)
    TF = islocalmax(finterp(x,utop(ii,:),xref),'MinProminence',min_prom);
    if sum(TF) ~= 0
        indices = 1:length(TF); I = indices(TF);
        peak_top(ii) = xref(I(end));
    end
    TF = islocalmax(finterp(x,uinter(ii,:),xref),'MinProminence',min_prom);
    if sum(TF) ~= 0
        indices = 1:length(TF); I = indices(TF);
        peak_inter(ii) = xref(I(end));
    end
end

%% 
if do_gifs
figure
set(gcf,'Position',[100, 100, 1500, 1000])
delete surf.gif
gif('surf.gif','nodither','DelayTime',mindelay);
delta_t = 0.1; tlast = -100;
for k = 1:length(t)
    if t(k) - tlast > delta_t
        tlast = t(k);
        surf(x,y(1:Ninter),squeeze(u(k,1:Ninter,:)),EdgeColor='none');
        hold on
        surf(x,y(Ninter:end),squeeze(u(k,Ninter:end,:)),EdgeColor='none',FaceAlpha=0.8);
        hold off
        colorbar
        zlim([-par.amp,par.amp])
        % view(0,90)
        gif
    end
end

%%
% plots of top of domain and interface
figure
delete interface_vs_top.gif
gif('interface_vs_top.gif','nodither','DelayTime',mindelay);
delta_t = 0.1; tlast = -100;
for k = 1:length(t)
    if t(k) - tlast > delta_t
        tlast = t(k);
        tiledlayout(2,1,'TileSpacing','Compact');
        nexttile
        plot(x,utop(k,:)); hold on;
        plot(peak_top(k),finterp(x,utop(k,:),peak_top(k)),"x"); hold off
        title("u(x,y,t) at the top of the domain");

        nexttile
        plot(x,uinter(k,:)); hold on;
        plot(peak_inter(k),finterp(x,uinter(k,:),peak_inter(k)),"x");
        title("u(x,y,t) at the interface of the domain");
        gif
    end
end
end

%% Heatmaps
figure
Nfigures = 5;
set(gcf,'Position',[100, 100, Nfigures*170, 160])
times = floor(linspace(1,length(t),Nfigures));
tiled_guy = tiledlayout(1,Nfigures,'TileSpacing','Compact','Padding','Compact');
for k = times
    nexttile
    surf(x,y(1:Ninter),squeeze(u(k,1:Ninter,:)),EdgeColor='none');
    hold on
    surf(x,y(Ninter:end),squeeze(u(k,Ninter:end,:)),EdgeColor='none',FaceAlpha=0.8);
    hold off
    view(0,90)
    xlim([x(1),x(end)])
    ylim([y(end),y(1)])
    title('t=' + string(t(k)))
end
%title(tiled_guy,"Vertical displacement u(x,z,t) over time")
xlabel(tiled_guy,"x")
ylabel(tiled_guy,"z")
cb = colorbar;
cb.Layout.Tile = 'east';

delete overhead.eps
print('overhead.eps','-depsc2','-r400');

%% 1D waves
figure;
set(gcf,'Position',[450 458 700 350])
formatSpec = '%.2f';
Nfigures = 5;
times = floor(linspace(1,length(t),Nfigures));
tiled_guy = tiledlayout(2,Nfigures,'TileSpacing','Compact','Padding','Compact');
for k = times
    nexttile
    plot(x,utop(k,:)); hold on
    title("t = " + num2str(t(k),formatSpec));
    if (k==1), ylabel("z = top of domain"), end
end

for k = times
    nexttile
    plot(x,uinter(k,:)); hold on
    if (k==1), ylabel("z = interface"), end
end
% title(tiled_guy,  "u(x,z,t) for z at the upper boundary and interface")
xlabel(tiled_guy, "x")
ylabel(tiled_guy, "u(x,z,t)")

delete u_inter_boundary.eps
print('u_inter_boundary.eps','-depsc2');

%% plot wave speed on the top and interface over time
tee = [t(2),t(end)];
c1 = par.speeds(1)*[1,1]; c2 = par.speeds(2)*[1,1];

figure; tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
set(gcf,'Position',[450 458 800 420])
nexttile
plot(t(2:end), abs(diff(peak_top)./diff(t')),tee,c1,tee,c2);
ylim([0,max(par.speeds)+0.3])
xlim(tee)
title("Wave speed at top over time")
legend("Wave speed","speed in top domain","speed in bottom domain",'Location','northwest')
nexttile
plot(t(2:end), abs(diff(peak_inter)./diff(t')),tee,c1,tee,c2);
ylim([0,max(par.speeds)+0.3])
xlim(tee)
title("Wave speed at interface over time")

delete wavespeeds.eps
print('wavespeeds.eps','-depsc2');
