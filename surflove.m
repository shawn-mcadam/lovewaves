clearvars;clc; close all
do_gifs = 0;

heavi = @(x,a,b) (a-b)*(x > 0)+b;
lin = @(y,a,b,L) a + (b-a)*y/L;
patchf = @(f,y,a,b,L) (y < 0).*f(0,a,b,L) + (y >= 0).*f(y,a,b,L);
% --------------- c_2 = 2c_1, explosion in bottom -----------------------
par1.rboom = 1/16;
par1.boom_depth = -0.1;
par1.amp = 1;
par1.initf = ...
    @(x,y) par1.amp*exp(-(x.^2 + (y-par1.boom_depth).^2)/par1.rboom^2);
par1.initfp = @(x,y) 0;
par1.speeds = [0.1, 0.2]; % top is face 1, bottom is face 2
par1.speedf = @(y) heavi(y, par1.speeds(1)^2, par1.speeds(2)^2);
par1.tfinal = 10;
par1.depth1 = 0.15; % depth of the interface
par1.depth2 = par1.speeds(2)*par1.tfinal - par1.boom_depth; % size of the second domain
par1.width  = 2*max(par1.speeds)*par1.tfinal+0.5;
par1.reltol = 1.0e-3; par1.abstol = 1.0e-6; % time tolerances
par1.step = 0.0035; % refinement in space
par1.nonlinearity = @(x) 0;

% --------------- c_2 = 2c_1, explosion in top ----------------------
par2.rboom = 1/16;
par2.boom_depth = 0.3;
par2.amp = 1;
par2.initf = ...
    @(x,y) par2.amp*exp(-(x.^2 + (y-par2.boom_depth).^2)/par2.rboom^2);
par2.initfp = @(x,y) 0;
par2.speeds = [0.1, 0.2]; % top is face 1, bottom is face 2
par2.speedf = @(y) heavi(y, par2.speeds(1)^2, par2.speeds(2)^2);
par2.tfinal = 17.5;
par2.depth1 = 0.5; % depth of the interface
par2.depth2 = par2.speeds(2)*par2.tfinal - (par2.speeds(1)/par2.speeds(2))*par2.boom_depth; % size of the second domain
par2.width  = 2*max(par2.speeds)*par2.tfinal;
par2.reltol = 1.0e-3; par2.abstol = 1.0e-6; % time tolerances
par2.step = 0.005; % refinement in space
par2.nonlinearity = @(x) 0;

% --------------- c_1 = 2c_2, explosion in bottom ----------------------
par3.rboom = 1/16;
par3.boom_depth = -0.3;
par3.amp = 1;
par3.initf = ...
    @(x,y) par3.amp*exp(-(x.^2 + (y-par3.boom_depth).^2)/par3.rboom^2);
par3.initfp = @(x,y) 0;
par3.speeds = [0.2, 0.1]; % top is face 1, bottom is face 2
par3.speedf = @(y) heavi(y, par3.speeds(1)^2, par3.speeds(2)^2);
par3.tfinal = 20;
par3.depth1 = 0.5; % depth of the interface
par3.depth2 = par3.speeds(2)*par3.tfinal - par3.boom_depth; % size of the second domain
par3.width  = 2*max(par3.speeds)*par3.tfinal;
par3.reltol = 1.0e-3; par3.abstol = 1.0e-6; % time tolerances
par3.step = 0.005; % refinement in space
par3.nonlinearity = @(x) 0;

% --------------- c_1 = 2c_2, explosion in top ----------------------
par4.rboom = 1/16;
par4.boom_depth = 0.3;
par4.amp = 1;
par4.initf = ...
    @(x,y) par4.amp*exp(-(x.^2 + (y-par4.boom_depth).^2)/par4.rboom^2);
par4.initfp = @(x,y) 0;
par4.speeds = [0.2, 0.1]; % top is face 1, bottom is face 2
par4.speedf = @(y) heavi(y, par4.speeds(1)^2, par4.speeds(2)^2);
par4.tfinal = 20;
par4.depth1 = 0.5; % depth of the interface
par4.depth2 = par4.speeds(2)*par4.tfinal - (par2.speeds(1)/par2.speeds(2))*par4.boom_depth; % size of the second domain
par4.width  = 2*max(par4.speeds)*par4.tfinal+0.5;
par4.reltol = 1.0e-3; par4.abstol = 1.0e-6; % time tolerances
par4.step = 0.005; % refinement in space
par4.nonlinearity = @(x) 0;

disp("Using the following parameters")
par = par2
tic
[t,x,y,u] = lovewave( ...
    par.speedf, ...
    par.tfinal, ...
    par.depth1, par.depth2, par.width, ...
    par.initf, par.initfp, ...
    par.reltol, par.abstol, ...
    par.step, ...
    par.nonlinearity);
toc
disp("Solved")

%% init important variables
zmax = par.amp;
Ninter = floor(par.depth1/par.step);
mindelay = 0.015;
utop   = u(:,1,:);
uinter = u(:,Ninter,:);
min_prom = 0.001;


% Track the maximum of a spline of u on the interface and top of the domain
finterp = @(x,y,xs) spline(x,y,xs);
xref = linspace(x(1),x(end),100*length(x)); % refined xinputs
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
        zlim([-par.amp/3,par.amp])
        clim([-par.amp/18,par.amp/14])
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
        ylim([-par.amp/6,par.amp/6])
        title("u(x,y,t) at the top of the domain");

        nexttile
        plot(x,uinter(k,:)); hold on;
        plot(peak_inter(k),finterp(x,uinter(k,:),peak_inter(k)),"x");
        title("u(x,y,t) at the interface of the domain");
        ylim([-par.amp/10,par.amp/6])
        gif
    end
end
end

%% Heatmaps
figure
formatSpec = '%.2f';
set(gcf,'Position',[450 458 650 450])
Nfigures = 6;
times = floor(linspace(1,length(t),Nfigures));
tiled_guy = tiledlayout(2,Nfigures/2,'TileSpacing','Compact','Padding','Compact');
for k = times
    nexttile
    surf(x,y(1:Ninter),squeeze(u(k,1:Ninter,:)),EdgeColor='none');
    hold on
    surf(x,y(Ninter:end),squeeze(u(k,Ninter:end,:)),EdgeColor='none',FaceAlpha=0.8);
    hold off
    view(0,90)
    xlim([x(1),x(end)])
    ylim([y(end),y(1)])
    colorbar
    title("t = " + num2str(t(k),formatSpec));
end
title(tiled_guy,"Vertical displacement u(x,z,t) over time")
xlabel(tiled_guy,"x")
ylabel(tiled_guy,"z")


delete overhead.jpg
print('overhead.eps','-depsc2','-r500');


%% 1D waves
%figure;
%formatSpec = '%.2f';
%Nfigures = 8;
%times = floor(linspace(1,length(t),Nfigures+1)); times(1)=[];
%set(gcf,'Position',[450 458 800 420])
%tiled_guy = tiledlayout(4,Nfigures/2,'TileSpacing','Compact','Padding','Compact');
%for k = times(1:Nfigures/2)
    %nexttile
    %plot(x,utop(k,:)); hold on
    %title("t = " + num2str(t(k),formatSpec));
    %ylim([-par.amp/10,par.amp/8])
    %if (k==times(1)), ylabel("z = top of domain"), end
%end
%
%for k = times(1:Nfigures/2)
    %nexttile
    %plot(x,uinter(k,:)); hold on
    %ylim([-par.amp/14,par.amp/10])
    %if (k==times(1)), ylabel("z = interface"), end
%end
%for k = times(Nfigures/2+1:end)
    %nexttile
    %plot(x,utop(k,:)); hold on
    %title("t = " + num2str(t(k),formatSpec));
    %ylim([-par.amp/10,par.amp/8])
    %if (k==times(Nfigures/2+1)), ylabel("z = top of domain"), end
%end
%
%for k = times(Nfigures/2+1:end)
    %nexttile
    %plot(x,uinter(k,:)); hold on
    %ylim([-par.amp/14,par.amp/10])
    %if (k==times(Nfigures/2+1)), ylabel("z = interface"), end
%end
% title(tiled_guy,  "u(x,z,t) for z at the upper boundary and interface")
%xlabel(tiled_guy, "x")
%ylabel(tiled_guy, "u(x,z,t)")
%
%delete u_inter_boundary.eps
%print('u_inter_boundary.eps','-depsc2');

%% plot wave speed on the top and interface over time
tee = [t(2),t(end)];
c1 = par.speeds(1)*[1,1]; c2 = par.speeds(2)*[1,1];

figure; tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(t(2:end), abs(diff(peak_top)./diff(t')),tee,c1,tee,c2);
ylim([0,max(par.speeds)+0.3])
title("Wave speed at top over time")
legend("Wave speed","speed in top domain","speed in bottom domain",'Location','northwest')
nexttile
plot(t(2:end), abs(diff(peak_inter)./diff(t')),tee,c1,tee,c2);
ylim([0,max(par.speeds)+0.3])
title("Wave speed at interface over time")

delete wavespeeds.eps
print('wavespeeds.eps','-depsc2');

