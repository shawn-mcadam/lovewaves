clearvars;clc; close all

heavi = @(x,a,b) (a-b)*(x > 0)+b;
lin = @(y,a,b,L) a + (b-a)*y/L;
patchf = @(f,y,a,b,L) (y < 0).*f(0,a,b,L) + (y >= 0).*f(y,a,b,L);
% --------------- first set of example parameters -----------------------
par1.rboom = 1/20;
par1.boom_depth = -0.7;
par1.hboom = 1;
par1.initf = ...
    @(x,y) par1.hboom*exp(-(x.^2 + (y-par1.boom_depth).^2)/par1.rboom);
par1.initfp = @(x,y) 0;
par1.speeds = [0.1, 0.2]; % top is face 1, bottom is face 2
par1.speedf = @(y) heavi(y, par1.speeds(1)^2, par1.speeds(2)^2);
par1.tfinal = 35;
par1.depth1 = 0.5; % depth of the interface
par1.depth2 = par1.speeds(2)*par1.tfinal - par1.boom_depth; % size of the second domain
par1.width  = 2*max(par1.speeds)*par1.tfinal;
par1.reltol = 1.0e-3; par1.abstol = 1.0e-6; % time tolerances
par1.step = 0.04; % refinement in space
par1.nonlinearity = @(x) 0;

% --------------- second set of example parameters ----------------------
lin_zero = @(x,a,b) (x < a).*(x - a) + (b < x).*(x - b);
par2.rboom = 1/10;
par2.boom_depth = 0.3;
par2.hboom = 1;
par2.initf = @(x,y) ...
    par2.hboom*exp(-(x.^2 + lin_zero(y,-par2.boom_depth,par2.boom_depth).^2)/par2.rboom);
par2.initfp = @(x,y) 0;
par2.speeds = [0.1,0.4];
par2.speedf = @(y) heavi(y, par2.speeds(2)^2, par2.speeds(1)^2);
par2.tfinal = 35;
par2.L = max(par2.speeds)*par2.tfinal;
par2.depth1 = 2;
par2.depth2 = 18;
par2.width  = 2*par2.L;
par2.reltol = 1.0e-3; par2.abstol = 1.0e-6;
par2.step = 0.1;
par2.nonlinearity = @(x) 0;

% --------------- symmetric set of example parameters ----------------------
par3.rboom = 1/10;
par3.boom_depth = 0.2;
par3.hboom = 1;
par3.initf = @(x,y) ...
    par3.hboom*exp(-(x.^2 + lin_zero(y,-par3.boom_depth,par3.boom_depth).^2)/par3.rboom);
par3.initfp = @(x,y) 0;
par3.speeds = [0.1,0.4];
par3.speedf = @(y) heavi(y, par3.speeds(2)^2, par3.speeds(1)^2);
par3.tfinal = 13;
par3.L = max(par3.speeds)*par3.tfinal;
par3.depth1 = 8; % depth of the interface
par3.depth2 = 8;
par3.width  = 2*par3.L;
par3.reltol = 1.0e-3; par3.abstol = 1.0e-6;
par3.step = 0.075;
par3.nonlinearity = @(x) 0;

% -------------- wave speed as a function ------------------------
par4 = par1;
par4.speeds = [0.1, 0.2];
par4.boom_depth = -0.5;
par4.initf = @(x,y) par4.hboom*exp(-(x.^2 + (y-par4.boom_depth).^2)/par4.rboom);
par4.initfp = @(x,y) 0;
par4.tfinal = 15;
par4.L = par4.speeds(2)*par4.tfinal;
par4.width  = 2*par4.L;
par4.depth1 = 0.75;
par4.depth2 = par4.L;
par4.speedf = @(y) patchf(lin,y,par4.speeds(2)^2,par4.speeds(1)^2,par4.depth1);
par4.step = 0.08;
par4.nonlinearity = @(x) 0;


% ----------- simulate love waves. Need a cylindrical domain... -----------------
par5.speeds = [0.5, 1]; % top is face 1, bottom is face 2
par5.k = 1;
par5.v = 0.6343;
par5.hboom = 1;
par5.depth1 = 1; % depth of the interface
par5.omega1 = par5.k*sqrt((par5.v/par5.speeds(1))^2-1);
par5.omega2 = par5.k*sqrt(1-(par5.v/par5.speeds(2))^2);
par5.depthf = @(z) (z > 0).*cos(par5.omega1*(z-par5.depth1)) + ...
    (z <= 0).*cos(par5.omega1*par5.depth1).*exp(par5.omega2*z);
par5.initf  = @(x,y) cos(par5.k*x).*par5.depthf(y);
par5.initfp = @(x,y) -par5.k*par5.v*sin(par5.k*x).*par5.depthf(y);
par5.speedf = @(y) heavi(y, par5.speeds(1)^2, par5.speeds(2)^2);
par5.tfinal = 7;
par5.depth2 = 4; % size of the second domain
par5.width  = 5*pi;
par5.reltol = 1.0e-3; par5.abstol = 1.0e-6; % time tolerances
par5.step = 0.04; % refinement in space
par5.nonlinearity = @(x) 0;


par = par1;
[t,x,y,u] = lovewave( ...
    par.speedf, ...
    par.tfinal, ...
    par.depth1, par.depth2, par.width, ...
    par.initf, par.initfp, ...
    par.reltol, par.abstol, ...
    par.step, ...
    par.nonlinearity);



%% init important variables
zmax = par.hboom;
Ninter = floor(par.depth1/par.step);
mindelay = 0.015;
utop   = u(:,1,:);
uinter = u(:,Ninter,:);

%% surf plot
figure
set(gcf,'Position',[100, 100, 1500, 1000])
delete surf.gif
gif('surf.gif','nodither','DelayTime',mindelay);
delta_t = 0.1; tlast = -100;
nrm = 0;
for k = 1:length(t)
    if t(k) - tlast > delta_t
        tlast = t(k);
        % surf(x,y(1:Ninter),squeeze(u(k,1:Ninter,:)),EdgeColor='none');
        % hold on
        % surf(x,y(Ninter:end-nrm),squeeze(u(k,Ninter:end-nrm,:)),EdgeColor='none',FaceAlpha=0.8);
        % hold off
        pcolor(x,y,squeeze(u(k,:,:)));
        colorbar
        %c = max(abs(min(u(k,:,:),[],"all")), abs(max(u(k,:,:),[],"all")));
        %clim([-c,c]);
        %zlim([-zmax/2,zmax]);
        gif
    end
end

%% plots of top of domain and interface
figure
delete interface_vs_top.gif
gif('interface_vs_top.gif','nodither','DelayTime',mindelay);
delta_t = 0.1; tlast = -100;
for k = 1:length(t)
    if t(k) - tlast > delta_t
        tlast = t(k);
        subplot(2,1,1);
        plot(x,utop(k,:));
        ylim([-zmax/2,zmax]);
        title("u(x,y,t) at the top of the domain");
        subplot(2,1,2);
        plot(x,uinter(k,:));
        ylim([-zmax/2,zmax]);
        title("u(x,y,t) at the interface of the domain");
        gif
    end
end

%%
% Track the maximum of a piecewise spline of the data
finterp = @(x,y,xs) spline(x,y,xs);
xref = linspace(x(1),x(end),100*length(x)); % refined xinputs

peak_top   = zeros(1,length(t));
peak_inter = zeros(1,length(t));
for ii = 1:length(t)
    TF = islocalmax(finterp(x,utop(ii,:),xref),'MinProminence',0.001);
    if sum(TF) ~= 0
        indices = 1:length(TF); I = indices(TF);
        peak_top(ii) = xref(I(end));
    end
    TF = islocalmax(finterp(x,uinter(ii,:),xref),'MinProminence',0.001);
    if sum(TF) ~= 0
        indices = 1:length(TF); I = indices(TF);
        peak_inter(ii) = xref(I(end));
    end
end

nrm = 0; % trimming last points because that tracks the wrong peak
tee = [t(2),t(end)];
c = par.speeds;

figure
subplot(2,1,1);
plot(t(2:end), diff(peak_top)./diff(t'),...
    tee, c(1)*[1,1], tee, c(2)*[1,1]);
ylim([0,max(c)+0.3])
title("Wave speed at top over time")
legend("Wave speed","speed in top domain","speed in bottom domain")
subplot(2,1,2);
plot(t(2:end-nrm), diff(peak_inter(1:end-nrm))./diff(t(1:end-nrm))',...
    tee, c(1)*[1,1], tee, c(2)*[1,1]);
ylim([0,max(c)+0.3])
title("Wave speed at interface over time")

% delete wavespeed_time.png
% exportgraphics(gcf,'wavespeed_time.png', 'Resolution', 1000)
