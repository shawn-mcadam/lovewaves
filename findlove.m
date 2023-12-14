clearvars; close all; clc
% parameters for speed = 0.05
par.reltol = 1.0e-6; par.abstol = 1.0e-6; % time tolerances

par.speed = 0.05;
par.depth1 = 1; % depth of the interface
par.rboom = 1/20; par.boom_depth = -0.7; par.hboom = 1;
par.tfinal = 1; par.step = 0.06; par.nonlin = 0;
par.initf = ...
    @(x,y) par.hboom*exp(-(x.^2 + (y-par.boom_depth).^2)/par.rboom);
p = [par, par, par];
p(1).step = 0.04; p(1).nonlin = 0; p(1).tfinal = 46; 
p(2).step = 0.04; p(2).nonlin = 1; p(2).tfinal = 44.5; 
p(3).step = 0.015; p(3).nonlin = 2; p(3).tfinal = 30; 

% parameters for speed = 0.2
par.speed = 0.2;
par.depth1 = 0.5;
par.boom_depth = -0.2;
par.initf = ...
    @(x,y) par.hboom*exp(-(x.^2 + (y-par.boom_depth).^2)/par.rboom);
par.step = 0.02;

p = [p par par par];
p(4).tfinal = 25; p(4).nonlin = 0;
p(5).tfinal = 20; p(5).nonlin = 0.1;
p(6).tfinal = 15; p(6).nonlin = 0.2;

heavi_step = @(x,a,b) (a-b)*(x > 0)+b;
for i=1:length(p)
    p(i).speeds = [0.1, p(i).speed];
    p(i).speedf = @(y) heavi_step(y, p(i).speeds(1)^2, p(i).speeds(2)^2);
    p(i).depth2 = p(i).speeds(2)*p(i).tfinal - p(i).boom_depth; % size of the second domain
    p(i).width  = 2*max(p(i).speeds)*p(i).tfinal;
    p(i).nonlinearity = @(y) heavi_step(y, p(i).nonlin, 0);
    
    tic
    [t,x,y,u] = lovewave( ...
        p(i).speedf, ...
        p(i).tfinal, ...
        p(i).depth1, p(i).depth2, p(i).width, ...
        p(i).initf, ...
        p(i).reltol, p(i).abstol, ...
        p(i).step, ...
        p(i).nonlinearity);
    toc

    utop   = u(:,1,:);
    uinter = u(:,floor(p(i).depth1/p(i).step),:);
    % Track the maximum of a piecewise spline of the data
    method = "spline"; minp = 0.001;
    xq = linspace(x(1),x(end),200*length(x)); % refined xinputs
    peak_top   = zeros(1,length(t)); peak_inter = zeros(1,length(t));
    for ii = 1:length(t)
        TF = islocalmax(interp1(x,utop(ii,:),xq,method),'MinProminence',minp);
        if sum(TF) ~= 0
            indices = 1:length(TF); I = indices(TF);
            peak_top(ii) = xq(I(end));
        end
        TF = islocalmax(interp1(x,uinter(ii,:),xq,method),'MinProminence',minp);
        if sum(TF) ~= 0
            indices = 1:length(TF); I = indices(TF);
            peak_inter(ii) = xq(I(end));
        end
    end
    
    nrm = 0; % trimming last points because that tracks the wrong peak
    tee = [t(2),t(end)];
    c = p(i).speeds;
    
    figure
    tfig = tiledlayout(2,1,'TileSpacing','Compact');
    title(tfig, "Ratio of speeds = " + string(c(1)/c(2)) + ". Nonlinearity factor = " + string(p(i).nonlin) + ".")
    xlabel(tfig, 'Time  (s)')
    ylabel(tfig, 'Speed (m/s)')
    
    nexttile
    plot(t(2:end), diff(peak_top)./diff(t'),...
        tee, c(1)*[1,1], tee, c(2)*[1,1]);
    ylim([0,max(c)+0.3])
    title("Outermost wave at the top of the domain")
    legend("Wave speed","speed in top domain","speed in bottom domain")
    
    nexttile
    plot(t(2:end-nrm), diff(peak_inter(1:end-nrm))./diff(t(1:end-nrm))',...
        tee, c(1)*[1,1], tee, c(2)*[1,1]);
    ylim([0,max(c)+0.3])
    title("Outermost wave at the interface of the domain")

    namej = "findlove/s" + string(p(i).speed) + "e" + string(p(i).nonlin);
    namej = strrep(namej,".","_")  + ".png";
    exportgraphics(gcf,namej, 'Resolution', 1000)
end
