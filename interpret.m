par1.name = "bottom_c1c2";
par1.depth1 = 0.15;
par1.step = 0.007;

par2.name = "top_c1c2";
par2.depth1 = 0.3;
par2.step = 0.0045;

par3.name = "bottom_c2c1";
par3.depth1 = 0.15;
par3.step = 0.0045;





par = par3;

Ninter = floor(par.depth1/par.step);
t = readmatrix(par.name + "/t.txt");
x = readmatrix(par.name + "/x.txt");
y = readmatrix(par.name + "/y.txt");

figure
Nfigures = 8;
times = floor(linspace(1,length(t),Nfigures));

formatSpec = '%.2f';
set(gcf,'Position',[450 458 1000 400]);
tiled_guy = tiledlayout(2,Nfigures/2,'TileSpacing','tight','Padding','tight');
for k = times
    u = load(par.name + "/u_t"+string(k)+".mat").bruh;
    nexttile
    surf(x,y(1:Ninter),u(1:Ninter,:),EdgeColor='none');
    hold on
    surf(x,y(Ninter:end),u(Ninter:end,:),EdgeColor='none',FaceAlpha=0.8);
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
print(par.name + '/overhead.jpg','-djpeg','-r1500');

