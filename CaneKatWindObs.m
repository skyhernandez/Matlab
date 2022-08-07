%
%% Copywrite (c) 2007 University of Miami
%

%**************************************************************************
%Program: CaneKatWindObs.m
%Project: Kat
%Author: Sky Hernandez
%Advisor: Dr. Brian Mapes & Derek Ortt
%Date: 02/07/07
%**************************************************************************

%
%% twenty*[zero/twelve]
%% 1$lat  2$lon 7$scfdir 8$scfwind 13$850dir 14$850wind 19$500dir 20$500wind
%% 25$250dir 26$250wind
%

load -ascii twentyfivezero2.txt % will load up rawindsonde/dropsonde file
load -ascii twentysixzero2.txt % will load up rawindsonde/dropsonde file
load -ascii twentysevenzero2.txt % will load up rawindsonde/dropsonde file
load -ascii coastline.dat % will load up coastline file

XX = 8
dir = 7
day = twentysevenzero2

theta = day(:,dir); vel= 1;

[x1,y1]=meshgrid(-97.5:.25:-76.5, 17.5:.25:35);
y1=day(:,1); x1=day(:,2);
v=cos(theta).*vel; u=sin(theta).*vel;

figure;
quiver(x1,y1,-u,-v,0.35, 'r');
hold on;
b=day(:,1); aa=day(:,2); c=day(:,XX);
[aai,bi]=meshgrid(-97.5:.25:-76.5, 17.5:.25:35);
c1=griddata(aa,b,c,aai,bi,'linear');
pcolor(aai,bi,c1);
shading interp;
colormap hsv
hold on;
a = plot (coastline(:,1), coastline (:,2), 'k')
set (a, 'linewidth', 1.2)
title ('scf winds 8/27/05/00')
ylabel ('latitude');
colorbar



eval (['print -djpeg scfnds8270500.jpg']);
