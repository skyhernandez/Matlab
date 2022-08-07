%
%% Copywrite (c) 2007 University of Miami
%

%**********************************************************************
%Program: CaneKatObs.m
%Project: Kat
%Author: Sky Hernandez
%Advisor: Dr. Brian Mapes & Derek Ortt
%Date: 02/07/07
%**********************************************************************

%
%% twenty*[zero/twelve]
%% 1$lat  2$lon 3$scfp  5$scftemp  6$scftd  7$scfdir  8$scfwind 10$850z  %% 11$850T 12$850td 13$850dir 14$850wind 16$500z  17$500T  18$500td 
%% 19$500dir 20$500wind 22$250z 23$250T 24$250td 25$250dir 26$250wind
%

load -ascii twentyfivezero2.txt% will load up rawindsonde/dropsonde file
load -ascii twentysixzero2.txt% will load up rawindsonde/dropsonde file
load -ascii twentysevenzero2.txt% will load up rawindsonde/dropsonde file
load -ascii coastline.dat % will load up coastline file

XX = 17;

y=twentyfivezero2(:,1); x=twentyfivezero2(:,2); z=twentyfivezero2(:,XX);
[xi,yi]=meshgrid(-97.5:.25:-76.5, 17.5:.25:35);
z1=griddata(x,y,z,xi,yi,'linear');
subplot(3,1,1);pcolor(xi,yi,z1);
shading interp;
hold on;

a = plot (coastline(:,1), coastline (:,2), 'k')
set (a, 'linewidth', 1.2)
title ('Title 8/25/05/00')
ylabel ('latitude');
colorbar

y=twentysixzero2(:,1); x=twentysixzero2(:,2); z=twentysixzero2(:,XX);
[xi,yi]=meshgrid(-97.5:.25:-76.5, 17.5:.25:35);
z1=griddata(x,y,z,xi,yi,'linear');
subplot(3,1,2);pcolor(xi,yi,z1);
shading interp;
hold on;

eval (['print -djpeg title825200500.jpg']);
 
