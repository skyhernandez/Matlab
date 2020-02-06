%
%% Copywrite (c) 2007 University Coporation for Atmospheric Research
%

%**********************************************************************
%Program:  mappingIndeces.m
%Project:  COSMIC CAPER
%Author:   MK Hernandez
%Debugger: Douglas Hunt
%Date:     6 July 2007
%**********************************************************************

%
%% Indeces text files
%% 1$indices 2$lat  3$lon 4$mycape  5$cape  6$mycin  7$cin  

load -ascii <indexfilename> % will load up rawindsonde and cosmic file
load -ascii coastline.dat % will load up global coastline file

XX = 4;

y = <indexfilename> (:,2); 
x = <indexfilename> (:,3); 
z = <indexfilename> (:,XX);
[xi,yi]=meshgrid(-180:.25:180, -90:.25:90);
z1 = griddata(x,y,z,xi,yi,'linear');
subplot(2,1,1); pcolor(xi,yi,z1);
shading interp;
hold on;
a = plot (coastline(:,1), coastline (:,2), 'k');
set (a, 'linewidth', 1.2);
title ('CAPE for 8/25/05/00');
ylabel ('latitude');
xlabel (‘longitude’);
colorbar

YY = 6;

y = <indexfilename> (:,2); 
x = <indexfilename> (:,3); 
z = <indexfilename> (:,YY);
[xi,yi]=meshgrid(-180:.25:180, -90:.25:90);
z2 = griddata(x,y,z,xi,yi,'linear');
subplot(2,1,2); pcolor(xi,yi,z2);
shading interp;
hold on;
a = plot (coastline(:,1), coastline (:,2), 'k')
set (a, 'linewidth', 1.2)
title (‘CIN for 8/25/05/00')
ylabel ('latitude');
xlabel (‘longitude’);
colorbar

eval (['print -djpeg title825200500.jpg']);

