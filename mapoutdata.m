%
%% Copywrite (c) 2007 University of Miami
%

%**************************************************************************
%Program: mapoutdata.m
%Project: Kat
%Author: MK Hernandez
%Advisor: Dr. Brian Mapes & Derek Ortt
%Date: 02/07/07
%**************************************************************************

%
%% twenty*[zero/twelve]
%% $lat  $lon $scfp  $scfz  $scftemp  $scftd  $scfdir  $scfwind $850  $850z  $850T  $850td 
%% $850dir $850wind $500  $500T  $500td $500dir $500wind $250  $250T  $250td $250dir $250wind
%

load -ascii twentyfourzero.txt % will load up rawindsonde file
load -ascii coastline.dat % will load up coastline file

contourf (twentyfourzero(:,1), twentyfourzero(:,2), twentyfourzero(:,4))
a = plot (coastline(:,1), coastline (:,2), 'K')
set (a, 'linewidth', 1.2)
title ('Surface pressure 8/24/05/00')
xlabel ('longitude'); ylabel ('latitude');
color bar

hold on

