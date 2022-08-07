%
%% Copywrite (c) 2007 University of Miami
%
%**********************************************************************
%Program: Katdifferencefields.m
%Project: Kat
%Author: Sky Hernandez
%Advisor: Dr. Brian Mapes & Derek Ortt
%Date: 02/07/07
%**********************************************************************
%
%% This will graph the coastline file first, and overlay the difference
%% of the (CMC-MM5).  NOTE: 28=sfc; 19=850; 16=700; 12=500; 6=250
%

    clear all; close all;

    load -ascii coastline.dat % will load up coastline file
    a=load('c2005082500.mat'); b=load('g2005082500.mat');
    c=load('c2005082506.mat'); d=load('g2005082506.mat');
    e=load('c2005082512.mat'); f=load('g2005082512.mat');
    g=load('c2005082518.mat'); h=load('g2005082518.mat');
        
        
% vars for CMC run
    y=(a.LATC(1:199,1:299)); x=(a.LONC(1:199,1:299)); 
    p1=(a.TPRESSURE(1:199,1:299,6)); p2=(c.TPRESSURE(1:199,1:299,6)); 
    p3=(e.TPRESSURE(1:199,1:299,6)); p4=(g.TPRESSURE(1:199,1:299,6)); 
  
   
% vars for MM5 run
    p8=(b.TPRESSURE(1:199,1:299,6)); p9=(d.TPRESSURE(1:199,1:299,6)); 
    p10=(f.TPRESSURE(1:199,1:299,6)); p11=(h.TPRESSURE(1:199,1:299,6)); 
    
       
% differences between the vars
    pa=p1-p8; pb=p2-p9; pc=p3-p10; pd=p4-p11; 

% % graphing of V vars
    subplot (2,2,1);
    PCOLOR(x,y,pa);
    TITLE('8/25/05/00Z');
    shading interp; colorbar;
    AXIS([-100, -70, 15, 35]);
    hold on;
    z = plot (coastline(:,1), coastline(:,2), 'k'); 
    set (z, 'linewidth', 1.2);
   
    subplot (2,2,2);
    PCOLOR(x,y,pb);
    TITLE('8/25/05/06Z');
    shading interp; colorbar;
    AXIS([-100, -70, 15, 35]);
    hold on;
    z = plot (coastline(:,1), coastline(:,2), 'k'); 
    set (z, 'linewidth', 1.2);

    subplot (2,2,3);    
    PCOLOR(x,y,pc);
    TITLE('8/25/05/12Z');
    shading interp; colorbar;
    AXIS([-100, -70, 15, 35]);
    hold on;
    z = plot (coastline(:,1), coastline(:,2), 'k'); 
    set (z, 'linewidth', 1.2);

    subplot (2,2,4);
    PCOLOR(x,y,pd);
    TITLE('8/25/05/18Z');
    shading interp; colorbar;
    AXIS([-100, -70, 15, 35]);
    hold on;
    z = plot (coastline(:,1), coastline(:,2), 'k'); 
    set (z, 'linewidth', 1.2);
   
% Print to file    
    eval (['print -djpeg diffTPRESSURE25082505.jpg']);
