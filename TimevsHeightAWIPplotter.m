%
%##################################################################
%##################################################################
%######                                                      ######
%######                    AWIPS Displayer                   ######
%######                                                      ######
%######                     Developed by                     ######
%######     Center for Analysis and Prediction of Storms     ######
%######                University of Oklahoma                ######
%######                                                      ######
%##################################################################
%##################################################################
%
%Program: AWIPSDisplayer.m
%Project: AWIPS
%Author: Sky Hernandez
%Date: 05/18/2012
%**********************************************************************

%%% 1$Time (mmddyyyyhhmm), 2$Tower observation height, 
%%% 3$Observed wind speed, 4$Differenced wind speed obs-fct 

%%% Loading files
%filename = sample.txt;
%fname = dir('*.txt');
AWIPSfilename = load ('diff_01312012_nttu.txt', '-ascii');
%date=sscanf(fname,'diff_%li')';
%disp(date)
nttulevs = 10;
%disp('Contents of AWIPSfilename:');
%whos -file AWIPSfilename

%%% Read in time for the x-dimension and height for the y-dimension
time=AWIPSfilename(:,1); height=AWIPSfilename(:,2); 
%%% Read in observed and differenced wind speed for countouring
observedspd   = AWIPSfilename(:,3);
differencespd = AWIPSfilename(:,4);
observeddir   = AWIPSfilename(:,5);
differencedir = AWIPSfilename(:,6);
forecastedspd    = differencespd+observedspd;
forecasteddir = differencedir+observeddir;
timesheight = length(observedspd)/nttulevs;
numpoints   = length(observedspd);



    for t = 1:timesheight
        for k = 1:nttulevs
            mat_observedspd(t,k)   = observedspd(nttulevs*(t-1)+k);
            mat_differencespd(t,k) = differencespd(nttulevs*(t-1)+k);
            mat_forecastedspd(t,k) = forecastedspd(nttulevs*(t-1)+k);
            mat_observeddir(t,k)   = observeddir(nttulevs*(t-1)+k);
            mat_differencedir(t,k) = differencedir(nttulevs*(t-1)+k);
            mat_forecasteddir(t,k) = forecasteddir(nttulevs*(t-1)+k);         
        end
    end
    for t = 1:timesheight
        vec_time(t) = time((t*nttulevs)-(nttulevs-1))-013120120000;
    end
    for k = 1:nttulevs
        vec_height(k)=height(k);
    end
   
  
%[xi,yi]=meshgrid(-97.5:.25:-76.5, 17.5:.25:35);
%z1=griddata(x,y,z,xi,yi,'linear');

figure;
subplot(3,1,1);pcolor(vec_time,vec_height,transpose(mat_observedspd)); 
hold on; contour(vec_time,vec_height,transpose(mat_observeddir));
shading interp;
hold on;
title ('observed');
xlabel ('Time');
ylabel ('Height');
colorbar;
caxis([0 25])

subplot(3,1,2);pcolor(vec_time,vec_height,transpose(mat_forecastedspd));
hold on; contour(vec_time,vec_height,transpose(mat_forecasteddir));
shading interp;
hold on;
title ('forecasted');
xlabel ('Time');
ylabel ('Height');
colorbar;
caxis([0 25])

subplot(3,1,3);pcolor(vec_time,vec_height,transpose(mat_differencespd));
hold on; contour(vec_time,vec_height,transpose(mat_differencedir));
shading interp;
hold on;
title ('difference');
xlabel ('Time');
ylabel ('Height');
colorbar;
caxis([-25 25])
eval ('print -djpeg ObsFctDiffplot.jpg');


figure;
x = -25:1:25;
hist(differencespd,x);
xlabel ('differencespd');
ylabel ('frequency');
eval ('print -djpeg histogramspd.jpg');

figure;
x = -360:1:0;
hist(differencedir,x);
xlabel ('differencedir');
ylabel ('frequency');
eval ('print -djpeg histogramdir.jpg');
%a = plot (coastline(:,1), coastline (:,2), 'k')
%set (a, 'linewidth', 1.2)
%title ('Title 8/25/05/00');
%xlabel ('Time');
%ylabel ('Height');
%colorbar

%y=twentysixzero2(:,1); x=twentysixzero2(:,2); z=twentysixzero2(:,XX);
%[xi,yi]=meshgrid(-97.5:.25:-76.5, 17.5:.25:35);
%z1=griddata(x,y,z,xi,yi,'linear');
%subplot(3,1,2);pcolor(xi,yi,z1);
%shading interp;
%hold on;


