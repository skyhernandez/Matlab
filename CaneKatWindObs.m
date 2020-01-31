%==================================================================================
%                              Lorenz butterfly equation
%==================================================================================
% This program solves the Lorenz (1963) equation using the forward in time scheme.
%     REFERENCE: http://en.wikipedia.org/wiki/Lorenz_attractor
% From the reference, the Lorenz equations are:
%         dx/dt = s(y-x)
%         dy/dt = x(r-z)-y
%         dz/dt = xy-bz
% These sets of formulas describe the behavior of the Lorenz oscillator, which is a 
% a 3d dynamical system that exhibits chaotic flow.  It is noted for its
% butterfly shape.  These equations are simplified weuations of convection
% rolls arising in our dynamical equations.  
%
% This effect also illustrates that our atmosphre may exhibit a variety
% quasi-periodic regimes that are subject to abrupt changes.
%
%----------------------------------------------------------------------------------
%
% @creator: Fuqing Zhang
% @updated: Michael Kevin Hernandez
% @version: 2.0.0
% @date   : 17 Sept. 2009 
%
%----------------------------------------------------------------------------------
% Version 2.0.0:
%   This version has the Lorenz system solved unter the Leap Frog, Runge
%   Kutta O(2) and O(3) Schemes.  
% 
%==================================================================================

%
%% key parameters 
%
  s = 10.; 
  r = 28.;
  b = 8/3;

%
%% time step
%

  dt = 0.002;
  T = 24;
  TT = T/dt+1-301;
  
%
%% Initializing x, y, & z for two different points to see how they diverge.
%

  x(1:TT)   =0; y(1:TT)   =0; z(1:TT)  =0;
  x1(1:TT)  =0; y1(1:TT)  =0; z1(1:TT) =0;
  x2(1:TT)  =0; y2(1:TT)  =0; z2(1:TT) =0;
  x21(1:TT) =0; y21(1:TT) =0; z21(1:TT)=0;
  x3(1:TT)  =0; y3(1:TT)  =0; z3(1:TT) =0;
  x31(1:TT) =0; y31(1:TT) =0; z31(1:TT)=0;
  x4(1:TT)  =0; y4(1:TT)  =0; z4(1:TT) =0;
  x41(1:TT) =0; y41(1:TT) =0; z41(1:TT)=0;

%
%% Initial conditions
%
%  x(1) = -4.398; y(1) = -3.8409; z(1) =23.2249;

  x(1)   = -10; y(1)   = 40; z(1)   = 40;
  x1(1)  = -10; y1(1)  = 40; z1(1)  = 40.001;
  
  x2(1)  = -10; y2(1)  = 40; z2(1)  = 40;
  x21(1) = -10; y21(1) = 40; z21(1) = 40.001;
  
  x3(1)  = -10; y3(1)  = 40; z3(1)  = 40;
  x31(1) = -10; y31(1) = 40; z31(1) = 40.001;
  
  x4(1)  = -10; y4(1)  = 40; z4(1)  = 40;
  x41(1) = -10; y41(1) = 40; z41(1) = 40.001;

%
%% Begin integration, old values are automatically replaced at the next time step
%

% Euler Scheme 
% Given by the originating code
% [x(t+1) = x(t) + dt*F(x,y,z,t)]

  for t = 1:TT-1
     x(t+1)  = x(t)  + dt * s * ( y(t)-x(t) );
     y(t+1)  = y(t)  + dt * ( r*x(t) - y(t) - x(t)*z(t));
     z(t+1)  = z(t)  + dt * ( x(t)*y(t) - b*z(t) );

     x1(t+1) = x1(t) + dt * s * ( y1(t)-x1(t) );
     y1(t+1) = y1(t) + dt * ( r*x1(t) - y1(t) - x1(t)*z1(t));
     z1(t+1) = z1(t) + dt * ( x1(t)*y1(t) - b*z1(t) );
  end

% Leap_frog (2nd order) 
% Given by <http://public.lanl.gov/balu/nadiga-06b.pdf>
%   and by <www.caam.rice.edu/~caam452/CAAM452Lecture4b.ppt>
% [x(t+1) = x(t-1) + 2*dt*F(x,y,z,t)]
  
  for t = 1:2
     x2(t+1)  = x2(t)  + dt * s * ( y2(t)-x2(t) );
     y2(t+1)  = y2(t)  + dt * ( r*x2(t) - y2(t) - x2(t)*z2(t));
     z2(t+1)  = z2(t)  + dt * ( x2(t)*y2(t) - b*z2(t) );

     x21(t+1) = x21(t) + dt * s * ( y21(t)-x21(t) );
     y21(t+1) = y21(t) + dt * ( r*x21(t) - y21(t) - x21(t)*z21(t));
     z21(t+1) = z21(t) + dt * ( x21(t)*y21(t) - b*z21(t) );
  end

  for t = 2:TT-1
     x2(t+1)  = x2(t-1)  + 2*dt * s * ( y2(t)-x2(t) );
     y2(t+1)  = y2(t-1)  + 2*dt * ( r*x2(t) - y2(t) - x2(t)*z2(t));
     z2(t+1)  = z2(t-1)  + 2*dt * ( x2(t)*y2(t) - b*z2(t) );

     x21(t+1) = x21(t-1) + 2*dt * s * ( y21(t)-x21(t) );
     y21(t+1) = y21(t-1) + 2*dt * ( r*x21(t) - y21(t) - x21(t)*z21(t));
     z21(t+1) = z21(t-1) + 2*dt * ( x21(t)*y21(t) - b*z21(t) );
  end
  
% Runge_Kutta (2nd order)
% Given by <http://www.calvin.edu/~scofield/courses/m231/materials/rungeKuttaFormulas.pdf>
%   and by <http://www.swarthmore.edu/NatSci/echeeve1/Ref/NumericInt/RK2.html>
% dx/dt = f(t,x,y), x(t=0) = xo
% dy/dt = g(t,x,y), y(t=0) = yo
%   x(t+1) = x(t) + dt/3(kx1 + 4kx2)
%   y(t+1) = y(t) + dt/3(ky1 + 4ky2)
%      kx1 = f(t,x,y)
%      ky1 = g(t,x,y)
%      kx2 = f(x(t) + dt*kx1*0.5, y(t) + dt*ky1*0.5)
%      ky2 = g(x(t) + dt*kx1*0.5, y(t) + dt*ky1*0.5)


  for t = 1:TT-1
      
     kx1      = s*(y3(t) - x3(t));     
     ky1      = (r*x3(t) - x3(t)*z3(t) - y3(t));
     kz1      = (x3(t)*y3(t) - b*z3(t));
     
     kx2      = s*((y3(t)+ 0.75*dt*ky1) - (x3(t)+ 0.75*dt*kx1)); 
     ky2      = (r*(x3(t)+ 0.75*dt*kx1) - (y3(t)+ 0.75*dt*ky1) - (x3(t)+ 0.75*dt*kx1)*(z3(t)+ 0.75*dt*kz1));
     kz2      = ((x3(t)+ 0.75*dt*kx1)*(y3(t)+ 0.75*dt*ky1) - b*(z3(t)+ 0.75*dt*kz1));
         
     x3(t+1)  = x3(t) + (2/3*kx2 + 1/3* kx1) * dt;
     y3(t+1)  = y3(t) + (2/3*ky2 + 1/3* ky1) * dt;
     z3(t+1)  = z3(t) + (2/3*kz2 + 1/3* kz1) * dt;   

     kx11     = s*(y3(t) - x3(t));     
     ky11     = (r*x3(t) - y3(t) - x3(t)*z3(t));
     kz11     = (x3(t)*y3(t) - b*z3(t));
     
     kx22     = s*((y31(t)+ 0.75*dt*ky11) - (x31(t)+ 0.75*dt*kx11)); 
     ky22     = (r*(x31(t)+ 0.75*dt*kx11) - (y31(t)+ 0.75*dt*ky11) - (x31(t)+ 0.75*dt*kx11)*(z3(t)+ 0.75*dt*kz1));
     kz22     = ((x31(t)+ 0.75*dt*kx11)*(y31(t)+ 0.75*dt*ky11) - b*(z31(t)+ 0.75*dt*kz11));
         
     x31(t+1) = x31(t) + (2/3*kx22 + 1/3* kx11) * dt;
     y31(t+1) = y31(t) + (2/3*ky22 + 1/3* ky11) * dt;
     z31(t+1) = z31(t) + (2/3*kz22 + 1/3* kz11) * dt;  

  end
    
% Runge_Kutta (3rd order)
% Given by <http://www.calvin.edu/~scofield/courses/m231/materials/rungeKuttaFormulas.pdf>
%   and by <http://people.revoledu.com/kardi/tutorial/ODE/Runge%20Kutta%203.htm>
% dx/dt = f(t,x,y), x(t=0) = xo
% dy/dt = g(t,x,y), y(t=0) = yo
%   x(t+1) = x(t) + dt/6(kx1 + 4kx2 + kx3)
%   y(t+1) = y(t) + dt/6(ky1 + 4ky2 + ky3)
%      kx1 = f(t,x,y)
%      ky1 = g(t,x,y)
%      kx2 = f(x(t) + dt*kx1*0.5, y(t) + dt*ky1*0.5)
%      ky2 = g(x(t) + dt*kx1*0.5, y(t) + dt*ky1*0.5)
%      kx3 = f(x(t) - dt*kx1 + 2*dt*kx2, y(t) - dt*ky1 + 2*dt*ky2)
%      ky3 = g(x(t) - dt*kx1 + 2*dt*kx2, y(t) - dt*ky1 + 2*dt*ky2)

  for t = 1:TT-1
     kx1      = s*(y4(t) - x4(t));     
     ky1      = (r*x4(t) - x4(t)*z4(t) - y4(t));
     kz1      = (x4(t)*y4(t) - b*z4(t));
     
     kx2      = s*((y4(t)+ 0.75*dt*ky1) - (x4(t)+ 0.75*dt*kx1)); 
     ky2      = (r*(x4(t)+ 0.75*dt*kx1) - (y4(t)+ 0.75*dt*ky1) - (x4(t)+ 0.75*dt*kx1)*(z4(t)+ 0.75*dt*kz1));
     kz2      = ((x4(t)+ 0.75*dt*kx1)*(y4(t)+ 0.75*dt*ky1) - b*(z4(t)+ 0.75*dt*kz1));
      
     kx3      = s*((y4(t)- dt*ky1 + 2*ky2*dt) - (x4(t)- dt*kx1 + 2*kx2*dt));
     ky3      = (r*(x4(t)- dt*kx1 + 2*kx2*dt) - (y4(t)- dt*ky1 + 2*ky2*dt) - (x4(t)- dt*kx1 + 2*kx2*dt)*(z4(t)- dt*kz1 + 2*kz2*dt));
     kz3      = ((x4(t)- dt*kx1 + 2*kx2*dt)*(y4(t)- dt*ky1  + 2*ky2*dt) - b*(z4(t)- dt*kz1 + 2*kz2*dt));
     
     x4(t+1)  = x4(t) + (2/3*kx2 + 1/6* kx1 + 1/6*kx3) * dt;
     y4(t+1)  = y4(t) + (2/3*ky2 + 1/6* ky1 + 1/6*ky3) * dt;
     z4(t+1)  = z4(t) + (2/3*kz2 + 1/6* kz1 + 1/6*kz3) * dt;   

     kx11     = s*(y4(t) - x4(t));     
     ky11     = (r*x4(t) - y4(t) - x4(t)*z4(t));
     kz11     = (x4(t)*y4(t) - b*z4(t));
     
     kx22     = s*((y41(t)+ 0.75*dt*ky11) - (x41(t)+ 0.75*dt*kx11)); 
     ky22     = (r*(x41(t)+ 0.75*dt*kx11) - (y41(t)+ 0.75*dt*ky11) - (x41(t)+ 0.75*dt*kx11)*(z4(t)+ 0.75*dt*kz1));
     kz22     = ((x41(t)+ 0.75*dt*kx11)*(y41(t)+ 0.75*dt*ky11) - b*(z41(t)+ 0.75*dt*kz11));
     
     kx33     = s*((y41(t)- dt*ky11 + 2*ky22*dt) - (x41(t)- dt*kx11 + 2*kx22*dt));
     ky33     = (r*(x41(t)- dt*kx11 + 2*kx22*dt) - (y41(t)- dt*ky11 + 2*ky22*dt) - (x41(t)- dt*kx11 + 2*kx22*dt)*(z41(t)- dt*kz11 + 2*kz22*dt));
     kz33     = ((x41(t)- dt*kx11 + 2*kx22*dt)*(y41(t)- dt*ky11  + 2*ky22*dt) - b*(z41(t)- dt*kz11 + 2*kz22*dt));
     
     x41(t+1) = x41(t) + (2/3*kx22 + 1/6* kx11 + 1/6*kx33) * dt;
     y41(t+1) = y41(t) + (2/3*ky22 + 1/6* ky11 + 1/6*ky33) * dt;
     z41(t+1) = z41(t) + (2/3*kz22 + 1/6* kz11 + 1/6*kz33) * dt; 
  end
       
 %
 %% Plotting the calculations
 %
  figure(1);clf 

  subplot(2,2,1)
  for t=1:20:TT
    plot3(x(t),y(t),z(t),'b:');hold on;
    plot3(x1(t),y1(t),z1(t),'r:');
  end
 
  plot3(x(1),y(1),z(1),'b*','Markersize',8); hold on;
  plot3(x1(1),y1(1),z1(1),'ro','Markersize',8); hold on;
  m=2;
  for t = 10:10:TT
      if t < 130 
         m=8;
  plot3(x(t),y(t),z(t),'b*','Markersize',m);
  plot3(x1(t),y1(t),z1(t),'ro','Markersize',m);
      elseif t > TT-10  
         m=12;
  plot3(x(t),y(t),z(t),'b*','Markersize',m);
  plot3(x1(t),y1(t),z1(t),'ro','Markersize',m);
      end
  end
  axis([-20 20 -20 40 0 50]);
  xlabel 'x';ylabel 'y'; zlabel 'z'; 
  title 'lorenz system: (Euler Scheme) '
  grid on;

  subplot(2,2,2)
  for t=1:20:TT
    plot3(x2(t),y2(t),z2(t),'b:');hold on;
    plot3(x21(t),y21(t),z21(t),'r:');
  end
 
  plot3(x2(1),y2(1),z2(1),'b*','Markersize',8); hold on;
  plot3(x21(1),y21(1),z21(1),'ro','Markersize',8); hold on;
  
  for t = 10:10:TT
      if t < 130 
         m=8;
  plot3(x2(t),y2(t),z2(t),'b*','Markersize',m);
  plot3(x2(t),y21(t),z21(t),'ro','Markersize',m);
      elseif t > TT-10  
         m=12;
  plot3(x2(t),y2(t),z2(t),'b*','Markersize',m);
  plot3(x21(t),y21(t),z21(t),'ro','Markersize',m);
      end
  end
  axis([-20 20 -20 40 0 50]);
  xlabel 'x';ylabel 'y'; zlabel 'z';
  title 'lorenz system: (Leap Frog O(2) Scheme)'
  grid on;
  
  subplot(2,2,3)
  for t=1:20:TT
    plot3(x3(t),y3(t),z3(t),'b:');hold on;
    plot3(x31(t),y31(t),z31(t),'r:');
  end
 
  plot3(x3(1),y3(1),z3(1),'b*','Markersize',8); hold on;
  plot3(x31(1),y31(1),z3(1),'ro','Markersize',8); hold on;
 
  for t = 10:10:TT
      if t < 130 
         m=8;
  plot3(x3(t),y3(t),z3(t),'b*','Markersize',m);
  plot3(x31(t),y31(t),z31(t),'ro','Markersize',m);
      elseif t > TT-10  
         m=12;
  plot3(x3(t),y3(t),z3(t),'b*','Markersize',m);
  plot3(x31(t),y31(t),z31(t),'ro','Markersize',m);
      end
  end
  axis([-20 20 -20 40 0 50]);
  xlabel 'x';ylabel 'y'; zlabel 'z';
  title 'lorenz system: (Runge Kutta O(2) Scheme)'
  grid on;

  subplot(2,2,4)
  for t=1:20:TT
    plot3(x4(t),y4(t),z4(t),'b:');hold on;
    plot3(x41(t),y41(t),z41(t),'r:');
  end
 
  plot3(x4(1),y4(1),z4(1),'b*','Markersize',8); hold on;
  plot3(x41(1),y41(1),z41(1),'ro','Markersize',8); hold on;
  
  for t = 10:10:TT
      if t < 130 
         m=8;
  plot3(x4(t),y4(t),z4(t),'b*','Markersize',m);
  plot3(x41(t),y41(t),z41(t),'ro','Markersize',m);
      elseif t > TT-10  
         m=12;
  plot3(x4(t),y4(t),z4(t),'b*','Markersize',m);
  plot3(x41(t),y41(t),z41(t),'ro','Markersize',m);
      end
  end
  axis([-20 20 -20 40 0 50]);
  xlabel 'x';ylabel 'y'; zlabel 'z';
  title 'lorenz system: (Runge Kutta O(3) Scheme)'
  grid on;

  eval ('print -djpeg Lorenz_butterfly.jpg');