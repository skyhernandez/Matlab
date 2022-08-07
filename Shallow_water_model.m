%==================================================================================
%                              Advection equation
%==================================================================================
% This program solves the advection with forcing numerically. From the 
% reference, the advection equation is:
%         du'/dt + U(du'/dx) +g(dh'/dt) = 0
%         dh'/dt + U(dh'/dx) +H(du'/dt) = 0
%
%----------------------------------------------------------------------------------
%
% @creator: Fuqing Zhang
% @updated: Skylar Hernandez
% @version: 2.0.0
% @date   : 2 Nov. 2009 
%
%----------------------------------------------------------------------------------
% Version 2.0.0:
%   This version has the 2d shallow water advection formula solved under
%   the Runge Kutta O(3) Schemes, using different initial perturbations.  This 
%   program includes the linear and nonlinear model.
% 
%==================================================================================

  clear all;         % clean up graphic and memory

%
%% Defintion of key parameters
% 

  
  L = 100000;         % The length of the domain in meters
  T = 3600  ;         % Total integration time in seconds
  g = 9.81  ;         % Constant acceleration towards the center of the earth
  u = 10    ;         % Total mean wind speed in m/s
  h = 10    ;         % Mean hieght in m

  
  dx = 2500 ;         % Grid spacing in meters
  dt = 1    ;         % Time step in seconds

  II = L/dx +1;       % Total grid points in x direction
  TT = T/dt +1;       % Total time steps of integration

%
%% Initialize the 2-D matrix of U and H with 0s
%

  U(TT, II)= 0.;
  H(TT, II)= 0.;
  Ustar(TT,II)  = 0.;
  U2star(TT,II) = 0.;
  Hstar(TT,II)  = 0.;
  H2star(TT,II) = 0.;

   
%
%% Define inititial conditions of U
%
 disp('*************************************************************************************')
 disp('*  Welcome to the 2-dimentional Advection forecast model using O(3) Runge Kutta     *')
 disp('*  Created by Michael Kevin Hernandez (C) 2009                                      *')
 disp('*************************************************************************************')


 IC = input('Type (1) sinusoidal or (2) Gaussian perturbation of hieght fields = ');                                                             
 
  for i = 1:1:II 
    if (IC == 1)
    %% sinusoidal perurbation of Hieght field
       H(1, i)      =  cos( 2*(i-1)*dx*pi/L );  
       Hstar(1, i)  =  cos( 2*(i-1)*dx*pi/L );
       H2star(1, i) =  cos( 2*(i-1)*dx*pi/L );
    elseif (IC == 2)    
    %% Gaussian perturbation of Hieght field with mean of 20 and SD of 3 and amplitude of 1.5.
       H(1, i)      = 1.5*exp(-1*(i-20)^2/(2*(3)^2));
       Hstar(1, i)  = 1.5*exp(-1*(i-20)^2/(2*(3)^2));
       H2star(1, i) = 1.5*exp(-1*(i-20)^2/(2*(3)^2));
    end
  end


 modeleqns = input('Type (1) Non Linearized or (2) Linearized model equations = ');                                                             

%
%% Discretize the PDE using Runge Kutta O(3) scheme.
%% Numerical integration or prediction.
%

if (modeleqns == 1)
  for t = 1:TT-1                                                       % start time integration
      for i = 2 : II-1                                                 % loop over grid points in space
                                                                     % solving with the discretized equation
                                                                     % solving with the discretized equation
%        Hprime(t+1,i)= Hprime(t,i)-(dt/(dx*2))*(U*(Hprime(t,i+1)-Hprime(t,i-1))+H*(Uprime(t,i+1)-Uprime(t,i-1))); 
%        Uprime(t+1,i)= Uprime(t,i)-(dt/(dx*2))*(U*(Uprime(t,i+1)-Uprime(t,i-1))+g*(Hprime(t,i+1)-Hprime(t,i-1)));                                                                          

      
      % Runge-Kutta O(3) from <http://cires.colorado.edu/science/groups/pielke/classes/at7500/wrfarw.pdf>
      % U*(t+1,i) = U'(t,i) - (cfl/3) *(0.5*(U'(t,i+1)+U'(t,i))-0.5*(U'(t,i-1)+U'(t,i)))
      % U**(t+1,i)= U'(t,i) - (cfl/2) *(0.5*(U*(t,i+1)+U*(t,i))-0.5*(U*(t,i-1)+U*(t,i)))
      % U(t+1,i)  = U'(t,i) - (cfl)   *(0.5*(U**(t,i+1)+U**(t,i))-0.5*(U**(t,i-1)+U**(t,i)))

          
          Hstar(t+1,i) = H(t,i) - (dt/(dx*3))*0.5*((u+((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))))*((H(t,i+1)+H(t,i))           -(H(t,i-1)+H(t,i)))          +(h+((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))))*((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))));         
          H2star(t+1,i)= H(t,i) - (dt/(dx*2))*0.5*((u+((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))))*((Hstar(t,i+1)+Hstar(t,i))   -(Hstar(t,i-1)+Hstar(t,i)))  +(h+((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))))*((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))));         
          H(t+1,i)     = H(t,i) - (dt/dx)    *0.5*((u+((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))))*((H2star(t,i+1)+H2star(t,i)) -(H2star(t,i-1)+H2star(t,i)))+(h+((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))))*((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))));
      
          Ustar(t+1,i) = U(t,i) - (dt/(dx*3))*0.5*((u+((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))))*((U(t,i+1)+U(t,i))           -(U(t,i-1)+U(t,i)))          +g*((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))));         
          U2star(t+1,i)= U(t,i) - (dt/(dx*2))*0.5*((u+((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))))*((Ustar(t,i+1)+Ustar(t,i))   -(Ustar(t,i-1)+Ustar(t,i)))  +g*((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))));         
          U(t+1,i)     = U(t,i) - (dt/dx)    *0.5*((u+((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))))*((U2star(t,i+1)+U2star(t,i)) -(U2star(t,i-1)+U2star(t,i)))+g*((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))));
      
      end
%
%% Periodic boundary condition is used which means U(t,0) = u(t, L) is true for all t
%

          Hstar(t+1,1) = H(t,1) - (dt/(dx*3))*0.5*((u+((U(t,2)+U(t,1))-(U(t,1)+U(t,II-1))))*((H(t,2)+H(t,1))           -(H(t,II-1)+H(t,1)))          +(h+((H(t,2)+H(t,1))-(H(t,1)+H(t,II-1))))*((U(t,1)+U(t,II-1))-(U(t,1)+U(t,II-1))));
          H2star(t+1,1)= H(t,1) - (dt/(dx*2))*0.5*((u+((U(t,2)+U(t,1))-(U(t,1)+U(t,II-1))))*((Hstar(t,2)+Hstar(t,1))   -(Hstar(t,II-1)+Hstar(t,1)))  +(h+((H(t,2)+H(t,1))-(H(t,1)+H(t,II-1))))*((U(t,1)+U(t,II-1))-(U(t,1)+U(t,II-1))));         
          H(t+1,1)     = H(t,1) - (dt/dx)    *0.5*((u+((U(t,2)+U(t,1))-(U(t,1)+U(t,II-1))))*((H2star(t,2)+H2star(t,1)) -(H2star(t,II-1)+H2star(t,1)))+(h+((H(t,2)+H(t,1))-(H(t,1)+H(t,II-1))))*((U(t,1)+U(t,II-1))-(U(t,1)+U(t,II-1))));
      
          Ustar(t+1,1) = U(t,1) - (dt/(dx*3))*0.5*((u+((U(t,2)+U(t,1))-(U(t,1)+U(t,II-1))))*((U(t,2)+U(t,1))           -(U(t,II-1)+U(t,1)))          +g*((H(t,2)+H(t,1)-(H(t,1)+H(t,II-1)))));         
          U2star(t+1,1)= U(t,1) - (dt/(dx*2))*0.5*((u+((U(t,2)+U(t,1))-(U(t,1)+U(t,II-1))))*((Ustar(t,2)+Ustar(t,1))   -(Ustar(t,II-1)+Ustar(t,1)))  +g*((H(t,2)+H(t,1)-(H(t,1)+H(t,II-1)))));         
          U(t+1,1)     = U(t,1) - (dt/dx)    *0.5*((u+((U(t,2)+U(t,1))-(U(t,1)+U(t,II-1))))*((U2star(t,2)+U2star(t,1)) -(U2star(t,II-1)+U2star(t,1)))+g*((H(t,2)+H(t,1)-(H(t,1)+H(t,II-1)))));
                      
          Hstar(t+1,II) = H(t+1,1);                                          
          H2star(t+1,II)= H(t+1,1);          
          H(t+1,II)     = H(t+1,1); 
      
          Ustar(t+1,II) = U(t+1,1);        
          U2star(t+1,II)= U(t+1,1);         
          U(t+1,II)     = U(t+1,1);     
  end                                                                 % end loop for space grid points
elseif(modeleqns==2)
  for t = 1:TT-1                                                       % start time integration
      for i = 2 : II-1                                                 % loop over grid points in space
                                                                     % solving with the discretized equation
                                                                     % solving with the discretized equation
%        Hprime(t+1,i)= Hprime(t,i)-(dt/(dx*2))*(U*(Hprime(t,i+1)-Hprime(t,i-1))+H*(Uprime(t,i+1)-Uprime(t,i-1))); 
%        Uprime(t+1,i)= Uprime(t,i)-(dt/(dx*2))*(U*(Uprime(t,i+1)-Uprime(t,i-1))+g*(Hprime(t,i+1)-Hprime(t,i-1)));                                                                          

      
      % Runge-Kutta O(3) from <http://cires.colorado.edu/science/groups/pielke/classes/at7500/wrfarw.pdf>
      % U*(t+1,i) = U'(t,i) - (cfl/3) *(0.5*(U'(t,i+1)+U'(t,i))-0.5*(U'(t,i-1)+U'(t,i)))
      % U**(t+1,i)= U'(t,i) - (cfl/2) *(0.5*(U*(t,i+1)+U*(t,i))-0.5*(U*(t,i-1)+U*(t,i)))
      % U(t+1,i)  = U'(t,i) - (cfl)   *(0.5*(U**(t,i+1)+U**(t,i))-0.5*(U**(t,i-1)+U**(t,i)))

          
          Hstar(t+1,i) = H(t,i) - (dt/(dx*3))*0.5*(u*((H(t,i+1)+H(t,i))           -(H(t,i-1)+H(t,i)))          +h*((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))));         
          H2star(t+1,i)= H(t,i) - (dt/(dx*2))*0.5*(u*((Hstar(t,i+1)+Hstar(t,i))   -(Hstar(t,i-1)+Hstar(t,i)))  +h*((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))));         
          H(t+1,i)     = H(t,i) - (dt/dx)    *0.5*(u*((H2star(t,i+1)+H2star(t,i)) -(H2star(t,i-1)+H2star(t,i)))+h*((U(t,i+1)+U(t,i))-(U(t,i)+U(t,i-1))));
      
          Ustar(t+1,i) = U(t,i) - (dt/(dx*3))*0.5*(u*((U(t,i+1)+U(t,i))           -(U(t,i-1)+U(t,i)))          +g*((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))));         
          U2star(t+1,i)= U(t,i) - (dt/(dx*2))*0.5*(u*((Ustar(t,i+1)+Ustar(t,i))   -(Ustar(t,i-1)+Ustar(t,i)))  +g*((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))));         
          U(t+1,i)     = U(t,i) - (dt/dx)    *0.5*(u*((U2star(t,i+1)+U2star(t,i)) -(U2star(t,i-1)+U2star(t,i)))+g*((H(t,i+1)+H(t,i))-(H(t,i)+H(t,i-1))));
      
      end
%
%% Periodic boundary condition is used which means U(t,0) = u(t, L) is true for all t
%

          Hstar(t+1,1) = H(t,1) - (dt/(dx*3))*0.5*(u*((H(t,2)+H(t,1))           -(H(t,II-1)+H(t,1)))          +h*((U(t,1)+U(t,II-1))-(U(t,1)+U(t,II-1))));
          H2star(t+1,1)= H(t,1) - (dt/(dx*2))*0.5*(u*((Hstar(t,2)+Hstar(t,1))   -(Hstar(t,II-1)+Hstar(t,1)))  +h*((U(t,1)+U(t,II-1))-(U(t,1)+U(t,II-1))));         
          H(t+1,1)     = H(t,1) - (dt/dx)    *0.5*(u*((H2star(t,2)+H2star(t,1)) -(H2star(t,II-1)+H2star(t,1)))+h*((U(t,1)+U(t,II-1))-(U(t,1)+U(t,II-1))));
      
          Ustar(t+1,1) = U(t,1) - (dt/(dx*3))*0.5*(u*((U(t,2)+U(t,1))           -(U(t,II-1)+U(t,1)))          +g*((H(t,2)+H(t,1)-(H(t,1)+H(t,II-1)))));         
          U2star(t+1,1)= U(t,1) - (dt/(dx*2))*0.5*(u*((Ustar(t,2)+Ustar(t,1))   -(Ustar(t,II-1)+Ustar(t,1)))  +g*((H(t,2)+H(t,1)-(H(t,1)+H(t,II-1)))));         
          U(t+1,1)     = U(t,1) - (dt/dx)    *0.5*(u*((U2star(t,2)+U2star(t,1)) -(U2star(t,II-1)+U2star(t,1)))+g*((H(t,2)+H(t,1)-(H(t,1)+H(t,II-1)))));
                      
          Hstar(t+1,II) = H(t+1,1);                                          
          H2star(t+1,II)= H(t+1,1);          
          H(t+1,II)     = H(t+1,1); 
      
          Ustar(t+1,II) = U(t+1,1);        
          U2star(t+1,II)= U(t+1,1);         
          U(t+1,II)     = U(t+1,1);     
  end                                                                 % end loop for space grid points
end 
      
%
%% Plot the initial conditions and final solution of U and H on the same graphics
%

  figure(1);clf; axis([0  40  -2  2]); hold on;
                 plot( 0:1:II-1, U(TT,:)); hold on;              % I.C.
                 plot( 0:1:II-1, H(TT,:),'r'); hold on;                  % Final solution
                 legend ('Wind Field','Hieght Field');
                 if (modeleqns == 1)
                     title('Non-linearized Forward in time, Centered in space Scheme');
                 elseif(modeleqns == 2)
                     title('Linearized Forward in time, Centered in space Scheme');
                 end

  PT = (TT-1)/8;                                                    % plot the solution at selected time steps with 
                                                                    % different panels

  figure(2);clf;
                 subplot(3,3,1); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(0*PT+1,:)); hold on;       % 0   (I.C.)
                 plot(0:1:II-1, H(0*PT+1,:),'r'); title 't = 0';
                 subplot(3,3,2); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(1*PT+1,:)); hold on;       % 1T/8  
                 plot(0:1:II-1, H(1*PT+1,:),'r'); title 't = 1t/8';
                 subplot(3,3,3); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(2*PT+1,:)); hold on;       % 2T/8
                 plot(0:1:II-1, H(2*PT+1,:),'r');title 't = 2t/8';
                 subplot(3,3,4); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(3*PT+1,:));  hold on;      % 3T/8
                 plot(0:1:II-1, H(3*PT+1,:),'r');title 't = 3t/8';
                 subplot(3,3,5); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(4*PT+1,:)); hold on;       % 4T/8
                 plot(0:1:II-1, H(4*PT+1,:),'r');title 't = 4t/8';
                 subplot(3,3,6); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(5*PT+1,:)); hold on;       % 5T/8
                 plot(0:1:II-1, H(5*PT+1,:),'r');title 't = 5t/8';
                 subplot(3,3,7); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(6*PT+1,:)); hold on;       % 6T/8
                 plot(0:1:II-1, H(6*PT+1,:),'r');title 't = 6t/8';
                 subplot(3,3,8); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(7*PT+1,:)); hold on;       % 7T/8
                 plot(0:1:II-1, H(7*PT+1,:),'r');title 't = 7t/8';
                 subplot(3,3,9); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(8*PT+1,:)); hold on;       % 8T/8
                 plot(0:1:II-1, H(8*PT+1,:),'r');title 't = 8t/8';
                 legend ('Wind','Height',4);  
