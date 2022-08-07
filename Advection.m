%==================================================================================
%                              Advection equation
%==================================================================================
% This program solves the advection without forcing numerically. From the 
% reference, the advection equation is:
%         du/dt + c(du/dx) = 0
%
%----------------------------------------------------------------------------------
%
% @creator: Fuqing Zhang
% @updated: Skylar Hernandez
% @version: 2.0.0
% @date   : 23 Sept. 2009 
%
%----------------------------------------------------------------------------------
% Version 2.0.0:
%   This version has the advection formula solved under the upstream in
%   space forward in time, Leap Frog, and Runge Kutta O(3) Schemes.  
% 
%==================================================================================

  clear all;         % clean up graphic and memory

%
%% Defintion of key parameters
% 
%% Compare the difference between these two schemes for dx=250, 2500 
%% and dt=50, 100, 200
%
  
  L = 100000;         % The length of the domain in meters
  T = 3600  ;         % Total integration time in seconds
  c = 10    ;         % Constant advection speed in m/s

  dx = 2500 ;         % Grid spacing in meters
  dt = 50 ;           % Time step in seconds
  cfl = c*dt/dx;      % Courant number, CFL criteria needs cfl < 1

  II = L/dx +1;       % Total grid points in x direction
  TT = T/dt +1;       % Total time steps of integration



%
%% Initialize the 2-D matrix of U with 0s
%

  U(TT, II)     = 0.;
  diff(TT,II)   = 0.;
  Ustar(TT,II)  = 0.;
  U2star(TT,II) = 0.;
  Usln(TT,II)   = 0.;
   
%
%% Define inititial conditions of U
%

 scheme = input('Type (1) Discretized, (2) Upstream, (3) Runge-Kutta, (4) Leap Frog = ');
    
  for i = 1:1:II 
      U(1, i)      = sin( 2*(i-1)*dx*pi/L );   % a sin wave with one full wavelength
      Ustar(1, i)  = sin( 2*(i-1)*dx*pi/L );
      U2star(1, i) = sin( 2*(i-1)*dx*pi/L );
  end
  

%
%% Periodic boundary condition is used which means U(t,0) = u(t, L) is true for all t
%
%% Discretize the PDE using forward in time and centered in space scheme; however,
%% this scheme is always unstable which means the wave amplitude will increase 
%% indefintely as time goes on even if cfl<1.
%


%
%% Numerical integration or prediction
%

  for t = 1:TT-1                                                     % start time integration
    for i = 2 : II-1                                                 % loop over grid points in space
      Usln(t+1,i) = sin( 2*((i-1)*dx-c*t*dt)*pi/L);  
      if (scheme == 1)
          U(t+1,i) = U(t,i) - (cfl/2) * ( U(t,i+1) - U(t,i-1)  );    % solving with the discretized equation
          diff(t+1,i) = sqrt((Usln(t+1,i)-U(t+1,i))^2);
      elseif (scheme == 2)
          U(t+1,i) = U(t,i) - (cfl)   * ( U(t,i)   - U(t,i-1)  );    %% Forward in time, upsrteam in space
          diff(t+1,i) = sqrt((Usln(t+1,i)-U(t+1,i))^2);
      elseif (scheme == 3)
          % Runge-Kutta O(3) from <http://cires.colorado.edu/science/groups/pielke/classes/at7500/wrfarw.pdf>
          % U*(t+1,i) = U(t,i) - (cfl/3) *(0.5*(U(t,i+1)+U(t,i))-0.5*(U(t,i-1)+U(t,i)))
          % U**(t+1,i)= U(t,i) - (cfl/2) *(0.5*(U*(t,i+1)+U*(t,i))-0.5*(U*(t,i-1)+U*(t,i)))
          % U(t+1,i)  = U(t,i) - (cfl)   *(0.5*(U**(t,i+1)+U*(t,i))-0.5**(U*(t,i-1)+U*(t,i)))
          Ustar(t+1,i) = U(t,i) - (cfl/3) *(0.5*(U(t,i+1)+U(t,i))-0.5*(U(t,i-1)+U(t,i)));         
          U2star(t+1,i)= U(t,i) - (cfl/2) *(0.5*(Ustar(t,i+1)+Ustar(t,i))-0.5*(Ustar(t,i-1)+Ustar(t,i)));
          U(t+1,i)     = U(t,i) - (cfl)   *(0.5*(U2star(t,i+1)+U2star(t,i))-0.5*(U2star(t,i-1)+U2star(t,i)));
          diff(t+1,i) = sqrt((Usln(t+1,i)-U(t+1,i))^2);
      else
      end
    end                                                              % end loop for space grid points

    if (scheme == 1)
      U(t+1,1) = U(t,1) - (cfl/2) * ( U(t,2)   - U(t,II-1) );        % with periodic B.C., U(t,2) is equivalent to U(t,II+1)
      U(t+1,II) = U(t+1, 1);                                         % apply periodic B.C. for the end point II;
    elseif (scheme == 2)
      U(t+1,1) =  U(t,1) - (cfl)  *  ( U(t,2)   - U(t,II-1) );       % with periodic B.C., U(t,2) is equivalent to U(t,II+1)
      U(t+1,II) = U(t+1, 1);                                         % apply periodic B.C. for the end point II;                                                                 %% Forward in time, upsrteam in space
    elseif (scheme == 3) 
      Ustar(t+1,1)  = U(t,i) - cfl/3*(0.5*(U(t,i+1)+U(t,i))-0.5*(U(t,i-1)+U(t,i)));        
      U2star(t+1,1) = U(t,i) - cfl/2*(0.5*(Ustar(t,i+1)+Ustar(t,i))-0.5*(Ustar(t,i-1)+Ustar(t,i)));
      U(t+1,1)      = U(t,i) - cfl*  (0.5*(U2star(t,i+1)+U2star(t,i))-0.5*(U2star(t,i-1)+U2star(t,i)));   
      Ustar(t+1,II) = U(t,i) - cfl/3*(0.5*(U(t,i+1)+U(t,i))-0.5*(U(t,i-1)+U(t,i)));   
      U2star(t+1,II)= U(t,i) - cfl/2*(0.5*(Ustar(t,i+1)+Ustar(t,i))-0.5*(Ustar(t,i-1)+Ustar(t,i)));
      U(t+1,II)     = U(t,i) - cfl*  (0.5*(U2star(t,i+1)+U2star(t,i))-0.5*(U2star(t,i-1)+U2star(t,i)));  
    else
    end
 
      Usln(t+1,1) = sin( 2*((i-1)*dx-c*t*dt)*pi/L);  
      Usln(t+1,II) = sin( 2*((i-1)*dx-c*t*dt)*pi/L);             % end loop for time integration
      diff(t+1,1) = sqrt((Usln(t+1,1)-U(t+1,1))^2);
      diff(t+1,II) = sqrt((Usln(t+1,II)-U(t+1,II))^2);
  end

if (scheme == 4)
U(1,i)   = sin( 2*(i-1)*dx*pi/L );                                %% IC to set the U(1,i)

  for t = 2:TT-1                                                  % start time integration
    for i = 2 : II-1                                              % loop over grid points in space
       U(t+1,i) = U(t-1,i) - (cfl) * ( U(t,i+1) - U(t,i-1)  );    %% Leap Frog
       diff(t+1,i) = sqrt((Usln(t+1,i)-U(t+1,i))^2);
    end

    
    U(t+1,1) = U(t-1,i) - (cfl) * ( U(t,i+1) - U(t,i-1)  );       % with periodic B.C., U(t,2) is equivalent to U(t,II+1)
    U(t+1,II) = U(t+1, 1);                                        % apply periodic B.C. for the end point II;
    diff(t+1,1) = sqrt((Usln(t+1,1)-U(t+1,1))^2);
    diff(t+1,II) = sqrt((Usln(t+1,II)-U(t+1,II))^2);
  end 
else
end
%
%% Plot the initial conditions and final solution of U on the same graphics
%

  figure(1);clf; axis([0  40  -2  2]); hold on;
                 plot( 0:1:II-1, U(1,:),'r' ); hold on;              % I.C.
                 plot( 0:1:II-1, U(TT,:)); hold on;                  % Final solution
                 legend ('Initial Condition','Final Solution');
                 if (scheme == 1)
                    title('Discretized Scheme');
                 elseif (scheme == 2)
                    title('Forward in Time & Upstream in Space Scheme');     
                 elseif (scheme == 3)
                    title('Runge-Kutta O(3) Scheme'); 
                 elseif (scheme == 4)
                    title('Leap Frog Scheme'); 
                 else
                 end

  PT = (TT-1)/8;                                                    % plot the solution at selected time steps with 
                                                                    % different panels

  figure(2);clf;
                 subplot(3,3,1); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(0*PT+1,:)); hold on;       % 0   (I.C.)
                 plot(0:1:II-1, U(0*PT+1,:),'r'); title 't = 0';
                 subplot(3,3,2); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(1*PT+1,:)); hold on;       % 1T/8  
                 plot(0:1:II-1, Usln(1*PT+1,:),'r'); title 't = 1t/8';
                 subplot(3,3,3); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(2*PT+1,:)); hold on;       % 2T/8
                 plot(0:1:II-1, Usln(2*PT+1,:),'r');title 't = 2t/8';
                 subplot(3,3,4); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(3*PT+1,:));  hold on;      % 3T/8
                 plot(0:1:II-1, Usln(3*PT+1,:),'r');title 't = 3t/8';
                 subplot(3,3,5); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(4*PT+1,:)); hold on;       % 4T/8
                 plot(0:1:II-1, Usln(4*PT+1,:),'r');title 't = 4t/8';
                 subplot(3,3,6); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(5*PT+1,:)); hold on;       % 5T/8
                 plot(0:1:II-1, Usln(5*PT+1,:),'r');title 't = 5t/8';
                 subplot(3,3,7); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(6*PT+1,:)); hold on;       % 6T/8
                 plot(0:1:II-1, Usln(6*PT+1,:),'r');title 't = 6t/8';
                 subplot(3,3,8); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(7*PT+1,:)); hold on;       % 7T/8
                 plot(0:1:II-1, Usln(7*PT+1,:),'r');title 't = 7t/8';
                 subplot(3,3,9); axis([0  40  -2  2]); hold on;
                 plot(0:1:II-1, U(8*PT+1,:)); hold on;       % 8T/8
                 plot(0:1:II-1, Usln(8*PT+1,:),'r');title 't = 8t/8';
                 legend ('Forecast','Truth',4);  
 
  figure(3);clf;
                 subplot(3,3,1); plot(0:1:II-1, diff(0*PT+1,:)); hold on;       % 0   (I.C.)
                 title 't = 0';
                 subplot(3,3,2); plot(0:1:II-1, diff(1*PT+1,:)); hold on;       % 1T/8  
                 title 't = 1t/8';
                 subplot(3,3,3); plot(0:1:II-1, diff(2*PT+1,:)); hold on;       % 2T/8
                 title 't = 2t/8';
                 subplot(3,3,4); plot(0:1:II-1, diff(3*PT+1,:));  hold on;      % 3T/8
                 title 't = 3t/8';
                 subplot(3,3,5); plot(0:1:II-1, diff(4*PT+1,:)); hold on;       % 4T/8
                 title 't = 4t/8';
                 subplot(3,3,6); plot(0:1:II-1, diff(5*PT+1,:)); hold on;       % 5T/8
                 title 't = 5t/8';
                 subplot(3,3,7); plot(0:1:II-1, diff(6*PT+1,:)); hold on;       % 6T/8
                 title 't = 6t/8';
                 subplot(3,3,8); plot(0:1:II-1, diff(7*PT+1,:)); hold on;       % 7T/8
                 title 't = 7t/8';
                 subplot(3,3,9); plot(0:1:II-1, diff(8*PT+1,:)); hold on;       % 8T/8
                 title 't = 8t/8';
                 legend ('RMS', 1); 

                 
%     (3) Calculate the root-mean-square difference between the numerical solution and analytical solution as a function of time and plotted in Figure 3.
