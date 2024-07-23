BME/EECS516
% MRI Project Template

% Other m-files required: ift2, ift, ft2, ft, blochsim_516
% Subfunctions: none
% MAT-files required: object18.mat

% Oct 2022; Last revision: Oct-30-2022

%% Select whether to load complex 2D object or create simple point object
complexobj = 0;
if complexobj
   % 2D Object for reconstruction
  load object22;
else
  % Single point object at (x,y,z) = (2,2,0) cm;
  % Point object has T1 of 1000 ms, T2 of 100 ms
  obj_x = 2;
  obj_y = 2;
  obj_z = 0;
  obj_T1 = 1000;
  obj_T2 = 100;
end

obj_n = length(obj_x); % Determine number of objects

%% Define simulation constants
% Physical constants
gambar = 42570;               % Gamma/2pi in kHz/T
gam = gambar*2*pi;            % Gamma in kiloradians/T

% Simulation values
dt = 0.05;                    % Time step for simulation, ms (50 us step size)
te = 10.0;                    % Echo time, ms
endtime = 15;                 % Total runtime of simulation, ms
time = [0:dt:endtime]';       % Vector containing each time step, ms (size #timepoints x 1)
totalTimepoints = length(time);          % Number of time points for simulation

% Initialize B vectors, the effective (x,y,z) applied magnetic field
% Vectors define applied magnetic field at time tp_n for object obj_n

bx = zeros([totalTimepoints obj_n]);
by = zeros([totalTimepoints obj_n]);
bz = zeros([totalTimepoints obj_n]);

% Define a 90 RF pulse
rf90pw = 2;                  % Pulse width in ms
sincper = rf90pw/4;          % in ms (this is the sinc stretch paramters
rf_timepoints = rf90pw/dt;   % Number of simulation steps for RF
rf_time = [-(rf_timepoints-1)/2:(rf_timepoints-1)/2]'.*dt; % Time vector for creating sinc, centered at 0
rf_shape = hanning(rf_timepoints).*sinc(rf_time./sincper); % Sinc waveform shape with hanning window, with amplitude 1

rf_amplitude90 = 0;          % REPLACE 0 with amplitude of the RF pulse here, in T

% Scale rf_shape by a_rf90 (amplitude), then fill the remainder of the time with zeros
b1_90 = rf_amplitude90.*[rf_shape; zeros([totalTimepoints-rf_timepoints 1])];

%% Create gradients
% Create gz
slThick = 1;             % Slick thickness in cm
gz1_a = 0;                % REPLACE 0 with amplitude of gz1 in T/cm
gz1_pw = rf90pw;           % Match the width of gz1 to the RF pulse
gz2_a = -gz1_a;
gz2_pw = rf90pw/2;

% Create gz with positve area gz1_a*gz1_pw, followed by negative area gz2_a*gz2pw
% gz step size is dt, with amplitude values in T/cm
gz =  (time < gz1_pw) .* gz1_a ...
       + (time >= gz1_pw).*(time < (gz1_pw+gz2_pw)) .* gz2_a;

% Create gx

% Create gy

%% Loop through slices and phase encodes
% for slice = 1:2  % slice loop
    % 
    % modify B1 here for other slices 
    % 

    % for pe = 1:npe % phase encode loop
        % disp(sprintf('PE %d of %d, Slice %2', pe, npe,slice))
        %
	    % assign gy amplitudes here 
	    %

	    % set beff here
	    % bx,by,bz should be of size npts x nobj and should describe the 
        % beff(t) in x,y,z direction in the *rotating* frame for each point
        % in space and time.  These b matrices should have units of Tesla
        
%         bx = ??
%         by = ??
%         bz = ??
        
        % Set the initial magnitization as a unit vector (0,0,1) for # obj
        m0 = [0 0 1]'*ones([1 obj_n]);
        
        % Do Bloch equation simulation
        [mx,my,mz] = blochsim_516(m0,bx,by,bz,obj_T1,obj_T2,dt);
        
        % mx,my,mz will have size totalTimepoints x obj_n and will describe M(t)  
        % in x,y,z direction in the rotating frame for each point in space

	    % Plot magnetization (comment out for parts 6-11)
            subplot(3,1,1)
            plot(time,mx);
            xlabel('time (ms)');
            ylabel('Mx');
            axis([0 15 -1 1]);
            title('put title here');
            
            subplot(3,1,2)
            plot(time,my);
            xlabel('time (ms)');
            ylabel('My');
            axis([0 15 -1 1]);
            
            subplot(3,1,3)
            plot(time,mz);
            xlabel('time (ms)');
            ylabel('Mz');
            axis([0 15 -1 1]);

	    % Do plots for part 6
	    % defined received (sampled) signal from mx,my,mz
        
        % now do plots
        % xpos = [-nread/2:nread/2-1]/nread*FOVx;
        % plot(xpos,abs(ift(sig)));

	%
	% place each phase encode line into the acquisition matrix
	% specifically, load M(kx,ky)
	%

    % end  % end of phase encode loop
    %
    % 
    % show images for parts 7-10
    % kxpos = ???  % vector of kx locations
    % kypos = ???  % vector of ky locations
    % imagesc(kxpos,kypos,real(M)); colormap gray; axis('image'); axis('xy')
    % xlabel('kx (units?)');
    % ylabel('ky (units?)');
    % title('put title here')
    % disp 'Press any key to continue...'; pause
    
    % imagesc(kxpos,kypos,imag(M)); colormap gray; axis('image'); axis('xy')
    % xlabel('kx (units?)');
    % ylabel('ky (units?)');
    % title('put title here')
    % disp 'Press any key to continue...'; pause

    % reconstruct the image here

    % xpos = ??  % vector of x locations
    % ypos = ??  % vector of y locations
    % imagesc(xpos,ypos,abs(im)); colormap gray; axis('image'); axis('xy')
    % xlabel('x (units?)');
    % ylabel('y (units?)');
    % title('put title here')

% end  % end of slice loop

