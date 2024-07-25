%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Daksh Maheshwari %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% EECS 516 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% MRI Project %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BME/EECS516
% Parts 7-11
% Oct 2022; Last revision: Oct-30-2022

plotdir = './generated_plots';
if ~exist(plotdir,'dir'), mkdir(plotdir); end
% variable to save plot (1 to save plot, 0 otherwise)
saveplot = 0;

% Single point object at (x,y,z) = (2,2,0) cm;
% Point object has T1 of 1000 ms, T2 of 100 ms
obj_x = 2;
obj_y = 2;
obj_z = 0;
obj_T1 = 1000;
obj_T2 = 100;

obj_n = length(obj_x); % Determine number of objects


%% Define simulation constants

% Physical constants
gambar = 42570; % Gamma/2pi in kHz/T
gam = gambar*2*pi; % Gamma in kiloradians/T

% Simulation values
dt = 0.05; % Time step for simulation, ms (50 us step size)
te = 10.0; % Echo time, ms
endtime = 15; % Total runtime of simulation, ms
time = [0:dt:endtime]'; % Vector containing each time step, ms (size #timepoints x 1)
totalTimepoints = length(time); % Number of time points for simulation

% Initialize B vectors, the effective (x,y,z) applied magnetic field
% Vectors define applied magnetic field at time tp_n for object obj_n
bx = zeros([totalTimepoints obj_n]);
by = zeros([totalTimepoints obj_n]);
bz = zeros([totalTimepoints obj_n]);

% Define a 90 RF pulse
rf90pw = 2; % Pulse width in ms
sincper = rf90pw/4; % in ms (this is the sinc stretch paramters
rf_timepoints = rf90pw/dt; % Number of simulation steps for RF
rf_time = [-(rf_timepoints-1)/2:(rf_timepoints-1)/2]'.*dt; % Time vector for creating sinc, centered at 0
rf_shape = hanning(rf_timepoints).*sinc(rf_time./sincper); % Sinc waveform shape with hanning window, with amplitude 1

% B1 amplitude alpha = gamma * A * T
alpha = pi/2; % flip angle 90 degrees
T = sum(rf_shape) * dt;
% T = sincper;
A = alpha/(gam * T);

rf_amplitude90 = A; % REPLACE 0 with amplitude of the RF pulse here, in T

% Scale rf_shape by a_rf90 (amplitude), then fill the remainder of the time with zeros
b1_90 = rf_amplitude90.*[rf_shape; zeros([totalTimepoints-rf_timepoints 1])];

%% Important values
Nx = 64;
Ny = 32;
Tread = 6.4; % in ms
Ty = 2; % in ms
Tx = 2; % in ms
FOVx = 10; % in cm
FOVy = 5; % in cm
Wkx = Nx/FOVx;
Wky = Ny/FOVy;
deltaz = 1; % in cm


%% Create gradients
% Create gz
slThick = 1; % Slick thickness in cm
gz1_a = 1 / (gambar * sincper * deltaz); % REPLACE 0 with amplitude of gz1 in T/cm
gz1_pw = rf90pw; % Match the width of gz1 to the RF pulse
gz2_a = -gz1_a;
gz2_pw = rf90pw/2;

% Create gz with positve area gz1_a*gz1_pw, followed by negative area gz2_a*gz2pw
% gz step size is dt, with amplitude values in T/cm
gz = (time < gz1_pw) .* gz1_a ...
 + (time >= gz1_pw).*(time < (gz1_pw+gz2_pw)) .* gz2_a;

% Create gx
%%%
% for symmetrical sample acquisition
gx1_a = -1*(Wkx/2) / (gambar * Tx);
gx2_a = Wkx/(gambar*Tread);
gx1_start = 5.8;
gx1_stop = 7.8;
gx2_stop = 14.2;
gx = (time >= gx1_start & time < gx1_stop) .* gx1_a ...
 + (time >= gx1_stop & time <= gx2_stop) .*gx2_a;
%%%
% Create gy
%%%
gy_max = Wky/(2*gambar*Ty);
gy_start = 5.8;
gy_stop = 7.8;
delta_gy = 1 / (gambar * FOVy * Ty);
gy_timevec = (time >= gy_start & time <= gy_stop);
% it is giving 33 length vector, making it 32
gy_a = -gy_max:delta_gy:(gy_max-delta_gy);
%%%


%% 7 & 8
M = zeros([Ny,Nx]);
npe = length(gy_a);
deltaw = gam.*gz*1;
%%%
for pe = 1:npe % phase encode loop
    disp(sprintf('PE %d of %d', pe, npe));
    % assign gy amplitudes here
    gy = gy_timevec .* gy_a(pe);

    % set beff here
    % bx,by,bz should be of size npts x nobj and should describe the
    % beff(t) in x,y,z direction in the *rotating* frame for each point
    % in space and time. These b matrices should have units of Tesla
    bx = b1_90;
    bz = gx .* obj_x + gy .* obj_y + gz .* obj_z;
    
    % Set the initial magnitization as a unit vector (0,0,1) for # obj
    m0 = [0 0 1]'*ones([1 obj_n]);
    
    % Do Bloch equation simulation
    [mx,my,mz] = blochsim_516(m0,bx,by,bz,obj_T1,obj_T2,dt);
    
    % mx,my,mz will have size totalTimepoints x obj_n and will describe M(t)
    % in x,y,z direction in the rotating frame for each point in space
    % Do plots for part 6
    mxy = mx + my*i;
    data_acq = mxy(find(time == gx1_stop): find(time == gx2_stop), :);
    sampled_mxy = data_acq(2:2:end, :);


% place each phase encode line into the acquisition matrix
% specifically, load M(kx,ky)
M(pe, :) = sum(sampled_mxy, 2).';
%
end % end of phase encode loop


f7a = figure('DefaultAxesFontSize',16);
f7a.WindowState = "maximized";
kxpos = [-Nx/2:Nx/2-1] / FOVx; % vector of kx locations
kypos = [-Ny/2:Ny/2-1] / FOVy; % vector of ky locations
imagesc(kxpos,kypos,real(M)); colormap gray; axis('image'); axis('xy');
xlabel('kx (cm^{-1})', 'FontSize', 16);
ylabel('ky (cm^{-1})', 'FontSize', 16);
title('7. real M(kx, ky)', 'FontSize', 21);

if saveplot == 1
    % saveas(f7a, './generated_plots/part_7_real.jpeg')
    pause(1);
    print(f7a, './generated_plots/part_7_real', '-dpng', '-r300');
end

f7b = figure('DefaultAxesFontSize',16);
f7b.WindowState = "maximized";
imagesc(kxpos,kypos,imag(M)); colormap gray; axis('image'); axis('xy');
xlabel('kx (cm^{-1})', 'FontSize', 16);
ylabel('ky (cm^{-1})', 'FontSize', 16);
title('7. imag M(kx, ky)', 'FontSize', 21);

if saveplot == 1
    % saveas(f7b, './generated_plots/part_7_imag.jpeg')
    pause(1);
    print(f7b, './generated_plots/part_7_imag', '-dpng', '-r300');
end

% reconstruct the image here
im = ift2(M);
f8 = figure('DefaultAxesFontSize',16);
f8.WindowState = "maximized";
xpos = [-Nx/2:Nx/2-1]/Nx * FOVx; % vector of x locations
ypos = [-Ny/2:Ny/2-1]/Ny * FOVy; % vector of y locations
imagesc(xpos,ypos,abs(im)); colormap gray; axis('image'); axis('xy');
xlabel('x (cm)', 'FontSize', 16);
ylabel('y (cm)', 'FontSize', 16);
title('8. reconstructed image abs(image(x,y))');

if saveplot == 1
    % saveas(f8, './generated_plots/part_8.jpeg');
    pause(1);
    print(f8, './generated_plots/part_8', '-dpng', '-r300');
end
%% 9
% setting point object at (2,4,0)
obj_y = 4;
M = zeros([Ny,Nx]);
npe = length(gy_a);
deltaw = gam.*gz*1;
%%%
for pe = 1:npe % phase encode loop
    disp(sprintf('PE %d of %d', pe, npe));
    % disp(sprintf('PE %d of %d, Slice %2', pe, npe,slice))
    %
    % assign gy amplitudes here
    gy = gy_timevec .* gy_a(pe);
    %
    % set beff here
    % bx,by,bz should be of size npts x nobj and should describe the
    % beff(t) in x,y,z direction in the *rotating* frame for each point
    % in space and time. These b matrices should have units of Tesla
    bx = b1_90;
    bz = gx .* obj_x + gy .* obj_y + gz .* obj_z;
    
    % Set the initial magnitization as a unit vector (0,0,1) for # obj
    m0 = [0 0 1]'*ones([1 obj_n]);
    
    % Do Bloch equation simulation
    [mx,my,mz] = blochsim_516(m0,bx,by,bz,obj_T1,obj_T2,dt);
    
    % mx,my,mz will have size totalTimepoints x obj_n and will describe M(t)
    % in x,y,z direction in the rotating frame for each point in space
    % Do plots for part 6
    mxy = mx + my*i;
    data_acq = mxy(find(time == gx1_stop): find(time == gx2_stop), :);
    sampled_mxy = data_acq(2:2:end, :);

   
% place each phase encode line into the acquisition matrix
% specifically, load M(kx,ky)
M(pe, :) = sum(sampled_mxy, 2).';
%
end % end of phase encode loop


f9a = figure('DefaultAxesFontSize',16);
f9a.WindowState = "maximized";
kxpos = [-Nx/2:Nx/2-1] / FOVx; % vector of kx locations
kypos = [-Ny/2:Ny/2-1] / FOVy; % vector of ky locations
imagesc(kxpos,kypos,real(M)); colormap gray; axis('image'); axis('xy');
xlabel('kx (cm^{-1})', 'FontSize', 16);
ylabel('ky (cm^{-1})', 'FontSize', 16);
title('9. real M(kx, ky)', 'FontSize', 21);
if saveplot == 1
    % saveas(f9a, './generated_plots/part_9_real.jpeg')
    pause(1);
    print(f9a, './generated_plots/part_9_real', '-dpng', '-r300');
end


f9b = figure('DefaultAxesFontSize',16);
f9b.WindowState = "maximized";
imagesc(kxpos,kypos,imag(M)); colormap gray; axis('image'); axis('xy');
xlabel('kx (cm^{-1})', 'FontSize', 16);
ylabel('ky (cm^{-1})', 'FontSize', 16);
title('9. imag M(kx, ky)', 'FontSize', 21);

if saveplot == 1
    % saveas(f9b, './generated_plots/part_9_imag.jpeg');
    pause(1);
    print(f9b, './generated_plots/part_9_imag', '-dpng', '-r300');
end

% reconstruct the image here
im = ift2(M);

f9c = figure('DefaultAxesFontSize',16);
f9c.WindowState = "maximized";
xpos = [-Nx/2:Nx/2-1]/Nx * FOVx; % vector of x locations
ypos = [-Ny/2:Ny/2-1]/Ny * FOVy; % vector of y locations
imagesc(xpos,ypos,abs(im)); colormap gray; axis('image'); axis('xy');
xlabel('x (cm)', 'FontSize', 16);
ylabel('y (cm)', 'FontSize', 16);
title('9. reconstructed image abs(image(x,y))');

if saveplot == 1
    % saveas(f9c, './generated_plots/part_9_image.jpeg')
    pause(1);
    print(f9c, './generated_plots/part_9_image', '-dpng', '-r300');
end


%% part 10-11 (loading object22 containing 1345 spins)
disp('loading object 22');
load object22;


%% Loop through slices and phase encodes
obj_n = length(obj_x); % Determine number of objects
% taking outer product for b1_90 to create b1_90 as 301x1345
b1_90 = b1_90*ones(1, obj_n);

% re-initializing bx,by,bz to match dims
% according to 1345 spins
bx = zeros([totalTimepoints obj_n]);
by = zeros([totalTimepoints obj_n]);
bz = zeros([totalTimepoints obj_n]);


for slice = 1:2 % slice loop
    M = zeros([Ny,Nx]);
    npe = length(gy_a);
    deltaw = gam.*gz*1;
    for pe = 1:npe % phase encode loop
        disp(sprintf('PE %d of %d', pe, npe));
        % disp(sprintf('PE %d of %d, Slice %2', pe, npe,slice))
        %
        % assign gy amplitudes here
        gy = gy_timevec .* gy_a(pe);
        %
        % set beff here
        % bx,by,bz should be of size npts x nobj and should describe the
        % beff(t) in x,y,z direction in the *rotating* frame for each point
        % in space and time. These b matrices should have units of Tesla
        if slice == 1
        bx = b1_90;
        else
        bx = b1_90 .* real(exp(-i*deltaw.*time));
        by = b1_90 .* imag(exp(-i*deltaw.*time));
        end
        bz = gx * obj_x + gy * obj_y + gz * obj_z;
        
        % Set the initial magnitization as a unit vector (0,0,1) for # obj
        m0 = [0 0 1]'*ones([1 obj_n]);
        
        % Do Bloch equation simulation
        [mx,my,mz] = blochsim_516(m0,bx,by,bz,obj_T1,obj_T2,dt);
        
        % mx,my,mz will have size totalTimepoints x obj_n and will describe M(t)
        % in x,y,z direction in the rotating frame for each point in space
        % Do plots for part 6
        mxy = mx + my*i;
        data_acq = mxy(find(time == gx1_stop): find(time == gx2_stop), :);
        sampled_mxy = data_acq(2:2:end, :);
    
        % place each phase encode line into the acquisition matrix
        % specifically, load M(kx,ky)
        M(pe, :) = sum(sampled_mxy, 2).';
        %
    end % end of phase encode loop


    f10_11a = figure('DefaultAxesFontSize',16);
    f10_11a.WindowState = "maximized";
    kxpos = [-Nx/2:Nx/2-1] / FOVx; % vector of kx locations
    kypos = [-Ny/2:Ny/2-1] / FOVy; % vector of ky locations
    imagesc(kxpos,kypos,abs(M)); colormap gray; axis('image'); axis('xy');
    xlabel('kx (cm^{-1})', 'FontSize', 16);
    ylabel('ky (cm^{-1})', 'FontSize', 16);
    if slice == 1
        title('10. abs(M(kx,ky))', 'FontSize', 21);
    
        if saveplot == 1
            % saveas(f10_11a, './generated_plots/part_10_abs(M).jpeg')
            pause(1);
            print(f10_11a, './generated_plots/part_10_abs(M)', '-dpng', '-r300');
        end
        
    else
        title('11. abs(M(kx,ky))', 'FontSize', 21);
        if saveplot == 1
            % saveas(f10_11a, './generated_plots/part_11_abs(M).jpeg')
            pause(1);
            print(f10_11a, './generated_plots/part_11_abs(M)', '-dpng', '-r300');
        end
    end

    % reconstruct the image here
    im = ift2(M);
    f10_11b = figure('DefaultAxesFontSize',16);
    f10_11b.WindowState = "maximized";
    xpos = [-Nx/2:Nx/2-1]/Nx * FOVx; % vector of x locations
    ypos = [-Ny/2:Ny/2-1]/Ny * FOVy; % vector of y locations
    imagesc(xpos,ypos,abs(im)); colormap gray; axis('image'); axis('xy');
    xlabel('x (cm)', 'FontSize', 16);
    ylabel('y (cm)', 'FontSize', 16);
    if slice == 1
        title('10. reconstructed image abs(image(x,y))', 'FontSize', 21);
        if saveplot == 1
            % saveas(f10_11b, './generated_plots/part_10_abs(image).jpeg');
            pause(1);
            print(f10_11b, './generated_plots/part_10_abs(image)', '-dpng', '-r300');
        end
    else
        title('11. reconstructed image abs(image(x,y))', 'FontSize', 21);
        if saveplot == 1
            % saveas(f10_11b, './generated_plots/part_11_abs(image).jpeg')
            pause(1);
            print(f10_11b, './generated_plots/part_11_abs(image)', '-dpng', '-r300');
        end
    end
end