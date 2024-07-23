%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Daksh Maheshwari %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% EECS 516 %%%%%%%%%%%%%%%%%%
%%%%%%% Ultrasound beamforming project %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%	This script requires the variables
%	Data [ntime, nelem]	RF signals from each array element
%	f0			Transducer center frequency (MHz)
%	fs			sampling frequency (MHz)
%	c			speed of sound (mm/usec)
%	dx			transducer element spacing (mm)
%
% Get needed variables
%
load data22;
f0 = 5;		    % MHz
fs = 20;	    % MHz
c = 1.54;	    % mm/us
lambda = c/f0;  % in mm
% 
% Output Array Parameters
%
% Determine the array spacing dx in mm
dx = lambda/2;
w0 = 2*pi*f0;

deltat=1/fs;
[ntime, nelem] = size(Data);		% # time samples, # array elements
disp(sprintf('f0=%g MHz, deltat=%g usec, dx=%g mm', f0, deltat, dx))
disp(sprintf('# of Time Samples=%g,  # of Array Elements=%g',ntime,nelem))

%
% --> QUESTION a. <--
% Make a "wavefield" plot of the raw Data
% Comment this out while you debug other parts of your program
%
showimage3(Data, 4)  % plotting the transpose here...time along x, tranducer x along y.
xlabel('element');
ylabel('time');
title('part (a): wavefield plot of the raw data');
disp 'hit key', pause

%
% --> QUESTION b. <--
% Compute the number of total beams and beam spacing
% I used variables called
%	nbeam = number of beams
%	sin_theta = vector of beam positions
%
a = (nelem-1) * dx/2;
delta_sin_theta = lambda/(2*a);

max_sin_theta = max(sin(pi/4));
min_sin_theta = min(sin(-pi/4));

nbeam = round((max_sin_theta - min_sin_theta) / delta_sin_theta);

disp(sprintf('nbeam = %g', nbeam))

sin_theta = min_sin_theta:delta_sin_theta:max_sin_theta;

t = (1:ntime) * deltat;     % time vector
r = (t*c) / 2;      % r vector
xn = linspace(-a, a, nelem);        % nth element w.r.t. center element.

%
% --> QUESTION f. <--
%
databb = zeros(ntime,nelem);	% Phase-Shifted Baseband data
lpf = abs(round(-ntime/2):round(ntime/2)-1) < round(ntime/4);

exp_demod = exp(1i*w0.*t).';

for ie = 1:nelem
    databb(:, ie) = ift(ft(Data(:, ie) .* exp_demod) .* lpf.');
end

% Plot and Show data 
%plot...
%showimage3(abs(databb), 4, 40) 
%disp 'hit key', pause

showimage3(real(databb), 4, 40) 
xlabel('time');
ylabel('element');
title('part (f): real baseband data');
disp 'hit key', pause

showimage3(imag(databb), 4, 40) 
xlabel('time');
ylabel('element');
title('part (f): imag baseband data');
disp 'hit key', pause

%
% Create room for answers and precompute information
%	This section merely makes matlab allocate memory before doing
%	any other processing. It is not a necessary section, but will
%	make Matlab run a little faster since it doesn't have to reallocate
%	space every time you add another beam
%
rsdata = zeros(ntime,nbeam);	            % r-sin(theta) data buffer
bbshift = zeros(ntime,nelem);	            % Phase-Shifted Baseband data

bbshift_rsdata = zeros(ntime, nbeam);       % r-sin(theta) for phase shifted baseband
bbshift_rsdata_han = zeros(ntime, nbeam);   % r-sin(theta) for phase shifted baseband with hanning filter
bbshift_rsdata_odd = zeros(ntime, nbeam);   % r-sin(theta) for phase shifted baseband with odd elements
shifted_pulse = zeros(ntime, nelem);        % shifted pulse by sampling delay tm_prime

%
% Repeat for every beam
%
for ib=1:nbeam
	disp(sprintf('Beam %d of %d', ib, nbeam))

	% --> QUESTION c, d, g, h <--
	% For the current beam, compute time delays (or delayed samples)
	% and phase rotations (if needed) for each channel (i.e. transducer 
    % element) as a function of range, as well as gain corrections.
        for ie=1:nelem
            tn_prime = -((xn(ie)*sin_theta(ib))/c) + ((xn(ie)^2)*(1-(sin_theta(ib))^2)) ./ (2*c*r.');
            m = round(tn_prime/deltat);
            tt = 1:ntime;
            shift_idx = tt' - m;

            for i = 1:length(shift_idx)
                if (shift_idx(i) <= 1800) && (shift_idx(i) > 0)
                    shifted_pulse(shift_idx(i), ie) = Data(i, ie);
                    bbshift(shift_idx(i), ie) = databb(i, ie);
                end
            end
            phase_rot = exp(-1i*w0.*tn_prime);
            bbshift(:, ie) = bbshift(:, ie) .* phase_rot;

        end

	% --> QUESTION e, i, j, k <--
	% Create data for r-sin(theta) buffer by computing the coherent
	% sum across the array using the delays & phase rotations from above
    rsdata(:, ib) = mean(shifted_pulse, 2);
    bbshift_rsdata(:, ib) = mean(bbshift, 2);
    han_window = hann(nelem);   % getting hanning window for nelem
    han_filter = repmat(han_window.', ntime, 1);    % building hanning filter
    bbshift_rsdata_han(:, ib) = mean(bbshift .* han_filter, 2);
    bbshift_rsdata_odd(:, ib) = mean(bbshift(:, (1:2:nbeam)), 2);
end



% Display the contents of the r-sin(theta) buffer as a gray scale image
% over a logarithmic scale of 40dB.  (change 40 to 20 to plot on a 20 dB
% scale)
%showimage3(abs(rsdata), 4, 40, axisx, ,axisy)


f=figure;
showimage3(abs(rsdata), 1, 40, 1, r(2)-r(1));
ylabel('r in mm');
xlabel('number of beams');
title('part (e): r-sin(theta)');

f=figure;
showimage3(abs(bbshift_rsdata), 1, 40, 1, r(2)-r(1));
ylabel('r in mm');
xlabel('number of beams');
title('part (i): r-sin(theta) for phase-shifted baseband data');

f=figure;
showimage3(abs(bbshift_rsdata_han), 1, 40, 1, r(2)-r(1));
ylabel('r in mm');
xlabel('number of beams');
title('part (j): r-sin(theta) for phase-shifted baseband data with hanning filter');

f=figure;
showimage3(abs(bbshift_rsdata_odd), 1, 40, 1, r(2)-r(1));
ylabel('r in mm');
xlabel('number of beams');
title('part (k): r-sin(theta) for phase-shifted baseband data with odd elements');


% part d and h
d_buffer = zeros(ntime, nelem);
h_buffer = zeros(ntime, nelem);

% for d and h, ib would nbeam/2 (theta = 0)

ib = round(nbeam/2);

for ie=1:nelem
    tn_prime = -((xn(ie)*sin_theta(ib))/c) + ((xn(ie)^2)*(1-(sin_theta(ib))^2))./(2*c*r.');
    m = round(tn_prime/deltat);
    tt = 1:ntime;
    shift_idx = tt' - m;

    for i=1:length(shift_idx)
        if (shift_idx(i) <= 1800) && (shift_idx(i) > 0)
            d_buffer(shift_idx(i), ie) = Data(i, ie);
            h_buffer(shift_idx(i), ie) = databb(i, ie);
        end
    end
    phase_rot = exp(-1i*w0.*tn_prime);
    h_buffer(:, ie) = h_buffer(:, ie) .* phase_rot;
end

f = figure();
tc1 = tiledlayout(2, 1);
nexttile(tc1);
hold on;
title('part (d): before shift pulse amplitude vs r distance travelled in time t');
plot(r, Data(:, 25));
plot(r, Data(:, 49));
plot(r, Data(:, 73));
legend('element 25', 'element 49', 'element 73');
xlabel('r');
ylabel('pulse amplitude');

nexttile(tc1);
hold on;
title('part (d): after shift pulse amplitude vs r distance travelled in time t');
plot(r, d_buffer(:, 25));
plot(r, d_buffer(:, 49));
plot(r, d_buffer(:, 73));
legend('element 25', 'element 49', 'element 73');
xlabel('r');
ylabel('pulse amplitude');


f = figure();
tc1 = tiledlayout(2, 1);
nexttile(tc1);
hold on;
title('part (h): before shift baseband vs r distance travelled in time t');
plot(r, databb(:, 25));
plot(r, databb(:, 49));
plot(r, databb(:, 73));
legend('element 25', 'element 49', 'element 73');
xlabel('r');
ylabel('pulse amplitude');

nexttile(tc1);
hold on;
title('part (h): after shift baseband vs r distance travelled in time t');
plot(r, h_buffer(:, 25));
plot(r, h_buffer(:, 49));
plot(r, h_buffer(:, 73));
legend('element 25', 'element 49', 'element 73');
xlabel('r');
ylabel('pulse amplitude');


