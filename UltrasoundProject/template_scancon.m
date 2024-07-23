%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Daksh Maheshwari %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% EECS 516 %%%%%%%%%%%%%%%%%%
%%%%%%% Ultrasound beamforming project %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%	template_scancon.m
%	template script for converting from r-sin(theta) data to x-y image
%	must load r-sin(theta) data: rsdata

% --> QUESTION l. <--
% Scan convert the r-sin(theta) buffer to produce a sector scan image.
% Use bilinear interpolation to compute the image values on the
% sector scan image grid.  Matlab's "interp2" function will help you
% do bilinear interpolation.

% compute values needed for interp2
nb_pixels = 512;
delx = 60/nb_pixels;
delz = 60/nb_pixels;
rmax = 69;

xz_range = linspace(-round(rmax/2), round(rmax/2), nb_pixels);
z_range = linspace(0, rmax, nb_pixels);
[XZ, Z] = meshgrid(xz_range, z_range);    % mesh in xz and x space

[ST, R] = meshgrid(sin_theta, r);

R_coords = sqrt(Z.^2 + XZ.^2);
ST_coords = XZ ./ R_coords;

% Create image w/ bilinear interpolation

%%%%%%%%%%%%%%%%%%% rsdata %%%%%%%%%%%%%%%%%%% 

im_rsdata = interp2(ST, R, abs(rsdata), ST_coords, R_coords, 'bilinear');
t = find(isnan(im_rsdata));
im_rsdata(t) = zeros(size(t));
% im(t) = zeros(size(t));

figure;
showimage3(im_rsdata, 1, 20, delx, delz); axis('image');    % Display 20dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (20 dB)');

figure;
showimage3(im_rsdata, 1, 40, delx, delz); axis('image');    % Display 40dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (40 dB)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% bbshift_rsdata %%%%%%%%%%%%%%%

im_bbshift_rsdata = interp2(ST, R, abs(bbshift_rsdata), ST_coords, R_coords, 'bilinear');
t = find(isnan(im_bbshift_rsdata));
im_bbshift_rsdata(t) = zeros(size(t));
% im(t) = zeros(size(t));

figure;
showimage3(im_bbshift_rsdata, 1, 20, delx, delz); axis('image');    % Display 20dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (20 dB) for phase-shifted baseband data (20 dB)');

figure;
showimage3(im_bbshift_rsdata, 1, 40, delx, delz); axis('image');    % Display 40dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (40 dB) for phase-shifted baseband data (40 dB)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% bbshift_rsdata_han %%%%%%%%%%%%%

im_bbshift_rsdata_han = interp2(ST, R, abs(bbshift_rsdata_han), ST_coords, R_coords, 'bilinear');
t = find(isnan(im_bbshift_rsdata));
im_bbshift_rsdata(t) = zeros(size(t));
% im(t) = zeros(size(t));

figure;
showimage3(im_bbshift_rsdata_han, 1, 20, delx, delz); axis('image');    % Display 20dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (20 dB) for phase-shifted baseband data with hanning filter (20 dB)');

figure;
showimage3(im_bbshift_rsdata_han, 1, 40, delx, delz); axis('image');    % Display 40dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (40 dB) for phase-shifted baseband data with hanning filter (40 dB)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% bbshift_rsdata_odd %%%%%%%%%%%%%

im_bbshift_rsdata_odd = interp2(ST, R, abs(bbshift_rsdata_odd), ST_coords, R_coords, 'bilinear');
t = find(isnan(im_bbshift_rsdata));
im_bbshift_rsdata(t) = zeros(size(t));
% im(t) = zeros(size(t));

figure;
showimage3(im_bbshift_rsdata_odd, 1, 20, delx, delz); axis('image');    % Display 20dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (20 dB) for phase-shifted baseband data with odd elements (20 dB)');

figure;
showimage3(im_bbshift_rsdata_odd, 1, 40, delx, delz); axis('image');    % Display 40dB scale image
xlabel('xz in mm');
ylabel('z in mm');
title('Scan of r-sin(theta) (40 dB) for phase-shifted baseband data with odd elements (40 dB)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do similarly all r-sin(theta) buffers

% % --> QUESTION m. <--
% % Use two images on a logarithmic scale to answer this question:
% % one on a 40dB scale, the other on a 20dB scale
% figure; showimage3(im, 1, 20, axisx, axisy); axis('image')		% Display 20 dB scale image
% figure; showimage3(im, 1, 40, axisx, axisy); axis('image')		% Display 40 dB scale image
% 
% % do similarly all r-sin(theta) buffers

