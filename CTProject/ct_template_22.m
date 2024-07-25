%
% template for 2022 BME/EECS 516 tomography project
%
% parameters for the 3 disk phantom
% x,y center, radius, 'amplitude' (e.g. attenuation coefficient)
%
% replace all ?'s
%
plotdir = './generated_plots';
if ~exist(plotdir,'dir'), mkdir(plotdir); end
saveplot = 1;
circ = [0 0 80 0.01; 30 30 20 0.02; -55 0 10 0.03];
nobj = size(circ,1);
%
% image parameters
%
FOV = 200; % mm
nx = 192; ny = 192;
dx = FOV/nx;
%
% geometry parameters
%
nr = nx; dr = dx; % # of radial samples, ray spacing
nrho = nr;
delta_rho = 1/FOV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 1 - determine number of angular samples
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = dr*[-nr/2:nr/2-1]'; % vector of R values
%%%
delta_theta = dx/(FOV/2);
na = pi/delta_theta;
na = 4 * round(na/4); % # of angular samples (make mult of 4)
%%%
ang = [0:(na-1)]'/na * pi;% angular sample positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 2 - compute sinogram for disk phantom
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1~=exist('sinogram1')
    % making sinogram for thick rays by averaging 5 infinitely thin
    % rays into a single thick ray - helps with sampling issues
    rf = dr/5*[-5*nr/2:5*nr/2-1]';
    rrf = rf(:,ones(1,na));
    sinogram1f = zeros(5*nr, na);
    for ii=1:size(circ,1)
        cx = circ(ii,1);
        cy = circ(ii,2);
        rad = circ(ii,3);
        amp = circ(ii,4);
        tau = cx * cos(ang) + cy * sin(ang);
        tau = tau(:,ones(1,5*nr))';
        t = find( (rrf-tau).^2 <= rad.^2 );
        if ii > 1, amp = amp - circ(1,4); end % small disks embedded in larger (subtract larger)
    sinogram1f(t) = sinogram1f(t)+amp*2*sqrt(rad^2-(rrf(t)-tau(t)).^2);
    end
    sinogram1ff = conv2(sinogram1f,ones([5 1]),'same')/5;
    sinogram1 = sinogram1ff(1:5:end,:);
end
%
% Output Image of Singram
%
f1 = figure('DefaultAxesFontSize',16);
f1.WindowState = "maximized";
imagesc(r,ang,sinogram1'); colormap('gray');
title('Part 2: Sinogram of Disk Phantom');
xlabel('R in (mm)');
ylabel('\theta in (rad)');
if saveplot == 1
    % saveas(f1, './generated_plots/part2.jpeg');
    pause(1);
    print(f1, './generated_plots/part2', '-dpng', '-r300');
end

% rename for use below
sinogram = sinogram1;

ang = [0:(na-1)]'/na * pi;
disp(sprintf('number of angles = %g', na));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 3 - Implement plain backprojection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implement backprogjection
bpimage = zeros(ny, nx);
for each_ang=1:na
    bpimage = bpimage + imrot3(repmat(sinogram(:, each_ang), 1, nr), ...
        ang(each_ang), 'bilinear');
end

%
% Display Image
%
f2 = figure('DefaultAxesFontSize',16);
f2.WindowState = "maximized";
imagesc(r,r,bpimage'); colormap('gray'); axis xy;
axis('image');
title('Part 3: Simple Backprojection Image');
xlabel('x in (mm)');
ylabel('y in (mm)');
if saveplot == 1
    % saveas(f2, './generated_plots/part3.jpeg');
    pause(1);
    print(f2, './generated_plots/part3', '-dpng', '-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 4 - Filter Projections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total 2nrho datapoints to be conformable
% with zeropadded sinogram
rho_zp = (-(nrho):(nrho-1)) .* delta_rho;

% making conformable dims
rho_zp = rho_zp.' * ones(1, na);

% padding sinogram
sinogram_zp = zeros(2*nr, na);
sinogram_zp(((2*nr)/4): ((2*nr)*3/4)-1, :) = sinogram;

% taking 1D FT of zero
% padded sinogram
sinogram_ft = ft(sinogram_zp);
sinogramfilt_zp = real(ift(sinogram_ft .* abs(rho_zp)));

% removing zeropadding
sinogramfilt = sinogramfilt_zp(((2*nr)/4): ((2*nr)*3/4)-1, :);

% sinogramfilt2 = ? % for Part 9
%
% Plot Filtered Sinogram
%
f3 = figure('DefaultAxesFontSize',16);
f3.WindowState = "maximized";
plot(r, sinogram(:,1)./max(sinogram(:,1)), '-',...
 r, sinogramfilt(:,1)./max(sinogramfilt(:,1)),':', 'LineWidth',1.0);
xlabel('R in (mm)');
ylabel('Radial projection value (Normalized by max value)');
legend('sinogram', 'filtered sinogram');
title('Part 4: Projection vs Filtered Projection at \theta = 0');
if saveplot == 1
    % saveas(f3, './generated_plots/part4.jpeg');
    pause(1);
    print(f3, './generated_plots/part4', '-dpng', '-r300');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 5 - Backproject the filtered sinogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagefbp = zeros(ny, nx);
for each_ang=1:na
    imagefbp = imagefbp + imrot3(repmat(sinogramfilt(:, each_ang), 1, nr), ...
        ang(each_ang), 'bilinear');
end

%
% Display Reconstructed Image with Negatives Set to Zero
%
f4 = figure('DefaultAxesFontSize',16);
f4.WindowState = "maximized";
imagesc(r,r,max(imagefbp,0)'); colormap('gray'); axis xy;
axis('image');
title('Part 5: FBP Reconstruction');
xlabel('x in (mm)');
ylabel('y in (mm)');
xlim([-100 100]);
ylim([-100 100]);
% yticks(-100:20:100)
% xticks(-100:20:100)
if saveplot == 1
    % saveas(f4, './generated_plots/part5.jpeg');
    pause(1);
    print(f4, './generated_plots/part5', '-dpng', '-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 6 - Generate Fourier Interpolation Image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = (-(nrho/2):(nrho/2-1)).*delta_rho;
[ang_mesh, rho_mesh] = meshgrid(ang, rho);
[v, u] = meshgrid(rho);
rho_intrp = sqrt(u.^2 + v.^2);
ang_intrp = atan2(v,u);

neg_ang_intrp_idx = ang_intrp<0;
ang_intrp(neg_ang_intrp_idx) = ang_intrp(neg_ang_intrp_idx) + pi;
rho_intrp(neg_ang_intrp_idx) = -1*rho_intrp(neg_ang_intrp_idx);

interp_sg = interp2(ang_mesh, rho_mesh, ft(sinogram), ang_intrp, rho_intrp, 'bilinear');
interp_sg(isnan(interp_sg)) = 0;
imagefi = ift2(interp_sg);

%
% Display Reconstructed Image
%
f5 = figure('DefaultAxesFontSize',16);
f5.WindowState = "maximized";
imagesc(r,r,abs(interp_sg)'); colormap('jet');
axis('image');
title('Part 6: F(u, v)');
set(gca, 'YDir', 'normal');
xlabel('x in (mm)');
ylabel('y in (mm)');
if saveplot == 1
    % saveas(f5, './generated_plots/part6_ft.jpeg')
    pause(1);
    print(f5, './generated_plots/part6_ft', '-dpng', '-r300');
end

f6 = figure('DefaultAxesFontSize',16);
f6.WindowState = "maximized";
imagesc(r,r,abs(imagefi)'); colormap('gray'); axis xy;
axis('image');
title('Part 6: FI Reconstruction');
xlabel('x in (mm)');
ylabel('y in (mm)');
if saveplot == 1
    % saveas(f6, './generated_plots/part6_image.jpeg');
    pause(1);
    print(f6, './generated_plots/part6_image', '-dpng', '-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 7 - Plot profiles through reconstructed images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = -dx * ([-ny/2:ny/2-1]);
linefbp = abs(imagefbp(nr:-1:1,nr/2));
linefi = abs(imagefi(nr:-1:1,nr/2));

f7 = figure('DefaultAxesFontSize',16);
f7.WindowState = "maximized";
plot(y, linefbp/max(linefbp),'-', y, linefi/max(linefi),'-', ...
    'LineWidth', 1.0);
xlabel('radial position (mm)');
ylabel('projection value');
title('Part 7: Profile y=0 for FI and FBP reconstruction method');
legend('FBP','FI');

if saveplot == 1
    % saveas(f7, './generated_plots/part7.jpeg');
    pause(1);
    print(f7, './generated_plots/part7', '-dpng', '-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 8 - 12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8
% subsampling sinogram
sinogram_subsample = sinogram(:,1:4:end);

rho_zp = (-(nrho):(nrho-1)).*delta_rho;
rho_zp = rho_zp.' .* ones(1,na/4);

% zero padding sinogram
sinogram_subsample_zp = zeros(2*nr, na/4);
sinogram_subsample_zp(((2*nr)/4): ((2*nr)*3/4)-1, :) = sinogram_subsample;

% taking 1D FT of zero
% padded sinogram
sinogram_subsample_ft = ft(sinogram_subsample_zp);
sinogramfilt_subsample_zp = real(ift(sinogram_subsample_ft .* abs(rho_zp)));
sinogramfilt1 = sinogramfilt_subsample_zp(((2*nr)/4): ((2*nr)*3/4)-1, :);

f8 = figure('DefaultAxesFontSize',16);
f8.WindowState = "maximized";
plot(r, sinogram_subsample(:,1)./max(sinogram_subsample(:,1)), '-',...
 r, sinogramfilt1(:,1)./max(sinogramfilt1(:,1)),':', 'LineWidth',1.0);
xlabel('R in (mm)');
ylabel('\theta in (rad)');

legend('subsampled sinogram', 'subsampled filtered sinogram');
title('subsampled sinogram projection vs subsampled sinogram filtered projection at \theta = 0');
imagefbp1 = zeros(ny, nx);

for each_ang=1:na/4
    imagefbp1 = imagefbp1 + imrot3(repmat(sinogramfilt1(:, each_ang), 1, nr), ...
        ang((each_ang-1)*4+1), 'bilinear');
end

f8 = figure('DefaultAxesFontSize',16);
f8.WindowState = "maximized";
imagesc(r,r,max(imagefbp1,0)'); colormap('gray'); axis xy;
axis('image');
title('Part 8: FBP Reconstruction from subsampled sinogram');
xlabel('x in (mm)');
ylabel('y in (mm)');
xlim([-100 100]);
ylim([-100 100]);

if saveplot == 1
    % saveas(f8, './generated_plots/part8.jpeg');
    pause(1);
    print(f8, './generated_plots/part8', '-dpng', '-r300');
end

%% 9
% subsampling sinogram
sinogram_subsample = sinogram(:,1:4:end);

rho_zp = (-(nrho):(nrho-1)).*delta_rho;
rho_zp = rho_zp.' .* ones(1,na/4);

% zero padding sinogram
sinogram_subsample_zp = zeros(2*nr, na/4);
sinogram_subsample_zp(((2*nr)/4): ((2*nr)*3/4)-1, :) = sinogram_subsample;

sinogram_subsample_ft = ft(sinogram_subsample_zp);

hanfilt = hanning(nr);
hanfilt = [zeros(1,nr/2), hanfilt.', zeros(1,nr/2)];
hanfilt = hanfilt.' .* ones(1,na/4);

sinogramfilt_subsample_zp = real(ift(sinogram_subsample_ft .* ...
 abs(rho_zp) .* hanfilt));
sinogramfilt2 = sinogramfilt_subsample_zp(((2*nr)/4): ((2*nr)*3/4)-1, :);

% figure;
% surf(abs(rho_padded) .* hanfilt);
figure(7)
plot(r, sinogram_subsample(:,1)./max(sinogram_subsample(:,1)), '-',...
 r, sinogramfilt2(:,1)./max(sinogramfilt2(:,1)),':', 'LineWidth', 1.0);
xlabel('R in (mm)');
ylabel('\theta in (rad)');
legend('subsampled sinogram', 'subsampled filtered sinogram');
title('subsampled sinogram projection vs subsampled sinogram filtered projection at \theta = 0');

imagefbp2 = zeros(ny, nx);
for each_ang=1:na/4
    imagefbp2 = imagefbp2 + imrot3(repmat(sinogramfilt2(:, each_ang), 1, nr), ...
    ang((each_ang-1)*4+1), 'bilinear');
end

f9 = figure('DefaultAxesFontSize',16);
f9.WindowState = "maximized";
imagesc(r,r,max(imagefbp2,0)'); colormap('gray'); axis xy;
axis('image');
title('Part 9: FBP Reconstruction from subsampled sinogram with hanning filter');
xlabel('x in (mm)');
ylabel('y in (mm)');
xlim([-100 100]);
ylim([-100 100]);

if saveplot == 1
    % saveas(f9, './generated_plots/part9.jpeg');
    pause(1);
    print(f9, './generated_plots/part9', '-dpng', '-r300');
end

%% 10
n0 = 10000;
proj1 = n0*exp(-sinogram1); % mean number of counts

% noisy sinogram
sinogram = log(n0./poissrnd(proj1));

% total 2nrho datapoints to be conformable
% with zeropadded sinogram
rho_zp = (-(nrho):(nrho-1)) .* delta_rho;

% making conformable dims
rho_zp = rho_zp.' * ones(1, na);

% padding sinogram
sinogram_noisy_zp = zeros(2*nr, na);
sinogram_noisy_zp(((2*nr)/4): ((2*nr)*3/4)-1, :) = sinogram;

% taking 1D FT of zero
% padded sinogram
sinogram_noisy_ft = ft(sinogram_noisy_zp);
sinogramfilt_noisy_zp = real(ift(sinogram_noisy_ft .* abs(rho_zp)));

% removing zeropadding
sinogramfilt3 = sinogramfilt_noisy_zp(((2*nr)/4): ((2*nr)*3/4)-1, :);

%
% Plot Filtered Sinogram
%
f10 = figure('DefaultAxesFontSize',16);
f10.WindowState = "maximized";
plot(r, sinogram(:,1)./max(sinogram(:,1)), '-',...
 r, sinogramfilt3(:,1)./max(sinogramfilt3(:,1)),':', 'LineWidth',1.0);
xlabel('R in (mm)');
ylabel('\theta in (rad)');
legend('sinogram', 'filtered sinogram');
title('Part 10: Noisy Projection vs Noisy Filtered Projection at \theta = 0; n_0 = 10000');

if saveplot == 1
    % saveas(f10, './generated_plots/part10_projection.jpeg');
    pause(1);
    print(f10, './generated_plots/part10_projection', '-dpng', '-r300');
end

% Reconstructing image
imagefbp3 = zeros(ny, nx);
for each_ang=1:na
    imagefbp3 = imagefbp3 + imrot3(repmat(sinogramfilt3(:, each_ang), 1, nr), ...
        ang(each_ang), 'bilinear');
end

% Display Reconstructed Image for Noisy Sinogram
f11 = figure('DefaultAxesFontSize',16);
f11 .WindowState = "maximized";
imagesc(r,r,max(imagefbp3,0)'); colormap('gray'); axis xy;
axis('image');
title('Part 10: FBP Reconstruction for Noisy Sinogram; n_0 = 10000');
xlabel('x in (mm)');
ylabel('y in (mm)');
xlim([-100 100]);
ylim([-100 100]);
% yticks(-100:20:100)
% xticks(-100:20:100)

if saveplot == 1
    % saveas(f11, './generated_plots/part10_recon_image.jpeg');
    pause(1);
    print(f11, './generated_plots/part10_recon_image', '-dpng', '-r300');
end

%% 11
n0 = 1000;
proj1 = n0*exp(-sinogram1); % mean number of counts

% noisy sinogram
sinogram = log(n0./poissrnd(proj1));

% total 2nrho datapoints to be conformable
% with zeropadded sinogram
rho_zp = (-(nrho):(nrho-1)) .* delta_rho;

% making conformable dims
rho_zp = rho_zp.' * ones(1, na);

% padding sinogram
sinogram_noisy_zp = zeros(2*nr, na);
sinogram_noisy_zp(((2*nr)/4): ((2*nr)*3/4)-1, :) = sinogram;

% taking 1D FT of zero
% padded sinogram
sinogram_noisy_ft = ft(sinogram_noisy_zp);
sinogramfilt_noisy_zp = real(ift(sinogram_noisy_ft .* abs(rho_zp)));

% removing zeropadding
sinogramfilt4 = sinogramfilt_noisy_zp(((2*nr)/4): ((2*nr)*3/4)-1, :);

%
% Plot Filtered Sinogram
%
f12 = figure('DefaultAxesFontSize',16);
f12.WindowState = "maximized";
plot(r, sinogram(:,1)./max(sinogram(:,1)), '-',...
 r, sinogramfilt4(:,1)./max(sinogramfilt4(:,1)),':', 'LineWidth',1.0);
xlabel('R in (mm)');
ylabel('\theta in (rad)');
legend('sinogram', 'filtered sinogram');
title('Part 11: Noisy Projection vs Noisy Filtered Projection at \theta = 0; n_0 = 1000');

if saveplot == 1
    % saveas(f12, './generated_plots/part11_projection.jpeg');
    pause(1);
    print(f12, './generated_plots/part11_projection', '-dpng', '-r300');
end

% Reconstructing image
imagefbp4 = zeros(ny, nx);
for each_ang=1:na
    imagefbp4 = imagefbp4 + imrot3(repmat(sinogramfilt4(:, each_ang), 1, nr), ...
        ang(each_ang), 'bilinear');
end

% Display Reconstructed Image for Noisy Sinogram
f12_1 = figure('DefaultAxesFontSize',16);
f12_1.WindowState = "maximized";
imagesc(r,r,max(imagefbp4,0)'); colormap('gray'); axis xy;
axis('image');
title('Part 11: FBP Reconstruction for Noisy Sinogram; n_0 = 1000');
xlabel('x in (mm)');
ylabel('y in (mm)');
xlim([-100 100]);
ylim([-100 100]);
% yticks(-100:20:100)
% xticks(-100:20:100)

if saveplot == 1
    % saveas(f12_1, './generated_plots/part11_recon_image.jpeg');
    pause(1);
    print(f12_1, './generated_plots/part11_recon_image', '-dpng', '-r300');
end

%% 12
n0 = 1000;
proj1 = n0*exp(-sinogram1); % mean number of counts

% noisy sinogram
sinogram = log(n0./poissrnd(proj1));

% total 2nrho datapoints to be conformable
% with zeropadded sinogram
rho_zp = (-(nrho):(nrho-1)) .* delta_rho;

% making conformable dims
rho_zp = rho_zp.' * ones(1, na);

hanfilt = hanning(nr);
hanfilt = [zeros(1,nr/2), hanfilt.', zeros(1,nr/2)];
hanfilt = hanfilt.' .* ones(1,na);

% padding sinogram
sinogram_noisy_zp = zeros(2*nr, na);
sinogram_noisy_zp(((2*nr)/4): ((2*nr)*3/4)-1, :) = sinogram;

% taking 1D FT of zero
% padded sinogram
sinogram_noisy_ft = ft(sinogram_noisy_zp);
sinogramfilt_noisy_zp = real(ift(sinogram_noisy_ft .* abs(rho_zp) .*hanfilt));

% removing zeropadding
sinogramfilt5 = sinogramfilt_noisy_zp(((2*nr)/4): ((2*nr)*3/4)-1, :);

%
% Plot Filtered Sinogram
%
f13 = figure('DefaultAxesFontSize',16);
f13.WindowState = "maximized";
plot(r, sinogram(:,1)./max(sinogram(:,1)), '-',...
 r, sinogramfilt5(:,1)./max(sinogramfilt5(:,1)),':', 'LineWidth',1.0);
xlabel('R in (mm)');
ylabel('\theta in (rad)');
legend('sinogram', 'filtered sinogram');
title(['Part 12: Noisy Projection vs Noisy Filtered Projection at \theta = 0;' ...
    ' n_0 = 1000 with hanning filtering']);

if saveplot == 1
    % saveas(f13, './generated_plots/part12_projection.jpeg');
    pause(1);
    print(f13, './generated_plots/part12_projection', '-dpng', '-r300');
end

% Reconstructing image
imagefbp5 = zeros(ny, nx);
for each_ang=1:na
    imagefbp5 = imagefbp5 + imrot3(repmat(sinogramfilt5(:, each_ang), 1, nr), ...
    ang(each_ang), 'bilinear');
end

% Display Reconstructed Image for Noisy Sinogram
f14 = figure('DefaultAxesFontSize',16);
f14 .WindowState = "maximized"; 
imagesc(r,r,max(imagefbp5,0)'); colormap('gray'); axis xy;
axis('image');
title('Part 12: FBP Reconstruction for Noisy Sinogram; n_0 = 1000 with hanning filter');
xlabel('x in (mm)');
ylabel('y in (mm)');
xlim([-100 100]);
ylim([-100 100]);
% yticks(-100:20:100)
% xticks(-100:20:100)

if saveplot == 1
    % saveas(f14, './generated_plots/part12_recon_image.jpeg');
    pause(1);
    print(f14, './generated_plots/part12_recon_image', '-dpng', '-r300');
end
%% 13

load mys22;

% different number of angles
[nr, na] = size(sinogram);

ang = [0:(na-1)]'/na * pi;
disp(sprintf('number of angles = %g', na));

f15 = figure('DefaultAxesFontSize',16);
f15.WindowState = "maximized";
imagesc(r,ang,sinogram'); colormap('gray');
title('Part 13: Mystery Sinogram (mys22.mat)');
xlabel('R in (mm)');
ylabel('\theta in (rad)');

if saveplot == 1
    % saveas(f15, './generated_plots/part13_sinogram.jpeg');
    pause(1);
    print(f15, './generated_plots/part13_sinogram', '-dpng', '-r300');

end

% total 2nrho datapoints to be conformable
% with zeropadded sinogram
rho_zp = (-(nrho):(nrho-1)) .* delta_rho;

% making conformable dims
rho_zp = rho_zp.' * ones(1, na);

% padding sinogram
mys_sinogram_zp = zeros(2*nr, na);
mys_sinogram_zp(((2*nr)/4): ((2*nr)*3/4)-1, :) = sinogram;

% taking 1D FT of zero
% padded sinogram
mys_sinogram_ft = ft(mys_sinogram_zp);
mys_sinogramfilt_zp = real(ift(mys_sinogram_ft .* abs(rho_zp)));

% removing zeropadding
mys_sinogramfilt = mys_sinogramfilt_zp(((2*nr)/4): ((2*nr)*3/4)-1, :);

%
% Plot Filtered Sinogram
%
f16 = figure('DefaultAxesFontSize',16);
f16.WindowState = "maximized";
plot(r, sinogram(:,1)./max(sinogram(:,1)), '-',...
 r, mys_sinogramfilt(:,1)./max(mys_sinogramfilt(:,1)),':', 'LineWidth',1.0);
xlabel('R in (mm)');
ylabel('\theta in (rad)');
legend('sinogram', 'filtered sinogram');
title('Part 13: Projection vs Filtered Projection at \theta = 0 for Mystery object');

if saveplot == 1
    % saveas(f16, './generated_plots/part13_proj.jpeg');
    pause(1);
    print(f16, './generated_plots/part13_proj', '-dpng', '-r300');
end

% Reconstructing image for Mystery object
imagefbp6 = zeros(ny, nx);
for each_ang=1:na
    imagefbp6 = imagefbp6 + imrot3(repmat(mys_sinogramfilt(:, each_ang), 1, nr), ...
    ang(each_ang), 'bilinear');
end

%
% Display Reconstructed Image with Negatives Set to Zero
%
f17 = figure('DefaultAxesFontSize',16);
f17.WindowState = "maximized";
imagesc(r,r,max(imagefbp6,0)'); colormap('gray'); axis xy;
axis('image');
title('Part 13: FBP Reconstruction for Mystery Object');
xlabel('x in (mm)');
ylabel('y in (mm)');
xlim([-100 100]);
ylim([-100 100]);
% yticks(-100:20:100)
% xticks(-100:20:100)
if saveplot == 1
    % saveas(f17, './generated_plots/part13_recon_image.jpeg');
    pause(1);
    print(f17, './generated_plots/part13_recon_image', '-dpng', '-r300');
end