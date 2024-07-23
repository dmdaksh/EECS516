%
%	template for 2022 BME/EECS 516 tomography project
%
%	parameters for the 3 disk phantom
%	x,y center, radius, 'amplitude' (e.g. attenuation coefficient)
%
%       replace all ?'s
%

circ = [0 0 80 0.01; 30 30 20 0.02; -55 0 10 0.03];
nobj = size(circ,1);

%
%	image parameters
%
FOV = 200; % mm
nx = 192; ny = 192;
dx = FOV/nx;		                 

%
%	geometry parameters
%

nr = nx;	dr = dx;		% # of radial samples, ray spacing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 1 - determine number of angular samples
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = dr*[-nr/2:nr/2-1]'; % vector of R values
na = ?		% # of angular samples (make mult of 4)
ang = [0:(na-1)]'/na * pi;	% angular sample positions

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
    if ii > 1, amp = amp - circ(1,4); end	% small disks embedded in larger (subtract larger)
    sinogram1f(t) = sinogram1f(t)+amp*2*sqrt(rad^2-(rrf(t)-tau(t)).^2);
  end
  sinogram1ff = conv2(sinogram1f,ones([5 1]),'same')/5;
  sinogram1 = sinogram1ff(1:5:end,:);
end

%
% Output Image of Singram 
%

  figure(1)
  imagesc(r,ang,sinogram1'); colormap('gray')
  title('Sinogram of Disk Phantom')
  xlabel('???')
  ylabel('???')

  % rename for use below
  sinogram = sinogram1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 10 - compute noisy sinogram for disk phantom
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n0 = 10000;
proj1 = n0*exp(-sinogram1); % mean number of counts
% sinogram = ?;  % for noisy data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 13 - Load mystery object singogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load mys22;
% different number of angles
[nr, na] = size(sinogram);



ang = [0:(na-1)]'/na * pi;
disp(sprintf('number of angles = %g', na))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 3 - Implement plain backprojection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% implement backprogjection
% bpimage = ?;

%
% Display Image
%

  figure(2)
  imagesc(r,r,bpimage'); colormap('gray')
  axis('image')
  title('Simple Backprojection Image')
  xlabel('???')
  ylabel('???')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 4 - Filter Projections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Filter sinogram
%

% sinogramfilt = ?
 

% sinogramfilt2 = ?  % for Part 9

%
% Plot Filtered Sinogram
%

  figure(3)
  plot(r, sinogram(:,1)./max(sinogram(:,1)), '-',...
       r, sinogramfilt(:,1)./max(sinogramfilt(:,1)),':');
  xlabel('???');
  ylabel('???');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 5 - Backproject the filtered sinogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  imagefbp = ?;

%
% Display Reconstructed Image with Negatives Set to Zero
%

  figure(4)
  imagesc(r,r,max(imagefbp,0)'); colormap('gray')
  axis('image')
  title('FBP Reconstruction')
  xlabel('???')
  ylabel('???')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 6 - Generate Fourier Interpolation Image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% imagefi = ?;

%
% Display Reconstructed Image 
%

  figure(5)
  imagesc(r,r,abs(imagefi)'); colormap('gray')
  axis('image')
  title('FI Reconstruction')
  xlabel('???')
  ylabel('???')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 7 - Plot profiles through reconstructed images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  y = -dx * ([-ny/2:ny/2-1]);
  linefbp = abs(imagefbp(nr:-1:1,nr/2));
  linefi = abs(imagefi(nr:-1:1,nr/2));
  figure(6)
  plot(y, linefbp/max(linefbp),'-', y, linefi/max(linefi),'-');
  xlabel('radial position (mm)');
  ylabel('projection value');
  legend('FBP','FI');

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 8 - 12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



