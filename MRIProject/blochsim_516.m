function [mx,my,mz] = blochsim_516(m0, bx, by, bz, T1, T2, dt)

% blochsim_516 - Bloch simulator for 516 MRI project
% Simulations rotations and relaxations for a given timecourse.

% Syntax:  [mx,my,mz] = function_name(m0, bx, by, bz, T1, T2, dt)
% Inputs:
%    m0        - [3,dim]     | Initial x,y,z magnitization
%    bx,by,bz  - [ntime,dim] | Effective X,Y,Z applied magnetic field (Tesla), 
%                              for rotating frame (no B0)
%    T1        - [dim]       | Spin-lattice relaxation time (ms)
%    T2        - [dim]       | spin-spin relaxation time (ms)
%	 dt	       - scalar		 | time interval (ms)
%    
%    
% Outputs:
%    mx,my,mz [ntime,dim]    | Simulated magnitization M(t) in x,y,z space

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none


% Oct 2018; Last revision: Oct-4-2018

% constants
gambar = 42.57e3;        % gamma/2pi in kHz/T
gam = gambar*2*pi;
% Put Beff into units rotations, T1 and T2 into losses
bx = bx.*(dt*gam);      % rotation angle/step
by = by.*(dt*gam);      % rotation angle/step
bz = bz.*(dt*gam);      % rotation angle/step

% Put relaxations into losses/recovery per step
T1 = (dt ./ T1(:)); 
T2 = (1 - dt ./ T2(:));

% size checks
if ~isreal(bx) | ~isreal(by) | ~isreal(bz) 
    bx = real(bx);
    by = real(by);
    bz = real(bz);
    warn('Warning: B field must be real valued - using only the real part');
end
if (size(bx) ~= size(by)) | (size(bx) ~= size(bz))
    error('Error: B vectors not the same length')
end
if (size(m0,2) ~= size(bx,2))
    error('Error: Initial magnetization not right size')
end
if (size(m0,2) ~= length(T1))
    error('Error: T1 vector not right size')
end
if (size(m0,2) ~= length(T2))
    error('Error: T2 vector not right size')
end

nstep = size(bx,1);
%
% Initialize outputs
mx = zeros(size(bx));
my = zeros(size(bx));
mz = zeros(size(bx));
mx(1,:) = m0(1,:);
my(1,:) = m0(2,:);
mz(1,:) = m0(3,:);

% stable bloch equation simulator: rotations are explicitly
% calculated and carried out on the magnetization vector 
for lp = 2:nstep
  B = [bx(lp-1,:); by(lp-1,:); bz(lp-1,:)]';  
  %	Compute sines & cosines of field angles:
  %	Theta = angle w.r.t positive z axis
  %	Phi   = angle w.r.t positive x axis
  %	Psi   = angle w.r.t transformed positive x axis
  %
  Bmag = sqrt(sum(B.^2,2));		% Magnitude of applied field
  Btrans = sqrt(B(:,1).^2 + B(:,2).^2);	% Magnitude of transverse applied field
  ct = ones(size(B,1),1);
  good = Bmag ~= 0;
  if any(good)
	ct(good) = B(good,3) ./ Bmag(good);	% cos(theta)
  end
  st = sqrt(1 - ct.^2);				% sin(theta) > 0

  cphi = ones(size(B,1),1);
  good = Btrans ~= 0;
  if any(good)
	cphi(good) = B(good,1) ./ Btrans(good);	% cos(phi)
  end
  sphi = sqrt(1 - cphi.^2) .* sign(B(:,2));	% sin(phi)

  cpsi = cos(Bmag);			% cos(psi)
  spsi = sin(Bmag);			% sin(psi)

  %
  %	Evolve
  %
  if any(Bmag ~= 0)
    Mx0 = mx(lp-1,:)';
    My0 = my(lp-1,:)';
    Mz0 = mz(lp-1,:)';
    
    Mx1 = cphi.*(ct.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
    + spsi.*(cphi.*My0-sphi.*Mx0))+st.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0))) ...
    - sphi.*(-spsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
    + cpsi.*(cphi.*My0-sphi.*Mx0));
    My1 = sphi.*(ct.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
    + spsi.*(cphi.*My0-sphi.*Mx0))+st.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0))) ...
    + cphi.*(-spsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
    + cpsi.*(cphi.*My0-sphi.*Mx0));
    Mz1 = ct.*(ct.*Mz0+st.*(sphi.*My0+cphi.*Mx0)) ...
    - st.*(cpsi.*(ct.*(sphi.*My0+cphi.*Mx0)-st.*Mz0) ...
    + spsi.*(cphi.*My0-sphi.*Mx0));
  else
    Mx1 = mx(lp-1,:)';
    My1 = my(lp-1,:)';
    Mz1 = mz(lp-1,:)';
  end
  % relaxation effects: "1" in Mz since Mo=1 by assumption
  mx(lp,:) = (Mx1 .* T2)';
  my(lp,:) = (My1 .* T2)';
  mz(lp,:) = (Mz1 + (1- Mz1).* T1)';

end % end loop through time