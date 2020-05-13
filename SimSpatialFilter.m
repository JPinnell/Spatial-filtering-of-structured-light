%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate spatial filtering of structured light
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J.Pinnell (2020)
% Based on Pinnell, J., Klug, A. and Forbes, A., "Spatial filtering of structured light" (2020)
% Simply run code, comment out sections and play around as desired
% Requires LG.m, HG.m, BG.m since MatLab in-built polynomials take long to evaluate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make coordinates (in mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = 1000; 
dx = 8e-3; 
x = dx.*(-H/2:H/2-1); 
[X,Y] = meshgrid(x,-x);
[Phi,R] = cart2pol(X,Y); % polar coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make structured light beam (the signal beam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = 1; % Gaussian waist radius (required for all beams)
% Gaussian beam
% Beam = exp(-(R.^2)./w0.^2);
% Laguerre-Gaussian beam
L = [3,-3]; P = [3,3]; weights = [1,1];
Beam = LG(R,Phi,P,L,weights,w0);
% Hermite-Gaussian beam
% N = [4]; M = [4]; weights = [1];
% Beam = HG(X,Y,N,M,weights,w0);
% Bessel-Gaussian beam
% L = 0; Kr = 10;
% Beam = BG(R,Phi,Kr,L,weights,w0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High frequency noise
a = 0.1; % amplitude of noise
k_n = 2.*pi.*[2,2]; % spatial frequencies of noise
Noise = a.*(sin(k_n(1).*X + k_n(2).*Y));
% Low frequency (Gaussian noise)
% a = 0.1; % amplitude of noise
% mu = 0; % mean of Gaussian noise
% sigma = 5; % standard deviation of Gaussian noise
% Noise =  a.*imnoise(ones(H),'gaussian',mu,sigma);
Beam_Noisy = Beam + Noise; 
% Note: in reality the noise is limited to the beam only 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make 2D optical DFT matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = 1000; % focal length of lens
lambda = 633e-6; % wavelength of light
Lk = lambda*f/dx; % side "length" at Fourier plane (inverse length)
dk = lambda*f/(H*dx); % sample "size" at Fourier plane (inverse length)
k = -Lk/2:dk:Lk/2-dk; % spatial frequency coordinate system
k = k./dk^2; % rescale
FTM = exp(-1i*2*pi/H).^(k'*x); % discrete optical Fourier transform matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 2; % mask width parameter
FourierBeam = FTM'*Beam*FTM; % performing 2D discrete Fourier transform
M = max(max(abs(FourierBeam))); % maximum amplitude of signal beam at Fourier plane
Mask = abs(FourierBeam) > M*exp(-t); % Eq. 19 in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spatialy filter noisy beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FourierBeam_Noisy = FTM'*Beam_Noisy*FTM; % Fourier transform noisy beam
Beam_Filtered = FourierBeam_Noisy.*Mask; % apply mask to noisy beam
Beam_Filtered = FTM'*Beam_Filtered*FTM; % Fourier transform back to original size

% plots
figure('color','w','units','points','position',[50 50 600 400]);
sx = 0.01; wx = 0.32; sy = 0.03; wy = 0.45; % positioning 
subplot('Position',[sx 2*sy+wy wx wy]); imagesc(abs(Beam).^2); axis image off; title('Ideal beam');
subplot('Position',[2*sx+wx 2*sy+wy wx wy]); imagesc(abs(Beam_Noisy).^2); axis image off; title('Noisy beam');
subplot('Position',[3*sx+2*wx 2*sy+wy wx wy]); imagesc(abs(FourierBeam_Noisy).^2); axis image off; title('Noisy beam at Fourier plane');
subplot('Position',[sx sx wx wy]); imagesc(Mask); axis image off; title('Filter mask');
subplot('Position',[2*sx+wx sx wx wy]); imagesc(abs(Beam_Filtered).^2); axis image off; title('Filtered beam');

% Correlation coefficient between beam amplitudes before and after filtering
fprintf('Initial correlation: %2.1f%% \n',corr2(abs(Beam),abs(Beam_Noisy))*100);
fprintf('Final correlation: %2.1f%% \n',corr2(abs(Beam),abs(Beam_Filtered))*100);