function [Mask] = GenFilterMask(Beam,X,Y,lambda,f,t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Mask] = GenFilterMask(Beam,lambda,f)
% v1 J.Pinnell 2020
% Based on Pinnell, J., Klug, A. and Forbes, A., "Spatial filtering of structured light" (2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the binary mask required to spatially filter 
% structured light. This is achieved by encoding the mask on a suitable 
% optical device (such as a spatial light modulator) at the Fourier plane. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X,Y    - x and y meshgrid coordinate system (can be SLM coordinates)
% Beam   - input structured light beam in the above coordinate system
% lambda - wavelength of the light
% f      - focal length of the Fourier lens
% t      - (optional) mask width parameter. Default value is 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mask - binary amplitude mask 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The coordinates must be a square grid with square pixels since it is 
% required for the optical DFT matrix. If required, it is straightforward 
% to adjust the mask afterwards by padding and using imresize.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Minimum working example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = 1000; dx = 8e-6; x = dx.*(-H/2:H/2-1); [X,Y] = meshgrid(x,-x); %coord
% w0 = 1e-3; Beam = exp(-(X.^2+Y.^2)./w0.^2); %Gaussian beam with 1mm width
% f = 1; lambda = 633e-9; %Parameters of the optical system
% Mask = GenFilterMask(Beam,X,Y,lambda,f,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preliminary checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,H] = size(X); % extract size
dx = abs(X(1,2) - X(1,1));
dy = abs(Y(2,1) - Y(1,1));
if V ~= H || dx ~= dy
    error('Coordinates must be square matrices with square pixels');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make 2D optical DFT matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = H*dx; % physical side length of image
x = X(1,:); % x coordinates
Lk = lambda*f/dx; % side "length" at Fourier plane (inverse length)
dk = lambda*f/L; % sample "size" at Fourier plane (inverse length)
k = -Lk/2:dk:Lk/2-dk; % spatial frequency coordinate system
k = k./dk^2; % rescale
FTM = exp(-1i*2*pi/H).^(k'*x); % discrete optical Fourier transform matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make mask by Fourier transforming Beam and thresholding amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FourierBeam = FTM'*Beam*FTM; % performing 2D discrete Fourier transform
M = max(max(abs(FourierBeam))); %maximum amplitude of beam at Fourier plane
if nargin == 4
    t = 2;
end
Mask = abs(FourierBeam) > M*exp(-t); % Eq. 19 in paper

end