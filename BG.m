function U = BG(R,Phi,Kr,L,weights,w0)
% This function computes a superposition of BG mode at the plane z=0
% R,Phi are coordinate matrices
% "weights" is a weight vector for the coefficients in the superposition
% E.g. generate mode U = BG_10^1 + BG_10^-1 then,
% x = -1:0.01:1, y = -1:0.01:1, [X,Y] = meshgrid(x,y); [Phi,R] = cart2pol(X,Y);
% Kr = [10,10], L = [1,-1], weights = [1,1]; w0 = 1;
% U = BG(R,Phi,Kr,L,weights,w0);

U = zeros(size(R)); %initialise electric field

for i = 1:length(weights)
    U = U + weights(i).*besselj(L(i),Kr(i).*R).*exp(-(R.^2)./(w0.^2)).*exp(1i.*L(i).*Phi);
end

end