
function S = L1Smooth(Im, alpha,scale,lambda, kappa)
%   Example
%   ==========
%   Im  = imread('pflower.jpg');
%   S  = L1Smooth(Im); % Default Parameters (alpha=0.02,scale=1.5,lambda = 2e-2, kappa = 2)
%   figure, imshow(Im), figure, imshow(S);

if ~exist('alpha','var')
    alpha = 0.02;
end
if ~exist('scale','var')
    scale = 1.5;
end

if ~exist('kappa','var')
    kappa = 2.0;
end
if ~exist('lambda','var')
    lambda = 2e-2;
end
S = im2double(Im);
betamax = 1e2;
fx = [1, -1];
fy = [1; -1];
[N,M,D] = size(Im);
sizeI2D = [N,M];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
Normin1 = fft2(S);
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
if D>1
    Denormin2 = repmat(Denormin2,[1,1,D]);
end
beta = 2*lambda;
% h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
% h = h/2.0
% mask_neg = h<0
% h = abs(h)
% mask_h_alpha = (h<alpha)
% h_alpha = ((h-alpha/alpha).^(scale))*alpha
% h(mask_h_alpha)=h_alpha(mask_h_alpha)
% h(mask_neg) = -h(mask_neg)
% 
% v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
% 
% v = v/2.0
% mask_neg = v<0
% v = abs(v)
% mask_v_alpha = (v<alpha)
% v_alpha = ((v-alpha/alpha).^(scale))*alpha
% v(mask_v_alpha)=v_alpha(mask_v_alpha)
% v(mask_neg) = -v(mask_neg)
% 
% h = h*2
% v = v*2
while beta < betamax
    Denormin   = 1 + beta*Denormin2;
    % h-v subproblem
    

    h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
    h = h/1.0
    mask_neg = h<0
    h = abs(h)
    mask_h_alpha = (h<alpha)
    h_alpha = (((h)/alpha).^(scale))*alpha
    
    h(mask_h_alpha)=h_alpha(mask_h_alpha)
    h(mask_neg) = -h(mask_neg)

    v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];

    v = v/1.0
    mask_neg = v<0
    v = abs(v)
    mask_v_alpha = (v<alpha)
    v_alpha = (((v)/alpha).^(scale))*alpha
    
    
    v(mask_v_alpha)=v_alpha(mask_v_alpha)
    v(mask_neg) = -v(mask_neg)

    h = h*1
    v = v*1
    
%     if D==1
%         t = (h.^2+v.^2)<=lambda/beta;
%     else
%         t = sum((h.^2+v.^2),3)<lambda/beta;
%         t = repmat(t,[1,1,D]);
%     end
%     h(t)=0; v(t)=0;
    % S subproblem
    Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
    Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
    FS = (Normin1 + beta*fft2(Normin2))./Denormin;
    S = real(ifft2(FS));
    beta = beta*kappa;
    fprintf('.');
    
end
fprintf('\n');
end