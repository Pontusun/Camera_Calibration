function [ Rs, ts ] = EstimateRt_linear( Hs, K )
% Linear parameter estimation of R and t
%
%   K: a camera calibration matrix. 3 x 3 matrix.
%
%   Hs: a homography from the world to images. 3 x 3 x N matrix, where N is 
%   the number of calibration images. 
%
%   Rs: rotation matrices. 3 x 3 x N matrix, where N is the number of calibration images. 
%
%   ts: translation vectors. 3 x 1 x N matrix, where N is the number of calibration images. 
%

%% Your code goes here.
n_img = size(Hs,3);
Rs = zeros(3, 3, n_img);
ts = zeros(3, 1, n_img);

K_inv = inv(K);

for i = 1:n_img
    z = sqrt(sum((K_inv*Hs(:,1,i)).^2));
    Rs(:,1,i) = K_inv*Hs(:,1,i)/z ;
    Rs(:,2,i) = K_inv*Hs(:,2,i)/z ;
    Rs(:,3,i) = cross(Rs(:,1,i), Rs(:,2,i)) ;
    ts(:,1,i) = K_inv*Hs(:,3,i)/z ;
end


end

