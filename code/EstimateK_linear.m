function [K, Hs] = EstimateK_linear(x, X)
% Linear parameter estimation of K
%
%   x:  2D points. n x 2 x N matrix, where n is the number of corners in
%   a checkerboard and N is the number of calibration images
%       
%   X:  3D points. n x 2 matrix, where n is the number of corners in a
%   checkerboard and assumes the points are on the Z=0 plane
%
%   imgs: calibration images. N x 1 cell, where N is the number of calibration images
%
%   K: a camera calibration matrix. 3 x 3 matrix.
%
%   Hs: a homography from the world to images. 3 x 3 x N matrix, where N is 
%   the number of calibration images. You can use est_homography.m to
%   estimate homographies.
%

%% Your code goes here
n_img = size(x,3);
V = zeros(2*n_img, 6);

for i = 1: n_img
    Hs(:,:,i) = est_homography(x(:,1,i),x(:,2,i),X(:,1),X(:,2));
    
    V(2*i-1,:) = get_v(Hs(:,:,i),1,1)' - get_v(Hs(:,:,i),2,2)';
    V(2*i,:) = get_v(Hs(:,:,i),1,2)';
end

[~, ~, temp] = svd(V);
B_col = temp(:,end);
B = [ B_col(1), B_col(2), B_col(4);
      B_col(2), B_col(3), B_col(5);
      B_col(4), B_col(5), B_col(6)];

py = (B(1,2)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
c = B(3,3) - ( B(1,3)^2+py*(B(1,2)*B(1,3)-B(1,1)*B(2,3)) )/B(1,1);
fy = sqrt( c*B(1,1)/(B(1,1)*B(2,2)-B(1,2)^2) );
fx = sqrt( c/B(1,1) );
s = -B(1,2)*fx^2*fy/c;
px = s*py/fy - B(1,3)*fx^2/c;
K = [fx, s, px;
      0, fy, py;
      0, 0, 1];

end

function v = get_v(H, i, j)
    
    v = [ H(1,i)*H(1,j), H(1,i)*H(2,j)+H(2,i)*H(1,j), H(2,i)*H(2,j), H(3,i)*H(1,j)+H(1,i)*H(3,j), H(3,i)*H(2,j)+H(2,i)*H(3,j),  H(3,i)*H(3,j)]';

end