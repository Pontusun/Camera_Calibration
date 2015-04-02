function [error, f] = GeoError(x, X, ks, K, Rs, ts)
% Measure a geometric error
%
%   x:  2D points. n x 2 x N matrix, where n is the number of corners in
%   a checkerboard and N is the number of calibration images
%       
%   X:  3D points. n x 2 matrix, where n is the number of corners in a
%   checkerboard and assumes the points are on the Z=0 plane
%
%   K: a camera calibration matrix. 3 x 3 matrix.
%
%   Rs: rotation matrices. 3 x 3 x N matrix, where N is the number of calibration images. 
%
%   ts: translation vectors. 3 x 1 x N matrix, where N is the number of calibration images. 
%
%   ks: radial distortion parameters. 2 x 1 matrix, where ks(1) = k_1 and
%   ks(2) = k_2.
%

%% Your code goes here
n_img = size(x,3);
n_points = size(x,1);

f = zeros(2,n_points,n_img);

%K(1,2) = 0;%?
px = K(1,3);
py = K(2,3);

for i = 1:n_img
    for j = 1:n_points
        ab_homo = [Rs(:,:,i), ts(:,1,i)] * [X(j,:), 0, 1]';
        a = ab_homo(1)/ab_homo(3);
        b = ab_homo(2)/ab_homo(3);
        r_sqr = a^2 + b^2 ;
        
        %ideal_homo = K*[Rs(:,:,i), ts(:,1,i)] * [X(j,:), 0, 1]';
        ideal_homo = K*[a; b; 1];
        u_ideal = ideal_homo(1)/ideal_homo(3);
        v_ideal = ideal_homo(2)/ideal_homo(3); 
        
        u_reproj = u_ideal + (u_ideal-px)*(ks(1)*r_sqr + ks(2)*r_sqr*r_sqr);
        v_reproj = v_ideal + (v_ideal-py)*(ks(1)*r_sqr + ks(2)*r_sqr*r_sqr);
        f(1,j,i) = x(j,1,i) - u_reproj;
        f(2,j,i) = x(j,2,i) - v_reproj;
    end
end

error = sum(sum(sum(f.^2)));

end