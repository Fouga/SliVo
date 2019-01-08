function [ alpha, phi, delta ] = planeParams( P )
%[ alpha, phi, d ] = planeParams( P ) for a plane P(1)x+P(2)y+P(3)z+P(4)=0,
%finds the angle of the normal vector from the vertical (alpha, [0 90] deg), the
%rotatation of the normal vector (phi, [0 360] deg) and the shift along the vertical
%dimention from the (0,0,0) point. alpha and phi are in degrees. Alpha is
%taken such that the normal vector 'looks upwards'. Phi is the angle from
%the X-axis in the counterclockwise direction (imagine X-axis looking to
%the right, and y-axis looking up in a plane).

if P(3)<0
    P = -P;
end
A = P(1);
B = P(2);
C = P(3);
D = P(4);
N = sqrt(A^2 + B^2 + C^2); % length of the normal vector

alpha = acos(C/N);
phi = atan2(B, A);
delta = -D/C;

end

