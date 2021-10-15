function [UT]=endstep(UT,U,stepnumber)

UT(:,:,stepnumber)=U; % 3D array for total displacement
