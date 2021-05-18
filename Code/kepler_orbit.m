function [dx] = kepler_orbit(~,x,mu)

% Differential equation to describe a Keplerian orbit (two-body problem)

% PROTOTYPE:
    % [dx] = kepler_orbit(~,x,mu);
    
% INPUT:
    % x: position and velocity vector [1x6]
    % mu: planetary constant of orbiting planet   [km^3/s^2]

% OUTPUT:
    % output1: []
    % output2: []
    % output3: []
    
% CONTRIBUTORS:
    % Elisa De Astis
    % Danilo Forte
    % Lorenzo Guerreschi
    % Sabrina Saban

% VERSIONS:
    % 18/12/2018: First version
    
    
% System to solve:
% x2'=-mu*x1/norm(x1)^3
% x1'=x2
% x2: velocity vector
% x1: position vector

r = sqrt(x(1)^2+x(2)^2+x(3)^2);
dx = [x(4:6)
    -(mu/(r)^3).*x(1:3)];
    
end