function [T,X] = orbit_plot(STATE,tspan,mu,plot_flag)

% Function to plot the orbit according to tspan

% PROTOTYPE:
   % [T,X] = orbit_plot(x,v,tspan,mu);
    
% INPUT:
   % x: position vector [1x3]             [km/s]
   % v: velocity vector [1x3]             [km/s]
   % tspan: time vector
   % (!! use linspace(0,t,n) with n = number of time intervals for closed orbit)
   % mu: planet gravitational parameter   [km^3/s^2]
   % plot: flag to use the plot section of the function ( 0 -> no plot, 1
   % -> plot)
% OUTPUT:
   % T: time vector
   % X: position and velocity according to time vector


% VERSIONS:
   % 18/12/2018: First version


options = odeset('RelTol', 1e-13,'AbsTol',1e-14 );

[T,X] = ode113(@(t,x) kepler_orbit(t,x,mu),tspan,STATE,options);

R_EA = 6371;
if plot_flag == 1
    TERRA = imread('EarthTexture.jpg','jpg');
    props.FaceColor = 'texture';
    props.EdgeColor = 'texture';
    props.FaceLighting = 'phong';
    props.Cdata = TERRA;
    Center = [0; 0; 0];
    [XX, YY, ZZ] = ellipsoid(0, 0, 0, R_EA, R_EA, R_EA, 40);
    surface(-XX, -YY, -ZZ, props);
    axis([-40000, 40000, -40000, 40000, -40000, 40000]);
    hold on
    axis equal

    plot3(X(:,1),X(:,2),X(:,3),'b')
end

end