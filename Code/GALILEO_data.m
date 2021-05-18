% Script to integrate and save NORAD TLE initial state of the GALILEO satellites 
% in the form of an SPICE SPK kernel using a Keplerian unpertrubed model
% up to an aribitrary epoch. The data is used to build up SPICE SPK
% kernel which are than validated against the initial data 

    
% INPUT:
   % SPK kernel containting an interval including the initial state of each
   % satellites 
% OUTPUT:
   % Text files containting epoch and state at each timestep for each
   % satellite


% VERSIONS:
   % 1/4/2021: First version


   
%% Loading the needed kernels
cspice_furnsh ( 'kernels\spk\TLE_GALILEO.bsp');
cspice_furnsh ( 'kernels\spk\2K25_GALILEO.bsp');
cspice_furnsh ( 'kernels\spk\de430.bsp');
cspice_furnsh ( 'kernels\lsk\naif0012.tls');


% Initial epoch
et0   = cspice_str2et('2020-DEC-10, 12:00:00  TDB');

% User imput final epoch
utctimEND = input ( 'Input ending UTC Time: ', 's' );

fprintf ( 'Converting ending UTC Time: %s\n', utctimEND )

etEND = cspice_str2et(utctimEND); %'2025-DEC-10, 12:00:00  TDB'

mu = 3.9860e+05; % Earth gravitational constant 

% Saving data flag
data = input ( 'Saving the data? ( 0 = no, 1 = y): ' );

if data == 1
    tspan = linspace(et0,etEND,1825*24); % 12 equally spaced point per orbit (24 a day) for 5 years
else
     tspan = [et0,etEND]; 
end

% Plot data flag
plot_flag = input ( 'Plottoing? ( 0 = no, 1 = y): ' );




for k = 1:24
    
    GALILEOID = -(650 + k);
    GALILEOID = char(string(GALILEOID));
    
    [state, ltime] = cspice_spkezr (GALILEOID, et0,'J2000', 'LT+S', '399');
    [T,STATE] = orbit_plot(state,tspan,mu,plot_flag);
    
    if data == 1
        
         state = [T,STATE];
         writematrix(state,GALILEOID);
    end
    
end

%% We have constructed the spk for the GALILEO constellation up to 2025. let's validate graphically the SPK

if data == 1
    etspan = linspace(etEND,etEND+86400*365*1.5,50000);
else
    etspan = T';
end


R_EA = 6371;
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
        
for k = 1:24
    
    GALILEOID = -(650 + k);
    GALILEOID = char(string(GALILEOID));
    
    [state, ltime] = cspice_spkezr (GALILEOID, etspan,'J2000', 'NONE', '399');
     

    plot3(state(1,:),state(2,:),state(3,:),'r')
    
    
end



%% 

cspice_kclear;