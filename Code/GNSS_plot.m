% Script to plot the four GNSS constellations and the  Heavy LNNSS

% Loading the kernels needed
METAKR = 'kernels\GNSS_visibility.tm';

cspice_furnsh ( METAKR );
cspice_furnsh('kernels/spk/8SV_H_LNNSS.bsp');

%% Selecting the time window to plot
utctimIN = '2020-DEC-10, 12:00:01  TDB';

fprintf ( 'Converting initial UTC Time: %s\n', utctimIN )

et0 = cspice_str2et(utctimIN); % 2020-DEC-10, 12:00:01  TDB


utctimEND = input ( 'Input ending UTC Time: ', 's' );

fprintf ( 'Converting ending UTC Time: %s\n', utctimEND )

etEND = cspice_str2et(utctimEND); % 2025-DEC-10, 12:00:00  TDB


etspan = linspace(et0,etEND,(etEND - et0)/600);

%% plotting

R_EA = cspice_bodvrd ( 'EARTH', 'RADII',3);
R_EA = mean(R_EA);

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
  
        
 % GPS constellation 

for k = 1:30

        GPSID = -(600 + k);
        GPSID = char(string(GPSID));

        [state, ltime] = cspice_spkezr (GPSID, etspan,'J2000', 'LT+S', '399');


        GPS(k) = plot3(state(1,:),state(2,:),state(3,:),'color','#EDB120','linewidth',1.2);
        
                 plot3(state(1,1),state(2,1),state(3,1),'o','color','#EDB120','linewidth',2);

end


 % GALILEO constellation 
        
for k = 1:24
    
    GALILEOID = -(650 + k);
    GALILEOID = char(string(GALILEOID));
    
    [state, ltime] = cspice_spkezr (GALILEOID, etspan,'J2000', 'LT+S', '399');
     

    GALILEO(k) = plot3(state(1,:),state(2,:),state(3,:),'color','#0072BD','linewidth',1.2);
    
                 plot3(state(1,1),state(2,1),state(3,1),'o','color','#0072BD','linewidth',2);
    
    
end


 % GLONASS constellation 

for k = 1:24
    
    GLONASSID = -(700 + k);
    GLONASSID = char(string(GLONASSID));
    
    [state, ltime] = cspice_spkezr (GLONASSID, etspan,'J2000', 'LT+S', '399');
     

    GLONASS(k) = plot3(state(1,:),state(2,:),state(3,:),'color','#7E2F8E','linewidth',1.2);
    
                 plot3(state(1,1),state(2,1),state(3,1),'o','color','#7E2F8E','linewidth',2);
    
    
end


 % BEIDOU III constellation 

for k = 1:30
    
    BEIDOUSID = -(750 + k);
    BEIDOUSID = char(string(BEIDOUSID));
    
    [state, ltime] = cspice_spkezr (BEIDOUSID, etspan,'J2000', 'LT+S', '399');
     

    BEIDOU(k) = plot3(state(1,:),state(2,:),state(3,:),'color','#D95319','linewidth',1.2);
    
                plot3(state(1,1),state(2,1),state(3,1),'o','color','#D95319','linewidth',2);
    
    
end


% Heavy LNNSS constellation

for k = 1:8
    
   LNNSSID = -(800 + k);
   LNNSSID = char(string(LNNSSID));
    
    [state, ltime] = cspice_spkezr (LNNSSID, etspan,'J2000', 'LT+S', '399');
     

    LNNSS(k) = plot3(state(1,:),state(2,:),state(3,:),'color','#77AC30','linewidth',1.2);
    
                plot3(state(1,1),state(2,1),state(3,1),'o','color','#77AC30','linewidth',2);
    
    
end




legend([GPS(1) GALILEO(1) GLONASS(1) BEIDOU(1) LNNSS(1)],{'GPS','GALILEO','GLONASS','BEIDOU III','LNNSS'})

cspice_kclear;
        