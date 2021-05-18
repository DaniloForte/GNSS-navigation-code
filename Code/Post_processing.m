%% Script to plot the results of the visibility analysis. 
% NOTE: THIS SCRIPT IS QUIETE IN AN HACK FORM. IT WORKS BUT SHOULD BE MODIFIED AND ADAPT MANUALLY
% FOR EACH USE CASE.


METAKR = 'kernels\GNSS_visibility.tm';

% Load the kernels that this program requires.  We
% will need a leapseconds kernel to convert input
% UTC time strings into ET.  We also will need
% SPK files for the GNSS constellations and the the solar system planets and satelittes.

cspice_furnsh ( METAKR );


% We also need the spk of the receiver trajectory

cspice_furnsh ('kernels\spk\MTO.bsp'); 
cspice_furnsh ('kernels\spk\LU_PATH.bsp'); 
cspice_furnsh ('kernels\spk\Halo_LUMIO.bsp'); 



% Plotting 

close all

epoch = GPS.epoch;

% Re-scaling of the time vector ( epoch(1) = 0)

time = epoch - epoch(1,1);


%% Plotting for each SV if it is tracked or not

% GPS

figure(1)

subplot(1,2,1)
hold on 
grid on

for k = 1:30
    
    plot(time/86400,GPS.acquired(:,k)*k,'o','color','#EDB120','linewidth',1.2);
       
end

V = axis;
V(2) =  time(end)/(24*3600) ;
V(3) = 0.05;
axis(V);

set(gca, 'FontSize', 18) ;  
title('GPS availability') ;
xlabel('time [days]', 'Interpreter', 'Latex')
ylabel('GPS SVs', 'Interpreter', 'Latex') 

% GALILEO

subplot(1,2,2)
hold on 
grid on

for k = 1:24
    
    plot(time/86400,GALILEO.acquired(:,k)*k,'o','color','#0072BD','linewidth',1.2);
       
end

V = axis;
V(2) =  time(end)/(24*3600);
V(3) = 0.05;
V(4) = 24;
axis(V);

set(gca, 'FontSize', 18) ;  
title('GALILEO availability') ;
xlabel('time [days]', 'Interpreter', 'Latex')
ylabel('Galileo SVs', 'Interpreter', 'Latex') 

% GLONASS

% subplot(2,2,3)
% hold on 
% grid on
% 
% for k = 1:24
%     
%     plot(time/86400,GLONASS.acquired(:,k)*k,'o','color','#7E2F8E','linewidth',1.2);
%        
% end
% 
% V = axis;
% V(2) =  time(end)/(24*3600);
% V(3) = 0.05;
% V(4) = 24;
% axis(V);
% 
% set(gca, 'FontSize', 18) ;  
% title('GLONASS availability') ;
% xlabel('time [days]', 'Interpreter', 'Latex')
% ylabel('GLONASS SVs', 'Interpreter', 'Latex') 
% 
% % BEIDOU III
% 
% subplot(2,2,4)
% hold on 
% grid on
% 
% for k = 1:30
%     
%     plot(time/86400,BEIDOU.acquired(:,k)*k,'o','color','#D95319','linewidth',1.2);
%        
% end
% 
% V = axis;
% V(2) =  time(end)/(24*3600) ;
% V(3) = 0.05;
% axis(V);
% 
% set(gca, 'FontSize', 18) ;  
% title('BEIDOU-III availability') ;
% xlabel('time [days]', 'Interpreter', 'Latex')
% ylabel('BeiDou-3 SVs', 'Interpreter', 'Latex') 


%% Plotting the total number of visible SVs


figure(2)

hold on
grid on

% Total tracked GPS SVs

% GPS.n_tracked(GPS.n_tracked >12 ) = 12;
plot_nGPS = plot(time/60,GPS.n_tracked,'color','#EDB120','linewidth',1.5);

% Total tracked GALILEO SVs

% GALILEO.n_tracked(GALILEO.n_tracked >12 ) = 12;
plot(time/60,GALILEO.n_tracked,'color','#0072BD','linewidth',1.5);

% GPS + GALILEO tracked SVs
plot(time/60,(GPS.n_tracked + GALILEO.n_tracked),'color','#7E2F8E','linewidth',1.5);


% 4 GNSS tracked SVs
% plot(time/60,(GPS.n_tracked + GALILEO.n_tracked + GLONASS.n_tracked + BEIDOU.n_tracked),'color','#D95319','linewidth',1.5);

% simple line at the 4 SVs threshold ( required for the PVT solution)

plot([time(1) time(end)]/60,[4,4],'g--','linewidth',1.5);

set(gca, 'FontSize', 20) ;  
title('Total SVs above threshold') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel('\# of SVs visible ', 'Interpreter', 'Latex') ;
legend('GPS','GALILEO','GPS + GALILEO','4 visible for PVT solution')

V = axis;
V(2) =  time(end)/(60);
V(4) = 30;
axis(V);

%% C_N0 

figure(3)
hold on
grid on

subplot(1,2,1)

% GPS

plot(time/60,GPS.C_N0(:,28),'color','#EDB120','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('GPS Block III-A-1 SV') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel('$ C/N_0$ [dB - Hz] ', 'Interpreter', 'Latex') ;

legend('BIII-A-1') ;

V = axis;
V(2) =  time(end)/(60);
axis(V);


subplot(1,2,2)

% GALILEO

plot(time/60,GALILEO.C_N0(:,12),'color','#0072BD','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel('$ C/N_0$ [dB - Hz] ', 'Interpreter', 'Latex') ;

legend('GSAT-0215') ;
V = axis;
V(2) =  time(end)/(60);
axis(V);

% subplot(4,1,3)
% 
% % GLONASS
% 
% plot(time/60,GLONASS.C_N0(:,8),'color','#7E2F8E','linewidth',1.5);
% 
% set(gca, 'FontSize', 18) ;  
% xlabel('time [min]', 'Interpreter', 'Latex')
% ylabel('$ C/N_0$ [dB - Hz] ', 'Interpreter', 'Latex') ;
% 
% legend('R 758') ;
% V = axis;
% V(2) =  time(end)/(60);
% axis(V);
% 
% subplot(4,1,4)
% 
% % BEIDOU III
% 
% plot(time/60,BEIDOU.C_N0(:,12),'color','#D95319','linewidth',1.5);
% 
% set(gca, 'FontSize', 18) ;  
% xlabel('time [min]', 'Interpreter', 'Latex')
% ylabel('$ C/N_0$ [dB - Hz] ', 'Interpreter', 'Latex') ;
% 
% legend('BDS M12') ;
% V = axis;
% V(2) =  time(end)/(60);
% axis(V);

%% Earth-RX distance with time

dist = zeros(1,numel(epoch));

R_E = cspice_bodvrd ( 'EARTH', 'RADII',3);
R_E = mean(R_E);

for s = 1 : numel(epoch)
    
    d = cspice_spkpos ( ID_RX, epoch(s)','J2000', 'NONE', '399');
    dist(s) = norm(d) - R_E;

end

figure(4)
hold on
grid on

plot(time/(24*3600),dist/R_E,'color','#0072BD','linewidth',2);

set(gca, 'FontSize', 20) ;  
title('Receiver Altitude') ;
xlabel('Time [days]', 'Interpreter', 'Latex')
ylabel(' \# Earth radii', 'Interpreter', 'Latex') ;

if ID_RX == -516
    
    legend('MTO')
    
end

V = axis;
V(2) =  time(end)/(24*3600);
axis(V);


%% GDOP

figure(5)
subplot(1,3,1)
hold on
grid on

plot(time/60,GPS.GDOP,'color','#EDB120','linewidth',1.5);

plot(time/60,GALILEO.GDOP,'color','#0072BD','linewidth',1.5);

plot(time/60,MULTI_GNSS.GPSandGAL_GDOP,'color','#7E2F8E','linewidth',1.5);

% plot(time/60,MULTI_GNSS.ALL_GNSS_GDOP,'color','#D95319','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Geometric dilution of precision (L5/E5a)') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS','GALILEO','GPS + GALILEO')

V = axis;
V(2) =  time(end)/(60);
if V(4) > 10000
   V(4) = 200000;
end
axis(V);

subplot(1,3,2)
hold on
grid on

plot(time/60,MULTI_GNSS.GPSandGAL_GDOP,'color','#7E2F8E','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS + GALILEO')

V = axis;
V(2) =  time(end)/(60) ;
V(4) = 5000;
axis(V);

% 4 GNSS

subplot(1,3,3)
hold on
grid on

% plot(time/60,MULTI_GNSS.ALL_GNSS_GDOP,'color','#D95319','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Geometric dilution of precision') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

%legend('4 GNSS')

V = axis;
V(2) =  time(end)/(60) ;
V(4) = 2500;
axis(V);


%% GDOP L5

figure(5)
subplot(1,2,1)
hold on
grid on

plot(time/60,GPS.GDOP,'color','#EDB120','linewidth',1.5);

plot(time/60,GALILEO.GDOP,'color','#0072BD','linewidth',1.5);

plot(time/60,MULTI_GNSS.GPSandGAL_GDOP,'color','#7E2F8E','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS','GALILEO','GPS + GALILEO')

V = axis;
V(2) =  time(end)/(60) ;
if V(4) > 1000
   V(4) = 200000;
end
axis(V);

subplot(1,2,2)
hold on
grid on

plot(time/60,MULTI_GNSS.GPSandGAL_GDOP,'color','#7E2F8E','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS + GALILEO')

V = axis;
V(2) =  time(end)/(60);
V(4) = 10000;
axis(V);


%% Doppler shift ( NOT QUIET INTERESTING)


% figure(6)
% hold on
% grid on
% 
% subplot(2,1,1)
% 
% plot(time/60,GPS.D_shift(:,15)/1000,'color','#EDB120','linewidth',1.5);
% 
% set(gca, 'FontSize', 18) ;  
% title('Selected SVs doppler shift') ;
% xlabel('time [min]', 'Interpreter', 'Latex')
% ylabel('frequency shift [kHz] ', 'Interpreter', 'Latex') ;
% 
% legend('BIIRM-6') ;
% 
% subplot(2,1,2)
% 
% plot(time/60,GALILEO.D_shift(:,12)/1000,'color','#0072BD','linewidth',1.5);
% 
% set(gca, 'FontSize', 18) ;  
% xlabel('time [min]', 'Interpreter', 'Latex')
% ylabel('frequency shift [kHz] ', 'Interpreter', 'Latex') ;
% 
% legend('GSAT-0215') ;


%% LNNSS

LNNSS = input ( 'Data including LNNSS? [Y/N]:', 's' );

if isequal(LNNSS,'Y')

%% Plotting for each SV if it is tracked or not

figure(7)

subplot(1,2,1)
hold on 
grid on

for k = 1:6
    
    plot(time/86400,LNNSS_6SV.acquired(:,k)*k,'o','color','#EDB120','linewidth',1.2);
       
end

V = axis;
V(3) = 0.05;
axis(V);

set(gca, 'FontSize', 18) ;  
title('6SV LNNSS availability') ;
xlabel('time [days]', 'Interpreter', 'Latex')
ylabel('SVs', 'Interpreter', 'Latex') 

% GALILEO

subplot(1,2,2)
hold on 
grid on

for k = 1:8
    
    plot(time/86400,LNNSS_8SV.acquired(:,k)*k,'o','color','#0072BD','linewidth',1.2);
       
end

V = axis;
V(3) = 0.05;
axis(V);

set(gca, 'FontSize', 18) ;  
title('8SV LNNSS availability') ;
xlabel('time [days]', 'Interpreter', 'Latex')
ylabel('SVs', 'Interpreter', 'Latex') 


%% Plotting the total number of visible SVs


figure(8)
subplot(1,2,1)
hold on
grid on

% Total tracked 6SV LNNSS SVs
% GPS + GALILEO tracked SVs
plot(time/60,(GPS.n_tracked + GALILEO.n_tracked),'color','#7E2F8E','linewidth',1.5);

plot(time/60,LNNSS_6SV.n_tracked,'color','#EDB120','linewidth',1.5);

plot(time/60,(GPS.n_tracked + GALILEO.n_tracked + LNNSS_6SV.n_tracked),'color','#D95319','linewidth',1.5);

% simple line at the 4 SVs threshold ( required for the PVT solution)

plot([time(1) time(end)]/60,[4,4],'g--','linewidth',1.5);

set(gca, 'FontSize', 16) ;  
title('Total SVs above threshold') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel('\# of SVs visible ', 'Interpreter', 'Latex') ;
legend('GPS+Galileo','6SVs LNNSS','GPS+Galileo+6SVs LNNSS','4 visible for PVT solution')

V = axis;
V(2) =  time(end)/(60);
V(4) = 15;
axis(V);


subplot(1,2,2)
hold on
grid on
% Total tracked 8SV MOON GALILEO SVs
plot(time/60,(GPS.n_tracked + GALILEO.n_tracked),'color','#7E2F8E','linewidth',1.5);

plot(time/60,LNNSS_8SV.n_tracked,'color','#0072BD','linewidth',1.5);

plot(time/60,(GPS.n_tracked + GALILEO.n_tracked + LNNSS_8SV.n_tracked),'color','#D95319','linewidth',1.5);

% simple line at the 4 SVs threshold ( required for the PVT solution)

plot([time(1) time(end)]/60,[4,4],'g--','linewidth',1.5);

set(gca, 'FontSize', 16) ;  
title('Total SVs above threshold') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel('\# of SVs visible ', 'Interpreter', 'Latex') ;
legend('GPS+Galileo','8SVs LNNSS','GPS+Galileo+8SVs LNNSS','4 visible for PVT solution')

V = axis;
V(2) =  time(end)/(60);
V(4) = 15;
axis(V);

%% C_N0 

figure(9)
hold on
grid on

%LNNSS SVs

plot(time/60,LNNSS_8SV.C_N0(:,2),'color','#EDB120','linewidth',1.5);

plot(time/60,LNNSS_8SV.C_N0(:,5),'color','#0072BD','linewidth',1.5);

plot([time(1) time(end)]/60,[26.5,26.5],'g--','linewidth',1.5);


set(gca, 'FontSize', 18) ;  
title('Selected SVs carrier to noise ratio') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel('$ C/N_0$ [dB - Hz] ', 'Interpreter', 'Latex') ;

legend('LNNSS SAT-2','LNNSS SAT-5','Threshold') ;

V = axis;
V(2) =  time(end)/(60);
V(3) = 25;
axis(V);

%% GDOP

figure(10)

subplot(1,3,1)
hold on
grid on

plot(time/60,LNNSS_6SV.GDOP,'color','#EDB120','linewidth',1.5);

plot(time/60,LNNSS_8SV.GDOP,'color','#0072BD','linewidth',1.5);


set(gca, 'FontSize', 16) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('LNNSS 6SVs','LNNSS 8SVs')

V = axis;
V(4) = 5000;
axis(V);

subplot(1,3,2)
hold on
grid on

plot(time/60,MULTI_GNSS.GPS_GAL_LNNSS6SV_GDOP,'color','#D95319','linewidth',1.5);


set(gca, 'FontSize', 16) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS+GALILEO+6SVs LNNSS')

V = axis;
V(4) = 5000;
axis(V);

subplot(1,3,3)
hold on
grid on

plot(time/60,MULTI_GNSS.GPS_GAL_LNNSS8SV_GDOP,'color','#D95319','linewidth',1.5);


set(gca, 'FontSize', 16) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS+GALILEO+8SVs LNNSS')

V = axis;
V(4) = 5000;
axis(V);

end

%% EPHEMERIDIS VISIBILITY

EV = input ( 'Ephemeridis visibility plot? [Y/N]:', 's' );

if isequal(EV,'Y')
    
%% Plotting for each SV if it is tracked or not

% GPS

figure(11)

subplot(1,2,1)
hold on 
grid on

for k = 1:30
    
    plot(time/86400,GPS.ephemeridis_tracking(:,k)*k,'o','color','#EDB120','linewidth',1.2);
       
end

V = axis;
V(3) = 0.05;
axis(V);

set(gca, 'FontSize', 18) ;  
title('GPS ephemeridis + tracking') ;
xlabel('time [days]', 'Interpreter', 'Latex')
ylabel('SVs', 'Interpreter', 'Latex') 

% GALILEO

subplot(1,2,2)
hold on 
grid on

for k = 1:24
    
    plot(time/86400,GALILEO.ephemeridis_tracking(:,k)*k,'o','color','#0072BD','linewidth',1.2);
       
end

V = axis;
V(3) = 0.05;
V(4) = 24;
axis(V);

set(gca, 'FontSize', 18) ;  
title('GALILEO ephemeridis + tracking') ;
xlabel('time [days]', 'Interpreter', 'Latex')
ylabel('SVs', 'Interpreter', 'Latex') 


%% Plotting the total number of visible SVs


figure(12)

hold on
grid on

% Total tracked GPS SVs

% GPS.n_tracked(GPS.n_tracked >12 ) = 12;
plot_nGPS = plot(time/60,GPS.ephemeridis_n_tracked,'color','#EDB120','linewidth',1.5);

% Total tracked GALILEO SVs

% GALILEO.n_tracked(GALILEO.n_tracked >12 ) = 12;
plot(time/60,GALILEO.ephemeridis_n_tracked,'color','#0072BD','linewidth',1.5);

% GPS + GALILEO tracked SVs
plot(time/60,(GPS.ephemeridis_n_tracked + GALILEO.ephemeridis_n_tracked),'color','#7E2F8E','linewidth',1.5);


% simple line at the 4 SVs threshold ( required for the PVT solution)

plot([time(1) time(end)]/60,[4,4],'g--','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Total SVs above threshold and with decoded ephemeridis') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel('\# of SVs visible ', 'Interpreter', 'Latex') ;
legend('GPS','GALILEO','GPS + GALILEO','4 visible for PVT solution')

V = axis;
V(2) = time(end)/60;
V(4) = 25;
axis(V);

%% GDOP

figure(13)
subplot(1,2,1)
hold on
grid on

plot(time/60,GPS.ephemeridis_GDOP,'color','#EDB120','linewidth',1.5);

plot(time/60,GALILEO.ephemeridis_GDOP,'color','#0072BD','linewidth',1.5);

plot(time/60,MULTI_GNSS.ephemeridis_GPSandGAL_GDOP,'color','r','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS','GALILEO','GPS + GALILEO')

V = axis;

if V(4) > 20000
   V(4) = 50000;
end
axis(V);

subplot(1,2,2)
hold on
grid on

plot(time/60,MULTI_GNSS.ephemeridis_GPSandGAL_GDOP,'color','r','linewidth',1.5);

set(gca, 'FontSize', 18) ;  
title('Geometric dilution of precision ') ;
xlabel('time [min]', 'Interpreter', 'Latex')
ylabel(' GDOP ', 'Interpreter', 'Latex') ;

legend('GPS + GALILEO')

V = axis;
V(4) = 1e4;
axis(V);

    
end



%% Clearing the kernel pool
cspice_kclear;


