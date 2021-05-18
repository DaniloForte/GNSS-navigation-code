%% B-IIF Block

% mean gain pattern data for BII-F SV avaible 

BIIF = readmatrix('Block-IIF_txgain_gpsace.txt');

% Azimuth

phi = 0:1:359;

% Off-bore angle

theta = 16:1:90;
theta = [13;14;15;theta'];

BIIF_GAIN = zeros(77,360);

% The data avaible is down to 16 deg of off-bore angle. To complete the
% pattern mean azimuthal values are included from other sources

BIIF_GAIN(1,:) = 14.15;
BIIF_GAIN(2,:) = 13.61;
BIIF_GAIN(3,:) = 13;


figure(1)
hold on 
grid on

% Extracting data from the table

for k = 4:78
    
 point = BIIF(:,1) == theta(k);   
 
 BIIF_GAIN(k,:) = BIIF(point,3);
                            
end

% Plotting some azimuthal slices

for z = 10:10:359
    
plot(theta,BIIF_GAIN(:,1),'linewidth',1.5)   

plot(theta,BIIF_GAIN(:,z),'linewidth',1.5)   
end

set(gca, 'FontSize', 18) ;  
title('B-IIF L1 antenna gain pattern') ;
xlabel('off-bore angle [Deg]', 'Interpreter', 'Latex')
ylabel('Gain [dB]', 'Interpreter', 'Latex') ;


%% B-IIR

% data for four BIIR SV avaible ( modernized panel)

SVN47 = readtable('L1_BIIR_47.xlsx','range','A4:AK94');
SVN59 = readtable('L1_BIIR_59.xlsx','range','A4:AK94');
SVN60 = readtable('L1_BIIR_60.xlsx','range','A4:AK94');
SVN61 = readtable('L1_BIIR_61.xlsx','range','A4:AK94');

% Azimuth

phi = 0:10:350;

% Off-bore angle

theta = SVN47.(1);

BIIR_GAIN = zeros(91,36);

% Gain correction factor

GCF = 1.3; % [dB]

% Mean value of the directivity of the block in each azimuth slice avaible

figure(2)
hold on 
grid on

for k = 2:37
    
 BIIR_GAIN(:,k-1) = (SVN47.(k) + SVN59.(k) + SVN60.(k) + SVN61.(k))/4 - GCF;  
 
  plot(theta,BIIR_GAIN(:,k-1),'linewidth',1.5)      
                         
end

set(gca, 'FontSize', 18) ;  
title('B-IIR L1 antenna gain pattern') ;
xlabel('off-bore angle [Deg]', 'Interpreter', 'Latex')
ylabel('Gain [dB]', 'Interpreter', 'Latex') ;

%% B-IIR-M

% data for nine B-IIR-M SV avaible

SVN48 = readtable('L1_IIRM_48.xlsx','range','A4:AK94');
SVN49 = readtable('L1_IIRM_49.xlsx','range','A4:AK94');
SVN50 = readtable('L1_IIRM_50.xlsx','range','A4:AK94');
SVN51 = readtable('L1_IIRM_51.xlsx','range','A4:AK94');
SVN52 = readtable('L1_IIRM_52.xlsx','range','A4:AK94');
SVN53 = readtable('L1_IIRM_53.xlsx','range','A4:AK94');
SVN55 = readtable('L1_IIRM_55.xlsx','range','A4:AK94');
SVN57 = readtable('L1_IIRM_57.xlsx','range','A4:AK94');
SVN58 = readtable('L1_IIRM_58.xlsx','range','A4:AK94');


BIIRM_GAIN = zeros(91,36);

% Gain correction factor

GCF = 1.4; % [dB]

% Mean value of the directivity of the block in each azimuth slice avaible

figure(3)
hold on 
grid on

for k = 2:37
    
 BIIRM_GAIN(:,k-1) = (SVN48.(k) + SVN49.(k) + SVN50.(k) + SVN51.(k) + SVN52.(k) + ...
                             SVN53.(k)  + SVN55.(k) + SVN57.(k) + SVN58.(k))/9 - GCF;  
 
 plot(theta,BIIRM_GAIN(:,k-1),'linewidth',1.5)      
                         
end

set(gca, 'FontSize', 18) ;  
title('BIIR-M L1 antenna gain pattern') ;
xlabel('off-bore angle [Deg]', 'Interpreter', 'Latex')
ylabel('Gain [dB]', 'Interpreter', 'Latex') ;


