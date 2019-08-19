% AW1=RAC3leafMetaDriveLight2('LightZhuL6.txt',0,0);
% AW2=RAC3leafMetaDriveLight2('LightZhuL6.txt',0,1);
% AW3=RAC3leafMetaDriveLight2('LightZhuL6.txt',1,1);
% 
% AW4=RAC3leafMetaDriveLight2('LightZhuL7.txt',0,0);
% AW5=RAC3leafMetaDriveLight2('LightZhuL7.txt',0,1);
% AW6=RAC3leafMetaDriveLight2('LightZhuL7.txt',1,1);
% 
% 
% AW7=RAC3leafMetaDriveLight2('LightZhuL8.txt',0,0);
% AW8=RAC3leafMetaDriveLight2('LightZhuL8.txt',0,1);
% AW9=RAC3leafMetaDriveLight2('LightZhuL8.txt',1,1);
% 
% 
% AW10=RAC3leafMetaDriveLight2('LightZhuL9.txt',0,0);
% AW11=RAC3leafMetaDriveLight2('LightZhuL9.txt',0,1);
% AW12=RAC3leafMetaDriveLight2('LightZhuL9.txt',1,1);
% 
% AWL6=[AW1;AW2;AW3];
% AWL7=[AW4;AW5;AW6];
% AWL8=[AW7;AW8;AW9];
% AWL9=[AW10;AW11;AW12];
AW=RAC3leafMetaDriveLight2('LightZhuL5.txt',1,1);


 Lightd=importdata('LightZhuL5.txt');
 Lightdata=Lightd.data;
 lightM(:,1)=Lightdata(217:1220,1);
 lightM(:,2)=Lightdata(217:1220,2);
 figure;
 plot(lightM(:,1),lightM(:,2));
 xlabel('Time (min)');
 ylabel('PAR (\mumol m^-2 s^-^1')
 
  Lightdata2=importdata('bon18182.dat');
 %LightM0=Lightdata.data;
 global lightM;
 lightM2(:,1)=Lightdata2(624:1440,7)-6;
 lightM2(:,2)=Lightdata2(624:1440,31)*1E6/(2.35E5);
%   figure;
%  plot(lightM2(:,1),lightM2(:,2));
%  xlabel('Time (min)');
%  ylabel('PAR (\mumolm^-^-2 s^-^1')
 figure;
 subplot(1,2,1);  plot(lightM2(:,1),lightM2(:,2)); xlabel('Time (min)'); ylabel('PAR (\mumol m^-^2 s^-^1)');
 subplot(1,2,2); plot(lightM(:,1),lightM(:,2)); xlabel('Time (min)');ylabel('PAR (\mumol m^-^2 s^-^1)');
 
