
%cassava simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RAC3leafMetaDrive('bon18182.dat',0,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AandW=RAC3leafMetaDrive(inputfile,Pst,PRca)
%clear all;
Convert=1E6/(2.35E5);
STdata=importdata('cassavaP2.txt');
CasSTdata=STdata.data;
% Lightdata=importdata('bon18182.dat');
 Lightdata=importdata(inputfile);
 %LightM0=Lightdata.data;
 global lightM;
 lightM(:,1)=Lightdata(624:1440,7)-10.4;
 lightM(:,2)=Lightdata(624:1440,31)*Convert;
 for i=1:817
     if lightM(i,2)<0
         lightM(i,2)=0;
     end
 end
 lightM(:,3)=Lightdata(624:1440,39);
 lightM(:,4)=Lightdata(624:1440,41);
for i=1:11
i
simNo=i;
dVcmax=CasSTdata(simNo,1);
dJmax=CasSTdata(simNo,2);
dKd=CasSTdata(simNo,3);
dKi=CasSTdata(simNo,4);
dSlop=CasSTdata(simNo,5);
dInter=CasSTdata(simNo,6);

Begin = 1;
%fin = SYSInitial(Begin);
global RAInteg;
RAInteg=PRca;% 1=consider Rac; 0=Rubisco always actived 
global activase; 
activase = 0.002*3; 
global RuACT_EPS_com;
RuACT_EPS_com = 1;

%%%%%%%%%%%%%RuACT_EPS_com
global BallBerryInterceptC3;
global BallBerrySlopeC3
BallBerryInterceptC3=1.6*dInter;%0.008;%WY201804%Ball 1988
BallBerrySlopeC3=1.6*dSlop*100;%WY201804  %9.29;%Ball 1988%10.5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ki_Gs;
global kd_Gs;
global GsResponse;
ki_Gs=1/dKi;%0.9*60;%6.9*60*exp(1);%2*60;0.9 %WY201804 cassava
kd_Gs=1/dKd;%6.5*60*exp(1);%1*60;4.1
GsResponse=Pst; %%if GsResponse=0 Ball berry model no time dependent Gs response ; %%if GsResponse=1 time dependent Gs response, using ki_Gs and kd_Gs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Para_mata;
global PhotosynthesisType;
Para_mata=1;%%if Para_mata=1, C3 Metabolic model and Gs model integrated  if Para_mata=0 Farquhar mdoel and gs model
PhotosynthesisType=1;% 1:C3 model 2:C4 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leaf model parameters

global Air_CO2;
global Air_O2;

global WeatherWind;
global Radiation_PAR;
global Vcmax25;
global Jmax25;
global WeatherTemperature;
global WeatherRH;
global Radiation_NIR;
global Radiation_LW;
global PhiLeaf

 
WeatherTemperature=25;
Air_CO2=400;
Air_O2=210.0;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=2000/Convert*0.85*0.85;%10*i;
Radiation_NIR=0;
Radiation_LW=0;
Vcmax25=dVcmax;
Jmax25=dJmax;
PhiLeaf=-0;%Mpa Min=-2.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global O2_cond


PCO2=Air_CO2;
O2_cond=Air_O2/1000*1.26;
Plight=Radiation_PAR*Convert;
PTemp=WeatherTemperature;

global CI
CI=Air_CO2*0.7/(3 * 10^4);
global CO2_Env;
global CO2_cond;
global LI;
global Jmax;
global alfa;
global fc;
global Theta;
global beta;
CO2_Env=PCO2;
CO2_cond = CO2_Env/(3 * 10^4);%CO2_Env*0.7/(3 * 10^4);% for leaf model input is CI
LI=Plight/1000;
Jmax=dJmax/1000;
alfa=0.85;
fc=0.15;
Theta=0.7;
beta=0.7519;



global Tp;
Tp=PTemp;
Begin = 1;
global tglobal;     % The total running time
tglobal =13.6*3600;
global options1 
options1 = odeset('RelTol',1e-4);
time = tglobal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variables used for obtaining flux and concentration data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     

global PS_OLD_TIME;
global PS_TIME_N;
global PS_VEL;
PS_OLD_TIME = 0;
PS_TIME_N= 0;
PS_VEL = zeros(1,1);

global PR_OLD_TIME;
global PR_TIME_N;
global PR_VEL;
PR_OLD_TIME = 0;
PR_TIME_N = 1;
PR_VEL = zeros(1,1);

global SUCS_OLD_TIME;
global SUCS_TIME_N;
global SUCS_VEL;
global SUCS_CON;

SUCS_OLD_TIME = 0;
SUCS_TIME_N = 1;
SUCS_VEL = zeros(1,3);    % Clean memory
SUCS_CON = zeros(3,1);    % Clean memory

global TIME_N;
global OLD_TIME;
global Gs_VEL;
TIME_N=0;
OLD_TIME=0;
Gs_VEL=zeros(1,10);
global RuACT_OLD_TIME;
global RuACT_TIME_N;
global RuACT_VEL;
global RuACT_CON;

RuACT_OLD_TIME = 0;
RuACT_TIME_N = 1;

RuACT_VEL = zeros(1,3);    % Clean memory
RuACT_CON = zeros(3,1);    % Clean memory

%%%%%%%%%%%%%%%%%%%%%%%%
%   Initialation step %
%%%%%%%%%%%%%%%%%%%%%%%%

global GP; 
GP = 0;
% Begin = 1;
% PS_PRs = zeros(23,1);
%CMs = CM_Ini(Begin);
CMs = RAC3leafMetaIni(Begin);

%ModelComb = IModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

global PR_PS_com;             % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 1;

global ATPActive;
ATPActive = 0;

global PSPR_SUCS_com;        % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

global RedoxReg_RA_com;
RedoxReg_RA_com = 0;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step %
%%%%%%%%%%%%%%%%%%%%%%%%

CM_Param = 0;

SUCS_Param = zeros(2,1);
SUCS_Param(1) = 1;
SUCS_Param(2) = 1;

PS_PR_Param = 0;


[Tt,d] = ode15s(@RAC3leafMetaMB,[0:1:time],CMs, options1,PS_PR_Param, SUCS_Param);   
%[Tt,d] = ode23tb(@CM_mb,[0,time],CMs, options1,PS_PR_Param, SUCS_Param);   
global d_plot;
d_plot=d;

global Tt_plot;
Tt_plot = Tt;
global Result;
Result =[Tt,d]; 
Eb=d(:,3);
gm=0.7;
Sc=3*10^4;
%vinf=gm*Sc*10^(-3)*(Ci/(3 * 10^4)-MC_CO2);
%vH2Ob=Gbw*(Eb-Ea)/(Pressure / 1000.0)*10^6.0;
global Gs_VEL;

Gsw=1.6*Result(:,5);
ESaturation = 0.611 * exp(17.502 * Result(:,6)./ (Result(:,6) + 240.97));
Eb=Result(:,4);
Pressure=101325.0;
% % % % % figure;
% % % % % subplot(2,3,1); plot(Result(:,1)/60,Result(:,2),'k');title('Ci');ylim([0,Air_CO2]);
% % % % % subplot(2,3,2); plot(Result(:,1)/60,Result(:,6),'k');title('Tleaf');
% % % % % subplot(2,3,3); plot(Result(:,1)/60,Result(:,5),'k');title('Gs');
% % % % % subplot(2,3,4); plot(Gs_VEL(:,1)/60,Gs_VEL(:,2),'k.');hold on;plot(Result(:,1)/60,gm*Sc*(Result(:,2)/(3 * 10^4)-Result(:,9)),'r');title('A');ylim([0,40])
% % % % % subplot(2,3,5); plot(Gs_VEL(:,1)/60,Gs_VEL(:,9),'k.');hold on;plot(Result(:,1)/60,Gsw.*(ESaturation-Eb)./(Pressure / 1000.0)*10^6.0);title('Transpiration');
% % % % % subplot(2,3,6); plot(RuACT_CON(:,1)/60,RuACT_CON(:,5),'k.');title('RA percent');
% % % % % 
% % % % % % % 
% % % % % % % figure;
% % % % % % % subplot(1,2,1);plot(Result(:,1)/60,Result(:,5),'k');title('Gs');ylim([0,0.3]);
% % % % % % % subplot(1,2,2);plot(Result(:,1)/60,gm*Sc*(Result(:,2)/(3 * 10^4)-Result(:,9)),'k');title('A');ylim([-1,30]);
% % % % % % % figure;
% % % % % % % plot(RuACT_CON(:,1)/60,RuACT_CON(:,5),'k.');title('RubiscoAct percent');xlabel('Time (min)');ylabel('RubiscoAct (%)')
% % % % % 
% % % % % AR1(:,1)=Result(:,1)/60-30;
% % % % % AR1(:,2)=Result(:,2);
% % % % % AR1(:,3)=Result(:,6);
% % % % % AR1(:,4)=Result(:,5);
% % % % % AR2(:,1)=Gs_VEL(:,1)/60;
% % % % % AR2(:,2)=Gs_VEL(:,2);
% % % % % AR2(:,3)=Gs_VEL(:,9);
% % % % % AR3(:,1)=RuACT_CON(:,1)/60;
% % % % % AR3(:,2)=RuACT_CON(:,5);
% % % % % % 
% % % % % % 
% % % % % % figure;
% % % % % % subplot(2,3,1); plot(AR1(:,1),AR1(:,2),'k');hold on;plot(AR13(:,1),AR13(:,2),'b');plot(AR12(:,1),AR12(:,2),'r');title('Ci');ylim([0,Air_CO2]);xlabel('Time (min)');ylabel('Ci (\mubar)')
% % % % % % subplot(2,3,2); plot(AR1(:,1),AR1(:,3),'k');hold on;plot(AR13(:,1),AR13(:,3),'b');plot(AR12(:,1),AR12(:,3),'r');title('Tleaf');xlabel('Time (min)');ylabel('Tleaf (^oC)')
% % % % % % subplot(2,3,3); plot(AR1(:,1),AR1(:,4),'k');hold on;plot(AR13(:,1),AR13(:,4),'b');plot(AR12(:,1),AR12(:,4),'r');title('Gs');xlabel('Time (min)');ylabel('gs (mol m^-^2 s^-^1)');
% % % % % % subplot(2,3,4); plot(AR2(:,1),AR2(:,2),'k');hold on;plot(AR23(:,1),AR23(:,2),'b');plot(AR22(:,1),AR22(:,2),'r');%plot(Result(:,1)/60,gm*Sc*(Result(:,2)/(3 * 10^4)-Result(:,9)),'r');
% % % % % %   title('A');ylim([0,40]);xlabel('Time (min)'); ylabel('A (\mumol m^-^2 s^-^1)')
% % % % % % subplot(2,3,5); plot(AR2(:,1),AR2(:,3),'k');hold on;plot(AR23(:,1),AR23(:,3),'b');plot(AR22(:,1),AR22(:,3),'r');%plot(Result(:,1)/60,Gsw.*(ESaturation-Eb)./(Pressure / 1000.0)*10^6.0);
% % % % % % title('Transpiration');xlabel('Time (min)');ylabel('T (\mumol m^-^2 s^-^1)')
% % % % % % subplot(2,3,6); plot(AR3(:,1),AR3(:,2),'k');hold on;plot(AR33(:,1),AR33(:,2),'b');plot(AR32(:,1),AR32(:,2),'r');title('RubiscoAct percent');xlabel('Time (min)');ylabel('RubiscoAct (%)')
% % % % % 
% % % % % 
% % % % % 
% % % % % % % R1(:,1)=Result(:,1)/60;
% % % % % % % R1(:,2)=Result(:,2);
% % % % % % % R1(:,3)=Result(:,6);
% % % % % % % R1(:,4)=Result(:,5);
% % % % % % % 
% % % % % % % R2(:,1)=Gs_VEL(:,1)/60;
% % % % % % % R2(:,2)=Gs_VEL(:,2);
% % % % % % % R2(:,3)=Gs_VEL(:,9);
% % % % % % [row2,col2]=size(Result);
% % % % % % Le2=find(Result(:,1)<2400);
% % % % % % Nu2=length(Le2);
% % % % % % Le3=find(Result(:,1)<4200);
% % % % % % Nu3=length(Le3);
% % % % % % Le4=find(Result(:,1)<3300);
% % % % % % Nu4=length(Le4);
% % % % % % Le5=find(Result(:,1)<5100);
% % % % % % Nu5=length(Le5);
% % % % % % 
% % % % % % AandW(1,1)=Result(Nu2,7);
% % % % % % AandW(2,1)=Result(Nu4,7)-Result(Nu2,7);
% % % % % % AandW(3,1)=Result(Nu3,7)-Result(Nu4,7);
% % % % % % AandW(4,1)=Result(Nu5,7)-Result(Nu3,7);
% % % % % % AandW(5,1)=Result(row2,7)-Result(Nu5,7);
% % % % % % AandW(6,1)=Result(row2,7);
% % % % % % AandW(1,2)=Result(Nu2,8);
% % % % % % AandW(2,2)=Result(Nu4,8)-Result(Nu2,8);
% % % % % % AandW(3,2)=Result(Nu3,8)-Result(Nu4,8);
% % % % % % AandW(4,2)=Result(Nu5,8)-Result(Nu3,8);
% % % % % % AandW(5,2)=Result(row2,8)-Result(Nu5,8);
% % % % % % AandW(6,2)=Result(row2,8);
% % % % % 
% % % % % 
% % % % %  [row2,col2]=size(Result);
% % % % % % % % % % % % AandW(1,i)=Result(row2,7);
% % % % % % % % % % % % AandW(2,i)=Result(row2,8);
% % % % % 
% % % % % Le2=find(Result(:,1)<3900);
% % % % % Nu2=length(Le2);
% % % % % Le3=find(Result(:,1)<6300);
% % % % % Nu3=length(Le3);
% % % % % Le4=find(Result(:,1)<8700);
% % % % % Nu4=length(Le4);
% % % % % Le5=find(Result(:,1)<11100);
% % % % % Nu5=length(Le5);
% % % % % 
% % % % % AandW(1,i)=Result(Nu3,7)-Result(Nu2,7);
% % % % % AandW(2,i)=Result(Nu4,7)-Result(Nu3,7);
% % % % % AandW(3,i)=Result(Nu5,7)-Result(Nu4,7);
% % % % % AandW(4,i)=Result(row2,7)-Result(Nu5,7);
% % % % % AandW(5,i)=Result(Nu3,8)-Result(Nu2,8);
% % % % % AandW(6,i)=Result(Nu4,8)-Result(Nu3,8);
% % % % % AandW(7,i)=Result(Nu5,8)-Result(Nu4,8);
% % % % % AandW(8,i)=Result(row2,8)-Result(Nu5,8);
% % % % % 
% % % % % AandW(9,i)=Result(Nu2+600,7)-Result(Nu2,7);
% % % % % AandW(10,i)=Result(Nu2+1200,7)-Result(Nu2+600,7);
% % % % % AandW(11,i)=Result(Nu2+1800,7)-Result(Nu2+1200,7);
% % % % % AandW(12,i)=Result(Nu2+2400,7)-Result(Nu2+1800,7);
% % % % % 
% % % % % AandW(13,i)=Result(Nu2+600,8)-Result(Nu2,8);
% % % % % AandW(14,i)=Result(Nu2+1200,8)-Result(Nu2+600,8);
% % % % % AandW(15,i)=Result(Nu2+1800,8)-Result(Nu2+1200,8);
% % % % % AandW(16,i)=Result(Nu2+2400,8)-Result(Nu2+1800,8);
% % % % % 
% % % % % AandW(17,i)=Result(Nu3+600,7)-Result(Nu3,7);
% % % % % AandW(18,i)=Result(Nu3+1200,7)-Result(Nu3+600,7);
% % % % % AandW(19,i)=Result(Nu3+1800,7)-Result(Nu3+1200,7);
% % % % % AandW(20,i)=Result(Nu3+2400,7)-Result(Nu3+1800,7);
% % % % % 
% % % % % AandW(21,i)=Result(Nu3+600,8)-Result(Nu3,8);
% % % % % AandW(22,i)=Result(Nu3+1200,8)-Result(Nu3+600,8);
% % % % % AandW(23,i)=Result(Nu3+1800,8)-Result(Nu3+1200,8);
% % % % % AandW(24,i)=Result(Nu3+2400,8)-Result(Nu3+1800,8);
% % % % % 
% % % % % AandW(25,i)=Result(Nu4+600,7)-Result(Nu4,7);
% % % % % AandW(26,i)=Result(Nu4+1200,7)-Result(Nu4+600,7);
% % % % % AandW(27,i)=Result(Nu4+1800,7)-Result(Nu4+1200,7);
% % % % % AandW(28,i)=Result(Nu4+2400,7)-Result(Nu4+1800,7);
% % % % % 
% % % % % AandW(29,i)=Result(Nu4+600,8)-Result(Nu4,8);
% % % % % AandW(30,i)=Result(Nu4+1200,8)-Result(Nu4+600,8);
% % % % % AandW(31,i)=Result(Nu4+1800,8)-Result(Nu4+1200,8);
% % % % % AandW(32,i)=Result(Nu4+2400,8)-Result(Nu4+1800,8);
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % AandW=AandW';
% 
% [row,col]=size(RuACT_CON);
% Leng=find(RuACT_CON(:,1)>=1800);
% Num=length(Leng);
% At=zeros(Num,2);
% Numn=row-Num+1;
% At(:,1)=(RuACT_CON(Numn:row,1)-1800)/60;
% At(:,2)=RuACT_CON(Numn:row,5);

% figure;
% plot(At(:,1),At(:,2));

% figure;
% subplot(2,3,1); plot(Result1(:,1)/60,Result1(:,2),'k');hold on;plot(Result2(:,1)/60,Result2(:,2),'b');plot(Result(:,1)/60,Result(:,2),'r');title('Ci');ylim([0,Air_CO2]);
% subplot(2,3,2); plot(Result1(:,1)/60,Result1(:,6),'k');hold on;plot(Result2(:,1)/60,Result2(:,6),'b');plot(Result(:,1)/60,Result(:,6),'r');title('Tleaf');
% subplot(2,3,3); plot(Result1(:,1)/60,Result1(:,5),'k');hold on;plot(Result2(:,1)/60,Result2(:,5),'b');plot(Result(:,1)/60,Result(:,5),'r');title('Gs');
% subplot(2,3,4); plot(Gs_VEL1(:,1)/60,Gs_VEL1(:,2),'k.');hold on;plot(Gs_VEL2(:,1)/60,Gs_VEL2(:,2),'b.');plot(Gs_VEL(:,1)/60,Gs_VEL(:,2),'r.');title('A');ylim([0,40])
% subplot(2,3,5); plot(Gs_VEL1(:,1)/60,Gs_VEL1(:,9),'k.');hold on; plot(Gs_VEL2(:,1)/60,Gs_VEL2(:,9),'b.');plot(Gs_VEL(:,1)/60,Gs_VEL(:,9),'r.');title('Transpiration');
% subplot(2,3,6); plot(RuACT_CON1(:,1)/60,RuACT_CON1(:,5),'k.');hold on;plot(RuACT_CON2(:,1)/60,RuACT_CON2(:,5),'b.');plot(RuACT_CON(:,1)/60,RuACT_CON(:,5),'r.');title('Actived Rubisco');

% figure;
% plot(Result(:,1)/60,Result(:,18),'k')
% figure;
% subplot(2,3,1); plot(RuACT_VEL(:,1),RuACT_VEL(:,2),'k.');title('v1');
% subplot(2,3,2); plot(RuACT_VEL(:,1),RuACT_VEL(:,3),'k.');title('vn1');
% subplot(2,3,3); plot(RuACT_VEL(:,1),RuACT_VEL(:,4),'k.');title('v7');
% subplot(2,3,4); plot(RuACT_VEL(:,1),RuACT_VEL(:,5),'k.');title('vn7');
% subplot(2,3,5); plot(RuACT_VEL(:,1),RuACT_VEL(:,6),'k.');title('v6_1');
% subplot(2,3,6); plot(RuACT_VEL(:,1),RuACT_VEL(:,7),'k.');title('v6_2');
% 
% figure;
% subplot(2,2,1); plot(Result(:,1)/60,Result(:,45),'k');title('ER');
% subplot(2,2,2); plot(Result(:,1)/60,Result(:,46),'k');title('EAF');
% subplot(2,2,3); plot(Result(:,1)/60,Result(:,47),'k');title('ECMR');
% subplot(2,2,4); plot(Result(:,1)/60,Result(:,48),'k');title('RuBP');
% figure;
% subplot(2,2,1); plot(Result(:,1)/60,Result(:,45)+Result(:,46)+Result(:,47),'k');title('Et');
% figure;
% subplot(2,3,1); plot(Result(:,1)/60,Result(:,2),'k');hold on; plot(Result0(:,1)/60,Result0(:,2),'r');title('Ci');ylim([0,Air_CO2]);
% subplot(2,3,2); plot(Result(:,1)/60,Result(:,6),'k');hold on; plot(Result0(:,1)/60,Result0(:,6),'r');title('Tleaf');
% subplot(2,3,3); plot(Result(:,1)/60,Result(:,5),'k');hold on; plot(Result0(:,1)/60,Result0(:,5),'r');title('Gs');
% subplot(2,3,4); plot(Gs_VEL(:,1)/60,Gs_VEL(:,2),'k.');hold on;plot(Gs_VEL0(:,1)/60,Gs_VEL0(:,2),'r.');title('A');ylim([0,40])
% subplot(2,3,5); plot(Gs_VEL(:,1)/60,Gs_VEL(:,9),'k.');hold on;plot(Gs_VEL0(:,1)/60,Gs_VEL0(:,9),'r.');title('Transpiration');

[row2,col2]=size(Result);
AandW(i,1)=Result(row2,7)';
AandW(i,2)=Result(row2,8)';
%%%%%%%%%%%%%%%%%%%%%%%
%   output  step. Notice if the graph needs to be displayed, then, the following line should not be commented out     %
%%%%%%%%%%%%%%%%%%%%%%%
%AMeM=CM_OUT(Tt,d);

PSPR_SUCS_com = 0;
IModelCom;
end
end
