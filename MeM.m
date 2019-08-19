%function AMeM = MeM(inputfile,PCO2,Plight,PTemp,PGRNC,PGRNT)
% MeM('MeM_input.txt',360,2000,25,0,0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ki_Gs;
global kd_Gs;
global GsResponse;
ki_Gs=0.9*60;%2*60;
kd_Gs=4.1*60;%1*60;
GsResponse=1; %%if GsResponse=0 Ball berry model no time dependent Gs response ; %%if GsResponse=1 time dependent Gs response, using ki_Gs and kd_Gs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Para_mata;
global PhotosynthesisType;
Para_mata=1;%%if Para_mata=1, C3 Metabolic model and Gs model integrated  if Para_mata=0 Farquhar mdoel and gs model
PhotosynthesisType=1;% 1:C3 model 2:C4 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leaf model parameters

global Air_CO2;
global Air_O2;
global WeatherRH;
global WeatherWind;
global Radiation_PAR;
global Vcmax25;
global Jmax25;
global WeatherTemperature;

global Radiation_NIR;
global Radiation_LW;
global PhiLeaf

Convert=1E6/(2.35E5); 
WeatherTemperature=25;
Air_CO2=400;
Air_O2=210.0;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=2000/Convert;%10*i;
Radiation_NIR=0;
Radiation_LW=0;
Vcmax25=100;
Jmax25=200;
PhiLeaf=-0;%Mpa Min=-2.3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global O2_cond
global GRNC;
global GRNT;
global Vfactor;

PCO2=Air_CO2;
O2_cond=Air_O2/1000*1.26;
Plight=Radiation_PAR*Convert;
PTemp=WeatherTemperature;
PGRNC=0;
PGRNT=0;
indata1= importdata('MeM_input.txt');
%indata2= importdata('../Input/MeM_static.txt');
GRN_data=indata1.data; 
%Env_data=indata2.data; 
Vfactor=GRN_data(:,10);%Vfactor=GRN_data(:,10);
GRNC=PGRNC;
GRNT=PGRNT;
global CI
CI=Air_CO2*0.4/(3 * 10^4);
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
Jmax=0.180;
alfa=0.85;
fc=0.15;
Theta=0.7;
beta=0.7519;



global Tp;
Tp=PTemp;
Begin = 1;
global tglobal;     % The total running time
tglobal = 1800;
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
%%%%%%%%%%%%%%%%%%%%%%%%
%   Initialation step %
%%%%%%%%%%%%%%%%%%%%%%%%

global GP; 
GP = 0;
% Begin = 1;
% PS_PRs = zeros(23,1);
%CMs = CM_Ini(Begin);
CMs = C3leafMetaIni(Begin);

ModelComb = IModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

global PR_PS_com;             % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 1;

global ATPActive;
ATPActive = 0;

global PSPR_SUCS_com;        % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step %
%%%%%%%%%%%%%%%%%%%%%%%%

CM_Param = 0;

SUCS_Param = zeros(2,1);
SUCS_Param(1) = 1;
SUCS_Param(2) = 1;

PS_PR_Param = 0;


[Tt,d] = ode15s(@CM_mb,[0,time],CMs, options1,PS_PR_Param, SUCS_Param);   
%[Tt,d] = ode23tb(@CM_mb,[0,time],CMs, options1,PS_PR_Param, SUCS_Param);   
global d_plot;
d_plot=d;

global Tt_plot;
Tt_plot = Tt;
global Result;

Result =[Tt,d]; 
% Eb=d(:,3);
% vinf=gm*Sc*10^(-3)*(Ci/(3 * 10^4)-MC_CO2);
% vH2Ob=Gbw*(Eb-Ea)/(Pressure / 1000.0)*10^6.0;
%global Gs_VEL;
figure;
subplot(2,3,1); plot(Result(:,1)/60,Result(:,2),'k');title('Ci');ylim([0,Air_CO2]);
subplot(2,3,2); plot(Result(:,1)/60,Result(:,6),'k');title('Tleaf');
subplot(2,3,3); plot(Result(:,1)/60,Result(:,5),'k');title('Gs');
subplot(2,3,4); plot(Gs_VEL(:,1)/60,Gs_VEL(:,2),'k.');title('A');ylim([0,40])
subplot(2,3,5); plot(Gs_VEL(:,1)/60,Gs_VEL(:,9),'k.');title('Transpiration');
figure;
plot(Result(:,1)/60,Result(:,18),'k')
% figure;
% subplot(2,3,1); plot(Result(:,1)/60,Result(:,2),'k');hold on; plot(Result0(:,1)/60,Result0(:,2),'r');title('Ci');ylim([0,Air_CO2]);
% subplot(2,3,2); plot(Result(:,1)/60,Result(:,6),'k');hold on; plot(Result0(:,1)/60,Result0(:,6),'r');title('Tleaf');
% subplot(2,3,3); plot(Result(:,1)/60,Result(:,5),'k');hold on; plot(Result0(:,1)/60,Result0(:,5),'r');title('Gs');
% subplot(2,3,4); plot(Gs_VEL(:,1)/60,Gs_VEL(:,2),'k.');hold on;plot(Gs_VEL0(:,1)/60,Gs_VEL0(:,2),'r.');title('A');ylim([0,40])
% subplot(2,3,5); plot(Gs_VEL(:,1)/60,Gs_VEL(:,9),'k.');hold on;plot(Gs_VEL0(:,1)/60,Gs_VEL0(:,9),'r.');title('Transpiration');


%%%%%%%%%%%%%%%%%%%%%%%
%   output  step. Notice if the graph needs to be displayed, then, the following line should not be commented out     %
%%%%%%%%%%%%%%%%%%%%%%%
%AMeM=CM_OUT(Tt,d);

PSPR_SUCS_com = 0;
IModelCom;
