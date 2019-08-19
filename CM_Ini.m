%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright ? 2007

%   This file is part of CarbonMetabolism.

%    CarbonMetabolism is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.

%    CarbonMetabolism is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    

% This function initialize different parameters used in the model CM

function CMs = CM_Ini(Begin)
global Tp;


global PS_C_CA;             %   Global constant for the total adenylates
global PS_C_CP;             %   Global constant for the total phosphate
global PS_C_CN;             %   Global constant for the total NADP+NADPH

PS_C_CP	= 25;   
PS_C_CA	=1.5;          
PS_C_CN	=1;

global NADHc;
global NADc;
global GLUc;
global KGc;

NADHc = 0.47;        
NADc = 0.4;        
GLUc=24;
KGc=0.4;

global SUCS_Pool;	
ATc =  1;      % mM
UTc =  1.5;    % mM
PTc = 15;      % mM  
	
SUCS_Pool	(	1	)	=	ATc;
SUCS_Pool	(	2	)	=	UTc;
SUCS_Pool	(	3	)	=	PTc;

RuBP	= 2.000;
PGA	    =2.400;
DPGA	=0.0011;
T3P	    =0.5;
NONE	=0;
FBP	    =0.670;
E4P	    =0.050;
S7P	    =2.000;
SBP	    =0.300;
ATP	    =0.68;         
NADPH	=0.21;         
HexP    = 2.2;
PenP    = 0.25;
CO2	    = 0.012;        
O2	    = 0.264; 

SERc= 7.5;                  % Serine in cytosol; 7.5 original value
GLYc = 1.8;                 % Glycine in cytosol; 1.8 original vlaue
%PGA = 4.3;                  % PGA in chloroplast;4.3 is the original value;
GOAc = 0.028;              % Glyoxylate in cytosol; 0.028; EXPERIMENTAL DATA;
GCAc = 0.36;                   % See the note for GCA.
GCA = 0.36;                    % Derived from radioactive labelling experiment; assuem equal concenatration 
                               % inside and outshide chloroplast
PGCA= 0.0029;                % Phosphoglycolate in chloroplast derived based on the Km112; orignal value is : 0.0029; 
GCEA =0.1812;                  % Glycerate in chloroplast; derived based on V113
GCEAc = 0.1812;                 % Glycerate in cytosol; assume at equilibrium with GCEA initially.
HPRc = 0.0035;                % HydroxylPyruvate; derived from equation 123;

T3Pc	=	2.3	;
FBPc	=	2	;
HexPc	=	5.8	;
F26BPc	=	7.8*10^(-6)	;
ATPc	=	0.35	;          
ADPc	=	0.65	;
UDPGc	=	0.57	;
UTPc	=	0.75	;       
SUCP	=	0;
SUC	=		0;
PGAc =  0; 


CMs = zeros(35,1);
CMs(1)	=RuBP;
CMs(2)	=PGA;
CMs(3)	=DPGA;
CMs(4)	=T3P;
CMs(5)	=FBP;
CMs(6)	=E4P;
CMs(7)	=S7P;
CMs(8)	=SBP;
CMs(9)	=ATP;
CMs(10)	=NADPH;
CMs(11)	= CO2;
CMs(12)	= O2;
CMs(13) = HexP;
CMs(14) = PenP;
CMs(15) = GCEA;
CMs(16) = GCA;
CMs(17) = PGCA;
CMs(18) = GCAc;
CMs(19) = GOAc;
CMs(20) = SERc;
CMs(21) = GLYc;
CMs(22) = HPRc;
CMs(23) = GCEAc;
CMs(24)= T3Pc	;
CMs(25)= FBPc	;
CMs(26)= HexPc	;
CMs(27)= F26BPc	;
CMs(28)= ATPc	;
CMs(29)= ADPc	;
CMs(30)= 0	;
CMs(31)= UDPGc	;
CMs(32)= UTPc	;
CMs(33)= SUCP	;
CMs(34)= SUC	;
CMs(35)= PGAc	;
       

% Initialize the Vmax for different reactions

global	V1	;
global	V2	;
global	V3	;
global	V5	;
global	V6	;
global	V7	;
global	V8	;
global	V9	;
global	V10	;
global	V13	;
global	V16	;
global	V23	;
global	V31	;
global	V32	;
global	V33	;

global GP; 
global Vcmax25;
fsoy=Vcmax25*0.95/30/2.913930914	;%WY201811
if GP == 0
	V1	=	2.913930914*fsoy/33;
	V2	=	30.1408176*fsoy/33;
	V3	=	4.039482839*fsoy/33	;
	V5	=	1.218890457*fsoy/33	;
	V6	=	0.726259575*fsoy/33;
	V7	=	3.122215944*fsoy/33	;
	V8	=	1.218890364	*fsoy/33;
	V9	=	0.324190582*fsoy/33	;
	V10	=	3.122068837*fsoy/33	;
	V13	=	10.83475664*fsoy/33	;
	V23	=	0.26684349/33;
    V16 = 5/33; 
end
V31		=	3.73/3/33;   % 1.05 *SC  *1.0 ;	%	(Lilley, Chon, Mosbach & Heldt, 1977b)	31	Phosphate translocator	DHAPi<->DHAPo   1.05 defulat
V32		=	3.73/3/33;   %1.05 *SC *1.0;	    %	(Lilley et al., 1977b)	32	Phosphate translocator	PGAi<->PGAo 1.05 default
V33		=	3.73/3/33;   %1.05 *SC * 1.0;	    %	(Lilley et al., 1977b)	33	Phosphate translocator	GAPi<->GAPo 1.05 default
%%%%%%%%%%%%%soy


global V31_ps2ca;
global V32_ps2ca;
global V33_ps2ca;
V31_ps2ca = V31;
V32_ps2ca = V32;
V33_ps2ca = V33;

% To set global information for different reactions

% Reaction: 111: RUBP+O2<-->PGlycolate + PGA
global V111;
global gp2V111;
V111 = gp2V111;
% Reaction: 112: PGlycolate-->Pi+Glycolate;
global V112;       
% Reaction 113  : Gcea+ATP<-->ADP + PGA
global V113;
% Reactoin 121; Glycolate +O2<-->H2O2+Glyoxylate
global V121;
% Reaction 122  : Glyoxylate + Serine<--> Hydoxypyruvate + Glycine;
global V122;
% Reaction 123: HydroxylPyruvate + NAD <--> NADH + Glycerate
global V123;
% Reaction 124: Glyoxylate + Glu  <--> KG + Glycine;
global V124;
% Reaction 131: NAD+Glycine <--> CO2+ NADH + NH3
global V131;
% The consant for calculating the glycerate uptake.
global V1T;
global V2T;

if GP == 0
    V111 = V1 * 0.22/33; 
	V112	=	52.41992121/33	;
	V113	=	5.715787563/33;
	V121	=	1.456108923/33;
	V122	=	3.306190845/33;
	V123	=	10.00978112/33;
	V124	=	2.745819515/33;
	V131	=	2.494745448/33;
end
V1T = 5/33;
% The constant for calculating the glycolate uptake
V2T = 6/33;      % The original value is 0.32.


% The following calculate the total concentration of different enzymes. 

global 	V51	;
global 	V52	;
global 	V55	;
global 	V56	;
global 	V57	;
global 	V58	;
global 	V59	;
global  V60;
global  V61;
global	V62;	


% Unit: mmol l-1 s-1;
if GP == 0
	V51	=	0.107376831/33;
	V52	=	0.063979048/33;
	V55	=	0.115403205/33;
	V56	=	0.055503446/33;
	V57	=	0.55503446/33;
	V58	=	0.016819226/33;
end

V59	=	V58 * 6/33	;	%	F6P + ATP --ADP + F26BP
V60	=	6.1/33;	%	ATP+UDP --UTP + ADP
V61	=	10000/33;	%	POPO --2PO
V62	=	2/33;	%	SUC Sink        0.2 works.

end