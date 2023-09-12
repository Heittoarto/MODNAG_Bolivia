%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for condensation model MODNAG                                                    %
% (Model for oligomerization and decomposition in nanoparticle growth)      %
%                                                                                                                   %
% By Arto Yli-Heitto, University of Eastern Finland, September 2020                %
% arto.heitto@uef.fi                                                                                        %
%                                                                                                                    %
% This routine calculates condensation of vapors on small monodisperse        %
% aerosol particles  and oligomerization and decomposition.                          %
%                                                                                                                    %
%                                                                                                                   %
% Vapor concentrations, RH and temperature are assumed to be constant     %
% during the simulation                                                                                 %
%                                                                                                                   %
% Modifications:                                                                                             %
%  Modified to work with Bolivia data AH 2021 
% Oligomerization and decompostiton excluded from the code AH 2021 %
%-------------------------------------------------------------------------             %
%                                                                                                                    %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version is for a system with 9 condensing compounds.
% Oligomerization and decomposition are omitted                 %
%                                                                                                                %
% The order of compounds for calculations is:                                              %
%   1 = water                                                                                               %
%   2 = ammonia                                                                                         %
%   3 = sulfuric acid                                                                                      %
%   4-9 = organics by VBS 1 - -4, no reactant                                              %
%                                                                    %
%                                                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input files that are needed are:
% -Time vector of measurement times in matlab time  format
% -Concentrations of organic compounds (#/cm^3) in matrix with one column for each model compound
%     and rows equal to measurement times
% -Molar masses of organic models compounds (g/mol) in similar matrix as concentrations
% -Diffusion coefficients of organic models compounds in similar matrix as concentrations
% -Saturation vapor pressures of organic models compounds (mug/m^3) in similar matrix as concentrations
% -Measured sulfuric acid concentrations (#/cm^3) in a vector with similar length as time vector 
% -Measured ambient temperature (K) in a vector with similar length as time vector 
% -Measured ambient RH in a vector with similar length as time vector 
% -Measured ambient pressure (Pa) in a vector with similar length as time vector 
% -Starting times of the modelled NPF events in matlab time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global kb R Temp Press Avogado M_vec sigmai rhoi pvi_vec p_sat whats psat_w_vec sigma tims
global mass_water alpham mass_NH3 p_event D_vec

whats=[0 3 -3 10 10 10 10 10 10]; % define the model compounds for MODNAG simulation.
% 0=water
% 1=amine
% -1=organic monoacid
% -2=organic diacid
% -3=H2SO4
% 3=ammonia


%Read the inputs
tims=load('Times.dat'); %Time vector
org_conc=load('Concs.dat'); %Concnetrations of organic model compounds 
org_M=load('Mmass.dat');% Molar masses of organic model compounds
Temp=load('Temp.dat'); %Ambient temperature
p_sat=load('pL.dat'); %Saturation vapor pressures of organic model compounds
Press=load('press.dat'); %Ambient pressure
H2SO4_conc=load('H2SO4.dat'); % Sulfuric acid concentration
Ammoconc=H2SO4_conc;
RH_vec=load('RH.dat');
Diffu=load('Diffu.dat');

org_M(org_conc==0)=1;
Diffu(Diffu==0)=1;

concs=[RH_vec Ammoconc H2SO4_conc org_conc];%

acids=find(whats<0); %Indexes for acids (H2SO4 and organic acids) in eaim
bases=find(whats>0 & whats<10); %indexes for bases (ammonia and amine) in eaim
neutrals=find(whats>=10); %Indexes for neutrals in eaim

% starting times of the NPF events
events=load('Events.dat');
p_event=events(1);


%Gas constant (J K^-1 mol^-1)
R = 8.314472;
%Boltzmann constant (m^2 kg s^-2 K^-1)
kb = 1.3806503e-23;
%Avogadro's constant (mol^-1)
Avogado = 6.0221415e+23;


%Gas phase molecular concentrations
%Acids (Only H2SO4, molecules/m3) 
c_a=concs(:,acids)*1e6;
%Bases (Only Ammonia, molecules/m3)
c_b=concs(:,bases)*1e6;
%Neutrals (Organics, molecules/m3)
c_n=concs(:,neutrals)*1e6;
    

% Molar masses (kg/mol)
M_vec=zeros(size(org_M,1),length(whats));
M_vec(:,whats==0)=0.018;%Water
M_vec(:,whats==-3)=0.098; %H2SO4
M_vec(:,whats==-4)=0.063; %HNO3
M_vec(:,whats==3)=0.017; %NH3
M_vec(:,whats~=0 & whats~=-3 & whats~=-4 & whats~=3)=org_M*1e-3; %Organic acids and neutrals
                        
% Densities of the compounds (kg/m^3)
rhoi=1200*ones(1,length(whats));

%Surface tension (N/m)
sigmai=0.03*ones(1,length(whats));

%Diffusion coefficients (m2/s)
D_vec=zeros(size(org_M,1),length(whats));
D_vec(:,whats==0)=21.485e-6; %water
D_vec(:,whats==-3)=9.3803e-6; %H2SO4
D_vec(:,whats==-4)=1.18e-5; %HNO3
D_vec(:,whats==3)=19e-6; %NH3
D_vec(:,whats~=0 & whats~=-3 & whats~=-4 & whats~=3)=Diffu;

% Mass accommodation coefficients 
alpham=1.0*ones(1,length(whats));

%Converting concentrations (molecules/m^3) to pressures (Pa)
p_a = c_a.*kb.*repmat(Temp,1,size(c_a,2));
p_b = c_b.*kb.*repmat(Temp,1,size(c_b,2));
p_n = c_n.*kb.*repmat(Temp,1,size(c_n,2));

%Calculating ambient partial pressure of water
%Coefficients from 1-component model file models.gro
A_w = [77.34491296e+0; 7235.424651e+0; 0.82e+1; 0.0057113e+0; 0.e+0];
psat_w_vec = exp(A_w(1) - A_w(2)./Temp - A_w(3).*log(Temp) + A_w(4).*Temp + A_w(5).*Temp.^2);
p_w = RH_vec/100 .* psat_w_vec;

%Ambient partial pressures
pvi_vec=zeros(size(org_M,1),length(whats));
pvi_vec(:,whats==0)=p_w;
pvi_vec(:,whats>0 & whats<10)=p_b;
pvi_vec(:,whats>=10)=p_n;
pvi_vec(:,whats<0)=p_a;
    
% Initial composition 
Ni=[0 0 40 0 zeros(1,size(org_conc,2)-1)];

% Masses of each component in a particle
mpi = Ni.*M_vec(1,:) ./Avogado;

% Mass of a particle
mp = nansum(mpi);

% Mass fractions
Xmass = mpi./mp;

% Mole fractions
Xmole = Ni./nansum(Ni);

% Density of the particle (kg/m^3)
rho = (nansum(Xmass./rhoi) )^(-1);

%Initial radius of the particle
rp0 = (3*(mp/rho)/(4*pi)).^(1.0/3.0);


%Surface tension
sigma = nansum(Xmole.*sigmai);

% Inputs: only variables are the masses in the particulate phase
input = mpi;


%% Make the time vector for simulated data 
spacing_data= [-3 101 3601
                2 3600 99400
                500 3000 10000
                1 2 2];
starttimes=[];
for ii=1:size(spacing_data,2)
    if spacing_data(4,ii)==1 %number one here means, that times are logarhytmically spaced during this section
        starttimes=[starttimes logspace(spacing_data(1,ii), spacing_data(2,ii),spacing_data(3,ii))];
        if ii>1 && spacing_data(4,ii)==1
            starttimes(end-spacing_data(3,ii)+1)=[];
        end
    elseif spacing_data(4,ii)==2%number two here means, that times are linearly spaced during this section
        starttimes=[starttimes linspace(spacing_data(1,ii), spacing_data(2,ii),spacing_data(3,ii))];      
    end
end

%%

        mass_water=input(whats==0);
        mass_NH3=input(whats==3);
        
        
      options = odeset('RelTol',1E-6,'AbsTol',1E-33);
     [tout, output, exitflag1] = ode15s(@flx_MODNAG_Bolivia,starttimes,input(whats<0 | whats>=10),options,whats);
    
    output=[output(:,1)/0.098*0.017 output];
        output = flx_Bolivia_massWater(output,starttimes);
         
%%  
Xmass = output./repmat(sum(output,2),1,length(whats));
rho = ( sum(Xmass./repmat(rhoi,length(Xmass),1),2) ).^(-1);
rp =( 3/(4*pi)*sum(output,2)./rho).^(1/3);

model_result = [p_event+tout/60/60/24 rp Xmass output(:,1:length(whats))];
