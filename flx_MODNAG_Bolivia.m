function flx = flx_MODNAG_Bolivia(time, input,whats)

%This function is for calculating the mass flux according to equations
%where particle diffusion and vapor molecule dimensions are taken into
%account according to Lehtinen & Kulmala 2003 (ACP) and Nieminen et al.
%2010 (ACP).
global Temp R kb Avogado Press p_event
global p_sat psat_w_vec sigma
global pvi_vec tims
global M_vec rhoi D_vec alpham %sigmai NumbConc
global mass_water 
timein=time/24/60/60+p_event;

if timein<tims(end)
M=interp1(tims,M_vec,timein);
    D=interp1(tims,D_vec,timein);
    pvi=interp1(tims,pvi_vec,timein);
    P_sats=interp1(tims,p_sat,timein);
    psat_w=interp1(tims,psat_w_vec,timein);
    T=interp1(tims,Temp,timein);
    press=interp1(tims,Press,timein);
else
M=interp1(tims,M_vec,timein,'pchip');
    D=interp1(tims,D_vec,timein,'pchip');
    pvi=interp1(tims,pvi_vec,timein,'pchip');
    P_sats=interp1(tims,p_sat,timein,'pchip');
    psat_w=interp1(tims,psat_w_vec,timein,'pchip');
    T=interp1(tims,Temp,timein,'pchip');
    press=interp1(tims,Press,timein,'pchip');
end


organs=find(whats==-2 | whats==-1 | whats==1 | whats==2 | whats>=10);
ite_dif=0.02;
psat_H2SO4=0;
psat_NH3=3.32e-6;
psat=zeros(size(whats));
psat(whats==0)=psat_w;
psat(whats==-3)=psat_H2SO4;
psat(whats==3)=psat_NH3;
psat(organs)=P_sats*101325;

%Masses of each component
mpi=zeros(1,length(whats));
mpi(whats==0)=mass_water;
mpi(whats==3)=input(1)/M(whats==-3)*M(whats==3);
mpi(whats<0 | whats>=10)=input;

mpi(mpi<0)=0;
moles=mpi./M;

 %Mass fractions
    Xmass=mpi/sum(mpi);
    
% Particle density
rho = ( sum(Xmass./rhoi) )^(-1);

 %Particle radius
    rp = (3*(sum(mpi)/rho)/(4*pi)).^(1.0/3.0);
    
act=1;
Ke = exp(2.*M.*sigma./R./T./rp./rho);

still_not_close_enough=1;
loop_nr=0;

while still_not_close_enough==1
    loop_nr=loop_nr+1;
    water_moles=(pvi(whats==0)*sum(moles(whats~=0)))/(act*psat(whats==0)*Ke(whats==0)-pvi(whats==0));
    
    old_water_moles=moles(whats==0);
    moles(whats==0)=water_moles;
    
    mpi=moles.*M;
    
    %Mass fractions
    Xmass=mpi/sum(mpi);
    
    % Particle density
    rho = ( sum(Xmass./rhoi) )^(-1);
    
    %Particle radius
    rp = (3*(sum(mpi)/rho)/(4*pi)).^(1.0/3.0);

    % Kelvin's effect
    Ke = exp(2.*M.*sigma./R./T./rp./rho);
    
    water_dif=abs((water_moles-old_water_moles)./water_moles);
    
    if water_dif<=ite_dif
        still_not_close_enough=0;
    end
end

%Particle mass
mp = sum(mpi);

%Mole fractions
Xmole=moles/sum(moles);

p_eq=(act.*Xmole.*psat.*Ke);

dp = 2.0*rp;

%--- Transition regime corrections ----------------------------------------
%According to Lehtinen & Kulmala and Neiminen et al.

%Calculating diffusion coefficient of particle
Ts=120; %Sutherland's coefficient for water (Licht and Stechert, 1944)
T0=291.15; %Reference temperaure (K)
visc_air0=18.27e-6; %Gas viscosity of air in reference temperature (kg/m/s)
visc_air=visc_air0*(T0+Ts)/(T+Ts)*(T/T0)^(3/2); %Gas viscosity of air (kg/m/s)
% visc_air=1.57e-5; %Gas viscosity of air (kg/m/s)
M_air = 29.0e-3;  %Molar mass of air (kg/mol)
lambda_air = 2*visc_air/(press*(8*M_air/(pi*R*T))^(1/2));  %Mean free bath of air molecules
Kn_air = 2*lambda_air/dp;  %Knudsen number using air mean free bath
Cc = 1+Kn_air*(1.257+0.4*exp(-1.1/Kn_air)); %Cunningham slip correction factor
Diffp = kb*T*Cc/(3*pi*visc_air*dp);  %Diffusion coefficient of particle

%Mass of vapor molecules
mv = M./Avogado;

%Mean thermal speed
cv = (8*kb*T./(pi.*mv)).^(1/2);
cp = (8*kb*T/(pi*mp))^(1/2);

%mean free bath
lambda = 3.0.*(Diffp+D)./(cp.^2+cv.^2).^(1/2);  

%Diametr of vapor molecules
dv = (6.*mv./(pi.*rhoi)).^(1/3);

%Knudsen number
Kn = 2.*lambda./(dp+dv); 

%Transition regime correction factor for mass flux
beta = (1.0 + Kn)./(1.0 + (4./(3.0.*alpham) + 0.377).*Kn + 4.0./(3.0.*alpham).*Kn.^2.);


%--- Mass fluxes ----------------------------------------------------------

%Mass fluxes
p_eq(p_eq<1e-15)=0;
flx_mass = M .* 2*pi.*(dv + dp).*(D + Diffp).*beta.*(pvi-p_eq)./(R*T);
flx_mass(mpi <= 0.0 & (pvi-p_eq) < 0)=0;

% Mass fluxes are calculated only for acids and neutrals (bases and water are assumed to
% be in equilibrium)
flx = flx_mass(whats<0 | whats>=10)';
%--------------------------------------------------------------------------
end