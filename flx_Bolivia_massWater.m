%Calc_equil_water_AH 
function mpi=flx_Bolivia_massWater(Masses,time)

global Temp R  whats
global p_sat psat_w_vec sigma
global pvi_vec tims
global M_vec rhoi  %sigmai NumbConc
global mass_water p_event

timein=time/24/60/60+p_event;

    M=interp1(tims,M_vec,timein);
    pvi=interp1(tims,pvi_vec,timein);
    P_sats=interp1(tims,p_sat,timein);
    psat_w=interp1(tims,psat_w_vec,timein);
    T=interp1(tims,Temp,timein);
    
    if any(isnan(M))
        M(isnan(M))=interp1(tims,M_vec,timein(isnan(M(:,9))),'pchip');
        pvi(isnan(pvi))=interp1(tims,pvi_vec,timein(isnan(pvi(:,end))),'pchip');
        P_sats(isnan(P_sats))=interp1(tims,p_sat,timein(isnan(P_sats(:,end))),'pchip');
        psat_w(isnan(psat_w))=interp1(tims,psat_w_vec,timein(isnan(psat_w)),'pchip');
        T(isnan(T))=interp1(tims,Temp,timein(isnan(T)),'pchip');
    end
organs=find(whats==-2 | whats==-1 | whats==1 | whats==2 | whats>=10);
ite_dif=0.02;
psat_H2SO4=0;
psat_NH3=3.32e-6;
psat=zeros(size(M));
psat(:,whats==0)=psat_w;
psat(:,whats==-3)=psat_H2SO4;
psat(:,whats==3)=psat_NH3;
psat(:,organs)=P_sats*101325;

%Masses of each component
mpi=zeros(size(Masses,1),length(whats));
mpi(:,whats==0)=mass_water;
mpi(:,whats<0 | whats>=10 | whats==3)=Masses;

moles=mpi./M;

activ=1;
 Xmass=mpi(1,:)/sum(mpi(1,:));
 mp=sum(mpi(1,:));
% Particle density
rho = ( sum(Xmass./rhoi) ).^(-1);
       
% Particle radius        
rp = (3*(mp/rho)/(4*pi)).^(1.0/3.0);

% Kelvin's effect
Ke = exp(2.*M(1,:).*sigma./R./T(1)./rp./rho);

for ii=1:size(Masses,1)
    still_not_close_enough=1;
    loop_nr=0;
    while still_not_close_enough==1
        loop_nr=loop_nr+1;
        water_moles=(pvi(ii,whats==0)*sum(moles(ii,whats~=0)))/(activ*psat(ii,whats==0)*Ke(whats==0)-pvi(ii,whats==0));

        old_water_moles=moles(ii,whats==0);
        moles(ii,whats==0)=water_moles;

        if ii==size(Masses,1)
                mpi(ii,whats==0)=moles(ii,whats==0).*M(ii,whats==0);
            else
                mpi(ii:ii+1,whats==0)=moles(ii,whats==0).*M(ii,whats==0);
            end
        Xmass=mpi(ii,:)/sum(mpi(ii,:));

        % Particle density
        rho = ( sum(Xmass./rhoi) )^(-1);

        mp=sum(mpi(ii,:));
        % Particle radius        
        rp = (3*(mp/rho)/(4*pi)).^(1.0/3.0);

        % Kelvin's effect
        Ke = exp(2.*M(ii,:).*sigma./R./T(ii)./rp./rho);

        water_dif=abs((water_moles-old_water_moles)./water_moles);

        if water_dif<=ite_dif
            still_not_close_enough=0;
        end
    end
end
 