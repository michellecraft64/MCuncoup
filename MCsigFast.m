function [endCond]=MCsigFast(Iapp_M,sigm,ReqNums,icY0,flName)
%implementing 1 M.C. based on Li & Cleland, uses sum of soma/dendrite
%values for all currents, similar to Michelle's but with gDR_M=15
% !!! t is in ms  !!! return spike times (spTimes) after transient time
% meant for SHORT times to show spikes AND voltage (not 200,000ms or anything)
% using RK4 time-stepping
% Only saving total number of spikes
% similar to MCfast.m BUT adding white noise in voltage
% UNITS of sigm = [micro-Amps/sqrt(ms)]

rng('shuffle'); %seed random # gen

dt=0.05; %assuming equally spaced, in ms

F=9.64853329e4; %Coulombs per Mole
T=273.15+35; %Kelvin
RTovF=8314.472*T/F; %R*1000 so in mV
thick_M=1; %thickeness of membrane shell, 1micron for Mitral Cell

% PARAMS
C_MP=1.2; %micro-F/cm^2
Rm_GM=30; %kilo-Ohm/cm^2
%maximal conduct
gNa_M=120; %mS/cm^2
gNaP=0.42; 
gDR_M=15; %70; 
gA_M=10; %only in soma
gKS=84;
gCaL=0.85;
gKCa_M=5; %only in soma
%reversal Poten
El_GM=-60; %mV
Ek=-80;
Ena=45;

sgsq=sigm*sqrt(dt)/C_MP;  %for white noise

vltThres=5; %(mV) voltage threshold for spike
minTbs=floor(0.5/dt); %min time between spikes in ms
indLst=1;
transTime=5000; %in ms
transInd=round(transTime/dt)+1;

%ReqNums=500000; %required # of spikes

%OUTPUT
spTim=zeros(ReqNums,1); % in Hz, AFTER transient time removed

numSpks=0;

caCon_M=0.05; %only in soma, micrMol/l?
Eca_M=70; %dyn variable
% set InitCond
if(length(icY0)==12)
    vltM=icY0(1);
    xVr_M=icY0(2:end); %gating variables
    xiR1_M=icY0(2:end); %aux for RK4
    xiR2_M=icY0(2:end); %aux for RK4
    xiR3_M=icY0(2:end); %aux for RK4
    xiR4_M=icY0(2:end); %aux for RK4
    
else
    vltM=-79.75;
    ic_aux=[.0001085; .99998; .00961145; .998242; .000877; .85574; .000008; .99368; .6217; .0000039736; .999996];
    xVr_M=ic_aux; %gating variables
    xiR1_M=ic_aux; %aux for RK4
    xiR2_M=ic_aux; %aux for RK4
    xiR3_M=ic_aux; %aux for RK4
    xiR4_M=ic_aux; %aux for RK4
end

%time-loop
j=2;
while(numSpks<ReqNums) %keep going until 
    
    % --- step 1 ---
    rhsV1_M=1/Rm_GM*(vltM-El_GM)+gNa_M*xVr_M(1).^3.*xVr_M(2).*(vltM-Ena)+gNaP*minf(vltM).*(vltM-Ena)+...
        gA_M*xVr_M(3).*xVr_M(4).*(vltM-Ek)+gKS*xVr_M(5).*xVr_M(6).*(vltM-Ek)+...
        gCaL*xVr_M(7).*xVr_M(8).*(vltM-Eca_M)+gKCa_M*xVr_M(9).*(vltM-Ek)+gDR_M*xVr_M(10).^2*xVr_M(11).*(vltM-Ek);
    rhsV1_M=dt/C_MP*(Iapp_M-rhsV1_M);
    
    xiR1_M(1)=dt*(ana(vltM).*(1-xVr_M(1)) - bna(vltM).*xVr_M(1)); %I_Na
    xiR1_M(2)=dt*(ahna(vltM).*(1-xVr_M(2)) - bhna(vltM).*xVr_M(2)); %I_Na
    xiR1_M(3)=dt*3.3./atau(vltM).*(aminf(vltM)-xVr_M(3)); %I_A
    xiR1_M(4)=dt*3.3./hatau(vltM).*(ahinf(vltM)-xVr_M(4)); %I_A
    xiR1_M(5)=dt/10*(ksminf(vltM)-xVr_M(5)); %I_KS
    xiR1_M(6)=dt./kshtau(vltM).*(kshinf(vltM)-xVr_M(6));
    xiR1_M(7)=dt*(cala(vltM).*(1-xVr_M(7)) - calb(vltM).*xVr_M(7)); %I_CaL
    xiR1_M(8)=dt*(calha(vltM).*(1-xVr_M(8)) - calhb(vltM).*xVr_M(8));
    xiR1_M(9)=dt*(kcaa(vltM,caCon_M).*(1-xVr_M(9)) - 0.05*xVr_M(9));
    xiR1_M(10)=dt*(nss_dr(vltM)-xVr_M(10))./taun_dr(vltM);    %I_DR
    xiR1_M(11)=dt*(kss_dr(vltM)-xVr_M(11))./50;
    %update calcium concentration
    Ica_M=gCaL*xVr_M(7).*xVr_M(8).*(vltM-Eca_M);
    caCr1_M=dt*(-Ica_M/(2*F*thick_M) + (.05-caCon_M)/10 );
    
    % --- step 2 ---
    rhsV2_M=1/Rm_GM*(vltM+.5*rhsV1_M-El_GM)+gNa_M*(xVr_M(1)+.5*xiR1_M(1)).^3.*(xVr_M(2)+.5*xiR1_M(2)).*(vltM+.5*rhsV1_M-Ena)+...
        gNaP*minf(vltM+.5*rhsV1_M).*(vltM+.5*rhsV1_M-Ena)+gA_M*(xVr_M(3)+.5*xiR1_M(3)).*(xVr_M(4)+...
        .5*xiR1_M(4)).*(vltM+.5*rhsV1_M-Ek)+gKS*(xVr_M(5)+.5*xiR1_M(5)).*(xVr_M(6)+.5*xiR1_M(6)).*(vltM+.5*rhsV1_M-Ek)+...
        gCaL*(xVr_M(7)+.5*xiR1_M(7)).*(xVr_M(8)+.5*xiR1_M(8)).*(vltM+.5*rhsV1_M-Eca_M)+gKCa_M*(xVr_M(9)+...
        .5*xiR1_M(9)).*(vltM+.5*rhsV1_M-Ek)+gDR_M*(xVr_M(10)+.5*xiR1_M(10)).^2*(xVr_M(11)+.5*xiR1_M(11)).*(vltM+.5*rhsV1_M-Ek);
    rhsV2_M=dt/C_MP*(Iapp_M-rhsV2_M);
    
    xiR2_M(1)=dt*(ana(vltM+.5*rhsV1_M).*(1-(xVr_M(1)+.5*xiR1_M(1))) - bna(vltM+.5*rhsV1_M).*(xVr_M(1)+.5*xiR1_M(1)));
    xiR2_M(2)=dt*(ahna(vltM+.5*rhsV1_M).*(1-(xVr_M(2)+.5*xiR1_M(2))) - bhna(vltM+.5*rhsV1_M).*(xVr_M(2)+.5*xiR1_M(2)));
    xiR2_M(3)=dt*3.3./atau(vltM+.5*rhsV1_M).*(aminf(vltM+.5*rhsV1_M)-(xVr_M(3)+.5*xiR1_M(3))); %I_A
    xiR2_M(4)=dt*3.3./hatau(vltM+.5*rhsV1_M).*(ahinf(vltM+.5*rhsV1_M)-(xVr_M(4)+.5*xiR1_M(4))); %I_A
    xiR2_M(5)=dt/10*(ksminf(vltM+.5*rhsV1_M)-(xVr_M(5)+.5*xiR1_M(5))); %I_KS
    xiR2_M(6)=dt./kshtau(vltM+.5*rhsV1_M).*(kshinf(vltM+.5*rhsV1_M)-(xVr_M(6)+.5*xiR1_M(6)));
    xiR2_M(7)=dt*(cala(vltM+.5*rhsV1_M).*(1-(xVr_M(7)+.5*xiR1_M(7))) - calb(vltM+.5*rhsV1_M).*(xVr_M(7)+.5*xiR1_M(7))); %I_CaL
    xiR2_M(8)=dt*(calha(vltM+.5*rhsV1_M).*(1-(xVr_M(8)+.5*xiR1_M(8))) - calhb(vltM+.5*rhsV1_M).*(xVr_M(8)+.5*xiR1_M(8)));
    xiR2_M(9)=dt*(kcaa(vltM+.5*rhsV1_M,caCon_M+.5*caCr1_M).*(1-(xVr_M(9)+.5*xiR1_M(9))) - 0.05*(xVr_M(9)+.5*xiR1_M(9)));
    xiR2_M(10)=dt*(nss_dr(vltM+.5*rhsV1_M)-(xVr_M(10)+.5*xiR1_M(10)))./taun_dr(vltM+.5*rhsV1_M);    %I_DR
    xiR2_M(11)=dt*(kss_dr(vltM+.5*rhsV1_M)-(xVr_M(11)+.5*xiR1_M(11)))./50;
    %update calcium concentration
    Ica_M=gCaL*(xVr_M(7)+.5*xiR1_M(7)).*(xVr_M(8)+.5*xiR1_M(8)).*(vltM+.5*rhsV1_M-Eca_M);
    caCr2_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+.5*caCr1_M))/10 );
    
    % --- step 3 ---
    rhsV3_M=1/Rm_GM*(vltM+.5*rhsV2_M-El_GM)+gNa_M*(xVr_M(1)+.5*xiR2_M(1)).^3.*(xVr_M(2)+.5*xiR2_M(2)).*(vltM+.5*rhsV2_M-Ena)+...
        gNaP*minf(vltM+.5*rhsV2_M).*(vltM+.5*rhsV2_M-Ena)+gA_M*(xVr_M(3)+.5*xiR2_M(3)).*(xVr_M(4)+...
        .5*xiR2_M(4)).*(vltM+.5*rhsV2_M-Ek)+gKS*(xVr_M(5)+.5*xiR2_M(5)).*(xVr_M(6)+.5*xiR2_M(6)).*(vltM+.5*rhsV2_M-Ek)+...
        gCaL*(xVr_M(7)+.5*xiR2_M(7)).*(xVr_M(8)+.5*xiR2_M(8)).*(vltM+.5*rhsV2_M-Eca_M)+gKCa_M*(xVr_M(9)+...
        .5*xiR2_M(9)).*(vltM+.5*rhsV2_M-Ek)+gDR_M*(xVr_M(10)+.5*xiR2_M(10)).^2*(xVr_M(11)+.5*xiR2_M(11)).*(vltM+.5*rhsV2_M-Ek);
    rhsV3_M=dt/C_MP*(Iapp_M-rhsV3_M);
    
    xiR3_M(1)=dt*(ana(vltM+.5*rhsV2_M).*(1-(xVr_M(1)+.5*xiR2_M(1))) - bna(vltM+.5*rhsV2_M).*(xVr_M(1)+.5*xiR2_M(1)));
    xiR3_M(2)=dt*(ahna(vltM+.5*rhsV2_M).*(1-(xVr_M(2)+.5*xiR2_M(2))) - bhna(vltM+.5*rhsV2_M).*(xVr_M(2)+.5*xiR2_M(2)));
    xiR3_M(3)=dt*3.3./atau(vltM+.5*rhsV2_M).*(aminf(vltM+.5*rhsV2_M)-(xVr_M(3)+.5*xiR2_M(3))); %I_A
    xiR3_M(4)=dt*3.3./hatau(vltM+.5*rhsV2_M).*(ahinf(vltM+.5*rhsV2_M)-(xVr_M(4)+.5*xiR2_M(4))); %I_A
    xiR3_M(5)=dt/10*(ksminf(vltM+.5*rhsV2_M)-(xVr_M(5)+.5*xiR2_M(5))); %I_KS
    xiR3_M(6)=dt./kshtau(vltM+.5*rhsV2_M).*(kshinf(vltM+.5*rhsV2_M)-(xVr_M(6)+.5*xiR2_M(6)));
    xiR3_M(7)=dt*(cala(vltM+.5*rhsV2_M).*(1-(xVr_M(7)+.5*xiR2_M(7))) - calb(vltM+.5*rhsV2_M).*(xVr_M(7)+.5*xiR2_M(7))); %I_CaL
    xiR3_M(8)=dt*(calha(vltM+.5*rhsV2_M).*(1-(xVr_M(8)+.5*xiR2_M(8))) - calhb(vltM+.5*rhsV2_M).*(xVr_M(8)+.5*xiR2_M(8)));
    xiR3_M(9)=dt*(kcaa(vltM+.5*rhsV2_M,caCon_M+.5*caCr2_M).*(1-(xVr_M(9)+.5*xiR2_M(9))) - 0.05*(xVr_M(9)+.5*xiR2_M(9)));
    xiR3_M(10)=dt*(nss_dr(vltM+.5*rhsV2_M)-(xVr_M(10)+.5*xiR2_M(10)))./taun_dr(vltM+.5*rhsV2_M);    %I_DR
    xiR3_M(11)=dt*(kss_dr(vltM+.5*rhsV2_M)-(xVr_M(11)+.5*xiR2_M(11)))./50;
    %update calcium concentration
    Ica_M=gCaL*(xVr_M(7)+.5*xiR2_M(7)).*(xVr_M(8)+.5*xiR2_M(8)).*(vltM+.5*rhsV2_M-Eca_M);
    caCr3_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+.5*caCr2_M))/10 );
    
    % --- step 4 ---
    rhsV4_M=1/Rm_GM*(vltM+rhsV3_M-El_GM)+gNa_M*(xVr_M(1)+xiR3_M(1)).^3.*(xVr_M(2)+xiR3_M(2)).*(vltM+rhsV3_M-Ena)+...
        gNaP*minf(vltM+rhsV3_M).*(vltM+rhsV3_M-Ena)+gA_M*(xVr_M(3)+xiR3_M(3)).*(xVr_M(4)+...
        xiR3_M(4)).*(vltM+rhsV3_M-Ek)+gKS*(xVr_M(5)+xiR3_M(5)).*(xVr_M(6)+xiR3_M(6)).*(vltM+rhsV3_M-Ek)+...
        gCaL*(xVr_M(7)+xiR3_M(7)).*(xVr_M(8)+xiR3_M(8)).*(vltM+rhsV3_M-Eca_M)+gKCa_M*(xVr_M(9)+...
        xiR3_M(9)).*(vltM+rhsV3_M-Ek)+gDR_M*(xVr_M(10)+xiR3_M(10)).^2*(xVr_M(11)+xiR3_M(11)).*(vltM+rhsV3_M-Ek);
    rhsV4_M=dt/C_MP*(Iapp_M-rhsV4_M);
    
    xiR4_M(1)=dt*(ana(vltM+rhsV3_M).*(1-(xVr_M(1)+xiR3_M(1))) - bna(vltM+rhsV3_M).*(xVr_M(1)+xiR3_M(1)));
    xiR4_M(2)=dt*(ahna(vltM+rhsV3_M).*(1-(xVr_M(2)+xiR3_M(2))) - bhna(vltM+rhsV3_M).*(xVr_M(2)+xiR3_M(2)));
    xiR4_M(3)=dt*3.3./atau(vltM+rhsV3_M).*(aminf(vltM+rhsV3_M)-(xVr_M(3)+xiR3_M(3))); %I_A
    xiR4_M(4)=dt*3.3./hatau(vltM+rhsV3_M).*(ahinf(vltM+rhsV3_M)-(xVr_M(4)+xiR3_M(4))); %I_A
    xiR4_M(5)=dt/10*(ksminf(vltM+rhsV3_M)-(xVr_M(5)+xiR3_M(5))); %I_KS
    xiR4_M(6)=dt./kshtau(vltM+rhsV3_M).*(kshinf(vltM+rhsV3_M)-(xVr_M(6)+xiR3_M(6)));
    xiR4_M(7)=dt*(cala(vltM+rhsV3_M).*(1-(xVr_M(7)+xiR3_M(7))) - calb(vltM+rhsV3_M).*(xVr_M(7)+xiR3_M(7))); %I_CaL
    xiR4_M(8)=dt*(calha(vltM+rhsV3_M).*(1-(xVr_M(8)+xiR3_M(8))) - calhb(vltM+rhsV3_M).*(xVr_M(8)+xiR3_M(8)));
    xiR4_M(9)=dt*(kcaa(vltM+rhsV3_M,caCon_M+caCr3_M).*(1-(xVr_M(9)+xiR3_M(9))) - 0.05*(xVr_M(9)+xiR3_M(9)));
    xiR4_M(10)=dt*(nss_dr(vltM+rhsV3_M)-(xVr_M(10)+xiR3_M(10)))./taun_dr(vltM+rhsV3_M);    %I_DR
    xiR4_M(11)=dt*(kss_dr(vltM+rhsV3_M)-(xVr_M(11)+xiR3_M(11)))./50;
    %update calcium concentration
    Ica_M=gCaL*(xVr_M(7)+xiR3_M(7)).*(xVr_M(8)+xiR3_M(8)).*(vltM+rhsV3_M-Eca_M);
    caCr4_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+caCr3_M))/10 );
    
    % FINAL Step... update voltage
    vltM=vltM+1/6*(rhsV1_M+2*rhsV2_M+2*rhsV3_M+rhsV4_M)+sgsq*randn;
    %update gating variables
    xVr_M=xVr_M+1/6*(xiR1_M+2*xiR2_M+2*xiR3_M+xiR4_M);
    
    %update calcium concentration
    %Ica=gCaL*xVr(8).*xVr(9).*(V(:,j)-Eca);
    caCon_M=caCon_M+1/6*(caCr1_M+2*caCr2_M+2*caCr3_M+caCr4_M);
    %update calcium reversal potential, Nernst eqn
    Eca_M=RTovF/2*log(10/caCon_M); %get 70mV when caCon=0.05
    
    %record spikes
    if( ((j-indLst)>minTbs) && vltM>vltThres )
        indLst=j;
        numSpks=numSpks+1;
        spTim(numSpks)=dt*j/1000-transTime/1000; %in sec
    end
    if(j==transInd) %reset after transient
        spTim=zeros(ReqNums,1); %reset spTim
        numSpks=0;
    end
    %increase j
    j=j+1;
end

endCond=[vltM; xVr_M];


fileID=fopen([pwd,'/1dats/',flName,'.dat'],'a');
fprintf(fileID,'%d\n',spTim); %b/c goes along cols 
fclose(fileID);

% --- Aux functions ---
    %Ina m alph & bet
    function a_na=ana(v)
        a_na=.32*(v+45)./(1-exp(-(v+45)./4));
    end
    function b_na=bna(v)
        b_na=-.28*(v+18)./(1-exp((v+18)./5));
    end
    %Ina h alph & bet
    function a_hna=ahna(v)
        a_hna=.128./exp((v+41)./18);
    end
    function b_hna=bhna(v)
        b_hna=4./(1+exp(-(v+18)./5));
    end
    %InaP
    function m_inf=minf(v)
        m_inf=1./(1+exp(-(v+50)./5));
    end
    %I_A
    function a_minf=aminf(v)
       a_minf=1./(1+exp(-(v-17.5)./14));
    end
    function a_tau=atau(v)
        a_tau=25*exp((v+45)./13.3)./(1+exp((v+45)./10));
    end
    function a_hinf=ahinf(v)
        a_hinf=1./(1+exp((v+41.7)./6));
    end
    function ha_tau=hatau(v)
        ha_tau=55.5*exp((v+70)./5.1)./(1+exp((v+70)./5));
    end
    %I_KS (tau=10)
    function ks_minf=ksminf(v)
        ks_minf=1./(1+exp(-(v+34)./6.5));
    end
    function ks_hinf=kshinf(v)
        ks_hinf=1./(1+exp((v+68)./6.6));
    end
    function ks_htau=kshtau(v)
        ks_htau=200+330./(1+exp(-(v+71.6)./6.85));
    end
    %I_CaL
    function cal_a=cala(v)
        cal_a=7.5./(1+exp(-(v-13)./7));
    end
    function cal_b=calb(v)
        cal_b=1.65./(1+exp((v-14)./4));
    end
    function cal_ha=calha(v)
        cal_ha=.0068./(1+exp((v+30)./12));
    end
    function cal_hb=calhb(v)
        cal_hb=.06./(1+exp(-v./11));
    end
    %I_KCA, beta_m=.05
    function kca_a=kcaa(v,CA)
        kca_a=-500*exp((v-65)./27).*(.015-CA)./(1-exp(-(CA-.015)./.0013));
    end
    %I_DR fast delayed rectifier; parms from getKssParm.m, stored in parm_fstDR.mat
    function xx=kss_dr(v)
        xx=.43315*(1+tanh(-(v+13.925)./13.0215))+.1337;
    end

    function xx=nss_dr(v)
        xx=((v+100)./150).^8.5849./(0.5747^8.5849+((v+100)./150).^8.5849);
    end
    function xx=taun_dr(v)
        xx=1./(.27654./exp((v+29.9998)./66.3783)+2.89./(1+exp(-(v-19.0524)./12.8786)));
    end

end %end of main function