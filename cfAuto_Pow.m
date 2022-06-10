% script to compare Auto & Power Specdt

load dFI
numDats=86-5+1;
flNameB='sps';
nzInd=(5:86)'; %indices of non-zero firing
frate0=frate;

load dFIsgmPt5.mat %!!!! assuming SAME Iap_v!!!!
frateP5=frate; clear frate;
numDatsP5=30-4+1;
flNameBP5='SgPtFv';
nzIndP5=(4:30)';
ind_diffP5=1; %add this much onto jInd to get corresponding sgm1 Iap_v

load dFIsgm1.mat %!!!! assuming SAME Iap_v!!!!
frate1=frate; clear frate;
numDats1=86-3+1;
flNameB1='SgOne';
nzInd1=(3:86)';
ind_diff1=3; %add this much onto jInd to get corresponding sgm1 Iap_v


load dFIsgm1pt5.mat 
frate1p5=frate; clear frate;
numDats1p5=31-3+1;
flNameB1p5='SgOnPtFv';
nzInd1p5=(3:31)'; %DIFF Iap_v2
ind_diff1p5=3; %add this much onto jInd to get corresponding sgm1 Iap_v

maxFreq=120; %in Hz
numSpks=100000; %numspikes for Power Spect
frqAll=(0:maxFreq)';  %freq for Pow-Spect for all
%% 

% %ndex of no-noise 
% jInd=2; %10--12, 8/9 are much weaker, COMPLETE set: (1:26) b/c only went to Iap_v(30)=144 in SgPtFv
% dt=0.001; %time units (sec) for Autoc
% dtP=.001; %time units (sec) for PowerSpect

% jInd=12;
% dt=0.001; %time units (sec) for Autoc
% dtP=.00005; %time units (sec) for PowerSpect

jInd=22;
dt=0.001; %time units (sec) for Autoc
dtP=.00005; %time units (sec) for PowerSpect

disp(['Curr input Iap=',num2str(Iap_v(jInd+4))]);

flName=[flNameB,num2str(jInd)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
% get Autoc
ySp=sparse(zeros(ceil(spt(numSpks)/dt)+1,1));
indOnz=round(spt(1:numSpks)./dt)+1;
ySp(indOnz)=1;
acfSpks0=autocorr(full(ySp),'NumLags',round(.1/dt));
taut0=(0:length(acfSpks0)-1)'*dt;
% get PowSpec
ySp=sparse(zeros(ceil(spt(numSpks)/dtP)+1,1));
indOnz=round(spt(1:numSpks)./dtP)+1;
ySp(indOnz)=1;
n=length(ySp); 
frq0=(0:n-1)*1/dtP/n;
numPS=round(maxFreq*dtP*n)+1; %up to 120Hz only
frq0=frq0(1:numPS); %truncate
powSpct0=abs(fft(full(ySp))).^2/n;
powSpct0=powSpct0(1:numPS);
%ps0=interp1(frq0,powSpct0,frqAll,'pchip');


% -- repeat with sigm=0.5
flName=[flNameBP5,num2str(jInd+ind_diffP5)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
% get Autoc
ySp=sparse(zeros(ceil(spt(numSpks)/dt)+1,1));
indOnz=round(spt(1:numSpks)./dt)+1;
ySp(indOnz)=1;
acfSpksP5=autocorr(full(ySp),'NumLags',round(.1/dt));
tautP5=(0:length(acfSpksP5)-1)'*dt;
% get PowSpec
ySp=sparse(zeros(ceil(spt(numSpks)/dtP)+1,1));
indOnz=round(spt(1:numSpks)./dtP)+1;
ySp(indOnz)=1;
n=length(ySp); 
frqP5=(0:n-1)*1/dtP/n;
numPS=round(maxFreq*dtP*n)+1; %up to 120Hz only
frqP5=frqP5(1:numPS); %truncate
powSpctP5=abs(fft(full(ySp))).^2/n;
powSpctP5=powSpctP5(1:numPS);
psP5=interp1(frqP5,powSpctP5,frqAll,'pchip');

% -- repeat with sigm=1
flName=[flNameB1,num2str(jInd+ind_diff1)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
% get Autoc
ySp=sparse(zeros(ceil(spt(numSpks)/dt)+1,1));
indOnz=round(spt(1:numSpks)./dt)+1;
ySp(indOnz)=1;
acfSpks1=autocorr(full(ySp),'NumLags',round(.1/dt));
taut1=(0:length(acfSpks1)-1)'*dt;
% get PowSpec
ySp=sparse(zeros(ceil(spt(numSpks)/dtP)+1,1));
indOnz=round(spt(1:numSpks)./dtP)+1;
ySp(indOnz)=1;
n=length(ySp); 
frq1=(0:n-1)*1/dtP/n;
numPS=round(maxFreq*dtP*n)+1; %up to 120Hz only
frq1=frq1(1:numPS); %truncate
powSpct1=abs(fft(full(ySp),n)).^2/n;
powSpct1=powSpct1(1:numPS);
ps1=interp1(frq1,powSpct1,frqAll,'pchip');

% -- repeat with sigm=1.5
flName=[flNameB1p5,num2str(jInd+ind_diff1p5)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
% get Autoc
ySp=sparse(zeros(ceil(spt(numSpks)/dt)+1,1));
indOnz=round(spt(1:numSpks)./dt)+1;
ySp(indOnz)=1;
acfSpks1p5=autocorr(full(ySp),'NumLags',round(.1/dt));
taut1p5=(0:length(acfSpks1p5)-1)'*dt;
% get PowSpec
ySp=sparse(zeros(ceil(spt(numSpks)/dtP)+1,1));
indOnz=round(spt(1:numSpks)./dtP)+1;
ySp(indOnz)=1;
n=length(ySp); 
frq1p5=(0:n-1)*1/dtP/n;
numPS=round(maxFreq*dtP*n)+1; %up to 120Hz only
frq1p5=frq1p5(1:numPS); %truncate
powSpct1p5=abs(fft(full(ySp),n)).^2/n;
powSpct1p5=powSpct1p5(1:numPS);
ps1p5=interp1(frq1p5,powSpct1p5,frqAll,'pchip');


figure
box off
hold on
plot(taut0,acfSpks0,'color',[0 0 0],'LineWidth',1)
plot(tautP5,acfSpksP5,'color',[0 0.435294117647059   0.443137254901961],'LineWidth',1)
plot(taut1,acfSpks1,'color',[.5 0 0],'LineWidth',1)
plot(taut1p5,acfSpks1p5,'color',[1 0 0],'LineWidth',1)
set(gca,'FontSize',18)
xlabel('Time Lag (s)')
ylabel('Autocov')

figure
subplot(2,2,1)
hold on
plot(frq0,powSpct0./powSpct0(1),'color',[0 0 0],'LineWidth',1)
set(gca,'FontSize',18)
ylabel('Power Spectrum')
xlabel('Frequency (Hz)')
set(gca,'XLim',[1 maxFreq])
subplot(2,2,2)
hold on
plot(frqAll,psP5./psP5(1),'color',[0 0.435294117647059   0.443137254901961],'LineWidth',1)
set(gca,'XLim',[1 maxFreq])
subplot(2,2,3)
hold on
plot(frqAll,ps1./ps1(1),'color',[.5 0 0],'LineWidth',1)
set(gca,'XLim',[1 maxFreq])
subplot(2,2,4)
hold on
plot(frqAll,ps1p5./ps1p5(1),'color',[1 0 0],'LineWidth',1)
set(gca,'XLim',[1 maxFreq])

%figure('Renderer', 'Painters');
figure
hold on
plot(frq0,powSpct0./powSpct0(1),'color',[0 0 0],'LineWidth',1)
set(gca,'FontSize',18)
ylabel('Power Spectrum')
xlabel('Frequency (Hz)')
plot(frqAll,psP5./psP5(1),'color',[0 0.435294117647059   0.443137254901961],'LineWidth',1)
plot(frqAll,ps1./ps1(1),'color',[.5 0 0],'LineWidth',1)
plot(frqAll,ps1p5./ps1p5(1),'color',[1 0 0],'LineWidth',1)
set(gca,'XLim',[1 maxFreq])
set(gca,'YLim',[0 2e-4])

% figure
% box off
% loglog(frq0,powSpct0./powSpct0(1),'color',[0 0 1],'LineWidth',1)
% hold on
% loglog(frqP5,powSpctP5./powSpctP5(1),'color',[0 1 1],'LineWidth',2)
% loglog(frq1,powSpct1./powSpct1(1),'color',[0 0 0.5],'LineWidth',2)
% set(gca,'FontSize',18)
% ylabel('Power Spectrum')
% xlabel('Frequency (Hz)')
% figure
% box off
% loglog(frq0,powSpct0./powSpct0(1),'color',[0 0 1],'LineWidth',1)
% hold on
% loglog(frqAll,psP5./psP5(1),'color',[0 1 1],'LineWidth',2)
% loglog(frqAll,ps1./ps1(1),'color',[0 0 0.5],'LineWidth',2)
% set(gca,'FontSize',18)
% ylabel('Power Spectrum')
% xlabel('Frequency (Hz)')