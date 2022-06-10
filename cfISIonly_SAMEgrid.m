% script to compare ISI (no noise, with noise, etc.) ONLY
% but NOW making meshes all the same, per Rev2

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

%% 
%idnex of no-noise 
%jInd=2; 
%dti=.01;

jInd=12; 
dti=.001;

% jInd=22; 
% dti=0.001;

seriCorr1=zeros(4,1);

flName=[flNameB,num2str(jInd)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
isi_v=diff(spt);
    crM=corrcoef(isi_v(1:2:end-1),isi_v(2:2:end));
    seriCorr1(1)=crM(1,2);
%edges=(min(isi_v) -.5*dti : dti : max(isi_v)+eps)';
edges=(0.5*dti : dti : max(isi_v)+.005)';

[fisi0,edges]=histcounts(isi_v,edges,'normalization','pdf');
xv0=edges(1:end-1)+.5*(edges(3)-edges(2));

disp(['Curr input Iap=',num2str(Iap_v(jInd+4))]);
mn0=mean(isi_v)
sum(fisi0*xv0)*(xv0(2)-xv0(1))
std0=std(isi_v)
sqrt(sum(fisi0*xv0.^2)*(xv0(2)-xv0(1))-(sum(fisi0*xv0)*(xv0(2)-xv0(1)))^2)
CV0=std0/mean(isi_v)

% -- repeat with sigm=0.5
flName=[flNameBP5,num2str(jInd+ind_diffP5)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
isi_v=diff(spt);
    crM=corrcoef(isi_v(1:2:end-1),isi_v(2:2:end));
    seriCorr1(2)=crM(1,2);
%dti=mean(isi_v)/12.15;
%edges=(min(isi_v) -.5*dti : dti : max(isi_v)+eps)';
edges=(.5*dti : dti : max(isi_v)+eps)';
[fisiP5,edges]=histcounts(isi_v,edges,'normalization','pdf');
xvP5=edges(1:end-1)+.5*(edges(3)-edges(2));

mnP5=mean(isi_v)
sum(fisiP5*xvP5)*(xvP5(2)-xvP5(1))
stdP5=std(isi_v)
sqrt(sum(fisiP5*xvP5.^2)*(xvP5(2)-xvP5(1))-(sum(fisiP5*xvP5)*(xvP5(2)-xvP5(1)))^2)
CVP5=stdP5/mean(isi_v)

% -- repeat with sigm=1
flName=[flNameB1,num2str(jInd+ind_diff1)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
isi_v=diff(spt);
    crM=corrcoef(isi_v(1:2:end-1),isi_v(2:2:end));
    seriCorr1(3)=crM(1,2);
%dti=mean(isi_v)/10;
%edges=(min(isi_v) -.5*dti : dti : max(isi_v)+eps)';
edges=(0.5*dti : dti : max(isi_v)+eps)';
[fisi1,edges]=histcounts(isi_v,edges,'normalization','pdf');
xv1=edges(1:end-1)+.5*(edges(3)-edges(2));

mn1=mean(isi_v)
sum(fisi1*xv1)*(xv1(2)-xv1(1))
std1=std(isi_v)
sqrt(sum(fisi1*xv1.^2)*(xv1(2)-xv1(1))-(sum(fisi1*xv1)*(xv1(2)-xv1(1)))^2)
CV1=std1/mean(isi_v)


% -- repeat with sigm=1.5
flName=[flNameB1p5,num2str(jInd+ind_diff1p5)];
spt=importdata([pwd,'/1dats/',flName,'.dat']);
isi_v=diff(spt);
    crM=corrcoef(isi_v(1:2:end-1),isi_v(2:2:end));
    seriCorr1(4)=crM(1,2);
%dti=0.003;%mean(isi_v)/8;%(max(isi_v)-min(isi_v))/44;
%edges=(min(isi_v) -.5*dti : dti : max(isi_v)+eps)';
edges=(0.5*dti : dti : max(isi_v)+eps)';
[fisi1p5,edges]=histcounts(isi_v,edges,'normalization','pdf');
xv1p5=edges(1:end-1)+.5*(edges(3)-edges(2));

mn1p5=mean(isi_v)
sum(fisi1p5*xv1p5)*(xv1p5(2)-xv1p5(1))
std1p5=std(isi_v)
sqrt(sum(fisi1p5*xv1p5.^2)*(xv1p5(2)-xv1p5(1))-(sum(fisi1p5*xv1p5)*(xv1p5(2)-xv1p5(1)))^2)
CV1p5=std1p5/mean(isi_v)

figure
box off
hold on
plot(xv0,fisi0,'color',[0 0 1],'LineWidth',1)
plot(xvP5,fisiP5,'color',[0 1 1],'LineWidth',2)
plot(xv1,fisi1,'color',[0 0 0.5],'LineWidth',2)
plot(xv1p5,fisi1p5,'color',[0 .5 0],'LineWidth',2)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('ISI PDF')

%save(['dMC_isi_Iap',num2str(Iap_v(jInd+4))],'xv*','fisi*')

% %if you want to show std,CV for increased noise:
% x=[std0 stdP5 std1 std1p5]';
% y=[CV0 CVP5 CV1 CV1p5]';
% figure
% subplot(2,1,1)
% plot(x,'o-')
% subplot(2,1,2)
% plot(y,'o-')