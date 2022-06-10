% assuming simple x'= -f(x) + sig for ISI density
% mixture Gaussians w/ additative noise

% can get resu with addit noise!  but make c vary with noise level
% tried: i) fix c, vary sig1/sig2 together & offset (sig1=sig2-.2), 
% ii) fix c & sig1, vary sig2, iii) varying c with sig1=sig2 varying works!

sigSqv=(.2:.2:5)';
lenSg=length(sigSqv);

cc=jet(lenSg);

cc1=[0 86 255;51 0 51]./255;

dx=0.01;
xVar=(0:dx:12)';
lenX=length(xVar);

fDens1=zeros(lenX,lenSg); %single peak at wavg
fDens2=zeros(lenX,lenSg); %2 peaks at w1,w2
mnISI_1=zeros(lenSg,1);
mnISI_2=zeros(lenSg,1);
vrISI_1=zeros(lenSg,1);
vrISI_2=zeros(lenSg,1);

w1=4;
w2=8;
wavg=(w1+w2)/2;

Knst = min(xVar.^4/4-(w1+w2+wavg)/3*xVar.^3+(w1*w2+w1*wavg+w2*wavg)/2*xVar.^2-w1*w2*wavg*xVar);

f1=@(x)(x-wavg);
F1=@(x)((x-wavg).^2./2);


cSt=0.65;
cEn=0.1;

for j=1:lenSg
    
    fDens1(:,j)=exp(-2*F1(xVar)./sigSqv(j));
    fDens1(:,j)=fDens1(:,j)./(sum(fDens1(:,j))*dx); %normalize
    
    c=(cEn-cSt)/(lenSg-1)*j+cSt-(cEn-cSt)/(lenSg-1); %weight of mixtures
    %c=.5;
    F2a=@(x)( (x-w1).^2./(2*sigSqv(j)) );
    F2b=@(x)( (x-w2).^2./(2*sigSqv(j)) );
    
    fDens2(:,j)=exp(-F2a(xVar))*c + exp(-F2b(xVar))*(1-c);
    fDens2(:,j)=fDens2(:,j)./(sum(fDens2(:,j))*dx); %normalize
end

mnISI_1=sum(fDens1'*xVar,2)*dx;
mnISI_2=sum(fDens2'*xVar,2)*dx;
vrISI_1=sum( fDens1'*(xVar).^2, 2)*dx - mnISI_1.^2;
vrISI_2=sum( fDens2'*(xVar).^2, 2)*dx - mnISI_2.^2;
CV_1=sqrt(vrISI_1)./mnISI_1;
    CV_1(isnan(CV_1))=0;
CV_2=sqrt(vrISI_2)./mnISI_2;
    CV_2(isnan(CV_2))=0;


h1=figure;
hold on
set(gca,'FontSize',20)
h2=figure;
hold on
set(gca,'FontSize',20)
for j=1:lenSg
    figure(h1)
    plot(xVar,fDens1(:,j),'color',cc(j,:))
    figure(h2)
    plot(xVar,fDens2(:,j),'color',cc(j,:))
end

figure
hold on
plot(sqrt(sigSqv),sqrt(vrISI_2),'color',cc1(2,:),'LineWidth',1)
plot(sqrt(sigSqv),sqrt(vrISI_1),'color',cc1(1,:),'LineWidth',1)
set(gca,'FontSize',20)
ylabel('Std.')
set(gca,'XLim',[sqrt(sigSqv(1))-eps sqrt(sigSqv(end))+eps])
figure
hold on
plot(sqrt(sigSqv),CV_2,'color',cc1(2,:),'LineWidth',1)
plot(sqrt(sigSqv),CV_1,'color',cc1(1,:),'LineWidth',1)
set(gca,'FontSize',20)
ylabel('CV')
set(gca,'XLim',[sqrt(sigSqv(1))-eps sqrt(sigSqv(end))+eps])

% figure
% hold on
% plot(sqrt(sigSqv),mnISI_1,'r--')
% plot(sqrt(sigSqv),mnISI_2,'k--')
% set(gca,'FontSize',18)
% ylabel('Mean')