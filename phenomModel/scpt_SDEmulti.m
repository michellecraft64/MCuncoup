% assuming simple x'= -f(x) + sig for ISI density
% 2 peaks makes it harder for Var(ISI) to go up with more sigv
% related so scp_simpSDE.m BUT now have multip noise (S depend)

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
%F1=@(x)((x-wavg).^2./2);

f2=@(x)((x-w1).*(x-w2).*(x-wavg));
%F2=@(x)(x.^4/4-(w1+w2+wavg)/3*x.^3+(w1*w2+w1*wavg+w2*wavg)/2*x.^2-w1*w2*wavg*x-Knst);

fxd=zeros(size(xVar));
minSg=3.5;
fxd(xVar<=minSg)=1-.95/minSg*xVar(xVar<=minSg);
fxd(xVar>minSg)=.95/(xVar(end)-minSg)*xVar(xVar>minSg)+(1-.95/(xVar(end)-minSg)*xVar(end));
for j=1:lenSg
    newSig=fxd+sigSqv(j);
    integF1=cumsum(f1(xVar)./newSig)*dx;
    integF2=cumsum(f2(xVar)./newSig)*dx;
   fDens1(:,j)=exp(-2*integF1)./newSig;
        fDens1(:,j)=fDens1(:,j)./(sum(fDens1(:,j))*dx); %normalize
   fDens2(:,j)=exp(-2*integF2+2*Knst/mean(newSig))./newSig;
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


figure
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on
for j=1:lenSg
    subplot(2,1,1)
    plot(xVar,fDens1(:,j),'color',cc(j,:))
    subplot(2,1,2)
    plot(xVar,fDens2(:,j),'color',cc(j,:))
end

figure
hold on
plot(sqrt(sigSqv),sqrt(vrISI_2),'color',cc1(2,:),'LineWidth',1)
plot(sqrt(sigSqv),sqrt(vrISI_1),'color',cc1(1,:),'LineWidth',1)
set(gca,'FontSize',20)
ylabel('Std')
set(gca,'XLim',[sqrt(.2)-eps sqrt(5)+eps])

figure
hold on
plot(sqrt(sigSqv),CV_2,'color',cc1(2,:),'LineWidth',1)
plot(sqrt(sigSqv),CV_1,'color',cc1(1,:),'LineWidth',1)
set(gca,'FontSize',20)
ylabel('CV')
set(gca,'XLim',[sqrt(.2)-eps sqrt(5)+eps])


figure
hold on
for j=1:lenSg
    newSigv=fxd+sigSqv(j);
    plot(xVar,newSigv,'color',cc(j,:))
end
set(gca,'FontSize',20)
ylabel('SigSq (S)')