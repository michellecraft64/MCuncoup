%script to get spike times VARYING noise, fixed level of Iapp

Iapp=120;
flNameB='spIap120_sig';
% Iapp=128;
% flNameB='spIap128_sig';

sigm=(0.1:.1:1.6)';

ReqNums=500000;

flName=[flNameB,num2str(0)];
MCsigFast(Iapp,0,ReqNums,0,flName);

% j=1;
% flName=[flNameB,num2str(j)];
% [endCond]=MCsigFast(Iapp,sigm(j),ReqNums,0,flName);
% 
% %load lstIC
%  for j=2:length(sigm)
%      flName=[flNameB,num2str(j)];
%     [endCond]=MCsigFast(Iapp,sigm(j),ReqNums,endCond,flName);
%     save(['/Volumes/GoogleDrive/My Drive/CRCNS_19-22/synch_MC/lstIC'],'endCond');
% end