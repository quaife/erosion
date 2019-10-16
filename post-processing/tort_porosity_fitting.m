load tort_volumn100.mat;
load phi.mat;

T_volumn = T100_volumn;
phi=phi100;

%T=1+p(1-phi)
p1=(T_volumn(1:end)-1)*(1-phi(1:end))/...
    ((1-phi(1:end))'*(1-phi(1:end)));
S1=sum((T_volumn(1:end)-(1+p1*(1-phi(1:end))')).^2);

%T=1-p*ln(phi)
p2=-(T_volumn(1:end)-1)*log(phi(1:end))/...
    (log(phi(1:end))'*log(phi(1:end)));
S2=(T_volumn(1:end)-1+p2*log(phi(1:end))')*(T_volumn(1:end)-1+p2*log(phi(1:end))')'
RMSD2=sqrt(S2/275);

%T=[1-p(1-phi)]^2
a=sum((1-phi(1:end)).^4);
b=sum(3*(1-phi(1:end)).^3);
c=(3-T_volumn(1:end))*(1-phi(1:end)).^2;
d=-(T_volumn(1:end)-1)*(1-phi(1:end));
p3= roots([a b c d]);
S31=sum((T_volumn(1:end)-(1+p3(1)*(1-phi(1:end))').^2).^2);
S32=sum((T_volumn(1:end)-(1+p3(2)*(1-phi(1:end))').^2).^2);
S33=sum((T_volumn(1:end)-(1+p3(3)*(1-phi(1:end))').^2).^2);

%T=phi^-p(p=-x)
x=-0.245:-0.00001:-0.246;%-0.206:-0.00001:-0.208;
phi1=repmat(phi(1:end),size(x)).^repmat(x,size(phi(1:end)));
phi2=repmat(phi(1:end),size(x)).^repmat(2*x,size(phi(1:end)));

y=T_volumn(1:end).*log(phi(1:end))'*phi1...
    -log(phi(1:end))'*phi2;

p=93*(-0.00001)-0.245;

dsdp=T_volumn(1:end).*log(phi(1:end))'*(phi(1:end).^p)...
    -log(phi(1:end))'*(phi(1:end).^(2*p));

s=sum((T_volumn(45:end)-phi(45:end)'.^p).^2);
RMSD = sqrt(s/275);


plot(phi, T_volumn,'.'); hold on;
plot(phi(1:end),1-p2*log(phi(1:end)),'k--','LineWidth',1.5); hold on;
plot(phi(1:end),phi(1:end).^p,'r--');