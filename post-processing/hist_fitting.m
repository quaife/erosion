load phi.mat;
load gap_mean.mat;
load gap_var.mat;
load gap_100b.mat;


%Weibull distribution
t=[1 46 112 164 207 246];
for j=1:6
    l=1;
    d_tmp=0;
    for i=1:258
        if (d(i,t(j))> 0)
            d_tmp(l) = d(i,t(j));
            l=l+1;
        end
    end
    figure;

    m=d_mean(t(j));
    s2=d_var(t(j));

    a=s2/m^2; 
    f = @(k_tmp) [gamma(1+2/k_tmp)-(gamma(1+1/k_tmp))^2]/(gamma(1+1/k_tmp))^2-a;
    k=fzero(f,[0.1 10]);
    lambda = m/gamma(1+1/k);
    k1(j)=k;
    lambda1(j)=lambda;
    
    x1=round(max(d_tmp),2);


    pdf = @(x)(k/lambda)*((x/lambda).^(k-1)).*exp(-(x/lambda).^k);
    [N, edges]=histcounts(d_tmp,15);
    area = (edges(2:end)-edges(1:end-1))* N';%scaling constant
    x_bar=(edges(1:end-1)+edges(2:end))/2;
    f2 = @(alpha) sum((alpha*pdf(x_bar)-N).^2);
%     alpha_real=fzero(f2,[1 20]);
%     [h,x]=hist(d_tmp,15);
%     bar(x,h/sum(h)); hold on;
     x=0:0.001:edges(end)+1;
    h=histogram(d_tmp,15); hold on;
    plot(x,area*pdf(x),'k-','LineWidth',2); hold on;
    xlabel('pore size','FontSize',20);
    ylabel('count','FontSize',20);
    a = get(gca,'YTickLabel');
   set(gca,'YTickLabel',a,'fontsize',18)
   xlim([0 edges(end)*1.029]);
end

%100 bodies at time 100th unit
load gap_100b100.mat;

m100=sum(d)/size(d,2);
s2_100=sum((d-m100).^2)/size(d,2);

a100=s2_100/m100^2; 
f = @(k_tmp) [gamma(1+2/k_tmp)-(gamma(1+1/k_tmp))^2]/(gamma(1+1/k_tmp))^2-a100;
k100=fzero(f,[0.1 10]);
lambda100 = m100/gamma(1+1/k100);

pdf = @(x)(k100/lambda100)*((x/lambda100).^(k100-1)).*exp(-(x/lambda100).^k100);
[N, edges]=histcounts(d,15);
area = (edges(2:end)-edges(1:end-1))* N';%scaling constant

x=0:0.001:edges(end)+1;
h=histogram(d,15); hold on;
plot(x,area*pdf(x),'k-','LineWidth',2); hold on;
xlabel('pore size','FontSize',20);
ylabel('count','FontSize',20);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',15)
xlim([0 edges(end)*1.029]);


