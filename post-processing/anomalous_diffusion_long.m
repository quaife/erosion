dt=0.001;

ntracer=1000;
iter=10000;
xtra=zeros(iter,1);
ytra=zeros(iter,1);
number=[1 46 82 112 139 164 207 227 246];

for t=1:10;
    
    fid = fopen(['xtracer' num2str(number(t)) '.txt']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    xtracer = cac{1}(1:ntracer*iter,1);
    

    fid = fopen(['ytracer' num2str(number(t)) '.txt']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    ytracer = cac{1}(1:ntracer*iter,1);
    
    %Glue the pieces of streamlines to a long streamline.
    for i=1:ntracer
        k=0;
        dx=0;
        dy=0;
        time=0;
        for j=1:iter
            xtra(j)=xtracer(i+(j-1)*(ntracer));
            ytra(j)=ytracer(i+(j-1)*(ntracer));
            if (j>1 && abs(xtra(j)-xtra(j-1)) > 1.9 )
                k=k+1;
                time(k)=j;
                dx(k) = xtra(j-1)-xtra(j);
                dy(k) = ytra(j-1)-ytra(j);
            end
        end
        if(k > 1)
            for t1=1:k-1;
                xtra(time(t1):time(t1+1)-1)=xtra(time(t1):time(t1+1)-1)+dx(t1);
                ytra(time(t1):time(t1+1)-1)=ytra(time(t1):time(t1+1)-1)+dy(t1);
                dx(t1+1)=dx(t1)+dx(t1+1);
                dy(t1+1)=dy(t1)+dy(t1+1);       
            end
            xtra(time(k):end)=xtra(time(k):end)+dx(k);
            ytra(time(k):end)=ytra(time(k):end)+dy(k);
        elseif(k==1)
            xtra(time(k):end)=xtra(time(k):end)+dx(k);
            ytra(time(k):end)=ytra(time(k):end)+dy(k);
        end
        tracer(:,:,i)=[xtra  ytra]; 
    end
    % find the first moment
    x=tracer(:,1,:);
    x = reshape(x,size(x,1),size(x,3),1);
    y=tracer(:,2,:);
    y = reshape(y,size(y,1),size(y,3));    
    
%    fileID = fopen(['mean' num2str(number(t)) '_long.txt'],'w');
    for j=1:iter
        if j>1
            MSD1(j)= MSD1(j-1)+sum(sqrt((x(j,:)-x(j-1,:)).^2 ...
            +(y(j,:)-y(j-1,:)).^2));
        else
            MSD1=0;
        end
    end
    MSD1 = MSD1/1000;
%    fprintf(fileID,'%25.20f\n',MSD1);
%    fclose(fileID); 
  
%  find the second moment 
    fileID = fopen(['sigma' num2str(number(t)) '_long.txt'],'w');
    for j=1:iter-1
        if j>1
            MSD2(j)= sum((sum(sqrt((x(2:j+1,:)-x(1:j,:)).^2 ...
            +(y(2:j+1,:)-y(1:j,:)).^2),1)-ones(1,ntracer)*MSD1(j)).^2);
        else
            MSD2=0;
        end
    end
    MSD2 = sqrt(MSD2/1000);
    fprintf(fileID,'%25.20f\n',MSD2);
    fclose(fileID); 
%     figure;
%     loglog([dt:dt:(iter-2)*dt], MSD2(2:end),'.');
%     hold on;
end


% draw the second moments with the fitting curves. 
load phi100.mat;

dt=0.001;
load sigma1_long.txt; 
load sigma46_long.txt; load sigma112_long.txt; load sigma164_long.txt;
load sigma207_long.txt; load sigma246_long.txt; 
% t1=dt:dt:(size(sigma1_long,1)-1)*dt;
t2=dt:dt:(size(sigma246_long,1)-1)*dt;
% k1=[1:10 11:3:25 28:10:328 329:100:size(t1,2)];
% k11=[2:11 12:3:26 29:10:329 330:100:size(sigma1_long,1)];
% k2=[1:10 11:3:25 28:10:328 329:100:3329 3330:1000:size(t2,2)];
% k21=[2:11 12:3:26 29:10:329 330:100:3330 3331:1000:size(sigma128_long,1)];
kk = -3:4/75:0.96;
k22 = round(10.^kk/dt,0);
% kkt= -3:4/49:-0.3061;
% k12 = round(10.^kkt/dt,0);
q=polyfit(log10(t2(1000:end))',log10(sigma1_long(1001:end)),1);
q1=polyfit(log10(t2(1000:end))',log10(sigma46_long(1001:end)),1);
q2=polyfit(log10(t2(1000:end))',log10(sigma112_long(1001:end)),1);
q3=polyfit(log10(t2(1000:end))',log10(sigma164_long(1001:end)),1);
q4=polyfit(log10(t2(1000:end))',log10(sigma207_long(1001:end)),1);
t=300:1000;
t1=1:300;
t3=size(t2,2)-9000:size(t2,2);
q5=polyfit(log10(t2(t))',log10(sigma246_long(t+1)),1);
q51=polyfit(log10(t2(t1))',log10(sigma246_long(t1+1)),1);
q52=polyfit(log10(t2(t3))',log10(sigma246_long(t3+1)),1);

h1=loglog(t2,2*sqrt(5)*t2/15*10^6,'k--'); hold on;
set(h1,'Linewidth',1.5);
loglog(t2(k22),sigma246_long(k22+1)*10^5,'x'); hold on;
loglog(t2(k22),sigma207_long(k22+1)*10^4,'o'); hold on;
loglog(t2(k22),sigma164_long(k22+1)*10^3,'*'); hold on;
loglog(t2(k22),sigma112_long(k22+1)*10^2,'v'); hold on;
loglog(t2(k22),sigma46_long(k22+1)*10,'^'); hold on;
loglog(t2(k22),sigma1_long(k22+1),'s'); hold on;


legend( strcat('\phi=100%'),...
    strcat('\phi=', num2str(phi100(246)*100,'%3.2f'),'%'),...
    strcat('\phi=', num2str(phi100(207)*100,'%3.2f'),'%'),...
    strcat('\phi=', num2str(phi100(164)*100,'%3.2f'),'%'),...
    strcat('\phi=', num2str(phi100(112)*100,'%3.2f'),'%'),...
    strcat('\phi=', num2str(phi100(46)*100,'%3.2f'),'%'),...  
    strcat('\phi=', num2str(phi100(1)*100,'%3.2f'),'%'));


h7=loglog(t2(t),t2(t).^q5(1)*10^4.45,'k-.'); hold on;
set(h7,'Linewidth',1.5);
h71=loglog(t2(t1),t2(t1).^q51(1)*10^4.6,'k-.'); hold on;
set(h71,'Linewidth',1.5);
h72=loglog(t2(t3),t2(t3).^q52(1)*10^4.45,'k-.'); hold on;
set(h72,'Linewidth',1.5);
h6=loglog(t2(end-9000:end),t2(end-9000:end).^q4(1)*10^3.3,'k-.'); hold on;
set(h6,'Linewidth',1.5);
h5=loglog(t2(end-9000:end),t2(end-9000:end).^q3(1)*10^2.3,'k-.'); hold on;
set(h5,'Linewidth',1.5);
h4=loglog(t2(end-9000:end),t2(end-9000:end).^q2(1)*10^1.5,'k-.'); hold on;
set(h4,'Linewidth',1.5);
h3=loglog(t2(end-9000:end),t2(end-9000:end).^q1(1)*10^0.6,'k-.'); hold on;
set(h3,'Linewidth',1.5);
h2=loglog(t2(end-9000:end),t2(end-9000:end).^q(1)*10^-0.5,'k-.'); hold on;
set(h2,'Linewidth',1.5);


xlabel('time','FontSize',15);
ylabel('particle spreading','FontSize',15);
set(gca,'ytick',[])
figure;
for n=1:ntracer
    plot(tracer(:,1,n),tracer(:,2,n),'color',[0 n/ntracer 1-n/ntracer ]); hold on;
end
