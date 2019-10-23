load data100b/nbody.dat;
%centx = zeros(1,20);
%centy = zeros(1,20);
%d = zeros(462,1);
vanT = zeros(21,1);
t=1;

for k=1:275
    %figure
 %     k=267;
%     l=1;
    fid = fopen(['data100b/x' num2str(k) '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    x = cac{1}(:,1);

    fid = fopen(['data100b/y' num2str(k) '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    y = cac{1}(:,1);

    nb=nbody(k);
    ninner=256;
    centx = zeros(1,nb);
    centy = zeros(1,nb);
    for i=1:nb;
         fill(x(1+ninner*(i-1):ninner*i),y(1+ninner*(i-1):ninner*i),'k',...
           'FaceColor', [153/256 153/256 153/256],'EdgeColor', 'none','linewidth',0.1); hold on;
       centx(i) = sum(x(1+ninner*(i-1):ninner*i))/ninner;
       centy(i) = sum(y(1+ninner*(i-1):ninner*i))/ninner;
%          text(centx(i)-0.05,centy(i),num2str(i)); hold on;
       %center(1,i)= sum(x(1+ninner*(i-1):ninner*i))/ninner;
       %center(2,i)= sum(y(1+ninner*(i-1):ninner*i))/ninner;
    end
    if (k==1) || (nb < oldnb)
        vanT(t)=k;
        t=t+1;
        num=1;
        par=10;
        cx=-num:2*num/(par-1):num;
        cy=num:-2*num/(par-1):-num;
        cent=[cx ones(1,par-1)*num cy(2:par) -ones(1,par-1)*num; -ones(1,par)*num cx(2:par) ones(1,par-1)*num cy(2:par) ];
        DT = delaunayTriangulation([centx';cx'; ones(par-1,1); cy(2:par)'; -ones(par-1,1)],...
           [centy'; -ones(par,1); cx(2:par)'; ones(par-1,1); cy(2:par)']); 
%         DT = delaunayTriangulation(centx',centy');       
       %DT = delaunayTriangulation([center(1,:)';cx'; ones(9,1); cy(2:10)'; -ones(9,1)],...
       %    [center(2,:)'; -ones(10,1); cx(2:10)'; ones(9,1); cy(2:10)']);        
       %%triplot(DT); hold on
        CL = [];
        numDT = size(DT,1);
        for j = 1:numDT; 
            if(max(DT.ConnectivityList(j,:))<=nb); 
               CL = [CL;DT.ConnectivityList(j,:)]; 
            end; 
        end;
        pairs = [CL(:,1:2); CL(:,2:3); CL(:,[1 3])];
        pairs = unique(pairs,'rows');
% eliminate the first level of redundancy
% this doesn't eliminate pairs that show as permutations of one another
        for l = 1:size(pairs,1)
           pairs(l,:) = sort(pairs(l,:));
        end
% always order them with smaller vertex first
        pairs = unique(pairs,'rows');
        CL = pairs;
    end
% remove the redundancy
    oldnb = nb;
    numLink = size(CL,1);
    for n=1:numLink
%          i=min(CL(n,:));
%          A1 = repmat(x(1+ninner*(i-1):ninner*i),1,256);
%          A2 = repmat(y(1+ninner*(i-1):ninner*i),1,256);
%          for m=1:3
%             j=CL(n,m);
%             if j > i
%                B1= repmat(x(1+ninner*(j-1):ninner*j)',256,1);  
%                B2= repmat(y(1+ninner*(j-1):ninner*j)',256,1);
%                d(l,k)=min(min((A1-B1).^2+(A2-B2).^2));
%                l=l+1;
%             end
%          end
          i1 = CL(n,1);
          A1 = repmat(x(1+ninner*(i1-1):ninner*i1),1,256);
          A2 = repmat(y(1+ninner*(i1-1):ninner*i1),1,256);
          i2 = CL(n,2);
          B1 = repmat(x(1+ninner*(i2-1):ninner*i2)',256,1);
          B2 = repmat(y(1+ninner*(i2-1):ninner*i2)',256,1);
          d(n,k)=min(min(sqrt((A1-B1).^2+(A2-B2).^2)));
     end
%           %intersection checking     
%           [xi,yi]=polyxpoly([centx(i) centx(j)],...
%                [centy(i) centy(j)],x(1+ninner*(i-1):ninner*i),...
%                 y(1+ninner*(i-1):ninner*i));
%           [xj,yj]=polyxpoly([centx(i) centx(j)],...
%                 [centy(i) centy(j)],x(1+ninner*(j-1):ninner*j),...
%                 y(1+ninner*(j-1):ninner*j));
%            TF1=isempty(xi);
%            TF2=isempty(xj);
%            if TF1==1
%                xi=centx(i); yi=centy(i);
%            end
%            if TF2==1
%                xj=centx(j); yj=centy(j);
%            end

       %triplot(CL,centx',centy')
       %triplot(CL,center(1,:)',center(2,:)')
    axis equal;
    axis off;
   
end

%calculate gap mean and variance 
for j=1:266
    k=1;
    d_tmp=0;
    for i=1:258
        if (d(i,j)> 0)
            d_tmp(k) = d(i,j);
            k=k+1;
        end
    end
    d_mean(j)=sum(d_tmp)/(k-1);
    d_var(j)=sum((d_tmp-d_mean(j)).^2)/(k-1);
end

 %gaps mean vs porosity
 load phi.mat;
 load gap_mean.mat;
 t=[1 46 112 164 207 246];
 
 figure;
 plot(phi100(1:size(d_mean,2)),d_mean,'-','LineWidth',2); hold on;
 plot(phi100(t),d_mean(t),'r*');
 xlabel('\phi','FontSize',20);
 ylabel('mean of gaps','FontSize',15);
 axis([0.5 1 0.05 0.51])
 
 %gaps variance vs porosity
 load phi.mat;
 load gap_var.mat;
 t=[1 46 112 164 207 246];
 
 figure;
 plot(phi100(1:size(d_var,2)),d_var,'-','LineWidth',2); hold on;
 plot(phi100(t),d_var(t),'r*');
 xlabel('\phi','FontSize',20);
 ylabel('variance of gaps','FontSize',15);
 axis([0.5 1 0 0.101])
 
 %histgram at six porosities
t=[1 46 112 164 207 246];
load phi.mat;
load gap_100b.mat;


for j=1:6
    k=1;
    d_tmp=0;
    for i=1:258
        if (d(i,t(j))> 0)
            d_tmp(k) = d(i,t(j));
            k=k+1;
        end
    end
    figure;
   %pd = fitdist(d_tmp1','Rayleigh');
   x1=round(max(d_tmp),2);
%    x=0:0.01:x1;
%    h=histogram(d_tmp,x);
    h=histogram(d_tmp,15);
   %axis([0.01 0.08 0 10]);
   f1 = gcf;
   %f1.Position = [0 0 700 300];
   xlabel('gap distance','FontSize',20);
   ylabel('count','FontSize',20);
%  pdfEst = pdf(pd,x);
%  line(x,pdfEst)
    
end

 