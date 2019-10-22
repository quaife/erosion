load data100b/nbody.dat;
t=1;

figure
k=100;

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
%  text(centx(i)-0.05,centy(i),num2str(i)); hold on;
       %center(1,i)= sum(x(1+ninner*(i-1):ninner*i))/ninner;
       %center(2,i)= sum(y(1+ninner*(i-1):ninner*i))/ninner;
end
%    if (k==1) || (nb < oldnb)
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
%    end
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
          [d(n),ind(n)]=min(reshape(sqrt((A1-B1).^2+(A2-B2).^2),256^2,1));
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
   
for n=1:254;
    ind1=rem(ind(n),256);
    ind2=((ind(n)-ind1)/256)+1;
    x1tmp = x(ninner*(CL(n,1)-1)+ind1);
    y1tmp = y(ninner*(CL(n,1)-1)+ind1);
    x2tmp = x(ninner*(CL(n,2)-1)+ind2);
    y2tmp = y(ninner*(CL(n,2)-1)+ind2);
    plot([x1tmp x2tmp],[y1tmp, y2tmp],'r','LineWidth',0.8); hold on

end

figure
x=0:0.01:0.25;
%pd = fitdist(d_tmp1','Rayleigh');
h=histogram(d,15);
%axis([0.01 0.08 0 10]);
f1 = gcf;
%f1.Position = [0 0 700 300];
xlabel('gap distance','FontSize',15);
ylabel('counts','FontSize',15);
% pdfEst = pdf(pd,x);
% line(x,pdfEst)