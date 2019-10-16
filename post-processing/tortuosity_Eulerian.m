load xtar_fmmbary.dat;
load ytar_fmmbary.dat;

ntarget=500;
xx = reshape(xtar_fmmbary,ntarget,ntarget);
yy = reshape(ytar_fmmbary,ntarget,ntarget);

xgrid=xx(1,:);
ygrid=[-1 yy(:,1)' 1];

for k=1:275;
    fid = fopen(['utar' num2str(k) '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    u = cac{1}(:,1);
    u =[zeros(1,500); reshape(u,ntarget,ntarget); zeros(1,500)];

    fid = fopen(['vtar' num2str(k) '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    v = cac{1}(:,1);
    v = [zeros(1,500); reshape(v,ntarget,ntarget); zeros(1,500)];


    F=sqrt(u.^2+v.^2);

    T1=trapz(ygrid,trapz(xgrid,F,2));
    T2=trapz(ygrid,trapz(xgrid,u,2));

    T100_volumn(k)=T1/T2;
end

% k=[1 46 82 112 139 164 186 207 227 246];
p2=plot(phi100,T100_volumn,'.'); hold on;
% p1=plot(phi100(k),T_geo,'r*','MarkerSize',8); hold on;
% p1=plot(phi100(100),T_geo100,'rs','MarkerSize',8);
%axis([0.35 1 1 1.22])
xlabel('\phi','FontSize',20);ylabel('T','FontSize',20)
%legend(p1,'Lagrangian approach',p2,'Eulerian approach');

figure;
plot(phi, 1-0.2216*log(phi),'.'); hold on;
plot(phi, 1+0.3533*(1-(phi)),'.'); hold on;
axis([0.01 1 0.01 1]);


