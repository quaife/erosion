load data100b/xtar_fmmbary.dat;
load data100b/ytar_fmmbary.dat;

load data100b/nbody.dat;
ntarget=500;
xx = reshape(xtar_fmmbary,ntarget,ntarget);
yy = reshape(ytar_fmmbary,ntarget,ntarget);


%load the data (velocity, vorticity, pressure drop) and
% drop the pressure on velocity and vorticity field
for k=1:275
%load(['x' num2str(k) '.dat']);
%load(['y' num2str(k) '.dat']);
%load(['vort_tar' num2str(k) '.dat']);
%load(['drop_press' num2str(k) '.dat']);
%zzu = reshape(utar_fmmbary,ntarget,ntarget);
%zzv = reshape(vtar_fmmbary,ntarget,ntarget);
%figure;
%surf(xx,yy,zzu));hold on;
%view(2) 
%shading interp 
%colorbar

%surf(xx,yy,zzv);hold on;
%view(2) 
%shading interp 
%colorbar

fid = fopen(['data100b/utar' num2str(k) '.dat']);
cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
fclose( fid );
u = cac{1}(:,1);

fid = fopen(['data100b/vtar' num2str(k) '.dat']);
cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
fclose( fid );
v = cac{1}(:,1);

 fid = fopen(['vort_tar' num2str(k) '.dat']);
 cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
 fclose( fid );
 vort = cac{1}(:,1);

fid = fopen(['data100b/drop_press' num2str(k) '.dat']);
cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
fclose( fid );
drop_press = cac{1}(1:1000,1);

%dy=0.99*2/500
pminus = sum(drop_press(23:500-23))/(477-23); 
pplus = sum(drop_press(523:1000-23))/455; 
p=pminus-pplus;


u_tar(:,:,k) = reshape(u/p,ntarget,ntarget);
v_tar(:,:,k) = reshape(v/p,ntarget,ntarget);
vor(:,:,k) = reshape(vort/p,ntarget,ntarget);
% 
% fid = fopen(['press_tar' num2str(k) '.dat']);
% cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
% fclose( fid );
% press = cac{1}(:,1);
% press_tar(:,:,k) = reshape(press/p,ntarget,ntarget);


end


%Flow rate U
load data100b/nbody.dat;
load pressure_drop100b.mat;

n= size(nbody,1);
time=linspace(0,1,n)';
semilogy(time, p(end)*p.^(-1)); hold on
xlabel('normalized time','FontSize',15); 
ylabel('flow rate U','FontSize',15);
axis([0 1 0.0003 1]);
%semilogy(time, 0.0001*exp(log(10000)*time));



% the video of vorticity

v = VideoWriter('vort_dense_100b_hr_complete.avi');
open(v);

for k=1:275
     k=200;
     figure;
%   surf(xx,yy,log10(u_tar(:,:,k).^2+v_tar(:,:,k).^2)/2);hold on;
    surf(xx,yy,vor(:,:,k));hold on;
    %surf(xx,yy,press_tar(:,:,k));hold on;
    %contour(xx,yy,press_tar(:,:,k));hold on;
    view(2) 
    shading interp 
    c=colorbar
    caxis([-0.2 0.2]);
    set(c,'FontSize',19);
    %title(['t=',num2str(k*0.0001)]);
    
    
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
    %for 
    for i=1:nb
       fill(x(1+ninner*(i-1):ninner*i),y(1+ninner*(i-1):ninner*i),'k',...
         'FaceColor', [153/256 153/256 153/256],'EdgeColor', 'none','linewidth',0.1); hold on;
%          fill3(x(1+ninner*(i-1):ninner*i),y(1+ninner*(i-1):ninner*i),...
%              ones(ninner,1),[153/256 153/256 153/256],'EdgeColor', 'none'); hold on;
%        centx = sum(x(1+ninner*(i-1):ninner*i))/ninner;
%        centy = sum(y(1+ninner*(i-1):ninner*i))/ninner;
%        text(centx-0.035,centy,num2str(i)); hold on;
    end
    axis equal;
    axis off;
    %print(gcf,'foo.png','-dpng','-r300'); 
    f = gcf;
    f.Position = [0 0 1200 900];
    if k==1
     set(gca,'units','pixels'); % set the axes units to pixels
     xpic = get(gca,'position'); % get the position of the axes
     set(gcf,'units','pixels'); % set the figure units to pixels
     ypic = get(f,'position'); % get the figure position
    end
    set(f,'position',[ypic(1) ypic(2) xpic(3) xpic(4)]);% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0.02 0.02 0.9 0.95]);
    F(k) = getframe(f);
    writeVideo(v,F(k));
    clf
end
close(v);

clf;
plot(vor(1,:),'b');hold on;
plot(vor(end,:),'r'); pause



%fid = fopen(['x' num2str(k) '.dat']);
%cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
%fclose( fid );
%x = cac{1}(:,1);

%fid = fopen(['y' num2str(k) '.dat']);
%cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
%fclose( fid );
%y = cac{1}(:,1);


%nb=nbody(k);
%ninner=256;
%for i=1:nb
%fill(x(1+ninner*(i-1):ninner*i),y(1+ninner*(i-1):ninner*i),'k',...
%    'FaceColor', [153/256 153/256 153/256],'EdgeColor', 'k','linewidth',0.1); hold on;
%end
%axis equal;
%axis off;
