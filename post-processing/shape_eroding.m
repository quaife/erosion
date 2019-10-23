nstep=1001;
ninner=2^8;
figure;
for k = 1:100:nstep; 
    fid = fopen(['geom' num2str(k,'%4.4i') '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',1, 'CollectOutput',true );
    fclose( fid );
    geom = cac{1}(:,1);
    ntar= geom(2);  
    nbody= geom(3);
    fid = fopen(['stress' num2str(k,'%4.4i') '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',1, 'CollectOutput',true );
    fclose( fid );
    stress = cac{1}(1:ninner,1);    
    stress=[ stress; stress(1)];
    %fprintf(fileID,'%i\n',nbody);
      %figure;
      %rectx=[-2 2 2 -2 -2];
      %recty=[-1 -1 1 1 -1];
      %fill(rectx,recty,'k','FaceColor','b','linewidth',1,'Facealpha',0.1,'LineStyle','none'); hold on;
      %plot([-2 2],[1 1],'k-','LineWidth',2); hold on;     
      %plot([-2 2],[-1 -1],'k-','LineWidth',2); hold on;    
      for n=1:nbody;
          x1=[ geom((n-1)*(3*ntar+3)+7+ntar:(n-1)*(3*ntar+3)+6+2*ntar) ;geom((n-1)*(3*ntar+3)+7+ntar)];
          y1=[ geom((n-1)*(3*ntar+3)+7+2*ntar:(n-1)*(3*ntar+3)+6+3*ntar); geom((n-1)*(3*ntar+3)+7+2*ntar)];
           %h=fill(x1,y1,'k',...
          %'FaceColor', [153/256 153/256 153/256],'EdgeColor', 'k','linewidth',0.1); hold on;
          color_line3(x1',y1',zeros(size(x1')),log10(stress)');hold on;
      end
      axis equal;
      axis off;
  
end
colormap jet
c=colorbar 
caxis([-4 1])
set(c,'FontSize',15);
