load nbody.dat;
n=size(nbody,1);
phi100=zeros(size(nbody));
for k=1:275
    %input the gemo data
    fid = fopen(['x' num2str(k) '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    x = cac{1}(:,1);

    fid = fopen(['y' num2str(k) '.dat']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    y = cac{1}(:,1);

    %calculate the bodies area
    nb=nbody(k);
    ninner=256;
    for i=1:nb
       phi100(k)=phi100(k)+polyarea(x(1+ninner*(i-1):ninner*i),y(1+ninner*(i-1):ninner*i));
    end
    
    %find the porosity in [-1,1]x[-1,1] (area([-1,1]x[-1,1])=4)
    phi100(k)= (4-phi100(k))/4;
end

for k=1:274
    dphi(k) = phi(k+1)-phi(k);
end

%save phi100 as phi.amt
load phi.mat;
x=linspace(0,1,n);
%area fraction
figure
plot(x,1-phi100);
xlabel('normalized time','FontSize',15); 
ylabel('area fraction','FontSize',15);
axis([0 1 0 0.5]);
axis equal;

%eroding speed
figure
plot(x(1:274),dphi*40000)
xlabel('normalized time','FontSize',15); 
ylabel('erdoing speed','FontSize',15);

