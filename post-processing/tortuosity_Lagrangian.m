
ntracer=1000;
iter=1000;
dx=(1-(-1))/ntracer;
L=2;

k=[1 46 82 112 139 164 186 207 227 246];
%k=100;
for i=1:10;   
    fid = fopen(['xtracer' num2str(k(i)) '.txt']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    x = cac{1}(1:ntracer*iter,1);
    x = reshape(x,ntracer,iter);

    fid = fopen(['ytracer' num2str(k(i)) '.txt']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    y = cac{1}(1:ntracer*iter,1);
    y = reshape(y,ntracer,iter);
    
    fid = fopen(['utracer' num2str(k(i)) '.txt']);
    cac = textscan(fid,'%f%f%f', 'Headerlines',0, 'CollectOutput',true );
    fclose( fid );
    u = cac{1}(1:ntracer,1);
    
    tmp1=(x(:,2:iter)-x(:,1:iter-1)).^2+(y(:,2:iter)-y(:,1:iter-1)).^2;
    
    leng=sum(sqrt(tmp1),2);
    
    %check the j-th tracer if it is at the section x=1
    %if not, use the arclength of (j-1)th tracer to 
    %replace it.
    for j=1:ntracer
        if(x(j,iter) > 0.99)
            leng_new(j)=leng(j)/L;
            u_new(j)=u(j);
            num=j;
        else
            leng_new(j)=leng(num)/L;
            u_new(j)=u(j);
        end
    end
    
    udx=sum(u_new*dx);
    T_geo100(i)=sum(u_new.*(leng_new*dx))/udx;
      
    
end  