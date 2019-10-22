load triangulation.mat

N = 256;
nb = 100;

x = reshape(x,N,nb);
y = reshape(y,N,nb);

xcenter = zeros(nb,1);
ycenter = zeros(nb,1);

for k = 1:nb
  xcenter(k) = mean(x(:,k));
  ycenter(k) = mean(y(:,k));
end

clf;
hold on;
plot(x,y,'r')
plot(xcenter,ycenter,'b.')
%plot([xcenter(1) xcenter(2)],[ycenter(1) ycenter(2)],'b')

pairs = [CL(:,1:2); CL(:,2:3); CL(:,[1 3])];

if 1
pairs = unique(pairs,'rows');
% eliminate the first level of redundancy
% this doesn't eliminate pairs that show as permutations of one another

for k = 1:size(pairs,1)
  pairs(k,:) = sort(pairs(k,:));
end
% always order them with smaller vertex first
pairs = unique(pairs,'rows');
% remove the redundancy
end

for k = 1:size(pairs,1)
  plot([xcenter(pairs(k,1)) xcenter(pairs(k,2))],...
       [ycenter(pairs(k,1)) ycenter(pairs(k,2))],'b')
end


%for k = 1:size(pairs,1)
%  ind1 = pairs(k,1);
%  ind2 = pairs(k,2);
%
%  dist2 =  ...
%  (repmat(x(:,ind1),1,N) - repmat(x(:,ind2)',N,1)).^2 + ...
%  (repmat(y(:,ind1),1,N) - repmat(y(:,ind2)',N,1)).^2;
%  clf
%  surf(dist2)
%  pause
%
%
%end

