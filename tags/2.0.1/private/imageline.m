function ixl = imageline(siz,ixs,ixe)

% create a line on an image similar to Bresenhams algorithm
%
% adopted from Chris Luengo 
% http://www.cb.uu.se/~cris/blog/index.php/archives/400#more-400
%
% 

[p1(1),p1(2)] = ind2sub(siz,ixs);
[p2(1),p2(2)] = ind2sub(siz,ixe);

p = p1;
d = p2-p1;
N = max(abs(d));
s = d/N;

subs = zeros(N,2);
subs(1,:) = p1;

for ii=2:N+1
   p = p+s;
   subs(ii,:) = round(p);
end

ixl = sub2ind(siz,subs(:,1),subs(:,2));