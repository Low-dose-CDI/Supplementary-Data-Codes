function out = makeSquareMask(radius,imgSize,cenx,ceny)

cent=floor(imgSize/2)+1; 

if nargin < 3
   ceny=cent;
   cenx=cent;
end

if radius <= 1 
   radius=fix(radius*imgSize/2);
end

out = zeros(imgSize,imgSize);
nc = floor(imgSize/2)+1;
out(nc-radius:nc+radius-1,nc-radius:nc+radius-1)=1;

out=circshift(out,[-(cent-ceny) -(cent-cenx)]);