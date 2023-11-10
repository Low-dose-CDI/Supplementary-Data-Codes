%Function for alignment

% Img1IN is the reference image
% Img2IN will be align based on Img1IN and saved on ImgOUT

function [ImgOUT]=align2D(Img1IN,Img2IN,showim)

if nargin<3
    showim=0;
end

%Img1IN = Img1IN/max(max(max(Img1IN)));
%Img2IN = Img2IN/max(max(max(Img2IN)));

error = Inf;
L = size(Img1IN);

Img = Img2IN;
C = fftshift(ifftn(fftn(Img1IN/max(max(max(Img1IN)))).*conj(fftn(Img)))); % cross correlation
ind = find(C == max(max(max(C))),1,'first'); % maximum cross correlation
[Y,X] = ind2sub(size(C),ind);
X = X - floor(L(2)/2) - 1;
Y = Y - floor(L(1)/2) - 1;
ImgTemp = circshift(Img,[Y X]);

errorTemp = nansum(nansum(nansum(abs(ImgTemp-Img1IN)))); %error


if errorTemp < error
    ImgOUT = ImgTemp; %save better align
    error = errorTemp;
end

Img = rot90(Img);   
% C = fftshift(ifftn(fftn(Img1IN/max(max(max(Img1IN)))).*conj(fftn(Img)))); % cross correlation
% ind = find(C == max(max(max(C))),1,'first'); % maximum cross correlation
% [Y,X] = ind2sub(size(C),ind);
% X = X - floor(L(2)/2) - 1;
% Y = Y - floor(L(1)/2) - 1;
% ImgTemp = circshift(Img,[Y X]);
% 
% errorTemp = nansum(nansum(nansum(abs(ImgTemp-Img1IN)))); %error
% 
% if errorTemp < error
%     ImgOUT = ImgTemp; %save better align
%     error = errorTemp;
% 
% end

Img = rot90(Img); 
C = fftshift(ifftn(fftn(Img1IN/max(max(max(Img1IN)))).*conj(fftn(Img))));% cross correlation
ind = find(C == max(max(max(C))),1,'first'); % maximum cross correlation
[Y,X] = ind2sub(size(C),ind);
X = X - floor(L(2)/2) - 1;
Y = Y - floor(L(1)/2) - 1;
ImgTemp = circshift(Img,[Y X]);

errorTemp = nansum(nansum(nansum(abs(ImgTemp-Img1IN)))); %error



if errorTemp < error
    ImgOUT = ImgTemp; %save better align
    error = errorTemp;
end

if showim~=0
figure(1)
subplot(1,2,1), imagesc(Img1IN), axis image, title('Model');
subplot(1,2,2), imagesc(ImgOUT), axis image, title('Aligned');
end
