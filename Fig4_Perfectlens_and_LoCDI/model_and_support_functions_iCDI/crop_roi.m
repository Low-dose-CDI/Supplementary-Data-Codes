function [cropped_im] = crop_roi(image, crop_size,centrex,centrey)
% Crop an image to a specified crop_size, centered around centrex, centrey
%   Inputs: 
%       image - image to be cropped
%       crop_size - size of cropped image, can specify [y_dim, x_dim]
%       centrex - x center of cropped image position
%       centrey - y center of cropped image position
%   Outputs:
%       cropped_im - cropped image

if length(crop_size) ~= 1
   crop_size_x = crop_size(2);
   crop_size_y = crop_size(1);
else
   crop_size_x = crop_size;
   crop_size_y = crop_size;
end

if nargin<3 || nargin<4
    centrey=floor(length(image)/2)+1;
    centrex=floor(length(image)/2)+1;
end

bigy = size(image,1);
bigx = size(image,2);

ycent = floor(bigy/2)+1;
xcent = floor(bigx/2)+1;

half_crop_size_x = floor(crop_size_x/2);
half_crop_size_y = floor(crop_size_y/2);

if mod(crop_size,2) == 0
    cropped_im = image(centrey - half_crop_size_y:centrey + (half_crop_size_y - 1),...
    centrex - half_crop_size_x:centrex + (half_crop_size_x - 1), :);
else
    cropped_im = image(centrey - half_crop_size_y:centrey + (half_crop_size_y),...
    centrex - half_crop_size_x:centrex + (half_crop_size_x), :);
end