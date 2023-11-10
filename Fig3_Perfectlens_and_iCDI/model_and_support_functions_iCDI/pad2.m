function out = pad2(in, desired_array_size)
% in                    -- input image to be padded
% desired_array_size    -- [y_size, x_size]. Desired y and x sizes to be padded to
 
[sy, sx] = size(in);
padding = ceil((desired_array_size - size(in))./2);
 
% pad in y
if mod(sy, 2) ~= 0 && mod(desired_array_size(1), 2) ~= 0 % odd-odd
    out = padarray(in, [padding(1), 0]);
elseif mod(sy, 2) ~= 0 && mod(desired_array_size(1), 2) == 0 % odd-even
    out = padarray(in, [padding(1), 0]);
    out(end,:) = [];
elseif mod(sy, 2) == 0 && mod(desired_array_size(1), 2) ~= 0 % even-odd 
    out = padarray(in, [padding(1), 0]);
    out(1,:) = [];
else % even-even
    out = padarray(in, [padding(1), 0]);
end
 
% pad in x
if mod(sx, 2) ~= 0 && mod(desired_array_size(2), 2) ~= 0 % odd-odd
    out = padarray(out, [0, padding(2)]);
elseif mod(sx, 2) ~= 0 && mod(desired_array_size(2), 2) == 0 % odd-even
    out = padarray(out, [0, padding(2)]);
    out(:,end) = [];
elseif mod(sx, 2) == 0 && mod(desired_array_size(2), 2) ~= 0 % even-odd 
    out = padarray(out, [0, padding(2)]);
    out(:,1) = [];
else % even-even
    out = padarray(out, [0, padding(2)]);
end

