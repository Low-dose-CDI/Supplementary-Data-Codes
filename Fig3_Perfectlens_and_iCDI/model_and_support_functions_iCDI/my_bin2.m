function out = bin(obj,bin_factor)
for ii=1:size(obj,1)/bin_factor
    for jj=1:size(obj,2)/bin_factor
        out(ii,jj)=sum(sum(obj((ii-1)*bin_factor+1:ii*bin_factor,(jj-1)*bin_factor+1:jj*bin_factor,1,1)));
    end
end
out=out/(bin_factor*bin_factor);
end

