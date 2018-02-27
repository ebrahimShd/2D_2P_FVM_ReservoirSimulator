function [ permTensor ] = perm2Tensor( permx,permy,permz )
%PERM2TENSOR Summary of this function goes here
%   Detailed explanation goes here
    permTensor = NaN(3,3,length(permx));
    for i=1:length(permx)
        permTensor(:,:,i)=[permx(i) 0 0;0 permy(i) 0;0 0 permz(i)];
    end
end

