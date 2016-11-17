function[samples] = sample_m2d2(m, N, prob)
% THIS DOES NOT WORK
% See sample_m2d.m
%
% Draw m samples from frequencies 
%   (n,k) \in {-N/2+1,...,N/2}x{-N/2+1,...,N/2}
% according to the distribution prob.
%
% ex. 
%   N = 128; C = pn2dC(N);
%   abab = sample_m2d2(round(N^2/10), N, @(n,k) pn2d(n,k,N,C)); 
%   figure, imagesc(fftshift(abab)); 
%   colormap('gray'); axis equal, axis tight, axis off
%   title('Sampling map')
%
% Returns NxN matrix of 0s and 1s.
%
  % Create storage.
    
    samples = zeros(N, N);
    for mm = 1:m
        p = 0;
        rand_num = rand;
        for ii = -N/2+1:N/2
            for jj = -N/2+1:N/2
                p = p+prob(ii, jj);
                if rand_num < p
                    samples(ii+N/2+1, jj+N/2+1) = 1;
                end
            end
        end
    end
    samples = fftshift(samples);
end