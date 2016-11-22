% compile
mex patchmatchMex.cpp patchmatch_simple.cc img.cc -DDONT_USE_MAIN -O
% test
a = double(imread('a.png'));
b = double(imread('b.png'));
[n,c] = patchmatchMex(a,b,5,0,30);
figure; 
subplot(2,2,1); imagesc(uint8(a)); title('a');
subplot(2,2,2); imagesc(uint8(b)); title('b');
subplot(2,2,3); imagesc(n(:,:,1)); title('dx');
subplot(2,2,4); imagesc(n(:,:,2)); title('dy');
% show help
patchmatchMex
