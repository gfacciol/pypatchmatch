% compile
mex patchmatchMex.cpp patchmatch_simple.cc img.cc -DDONT_USE_MAIN -O
% show help
patchmatchMex
% test
a = double(imread('a.png'));
b = double(imread('b.png'));
[n,c] = patchmatchMex(a,b);
