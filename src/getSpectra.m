function [x,X,n,Fs] = getSpectra(m4aFile,srcPath,m4aPath,m4areadPath)
%GETSPECTRA Summary of this function goes here
%{
m4aPath: path of .m4a file directory.
m4areadPath: path of m4aread package.
x: time spectrum of noise.
X: unshifted frequency spectrum of noise.
n: length of time spectrum.
Fs: continuous time sampling rate (Hz).
%}

cd(m4areadPath);

[x,Fs] = m4aread(fullfile(m4aPath,m4aFile.name));
x = x(:,1);
n = length(x);
X = fft(x)/n;

cd(srcPath);

end

