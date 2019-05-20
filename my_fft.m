% https://www.mathworks.com/help/matlab/fourier-analysis-and-filtering.html?s_tid=CRUX_lftnav
% esta funcao my_fft pode ser útil para o projecto

function [f, sf]=my_fft(st,fs)
 
 sf=fft(st);
 
 sf=fftshift(sf)./(length(st)-1);
 
 f=[-fs/2:fs/(length(sf)-1):fs/2];