function y = fourier_func(t,a0,a,b,T)
% Return periodic signal using Fourier coefficients

y = a0;
for i=1:length(a)
  y = y + a(i)*cos(2*pi*i/T*t)+b(i)*sin(2*pi*i/T*t);
end

end
