using ChirpSignal

fL = 0
fH = 100
s = chirp(1.0, 16000, fL, fH, method="linear");

n = 0.5randn(length(s));
x = s + n;

y = emd(x, 17, minsnr=0, maxiters=3);
N = length(y)
for i = 1:N
    plot(x)
    plot!(y[i], leg=nothing)
    gui();sleep(0.5)
end

