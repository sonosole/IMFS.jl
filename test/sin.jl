using Plots

fs = 16000
t = collect(0 : 1/fs : 10000/fs);
u = @. exp(-(t-0.3)^2/1e-3) + 1;

n = 0.5randn(length(t)) .* u;
s = sin.(2*π*50*t) .* u;
x = s + n;
y = emd(x, 17, minsnr=0, maxiters=3);
N = length(y)
for i = 1:N
    plot(x)
    plot!(y[i], leg=nothing)
    gui();sleep(0.5)
end

