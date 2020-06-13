using Pkg
Pkg.add("QuadGK")
using QuadGK
using Plots

f(x) = x^2

print(quadgk(f, 2,3))

table = zeros(2, 3, 5)
for k in 1:5, j in 1:3, i in 1:2
    table[i, j, k] = i*j*k - i*j
end

print(table)

global x = pi/4
global y = 3*pi/4

for i in x:y
    f[i] = sin(i)
end

plot(f, label="hi")
