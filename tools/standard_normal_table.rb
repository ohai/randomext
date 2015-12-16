include Math

R = 3.442619855899
V = 9.91256303526217e-3
K = 7
M = 64
N = (1<<K)

w = Array.new(N)
k = Array.new(N)
f = Array.new(N)

def pow2(r)
  (1 << r).to_f
end

def sn(x)
  exp(-x*x/2)
end

w[N-1] = V*exp(R*R/2)/pow2(M-K-1);
w[N-2] = R/pow2(M-K-1);
k[N-1] = (R/w[N-1]).floor
f[N-1] = sn(R);
xi = R;

(N-2).downto(1) do |i|
  xi = sqrt(-2*log(sn(xi)+V/xi));
  w[i-1] = xi/pow2(M-K-1);
  k[i] = (xi/w[i]).floor;
  f[i] = sn(xi);
end

k[0] = 0;
f[0] = 1.0;

def print_array(type, name, ary)
  puts "static #{type} #{name}[N] = {"
  ary.each do |val|
    puts "  #{val},"
  end
  puts "};"
end

print_array("double", "w", w)
puts
print_array("uint64_t", "k", k)
puts
print_array("double", "f", f)
