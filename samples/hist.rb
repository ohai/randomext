require 'randomext'
require 'tempfile'
require 'benchmark'
require 'bigdecimal'

include Math

module Enumerable
  def sum
    inject(0){|r, i| r+i }
  end
end

class Numeric
  def sign
    return 1 if self > 0
    return -1 if self < 0
    return 0
  end
end

def histogram(num_bins, num_samples, min, max, gen)
  binsize = (max - min).to_f / num_bins
  bins = Array.new(num_bins, 0)
  
  num_samples.times do
    n = (gen[] - min).div(binsize)
    if 0 <= n && n < num_bins
      bins[n] += 1
    end
  end

  return bins.map.with_index{|k, n|
    [binsize*(n+0.5)+min, k/(binsize*num_samples.to_f)]
  }
end

def plot(filename, data, curve, curve_type)
  hist_tmp = Tempfile.new("hist")
  curve_tmp = Tempfile.new("curve")
  
  data.each{|k, v| hist_tmp.puts "#{k} #{v}" }
  hist_tmp.close(false)
  curve.each{|k, v| curve_tmp.puts "#{k} #{v}" }
  curve_tmp.close(false)
  
  IO.popen("gnuplot", "w") do |process|
    process.puts("set term png")
    process.puts("set output \"#{filename}\"")
    process.puts("set nokey")
    process.puts("plot \"#{hist_tmp.path}\" w boxes, \"#{curve_tmp.path}\" w #{curve_type}")
    process.flush
  end
ensure
  hist_tmp.close(true)
  curve_tmp.close(true)
end

def match?(name)
  return true if ARGV.empty?
  ARGV.any?{|pattern| File.fnmatch?(pattern, name) }
end

def draw_histogram(name, num_bins, num_samples, min, max, reporter, gen, func=nil)
  return unless match?(name)
  
  hist = nil
  reporter.report("#{name}:") do
    hist = histogram(num_bins, num_samples, min, max, gen)
  end
  
  if func
    curve = hist.map{|x, _| [x, func[x]] }
  else
    curve = []
  end

  plot("#{name}.png", hist, curve, "lines")
end

def draw_disc_histogram(name, num_samples, reporter, gen, func)
  return unless match?(name)
  
  h = Hash.new(0)
  reporter.report("#{name}:") do
    num_samples.times{ h[gen[]] += 1 }
  end
  hist = h.map{|k, v| [k, v.to_f/num_samples] }.sort_by{|k, v| k }
  curve = hist.map{|x, _| [x, func[x]] }

  plot("#{name}.png", hist, curve, "points")
end

module Distribution
  module_function
  def normal(x, average = 0, sd = 1.0)
    exp(-(x - average)**2/(2 * sd**2))/sqrt(2*PI*sd**2)
  end
  
  def lognormal(x, mu=0.0, sigma=1.0)
    exp(-(log(x)-mu)**2/(2*sigma**2))/(x*sqrt(2*PI)*sigma)
  end

  def cauthy(x, mu, theta)
    theta/(theta**2 + (x-mu)**2)/PI
  end

  def levy(x, mu, theta)
    sqrt(theta/(2*PI))*(x-mu)**(-1.5)*exp(-theta/(2*(x-mu)))
  end

  def exponential(x, theta)
    exp(-x/theta)/theta
  end

  def laplace(x, mu, theta)
    exp(-(x-mu).abs/theta)/(2*theta)
  end

  def rayleigh(x, sigma)
    x/sigma**2*exp(-x**2/(2*sigma**2))
  end

  def weibull(x, g, mu)
    g/mu*(x/mu)**(g-1)*exp(-(x/mu)**g)
  end

  def gumbel(x, mu, nu)
    1/nu*exp(-(x-mu)/nu)*exp(-exp(-(x-mu)/nu))
  end
  
  def gamma(x, alpha, beta)
    beta**(-alpha)*x**(alpha-1)*exp(-x/beta)/Math.gamma(alpha)
  end

  def beta(x, alpha, beta)
    b = Math.gamma(alpha)*Math.gamma(beta)/Math.gamma(alpha+beta)
    x**(alpha-1)*(1-x)**(beta-1)/b
  end

  def power(x, gamma, a, b)
    gamma*(x-a)**(gamma-1)/(b-a)**gamma
  end
  
  def chi_square(x, r)
    x**(r/2.0-1)*exp(-x/2.0)/(2**(r/2.0)*Math.gamma(r/2.0))
  end
  
  def F(x, r1, r2)
    r1 = r1.to_f
    r2 = r2.to_f
    (r1/r2)**(r1/2)*x**(r1/2-1)/(B(r1/2, r2/2)*(1+r1/r2*x)**((r1+r2)/2))
  end

  def B(a, b)
    Math.gamma(a)*Math.gamma(b)/Math.gamma(a+b)
  end

  def t(x, r)
    r = r.to_f
    Math.gamma((r+1)/2)/(sqrt(PI*r)*Math.gamma(r/2)*(1+x**2/r)**((r+1)/2))
  end

  def wald(x, mu, lambda)
    Math.sqrt(lambda/(2*PI*x**3))*exp(-lambda*(x-mu)**2/(2*mu**2*x))
  end
  
  def pareto(x, a, b)
    a*b**a/x**(a+1)
  end
  
  def logistic(x, mu, theta)
    1.0/(4*theta*(cosh((x-mu)/(2*theta))**2))
  end

  def vonmises(x, mu, kappa)
    i0 = 1 + (1..100).inject(0){|r, i|
      l = -2*sumlog(1, i) + 2*i*log(kappa/2)
      r + exp(l)
    } 
    exp(kappa*cos(x-mu))/(2*PI*i0)
  end

  def non_central_chi_square(x, r, l)
    r = r.to_f; l = l.to_f
    u = (0..70).map do |j|
      w = j*log(l.abs/2) - l/2 + (r/2+j-1)*log(x) - x/2 -
        sumlog(1, j) - (r/2+j)*log(2) -lgamma(r/2+j).first
      (l.sign)**j * exp(w)
    end
    u.sum
  end
  
  def non_central_t(x, r, lambda)
    r = r.to_f; lambda = lambda.to_f
    u = (0..70).map do |j|
      w = j*log(lambda.abs) - lambda**2/2 - sumlog(1, j) +
        j/2.0 *log(2.0/r) + lgamma((r+j+1)/2.0).first + j * log(x.abs) -
        0.5*log(PI*r) - lgamma(r/2.0).first -
        (r+j+1)/2.0*log(1 + (x**2.0)/r)
      (x.sign*lambda.sign)**j*exp(w)
    end
    u.sum
  end

  def planck(x, a, b)
    (b**(-a-1) * x**a)/(Math.gamma(a+1) * Zeta(a+1) * (exp(x/b)-1))
  end
  
  def combination(n ,r)
    r = n - r if n/2 < r
    ret = 1
    n.downto(n-r+1) do |k|
      ret = ret*k/(n-k+1)
    end
    ret
  end
  
  def binomial(x, n, p)
    ret = 1.0
    q = 1.0-p
    n.downto(n-x+1) do |k|
      ret = ret*k/(n-k+1)*p
      ret *= q if n-k+1 <= n-x
    end
    (n-2*x).times{ ret *= q }
    return ret
  end

  def geometric(x, theta)
    theta*(1-theta)**(x-1)
  end
  
  def poisson(x, lambda)
    exp(-(1..x).inject(0){|r, i| r+log(i)} - lambda + x*log(lambda))
  end

  def hypergeometric(x, n, m, k)
    exp(logcombination(m, x) + logcombination(n-m, k-x) - logcombination(n, k))
  end
  
  def sumlog(from, to)
    (from .. to).inject(0.0){|r, i| r + Math.log(i)}
  end
  
  def logcombination(n, m)
    sumlog(n-m+1, n) - sumlog(1, m)
  end

  def negative_binomial(x, r, theta)
    exp(lgamma(r+x)[0] - lgamma(r)[0] - sumlog(1, x) +
        r*log(theta) + x*log(1-theta))
  end

  def logseries(x, theta)
    -1/log(1-theta)*theta**x/x
  end
  
  def zipf_mandelbrot(x, n, q=0.0, s=1.0)
    sum = (1..n).inject(0.0){|r, i| r + 1.0/(i+q)**s }
    1.0 / (x+q)**s / sum
  end

  def zeta(x, s)
    x**(-s)/Zeta(s)
  end

  # Riemann's zeta function
  def Zeta(s)
    (1..100).inject(0){|r, n| r + n**(-s) }
  end
end

rng = Random.new

Benchmark.bm(14) do |reporter|
  draw_histogram("snormal", 80, 100000, -6.0, 6.0, reporter,
                 rng.method(:standard_normal), Distribution.method(:normal))
  
  
  draw_histogram("normal", 80, 100000, -3.0, 9.0, reporter,
                 proc{ rng.normal(3.0, 1.7) },
                 proc{|x| Distribution.normal(x, 3.0, 1.7) })

  draw_histogram("lognormal", 100, 100000, 0.0, 8.0, reporter,
                 rng.method(:lognormal), Distribution.method(:lognormal))

  draw_histogram("cauthy", 100, 100000, -40.0, 46.0, reporter,
                 proc{ rng.cauthy(3.0, 1.5) },
                 proc{|x| Distribution.cauthy(x, 3.0, 1.5) })
  draw_histogram("levy", 200, 100000, 0.0, 40.0, reporter,
                 proc{ rng.levy(0.8, 1.2) },
                 proc{|x| Distribution.levy(x, 0.8, 1.2)})
  draw_histogram("exponential", 100, 100000, 0.0, 10.0, reporter,
                 proc{ rng.exponential(1.3) },
                 proc{|x| Distribution.exponential(x, 1.3) })
  draw_histogram("standard_exponential", 100, 100000, 0.0, 10.0, reporter,
                 proc{ rng.standard_exponential },
                 proc{|x| Distribution.exponential(x, 1.0) })
  draw_histogram("laplace", 100, 100000, -6.0, 6.0, reporter,
                 proc{ rng.laplace(0, 1) },
                 proc{|x| Distribution.laplace(x, 0, 1)})
  draw_histogram("rayleigh", 100, 100000, 0.0, 5.0, reporter,
                 proc{ rng.rayleigh(1.2) },
                 proc{|x| Distribution.rayleigh(x, 1.2) })
  draw_histogram("weibull", 100, 100000, 0.0, 5.0, reporter,
                 proc{ rng.weibull(3.0, 1) },
                 proc{|x| Distribution.weibull(x, 3.0, 1) })
  draw_histogram("gumbel", 100, 100000, -4.0, 8.0, reporter,
                 proc{ rng.gumbel(0, 1) },
                 proc{|x| Distribution.gumbel(x, 0, 1) })
  draw_histogram("gamma", 100, 100000, 0.0, 10.0, reporter,
                 proc{ rng.gamma(2.0, 1.0) },
                 proc{|x| Distribution.gamma(x, 2.0, 1.0) })
  draw_histogram("gamma2", 200, 100000, 0.0, 4.0, reporter,
                 proc{ rng.gamma(0.4, 1.0) },
                 proc{|x| Distribution.gamma(x, 0.4, 1.0) })
  draw_histogram("beta", 100, 100000, 0.0, 1.0, reporter,
                 proc{ rng.beta(4.3, 7.2) },
                 proc{|x| Distribution.beta(x, 4.3, 7.2) })
  draw_histogram("beta2", 100, 100000, 0.0, 1.0, reporter,
                 proc{ rng.beta(0.7, 0.42) },
                 proc{|x| Distribution.beta(x, 0.7, 0.42) })
  draw_histogram("power-4", 100, 100000, 0.0, 1.0, reporter,
                 proc{ rng.power(4.0, 0.0, 1.0) },
                 proc{|x| Distribution.power(x, 4.0, 0.0, 1.0) })
  draw_histogram("power-0.5", 100, 100000, 0.0, 1.0, reporter,
                 proc{ rng.power(0.5, 0.0, 1.0) },
                 proc{|x| Distribution.power(x, 0.5, 0.0, 1.0) })
  draw_histogram("chi_square-1",100, 100000, 0.0, 20.0, reporter,
                 proc{ rng.chi_square(1) },
                 proc{|x| Distribution.chi_square(x, 1) })
  draw_histogram("chi_square-5",100, 100000, 0.0, 20.0, reporter,
                 proc{ rng.chi_square(5) },
                 proc{|x| Distribution.chi_square(x, 5) })
  draw_histogram("F-2-8", 100, 100000, 0.0, 10.0, reporter,
                 proc{ rng.F(2, 8) },
                 proc{|x| Distribution.F(x, 2, 8) })
  draw_histogram("F-20-45", 100, 100000, 0.0, 3.5, reporter,
                 proc{ rng.F(20, 45) },
                 proc{|x| Distribution.F(x, 20, 45) })
  draw_histogram("t-1", 100, 100000, -8.0, 8.0, reporter,
                 proc{ rng.t(1) },
                 proc{|x| Distribution.t(x, 1) })
  draw_histogram("t-2", 100, 100000, -8.0, 8.0, reporter,
                 proc{ rng.t(2) },
                 proc{|x| Distribution.t(x, 2) })
  draw_histogram("t-5", 100, 100000, -8.0, 8.0, reporter,
                 proc{ rng.t(5) },
                 proc{|x| Distribution.t(x, 5) })
  draw_histogram("t-20", 100, 100000, -8.0, 8.0, reporter,
                 proc{ rng.t(20) },
                 proc{|x| Distribution.t(x, 20) })
  draw_histogram("wald", 100, 100000, 0.0, 3.0, reporter,
                 proc{ rng.wald(2.0, 1.2) },
                 proc{|x| Distribution.wald(x, 2.0, 1.2) })
  draw_histogram("pareto1-4.5", 100, 100000, 1.0, 4.0, reporter,
                 proc{ rng.pareto(4.5, 1) },
                 proc{|x| Distribution.pareto(x, 4.5, 1) })
  
  draw_histogram("logistic", 100, 100000, -10, 10, reporter,
                 proc{ rng.logistic(0.8, 1.2) },
                 proc{|x| Distribution.logistic(x, 0.8, 1.2) })
  
  draw_histogram("vonmises-0.2", 80, 100000, -PI, PI, reporter,
                 proc{ rng.vonmises(0.0, 0.2) },
                 proc{|x| Distribution.vonmises(x, 0.0, 0.2) })
  draw_histogram("vonmises-16", 100, 100000, -PI, PI, reporter,
                 proc{ rng.vonmises(0.0, 16) },
                 proc{|x| Distribution.vonmises(x, 0.0, 16) })

  draw_histogram("non_central_chi_square", 100, 100000, 0.0, 24.0, reporter,
                 proc{ rng.non_central_chi_square(2, 1.8) },
                 proc{|x| Distribution.non_central_chi_square(x, 2, 1.8) })
  
  draw_histogram("non_central_t-10-3", 100, 100000, -7.0, 13.0, reporter,
                 proc{ rng.non_central_t(10, 3.0) },
                 proc{|x| Distribution.non_central_t(x, 10, 3.0) })
  draw_histogram("non_central_t-1-1.2", 100, 100000, -7.0, 12.0, reporter,
                 proc{ rng.non_central_t(1, 1.2) },
                 proc{|x| Distribution.non_central_t(x, 1, 1.2) })
  draw_histogram("non_central_t-2-1.2", 100, 100000, -7.0, 12.0, reporter,
                 proc{ rng.non_central_t(2, 1.2) },
                 proc{|x| Distribution.non_central_t(x, 2, 1.2) })
  draw_histogram("non_central_t-2--1.2", 100, 100000, -12.0, 7.0, reporter,
                 proc{ rng.non_central_t(2, -1.2) },
                 proc{|x| Distribution.non_central_t(x, 2, -1.2) })

  draw_histogram("planck", 100, 100000, 0.0, 12.0, reporter,
                 proc{ rng.planck(3, 1) },
                 proc{|x| Distribution.planck(x, 3, 1) })
  
  draw_disc_histogram("bernoulli", 100000, reporter,
                      proc{ rng.bernoulli(0.65) },
                      proc{|x| x == 0 ? 0.35 : 0.65 })

  draw_disc_histogram("geometric", 100000, reporter,
                      proc{ rng.geometric(0.1) },
                      proc{|x| Distribution.geometric(x, 0.1) })
  
  draw_disc_histogram("binomial1-20", 100000, reporter,
                      proc{ rng.binomial(20, 0.45) },
                      proc{|x| Distribution.binomial(x, 20, 0.45) })
  draw_disc_histogram("binomial1-200", 100000, reporter,
                      proc{ rng.binomial(200, 0.65) },
                      proc{|x| Distribution.binomial(x, 200, 0.65) })

  binomial = Random::Binomial.new(rng, 20, 0.45);
  draw_disc_histogram("binomial2-20", 100000, reporter,
                      proc{ binomial.rand },
                      proc{|x| Distribution.binomial(x, 20, 0.45) })
  binomial = Random::Binomial.new(rng, 200, 0.65);
  draw_disc_histogram("binomial2-200", 100000, reporter,
                      proc{ binomial.rand },
                      proc{|x| Distribution.binomial(x, 200, 0.65) })

  draw_disc_histogram("poisson-5", 100000, reporter,
                      proc{ rng.poisson(5.0) },
                      proc{|x| Distribution.poisson(x, 5.0) })
  draw_disc_histogram("poisson-50", 100000, reporter,
                      proc{ rng.poisson(50.0) },
                      proc{|x| Distribution.poisson(x, 50.0) })

  draw_disc_histogram("hypergeometric-50-14-10", 100000, reporter,
                      proc{ rng.hypergeometric(50, 14, 10) },
                      proc{|x| Distribution.hypergeometric(x, 50, 14, 10) })
  draw_disc_histogram("hypergeometric-500-140-100", 100000, reporter,
                      proc{ rng.hypergeometric(500, 140, 100) },
                      proc{|x| Distribution.hypergeometric(x, 500, 140, 100) })
  draw_disc_histogram("hypergeometric-5000-1400-1000", 100000, reporter,
                      proc{ rng.hypergeometric(5000, 1400, 1000) },
                      proc{|x| Distribution.hypergeometric(x, 5000, 1400, 1000) })
  
  draw_disc_histogram("negative_binomial-5-0.2", 100000, reporter,
                      proc{ rng.negative_binomial(5, 0.2) },
                      proc{|x| Distribution.negative_binomial(x, 5, 0.2) })
  draw_disc_histogram("negative_binomial-10-0.2", 100000, reporter,
                      proc{ rng.negative_binomial(10, 0.2) },
                      proc{|x| Distribution.negative_binomial(x, 10, 0.2) })
  draw_disc_histogram("negative_binomial-0.6-0.5", 100000, reporter,
                      proc{ rng.negative_binomial(0.6, 0.5) },
                      proc{|x| Distribution.negative_binomial(x, 0.6, 0.5) })
  draw_disc_histogram("logseries-0.6", 100000, reporter,
                      proc{ rng.logseries(0.6) },
                      proc{|x| Distribution.logseries(x, 0.6) })
  draw_disc_histogram("logseries-0.95", 100000, reporter,
                      proc{ rng.logseries(0.95) },
                      proc{|x| Distribution.logseries(x, 0.95) })
  draw_disc_histogram("zipf_mandelbrot", 100000, reporter,
                      proc{ rng.zipf_mandelbrot(20, 2.1, 1.9) },
                      proc{|x| Distribution.zipf_mandelbrot(x, 20, 2.1, 1.9) })
  draw_histogram("zeta", 20, 100000, 0.5, 20.5, reporter,
                 proc{ rng.zeta(2) },
                 proc{|x| Distribution.zeta(x, 2) })
end

