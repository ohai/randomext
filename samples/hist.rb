require 'randomext'
require 'tempfile'
require 'benchmark'
require 'bigdecimal'

include Math


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

def draw_histogram(name, num_bins, num_samples, min, max, reporter, gen, func=nil)
  return if !ARGV.empty? && !ARGV.include?(name)
  
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
  return if !ARGV.empty? && !ARGV.include?(name)
  
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

  def gamma(x, alpha, beta)
    beta**(-alpha)*x**(alpha-1)*exp(-x/beta)/Math.gamma(alpha)
  end

  def beta(x, alpha, beta)
    b = Math.gamma(alpha)*Math.gamma(beta)/Math.gamma(alpha+beta)
    x**(alpha-1)*(1-x)**(beta-1)/b
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
    combination(n, x) * BigDecimal(p,14)**x * BigDecimal(1-p,14)**(n-x)
  end

  def poisson(x, lambda)
    exp(-(1..x).inject(0){|r, i| r+log(i)} - lambda + x*log(lambda))
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
  draw_disc_histogram("bernoulli", 100000, reporter,
                      proc{ rng.bernoulli(0.65) },
                      proc{|x| x == 0 ? 0.35 : 0.65 })

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
end


