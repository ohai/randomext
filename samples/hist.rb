require 'randomext'
require 'tempfile'

include Math


def histogram(num_bins, num_samples, min, max, gen)
  binsize = (max - min).to_f / num_bins
  bins = Array.new(num_bins, 0)
  datasize = 0
  
  num_samples.times do
    n = (gen[] - min).div(binsize)
    if 0 <= n && n < num_bins
      datasize += 1
      bins[n] += 1
    end
  end

  return bins.map.with_index{|k, n|
    [binsize*(n+0.5)+min, k/(binsize*datasize.to_f)]
  }
end

def plot(filename, data, curve)
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
    process.puts("plot \"#{hist_tmp.path}\" w boxes, \"#{curve_tmp.path}\" w l")
    process.flush
  end
ensure
  hist_tmp.close(true)
  curve_tmp.close(true)
end

def draw_histogram(filename, num_bins, num_samples, min, max, gen, func=nil)
  hist = histogram(num_bins, num_samples, min, max, gen)
  if func
    curve = hist.map{|x, _| [x, func[x]] }
  else
    curve = []
  end

  plot(filename, hist, curve)
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
end

rng = Random.new
draw_histogram("snormal.png", 80, 100000, -6.0, 6.0,
               rng.method(:standard_normal), Distribution.method(:normal))
               
               
draw_histogram("normal.png", 80, 100000, -3.0, 9.0,
               proc{ rng.normal(3.0, 1.7) },
               proc{|x| Distribution.normal(x, 3.0, 1.7) })

draw_histogram("lognormal.png", 100, 100000, 0.0, 8.0,
               rng.method(:lognormal), Distribution.method(:lognormal))

draw_histogram("cauthy.png", 100, 100000, -40.0, 46.0,
               proc{ rng.cauthy(3.0, 1.5) },
               proc{|x| Distribution.cauthy(x, 3.0, 1.5) })
