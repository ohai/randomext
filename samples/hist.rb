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

def normal_distribution(x, average = 0, sd = 1.0)
  exp(-(x - average)**2/(2 * sd**2))/sqrt(2*PI*sd**2)
end

rng = Random.new
draw_histogram("snormal.png", 80, 100000, -6.0, 6.0, rng.method(:standard_normal),
               method(:normal_distribution))
               
draw_histogram("normal.png", 80, 100000, -3.0, 9.0,
               proc{ rng.normal(3.0, 1.7) },
               proc{|x| normal_distribution(x, 3.0, 1.7) })
