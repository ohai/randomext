
class Random
  # Draw a random sample from standard normal distribution.
  #
  # Computed using Polar method.
  def standard_normal
    if @z
      r = @z; @z = nil; return r
    end
    
    begin
      v1 = self.rand(-1.0 ... 1.0)
      v2 = self.rand(-1.0 ... 1.0)
      v = v1**2 + v2**2
    end until 0.0 < v && v < 1.0
    w = Math.sqrt(-2.0*Math.log(v)/v)
    @z = v1 * w
    return v2 * w
  end
end
