
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

  # Draw a random sample from normal(Gaussian) distribution.
  #
  # @param [Float] mean mean/average
  # @param [Float] sd Standard deviarion
  def normal(mean=0.0, sd=1.0)
    mean + standard_normal()*sd
  end

  # Draw a random sample from a log normal distribution.
  #
  # The lognormal distribution with parameters mu and sigma
  # is defined
  #   1/sqrt(2*PI*sigma**2)*exp(-(log(x)-mu)**2/(2*sigma**2))
  # @param [Float] mu mean in the normal distribution
  # @param [Float] sigma standard deviarion in the normal distribution
  def lognormal(mu=0.0, sigma=1.0)
    Math.exp(normal(mu, sigma))
  end
end
