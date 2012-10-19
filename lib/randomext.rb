
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
  # @param [Float] mean mean
  # @param [Float] sd standard deviarion
  def normal(mean=0.0, sd=1.0)
    mean + standard_normal()*sd
  end

  # Draw a random sample from a log normal distribution.
  #
  # The lognormal distribution with parameters mu and sigma
  # is defined:
  #   1/sqrt(2*PI*sigma**2)*exp(-(log(x)-mu)**2/(2*sigma**2))
  # @param [Float] mu the mean in a normal distribution
  # @param [Float] sigma the standard deviarion in a normal distribution
  def lognormal(mu=0.0, sigma=1.0)
    Math.exp(normal(mu, sigma))
  end

  # Draw a random sample from a Cauthy distribution.
  #
  # @param [Float] mu
  # @param [Float] theta
  def cauthy(mu, theta)
    mu + theta*standard_cauthy()
  end

  # Draw a random sample from the standard Cauthy distribution.
  #
  # Computed using Polar method from the standard normal distribution.
  def standard_cauthy
    y1 = standard_normal()
    begin; y2 = standard_normal(); end until y2 != 0.0
    return y1/y2
  end
end
