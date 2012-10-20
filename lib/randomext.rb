require 'randomext_native'

class Random
  # Draw a random sample from a normal(Gaussian) distribution.
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
  # @param [Float] loc location parameter
  # @param [Float] scale scale parameter
  def cauthy(loc, scale)
    loc + scale*standard_cauthy()
  end

  # Draw a random sample from the standard Cauthy distribution.
  #
  # Computed using Polar method from the standard normal distribution.
  def standard_cauthy
    y1 = standard_normal()
    begin; y2 = standard_normal(); end until y2 != 0.0
    return y1/y2
  end

  # Draw a random sample from a Levy distribution.
  #
  # @param loc location parameter
  # @param scale scale parameter
  def levy(loc=0.0, scale=1.0)
    begin z = standard_normal.abs; end until z > 0
    loc + scale/z**2
  end

  # Draw a random sample from a exponential distribution.
  #
  # Inverse function method is used.
  # @param scale scale parameter
  def exponential(scale=1.0)
    -scale * Math.log(1-rand)
  end

  # Draw a random sample from a gamma distribution
  #
  # @param shape shape parameter
  # @param scale scale parameter
  def gamma(shape, scale=1.0)
    case
    when shape <= 0.0
      raise ArgumentError, "Random#gamma: shape parameter should be positive"
    when shape > 1.0
      scale * _gamma(shape)
    when shape == 1.0
      exponential(scale)
    when shape < 1.0
      scale*_gamma(shape+1)*rand_open_interval**(1.0/shape)
    end
  end

  # Draw a sample from a beta distribution.
  def beta(alpha, beta)
    y1 = gamma(alpha); y2 = gamma(beta)
    y1/(y1+y2)
  end
  
  # Draw a sample from the uniform distribution on (0, 1)
  def rand_open_interval
    begin; x = rand; end until x != 0.0
    x
  end
end
