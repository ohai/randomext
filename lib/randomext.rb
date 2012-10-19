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
