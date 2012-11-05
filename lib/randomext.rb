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
  # @param [Float] loc location parameter
  # @param [Float] scale scale parameter
  def levy(loc=0.0, scale=1.0)
    begin z = standard_normal.abs; end until z > 0
    loc + scale/z**2
  end

  # Draw a random sample from a exponential distribution.
  #
  # Inverse function method is used.
  # @param [Float] scale scale parameter
  def exponential(scale=1.0)
    if scale < 0.0
      raise ArgumentError, "Random#exponential: scale parameter must be positive"
    end
    scale * standard_exponential
  end

  # Draw a random sample from a gamma distribution
  #
  # @param [Float] shape shape parameter
  # @param [Float] scale scale parameter
  def gamma(shape, scale=1.0)
    if scale <= 0.0
      raise ArgumentError, "Random#gamma: scale parameter must be positive"
    end
    
    case
    when shape <= 0.0
      raise ArgumentError, "Random#gamma: shape parameter must be positive"
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

  # Draws a random sample from a power function distribution
  #
  # @param [Float] shape shape parameter
  # @param [Float] a lower boundary parameter
  # @param [Float] b upper boundary parameter
  def power(shape, a, b)
    if shape <= 0 || a >= b
      raise ArgumentError, "Random#power: shape must be positive, and b should be greater than a"
    end
    
    a + (b-a)*(rand_open_interval**(1/shape))
  end
  # Draw a random sample from a chi_square distribution.
  #
  # @param [Integer] r degree of freedom
  def chi_square(r)
    if r == 1
      standard_normal ** 2
    elsif r > 1
      gamma(r*0.5, 2)
    else
      raise ArgumentError, "Random#chi_square:r (degree of freedom) must be >= 1"
    end
  end

  # Draws a random sample from a F distribution.
  #
  # @param [Integer] r1 degree of freedom
  # @param [Integer] r2 degree of freedom
  def F(r1, r2)
    f = r2 / r1.to_f
    f*chi_square(r1)/chi_square(r2)
  end

  # Draws a random sample from a t distribution.
  #
  # @param [Integer] r degree of freedom
  def t(r)
    if r ==1
      standard_cauthy
    elsif r == 2
      standard_normal/Math.sqrt(exponential(1))
    elsif r > 2
      rdiv2 = r/2.0
      Math.sqrt(rdiv2)*standard_normal/Math.sqrt(_gamma(rdiv2))
    else
      raise ArgumentError, "Random#t: r (degree of freedom) must be >= 1"
    end
  end

  # Draws a random sample from a Pareto distribution.
  #
  # The probabilistic mass function for the Pareto distribution
  # with parameters a and b is defined as:
  #   p(x) = a*b**a/x**(a+1)
  #
  # @param [Float] a shape parameter
  # @param [Float] b scale parameter
  def pareto(a, b=1.0)
    if a <= 0 || b <= 0
      raise ArgumentError, "Random#pareto: parameters a and b must be positive"
    end
    b * (1.0 - rand)**(-1/a)
  end
  
  # Draws a random sample from a logistic distribution.
  #
  # @param [Float] mu the location parameter
  # @param [Float] theta the scale parameter
  def logistic(mu, theta)
    u = rand_open_interval
    mu + theta*log(u/(1-u))
  end
  
  # Draw a sample from Bernoulli distribution.
  #
  # @param [Float] p the probability returning 1 
  def bernoulli(p)
    (rand < p) ? 1 : 0
  end

  # Draws a random sample from a geometric distribution.
  #
  # @param [Float] theta the probability of sucess
  def geometric(theta)
    if theta <= 0.0 || theta >= 1.0
      raise ArgumentError, "Random#geometric: theta should be in (0, 1)"
    end
    
    d= -1/(Math.log(1-theta))
    (d * standard_exponential).floor + 1
  end

  # Draws a random sample from a negative binomial distribution.
  #
  # @param [Float] r positive parameter
  # @param [Float] theta parameter in (0, 1)
  def negative_binomial(r, theta)
    if r <= 0.0
      raise ArgumentError, "Random#negative_binomial: r must be positive"
    end
    if theta <= 0.0 && theta >= 1.0
      raise ArgumentError, "Random#negative_binomial: theta must be in (0, 1)"
    end
    poisson(gamma(r, 1/theta - 1))
  end
  
  # Draw a sample from the uniform distribution on (0, 1)
  def rand_open_interval
    begin; x = rand; end until x != 0.0
    x
  end
end
