require 'randomext_native'

class Random
  # Draws a random sample from a normal(Gaussian) distribution.
  #
  # @param [Float] mean mean
  # @param [Float] sd standard deviarion
  # @return [Float] a random sample
  def normal(mean=0.0, sd=1.0)
    mean + standard_normal()*sd
  end

  # Draws a random sample from a log normal distribution.
  #
  # The lognormal distribution with parameters mu and sigma
  # is defined:
  #   1/sqrt(2*PI*sigma**2)*exp(-(log(x)-mu)**2/(2*sigma**2))
  # @param [Float] mu the mean in a normal distribution
  # @param [Float] sigma the standard deviarion in a normal distribution
  # @return [Float] a random sample in (0, INFINITY)
  def lognormal(mu=0.0, sigma=1.0)
    Math.exp(normal(mu, sigma))
  end

  # Draws a random sample from a Cauthy distribution.
  #
  # @param [Float] loc location parameter
  # @param [Float] scale scale parameter
  # @return [Float] a random sample
  def cauthy(loc, scale)
    loc + scale*standard_cauthy()
  end

  # Draws a random sample from the standard Cauthy distribution.
  #
  # Computed using Polar method from the standard normal distribution.
  # @return [Float] a random sample
  def standard_cauthy
    y1 = standard_normal()
    begin; y2 = standard_normal(); end until y2 != 0.0
    return y1/y2
  end

  # Draws a random sample from a Levy distribution.
  #
  # @param [Float] loc location parameter
  # @param [Float] scale scale parameter
  # @return [Float] a random sample
  def levy(loc=0.0, scale=1.0)
    begin z = standard_normal.abs; end until z > 0
    loc + scale/z**2
  end

  # Draws a random sample from a exponential distribution.
  #
  # Inverse function method is used.
  # @param [Float] scale scale parameter (scale > 0)
  # @return [Float] a random sample
  def exponential(scale=1.0)
    if scale < 0.0
      raise ArgumentError, "Random#exponential: scale parameter must be positive"
    end
    scale * standard_exponential
  end

  # Draws a random sample from a Laplace distribution
  #
  # @param [Float] loc location parameter
  # @param [Float] scale scale parameter
  # @return [Float] a random sample
  def laplace(loc=0.0, scale=1.0)
    sign = rand(2) == 1 ? 1 : -1
    loc + sign*scale*standard_exponential
  end

  # Draws a random sample from a Rayleigh distribution
  #
  # @param [Float] sigma scale parameter
  # @return [Float] a random sample
  def rayleigh(sigma=1.0)
    sigma*Math.sqrt(2*standard_exponential)
  end

  # Draws a random sample from a Weibull distribution
  #
  # @param [Float] g shape parameter (g > 0.0)
  # @param [Float] mu scale parameter 
  # @return [Float] a random sample
  def weibull(g, mu=1.0)
    if g <= 0
      raise ArgumentError, "Random#weibull: shape parameter must be positive"
    end
    mu * standard_exponential**(1.0/g)
  end

  # Draws a random sample from a Gumbel distribution
  #
  # @param [Float] loc location parameter
  # @param [Float] scale scale parameter
  # @return [Float] a random sample
  def gumbel(loc=0.0, scale=1.0)
    loc - scale * Math.log(standard_exponential)
  end
  
  # Draws a random sample from a gamma distribution
  #
  # @param [Float] shape shape parameter (shape > 0.0)
  # @param [Float] scale scale parameter (scale > 0.0)
  # @return [Float] a random sample
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

  # Draws a random sample from a beta distribution.
  #
  # @param [Float] alpha a shape parameter (alpha > 0.0)
  # @param [Float] beta another shape parameter (beta > 0.0)
  # @return [Float] a random sample
  def beta(alpha, beta)
    y1 = gamma(alpha); y2 = gamma(beta)
    y1/(y1+y2)
  end

  # Draws a random sample from a power function distribution
  #
  # @param [Float] shape shape parameter (shape > 0.0)
  # @param [Float] a lower boundary parameter
  # @param [Float] b upper boundary parameter (a < b)
  # @return [Float] a random sample
  def power(shape, a, b)
    if shape <= 0 || a >= b
      raise ArgumentError, "Random#power: shape must be positive, and b should be greater than a"
    end
    
    a + (b-a)*(rand_open_interval**(1/shape))
  end
  # Draw a random sample from a chi_square distribution.
  #
  # @param [Integer] r degree of freedom (r >= 1)
  # @return [Float] a random sample
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
  # @param [Integer] r1 degree of freedom (r1 >= 1)
  # @param [Integer] r2 degree of freedom (r2 >= 1)
  # @return [Float] a random sample
  def F(r1, r2)
    f = r2 / r1.to_f
    f*chi_square(r1)/chi_square(r2)
  end

  # Draws a random sample from a t distribution.
  #  
  # @param [Integer] r degree of freedom (r >= 1)
  # @return [Float] a random sample
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

  # Draws a random sample from a wald distribution.
  #
  # A wald distribution is also called an inverse Gaussian distribution.
  #
  # @param [Float] mean mean
  # @param [Float] shape shape parameter (shape > 0.0)
  # @return [Float] a random sample in (0, INFINITY)
  def wald(mean, shape)
    if shape <= 0.0
      raise ArgumentError, "Random#wald: shape parameter must be positive"
    end
    p = mean**2
    q = p/(2*shape)
    z = standard_normal
    return mean if z == 0.0
    v = mean + q*z**2
    x1 = v + Math.sqrt(v**2-p)
    return x1 if rand*(x1 + mean) <= mean
    return p/x1
  end
  
  # Draws a random sample from a Pareto distribution.
  #
  # The probabilistic mass function for the distribution is defined as:
  #   p(x) = a*b**a/x**(a+1)
  #
  # @param [Float] a shape parameter (a > 0.0)
  # @param [Float] b scale parameter (b > 0.0)
  # @return [Float] a random sample in [b, INFINITY)
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
  # @return [Float] a random sample
  def logistic(mu, theta)
    u = rand_open_interval
    mu + theta*log(u/(1-u))
  end

  def non_central_t(r, lambda)
    if lambda == 0.0
      raise ArgumentError, "Random#non_central_t: lambda must not be 0"
    end

    if r == 1
      z = standard_normal + lambda
      w = standard_normal.abs
      z/w
    elsif r == 2
      z = standard_normal + lambda
      w = standard_exponential
      z/Math.sqrt(w)
    elsif r > 2
      d = Math.sqrt(r/2.0)
      z = standard_normal + lambda
      w = _gamma(r/2.0)
      d*z/Math.sqrt(w)
    else
      raise ArgumentError, "Random#non_central_t: r must be positive"
    end
  end
  
  # Draw a random sample from a Bernoulli distribution.
  #
  # @param [Float] p the probability returning 1
  # @return [Integer] a random sample, 0 or 1
  def bernoulli(p)
    (rand < p) ? 1 : 0
  end

  # Draws a random sample from a geometric distribution.
  #
  # @param [Float] theta the probability of sucess (0 < theta < 1)
  # @return [Integer] a random sample in [1, INFINITY)
  def geometric(theta)
    if theta <= 0.0 || theta >= 1.0
      raise ArgumentError, "Random#geometric: theta should be in (0, 1)"
    end
    
    d= -1/(Math.log(1-theta))
    (d * standard_exponential).floor + 1
  end

  # Draws a random sample from a negative binomial distribution.
  #
  # @param [Float] r s parameter (0 < r)
  # @param [Float] theta a parameter (0 < theta < 1)
  # @return [Integer] a random sample in [0, INFINITY)
  def negative_binomial(r, theta)
    if r <= 0.0
      raise ArgumentError, "Random#negative_binomial: r must be positive"
    end
    if theta <= 0.0 && theta >= 1.0
      raise ArgumentError, "Random#negative_binomial: theta must be in (0, 1)"
    end
    poisson(gamma(r, 1/theta - 1))
  end

  # Draws a random sample from a log series distribution.
  #
  # @param [Float] theta a parameter (0 < theta < 1)
  # @return [Integer] a random sample in [0, INFINITY)
  def logseries(theta)
    if theta <= 0 || 1 <= theta
      raise ArgumentError, "Random#logseries: theta must be in (0, 1)"
    end
    q = 1 - theta
    v = rand_open_interval
    if v >= theta
      1
    else
      u = rand_open_interval
      (log(v)/log(1-q**u)).ceil
    end
  end
  
  # Draw a sample from the uniform distribution on (0, 1)
  #
  # @return [Float] a random sample in (0, 1)
  def rand_open_interval
    begin; x = rand; end until x != 0.0
    x
  end
end
