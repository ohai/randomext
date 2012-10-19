
class Random
  # Draw a random sample from standard normal distribution.
  #
  # Computed using Polar method.
  def standard_normal
    raise
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

  module ZigguratStandardNormal
    module_function
    def f(x); Math.exp(-x**2/2); end
    
    R = 3.442619855899
    V = 9.91256303526217e-3
    K = 7
    M = 64
    N = 2**K
    
    W = Array.new(N)
    KK = Array.new(N)
    F = Array.new(N+1)
    X = Array.new(N)
    
    W[N-1] = V*Math.exp(R**2/2)/(2**(M-K-1))
    W[N-2] = R/(2**(M-K-1))
    KK[N-1] = (R/W[N-1]).ceil
    F[N-1] = Math.exp(-R**2/2)
    X[N] = V*f(R)
    X[N-1] = R
    (N-2).downto(1) do |i|
      X[i] = Math.sqrt(-2*Math.log(f(X[i+1])+V/X[i+1]))
      W[i-1] = X[i]/(2**(M-K-1))
      KK[i] = (X[i]/W[i]).ceil
      F[i] = Math.exp(-X[i]**2/2)
    end
    KK[0] = 0
    F[0] = 1
    
    def generate(rng)
      loop do 
        sign = (rng.rand(2) == 0) ? 1 : -1
        i = rng.rand(N)
        u = rng.rand(2**(M-K-1))
        return sign*u*W[i] if u < KK[i]
        return sign*sample_from_tail(rng) if i == N-1
        ux = u * W[i]
        f = Math.exp(-ux**2/2)
        return sign*ux if rng.rand*(F[i] - F[i+1]) <= f - F[i+1]
      end
    end

    def sample_from_tail(rng)
      loop do
        x = Math.sqrt(R**2-2*Math.log(1-rng.rand))
        return x if x*rng.rand <= R
      end
    end
  end

  def standard_normal
    ZigguratStandardNormal.generate(self)
  end
end
