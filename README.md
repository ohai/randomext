# randomext
The library extends class Random in the Ruby standard library.

## Overview
The Random class in the Ruby standard library supports only 
random sampling from discrete/continuous uniform distribution.

This library provides random sampling methods from 
many kinds of probability distributions such as:

* normal (Gaussian)
* lognormal
* Cauthy
* levy
* exponential
* Laplace
* Rayleigh
* Weibull
* Gumbel
* gamma
* beta
* power
* chi_square
* F
* t
* Wald (inverse Gaussian)
* Pareto
* logistic
* von Mises
* Non-Central Chi-Square
* Non-Central t
* Planck
* Bernoulli
* binomial
* Poisson
* geometric
* negative binomial
* log series
* Zipf-Mandelbrot
* zeta

## Usage
To use this library, you need to install randomext gem:

    gem install randomext

And write 

    require 'randomext'

in your ruby script, then you can use some additional methods in Random class.

## References
Almost all algorithms are based on:
四辻哲章, "計算機シミュレーションのための確率分布乱数生成法", プレアデス出版 (2010)

I examine numpy to select nice distributions.