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
* Chi-square
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

## Example
Create an sequence of random numbers from the Gaussian distribution with
the mean is 0.0 and S.D. is 2.0.

    require 'randomext'
    random_numbers = Array.new(100){ Random::DEFAULT.normal(0.0, 2.0) }


## References
Almost all algorithms are based on:
四辻哲章, "計算機シミュレーションのための確率分布乱数生成法", プレアデス出版 (2010)

I examine numpy to select nice distributions.

## URLs
* [Web site](http://www.kmc.gr.jp/~ohai/randomext/)
* [Repository and issue tracker](https://bitbucket.org/ohai/randomext)

## Author
Ippei Obayashi <ohai@kmc.gr.jp>

## Copyright
Copyright (c) 2012, Ippei Obayashi
All rights reserved.

The software is licensed under the BSD 2-Clause License.
Please see the {file:LICENSE} for more information.
