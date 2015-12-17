# -*- Ruby -*-

Gem::Specification.new do |spec|
  spec.name = "randomext"
  spec.version = "0.1.1"
  spec.summary = "This library extend Random class of Ruby standard library"
  spec.author = "Ippei Obayashi"
  spec.email = "ohai@kmc.gr.jp"
  spec.description = <<EOS
This library extends class Random in the Ruby standard library.
The Random class in the Ruby standard library supports only 
random sampling from discrete/continuous uniform distribution.

This library provides random sampling methods from 
many kinds of probability distributions such as normal, gamma,
beta, chi_square, t, F, binomial, Poisson, and many other
distributions.
EOS
  spec.files =
    Dir.glob("lib/**/*.rb") +
    Dir.glob("tests/*.rb") +
    ["ext/extconf.rb"] +
    Dir.glob("ext/*.[ch]")
  spec.extensions << 'ext/extconf.rb'
  spec.has_rdoc = false
  spec.homepage = "http://www.kmc.gr.jp/~ohai/randomext/"
  spec.license = "BSD-2-Clause"
end
