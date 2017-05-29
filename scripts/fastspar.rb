class Fastspar < Formula
  desc "A fast c++ implementation of SparCC"
  homepage "https://github.com/scwatts/fastspar"
  url "https://github.com/scwatts/fastspar/archive/v0.0.4.tar.gz"
  sha256 "64db95f37251031cfb757590477dc4d7545fa213ee08f0dc9560e72490cfb5f0"
  head "https://github.com/scwatts/fastspar.git"

  depends_on "armadillo" => :run
  depends_on "gsl" => :run
  depends_on "gnu-getopt" => :run
  depends_on "openblas" => :run

  needs :openmp

  def install
    system "./configure", "CC=gcc-6",
                          "CXX=g++-6",
                          "LDFLAGS=-L#{HOMEBREW_PREFIX}/lib/openblas/",
                          "--prefix=#{prefix}"
    system "make", "-j", "install"
  end

  test do
    system "false"
  end
end
