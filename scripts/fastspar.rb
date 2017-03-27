class Fastspar < Formula
  desc "A fast c++ implementation of SparCC"
  homepage "https://github.com/scwatts/fastspar"
  url "https://github.com/scwatts/fastspar/archive/v0.0.3.tar.gz"
  sha256 "5b9aab1e3a4e9935868a9b43cd93095d20daf37f9da569b51c4374af4cf24564"
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
