class Fastspar < Formula
  desc "A fast c++ implementation of SparCC"
  homepage "https://github.com/scwatts/fastspar"
  url "https://github.com/scwatts/fastspar/archive/v0.0.3.tar.gz"
  sha256 "7dc39fd5c40132f17002d61219231e8d1a950cb8592e418c2464f442c26109f1"
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
