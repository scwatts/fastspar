class Fastspar < Formula
  desc " Rapid and scalable correlation estimation for compositional data"
  homepage "https://github.com/scwatts/fastspar"
  url "https://github.com/scwatts/fastspar/archive/v0.0.7.tar.gz"
  sha256 "af401ac1bc2894aa8aff136aa42066b8b309a40b1bca1bfe6ec119d46f5c7ddf"
  head "https://github.com/scwatts/fastspar.git"

  depends_on "autoconf" => :build
  depends_on "autoconf-archive" => :build
  depends_on "automake" => :build
  depends_on "armadillo"
  depends_on "gsl"
  depends_on "gnu-getopt"
  depends_on "openblas"

  needs :openmp

  def install
    system "./autogen.sh"
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
