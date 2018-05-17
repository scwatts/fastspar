class Fastspar < Formula
  desc " Rapid and scalable correlation estimation for compositional data"
  homepage "https://github.com/scwatts/fastspar"
  url "https://github.com/scwatts/fastspar/archive/v0.0.6.tar.gz"
  sha256 "07dc76479d4a8aa30ba6d69e20b50241bfbbcdd2f8fd2858b9f5b594bc2cba98"
  head "https://github.com/scwatts/fastspar.git"

  depends_on "autoconf" => :build
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
