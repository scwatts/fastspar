class Fastspar < Formula
  desc " Rapid and scalable correlation estimation for compositional data"
  homepage "https://github.com/scwatts/fastspar"
  url "https://github.com/scwatts/fastspar/archive/v0.0.5.tar.gz"
  sha256 "c6fafdff151fa3dcf430fb4c593d392ac27f721807ac4f49d3412c96d25a8389"
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
