version       = "0.0.1"
author        = "Andreas Wilm and others"
description   = "Say Li's stuff"
license       = "MIT"

requires "nim >= 0.19"
requires "cligen >= 0.9.16"
requires "zip >= 0.2.1"
requires "yaml >= 0.11.0"

srcDir = "src"

bin = @["saylidemux"]

skipDirs = @["tests"]
skipExt = @["nim"]

task test, "run tests":
  withDir "tests":
    exec "nim c --lineDir:on --debuginfo -r all"
