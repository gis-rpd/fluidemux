## Custom FastQ demultiplexer developed for analyzing shRNA libraries
## containing a special construct
##
## - Author: Andreas Wilm <wilma@gis.a-star.edu.sg>
## - License: The MIT License


import cligen
import fastq
import tables
#import yaml/serialization, streams
import strutils
import zip/gzipfiles


# from https://bitbucket.org/andreas-wilm/decontanimate/
type NucAlphabet* = enum
  ## Simple enum for nucleid acids alphabets
  DNA, RNA


proc complement*(nucleicAcid: string, toAlphabet: NucAlphabet): string =
  ## Returns complement of a nucleic acid. Case aware.
  ## IUPAC ambiguity characters will raise ValueError.
  result = newString(nucleicAcid.len)# FIXME would it be faster to copy input and change in place?
  let ignChars = ['-', 'N', 'n']
  var trans = {'C': 'G', 'c': 'g',
                'G': 'C', 'g': 'c',
                'U': 'A', 'u': 'a',
                'T': 'A', 't': 'a'
              }.toTable
  if toAlphabet == DNA:
      trans['A'] = 'T'
      trans['a'] = 't'
  elif toAlphabet == RNA:
      trans['A'] = 'U'
      trans['a'] = 'u'
  else:
    raise newException(ValueError, "Unknown Alphabet " & $toAlphabet)

  for i, nuc in nucleicAcid:
    if nuc in ignChars:
      result[i] = nuc
    else:
      doAssert trans.hasKey(nuc)
      result[i] = trans[nuc]


proc reversed(s: string): string =
  ## Return reversed copy of input string.
  ## From https://www.rosettacode.org/wiki/Reverse_a_string
  result = newString(s.len)
  for i, c in s:
    result[s.high - i] = c


proc parseBarcodeYaml(barcodeYaml: string): Table[string, string] =
  ## parse barcode to name mapping from given yaml file
  result = initTable[string, string]()
  #type bc2name = object
  #  bc: string
  #  name: string
  #var bc2nameList: seq[bc2name]
  #var s = newFileStream(barcodeYaml)
  #load(s, bc2nameList)
  #s.close()
  #for x in bc2nameList:
  #  result[x.bc] = x.name
  for line in lines barcodeYaml:
    let fields = line.split(':', 1)
    let bc = fields[0].strip()
    let name = fields[1].strip()
    doAssert len(bc) == 6 and name.startsWith("ROW")
    result[bc] = name


proc main(fastq1: string, fastq2: string, barcodeYaml: string, outpref: string) =
  var bcMap = parseBarcodeYaml(barcodeYaml)
  var bcClassCounts = initCountTable[string]()
  var bcCounts = initCountTable[string]()
  var validBarcodes: seq[string]
  var fastqOutStreams = initTable[string, (GzFileStream, GzFileStream)]()
  for k in keys bcMap:
    validBarcodes.add(k)

  echo "DEBUG ", bcMap;

  for r1, r2 in parsePairedEnd(fastq1, fastq2):
    var bc1 = r1.sequence[0..5]
    var bc2 = r2.sequence[r2.len-6..<r2.len]
    var class = ""
    var bc = ""
    bc1 = complement(reversed(bc1), DNA)

    # classify barcodes 
    if bc1 in validBarcodes and bc2 in validBarcodes:
      if bc1 == bc2:
        class = "both-valid"
        bc = bc1
      else:
        class = "both-valid-but-different"
    elif bc1 in validBarcodes or bc2 in validBarcodes:
      class = "one-valid"
      if bc1 in validBarcodes:
        bc = bc1
      elif bc2 in validBarcodes:
        bc = bc2
      else:
       raise newException(ValueError, "Internal error")
    elif "N" in bc1 or "N" in bc2:
      class = "N-in-either"
    else:
      class = "invalid-no-Ns"
    if bcClassCounts.hasKey(class) == false:
      bcClassCounts[class] = 0
    bcClassCounts.inc(class)

    # write only the ones for which we have exactly one or two matches
    if class == "both-valid" or class == "one-valid":
      var bcName: string
      bcName = bcMap[bc]
      assert bcName != ""
      if bcCounts.hasKey(bcName) == false:
        bcCounts[bcName] = 0
      bcCounts.inc(bcName)

      var clippedR1 = r1[6..<r1.len]
      var clippedR2 = r2[0..<r2.len-6]
      assert len(clippedR1)+6 == len(r1)
      if fastqOutStreams.haskey(bcname) == false:
        let fq1Stream = newGZFileStream(outpref & bcname & ".R1.fastq.gz", fmWrite)
        let fq2Stream = newGZFileStream(outpref & bcname & ".R2.fastq.gz", fmWrite)
        fastqOutStreams[bcname] = (fq1Stream, fq2Stream)

      let (fq1Stream, fq2Stream) = fastqOutStreams[bcname]
      writeRecord(fq1Stream, clippedR1)
      writeRecord(fq2Stream, clippedR2)

  for k in keys bcClassCounts:
    echo "Pairs in class ", k, ": ", bcClassCounts[k]
  for k in keys bcCounts:
    echo "Pairs for barcode ", k, ": ", bcCounts[k]

  # close all streams
  for bcName in keys fastqOutStreams:
    let (fq1Stream, fq2Stream) = fastqOutStreams[bcname]
    fq1Stream.close
    fq2Stream.close


when isMainModule:
  import cligen
  dispatch(main, 
    help = {
      "fastq1":  "FastQ1",
      "fastq2" : "FastQ2",
      "barcodeYaml": "YAML file with barcodes to name mapping",
      "outpref": "Output prefix"
    },
    short = {
          "fastq1": '1',
          "fastq2": '2'})
