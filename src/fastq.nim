## Nim routines for processing FastQ files

## Adopted from Jonathan Badger's nimbioseq: https://github.com/jhbadger/nimbioseq
## and https://bitbucket.org/andreas-wilm/decontanimate
##
#import sequtils, strutils, math, tables, osproc, streams
import zip/gzipfiles
import sequtils
import strutils


type Record* = object
  ## This type represents a genetic sequence with optional quality
  id*: string
  description*: string
  qualities*: string
  sequence*: string


proc len*(rec: Record): int =
  return len(rec.sequence)


proc `$`*(rec: Record): string =
  ## Returns formatted string of FastQ record ``rec``.
  ## Borrowed from `jhbadger's bioseq <https://github.com/jhbadger/nimbioseq.git>`_
  var header = "@" & rec.id
  if rec.description != "":
    header = header & " " & rec.description
  header & "\n" & rec.sequence & "\n+\n" & rec.qualities


proc `[]`*(srcRecord: Record, x: HSlice): Record =
  result = srcRecord
  result.qualities = result.qualities[x]
  result.sequence = result.sequence[x]


proc parseNextRecord*(fqStream: GZFileStream): Record =
  var line: string

  if readLine(fqStream, line) == false:
    return
  var fields = split(line, ' ', 1)
  if len(fields) > 1:
    (result.id, result.description) = fields
  else:
    (result.id, result.description) = (fields[0], "")
  assert result.id[0] == '@'
  result.id = result.id[1..<len(result.id)]

  if readLine(fqStream, line) == false:
    raise newException(ValueError, "Premature end of file")
  result.sequence = line

  if readLine(fqStream, line) == false:
    raise newException(ValueError, "Premature end of file")
  doAssert line[0] == '+'

  if readLine(fqStream, line) == false:
    raise newException(ValueError, "Premature end of file")
  result.qualities = line
  doAssert len(result.qualities) == len(result.sequence)


# FIXME would be easier if there was zip support for iterators
iterator parsePairedEnd*(fastqGz1: string, fastqGz2: string): (Record, Record) =
  var rec1, rec2: Record
  var fqStream1 = newGZFileStream(fastqGz1, fmRead)
  var fqStream2 = newGZFileStream(fastqGz2, fmRead)
  var emptyRec: Record

  while true:
    rec1 = parseNextRecord(fqStream1)
    rec2 = parseNextRecord(fqStream2)
    if rec1 == emptyRec:
      doAssert rec2 == emptyRec
      break
    else:
      doAssert rec2 != emptyRec
    yield (rec1, rec2)


proc writeRecord*[T](fh: T, rec: Record) =
  fh.write($rec & "\n")

