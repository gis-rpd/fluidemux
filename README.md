# rpd-fluidemux

Project specific FastQ demultiplexer. Of limited general use, but existing code can be easily
changed according to needs. Contains plenty of hardcoded assumptions, e.g. barcodes are of length 6
and names start with "ROW".  Provided barcodes come from the [Fluidigm C1 mRNA Seq HT Demultiplex
Script](https://www.fluidigm.com/software).

## Building

Install nim. Then run `nimble build`

## Usage

```
$ ./fluidemux  -h
Usage:
  main [required&optional-params]
  Options(opt-arg sep :|=|spc):
  -h, --help                             write this help to stdout
  -1=, --fastq1=       string  REQUIRED  FastQ1
  -2=, --fastq2=       string  REQUIRED  FastQ2
  -b=, --barcodeYaml=  string  REQUIRED  YAML file with barcodes to name mapping
  -o=, --outpref=      string  REQUIRED  Output prefix
```


