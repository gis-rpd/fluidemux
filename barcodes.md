# Fluidigm Barcodes

- See https://bitbucket.org/andreas-wilm/mrnaseqht_demux/src/master/mRNASeqHT_demultiplex.pl
- From File: mRNASeqHT_demultiplex, Date: 2015-08-19, Version: 1.0.2


CACGTA	=>	'ROW01'	,
	CTCACA	=>	'ROW02'	,
	TGCATC	=>	'ROW03'	,
	TCAGAC	=>	'ROW04'	,
	CGATGT	=>	'ROW05'	,
	TACTGC	=>	'ROW06'	,
	ATGCTC	=>	'ROW07'	,
	CATCTG	=>	'ROW08'	,
	GACTCA	=>	'ROW09'	,
	AGATCG	=>	'ROW10'	,
	ATCAGC	=>	'ROW11'	,
	GCTACA	=>	'ROW12'	,
	CAGATC	=>	'ROW13'	,
	CACAGT	=>	'ROW14'	,
	TACGAG	=>	'ROW15'	,
	CGACTA	=>	'ROW16'	,
	GCATCT	=>	'ROW17'	,
	AGCACT	=>	'ROW18'	,
	ACACTG	=>	'ROW19'	,
	CGTAGA	=>	'ROW20'	,
	ACTCGA	=>	'ROW21'	,
	ACATGC	=>	'ROW22'	,
	CGTCAT	=>	'ROW23'	,
	TGTACG	=>	'ROW24'	,
	GCAGTA	=>	'ROW25'	,
	TCACGT	=>	'ROW26'	,
	ACGTCA	=>	'ROW27'	,
	CTCGAT	=>	'ROW28'	,
	ATCGTG	=>	'ROW29'	,
	GCTGAT	=>	'ROW30'	,
	GTCTAC	=>	'ROW31'	,
	CATGCT	=>	'ROW32'	,
	TAGCAC	=>	'ROW33'	,
	GTGCAT	=>	'ROW34'	,
	TAGTCG	=>	'ROW35'	,
	TCTCAG	=>	'ROW36'	,
	CTAGTC	=>	'ROW37'	,
	TCTAGC	=>	'ROW38'	,
	ATGACG	=>	'ROW39'	,
	GAGCTA	=>	'ROW40'	

[A[A
	 awk '/=/ {print $1}' barcodes.md  > barcodes.txt
