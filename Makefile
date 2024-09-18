CHR1 = ../data/t2t-chr1.txt
#CHR1 = /xdisk/twheeler/data/genomes/human/t2t/chr1.txt
SEQ1 = tests/inputs/1.fa
SEQ2 = tests/inputs/2.fa
SEQ3 = tests/inputs/2.fa
SUFR   = ./target/release/sufr

CREATE_DEBUG = cargo run -- create --log debug

CREATE = cargo run -- create --log info
#CREATE = $(SUFR) create --log info

#READ   = cargo run -- check
#CHECK  = cargo run -- read

#CHECK  = $(SUFR) read
#READ   = $(SUFR) check

perf:
	perf record --call-graph dwarf $(SUFR) create -t 16 --log info $(CHR1)

s1:
	$(CREATE_DEBUG) $(SEQ1) -n 2 --check

s2:
	$(CREATE_DEBUG) $(SEQ2) --check

s3:
	$(CREATE_DEBUG) $(SEQ3) --check

ecoli:
	$(CREATE) ../data/ecoli.fa -o ecoli.sa --check

chr1:
	$(CREATE) ../data/t2t-chr1.fa -n 64

check-t2t-chr1:
	$(READ) -s ../data/t2t-chr1.txt -a t2t-chr1.sa

create-chr1:
	$(CREATE) ../data/chr1.txt -o t2t-chr1.sa --ignore-start-n

check-chr1:
	$(READ) -s ../data/chr1.txt -a t2t-chr1.sa

valcache:
	valgrind --tool=cachegrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info

valcall:
	valgrind --tool=callrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info
