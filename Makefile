CHR1 = ../data/chr1.fa
SEQ1 = sufr/tests/inputs/1.fa
SEQ2 = sufr/tests/inputs/2.fa
SEQ3 = sufr/tests/inputs/3.fa
SUFR = ./target/release/sufr

#CREATE = cargo run -- create --log info
CREATE = $(SUFR) create --log info
CREATE_DEBUG = cargo run -- create --log debug

#READ   = cargo run -- check
#CHECK  = cargo run -- read
#CHECK  = $(SUFR) read
#READ   = $(SUFR) check

perf:
	perf record --call-graph dwarf $(SUFR) create -t 16 --log info $(CHR1)

s1:
	$(CREATE_DEBUG) $(SEQ1) -n 2 --dna --check

s2:
	$(CREATE_DEBUG) $(SEQ2) --dna --check

s3:
	$(CREATE_DEBUG) $(SEQ3) --dna --check

ecoli:
	$(CREATE) ../data/ecoli.fa --dna --check

chr1:
	$(CREATE) --dna -n 64 $(CHR1)

valcache:
	valgrind --tool=cachegrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info

valcall:
	valgrind --tool=callrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info
