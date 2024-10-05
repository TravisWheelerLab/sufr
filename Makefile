CHR1  = ../data/chr1.fa
HUMAN = ../data/hs1.fa.gz
SEQ1  = sufr/tests/inputs/1.fa
SEQ2  = sufr/tests/inputs/2.fa
SEQ3  = sufr/tests/inputs/3.fa
SUFR  = ./target/release/sufr

#CREATE = cargo run -- create --log info
CREATE_RELEASE = /usr/bin/time -l $(SUFR) --log info create 
CREATE_DEBUG = cargo run -- --log debug create 

SEARCH = cargo run -- search

#READ   = cargo run -- check
#CHECK  = cargo run -- read
#CHECK  = $(SUFR) read
#READ   = $(SUFR) check

perf:
	perf record --call-graph dwarf $(SUFR) create -t 16 --log info $(CHR1)

r1:
	$(SEARCH) C 1.sufr

s1:
	$(CREATE_DEBUG) $(SEQ1) -n 2 --dna 

s2:
	$(CREATE_DEBUG) $(SEQ2) -n 3 --dna 

s3:
	$(CREATE_DEBUG) $(SEQ3) --dna 

ecoli:
	$(CREATE_RELEASE) --dna ../data/ecoli.fa 

chr1:
	$(CREATE_RELEASE) --dna -n 64 $(CHR1)

human:
	$(CREATE_RELEASE) -t 16 --dna -n 8000 $(HUMAN)

valcache:
	valgrind --tool=cachegrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info

valcall:
	valgrind --tool=callrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info
