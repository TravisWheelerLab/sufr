CHR1  = ../data/chr1.fa
HUMAN = ../data/hs1.fa
ELEGANS = ../data/elegans.fa
SEQ1  = sufr/tests/inputs/1.fa
SEQ2  = sufr/tests/inputs/2.fa
SEQ3  = sufr/tests/inputs/3.fa
ABBA  = libsufr/tests/inputs/abba.fa
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

testfiles: s1 s1s s1n s2 s2d s2n s2s s2ns s3 abba
	cp [1-3]*.sufr sufr/tests/expected/
	cp abba.sufr libsufr/tests/inputs/

s1:
	$(CREATE_DEBUG) $(SEQ1) -n 2 --dna -o 1.sufr

s1s:
	$(CREATE_DEBUG) $(SEQ1) -n 2 --dna -o 1s.sufr --ignore-softmask

s1n:
	$(CREATE_DEBUG) $(SEQ1) -n 2 --dna --allow-ambiguity -o 1n.sufr

s2:
	$(CREATE_DEBUG) $(SEQ2) -n 3 --dna 
	
s2d:
	$(CREATE_DEBUG) $(SEQ2) -n 2 --dna --sequence-delimiter 'N' -o 2d.sufr

s2n:
	$(CREATE_DEBUG) $(SEQ2) -n 3 --dna --allow-ambiguity -o 2n.sufr

s2s:
	$(CREATE_DEBUG) $(SEQ2) -n 3 --dna --ignore-softmask -o 2s.sufr

s2ns:
	$(CREATE_DEBUG) $(SEQ2) -n 3 --dna --allow-ambiguity --ignore-softmask -o 2ns.sufr

s3:
	$(CREATE_DEBUG) $(SEQ3) --dna

abba:
	$(CREATE_DEBUG) $(ABBA)

elegans:
	$(CREATE_RELEASE) --dna -m 12 $(ELEGANS)

ecoli:
	$(CREATE_RELEASE) --dna -m 16 ../data/ecoli.fa 

chr1:
	$(CREATE_RELEASE) --dna -n 64 $(CHR1)

human:
	$(CREATE_RELEASE) -t 8 --dna -n 800 $(HUMAN)

valcache:
	valgrind --tool=cachegrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info

valcall:
	valgrind --tool=callrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info
