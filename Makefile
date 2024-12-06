CHR1    = ../data/chr1.fa
HUMAN   = ../data/hs1.fa
ELEGANS = ../data/elegans.fa
SEQ1    = data/inputs/1.fa
SEQ2    = data/inputs/2.fa
SEQ3    = data/inputs/3.fa
SEQA    = data/inputs/a.fa
SMOL    = data/inputs/smol.fa
ABBA    = data/inputs/abba.fa
SUFR    = ./target/release/sufr

CREATE_RELEASE = /usr/bin/time -l $(SUFR) --log info create 
CREATE_DEBUG = cargo run -- --log debug create 

perf:
	perf record --call-graph dwarf $(SUFR) create -t 16 --log info $(CHR1)

r1:
	$(SEARCH) C 1.sufr

testfiles: s1 s1s s1n s2 s2d s2n s2s s2ns s3 abba
	cp [1-3]*.sufr data/expected/
	cp abba.sufr data/inputs/

smol:
	$(CREATE_DEBUG) $(SMOL) --dna -m 4

smol-mask:
	$(CREATE_DEBUG) $(SMOL) --dna --seed-mask 1001 -o smol.mask.sufr

masked:
	$(CREATE_DEBUG) $(SEQA) --dna --seed-mask 1001 -o a-m1011.sufr

s1:
	$(CREATE_DEBUG) $(SEQ1) -n 1 --dna -o 1.sufr

s1m:
	$(CREATE_DEBUG) $(SEQ1) -n 1 --dna --seed-mask 1001 -o 1m.sufr

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

chr1-mask:
	$(CREATE_RELEASE) --dna -n 64 -s 111010010100110111 -o chr1-mask.sufr $(CHR1) 

human:
	$(CREATE_RELEASE) -t 8 --dna -n 800 $(HUMAN)

humask:
	$(CREATE_RELEASE) -t 8 --dna -n 800 -s 111010010100110111 $(HUMAN)

valcache:
	valgrind --tool=cachegrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info

valcall:
	valgrind --tool=callrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info
