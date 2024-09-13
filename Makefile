CHR1 = ../data/t2t-chr1.txt
#CHR1 = /xdisk/twheeler/data/genomes/human/t2t/chr1.txt

CHECK  = cargo run -- read
CREATE = cargo run -- create --log info
READ   = cargo run -- check

#SUFR   = ./target/release/sufr
#CHECK  = $(SUFR) read
#CREATE = $(SUFR) create --log info
#READ   = $(SUFR) check

perf:
	perf record --call-graph dwarf $(SUFR) create -t 16 --log info $(CHR1)

create-s1:
	#$(CREATE) tests/inputs/seq1.txt -o seq1.sa
	cargo run -- create --log debug tests/inputs/seq1.txt -o seq1.sa -n 3

read-s1:
	$(READ) -s tests/inputs/seq1.txt -a seq1.sa -e 1-100

check-s1:
	$(READ) -s tests/inputs/seq1.txt -a seq1.sa

create-s2:
	$(CREATE) tests/inputs/seq2.txt -o seq2.sa

check-s2:
	$(READ) -s tests/inputs/seq2.txt -a seq2.sa

read-s2:
	$(READ) -s tests/inputs/seq2.txt -a seq2.sa -e 1-100

create-ecoli:
	$(CREATE) ../data/ecoli.txt -o ecoli.sa

check-ecoli:
	$(READ) -s ../data/ecoli.txt -a ecoli.sa

create-t2t-chr1:
	$(CREATE) ../data/t2t-chr1.txt -o t2t-chr1.sa -n 64

check-t2t-chr1:
	$(READ) -s ../data/t2t-chr1.txt -a t2t-chr1.sa

create-chr1:
	$(CREATE) ../data/chr1.txt -o t2t-chr1.sa --ignore-start-n

check-chr1:
	$(READ) -s ../data/chr1.txt -a t2t-chr1.sa

s1: create-s1 check-s1

s2: create-s2 check-s2

ecoli: create-ecoli check-ecoli

t2t-chr1: create-t2t-chr1 check-t2t-chr1

chr1: create-chr1 check-chr1

valcache:
	valgrind --tool=cachegrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info

valcall:
	valgrind --tool=callrind ./target/release/sufr create ../data/chr1.fa --ignore-start-n -o chr1.sa --log info
