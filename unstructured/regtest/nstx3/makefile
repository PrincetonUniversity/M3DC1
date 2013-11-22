green='\e[0;32m'
NC='\e[0m'

all : nstx3

nstx3 : C1ke
	@$(PY) ../compare.py; test "$$?" -eq 0
	@echo -e ${green}Success!${NC}

reset : C1ke
	cp C1ke C1ke.base

C1ke : 
	$(M3DC1)-opt-25/m3dc1 >& C1stdout

clean :
	rm -fr *.h5 C1stdout C1ke fort* *~ profilesdb-* testout