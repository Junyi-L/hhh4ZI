## Run demo scripts and compare output to previously saved results
R := R
DEMOSCRIPTS := $(wildcard *.R)
ROUTFILES := ${DEMOSCRIPTS:.R=.Rout}

all : ${ROUTFILES}

${ROUTFILES} :
%.Rout : %.R
	LANGUAGE=en ${R} -f "$<" > "$@" 2>&1
	@if [ -f Rplots.pdf ]; then mv Rplots.pdf ${<:.R=.pdf}; fi

diff : $(wildcard *.save)
	$(foreach outfile,$(wildcard $(basename $^)),${R} CMD Rdiff ${outfile} ${outfile}.save;)

clean :
	rm -f ${DEMOSCRIPTS:.R=.pdf} ${DEMOSCRIPTS:.R=.Rout}

.PHONY: all clean diff
