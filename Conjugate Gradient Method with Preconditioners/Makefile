        CC = gcc

    CFLAGS = -O0
    LFLAGS = -lm

      PROG = cgSolver
      MODULES = sislin utils $(PROG)
      OBJS = $(addsuffix .o,$(MODULES))
      SRCS = $(addsuffix .c,$(MODULES)) $(addsuffix .h,$(MODULES))

# Lista de arquivos para distribuição
DISTFILES = *.c *.h Makefile LEIAME
DISTDIR = os24-k24

.PHONY: clean purge dist all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

debug:   CFLAGS+=-D__DEBUG__
debug: $(PROG)

clean:
	@echo "Limpando ....."
	@rm -rf core *~ *.bak
	@rm -f a.out *.o $(PROG)


dist: clean
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@ln -s . $(DISTDIR)
	@tar -chzvf $(DISTDIR).tgz $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)
