
CC = gcc -O3 -g -Wall -Wno-unused-function -Wno-stringop-truncation

CC_CLANG = clang-12 -O3 -g -Wall -Wno-unused-function -Wnewline-eof

CC_DBG = gcc -O0 -g -pg -Wall -Wno-unused-function -Wno-unused-variable -Wno-stringop-truncation
LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c, %, $(srcs))
bins_clang = $(patsubst %.c, %_clang, $(srcs))
bins_dbg = $(patsubst %.c, %_d, $(srcs))
deps = Makefile
includes = -I../../../../common/c

bins: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $(includes) $< -o $@ $(LM)

$(bins_clang): %_clang : %.c $(deps)
	$(CC_CLANG) $(includes) $< -o $@ $(LM)

$(bins_dbg): %_d : %.c $(deps)
	$(CC_DBG) $(includes) $< -o $@ $(LM)

clean::
	$(RM) -rf $(bins) $(bins_clang) $(bins_dbg) .*.swp *~

