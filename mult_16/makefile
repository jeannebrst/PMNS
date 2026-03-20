CC = gcc
CFLAGS = -O3 -funroll-loops -mfpu=neon -mcpu=cortex-a72
LIBS = -lmphell -lgmp

SRC_COMMON = mult_16.c ../util.c

SRC_MAIN = main_16.c $(SRC_COMMON)
SRC_PERFS = main_16_perfs.c $(SRC_COMMON)

OUT_MAIN = mult_16
OUT_PERFS = mult_16_perfs

# Compilation par défaut:
all: $(OUT_MAIN) $(OUT_PERFS)

# Programme normal:
$(OUT_MAIN): $(SRC_MAIN)
	$(CC) $(CFLAGS) $(SRC_MAIN) $(LIBS) -o $(OUT_MAIN)

# Programme benchmark:
$(OUT_PERFS): $(SRC_PERFS)
	$(CC) $(CFLAGS) $(SRC_PERFS) $(LIBS) -o $(OUT_PERFS)

# Run interactif:
run: $(OUT_MAIN)
	./$(OUT_MAIN)

# Run interactif:
perfs: $(OUT_PERFS)
	./$(OUT_PERFS)

# Nettoyage:
clean:
	rm -f $(OUT_MAIN) $(OUT_PERFS)