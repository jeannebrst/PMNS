SIZES = 14 15 16

# Compilation par défaut:
all:
	@echo "Usage:"
	@echo "make <size>			# Ex: make 14"
	@echo "make perfs_<size>	# Ex: make perfs_14"

# Lancer main_X:
$(SIZES):
	@cd mult_$@ && $(MAKE) run

# Lancer main_X_perfs:
perfs: ; @echo "Please run: make perfs_<size>"

perfs_%:
	@size=$*; \
	if [ "$$size" != 14 ] && [ "$$size" != 15 ] && [ "$$size" != 16 ]; then \
		echo "Error: size must be 14, 15 or 16"; exit1; \
	fi; \
	cd mult_$${size} && $(MAKE) perfs

clean:
	@cd mult_14/ && make clean && cd ..
	@cd mult_15/ && make clean && cd ..
	@cd mult_16/ && make clean && cd ..