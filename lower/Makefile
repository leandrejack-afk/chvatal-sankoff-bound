CC ?= clang
CFLAGS ?= -O3 -march=native
all: klow klow_plain kedit kedit_plain check/lowcheck
klow: klow.c ; $(CC) $(CFLAGS) -o $@ $< -lpthread
klow_plain: klow_plain.c ; $(CC) $(CFLAGS) -o $@ $< -lpthread
kedit: kedit.c ; $(CC) $(CFLAGS) -o $@ $< -lpthread
kedit_plain: kedit_plain.c ; $(CC) $(CFLAGS) -o $@ $<
check/lowcheck: check/lowcheck.c ; $(CC) $(CFLAGS) -o $@ $< -lpthread
clean: ; rm -f klow klow_plain kedit kedit_plain check/lowcheck
.PHONY: all clean
