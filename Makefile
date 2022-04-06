# Makefile to compile SOC_data_sim.c file
# Contributors: Andrey Yakymovych

CC = gcc
CFLAGS  = -O1 -Wall -Werror -lm -Wno-unused-result

TARGET = SOC_data_sim

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -g -o $(TARGET) $(TARGET).c

clean:
	$(RM) $(TARGET)

