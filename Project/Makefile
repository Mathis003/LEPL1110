CC = gcc
CFLAGS = -Wall -Wextra -std=c99
SRC = main.c
TARGET = main

.PHONY: all clean

all: $(TARGET)

# Build the program
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $^ -o $@

# Not use it if you want to run the program with options => ./main [options]
run: all
	./$(TARGET)

# Clean the program
clean:
	rm -f $(TARGET)