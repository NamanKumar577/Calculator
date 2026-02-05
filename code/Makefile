MCU = atmega328p
F_CPU = 16000000UL
FORMAT = ihex
TARGET = main
SRC = main.c
OPT = 2

CC = avr-gcc
OBJCOPY = avr-objcopy
OBJDUMP = avr-objdump
NM = avr-nm
SIZE = avr-size
REMOVE = rm -f

CSTANDARD = -std=gnu99
CWARN = -Wall -Wstrict-prototypes
CFLAGS = -O$(OPT) $(CWARN) $(CSTANDARD)
CPPFLAGS = -DF_CPU=$(F_CPU)

ALL_CFLAGS = -mmcu=$(MCU) -I. $(CFLAGS) $(CPPFLAGS)
LDFLAGS = -mmcu=$(MCU)

OBJ = $(SRC:.c=.o)

all: $(TARGET).hex $(TARGET).eep

$(TARGET).elf: $(OBJ)
	$(CC) $(ALL_CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) -c $(ALL_CFLAGS) $< -o $@

$(TARGET).hex: $(TARGET).elf
	$(OBJCOPY) -O $(FORMAT) -R .eeprom $< $@

$(TARGET).eep: $(TARGET).elf
	-$(OBJCOPY) -j .eeprom --set-section-flags=.eeprom="alloc,load" \
	--change-section-lma .eeprom=0 -O $(FORMAT) $< $@

clean:
	$(REMOVE) $(TARGET).hex $(TARGET).eep $(TARGET).elf $(OBJ)

# Optional upload target (only run when board is connected)
AVRDUDE_PROGRAMMER = arduino
AVRDUDE_PORT = /dev/ttyACM0
AVRDUDE_BUDRATE = 115200
AVRDUDE = avrdude

flash: $(TARGET).hex
	$(AVRDUDE) -p $(MCU) -P $(AVRDUDE_PORT) -c $(AVRDUDE_PROGRAMMER) -b $(AVRDUDE_BUDRATE) \
	-U flash:w:$(TARGET).hex

.PHONY: all clean flash
