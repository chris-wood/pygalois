#include "galois.h"
#include <stdio.h>

int main()
{
	uint16_t x = 0x142a;
	uint16_t y = 0x521f;
	QR result;
	printf("%x + %x = %x\n", x, y, g_add(x, y));
	printf("%x - %x = %x\n", x, y, g_sub(x, y));
	printf("%x * %x = %x\n", x, y, g_mult16(x, y));
	result = g_div16(x, y);
	printf("%x / %x = %x, %x\n", x, y, result.q, result.r);

	return 0;
}
