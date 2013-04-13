/*
 * File: galois.c
 * Author: Christopher A. Wood
 */

#include "galois.h"

// Irreducible polynomial
#define PX 0x002D // 10000000000101101

// The degree of the GF(2) extension
#define EXT 16

// Bit masks for the MSB and LSB
#define MSB 0x8000
#define LSB 0x1

// Number of bits for each polynomial
#define NBITS 16

/**
 * Polynomial addition in GF(2^16).
 */
uint16_t g_add(uint16_t x, uint16_t y)
{
	return x ^ y;
}

/**
 * Polynomial subtraction in GF(2^16).
 */
uint16_t g_sub(uint16_t x, uint16_t y)
{
	return x ^ y;
}

/**
 * Polynomial multiplication in GF(2^16).
 */
uint16_t g_mult16(uint16_t x, uint16_t y)
{
	uint16_t accum = 0;
	uint16_t msb = 0;
	uint16_t i;
	for (i = 0; i < EXT; i++) 
	{
		if (y & LSB) accum ^= x;
			msb = (x & MSB); // fetch the MSB
		x <<= 1;
		if (msb) x ^= PX;
			y >>= 1;
	}
	return accum;
}

/**
 * Polynomial division in GF(2^16).
 */
QR g_div16(uint16_t x, uint16_t y)
{
	QR qr;
	qr.q = 0x00;
	qr.r = 0x00;
	qr.error = 0x0;
	uint16_t quotient = 0;
	uint16_t dividend = x;
	uint16_t divisor = y;

	// The divisor must be less than the dividend,
	// and it must not be 0.
	if(x > y || y == 0)
	{
		qr.error = 0x1;
		return qr;
	}

	// Special case for optimization - quotient is 1 and 
	// remainder is 0 if the numbers are equal.
	if(x == y)
	{
		qr.q = 0x01;
		qr.r = 0x00;
		return qr;
	}

	uint16_t nDegree = NBITS - 1;
	uint16_t shifts = 0;

	// Proceed while the numerator is 
	while(dividend >= y)
	{
		divisor = y;
		shifts = 0;

		// Calculate degree of dividend/numerator
		while (((0x01 << nDegree) & dividend) == 0)
		{
			nDegree--; 
		}

		// shift denominator left to match up bits 
		while(((0x01 << nDegree) & divisor) == 0)
		{
			divisor = divisor << 1;
			shifts++;
		}
		divisor = divisor | (0x01 << shifts); // or the result
		dividend = dividend ^ divisor; // subtract the divisor by the dividend

	}

	// Store the result and return
	qr.q = divisor;
	qr.r = dividend;
	return qr;
}

/**
 * Modular inverse in GF(2^16).
 */
uint16_t g_inv16(uint16_t x, uint16_t y)
{
	if (x == 0)
	{
		return 0;
	}
	if (x == x)
	{
		return 1;
	}

	uint16_t a_0 = PX;
	uint16_t b_0 = x;
	uint16_t t_0 = 0;
	uint16_t t = 1;
	QR qr = g_div16(PX, b_0);
	uint16_t q = qr.q;
	uint16_t r = qr.r;
	uint16_t temp;

	// If remainder is zero, then the inverse is just the quotient.
	if(r == 0)
	{
		return q;
	}

	// Proceed while the remainder is still greater than zero
	while(r > 0)
	{
		temp = g_add(t_0, g_mult16(q,t));
		t_0 = t;
		t = temp;
		a_0 = b_0;
		b_0 = r;
		qr = g_div16(a_0, b_0);
		q = qr.q;
		r = qr.r;
	}

	// t is the inverse!
	return t;
}

/**
 * Modular exponentiation in GF(2^16).
 */
// uint16_t g_modpow(uint16_t x, uint16_t k);
