/*
 * galois.h
 *
 *  Created on: Mar 29, 2013
 *      Author: caw4567
 */

#ifndef GALOIS_H_
#define GALOIS_H_

#include <stdint.h>

// Quotient and remainder struct
typedef struct 
{
	uint16_t q;
	uint16_t r;
	uint8_t error;
} QR;

/**
 * Polynomial addition in GF(2^n).
 */
uint16_t g_add(uint16_t x, uint16_t y);

/**
 * Polynomial subtraction in GF(2^n).
 */
uint16_t g_sub(uint16_t x, uint16_t y);

/**
 * Polynomial multiplication in GF(2^16).
 */
uint16_t g_mult16(uint16_t x, uint16_t y);

/**
 * Polynomial division in GF(2^16).
 */
QR g_div16(uint16_t x, uint16_t y);

/**
 * Modular inverse in GF(2^16).
 */
uint16_t g_inv16(uint16_t x, uint16_t y);

/**
 * Modular exponentiation in GF(2^16).
 */
// uint16_t g_modpow(uint16_t x, uint16_t k);

#endif /* GALOIS_H_ */