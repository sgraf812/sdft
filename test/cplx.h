#pragma once

typedef struct {
	double real;
	double imag;
} cplx;

const static cplx cplx_zero = { 0, 0 };

static cplx cplx_mult(cplx* a, cplx* b) 
{
	cplx ret;
	ret.real = a->real * b->real - a->imag * b->imag;
	ret.imag = a->real * b->imag + a->imag * b->real;
	return ret;
}

static cplx cplx_add(cplx* a, cplx* b) 
{
	cplx ret;
	ret.real = a->real + b->real;
	ret.imag = a->imag + b->imag;
	return ret;
}

static cplx cplx_sub(cplx* a, cplx* b) 
{
	cplx ret;
	ret.real = a->real - b->real;
	ret.imag = a->imag - b->imag;
	return ret;
}

static double cplx_abs(cplx* a) 
{
	return sqrt(a->real*a->real + a->imag * a->imag);
}

static int cplx_equal(cplx* a, cplx* b) 
{
	return a->real == b->real && a->imag == b->imag;
}