#pragma once

typedef struct {
    double real;
    double imag;
} my_complex;

const static my_complex my_complex_zero = {0, 0};

static my_complex my_complex_mult(my_complex *a, my_complex *b)
{
    my_complex ret;
    ret.real = a->real * b->real - a->imag * b->imag;
    ret.imag = a->real * b->imag + a->imag * b->real;
    return ret;
}

static my_complex my_complex_add(my_complex *a, my_complex *b)
{
    my_complex ret;
    ret.real = a->real + b->real;
    ret.imag = a->imag + b->imag;
    return ret;
}

static my_complex my_complex_sub(my_complex *a, my_complex *b)
{
    my_complex ret;
    ret.real = a->real - b->real;
    ret.imag = a->imag - b->imag;
    return ret;
}

static double my_complex_abs(my_complex *a)
{
    return sqrt(a->real * a->real + a->imag * a->imag);
}

static int my_complex_equal(my_complex *a, my_complex *b)
{
    return a->real == b->real && a->imag == b->imag;
}