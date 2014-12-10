#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sdft/sdft.h"
#include "minunit.h"

int tests_run = 0;

#define cplx double complex

void dft(cplx *signal, cplx *spec, size_t N)
{
    for (int i = 0; i < N; ++i) {
        spec[i] = 0;
        for (int j = 0; j < N; ++j) {
            cplx a = -2.0 * M_PI * i * j / N * I;
            spec[i] += signal[j] * cexp(a);
        }
    }
}

static const size_t N = 16;

char *compare_sdft_to_dft(cplx (*new_signal)[N], enum sdft_SignalTraits traits)
{
    // 1. Allocate the buffers of the complex number type to use.
    cplx signal_buffer[N];
    cplx spec_buffer[N];
    cplx phase_buffer[N];
    memset(signal_buffer, 0, N * sizeof(cplx));
    memset(spec_buffer, 0, N * sizeof(cplx));

    // 2. Allocate and initialize the sdft_State with the help of the first two library functions
    struct sdft_State *s = malloc(sdft_size_of_state());
    sdft_init(
            s,             // pointer to the allocated state struct
            SDFT_DOUBLE,   // floating point type to use for the complex computations
            signal_buffer, // user allocated and initialized buffer containing the signal, length at least N
            spec_buffer,   // user allocated and initialized buffer containing the spectrum matching the signal
                           // its length depends on the supplied traits value, at least N/2 or at least N.
                           // Here, we always allocate N to be safe.
            phase_buffer,  // user allocated buffer used and initialized internal by the algorithm.
            N,             // length of the signal / number of samples
            traits);       // signal traits of the intended usage. See the docstrings.

    // 3. Push through all new N values of our incoming signal.
    //    This gets the actual work done and computes a new spectrum for each sample.
    //    Throughout the rest of the function, we can access the spectrum by spec_buffer.
    for (int i = 0; i < N; ++i) {
        sdft_push_next_sample(s, *new_signal + i);
    }

    // Steps 4- are rather optional.

    // 4. When we want to inspect what the actual signal values are, we have to unshift the
    //    internal representation first.
    sdft_unshift_signal(s);

    // 5. Now we can also access signal_buffer for the stored signal.
    //    In this case, since we pushed all N values of new_signal, the signal_buffer must contain
    //    the same samples as new_signal.
    for (int i = 0; i < N; ++i) {
        MU_ASSERT("signal equals new_signal", signal_buffer[i] == (*new_signal)[i]);
    }

    cplx expected_spec[N];
    dft(*new_signal, expected_spec, N);

    // 6. Here, we compare the spectrum computed by the SDFT to that of a classic DFT.
    size_t n_bins = traits == SDFT_REAL_AND_IMAG ? N : N / 2;
    for (int i = 0; i < n_bins; ++i) {
        MU_ASSERT("spectrum equal to that of the dft", cabs(spec_buffer[i] - expected_spec[i]) < 0.001);
    }

    // 7. Finally, we have to be sure to free all allocated buffers.
    //    In this case, we allocated all buffers on the stack, the only exception being
    //    the sdft_State variable, which was malloc'ed. (Technically, we could use VLA to allocate it
    //    on the stack, too).
    free(s);

    return 0;
}

char *test_build_dft_from_zeros()
{
    cplx new_signal[N] = {
            51, 2, 42 + 5 * I, 0.2 + 0.5 * I,
            1, 765, 34, 2903,
            4096 + 256 * I, 5334 * I, 3, 6,
            4, 1, 74 * I, 79 + 74.5 * I
    };

    return compare_sdft_to_dft(&new_signal, SDFT_REAL_AND_IMAG);
}

char *test_real_signal()
{
    cplx new_signal[N] = {
            51, 2, 42, 0.2,
            1, 765, 34, 2903,
            4096, 5334, 3, 6,
            4, 1, 74, 79
    };

    return compare_sdft_to_dft(&new_signal, SDFT_REAL_ONLY);
}

char *test_imag_signal()
{
    cplx new_signal[N] = {
            51 * I, 2 * I, 42 * I, 0.2 * I,
            1 * I, 765 * I, 34 * I, 2903 * I,
            4096 * I, 5334 * I, 3 * I, 6 * I,
            4 * I, 1 * I, 74 * I, 79 * I
    };

    return compare_sdft_to_dft(&new_signal, SDFT_IMAG_ONLY);
}

char *test_suite(void)
{
    MU_RUN_TEST(test_build_dft_from_zeros);
    MU_RUN_TEST(test_real_signal);
    MU_RUN_TEST(test_imag_signal);
    return 0;
}

int main(int ac, char **av)
{
    char *result = test_suite();

    printf("number of tests run: %d\n", tests_run);
    if (result) printf("FAIL: %s\n", result);
    return 0 != result;
}