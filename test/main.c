#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sdft/sdft.h"
#include "minunit.h"
#include "my_complex.h"

int tests_run = 0;

void dft(my_complex *signal, my_complex *spec, size_t N)
{
    for (size_t k = 0; k < N; ++k) {
        spec[k] = my_complex_zero;
        for (size_t j = 0; j < N; ++j) {
            const double double_pi = 2 * 3.141592653589793238462643383279502884;
            size_t i = (k * j) % N;
            double angle = -double_pi * i / N;
            // spec[i] += signal[j] * cexp(angle*I);
            my_complex tmp = {cos(angle), sin(angle)};
            tmp = my_complex_mult(signal + j, &tmp);
            spec[k] = my_complex_add(spec + k, &tmp);
        }
    }
}

char *compare_sdft_to_dft(struct sdft_State *s, my_complex *signal, size_t signal_length,
        enum sdft_SignalTraits traits, size_t window_size)
{
    // 3. Push through all values of our incoming signal.
    //    This gets the actual work done and computes a new spectrum for each sample.
    for (size_t i = 0; i < signal_length; ++i) {
        sdft_push_next_sample(s, signal + i);
    }

    // Steps 4-6 are rather optional.

    // 4. When we want to inspect what the samples in the window are, we also have to unshift the
    //    internal representation.
    my_complex *window = sdft_unshift_and_get_window(s);

    // 5. Now we can also access window for the stored samples.
    //    In this case, since we pushed all signal_length values of signal, the window_buffer must contain
    //    the last window_size samples.
    size_t signal_offset = signal_length - window_size;
    for (size_t i = 0; i < window_size; ++i) {
        MU_ASSERT("window values don't equal signal", my_complex_equal(window + i, signal + signal_offset + i));
    }

    // 6. Here, we compare the spectrum computed by the SDFT to that of a classic DFT.
    my_complex *expected_spec = malloc(sizeof(my_complex) * window_size);
    dft(signal + signal_offset, expected_spec, window_size);
    my_complex *actual_spec = sdft_get_spectrum(s);
    size_t n_bins = traits == SDFT_REAL_AND_IMAG ? window_size : window_size / 2;
    for (size_t i = 0; i < n_bins; ++i) {
        my_complex delta = my_complex_sub(actual_spec + i, expected_spec + i);
        MU_ASSERT("spectrum isn't equal to that of the dft", my_complex_abs(&delta) < 0.001);
    }

    // We allocated, so we have to free
    free(expected_spec);

    return 0;
}

char *simple_sdft(my_complex *signal, size_t signal_length, enum sdft_SignalTraits traits, size_t window_size)
{
    // 1. Allocate the buffers of the complex number type to use.
    //    Here, we use a custom typedef of just two doubles (double[2]).
    my_complex *window_buffer = malloc(sizeof(my_complex) * window_size);
    my_complex *spec_buffer = malloc(sizeof(my_complex) * window_size);
    my_complex *phase_buffer = malloc(sizeof(my_complex) * window_size);
    memset(window_buffer, 0, window_size * sizeof(my_complex));
    memset(spec_buffer, 0, window_size * sizeof(my_complex));

    // 2. Allocate and initialize the sdft_State with the help of the first two library functions
    struct sdft_State *s = malloc(sdft_size_of_state());
    sdft_init_from_buffers(
            s,             // pointer to the allocated state struct
            SDFT_DOUBLE,   // floating point type to use for the complex computations
            window_buffer, // user allocated and initialized buffer containing the signal, length at least N
            spec_buffer,   // user allocated and initialized buffer containing the spectrum matching the signal
            // its length depends on the supplied traits value, at least N/2 or at least N.
            // Here, we always allocate N to be safe.
            phase_buffer,  // user allocated buffer used and initialized internal by the algorithm.
            window_size,             // length of the signal / number of samples
            traits);       // signal traits of the intended usage. See the docstrings.

    // 3.-6. See compare_sdft_dft.
    char *msg = compare_sdft_to_dft(s, signal, signal_length, traits, window_size);

    // 7. Finally, we have to be sure to free all allocated buffers.
    //    In this case, we allocated all buffers on the stack, the only exception being
    //    the sdft_State variable, which was malloc'ed. (Technically, we could use a C99 VLA to allocate it
    //    on the stack, too).
    free(s);
    free(window_buffer);
    free(spec_buffer);
    free(phase_buffer);

    return msg;
}

char *combined_sdft(my_complex *signal, size_t signal_length, enum sdft_SignalTraits traits, size_t window_size)
{
    //  1. Here, because we want to combine 2 sdfts runs, we allocate buffers twice as big.
    my_complex *window_buffer = malloc(sizeof(my_complex) * 2 * window_size);
    my_complex *spec_buffer = malloc(sizeof(my_complex) * 2 * window_size);
    my_complex *phase_buffer = malloc(sizeof(my_complex) * 2 * window_size);
    memset(window_buffer, 0, 2 * window_size * sizeof(my_complex));
    memset(spec_buffer, 0, 2 * window_size * sizeof(my_complex));

    // 2. Allocate and initialize the sdft_State with the help of the first two library functions
    struct sdft_State *fst = malloc(sdft_size_of_state());
    struct sdft_State *snd = malloc(sdft_size_of_state());
    struct sdft_State *combined = malloc(sdft_size_of_state());

    // Initialize the two sub struct in the same manner. The buffers of the second state are offset by window_size.
    sdft_init_from_buffers(fst, SDFT_DOUBLE, window_buffer, spec_buffer, phase_buffer, window_size, traits);
    sdft_init_from_buffers(snd, SDFT_DOUBLE, window_buffer + window_size, spec_buffer + window_size,
            phase_buffer + window_size, window_size, traits);

    // Now we combine the fst and the snd state.
    sdft_init_combine(combined, fst, snd);

    char *msg = compare_sdft_to_dft(combined, signal, signal_length, traits, window_size);

    free(fst);
    free(snd);
    free(combined);
    free(window_buffer);
    free(spec_buffer);
    free(phase_buffer);

    return msg;
}

char *run_all_combinations(my_complex *signal, size_t signal_length, enum sdft_SignalTraits traits)
{
    // runs all combinations of test functions (simple, combined) and window_sizes (1..signal_length)

    typedef char *(*dft_function)(my_complex *, size_t, enum sdft_SignalTraits, size_t);
    dft_function test_functions[2] = {&simple_sdft, &combined_sdft};

    char *msg;
    for (size_t window_size = 1; window_size < signal_length; ++window_size) {
        for (size_t i = 0; i < 2; ++i) {
            if ((msg = test_functions[i](signal, signal_length, traits, window_size))) {
                return msg;
            }
            tests_run++;
        }
    }
    return 0;
}

char *test_mixed_signal()
{
    my_complex signal[] = {
            {51, 0}, {2, 0}, {42, 5}, {0.2, 0.5},
            {1, 0}, {765, 0}, {34, 0}, {2903, 0},
            {4096, 256}, {0, 5334}, {3, 0}, {6, 0},
            {4, 0}, {1, 0}, {0, 74}, {79, 74.5}
    };

    return run_all_combinations(signal, 16, SDFT_REAL_AND_IMAG);
}

char *test_real_signal()
{
    my_complex signal[] = {
            {51, 0}, {2, 0}, {42, 0}, {0.2, 0},
            {1, 0}, {765, 0}, {34, 0}, {2903, 0},
            {4096, 0}, {5334, 0}, {3, 0}, {6, 0},
            {4, 0}, {1, 0}, {74, 0}, {79, 0}
    };

    return run_all_combinations(signal, 16, SDFT_REAL_ONLY);
}

char *test_imag_signal()
{
    my_complex signal[] = {
            {0, 51}, {0, 2}, {0, 42}, {0, 0.2},
            {0, 1}, {0, 765}, {0, 34}, {0, 2903},
            {0, 4096}, {0, 5334}, {0, 3}, {0, 6},
            {0, 4}, {0, 1}, {0, 74}, {0, 79}
    };

    return run_all_combinations(signal, 16, SDFT_IMAG_ONLY);
}

#include "actual_signal.h"

char *test_actual_signal()
{
    const size_t N = 512;
    my_complex signal[512];
    for (size_t i = 0; i < N; ++i) {
        signal[i].real = actual_signal[i];
        signal[i].imag = 0;
    }

    return run_all_combinations(signal, N, SDFT_REAL_AND_IMAG);
}

char *test_suite(void)
{
    MU_RUN_TESTS(test_mixed_signal);
    MU_RUN_TESTS(test_real_signal);
    MU_RUN_TESTS(test_imag_signal);
    MU_RUN_TESTS(test_actual_signal);
    return 0;
}

int main(int ac, char **av)
{
    char *result = test_suite();

    printf("number of tests run: %d\n", tests_run);
    if (result) printf("FAIL: %s\n", result);
    return 0 != result;
}