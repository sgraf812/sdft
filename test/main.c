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

char *compare_sdft_to_dft(my_complex *new_signal, enum sdft_SignalTraits traits, size_t N)
{
    // 1. Allocate the buffers of the complex number type to use.
    //    Here, we use a custom typedef of just two doubles (double[2]).
    my_complex *signal_buffer = malloc(sizeof(my_complex) * N);
    my_complex *spec_buffer = malloc(sizeof(my_complex) * N);
    my_complex *phase_buffer = malloc(sizeof(my_complex) * N);
    memset(signal_buffer, 0, N * sizeof(my_complex));
    memset(spec_buffer, 0, N * sizeof(my_complex));

    // 2. Allocate and initialize the sdft_State with the help of the first two library functions
    struct sdft_State *s = malloc(sdft_size_of_state());
    sdft_init_from_buffers(
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
    for (size_t i = 0; i < N; ++i) {
        sdft_push_next_sample(s, new_signal + i);
    }

    // Steps 4-6 are rather optional.

    // 4. When we want to inspect what the actual signal values are, we have to unshift the
    //    internal representation first.
    my_complex *sig = sdft_unshift_and_get_signal(s);

    // 5. Now we can also access signal_buffer for the stored signal.
    //    In this case, since we pushed all N values of new_signal, the signal_buffer must contain
    //    the same samples as new_signal.
    for (size_t i = 0; i < N; ++i) {
        MU_ASSERT("signal equals new_signal", my_complex_equal(sig + i, new_signal + i));
    }

    my_complex *expected_spec = malloc(sizeof(my_complex) * N);
    dft(new_signal, expected_spec, N);

    // 6. Here, we compare the spectrum computed by the SDFT to that of a classic DFT.
    my_complex *actual_spec = sdft_get_spectrum(s);
    size_t n_bins = traits == SDFT_REAL_AND_IMAG ? N : N / 2;
    for (size_t i = 0; i < n_bins; ++i) {
        my_complex delta = my_complex_sub(actual_spec + i, expected_spec + i);
        MU_ASSERT("spectrum equal to that of the dft", my_complex_abs(&delta) < 0.001);
    }

    double *abs_spec = malloc(sizeof(double) * N);
    for (size_t i = 0; i < N; i++) {
        abs_spec[i] = my_complex_abs(spec_buffer + i);
    }

    // 7. Finally, we have to be sure to free all allocated buffers.
    //    In this case, we allocated all buffers on the stack, the only exception being
    //    the sdft_State variable, which was malloc'ed. (Technically, we could use VLA to allocate it
    //    on the stack, too).
    free(s);
    free(signal_buffer);
    free(spec_buffer);
    free(phase_buffer);
    free(expected_spec);
    free(abs_spec);

    return 0;
}

char *compare_sdft_to_dft_combined(my_complex *new_signal, enum sdft_SignalTraits traits, size_t N)
{
    //  1. Here, because we want to combine 2 sdfts runs, we allocate buffers twice as big.
    my_complex *signal_buffer = malloc(sizeof(my_complex) * 2 * N);
    my_complex *spec_buffer = malloc(sizeof(my_complex) * 2 * N);
    my_complex *phase_buffer = malloc(sizeof(my_complex) * 2 * N);
    memset(signal_buffer, 0, 2 * N * sizeof(my_complex));
    memset(spec_buffer, 0, 2 * N * sizeof(my_complex));

    // 2. Allocate and initialize the sdft_State with the help of the first two library functions
    struct sdft_State *fst = malloc(sdft_size_of_state());
    struct sdft_State *snd = malloc(sdft_size_of_state());
    struct sdft_State *combined = malloc(sdft_size_of_state());

    // Initialize the two sub struct in the same manner. The buffers of the second state are offset by N.
    sdft_init_from_buffers(fst, SDFT_DOUBLE, signal_buffer, spec_buffer, phase_buffer, N, traits);
    sdft_init_from_buffers(snd, SDFT_DOUBLE, signal_buffer + N, spec_buffer + N, phase_buffer + N, N, traits);

    // Now we combine the fst and the snd state.
    sdft_init_combine(combined, fst, snd);

    // 3. Push through all new N values of our incoming signal.
    //    To make this more interesting (and to actually test that the combination works),
    //    We will push through the values 10 times.
    my_complex *expected_spec = malloc(sizeof(my_complex) * N);
    dft(new_signal, expected_spec, N);
    for (size_t i = 0; i < 10 * N; ++i) {
		if (i >= N) {
			// We can now assert that we indeed always have a valid signal.
			my_complex *actual_sig = sdft_unshift_and_get_signal(combined);
			for (size_t j = 0; j < N; ++j) {
				int equal = my_complex_equal(actual_sig + j, new_signal + (i + j) % N);
				MU_ASSERT("actual_signal doesn't match new_signal", equal);
			}

			if (i % N == 0) {
				my_complex *actual_spec = sdft_get_spectrum(combined);
				size_t n_bins = traits == SDFT_REAL_AND_IMAG ? N : N / 2;
				for (size_t j = 0; j < n_bins; ++j) {
					my_complex delta = my_complex_sub(actual_spec + j, expected_spec + j);
					int equal = my_complex_abs(&delta) < 0.001;
					MU_ASSERT("spectrum isn't equal to that of the dft", equal);
				}
			} 
		}

        sdft_push_next_sample(combined, new_signal + i%N);
    }

    double *abs_spec = malloc(sizeof(double) * N);
    for (size_t i = 0; i < N; i++) {
        abs_spec[i] = my_complex_abs(expected_spec + i);
    }

    // Finally, we have to be sure to free all allocated buffers.
    free(fst);
    free(snd);
    free(combined);
    free(signal_buffer);
    free(spec_buffer);
    free(phase_buffer);
    free(expected_spec);
    free(abs_spec);

    return 0;
}

char *test_build_dft_from_zeros()
{
    my_complex new_signal[] = {
            {51, 0}, {2, 0}, {42, 5}, {0.2, 0.5},
            {1, 0}, {765, 0}, {34, 0}, {2903, 0},
            {4096, 256}, {0, 5334}, {3, 0}, {6, 0},
            {4, 0}, {1, 0}, {0, 74}, {79, 74.5}
    };

    char *msg;
    if ((msg = compare_sdft_to_dft(new_signal, SDFT_REAL_AND_IMAG, 16))) {
        return msg;
    }
    return compare_sdft_to_dft_combined(new_signal, SDFT_REAL_AND_IMAG, 16);
}

char *test_real_signal()
{
    my_complex new_signal[] = {
            {51, 0}, {2, 0}, {42, 0}, {0.2, 0},
            {1, 0}, {765, 0}, {34, 0}, {2903, 0},
            {4096, 0}, {5334, 0}, {3, 0}, {6, 0},
            {4, 0}, {1, 0}, {74, 0}, {79, 0}
    };

    char *msg;
    if ((msg = compare_sdft_to_dft(new_signal, SDFT_REAL_ONLY, 16))) {
        return msg;
    }
    return compare_sdft_to_dft_combined(new_signal, SDFT_REAL_ONLY, 16);
}

char *test_imag_signal()
{
    my_complex new_signal[] = {
            {0, 51}, {0, 2}, {0, 42}, {0, 0.2},
            {0, 1}, {0, 765}, {0, 34}, {0, 2903},
            {0, 4096}, {0, 5334}, {0, 3}, {0, 6},
            {0, 4}, {0, 1}, {0, 74}, {0, 79}
    };

    char *msg;
    if ((msg = compare_sdft_to_dft(new_signal, SDFT_IMAG_ONLY, 16))) {
        return msg;
    }
    return compare_sdft_to_dft_combined(new_signal, SDFT_IMAG_ONLY, 16);
}

#include "actual_signal.h"

char *test_actual_signal()
{
    const size_t N = 512;
    my_complex new_signal[512];
    for (size_t i = 0; i < N; ++i) {
        new_signal[i].real = actual_signal[i];
        new_signal[i].imag = 0;
    }

    char *msg;
    if ((msg = compare_sdft_to_dft(new_signal, SDFT_REAL_ONLY, 16))) {
        return msg;
    }
    return compare_sdft_to_dft_combined(new_signal, SDFT_REAL_ONLY, 16);
}

char *test_suite(void)
{
    MU_RUN_TEST(test_build_dft_from_zeros);
    MU_RUN_TEST(test_real_signal);
    MU_RUN_TEST(test_imag_signal);
    MU_RUN_TEST(test_actual_signal);
    return 0;
}

int main(int ac, char **av)
{
    char *result = test_suite();

    printf("number of tests run: %d\n", tests_run);
    if (result) printf("FAIL: %s\n", result);
    return 0 != result;
}