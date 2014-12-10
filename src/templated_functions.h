#ifndef APPEND_TYPE_SUFFIX
#error APPEND_TYPE_SUFFIX has to be defined before this header can be included.
#endif

#include <math.h>
#include <assert.h>

#define CPLX APPEND_TYPE_SUFFIX(cplx)

static void APPEND_TYPE_SUFFIX(impl_push_next_sample)(struct sdft_State *s, void *next_sample)
{
    assert(s->signal_index >= 0 && s->signal_index < s->number_of_samples);

    CPLX ns = *(CPLX *) next_sample;
    assert(s->signal_traits == SDFT_REAL_AND_IMAG
            || s->signal_traits == SDFT_REAL_ONLY && APPEND_TYPE_SUFFIX(cimag)(ns) == 0
            || s->signal_traits == SDFT_IMAG_ONLY && APPEND_TYPE_SUFFIX(creal)(ns) == 0);

    CPLX *signal = s->signal;
    CPLX delta = ns - signal[s->signal_index];

    size_t n_samples = s->number_of_samples;
    size_t n_bins = number_of_bins(s);
    CPLX *po = s->phase_offsets;
    CPLX *spec = s->spectrum;
    for (size_t i = 0; i < n_bins; ++i) {
        spec[i] = (spec[i] + delta) * po[i];
    }

    signal[s->signal_index] += delta;
    if (++s->signal_index == n_samples)
        s->signal_index = 0;
}

void APPEND_TYPE_SUFFIX(impl_unshift_signal)(struct sdft_State *s)
{
    size_t n = s->number_of_samples;

    if (s->signal_index == 0) {
        // we are lucky, there's nothing to shift
        return;
    }

    // instead of taking the signal index as offset, we take n - signal_index
    // as offset to follow a cycle through the array.
    assert(n > s->signal_index);
    size_t ofs = s->signal_index;

    CPLX *signal = s->signal;
    CPLX first = signal[0];
    size_t cur = 0;
    for (int i = 0; i < n - 1; ++i) {
        size_t next = (cur + ofs) % n;
        signal[cur] = signal[next];
        cur = next;
    }

    assert((cur + ofs) % n == 0); // or: next == 0
    signal[cur] = first;
}

static void APPEND_TYPE_SUFFIX(impl_init)(
        struct sdft_State *s,
        void *signal,
        void *spectrum,
        void *buffer,
        size_t number_of_samples,
        enum sdft_SignalTraits signal_traits)
{
    s->signal = signal;
    s->spectrum = spectrum;
    s->phase_offsets = buffer;
    s->signal_index = 0;
    s->number_of_samples = number_of_samples;
    s->signal_traits = signal_traits;
    s->push_next_sample = &APPEND_TYPE_SUFFIX(impl_push_next_sample);
    s->unshift_signal = &APPEND_TYPE_SUFFIX(impl_unshift_signal);

    // check for violations of the signal traits
    for (int i = 0; i < number_of_samples; ++i) {
        assert(s->signal_traits == SDFT_REAL_AND_IMAG
                || s->signal_traits == SDFT_REAL_ONLY && APPEND_TYPE_SUFFIX(cimag)(s->signal[i]) == 0
                || s->signal_traits == SDFT_IMAG_ONLY && APPEND_TYPE_SUFFIX(creal)(s->signal[i]) == 0);
    }

    // generate the phase offsets
    CPLX *offsets = buffer;
    for (size_t i = 0; i < number_of_samples; i++) {
        APPEND_TYPE_SUFFIX(cplx) angle = 2 * M_PI * i / number_of_samples * I;
        offsets[i] = APPEND_TYPE_SUFFIX(cexp)(angle);
    }
}

#undef CPLX