#include <cassert>

#include <complex>

#include "sdft/sdft.h"

//
// Forward declarations for precision dependent implementations
//

template<typename Float>
static void impl_init(
        struct sdft_State *s,
        void *signal,
        void *spectrum,
        void *buffer,
        size_t number_of_samples,
        enum sdft_SignalTraits signal_traits);

template<typename Float>
void impl_push_next_sample(struct sdft_State *s, void *next_sample);

template<typename Float>
void impl_unshift_signal(struct sdft_State *s);

//
// Exported state struct an size_of_state, init, push_next_sample and unshift_signal functions
//

struct sdft_State {
    void *signal;
    void *spectrum;
    void *phase_offsets;
    size_t signal_index;
    size_t number_of_samples;
    void (*push_next_sample)(struct sdft_State*, void *next_sample);
    void (*unshift_signal)(struct sdft_State*);
    enum sdft_SignalTraits signal_traits;
};

size_t sdft_size_of_state()
{
    return sizeof(struct sdft_State);
}

void sdft_init(
        struct sdft_State *s,
        enum sdft_FloatPrecision precision,
        void *signal,
        void *spectrum,
        void *buffer,
        size_t number_of_samples,
        enum sdft_SignalTraits signal_traits)
{
    switch (precision) {
        case SDFT_SINGLE:
            impl_init<float>(s, signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
        case SDFT_DOUBLE:
            impl_init<double>(s, signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
        case SDFT_LONG_DOUBLE:
            impl_init<long double>(s, signal, spectrum, buffer, number_of_samples, signal_traits);
            break;
    }
}

void sdft_push_next_sample(struct sdft_State *s, void *next_sample)
{
    s->push_next_sample(s, next_sample);
}

void sdft_unshift_signal(struct sdft_State *s)
{
    s->unshift_signal(s);
}

//
// Templated implementations of the 3 precision dependent functions
//

template<typename Float>
static void impl_init(
        struct sdft_State *s,
        void *signal,
        void *spectrum,
        void *buffer,
        size_t number_of_samples,
        enum sdft_SignalTraits signal_traits)
{
    typedef std::complex<Float> cplx;
    s->signal = signal;
    s->spectrum = spectrum;
    s->phase_offsets = buffer;
    s->signal_index = 0;
    s->number_of_samples = number_of_samples;
    s->signal_traits = signal_traits;
    s->push_next_sample = &impl_push_next_sample<Float>;
    s->unshift_signal = &impl_unshift_signal<Float>;

    cplx *sig = (cplx *) signal;
    cplx *offsets = (cplx *) buffer;

    // check for violations of the signal traits
    for (size_t i = 0; i < number_of_samples; ++i) {
        assert(s->signal_traits == SDFT_REAL_AND_IMAG
                || s->signal_traits == SDFT_REAL_ONLY && imag(sig[i]) == 0
                || s->signal_traits == SDFT_IMAG_ONLY && real(sig[i]) == 0);
    }

    // generate the phase offsets
	const Float double_pi = static_cast<Float>(2*3.141592653589793238462643383279502884); // Enough precision for everyone!
    for (size_t i = 0; i < number_of_samples; i++) {
        cplx angle(0, double_pi * i / number_of_samples);
        offsets[i] = exp(angle);
    };
}

template<typename Float>
void impl_push_next_sample(struct sdft_State *s, void *next_sample)
{
    typedef std::complex<Float> cplx;
    assert(s->signal_index >= 0 && s->signal_index < s->number_of_samples);

    cplx *signal = (cplx *) s->signal;
    cplx *spec = (cplx *) s->spectrum;
    cplx *po = (cplx *) s->phase_offsets;

    cplx ns = *(cplx *) next_sample;

    assert(s->signal_traits == SDFT_REAL_AND_IMAG
            || s->signal_traits == SDFT_REAL_ONLY && imag(ns) == 0
            || s->signal_traits == SDFT_IMAG_ONLY && real(ns) == 0);

    cplx delta = ns - signal[s->signal_index];

    size_t n_samples = s->number_of_samples;
    size_t n_bins = s->signal_traits == SDFT_REAL_AND_IMAG
            ? s->number_of_samples
            : s->number_of_samples / 2; // only first half of spectrum relevant

    for (size_t i = 0; i < n_bins; ++i) {
        spec[i] = (spec[i] + delta) * po[i];
    }

    signal[s->signal_index] += delta;
    if (++s->signal_index == n_samples)
        s->signal_index = 0;
}

template<typename Float>
void impl_unshift_signal(struct sdft_State *s)
{
    typedef std::complex<Float> cplx;
    cplx *signal = (cplx *) s->signal;
    size_t n = s->number_of_samples;

    if (s->signal_index == 0) {
        // we are lucky, there's nothing to shift
        return;
    }

    // instead of taking the signal index as offset, we take n - signal_index
    // as offset to follow a cycle through the array.
    assert(n > s->signal_index);
    size_t ofs = s->signal_index;

    cplx first = signal[0];
    size_t cur = 0;
    for (size_t i = 0; i < n - 1; ++i) {
        size_t next = (cur + ofs) % n;
        signal[cur] = signal[next];
        cur = next;
    }

    assert((cur + ofs) % n == 0); // or: next == 0
    signal[cur] = first;
}
