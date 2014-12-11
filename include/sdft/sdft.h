#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
* \brief An opaque state struct which is internally used for preservation of state
* between sdft calls.
*/
struct sdft_State;

/**
* \brief Each value describes an available floating point precision the SDFT can be run with.
*/
enum sdft_FloatPrecision {
    SDFT_SINGLE,
    SDFT_DOUBLE,
    SDFT_LONG_DOUBLE
};

/**
* \brief An enumeration describing the traits of the signal wrt. real and imaginary parts.
*
* Is used as a parameter to sdft_init. If it is clear in beforehand that the signal is purely
* imaginary or purely real, sdft_push_next_sample can perform 50% less work because of the
* symmetry of the spectrum.
*/
enum sdft_SignalTraits {
    SDFT_REAL_AND_IMAG,
    SDFT_REAL_ONLY,
    SDFT_IMAG_ONLY
};

/**
* \brief Returns the size of the sdft_State struct which has to be allocated by the
*        user of the library.
*/
size_t sdft_size_of_state();

/**
* \brief Initializes the sdft_State with user allocated buffers and the guaranteed signal traits.
*
* As the consequence of guaranteeing no internal heap allocation, it is the user's task to provide
* buffers for the signal, spectrum and phase_offset. Each of the buffers will contain complex numbers,
* which are binary compatible (in fact are equal) to C99's complex number types, handled as an array of 2
* floating point numbers each with appropriate precision.
*
* \param state the state which is to be initialized.
* \param precision the precision to use in floating point operations and buffers. Affects the buffer sizes.
* \param signal the buffer containing the initial signal. At least number_of_samples complex elements.
* \param spectrum the buffer containing the initial spectrum. At least number_of_samples
*        (signal_traits == SDFT_REAL_AND_IMAG) or number_of_samples (signal_traits != SDFT_REAL_AND_IMAG)
*        complex elements, depending on the signal_trait which is guaranteed by the user.
* \param phase_offsets a buffer for internal use whose content will be overwritten. At least number_of_samples
*        complex elements.
* \param number_of_samples the number of samples in signal and also the number of bins in spectrum.
* \param signal_traits can pass guarantees to the SDFT about the signal which can be exploited.
*        In particular, purely real and purely imaginary signals can be represented by a spectrum of half the
*        sample length because of the implied symmetries in the frequency domain.
*        By specifying either SDFT_REAL_ONLY or SDFT_IMAG_ONLY, only half of the spectrum gets updated in each
*        time step. As such, the upper half of spectrum does not contain usable values, but can be reconstructed
*        by the user exploiting the symmetries.
*        Of course, SDFT_REAL_AND_IMAG will always be the conservative but correct choice and should be used if
*        not certain of the traits of the signal.
* \see sdft_State, sdft_SignalTraits, sdft_float
*/
void sdft_init(
        struct sdft_State *state,
        enum sdft_FloatPrecision precision,
        void *signal,
        void *spectrum,
        void *phase_offsets,
        size_t number_of_samples,
        enum sdft_SignalTraits signal_traits);

/**
* \brief Pushes a fresh sample through the sDFT and updates the spectrum, which is immediately usable.
*
* The internal signal buffer does not necessarily contain the actual signal in the right order.
* The initial ordering of the signal buffer can be restored with sdft_unshift_signal, which
* is NOT needed for sdft_push_next_sample.
*
* \param state the internal state, which has to be initialized with sdft_init prior to usage.
* \param next_sample a pointer the next sample, which is assumed to be a single complex number in the format
*        described in sdft_init.
*
* Runtime: O(number_of_samples)
*/
void sdft_push_next_sample(struct sdft_State *state, void *next_sample);

/**
* \brief Restores the initial ordering of the signal buffer after having called sdft_push_next_sample one or more times.
*
* Should only be done for inspecting the stored signal values.
*
* Runtime: O(number_of_samples)
*/
void sdft_unshift_signal(struct sdft_State *state);

#ifdef __cplusplus
};
#endif