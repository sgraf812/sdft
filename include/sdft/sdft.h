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
* \brief Error codes returned by the functions of this library.
*/
enum sdft_Error {
    /**
    * The function executed without producing errors.
    */
    SDFT_NO_ERROR,
    /**
    * The passed window size was too small (e.g. < 1) to be useful.
    */
    SDFT_WINDOW_TOO_SHORT,
    /**
    * One of the samples of the signal didn't match the desired signal trait.
    */
    SDFT_SIGNAL_TRAIT_VIOLATION,
    /**
    * The two sdft_State structs passed to sdft_init_combine were not combinable due to
    * different configuration parameters (e.g. different window sizes or signal traits).
    */
    SDFT_NOT_COMBINABLE,
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
* \param window the sliding window buffer containing the initial signal. At least number_of_samples complex elements.
* \param spectrum the buffer containing the initial spectrum. At least number_of_samples
*        (signal_traits == SDFT_REAL_AND_IMAG) or number_of_samples (signal_traits != SDFT_REAL_AND_IMAG)
*        complex elements, depending on the signal_trait which is guaranteed by the user.
* \param phase_offsets a buffer for internal use whose content will be overwritten. At least number_of_samples
*        complex elements.
* \param window_size the number of samples in the sliding window buffer and also the number of bins in the spectrum.
* \param signal_traits can pass guarantees to the SDFT about the signal which can be exploited.
*        In particular, purely real and purely imaginary signals can be represented by a spectrum of half the
*        sample length because of the implied symmetries in the frequency domain.
*        By specifying either SDFT_REAL_ONLY or SDFT_IMAG_ONLY, only half of the spectrum gets updated in each
*        time step. As such, the upper half of spectrum does not contain usable values, but can be reconstructed
*        by the user exploiting the symmetries.
*        Of course, SDFT_REAL_AND_IMAG will always be the conservative but correct choice and should be used if
*        not certain of the traits of the signal.
* \returns an error code indicating success or failure.
*          SDFT_WINDOW_TOO_SHORT: The window_size was too small (e.g. < 1).
*          SDFT_SIGNAL_TRAIT_VIOLATION: The initial values in the window buffer didn't match the desired signal trait.
* \see sdft_State, sdft_SignalTraits, sdft_float
*
* Runtime: O(window_size)
*/
enum sdft_Error sdft_init_from_buffers(
        struct sdft_State *state,
        enum sdft_FloatPrecision precision,
        void *window,
        void *spectrum,
        void *phase_offsets,
        size_t window_size,
        enum sdft_SignalTraits signal_traits);

/**
* \brief Combines two combinable sdft_State structs (the buffers of which must not overlap) for vastly increased
*        numberical stability.
*
* Two states are combinable iff their floating point precision, signal trait and window sizes match.
* Upon successful initialization, the first struct is assumed to contain the initial values and the buffers
* of the second struct are reset.
*
* The two states are combined in such a way that as soon as one buffer reaches a steady-state (e.g. the window is
* filled with the last window_size samples of the signal), the other buffer is reset. Thus, no state will ever aggregate
* errors for more than 2*window_size samples. Some ASCII art depicting the number of samples in the windows of each
* state:
*
*
*               reset first
*                    |
*              /-----+----\
*              |          |
*              v          v
*
*         -----     /-----     /-
*         |   |    /     |    /
*         |   |   /      |   /
*         |   |  /       |  /
*         |   | /        | /
* first:  -----------------------
*
*              /-----     /-----
*             /     |    /     |
*            /      |   /      |
*           /       |  /       |
*          /        | /        |
* second: -----------------------
*
*                    ʌ          ʌ
*                    |          |
*                    \-----+----/
*                          |
*                    reset second
*
* At all times we get have a valid window stored in either one of the states. The reset forces the spectrum and
* the window to be in sync (e.g. both are zero'ed out), in this way numerical error's can't propagate endlessly.
*
* \param state the allocated sdft_State struct which is initialized by this function.
* \param first the first sdft_State struct. Must be combinable with second.
* \param second the second sdft_State struct. Must be combinable with first.
*
* \returns an error code indicating success or failure.
*          SDFT_NOT_COMBINABLE: second and first were not combinable because of a mismatch in precision, signal traits
*                               or window sizes.
*
* Runtime: O(1)
*/
enum sdft_Error sdft_init_combine(
        struct sdft_State *state,
        struct sdft_State *first,
        struct sdft_State *second);

/**
* \brief Pushes a fresh sample through the sDFT and updates the spectrum, which is immediately usable.
*
* The internal window buffer(s) does not necessarily contain the actual samples in the right order.
* The ordering of the window can be restored with sdft_unshift_and_get_window, which
* does NOT need to be called prior to sdft_push_next_sample.
*
* \param state the internal state, which has to be initialized with sdft_init prior to usage.
* \param next_sample a pointer the next sample, which is assumed to be a single complex number in the format
*        described in sdft_init.
*
* Runtime: O(window_size)
*/
enum sdft_Error sdft_push_next_sample(struct sdft_State *state, void *next_sample);

/**
* \brief Returns a pointer to the current spectrum buffer.
*
* For same states, the result of this function may change after invocations of sdft_push_next_sample if state
* is initialized as combined.
*/
void *sdft_get_spectrum(struct sdft_State *state);

/**
* \brief Restores the temporal ordering of the window buffer after having called sdft_push_next_sample one or more times
*        and returns it.
*
* The returned window pointer may change between invocations if sdft_push_next_sample was called!
*
* \returns The current window buffer with its samples ordered temporally.
*
* Runtime: O(window_size)
*/
void *sdft_unshift_and_get_window(struct sdft_State *state);

#ifdef __cplusplus
};
#endif