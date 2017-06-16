//fftwcpp_templates.h
//Marc Brooker mbrooker@rrsg.ee.uct.ac.za
//12 June 2006

//Preprocessor magic to link to the right fftw precision
//C++ experts consider preprocessor magic to be bad. I used it here because the alternative was a mess of templates.
#if RS_FLOAT_LONG_DOUBLE == 1
#define fftw_ fftwf_
#elif RS_FLOAT_FLOAT == 1
#define fftw_ fftwl_
#endif
