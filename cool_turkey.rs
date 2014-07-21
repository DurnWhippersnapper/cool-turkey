extern crate num;
use num::complex::Complex;
use std::num::Zero;
use std::num::One;
use std::num::FromPrimitive;

/*
trait MultipliesComplex<A>
{
    fn complex_mult<A>(&self, rhs: &Complex<A>) -> Complex<A>;
}

impl<A: Clone + Num> MultipliesComplex<A> for Complex<A>
{
    fn complex_mult<A>(&self, rhs: &Complex<A>)
    {
        self * rhs
    }
}
*/

/*
    Regular O(n^2) DFT calculation
*/
pub fn dft<T: Num + Clone + FloatMath + FromPrimitive, U: Mul<Complex<T>, Complex<T>> >(signal: &[U], spectrum: &mut [Complex<T>])
{
    for (k, spec_bin) in spectrum.mut_iter().enumerate()
    {
        let mut sum: Complex<T> = Zero::zero();
        let mut angle: T = Zero::zero();
        // TODO figure out how to make this generic
        //let rad_per_sample = k * Float::two_pi() / signal.len();
        let rad_per_sample : T = One::one();
        for x in signal.iter()
        {
            let one: T = One::one();
            let twiddle: Complex<T> = Complex::from_polar(&one, &angle);
            sum = sum + (*x) * twiddle;
            angle = angle - rad_per_sample;
        }
        *spec_bin = sum;
    }
}

/*
    Testing
*/
fn compare_vectors(vec1: &[Complex<f32>], vec2: &[Complex<f32>]) -> bool
{
    let mut sse = 0f32;
    for (&a, &b) in vec1.iter().zip(vec2.iter())
    {
        sse = sse + (a - b).norm();
    }
    return sse < 1f32;
}

#[test]
fn dft_test()
{

    let signal = vec![Complex{re: 1f32, im: 0f32},
                      Complex{re:-1f32, im: 0f32}];
    let mut spectrum = signal.to_owned();
    dft(signal.as_slice(), spectrum.as_mut_slice());
    assert!(compare_vectors(spectrum.as_slice(), vec![Complex{re: 0f32, im: 0f32}, 
                                                      Complex{re: 2f32, im: 0f32}].as_slice()));

    let signal = vec![Complex{re: 1f32, im: 1f32},
                      Complex{re: 2f32, im:-3f32},
                      Complex{re:-1f32, im: 4f32}];
    let mut spectrum = signal.to_owned();
    dft(signal.as_slice(), spectrum.as_mut_slice());
    assert!(compare_vectors(spectrum.as_slice(), vec![Complex{re: 2f32, im: 2f32}, 
                                                      Complex{re: -5.562177f32, im: -2.098076f32},
                                                      Complex{re: 6.562178f32, im: 3.09807f32}].as_slice()));

    let signal = vec![Complex{re: 0f32, im: 1f32},
                      Complex{re: 2.5f32, im:-3f32},
                      Complex{re:-1f32, im: -1f32},
                      Complex{re: 4f32, im: 0f32}];
    let mut spectrum = signal.to_owned();
    dft(signal.as_slice(), spectrum.as_mut_slice());
    assert!(compare_vectors(spectrum.as_slice(), vec![Complex{re: 5.5f32, im: -3f32}, 
                                                      Complex{re: -2f32, im: 3.5f32},
                                                      Complex{re: -7.5f32, im: 3f32},
                                                      Complex{re: 4f32, im: 0.5f32}].as_slice()));

}
