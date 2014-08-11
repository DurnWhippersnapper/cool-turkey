extern crate num;
use num::complex::Complex;
use std::num::Zero;

/*
    Regular O(n^2) DFT calculation
*/
pub fn dft(signal: &[Complex<f32>], spectrum: &mut [Complex<f32>])
{
    for (k, spec_bin) in spectrum.mut_iter().enumerate()
    {
        let mut sum: Complex<f32> = Zero::zero();
        for (i, &x) in signal.iter().enumerate()
        {
            let angle = -1f32 * (i as f32) * (k as f32) * Float::two_pi() / (signal.len() as f32);
            let twiddle: Complex<f32> = Complex::from_polar(&1f32, &angle);
            sum = sum + twiddle * x;
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
