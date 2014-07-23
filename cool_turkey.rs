extern crate num;
use num::complex::Complex;
use std::num::Zero;
use std::iter::Repeat;

fn factor(n: uint) -> Vec<uint>
{
    let mut factors: Vec<uint> = Vec::new();
    let mut next = n;
    while next > 1
    {
        for div in Repeat::new(2u).take(1).chain(std::iter::range_step_inclusive(3u, next, 2))
        {
            if next % div == 0
            {
                next = next / div;
                factors.push(div);
                break;
            }
        }
    }
    return factors;
}

pub fn cooley_tukey(signal: &[Complex<f32>], spectrum: &mut [Complex<f32>])
{
    let factors = factor(signal.len());
    cooley_tukey_work(signal, spectrum, 1, factors.as_slice());
}

fn multiply_by_twiddles(xs: &mut [Complex<f32>], stride: uint, n1: uint, n2: uint)
{
    for k2 in std::iter::range(0, n2)
    {
        for k1 in std::iter::range(0, n1)
        {
            let angle: f32 = (-1i as f32) * ((k2 * k1) as f32) * Float::two_pi() / ((n1 * n2) as f32);
            let twiddle: Complex<f32> = Complex::from_polar(&1f32, &angle);
            let idx = (k2 * n1 + k1) * stride;
            xs[idx] = xs[idx] * twiddle;
        }
    }
}

fn cooley_tukey_work(signal: &[Complex<f32>], spectrum: &mut [Complex<f32>], stride: uint, factors: &[uint])
{
    if factors.len() == 1
    {
        cooley_tukey_base(signal, spectrum, stride, stride);
    }
    else
    {
        let n1 = factors[0];
        let n2 = signal.len() / n1;
        for i in range(0, n1)
        {
            cooley_tukey_work(signal.slice_from(i), spectrum.mut_slice_from(i), stride * n1, factors.slice_from(1));
        }

        multiply_by_twiddles(spectrum, stride, n1, n2);

        let spectrum_copy = spectrum.to_owned();
        let spectrum_chunks = spectrum_copy.as_slice().chunks(stride * n1);

        for (i, chunk) in spectrum_chunks.enumerate()
        {
            println!("{0}, {1}", i, chunk);
            cooley_tukey_base(chunk, spectrum.mut_slice_from(i), stride, n2 * stride);
        }

    }
}

fn cooley_tukey_base(signal: &[Complex<f32>], spectrum: &mut [Complex<f32>], signal_stride: uint, spectrum_stride: uint)
{
    let idxs = std::iter::range_step(0, spectrum.len(), spectrum_stride);
    for (k, idx) in idxs.enumerate()
    {
        let mut sum: Complex<f32> = Zero::zero();
        let mut angle = 0f32;
        // TODO get rid of this ceiling crap
        let rad_per_sample = (k as f32) * Float::two_pi() / ((signal.len() as f32) / (signal_stride as f32)).ceil();
        for signal_idx in std::iter::range_step(0, signal.len(), signal_stride)
        {
            let twiddle: Complex<f32> = Complex::from_polar(&1f32, &angle);
            sum = sum + twiddle * signal[signal_idx];
            angle = angle - rad_per_sample;
        }
        spectrum[idx] = sum;
    }
}

/*
    Regular O(n^2) DFT calculation
*/
pub fn dft(signal: &[Complex<f32>], spectrum: &mut [Complex<f32>])
{
    for (k, spec_bin) in spectrum.mut_iter().enumerate()
    {
        let mut sum: Complex<f32> = Zero::zero();
        let mut angle = 0f32;
        let rad_per_sample = (k as f32) * Float::two_pi() / (signal.len() as f32);
        for &x in signal.iter()
        {
            let twiddle: Complex<f32> = Complex::from_polar(&1f32, &angle);
            sum = sum + twiddle * x;
            angle = angle - rad_per_sample;
        }
        *spec_bin = sum;
    }
}

/*
    Testing
*/
fn sse_vectors(vec1: &[Complex<f32>], vec2: &[Complex<f32>]) -> f32
{
    let mut sse = 0f32;
    for (&a, &b) in vec1.iter().zip(vec2.iter())
    {
        sse = sse + (a - b).norm();
    }
    return sse;
}

#[test]
fn dft_test()
{

    let signal = vec![Complex{re: 1f32, im: 0f32},
                      Complex{re:-1f32, im: 0f32}];
    let mut spectrum = signal.to_owned();
    dft(signal.as_slice(), spectrum.as_mut_slice());
    assert!(1f32 > sse_vectors(spectrum.as_slice(), vec![Complex{re: 0f32, im: 0f32}, 
                                                      Complex{re: 2f32, im: 0f32}].as_slice()));

    let signal = vec![Complex{re: 1f32, im: 1f32},
                      Complex{re: 2f32, im:-3f32},
                      Complex{re:-1f32, im: 4f32}];
    let mut spectrum = signal.to_owned();
    dft(signal.as_slice(), spectrum.as_mut_slice());
    assert!(1f32 > sse_vectors(spectrum.as_slice(), vec![Complex{re: 2f32, im: 2f32}, 
                                                      Complex{re: -5.562177f32, im: -2.098076f32},
                                                      Complex{re: 6.562178f32, im: 3.09807f32}].as_slice()));

    let signal = vec![Complex{re: 1f32, im: 1f32},
                      Complex{re: 2f32, im: 2f32},
                      Complex{re: 3f32, im: 3f32},
                      Complex{re: 4f32, im: 4f32},
                      Complex{re: 5f32, im: 5f32},
                      Complex{re: 6f32, im: 6f32}];
    let mut spectrum = signal.to_owned();
    dft(signal.as_slice(), spectrum.as_mut_slice());
    assert!(1f32 > sse_vectors(spectrum.as_slice(), vec![Complex{re: 21f32, im: 21f32}, 
                                                      Complex{re: -8.16f32, im: 2.16f32},
                                                      Complex{re: -4.76f32, im: -1.24f32},
                                                      Complex{re: -3f32, im: -3f32},
                                                      Complex{re: -1.24f32, im: -4.76f32},
                                                      Complex{re: 2.16f32, im: -8.16f32}].as_slice()));

    let signal = vec![Complex{re: 0f32, im: 1f32},
                      Complex{re: 2.5f32, im:-3f32},
                      Complex{re:-1f32, im: -1f32},
                      Complex{re: 4f32, im: 0f32}];
    let mut spectrum = signal.to_owned();
    dft(signal.as_slice(), spectrum.as_mut_slice());
    assert!(1f32 > sse_vectors(spectrum.as_slice(), vec![Complex{re: 5.5f32, im: -3f32}, 
                                                      Complex{re: -2f32, im: 3.5f32},
                                                      Complex{re: -7.5f32, im: 3f32},
                                                      Complex{re: 4f32, im: 0.5f32}].as_slice()));

}

#[test]
fn cooley_tukey_test_4()
{
    let signal = vec![Complex{re: 0f32, im: 1f32},
                      Complex{re: 2.5f32, im:-3f32},
                      Complex{re:-1f32, im: -1f32},
                      Complex{re: 4f32, im: 0f32}];
    let mut spectrum_dft = signal.to_owned();
    let mut spectrum_ct = signal.to_owned();
    dft(signal.as_slice(), spectrum_dft.as_mut_slice());
    cooley_tukey(signal.as_slice(), spectrum_ct.as_mut_slice());
    assert!(1f32 > sse_vectors(spectrum_dft.as_slice(),
                               spectrum_ct.as_slice()));
}

#[test]
fn cooley_tukey_test_6()
{
    let signal = vec![Complex{re: 1f32, im: 1f32},
                      Complex{re: 2f32, im: 2f32},
                      Complex{re: 3f32, im: 3f32},
                      Complex{re: 4f32, im: 4f32},
                      Complex{re: 5f32, im: 5f32},
                      Complex{re: 6f32, im: 6f32}];
    let mut spectrum_dft = signal.to_owned();
    let mut spectrum_ct = signal.to_owned();
    dft(signal.as_slice(), spectrum_dft.as_mut_slice());
    cooley_tukey(signal.as_slice(), spectrum_ct.as_mut_slice());
    assert!(1f32 > sse_vectors(spectrum_dft.as_slice(),
                               spectrum_ct.as_slice()));
}

#[test]
fn cooley_tukey_test_nofactor()
{
    let signal = vec![Complex{re: 0f32, im: 1f32},
                      Complex{re:-1f32, im:-1f32},
                      Complex{re: 4f32, im: 0f32}];
    let mut spectrum_dft = signal.to_owned();
    let mut spectrum_ct = signal.to_owned();
    dft(signal.as_slice(), spectrum_dft.as_mut_slice());
    cooley_tukey(signal.as_slice(), spectrum_ct.as_mut_slice());
    assert_eq!(spectrum_dft, spectrum_ct);
}
