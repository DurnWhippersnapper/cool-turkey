extern crate num;
use num::complex::Complex;
use std::num::Zero;
use std::iter::Repeat;

//TODO we can do better than this
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
    let mut temp: Vec<Complex<f32>> = Vec::new();
    temp.grow(signal.len(), Zero::zero());
    cooley_tukey_work(signal, spectrum, temp.as_mut_slice(), 1, factors.as_slice());
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

fn cooley_tukey_work(signal: &[Complex<f32>], spectrum: &mut [Complex<f32>], temp: &mut [Complex<f32>], stride: uint, factors: &[uint])
{
    //TODO do fancy pattern matching on factors
    if factors.len() == 1
    {
        cooley_tukey_base(signal, spectrum, stride, stride);
    }
    else
    {
        let n1 = factors[0];
        //TODO get rid of the ceiling stuff
        let n2 = ((signal.len() as f32) / (n1 as f32) / (stride as f32)).ceil() as uint;
        for i in range(0, n1)
        {
            cooley_tukey_work(signal.slice_from(i * stride), temp.slice_from_mut(i * stride), spectrum.slice_from_mut(i * stride), stride * n1, factors.slice_from(1));
        }

        multiply_by_twiddles(temp, stride, n1, n2);

        for (i, chunk) in temp.chunks(stride * n1).enumerate()
        {
            cooley_tukey_base(chunk, spectrum.slice_from_mut(i * stride), stride, n2 * stride);
        }

    }
}

//TODO make this the function as dft() with use of iterators and generics
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
    for (k, spec_bin) in spectrum.iter_mut().enumerate()
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
