#![feature(macro_rules)]
use std::kinds::marker;
use std::mem;
use std::mem::transmute;
use std::slice::Items;

#[cfg(test)]
mod stride_items_test;

pub struct StrideItems<'a, T>
{
    ptr: *const T,
    end: *const T,
    stride: uint,
    marker: marker::ContravariantLifetime<'a>
}

pub struct MutStrideItems<'a, T>
{
    ptr: *mut T,
    end: *mut T,
    stride: uint,
    marker: marker::ContravariantLifetime<'a>,
    marker2: marker::NoCopy
}

// The shared definition of the 'StrideItems' and 'MutStrideItems' iterators
macro_rules! iterator{
    (struct $name:ident -> $ptr:ty, $elem:ty) => {
        impl<'a, T> Iterator<$elem> for $name<'a, T> {
            #[inline]
            fn next(&mut self) -> Option<$elem> {
                // could be implemented with slices, but this avoids bounds checks
                unsafe {
                    if self.ptr == self.end {
                        None
                    } else {
                        if mem::size_of::<T>() == 0 {
                            // purposefully don't use 'ptr.offset' because for
                            // vectors with 0-size elements this would return the
                            // same pointer.
                            self.ptr = transmute(self.ptr as uint + self.stride);

                            // Use a non-null pointer value
                            Some(transmute(1u))
                        } else {
                            let old = self.ptr;
                            self.ptr = self.ptr.offset(self.stride as int);

                            Some(transmute(old))
                        }
                    }
                }
            }

            #[inline]
            fn size_hint(&self) -> (uint, Option<uint>) {
                let diff = (self.end as uint) - (self.ptr as uint);
                let size = mem::size_of::<T>();
                let exact = diff / (self.stride * (if size == 0 {1} else {size}));
                (exact, Some(exact))
            }
        }

        impl<'a, T> DoubleEndedIterator<$elem> for $name<'a, T> {
            #[inline]
            fn next_back(&mut self) -> Option<$elem> {
                // could be implemented with slices, but this avoids bounds checks
                unsafe {
                    if self.end <= self.ptr {
                        None
                    } else {
                        if mem::size_of::<T>() == 0 {
                            // See above for why 'ptr.offset' isn't used
                            self.end = transmute(self.end as uint - self.stride);

                            // Use a non-null pointer value
                            Some(transmute(1u))
                        } else {
                            self.end = self.end.offset(-1 * self.stride as int);

                            Some(transmute(self.end))
                        }
                    }
                }
            }
        }
    }
}

iterator!{struct StrideItems -> *const T, &'a T}
iterator!{struct MutStrideItems -> *mut T, &'a mut T}
impl<'a, T> ExactSize<&'a T> for StrideItems<'a, T>{}
impl<'a, T> ExactSize<&'a mut T> for MutStrideItems<'a, T>{}

// Composing Stride Iterators
fn add_stride<'a, T>(orig: & StrideItems<'a, T>, stride: uint) -> StrideItems<'a, T>
{
    StrideItems{ptr: orig.ptr,
                end: orig.end,
                marker: orig.marker,
                stride: orig.stride * stride}
}

fn add_mut_stride<'a, T>(orig: & MutStrideItems<'a, T>, stride: uint) -> MutStrideItems<'a, T>
{
    MutStrideItems{ptr: orig.ptr,
                   end: orig.end,
                   marker: orig.marker,
                   marker2: marker::NoCopy,
                   stride: orig.stride * stride}
}

// Creating StrideItems from Slices
// should we take an offset as well?
fn stride<'a, T>(items: &[T], stride: uint) -> StrideItems<'a, T>
{
    unsafe
    {
        let p = items.as_ptr();
        if mem::size_of::<T>() == 0 {
            StrideItems{ptr: p,
                        end: (p as uint + items.len()) as *const T,
                        stride: stride,
                        marker: marker::ContravariantLifetime::<'a>}
        } else {
            StrideItems{ptr: p,
                        end: p.offset(items.len() as int),
                        stride: stride,
                        marker: marker::ContravariantLifetime::<'a>}
        }
    }
}

fn mut_stride<'a, T>(items: &mut [T], stride: uint) -> MutStrideItems<'a, T>
{
    unsafe
    {
        let p = items.as_mut_ptr();
        if mem::size_of::<T>() == 0 {
            MutStrideItems{ptr: p,
                           end: (p as uint + items.len()) as *mut T,
                           stride: stride,
                           marker: marker::ContravariantLifetime::<'a>,
                           marker2: marker::NoCopy}
        } else {
            MutStrideItems{ptr: p,
                           end: p.offset(items.len() as int),
                           stride: stride,
                           marker: marker::ContravariantLifetime::<'a>,
                           marker2: marker::NoCopy}
        }
    }
}
