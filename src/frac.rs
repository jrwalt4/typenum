//! Type-level fractions.

use crate::{
    consts::*,
    marker_traits::{NonZero, Unsigned},
    operator_aliases::{Gcf, Lcm, PartialQuot, Prod, Sum},
    type_operators::{Gcd, Lcd, PartialDiv},
};

use core::ops::{Add, Mul};

/// Type-level unsigned fraction.
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug, Default)]
#[cfg_attr(feature = "scale_info", derive(scale_info::TypeInfo))]
pub struct UFrac<N: Unsigned + NonZero, D: Unsigned + NonZero> {
    pub(crate) n: N,
    pub(crate) d: D,
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> UFrac<N, D> {
    /// Instantiates a singleton representing this unsigned fraction.
    #[inline]
    pub fn new() -> UFrac<N, D> {
        UFrac::default()
    }
}

/// Unsigned zero fraction.
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug, Default)]
#[cfg_attr(feature = "scale_info", derive(scale_info::TypeInfo))]
pub struct UF0;

impl UF0 {
    /// Instantiates a singleton representing zero fraction.
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }
}

/// Marker trait for unsigned rational numbers.
pub trait UnsignedRational {
    /// The reduced numerator of a rational number.
    type Numer: Unsigned;

    /// The reduced denominator of a rational number.
    ///
    /// Could possibly allow for Zero denominator in the future
    /// to support +/- Infinity.
    type Denom: Unsigned + NonZero;

    /// The 32-bit floating point representation of this rational number.
    const F32: f32;

    /// The 64-bit floating point representation of this rational number.
    const F64: f64;

    #[allow(missing_docs)]
    #[inline]
    fn to_f32() -> f32 {
        Self::F32
    }

    #[allow(missing_docs)]
    #[inline]
    fn to_f64() -> f64 {
        Self::F64
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> UnsignedRational for UFrac<N, D> {
    type Numer = N;
    type Denom = D;

    const F32: f32 = <N as Unsigned>::I32 as f32 / <D as Unsigned>::I32 as f32;

    const F64: f64 = <N as Unsigned>::I64 as f64 / <D as Unsigned>::I64 as f64;
}

impl UnsignedRational for UF0 {
    type Numer = U0;
    type Denom = U1;

    const F32: f32 = 0.0;

    const F64: f64 = 0.0;
}

// ---------------------------------------------------------------------------------------
// Add

impl Add<UF0> for UF0 {
    type Output = UF0;

    #[inline]
    fn add(self, _rhs: UF0) -> Self::Output {
        self
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Add<UF0> for UFrac<N, D> {
    type Output = Self;

    #[inline]
    fn add(self, _rhs: UF0) -> Self::Output {
        self
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Add<UFrac<N, D>> for UF0 {
    type Output = UFrac<N, D>;

    #[inline]
    fn add(self, rhs: UFrac<N, D>) -> Self::Output {
        rhs
    }
}

impl<Nl, Dl, Nr, Dr> Add<UFrac<Nr, Dr>> for UFrac<Nl, Dl>
where
    Nl: Unsigned + NonZero,
    Dl: Unsigned + NonZero + Gcd<Dr>,
    Nr: Unsigned + NonZero,
    Dr: Unsigned + NonZero,
    Dl: Lcd<Dr>,
    Gcf<Dl, Dr>: Unsigned + NonZero,
    Lcm<Dl, Dr>: Unsigned + NonZero,
    Nl: Mul<Dr>,
    Prod<Nl, Dr>: PartialDiv<Gcf<Dl, Dr>>,
    Nr: Mul<Dl>,
    Prod<Nr, Dl>: PartialDiv<Gcf<Dl, Dr>>,
    PartialQuot<Prod<Nl, Dr>, Gcf<Dl, Dr>>: Add<PartialQuot<Prod<Nr, Dl>, Gcf<Dl, Dr>>>,
    Sum<PartialQuot<Prod<Nl, Dr>, Gcf<Dl, Dr>>, PartialQuot<Prod<Nr, Dl>, Gcf<Dl, Dr>>>:
        Unsigned + NonZero,
{
    type Output = UFrac<
        Sum<PartialQuot<Prod<Nl, Dr>, Gcf<Dl, Dr>>, PartialQuot<Prod<Nr, Dl>, Gcf<Dl, Dr>>>,
        Lcm<Dl, Dr>,
    >;
    #[inline]
    fn add(self, _rhs: UFrac<Nr, Dr>) -> Self::Output {
        UFrac::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        consts::*,
        frac::{UFrac, UnsignedRational, UF0},
        operator_aliases::Sum,
    };

    #[test]
    fn unsigned_rational() {
        assert_eq!(0.5_f32, UFrac::<U1, U2>::F32);
    }

    #[test]
    fn ufrac_add() {
        assert_eq!(0.5_f32, Sum::<UFrac::<U1, U2>, UF0>::F32);
        assert_eq!(0.5_f32, Sum::<UF0, UFrac::<U1, U2>>::F32);
    }
}
