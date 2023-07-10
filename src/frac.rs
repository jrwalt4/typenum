//! Type-level fractions.

use crate::{
    consts::*,
    marker_traits::{Bit, NonZero, Unsigned},
    operator_aliases::{Diff, Gcf, Lcm, PartialQuot, Prod, Sum},
    private::{PrivateSub, PrivateSubOut, Trim, TrimOut},
    type_operators::{Gcd, Lcd, PartialDiv},
    uint::{UInt, UTerm},
};

use core::ops::{Add, Mul, Sub};

/// Type-level unsigned fraction.
///
/// `N` should never be `NonZero`, but is allowed in signature so that UFrac<UTerm, D>
/// can be allowed in intermediate calculation steps.
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug, Default)]
#[cfg_attr(feature = "scale_info", derive(scale_info::TypeInfo))]
pub struct UFrac<N: Unsigned, D: Unsigned + NonZero> {
    pub(crate) n: N,
    pub(crate) d: D,
}

impl<N: Unsigned, D: Unsigned + NonZero> UFrac<N, D> {
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

// `N` must be `NonZero`, otherwise it's UF0.
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
// Trim (a.k.a 'Simplify')

impl Trim for UF0 {
    type Output = UF0;

    #[inline]
    fn trim(self) -> Self::Output {
        self
    }
}

impl<D: Unsigned + NonZero> Trim for UFrac<UTerm, D> {
    type Output = UF0;

    #[inline]
    fn trim(self) -> Self::Output {
        UF0
    }
}

impl<U: Unsigned, B: Bit, D: Unsigned + NonZero> Trim for UFrac<UInt<U, B>, D>
where
    UInt<U, B>: Gcd<D>,
    UInt<U, B>: PartialDiv<Gcf<UInt<U, B>, D>>,
    PartialQuot<UInt<U, B>, Gcf<UInt<U, B>, D>>: Unsigned,
    D: PartialDiv<Gcf<UInt<U, B>, D>>,
    PartialQuot<D, Gcf<UInt<U, B>, D>>: Unsigned + NonZero,
{
    type Output =
        UFrac<PartialQuot<UInt<U, B>, Gcf<UInt<U, B>, D>>, PartialQuot<D, Gcf<UInt<U, B>, D>>>;

    #[inline]
    fn trim(self) -> Self::Output {
        UFrac::new()
    }
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

// ---------------------------------------------------------------------------------------
// Sub

impl Sub for UF0 {
    type Output = UF0;

    #[inline]
    fn sub(self, _rhs: Self) -> Self::Output {
        self
    }
}

impl<N: Unsigned, D: Unsigned + NonZero> Sub<UF0> for UFrac<N, D> {
    type Output = Self;

    #[inline]
    fn sub(self, _rhs: UF0) -> Self::Output {
        self
    }
}

impl<Nl, Dl, Nr, Dr> Sub<UFrac<Nr, Dr>> for UFrac<Nl, Dl>
where
    Nl: Unsigned + NonZero,
    Dl: Unsigned + NonZero,
    Nr: Unsigned + NonZero,
    Dr: Unsigned + NonZero,
    UFrac<Nl, Dl>: PrivateSub<UFrac<Nr, Dr>>,
    PrivateSubOut<UFrac<Nl, Dl>, UFrac<Nr, Dr>>: Trim,
{
    type Output = TrimOut<PrivateSubOut<UFrac<Nl, Dl>, UFrac<Nr, Dr>>>;

    #[inline]
    fn sub(self, rhs: UFrac<Nr, Dr>) -> Self::Output {
        self.private_sub(rhs).trim()
    }
}

impl<Nl, Dl, Nr, Dr> PrivateSub<UFrac<Nr, Dr>> for UFrac<Nl, Dl>
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
    PartialQuot<Prod<Nl, Dr>, Gcf<Dl, Dr>>: Sub<PartialQuot<Prod<Nr, Dl>, Gcf<Dl, Dr>>>,
    Diff<PartialQuot<Prod<Nl, Dr>, Gcf<Dl, Dr>>, PartialQuot<Prod<Nr, Dl>, Gcf<Dl, Dr>>>: Unsigned,
{
    type Output = UFrac<
        Diff<PartialQuot<Prod<Nl, Dr>, Gcf<Dl, Dr>>, PartialQuot<Prod<Nr, Dl>, Gcf<Dl, Dr>>>,
        Lcm<Dl, Dr>,
    >;
    #[inline]
    fn private_sub(self, _rhs: UFrac<Nr, Dr>) -> Self::Output {
        UFrac::new()
    }
}

// ---------------------------------------------------------------------------------------
// Mul

impl Mul for UF0 {
    type Output = UF0;

    #[inline]
    fn mul(self, _rhs: Self) -> Self::Output {
        self
    }
}

impl<N: Unsigned, D: Unsigned + NonZero> Mul<UF0> for UFrac<N, D> {
    type Output = UF0;

    #[inline]
    fn mul(self, rhs: UF0) -> Self::Output {
        rhs
    }
}

impl<N: Unsigned, D: Unsigned + NonZero> Mul<UFrac<N, D>> for UF0 {
    type Output = UF0;

    #[inline]
    fn mul(self, _rhs: UFrac<N, D>) -> Self::Output {
        self
    }
}

impl<Nl, Dl, Nr, Dr> Mul<UFrac<Nr, Dr>> for UFrac<Nl, Dl>
where
    Nl: Unsigned + Mul<Nr>,
    Dl: Unsigned + NonZero + Mul<Dr>,
    Nr: Unsigned,
    Dr: Unsigned + NonZero,
    Prod<Nl, Nr>: Unsigned,
    Prod<Dl, Dr>: Unsigned + NonZero,
    UFrac<Prod<Nl, Nr>, Prod<Dl, Dr>>: Trim,
{
    type Output = TrimOut<UFrac<Prod<Nl, Nr>, Prod<Dl, Dr>>>;

    #[inline]
    fn mul(self, _rhs: UFrac<Nr, Dr>) -> Self::Output {
        Trim::trim(UFrac::<Prod<Nl, Nr>, Prod<Dl, Dr>>::default())
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        assert_type_eq,
        consts::*,
        frac::{UFrac, UnsignedRational, UF0},
        operator_aliases::{Diff, Prod, Sum},
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

    #[test]
    fn ufrac_sub() {
        assert_eq!(0.25_f32, Diff::<UFrac::<U1, U2>, UFrac::<U1, U4>>::F32);
        assert_eq!(0.0_f32, Diff::<UFrac::<U1, U2>, UFrac::<U1, U2>>::F32);
        assert_eq!(0.5_f32, Diff::<UFrac::<U1, U2>, UF0>::F32);
        assert_eq!(0.0_f32, Diff::<UF0, UF0>::F32);
    }

    #[test]
    fn ufrac_mul() {
        assert_eq!(0.25_f32, Prod::<UFrac::<U1, U2>, UFrac::<U1, U2>>::F32);
        assert_eq!(0.125_f32, Prod::<UFrac::<U1, U2>, UFrac::<U1, U4>>::F32);
        assert_eq!(0.1_f32, Prod::<UFrac::<U2, U5>, UFrac::<U1, U4>>::F32);
        assert_type_eq!(UFrac<U1, U10>, Prod<UFrac<U2, U5>, UFrac::<U1, U4>>);
        assert_eq!(0_f32, Prod::<UFrac::<U1, U2>, UF0>::F32);
        assert_eq!(0_f32, Prod::<UF0, UFrac::<U1, U2>>::F32);
    }
}
