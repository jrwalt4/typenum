//! Type-level fractions.

use crate::{
    consts::*,
    int::{NInt, PInt, Z0},
    marker_traits::{Bit, Integer, NonZero, Unsigned},
    operator_aliases::{Compare, Diff, Gcf, Lcm, PartialQuot, Prod, Sum},
    private::{
        Internal, InternalMarker, PrivateRationalAdd, PrivateRationalAddOut, PrivateSub,
        PrivateSubOut, Trim, TrimOut,
    },
    sealed::Sealed,
    type_operators::{Cmp, Gcd, Lcd, PartialDiv},
    uint::{UInt, UTerm},
    Equal, Greater, Less,
};

use core::ops::{Add, Div, Mul, Sub};

/// Type-level unsigned fraction.
///
/// `N` should never be `NonZero`, but is allowed in signature so that UFrac<UTerm, D>
/// can be allowed in intermediate calculation steps.
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug, Default)]
#[cfg_attr(feature = "scale_info", derive(scale_info::TypeInfo))]
pub struct UFrac<N: Unsigned, D: Unsigned + NonZero = U1> {
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
pub trait UnsignedRational: Sealed + Copy + Default {
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

/// Signed rational
pub trait Rational {
    /// Signed numerator
    type Numer: Integer;

    /// Denomenator
    type Denom: Unsigned + NonZero;

    #[allow(missing_docs)]
    const F32: f32;

    #[allow(missing_docs)]
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

/// Positive fraction.
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug, Default)]
#[cfg_attr(feature = "scale_info", derive(scale_info::TypeInfo))]
pub struct PFrac<U: UnsignedRational>(U);

impl<U: UnsignedRational> Rational for PFrac<U>
where
    U::Numer: NonZero,
{
    type Numer = PInt<U::Numer>;
    type Denom = U::Denom;
    const F32: f32 = U::F32;
    const F64: f64 = U::F64;
}

/// Negative fraction.
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug, Default)]
#[cfg_attr(feature = "scale_info", derive(scale_info::TypeInfo))]
pub struct NFrac<U: UnsignedRational>(U);

impl<U: UnsignedRational> Rational for NFrac<U>
where
    U::Numer: NonZero,
{
    type Numer = NInt<U::Numer>;
    type Denom = U::Denom;
    const F32: f32 = -1.0 * U::F32;
    const F64: f64 = -1.0 * U::F64;
}

/// Signed zero fraction.
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug, Default)]
#[cfg_attr(feature = "scale_info", derive(scale_info::TypeInfo))]
pub struct F0;

impl Rational for F0 {
    type Numer = Z0;
    type Denom = U1;
    const F32: f32 = 0.0_f32;
    const F64: f64 = 0.0_f64;
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

impl Trim for F0 {
    type Output = F0;

    #[inline]
    fn trim(self) -> Self::Output {
        self
    }
}

impl<U: UnsignedRational + Trim> Trim for PFrac<U>
where
    TrimOut<U>: UnsignedRational,
{
    type Output = PFrac<TrimOut<U>>;

    #[inline]
    fn trim(self) -> Self::Output {
        PFrac::default()
    }
}

impl<U: UnsignedRational + Trim> Trim for NFrac<U>
where
    TrimOut<U>: UnsignedRational,
{
    type Output = NFrac<TrimOut<U>>;

    #[inline]
    fn trim(self) -> Self::Output {
        NFrac::default()
    }
}

// ---------------------------------------------------------------------------------------
// Cmp

impl Cmp for UF0 {
    type Output = Equal;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &Self) -> Self::Output {
        Equal
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Cmp<UFrac<N, D>> for UF0 {
    type Output = Less;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &UFrac<N, D>) -> Self::Output {
        Less
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Cmp<UF0> for UFrac<N, D> {
    type Output = Greater;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &UF0) -> Self::Output {
        Greater
    }
}

impl<Nl, Dl, Nr, Dr> Cmp<UFrac<Nr, Dr>> for UFrac<Nl, Dl>
where
    Nl: Unsigned + NonZero + Mul<Dr>,
    Dl: Unsigned + NonZero,
    Nr: Unsigned + NonZero + Mul<Dl>,
    Dr: Unsigned + NonZero,
    Prod<Nl, Dr>: Cmp<Prod<Nr, Dl>>,
    Compare<Prod<Nl, Dr>, Prod<Nr, Dl>>: Default,
{
    type Output = Compare<Prod<Nl, Dr>, Prod<Nr, Dl>>;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &UFrac<Nr, Dr>) -> Self::Output {
        Default::default()
    }
}

impl Cmp for F0 {
    type Output = Equal;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &Self) -> Self::Output {
        Equal
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Cmp<PFrac<UFrac<N, D>>> for F0 {
    type Output = Less;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &PFrac<UFrac<N, D>>) -> Self::Output {
        Less
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Cmp<NFrac<UFrac<N, D>>> for F0 {
    type Output = Greater;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &NFrac<UFrac<N, D>>) -> Self::Output {
        Greater
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Cmp<F0> for PFrac<UFrac<N, D>> {
    type Output = Greater;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &F0) -> Self::Output {
        Greater
    }
}

impl<N: Unsigned + NonZero, D: Unsigned + NonZero> Cmp<F0> for NFrac<UFrac<N, D>> {
    type Output = Less;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &F0) -> Self::Output {
        Less
    }
}

impl<Ul: UnsignedRational + Cmp<Ur>, Ur: UnsignedRational> Cmp<PFrac<Ur>> for PFrac<Ul> {
    type Output = Compare<Ul, Ur>;

    #[inline]
    fn compare<IM: InternalMarker>(&self, rhs: &PFrac<Ur>) -> Self::Output {
        self.0.compare::<IM>(&rhs.0)
    }
}

impl<Ul: UnsignedRational, Ur: UnsignedRational> Cmp<NFrac<Ur>> for PFrac<Ul> {
    type Output = Greater;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &NFrac<Ur>) -> Self::Output {
        Greater
    }
}

impl<Ul: UnsignedRational, Ur: UnsignedRational> Cmp<PFrac<Ur>> for NFrac<Ul> {
    type Output = Less;

    #[inline]
    fn compare<IM: InternalMarker>(&self, _: &PFrac<Ur>) -> Self::Output {
        Less
    }
}

impl<Ul: UnsignedRational, Ur: UnsignedRational + Cmp<Ul>> Cmp<NFrac<Ur>> for NFrac<Ul> {
    type Output = Compare<Ur, Ul>;

    #[inline]
    fn compare<IM: InternalMarker>(&self, rhs: &NFrac<Ur>) -> Self::Output {
        rhs.0.compare::<IM>(&self.0)
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

impl Add for F0 {
    type Output = F0;

    #[inline]
    fn add(self, _rhs: F0) -> Self::Output {
        self
    }
}

impl<U: UnsignedRational> Add<PFrac<U>> for F0 {
    type Output = PFrac<U>;

    #[inline]
    fn add(self, rhs: PFrac<U>) -> Self::Output {
        rhs
    }
}

impl<U: UnsignedRational> Add<NFrac<U>> for F0 {
    type Output = NFrac<U>;

    #[inline]
    fn add(self, rhs: NFrac<U>) -> Self::Output {
        rhs
    }
}

impl<U: UnsignedRational> Add<F0> for PFrac<U> {
    type Output = PFrac<U>;

    #[inline]
    fn add(self, _: F0) -> Self::Output {
        self
    }
}

impl<U: UnsignedRational> Add<F0> for NFrac<U> {
    type Output = NFrac<U>;

    #[inline]
    fn add(self, _: F0) -> Self::Output {
        self
    }
}

impl<Ul, Ur> Add<PFrac<Ur>> for PFrac<Ul>
where
    Ul: UnsignedRational + Add<Ur>,
    Ur: UnsignedRational,
    Sum<Ul, Ur>: UnsignedRational,
{
    type Output = PFrac<Sum<Ul, Ur>>;

    #[inline]
    fn add(self, rhs: PFrac<Ur>) -> Self::Output {
        PFrac(self.0 + rhs.0)
    }
}

impl<Ul, Ur> Add<NFrac<Ur>> for NFrac<Ul>
where
    Ul: UnsignedRational + Add<Ur>,
    Ur: UnsignedRational,
    Sum<Ul, Ur>: UnsignedRational,
{
    type Output = NFrac<Sum<Ul, Ur>>;

    #[inline]
    fn add(self, rhs: NFrac<Ur>) -> Self::Output {
        NFrac(self.0 + rhs.0)
    }
}

impl<Ul: UnsignedRational, Ur: UnsignedRational> Add<NFrac<Ur>> for PFrac<Ul>
where
    Ul: Cmp<Ur> + PrivateRationalAdd<Compare<Ul, Ur>, Ur>,
{
    type Output = PrivateRationalAddOut<Ul, Compare<Ul, Ur>, Ur>;

    #[inline]
    fn add(self, rhs: NFrac<Ur>) -> Self::Output {
        let lhs = self.0;
        let rhs = rhs.0;
        let lhs_cmp_rhs = lhs.compare::<Internal>(&rhs);
        lhs.private_rational_add(lhs_cmp_rhs, rhs)
    }
}

impl<Ul: UnsignedRational, Ur: UnsignedRational> Add<PFrac<Ur>> for NFrac<Ul>
where
    Ur: Cmp<Ul> + PrivateRationalAdd<Compare<Ur, Ul>, Ul>,
{
    type Output = PrivateRationalAddOut<Ur, Compare<Ur, Ul>, Ul>;

    #[inline]
    fn add(self, rhs: PFrac<Ur>) -> Self::Output {
        let lhs = self.0;
        let rhs = rhs.0;
        let rhs_cmp_lhs = rhs.compare::<Internal>(&lhs);
        rhs.private_rational_add(rhs_cmp_lhs, lhs)
    }
}

/// `P + N = 0` where `P == N`
impl<N: UnsignedRational, P: UnsignedRational> PrivateRationalAdd<Equal, N> for P {
    type Output = F0;

    #[inline]
    fn private_rational_add(self, _: Equal, _: N) -> Self::Output {
        F0
    }
}

/// `P + N = Positive` where `P > N`
impl<N: UnsignedRational, P: UnsignedRational> PrivateRationalAdd<Greater, N> for P
where
    P: Sub<N>,
    Diff<P, N>: UnsignedRational,
{
    type Output = PFrac<Diff<P, N>>;

    #[inline]
    fn private_rational_add(self, _: Greater, n: N) -> Self::Output {
        PFrac(self - n)
    }
}

/// `P + N = Negative` where `P < N`
impl<N: UnsignedRational, P: UnsignedRational> PrivateRationalAdd<Less, N> for P
where
    N: Sub<P>,
    Diff<N, P>: UnsignedRational,
{
    type Output = NFrac<Diff<N, P>>;

    #[inline]
    fn private_rational_add(self, _: Less, n: N) -> Self::Output {
        NFrac(n - self)
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

// ---------------------------------------------------------------------------------------
// Div

impl<N: Unsigned, D: Unsigned + NonZero> Div<UFrac<N, D>> for UF0 {
    type Output = UF0;

    #[inline]
    fn div(self, _rhs: UFrac<N, D>) -> Self::Output {
        self
    }
}

impl<Nl, Dl, Nr, Dr> Div<UFrac<Nr, Dr>> for UFrac<Nl, Dl>
where
    Nl: Unsigned,
    Dl: Unsigned + NonZero,
    Nr: Unsigned + NonZero,
    Dr: Unsigned + NonZero,
    Self: Mul<UFrac<Dr, Nr>>,
{
    type Output = Prod<Self, UFrac<Dr, Nr>>;

    #[inline]
    fn div(self, _rhs: UFrac<Nr, Dr>) -> Self::Output {
        self * UFrac::<Dr, Nr>::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        assert_type_eq,
        consts::*,
        frac::{NFrac, PFrac, Rational, UFrac, UnsignedRational, F0, UF0},
        operator_aliases::{Diff, Prod, Quot, Sum},
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

    #[test]
    fn ufrac_div() {
        assert_eq!(1_f32, Quot::<UFrac::<U1, U2>, UFrac::<U1, U2>>::F32);
        assert_eq!(0.25_f32, Quot::<UFrac::<U1, U2>, UFrac::<U2>>::F32);
        assert_eq!(0.8_f32, Quot::<UFrac::<U2, U5>, UFrac::<U1, U2>>::F32);
        assert_type_eq!(UFrac<U4, U5>, Quot<UFrac<U2, U5>, UFrac::<U2, U4>>);
        assert_eq!(0_f32, Quot::<UF0, UFrac::<U1, U2>>::F32);
    }

    #[test]
    fn frac_add() {
        assert_eq!(0_f32, Sum::<F0, F0>::F32);
        assert_eq!(0.5_f32, Sum::<F0, PFrac<UFrac<U1, U2>>>::F32);
        assert_eq!(0.5_f32, Sum::<PFrac<UFrac<U1, U2>>, F0>::F32);
        assert_eq!(
            1_f32,
            Sum::<PFrac<UFrac<U1, U2>>, PFrac<UFrac<U1, U2>>>::F32
        );
        assert_eq!(
            -2_f32,
            Sum::<NFrac<UFrac<U3, U2>>, NFrac<UFrac<U1, U2>>>::F32
        );
        assert_eq!(
            0_f32,
            Sum::<PFrac<UFrac<U1, U2>>, NFrac<UFrac<U1, U2>>>::F32
        );
        assert_eq!(0.5_f32, Sum::<PFrac<UFrac<U3, U2>>, NFrac<UFrac<U1>>>::F32);
        assert_eq!(-0.5_f32, Sum::<PFrac<UFrac<U1, U2>>, NFrac<UFrac<U1>>>::F32);
    }
}
