
pub mod bit;

pub mod uint;
pub mod const_uints;

pub mod int;

pub mod __private;

pub trait Same<Rhs = Self> {
    /// `Output` should always be `Self`
    type Output;
}

pub trait Not {
    type Output;
}
pub trait Neg {
    type Output;
}
pub trait And<Rhs = Self> {
    type Output;
}
pub trait Or<Rhs = Self> {
    type Output;
}
pub trait Xor<Rhs = Self> {
    type Output;
}

pub trait Shl<Rhs = Self> {
    type Output;
}
pub trait Shr<Rhs = Self> {
    type Output;
}

pub trait Add<Rhs = Self> {
    type Output;
}
pub trait Sub<Rhs = Self> {
    type Output;
}
pub trait Mul<Rhs = Self> {
    type Output;
}
pub trait Div<Rhs = Self> {
    type Output;
}
pub trait Rem<Rhs = Self> {
    type Output;
}

/// Gives the size of a type number in bits as a `UInt`
pub trait SizeOf {
    type Output;
}
