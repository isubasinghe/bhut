//! Interpolated brightness on a floating-point-spaced grid; exact mode remains available.
use std::sync::OnceLock;

pub const MIN_BITS: u32 = 107 << 23; // 2^-20, below the smallest nonzero blurred count.
pub const SHIFT: u32 = 13; // 1024 intervals per power of two.
pub const TABLE_LEN: usize = (((159 << 23) - MIN_BITS) >> SHIFT) as usize + 1; // through 2^32

pub fn table() -> &'static [f32] {
    static TABLE: OnceLock<Vec<f32>> = OnceLock::new();
    TABLE.get_or_init(|| {
        (0..TABLE_LEN)
            .map(|i| {
                let total = f32::from_bits(MIN_BITS + ((i as u32) << SHIFT));
                (10.0 * total).ln_1p().powf(0.6)
            })
            .collect()
    })
}

pub struct Brightness {
    exposure: f32,
    denominator: f32,
    gain: f32,
    fast: bool,
    lookup: &'static [f32],
}
impl Brightness {
    pub fn new(exposure: f32, fast: bool) -> Self {
        let lookup = if fast { table() } else { &[] };
        let denominator = (10.0 * exposure).ln_1p();
        Self {
            exposure,
            denominator,
            gain: denominator.powf(-0.6),
            fast,
            lookup,
        }
    }
    pub fn value(&self, total: f32) -> f32 {
        if total >= self.exposure {
            return 1.0;
        }
        if self.fast && total.to_bits() >= MIN_BITS {
            let offset = total.to_bits() - MIN_BITS;
            let index = (offset >> SHIFT) as usize;
            if index + 1 < TABLE_LEN {
                let fraction = (offset & ((1 << SHIFT) - 1)) as f32 / (1 << SHIFT) as f32;
                let values = self.lookup;
                return ((values[index] + (values[index + 1] - values[index]) * fraction)
                    * self.gain)
                    .min(1.0);
            }
        }
        ((10.0 * total).ln_1p() / self.denominator)
            .min(1.0)
            .powf(0.6)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn interpolation_stays_below_one_rgb_quantization_step() {
        for exposure in [1.0, 10.0, 1000.0, 100000.0, 2.0_f32.powi(31)] {
            let exact = Brightness::new(exposure, false);
            let fast = Brightness::new(exposure, true);
            for i in 0..200000 {
                let total = 2.0_f32.powf(-30.0 + i as f32 * 62.0 / 199999.0);
                assert!(
                    (exact.value(total) - fast.value(total)).abs() * 255.0 < 0.01,
                    "exposure={exposure} total={total}"
                );
            }
            assert_eq!(fast.value(0.0), 0.0);
        }
    }
}
