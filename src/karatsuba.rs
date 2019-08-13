trait MulNaive<'a, Var, C: Coeff>{
    fn mul_naive(self, b: SeriesSlice<'a, Var, C>) -> Series<Var, C>;
}

const MIN_KARATSUBA_SIZE: usize = 8;

// dubious helpers trait that only serve to prevent obscure
// compiler errors (rust 1.36.0)
trait MulHelper<Var, C: Coeff> {
    fn add_poly_prod<'a, 'b>(
        &mut self,
        a: SeriesSlice<'a, Var, C>,
        b: SeriesSlice<'b, Var, C>,
    );
    fn add_poly_prod_naive<'a, 'b>(
        &mut self,
        a: SeriesSlice<'a, Var, C>,
        b: SeriesSlice<'b, Var, C>,
    );
    fn add_poly_prod_karatsuba<'a, 'b>(
        &mut self,
        a: SeriesSlice<'a, Var, C>,
        b: SeriesSlice<'b, Var, C>,
    );
}

impl<Var: Clone, C: Coeff> MulHelper<Var, C> for Series<Var, C>
where
    for<'d,'e> &'d C: Mul<&'e C, Output = C> + Add<&'e C, Output = C>,
    C: AddAssign + Clone,
    for <'d> Series<Var, C>: AddAssign<SeriesSlice<'d, Var, C>>
    + SubAssign<SeriesSlice<'d, Var, C>>
    + SubAssign<Series<Var, C>>,
{
    fn add_poly_prod<'a, 'b>(
        &mut self,
        a: SeriesSlice<'a, Var, C>,
        b: SeriesSlice<'b, Var, C>,
    ) {
        debug_assert_eq!(a.min_pow(), 0);
        debug_assert_eq!(b.min_pow(), 0);
        if self.coeffs.len() < MIN_KARATSUBA_SIZE {
            self.add_poly_prod_naive(a, b)
        }
        else {
            self.add_poly_prod_karatsuba(a, b)
        }
    }

    fn add_poly_prod_naive<'a, 'b>(
        &mut self,
        a: SeriesSlice<'a, Var, C>,
        b: SeriesSlice<'b, Var, C>,
    ) {
        debug_assert_eq!(self.min_pow(), 0);
        for (i, a) in a.iter() {
            for (j, b) in b.iter() {
                let idx = (i+j) as usize;
                self.coeffs[idx] += a*b;
            }
        }
    }

    fn add_poly_prod_karatsuba<'a, 'b>(
        &mut self,
        a: SeriesSlice<'a, Var, C>,
        b: SeriesSlice<'b, Var, C>,
    ) {
        let mid = ((min(a.coeffs.len(), b.coeffs.len()) + 1)/2) as isize;
        let (a_low, mut a_high) = a.split_at(mid);
        a_high.min_pow = 0;
        let (b_low, mut b_high) = b.split_at(mid);
        b_high.min_pow = 0;
        let mut t = Series::with_cutoff(
            self.var.clone(), 0, a_high.cutoff_pow(), vec![]
        );
        t += a_low;
        t += a_high;
        let mut u = Series::with_cutoff(
            self.var.clone(), 0, a_high.cutoff_pow(), vec![]
        );
        u += b_low;
        u += b_high;
        u.min_pow = mid;
        self.add_poly_prod(t.as_slice(), u.as_slice());
        let mut t = t;
        t.coeffs.clear();
        for _ in 0..2*mid {
            t.coeffs.push(C::from(0));
        }
        t.add_poly_prod(a_low, b_low);
        *self += t.as_slice();
        t.min_pow = mid;
        *self -= t.as_slice();
        for c in &mut t.coeffs {
            *c = C::from(0);
        }
        t.add_poly_prod(a_high, b_high);
        t.min_pow = 2*mid;
        *self += t.as_slice();
        t.min_pow = mid;
        *self -= t
    }
}
