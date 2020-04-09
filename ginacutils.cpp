/* Multiply a submatrix by the LCM of all of its denominators.
 * Divide it by the GCD of all of its numerators.
 * Apply normal() on each element.
 */
void
rescale_submatrix(matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc)
{
    vector<ex> n(nr*nc), d(nr*nc);
    ex mul = 0;
    ex div = 0;
    auto it_n = n.begin();
    auto it_d = d.begin();
    for (unsigned i = r; i < r + nr; i++) {
        for (unsigned j = c; j < c + nc; j++) {
            ex nd = m(i, j).normal().numer_denom();
            ex numer = nd.op(0);
            ex denom = nd.op(1);
            *it_n++ = numer;
            *it_d++ = denom;
            div = div.is_zero() ? numer : gcd(div, numer);
            //logd("lcm: {} {}", mul, denom);
            mul = mul.is_zero() ? denom : lcm(mul, denom);
        }
    }
    if (div.is_zero()) return;
    if (mul.is_zero()) return;
    // It would be tempting to exit here if div=mul=1, but we'd
    // then discard the normal() call results.
    it_n = n.begin();
    it_d = d.begin();
    for (unsigned i = r; i < r + nr; i++) {
        for (unsigned j = c; j < c + nc; j++) {
            ex nn, dd;
            bool ok1 = divide(*it_n++, div, nn);
            bool ok2 = divide(mul, *it_d++, dd);
            assert(ok1 && ok2);
            m(i, j) = nn*dd;
        }
    }
}

class matrix_hack : public matrix {
    public:
    void append_rows(const matrix &src);
    void append_row(const exvector &src);
    void resize(unsigned nrows);
    exvector &mvec();
};

void
matrix_hack::append_rows(const matrix &src)
{
    assert(col == src.cols());
    row += src.rows();
    m.insert(m.end(), ((const matrix_hack*)&src)->m.begin(), ((const matrix_hack*)&src)->m.end());
}

void
matrix_hack::append_row(const exvector &src)
{
    assert(col == src.size());
    row += 1;
    m.insert(m.end(), src.begin(), src.end());
}

void
matrix_hack::resize(unsigned nrows)
{
    row = nrows;
    m.resize(row*col);
}

exvector &
matrix_hack::mvec()
{
    ensure_if_modifiable();
    return m;
}

/* Transform a given matrix into upper-echelon form via Gauss
 * elimination.
 *
 * Requires O(N^3) GCD operations and about the same amount of
 * arithmetic ones for dense matrices, but about O(N*K) for
 * sparse matrices with K entries.
 */
void
echelon_form_gauss(matrix &m)
{
    unsigned nr = m.rows();
    unsigned nc = m.cols();
    exvector &mv = ((matrix_hack*)&m)->mvec();
    for (unsigned i = 0; i < nr*nc; i++) {
        mv[i] = mv[i].normal();
    }
    unsigned r0 = 0;
    unsigned c0 = 0;
    for (; (c0 < nc) && (r0 < nr - 1); c0++) {
        for (unsigned r = r0; r < nr; r++) {
            mv[r*nc + c0] = mv[r*nc + c0].normal();
        }
        // No normalization before is_zero() here, because
        // we maintain the matrix normalized throughout the
        // algorithm.
        unsigned pivot = r0;
        while ((pivot < nr) && mv[pivot*nc + c0].is_zero()) {
            pivot++;
        }
        if (pivot == nr) {
            // The whole column below r0:c0 is zero, let's skip
            // it.
            continue;
        }
        if (pivot > r0) {
            // Found a non-zero row somewhere below r0; let's
            // swap it in.
            for (unsigned c = c0; c < nc; c++) {
                mv[pivot*nc + c].swap(mv[r0*nc + c]);
            }
        }
        ex a = mv[r0*nc + c0];
        for (unsigned r = r0 + 1; r < nr; r++) {
            ex b = mv[r*nc + c0];
            if (!b.is_zero()) {
                ex k = b/a;
                mv[r*nc + c0] = 0;
                for (unsigned c = c0 + 1; c < nc; c++) {
                    mv[r*nc + c] = normal(mv[r*nc + c] - k*mv[r0*nc + c]);
                }
            }
        }
        r0++;
    }
    // Zero out the remaining rows (just in case).
    for (unsigned r = r0 + 1; r < nr; r++) {
        for (unsigned c = 0; c < nc; c++) {
            mv[r*nc + c] = 0;
        }
    }
}

/* A vector (sub-)space represented by a set of basis vectors.
 *
 * The basis vectors are stored as row vectors, but can be viewed
 * as either row or column vectors; hence the *_row and *_col
 * set of functions.
 */
struct vspace {
    matrix basis;
    vspace(unsigned n);
    vspace(const matrix &b);
    unsigned dim() const;
    unsigned length() const;
    matrix basis_col(unsigned i) const;
    matrix basis_row(unsigned i) const;
    const matrix basis_cols() const;
    const matrix &basis_rows() const;
    bool contains(const matrix &v) const;
    void add_rows(const matrix &v);
    void add_row(const exvector &v);
    void normalize();
};

vspace::vspace(unsigned n)
    : basis(0, n)
{
}

vspace::vspace(const matrix &b)
    : basis(b)
{
    for (unsigned i = 0; i < basis.rows(); i++) {
        rescale_submatrix(basis, i, 1, 0, basis.cols());
    }
    normalize();
}

void
vspace::add_rows(const matrix &v)
{
    ((matrix_hack*)&basis)->append_rows(v);
}

void
vspace::add_row(const exvector &v)
{
    ((matrix_hack*)&basis)->append_row(v);
}

void
vspace::normalize()
{
    echelon_form_gauss(basis);
    unsigned nrows = basis.rows();
    for (; nrows > 0; nrows--) {
        for (unsigned c = 0; c < basis.cols(); c++) {
            if (!basis(nrows - 1, c).normal().is_zero()) goto done;
        }
    }
done:;
    ((matrix_hack*)&basis)->resize(nrows);
    for (unsigned i = 0; i < basis.rows(); i++) {
        //rescale_submatrix(basis, i, 1, 0, basis.cols());
    }
}

unsigned
vspace::dim() const
{
    return basis.rows();
}

unsigned
vspace::length() const
{
    return basis.cols();
}

matrix
vspace::basis_col(unsigned i) const
{
    assert(i < basis.rows());
    matrix v(basis.cols(), 1);
    for (unsigned j = 0; j < basis.cols(); j++) {
        v.let_op(j) = basis(i, j);
    }
    return v;
}

matrix
vspace::basis_row(unsigned i) const
{
    assert(i < basis.rows());
    matrix v(1, basis.cols());
    for (unsigned j = 0; j < basis.cols(); j++) {
        v.let_op(j) = basis(i, j);
    }
    return v;
}

const matrix &
vspace::basis_rows() const
{
    return basis;
}

const matrix
vspace::basis_cols() const
{
    return basis.transpose();
}

bool
vspace::contains(const matrix &v) const
{
    //logd("BASIS:\n{}", basis);
    assert(v.nops() == basis.cols());
    matrix vv = v;
    rescale_submatrix(vv, 0, v.rows(), 0, v.cols());
    unsigned p = 0;
    // Division-free subtraction of basis vectors from v.
    for (unsigned i = 0; i < basis.rows(); i++, p++) {
        // Advance p to the first non-zero column of basis[i].
        for (;;) {
            // This assertion should only fail if the normalize()
            // was not called between add_rows() and contains().
            assert(p < basis.cols());
            if (!basis(i, p).is_zero()) break;
            vv.let_op(p) = normal(vv.op(p));
            // If vv has non-zero columns before p, it's not in
            // the basis.
            if (!vv.op(p).is_zero())
                return false;
            p++;
        }
        // Subtract basis[i] from vv, if vv[p] != 0.
        const ex &vv_p = vv.op(p);
        if (!vv_p.is_zero()) {
            const ex b_ip = basis(i, p);
            vv.let_op(p) = 0;
            for (unsigned j = p + 1; j < basis.cols(); j++) {
                vv.let_op(j) = normal(vv.op(j)*b_ip - basis(i, j)*vv_p);
            }
        }
    }
    for (unsigned i = p; i < basis.cols(); i++) {
        vv.let_op(i) = normal(vv.op(i));
        if (!vv.op(i).is_zero())
            return false;
    }
    return true;
}

/* Return a rectangular submatrix of a matrix.
 */
matrix
matrix_cut(const matrix &m, unsigned r, unsigned nr, unsigned c, unsigned nc)
{
    matrix res(nr, nc);
    for (unsigned i = 0; i < nr; i++) {
        for (unsigned j = 0; j < nc; j++) {
            res(i, j) = m(r + i, c + j);
        }
    }
    return res;
}

/* Iterate through terms of e, call yield(t) for each one.
 */
template <typename F> void
term_iter(const ex &e, F yield)
{
    if (is_a<add>(e)) {
        for (const auto &t : e) {
            yield(t);
        }
    } else {
        yield(e);
    }
}

/* Iterate through factors of e, call yield(f, k) for each
 * factor of the form f^k.
 *
 * Note that this function doesn't factor e itself, only iterates
 * through the factors already explicitly present.
 */
template <typename F> void
factor_iter(const ex &e, F yield)
{
    if (is_a<mul>(e)) {
        for (const auto &f : e) {
            if (is_a<power>(f)) {
                yield(f.op(0), ex_to<numeric>(f.op(1)).to_int());
            } else {
                yield(f, 1);
            }
        }
    } else {
        if (is_a<power>(e)) {
            yield(e.op(0), ex_to<numeric>(e.op(1)).to_int());
        } else {
            yield(e, 1);
        }
    }
}

/* Return a minor matrix obtained by crossing out a given row
 * and column.
 */
matrix
minor(const matrix &m, int row, int col)
{
    matrix minor(m.rows() - 1, m.cols() - 1);
    for (int i = 0; i < m.rows()-1; i++)
    for (int j = 0; j < m.cols()-1; j++) {
        minor(i, j) = m(i < row ? i : i+1, j < col ? j : j + 1);
    }
    return minor;
}

/* Return an adjugate matrix, det(A)*A^(-1).
 */
matrix
adjugate(const matrix &m)
{
    matrix cofactor(m.rows(), m.cols());
    matrix minor(m.rows() - 1, m.cols() - 1);
    for (int i = 0; i < m.rows(); i++)
    for (int j = 0; j < m.cols(); j++) {
        for (int ii = 0; ii < m.rows()-1; ii++)
        for (int jj = 0; jj < m.cols()-1; jj++) {
            minor(ii, jj) = m(ii < i ? ii : ii+1, jj < j ? jj : jj + 1);
        }
        cofactor(i, j) = (i + j) % 2 == 0 ?
            minor.determinant() : -minor.determinant();
    }
    return cofactor.transpose();
}

/* Return a given cell of an adjugate matrix.
 */
ex
adjugate(const matrix &m, int row, int col)
{
    return (row + col) % 2 == 0 ?
        minor(m, col, row).determinant() :
        -minor(m, col, row).determinant();
}

