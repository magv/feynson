static const char usagetext[] = R"(
Ss{NAME}
    Nm{feynson} -- a tool for Feynman integral symmetries.

Ss{SYNOPSYS}
    Nm{feynson} [Fl{options}] Cm{command} Ar{args} ...

Ss{DESCRIPTION}

Ss{COMMANDS}
    Cm{symmetrize} Ar{spec-file}
        Print a list of momenta substitutions that make symmetries
        between a list of integral families explicit.

        These momenta substitutions will make the set of
        denominators of two integral families exactly match (up
        to a reordering) if the families are isomorphic, and
        will make one a subset of the other if one family is
        isomorphic to a subsector of another family.

        The input specification file should be a list of three
        elements:
        1) a list of all integral families, with each family
           being a list of propagators (e.g. Ql{(l1+l2)^2});
        2) a list of all loop momenta;
        3) a list of external invariant substitution rules, each
           rule being a list of two elements: a scalar product
           and its substitution (e.g. Ql[{q^2, 1}] or Ql[{p1*p2, s12}]).

        Each family that can be mapped to (a subsector of) another
        is guaranteed to be mapped to the first possible family,
        in the order of specification.

        The families must come in the order of decreasing number
        of denominators.

    Cm{zero-sectors} [Fl{-s}] Ar{spec-file}
        Print a list of all zero sectors of a given integral
        family.

        The input specification file should be a list of four
        elements:
        1) a list of all propagator momenta (e.g. Ql{(l1-q)^2});
        2) a list of cut flags, Ql{0} for normal propagators, Ql{1}
           for cut propagators;
        3) a list of all loop momenta (e.g. Ql{l1});
        4) and a list of external invariant substitutions (e.g.
           Ql[{q^2, 1}]).

        For example: Ql[{ {(q-l)^2, l^2}, {0, 0}, {l}, {{q^2,1}} }].

        The output will be a list of zero sectors, each denoted
        by an integer s=2^{i_1-1} + ... + 2^{i_n-1}, where i_k
        are the indices of denominators that belong to this
        sector (counting from 1).

        If the Fl{-s} flag is given, the output will be shortened
        by only listing the topmost zero sectors: all the remaining
        zero sectors are their subsectors.

        Every sector that is missing a cut propagator of its
        supersectors will be reported as zero.

    Cm{ufx} Ar{spec-file}
        Print Feynman parametrization (U, F, X) of an integral
        defined by a set of propagators.

        The input specification file should be a list of three
        elements:
        1) a list of all propagators, e.g. Ql{(l1-q)^2};
        2) a list of all loop momenta, e.g. Ql{l1};
        3) and a list of external invariant substitutions, e.g.
           Ql[{q^2, 1}].

        For example: Ql[{ {(q-l)^2, l^2}, {l}, {{q^2,1}} }].

        The output will be a list of three items: the U polynomial,
        the F polynomial, and the list of Feynman parameter
        variables.

Ss{OPTIONS}
    Fl{-j} Ar{jobs}    Parallelize calculations using at most this many workers.
    Fl{-s}         Shorten the output (depending on the command).
    Fl{-q}         Print a more quiet log.
    Fl{-h}         Show this help message.
    Fl{-V}         Print version information.

Ss{ARGUMENTS}
    Ar{spec-file}  Filename of the input file, with Ql{-} meaning the standard input.

Ss{ENVIRONMENT}
    Ev{TMPDIR}     Temporary files will be created here.

Ss{AUTHORS}
    Vitaly Magerya <vitaly.magerya@tx97.net>
)";

#include <cln/real.h>
#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <assert.h>
#include <atomic>
#include <chrono>
#include <fstream>
#include <queue>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <tuple>
#include <unistd.h>

// The fact that a random macro defined in sys/sysmacros.h somehow
// finds its way here to mess up my day is one of the reasons
// people should stay away from C++ altogether. "A minor device
// number"?
#undef minor
#undef major

using namespace GiNaC;
using namespace std;

// Nauty defines the "set" type, and we don't need this sort of
// a conflict.
#define set Nauty_set
#include <nauty/nausparse.h>
#undef set

#include "blake2b.c"
#include "ginacutils.cpp"

static bool COLORS = !!isatty(STDOUT_FILENO);
static bool VERBOSE = true;
static int JOBS = 1;
static int WORKER = -1;

typedef vector<symbol> symvector;
struct hash_t { uint8_t hash[32]; };

static inline bool
operator ==(const hash_t &a, const hash_t &b)
{
    return memcmp((void*)&a.hash[0], (void*)&b.hash[0], sizeof(a.hash)) == 0;
}

static inline bool
operator <(const hash_t &a, const hash_t &b)
{
    return memcmp((void*)&a.hash[0], (void*)&b.hash[0], sizeof(a.hash)) < 0;
}

/* FORMATTING
 *
 * This needs to go before logging, so that log_format would be
 * able to see these definitions.
 */

ostream&
operator <<(ostream &o, const set<unsigned> &v)
{
    o << "set{";
    bool first = true;
    for (auto &&i : v) {
        if (first) {
            o << i;
            first = false;
        } else {
            o << ", " << i;
        }
    }
    o << "}";
    return o;
}

ostream&
operator <<(ostream &o, const hash_t &v)
{
    for (size_t i = 0; i < sizeof(v.hash); i++) {
        uint8_t dlo = v.hash[i] >> 4;
        uint8_t dhi = v.hash[i] & 0x0f;
        o << (char)(dhi < 10 ? '0' + dhi : 'a' - 10 + dhi);
        o << (char)(dlo < 10 ? '0' + dlo : 'a' - 10 + dlo);
    }
    return o;
}

template<typename T> ostream&
operator <<(ostream &o, const vector<T> &v)
{
    o << "{";
    for (size_t i = 0; i < v.size(); i++) {
        if (i != 0) o << ", ";
        o << v[i];
    }
    o << "}";
    return o;
}

ostream&
operator <<(ostream &o, const vector<uint8_t> &v)
{
    o << "{";
    for (size_t i = 0; i < v.size(); i++) {
        if (i != 0) o << ", ";
        o << (int)v[i];
    }
    o << "}";
    return o;
}

string
to_string(const ex& expr)
{
    ostringstream s;
    s << expr;
    return s.str();
}

/* LOGGING
 * ============================================================
 *
 * As with everything in C++, you can "optimize" this piece of
 * code to e.g. use compile-only format string parsing, minimize
 * number of created functions, and so on, but at the expense
 * of kilolines of code, and your own sanity lost in the fight
 * versus byzantine template rules. Please don't.
 */

static auto _log_starttime = chrono::steady_clock::now();
static auto _log_lasttime = chrono::steady_clock::now();
static int _log_depth = 0;

const char*
log_adv(const char *fmt)
{
    for (int i = 0; ; i++) {
        if (fmt[i] == '{') {
            cerr.write(fmt, i);
            return fmt + i + 2;
        }
        if (fmt[i] == 0) {
            cerr.write(fmt, i);
            return fmt + i;
        }
    }
}

/* This function is used to print objects into the log. Override it
 * for the data types you care about to modify their appearance.
 */
template<typename T> static inline void
log_format(ostream &o, const T &value)
{
    o << value;
}

void
log_print_start(const char *pre, const char *post)
{
    auto t = chrono::steady_clock::now();
    auto dt = chrono::duration_cast<chrono::duration<double>>(t - _log_starttime).count();
    cerr << pre;
    if ((JOBS > 1) && (WORKER >= 0)) cerr << "w" << WORKER + 1 << " ";
    cerr << std::fixed << std::setprecision(4) << dt << "s +";
    cerr << chrono::duration_cast<chrono::duration<double>>(t - _log_lasttime).count() << "s";
    for (int i = 0; i < _log_depth; i++) {
        cerr << " *";
    }
    cerr << post;
    _log_lasttime = t;
}

template<typename T> const char *
log_print_one(const char *fmt, const T &value)
{
    fmt = log_adv(fmt);
    log_format(cerr, value);
    return fmt;
}

void
log_print_end(const char *fmt)
{
    cerr << fmt;
    if (COLORS) cerr << "\033[0m";
    cerr << endl;
}

struct _sequencehack {
    template<typename ...Args>
    _sequencehack(Args &&...) {}
};

template<typename ...Args> static void
log_fmt(const char *pre, const char *post, const char *fmt, const Args &...args)
{
    log_print_start(pre, post);
    (void) _sequencehack {
        (fmt = log_print_one(fmt, args), 0)
        ...
    };
    log_print_end(fmt);
}

/* Log an debug message. These can be suppressed by setting
 * VERBOSE to false.
 */
template<typename... Args> static inline void
logd(const char *fmt, const Args &...args)
{
    if (VERBOSE) {
        log_fmt(COLORS ? "\033[2;37m[dbg " : "[dbg ", "] ", fmt, args...);
    }
}

/* Log an information message.
 */
template<typename... Args> static inline void
logi(const char *fmt, const Args &...args)
{
    log_fmt(COLORS ? "\033[32m[inf " : "[inf ", COLORS ? "]\033[0m " : "] ", fmt, args...);
}

/* Log a warning message.
 */
template<typename... Args> static inline void
logw(const char *fmt, const Args &...args)
{
    log_fmt(COLORS ? "\033[1;34m[wrn " : "[wrn ", "] ", fmt, args...);
}

/* Log an error message.
 */
template<typename... Args> static inline void
loge(const char *fmt, const Args &...args)
{
    log_fmt(COLORS ? "\033[1;31m[err " : "[err ", "] ", fmt, args...);
}

template<typename F>
struct _scopeexithack {
    _scopeexithack(F f) : f(f) {}
    ~_scopeexithack() { f(); }
    F f;
};

/* Place this macro as the start of a function, and you'll get
 * log entries every time this function is entered and exited.
 *
 * As a general rule, if you're adding logging statements to a
 * function, add LOGME to its start as well.
 */
#define LOGME \
    const auto &__log_func = __func__; \
    logd("> {}()", __log_func); \
    const auto __log_t0 = _log_lasttime; \
    _log_depth++; \
    auto __log_f = [&]{ \
        _log_depth--; \
        if (VERBOSE) { \
            auto t = chrono::steady_clock::now(); \
            auto dt = chrono::duration_cast<chrono::duration<double>>(t - __log_t0).count(); \
            logd("< {}(+{}s)",__log_func, dt); \
        } \
    }; \
    auto __log_s = _scopeexithack<decltype(__log_f)>(__log_f);

/* CONCURRENCY UTILS
 * =================
 */

/* Allocate size bytes of memory identically visible across the
 * forked workers.
 */
void *
shared_alloc(size_t size)
{
    // We could have used calloc() if FORK is false. Meh.
    void *mem = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    if (mem == MAP_FAILED) {
        loge("Failed to mmap() {} bytes of shared memory: {}", size, strerror(errno));
        exit(1);
    }
    return mem;
}

void
shared_free(void *memory, size_t size)
{
    if (munmap(memory, size) != 0) {
        logw("Failed to munmap() {} bytes of shared memory: {}", size, strerror(errno));
    }
}

/* FORK-JOIN MULTIPROCESSING
 * =========================
 *
 * Usage:
 *   FORK_BEGIN
 *      Worker code.
 *      Use WORKER variable (an int between 0 and JOBS-1) to
 *      determine which worker are you in.
 *      If FORK is false (or, equivalently, JOBS is 1), then no
 *      fork is done, and the worker code is executed in the
 *      main process.
 *   FORK_END
 *   Main process code (all the workes have joined at this point).
 *
 * Note: we must use fork-join instead of threads or OpenMP
 * because GiNaC is not thread-safe.
 */

#define FORK (JOBS > 1)

#define FORK_BEGIN \
    { \
        vector<pid_t> workerpids(JOBS); \
        if (FORK) logd("Forking {} workers", JOBS); \
        for (WORKER = 0; WORKER < JOBS; WORKER++) { \
            if (!FORK || ((workerpids[WORKER] = fork()) == 0)) { \
                {

#define FORK_END \
                } \
                if (FORK) exit(0); \
            } \
        } \
        WORKER = -1; \
        if (FORK) { \
            logd("Waiting for workers"); \
            int wst; \
            for (pid_t wpid; (wpid = wait(&wst)) > 0;) { \
                if (WIFEXITED(wst) && (WEXITSTATUS(wst) == 0)) { \
                    logd("Worker joined: {}", wpid); \
                } else { \
                    loge("Worker failed: {}, quitting", wpid); \
                    for (int i = 0; i < JOBS; i++) { \
                        kill(workerpids[i], SIGTERM); \
                    }; \
                    exit(1); \
                } \
            } \
        } \
    }

/* TEMPORARY FILES
 * ===============
 */

namespace tmpdir {
    static char dirname[2048];

    static void
    create(const char *prefix)
    {
        const char *tmp = getenv("TMPDIR");
        snprintf(dirname, sizeof(dirname), "%s/%s.%ld.XXXXXX",
                tmp == NULL ? "/tmp" : tmp, prefix, (long)getpid());
        if (mkdtemp(dirname) != dirname) {
            loge("Failed to create temporary directory: {}", strerror(errno));
            exit(1);
        }
        logd("Created temporary directory: {}", dirname);
    }

    static void
    remove()
    {
        logd("Removing temporary directory: {}", dirname);
        // There isn't a short way to run a command with a list
        // of arguments without worying about quoting. Fork,
        // exec, and then wait is the shortest...
        if (fork() == 0) {
            execl("/bin/rm", "/bin/rm", "-rf", dirname, NULL);
            exit(1);
        };
        wait(NULL);
    }

    static string
    filename(int suffix)
    {
        return string(dirname) + "/" + to_string(suffix);
    }
}

/* BIT UTILS
 * =========
 */

static unsigned
bitlength(unsigned x)
{
    unsigned n = 0;
    for (; x; x >>= 1) n++;
    return n;
}

static unsigned
bitcount(uint64_t x)
{
    unsigned n = 0;
    for (; x; x >>= 1) n += x&1;
    return n;
}

static int
bitposition(unsigned x, int bitn)
{
    int i = 0;
    for (int n = bitn + 1; x; x >>= 1, i++) {
        n -= x&1;
        if (n == 0) return i;
    }
    return -1;
}

static inline bool
is_subsector(uint64_t subsector, uint64_t sector)
{
    return (sector | subsector) == sector;
}

/* MISC UTILS
 * ==========
 */

ex
readfile(const char *filename, parser &reader)
{
    if (strcmp(filename, "-") == 0) {
        return reader(cin);
    } else {
        ifstream i(filename);
        if (!i) throw parse_error("the file was not found");
        return reader(i);
    }
}

symvector
symbolsequence(const char *prefix, int n)
{
    symvector x(n);
    for (int i = 0; i < n; i++) {
        ostringstream name;
        name << prefix << i + 1;
        x[i] = symbol(name.str());
    }
    return x;
}

ex
normalize_productrules(ex productrules)
{
    lst rules;
    for (unsigned i = 0; i < productrules.nops(); i++) {
        const ex &r = productrules.op(i);
        rules.append(r.op(0) == r.op(1));
    }
    return rules;
}

symvector
feynman_x(int ndens)
{
    return symbolsequence("x", ndens);
}

pair<ex, ex>
feynman_uf(const ex &denominators, const ex &loopmom, const ex &sprules, const symvector &X)
{
    LOGME;
    ex dens = denominators.expand().subs(sprules, subs_options::algebraic);
    unsigned N = dens.nops();
    unsigned L = loopmom.nops();
    matrix A(L, L);
    matrix B(L, 1);
    ex C;
    assert(denominators.nops() == X.size());
    for (unsigned d = 0; d < N; d++) {
        // Decompose D_i = a_ij l_i l_j + b_i l_i + c
        ex D = dens.op(d).expand();
        matrix a(L, L);
        matrix b(L, 1);
        bool has_a = false;
        bool has_b = false;
        for (unsigned i = 0; i < L; i++) {
            // l_i^2
            a(i, i) = D.coeff(loopmom.op(i), 2);
            has_a = has_a || (!a(i, i).is_zero());
            D -= a(i, i)*loopmom.op(i)*loopmom.op(i);
            // l_i^1
            ex ki = D.coeff(loopmom.op(i), 1);
            for (unsigned j = i + 1; j < L; j++) {
                ex k = ki.coeff(loopmom.op(j), 1);
                a(i, j) = a(j, i) = k/2;
                has_a = has_a || (!k.is_zero());
                D -= k*loopmom.op(i)*loopmom.op(j);
            }
            b(i, 0) = D.coeff(loopmom.op(i), 1);
            has_b = has_b || (!b(i, 0).is_zero());
            D -= (b(i, 0)*loopmom.op(i)).expand();
        }
        A = A.add(a.mul_scalar(X[d]));
        B = B.add(b.mul_scalar(X[d]/2));
        C += X[d]*D;
        if (!has_a) {
            if (has_b) {
                loge("feynman_uf(): denominator #{} (={}) is linear in the loop momenta", d+1, dens.op(d));
            } else {
                loge("feynman_uf(): denominator #{} (={}) is free of the loop momenta", d+1, dens.op(d));
            }
            exit(1);
        }
    }
    ex U = A.determinant().expand();
    ex F = C*U;
    for (unsigned i = 0; i < L; i++)
    for (unsigned j = 0; j < L; j++) {
        ex k = (B(i, 0)*B(j, 0)).expand().subs(sprules, subs_options::algebraic);
        if (!k.is_zero()) F -= adjugate(A, i, j)*k;
    }
    return make_pair(U, F.expand());
}

map<vector<int>, ex>
bracket(const ex &expr, const symvector &stemset)
{
    LOGME;
    map<vector<int>, ex> result;
    map<ex, int, ex_is_less> sym2id;
    for (unsigned i = 0; i < stemset.size(); i++) {
        sym2id[stemset[i]] = i;
    }
    term_iter(expr.expand(), [&](const ex &term) {
        vector<int> stemidx(stemset.size());
        exvector coef;
        factor_iter(term, [&](const ex &factor, int power) {
            auto it = sym2id.find(factor);
            if (it != sym2id.end()) {
                stemidx[it->second] = power;
            } else {
                if (power == 1) {
                    coef.push_back(factor);
                } else {
                    coef.push_back(pow(factor, power));
                }
            }
        });
        result[stemidx] += mul(coef);
    });
    return result;
}

bool
family_is_zero(const ex &G, const symvector &X)
{
    LOGME;
    symvector K = symbolsequence("k", X.size());
    ex eqn = G;
    for (unsigned i = 0; i < X.size(); i++) {
        eqn -= K[i]*X[i]*G.diff(X[i]);
    }
    matrix eqns(0, K.size() + 1);
    for (auto &&kv : bracket(eqn, X)) {
        uint64_t sector = 0;
        for (unsigned i = 0; i < kv.first.size(); i++) {
            if (kv.first[i] != 0) sector |= 1ul << i;
        }
        ex c0 = kv.second;
        exvector c(K.size() + 1);
        for (unsigned i = 0; i < K.size(); i++) {
            c[i] = c0.coeff(K[i]);
            c0 -= c[i]*K[i];
        }
        c[X.size()] = c0.expand();
        ((matrix_hack*)&eqns)->append_row(c);
    }
    return eqns.rank() <= K.size();
}

vector<bool>
zero_sectors(const ex &G, const symvector &X, uint64_t cutmask)
{
    LOGME;
    // Lee criterium: k_i x_i dG/dx_i == G
    symvector K = symbolsequence("k", X.size());
    ex eqn = G;
    for (unsigned i = 0; i < X.size(); i++) {
        eqn -= K[i]*X[i]*G.diff(X[i]);
    }
    // Split eqn into {x1^p1 ... xn^pn} * { c0 + k1*c1 + ... + kn*cn },
    // and save the c vectors.
    vector<pair<uint64_t, exvector>> eqns;
    for (auto &&kv : bracket(eqn, X)) {
        uint64_t sector = 0;
        for (unsigned i = 0; i < kv.first.size(); i++) {
            if (kv.first[i] != 0) sector |= 1ul << i;
        }
        ex c0 = kv.second;
        exvector c(K.size() + 1);
        for (unsigned i = 0; i < K.size(); i++) {
            c[i] = c0.coeff(K[i]);
            c0 -= c[i]*K[i];
        }
        c[X.size()] = c0.expand();
        eqns.push_back(make_pair(sector, c));
    }
    logd("Total equations: {}", eqns.size());
    // Solve the equation for each sector;
    uint64_t nsectors = 1UL << X.size();
    logd("Total sectors: {}", nsectors);
    uint8_t *zeros = (uint8_t*)shared_alloc(nsectors);
    uint8_t *nonzeros = (uint8_t*)shared_alloc(nsectors);
    zeros[0] = true;
    matrix seceqns(0, K.size() + 1);
    FORK_BEGIN;
    // Any order of traversal will work; bottom up seems to be
    // slightly faster.
    for (uint64_t sector = nsectors/JOBS*WORKER; sector < nsectors; sector++) {
        if (zeros[sector] || nonzeros[sector]) continue;
        if (~sector & cutmask) {
            // If sector has a 0 where cutmask has a 1.
            zeros[sector] = true;
            continue;
        }
        ((matrix_hack*)&seceqns)->resize(0);
        for (auto &&sc : eqns) {
            if (is_subsector(sc.first, sector)) {
                ((matrix_hack*)&seceqns)->append_row(sc.second);
            }
        }
        unsigned rank = seceqns.rank();
        unsigned ndens = bitcount(sector);
        if (rank <= ndens) {
            for (uint64_t s = sector; s != 0; s = (s - 1) & sector) {
                assert(is_subsector(s, sector));
                zeros[s] = true;
            }
        } else {
            for (uint64_t s = sector; s < nsectors; s = (s + 1) | sector) {
                assert(is_subsector(sector, s));
                nonzeros[s] = true;
            }
        }
    }
    FORK_END;
    vector<bool> result(zeros, zeros + nsectors);
    shared_free((void*)nonzeros, nsectors);
    shared_free((void*)zeros, nsectors);
    return result;
}

static vector<pair<vector<int>, int>>
subsector_bracket(const vector<pair<vector<int>, int>> &br, unsigned nx, unsigned mask)
{
    assert((br.size() == 0) || (br[0].first.size() == nx));
    int newnx = bitcount(mask);
    vector<pair<vector<int>, int>> result;
    for (auto &&stemcoef : br) {
        vector<int> stem(newnx);
        for (unsigned i = 0, j = 0; i < nx; i++) {
            if (mask & (1u << i)) {
                stem[j++] = stemcoef.first[i];
            } else {
                if (stemcoef.first[i] != 0) goto skip;
            }
        }
        result.push_back(make_pair(stem, stemcoef.second));
    skip:;
    }
    return result;
}

static vector<uint8_t>
canonical_variable_permutation(const vector<pair<vector<int>, int>> &polybr, unsigned nx)
{
    assert((polybr.size() == 0) || (polybr[0].first.size() == nx));
    // A list of stems and corresponding unique coefficient ids.
    set<int> coefset;
    for (auto &&stemcoef : polybr) {
        coefset.insert(stemcoef.second);
    }
    int maxexponent = 1;
    for (auto &&stemcoef : polybr) {
        for (int p : stemcoef.first) {
            maxexponent = max(maxexponent, p);
        }
    }
    // We shall map the polynomial to a graph, with
    // - a vertex x_i of color 0 for each x_i in X;
    // - a vertex stem_j of color 1+uniqid[coef_j] for each
    //   term stem_j*coef_j in the bracketed polynomial;
    // - an edge stem_k--x_i of color p for each x_i^p in stem_k.
    // This construction ensures that all permutations of X
    // that leave the polynomial unchanged correspond to graph
    // relabelings, and vice versa.
    // Note that edge colors are 1...n, because edge color 0
    // corresponds to "no edge".
    // Because nauty doesn't support edge colors directly, only
    // vertex colors, we shall adjust by:
    // - making multiple layers (copies) of vertices, each layer
    //   having a distinct set of vertex colors, with one line
    //   connecting all copies of the same vertex;
    // - replicating an edge of color c in layer l if c&(1<<l).
    int maxedgecolor = maxexponent;
    unsigned nlayers = bitlength(maxedgecolor);
    // Construct the graph.
    SG_DECL(g);
    g.nv = 0; g.nde = 0;
    unsigned nterms = polybr.size();
    SG_ALLOC(g, 1 + (nx + nterms)*nlayers, 32, "graph malloc");
#define BEGIN_VERTEX g.v[g.nv] = g.nde; g.d[g.nv] = 0;
#define ADD_EDGE(vertex) \
    while (g.nde > g.elen) {DYNREALLOC(int, g.e, g.elen, g.elen*2, "graph realloc");}; \
    g.e[g.nde++] = (vertex); \
    g.d[g.nv]++;
#define END_VERTEX g.nv++;
    // Vertices go in this order:
    //   x1(1) ... xN(1) x1(2) ... xN(L) ... x1(L) ... xN(L)
    //   BR1(1) ... BR1(L) ...
#define VERTEX_Xil(i,l) nx*(l) + (i)
#define VERTEX_Bkl(k,l) nx*nlayers + nlayers*(k) + (l)
    for (unsigned l = 0, lbit = 1; l < nlayers; l++, lbit <<= 1) {
        for (unsigned i = 0; i < nx; i++) {
            BEGIN_VERTEX;
            if (l > 0) { ADD_EDGE(VERTEX_Xil(i, l - 1)); }
            if (l < nlayers - 1) { ADD_EDGE(VERTEX_Xil(i, l + 1)); }
            unsigned k = 0;
            for (auto &&stemcoef : polybr) {
                if (stemcoef.first[i] & lbit) {
                    ADD_EDGE(VERTEX_Bkl(k, l));
                }
                k++;
            }
            END_VERTEX;
        }
    }
    unsigned k = 0;
    for (auto &&stemcoef : polybr) {
        for (unsigned l = 0, lbit = 1; l < nlayers; l++, lbit <<= 1) {
            BEGIN_VERTEX;
            if (l > 0) { ADD_EDGE(VERTEX_Bkl(k, l - 1)); }
            if (l < nlayers - 1) { ADD_EDGE(VERTEX_Bkl(k, l + 1)); }
            for (unsigned i = 0; i < nx; i++) {
                if (stemcoef.first[i] & lbit) {
                    ADD_EDGE(VERTEX_Xil(i, l));
                }
            }
            END_VERTEX;
        }
        k++;
    }
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    DYNALLOC1(int, lab, lab_sz, g.nv, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, g.nv, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, g.nv, "malloc");
    int partidx = 0;
#define BEGIN_PARTITION
#define ADD_VERTEX(vertex) lab[partidx] = (vertex); ptn[partidx] = 1; partidx++;
#define END_PARTITION ptn[partidx-1] = 0;
    /* The partition has this order:
     *   x1(1) ... xn(1) | x1(2) ... x1(2) | ...
     *   BR1_c1(1) ... BRn_c1(1) | BR1_c2(1) ... BRn_c2(1) | ...
     */
    for (unsigned l = 0; l < nlayers; l++) {
        BEGIN_PARTITION; // layer #l vertices
        for (unsigned i = 0; i < nx; i++) {
            ADD_VERTEX(VERTEX_Xil(i, l));
        }
        END_PARTITION;
    }
    for (unsigned l = 0; l < nlayers; l++) {
        // Walk coefficients in increasing order.
        for (int coef : coefset) {
            BEGIN_PARTITION; // layer #l, terms with coef #coef
            unsigned k = 0;
            for (auto &&stemcoef : polybr) {
                if (stemcoef.second == coef) {
                    ADD_VERTEX(VERTEX_Bkl(k, l));
                }
                k++;
            }
            END_PARTITION;
        }
    }
    assert(partidx == g.nv);
    SG_DECL(cg);
    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.defaultptn = FALSE;
    options.getcanon = TRUE;
    statsblk stats;
    sparsenauty(&g, lab, ptn, orbits, &options, &stats, &cg);
    SG_FREE(cg);
    SG_FREE(g);
    vector<uint8_t> result(nx);
    for (unsigned i = 0; i < nx; i++) {
        result[i] = lab[i];
    }
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
    return result;
}

ex
canonical_poly(const ex &poly, symvector X, vector<int> perm)
{
    lst sub;
    for (unsigned i = 0; i < X.size(); i++) {
        sub.append(X[perm[i]] == X[i]);
    }
    return poly.subs(sub);
}

hash_t
canonical_hash(const vector<pair<vector<int>, int>> &polybr, unsigned nx, vector<uint8_t> perm)
{
    assert((polybr.size() == 0) || (polybr[0].first.size() == nx));
    vector<vector<int>> expr;
    // { 1 0 2 }  ->  x1 * x3^2  /.  x1->x2 x2->x3 x3->x1  ->  x2 x1^2  ->  { 2 1 0 }
    for (unsigned i = 0; i < polybr.size(); i++) {
        vector<int> term(nx + 1);
        for (unsigned j = 0; j < nx; j++) {
            term[j] = polybr[i].first[perm[j]];
        }
        term[nx] = polybr[i].second;
        expr.push_back(term);
    }
    unsigned nterms = expr.size();
    sort(expr.begin(), expr.end());
    blake2b_state S;
    blake2b_init(&S, sizeof(hash_t));
    blake2b_update(&S, (void*)&nx, sizeof(nx));
    blake2b_update(&S, (void*)&nterms, sizeof(nterms));
    for (auto &&term : expr) {
        blake2b_update(&S, (void*)&term[0], sizeof(term[0])*term.size());
    }
    hash_t result;
    blake2b_final(&S, result.hash, sizeof(result.hash));
    return result;
}

vector<uint8_t>
canonical_variable_permutation(const ex &poly, const symvector &X)
{
    // A list of stems and corresponding unique coefficient ids.
    map<vector<int>, ex> br = bracket(poly.expand(), X);
    // Find unique coefficients and sort them.
    exset coefset;
    for (auto &&stemcoef : br) {
        coefset.insert(stemcoef.second);
    }
    // Assign each unique coefficient an index, in increasing
    // order.
    map<ex, int, ex_is_less> coef2uniqid;
    for (auto &&coef : coefset) {
        logd(" C{} = {}", coef2uniqid.size(), coef);
        coef2uniqid[coef] = coef2uniqid.size();
    }
    vector<pair<vector<int>, int>> polybr;
    for (auto &&stemcoef : br) {
        polybr.push_back(make_pair(stemcoef.first, coef2uniqid[stemcoef.second]));
    }
    return canonical_variable_permutation(polybr, X.size());
}

exvector
lincoefficients(const ex &expr, const symvector &vars)
{
    exvector coefs;
    coefs.reserve(vars.size() + 1);
    ex restofexpr = expr;
    for (unsigned i = 0; i < vars.size(); i++) {
        if (restofexpr.degree(vars[i]) > 1) {
            loge("lincoefficients: {} is not linear in {}", expr, vars[i]);
            exit(1);
        }
        ex c = restofexpr.coeff(vars[i]);
        coefs.push_back(c);
        restofexpr = restofexpr - c*vars[i];
    }
    coefs.push_back(restofexpr);
    return coefs;
}

ex
find_momenta_map(const exvector &src, const exvector &dst, const ex &loopmom)
{
    LOGME;
    assert(src.size() == dst.size());
    symvector newl = symbolsequence("$l", loopmom.nops());
    exmap newloopmommap;
    for (unsigned i = 0; i < loopmom.nops(); i++) {
        newloopmommap[loopmom.op(i)] = newl[i];
    }
    vector<vector<exvector>> eqnsets;
    for (unsigned i = 0; i < src.size(); i++) {
        ex eq = src[i].subs(newloopmommap) - dst[i];
        vector<exvector> eqns;
        factor_iter(factor(eq), [&](const ex &factor, int power) {
            (void)power;
            if (is_a<numeric>(factor)) return;
            eqns.push_back(lincoefficients(factor, newl));
        });
        sort(eqns.begin(), eqns.end());
        eqnsets.push_back(eqns);
    }
    sort(eqnsets.begin(), eqnsets.end());
    vector<unsigned> eqnum(eqnsets.size());
    unsigned eqset = 0;
    stack<vspace> mxstack;
    mxstack.push(vspace(loopmom.nops() + 1));
    while (eqset < eqnsets.size()) {
        auto mx = mxstack.top();
        // The first eqset equations are chosen and are compatible.
        // Choose the next one.
        mx.add_row(eqnsets[eqset][eqnum[eqset]]);
        mx.normalize();
        bool lastrowiszero = true;
        for (unsigned j = 0; j < mx.length() - 1; j++) {
            lastrowiszero = lastrowiszero && (mx.basis_rows()(mx.dim() - 1, j).is_zero());
        }
        if (lastrowiszero) {
            // The system is inconsistent, choose different
            // equation from this set.
            while (++eqnum[eqset] >= eqnsets[eqset].size()) {
                // No more equations at this level, roll back
                // to the previous set.
                eqnum[eqset] = 0;
                if (eqset == 0) {
                    // Nowhere to roll back. Fail.
                    loge("find_momenta_map(): no solutions at all; an external momenta symmetry?");
                    exit(1);
                }
                eqset--;
                mxstack.pop();
            }
        } else {
            // The system remains consistent, continue to the
            // next set.
            eqset++;
            mxstack.push(mx);
        }
    }
    const matrix &mx = mxstack.top().basis_rows();
    matrix M = matrix_cut(mx, 0, loopmom.nops(), 0, loopmom.nops());
    matrix C = matrix_cut(mx, 0, loopmom.nops(), loopmom.nops(), 1);
    matrix sol = M.inverse().mul(C).mul_scalar(-1);
    lst result;
    for (unsigned i = 0; i < loopmom.nops(); i++) {
        if (!loopmom.op(i).is_equal(sol(i, 0))) {
            result.append(lst{loopmom.op(i), sol(i, 0)});
        }
    }
    return result;
}

/* MAIN
 */

void
usage()
{
    const char *p = strchr(usagetext, '\n') + 1;
    for (;;) {
        const char *l1 = strchr(p + 2, '{');
        const char *l2 = strchr(p + 2, '[');
        if ((l1 == NULL) && (l2 == NULL)) break;
        const char *l = l1 == NULL ? l2 : l2 == NULL ? l1 : l1 < l2 ? l1 : l2;
        const char *r = strchr(l, (*l == '{') ? '}' : ']');
        if (r == NULL) break;
        const char *a = "", *b = "\033[0m";
        if (l[-2] == 'S' && l[-1] == 's') { a = "\033[1m"; goto found; }
        if (l[-2] == 'N' && l[-1] == 'm') { a = "\033[1;35m"; goto found; }
        if (l[-2] == 'F' && l[-1] == 'l') { a = "\033[33m"; goto found; }
        if (l[-2] == 'C' && l[-1] == 'm') { a = "\033[1m"; goto found; }
        if (l[-2] == 'A' && l[-1] == 'r') { a = "\033[32m"; goto found; }
        if (l[-2] == 'E' && l[-1] == 'v') { a = "\033[34m"; goto found; }
        if (l[-2] == 'Q' && l[-1] == 'l') { a = "\033[35m"; goto found; }
        cout.write(p, r + 1 -p);
        p = r + 1;
        continue;
found:
        cout.write(p, l - p - 2);
        if (COLORS) cout << a;
        cout.write(l + 1, r - l - 1);
        if (COLORS) cout << b;
        p = r + 1;
    }
    cout << p;
}

void
main_ufx(const char *specfile)
{
    parser reader;
    ex input = readfile(specfile, reader);
    ex denominators = input.op(0);
    ex loopmomenta = input.op(1);
    ex productrules = normalize_productrules(input.op(2));
    auto x = feynman_x(denominators.nops());
    auto uf = feynman_uf(denominators, loopmomenta, productrules, x);
    cout << "{\n " << uf.first << ",\n " << uf.second << ",\n " << x << "\n}" << endl;
}

void
main_zerosectors(const char *specfile, bool SHORT)
{
    parser reader;
    ex input = readfile(specfile, reader);
    ex denominators = input.op(0);
    ex cutflags = input.op(1);
    ex loopmomenta = input.op(2);
    ex productrules = normalize_productrules(input.op(3));
    auto x = feynman_x(denominators.nops());
    auto uf = feynman_uf(denominators, loopmomenta, productrules, x);
    uint64_t cutmask = 0;
    for (unsigned i = 0; i < cutflags.nops(); i++) {
        if (!cutflags.op(i).is_zero()) cutmask |= (1ul << i);
    }
    auto zeros = zero_sectors(uf.first + uf.second, x, cutmask);
    if (SHORT) {
        // Only keep the topmost-level sectors; hide all zero
        // sectors that are subsectors of other zero sectors.
        for (uint64_t sec = 0; sec < zeros.size(); sec++) {
            if (zeros[sec]) {
                for (uint64_t s = (sec - 1) & sec; s > 0; s = (s - 1) & sec) {
                    zeros[s] = false;
                }
                zeros[0] = false;
            }
        }
    }
    cout << "{";
    bool first = true;
    for (uint64_t sec = 0; sec < zeros.size(); sec++) {
        if (zeros[sec]) {
            if (first) {
                cout << "\n " << sec;
                first = false;
            } else {
                cout << ",\n " << sec;
            }
        }
    }
    cout << "\n}" << endl;
}

void
main_symmetrize(const char *specfile)
{
    parser reader;
    ex input = readfile(specfile, reader);
    ex families = input.op(0);
    ex loopmomenta = input.op(1);
    ex productrules = normalize_productrules(input.op(2));
    // Compute the set of all family sizes. We will only be
    // interested in sectors of these sizes.
    set<unsigned> familysizeset;
    unsigned maxfamilysize = 0;
    for (auto &&family : families) {
        unsigned ndens = family.nops();
        assert(ndens < 256);
        familysizeset.insert(ndens);
        maxfamilysize = max(maxfamilysize, ndens);
    }
    logd("Family size set: {}", familysizeset);
    // Compute G polynomials, their brackets, and the zero
    // sectors for each family.
    symvector x = feynman_x(maxfamilysize);
    map<ex, int, ex_is_less> coef2uniqid;
    map<int, vector<pair<vector<int>, int>>> gbrackets;
    if (FORK) {
        tmpdir::create("feynson");
    }
    FORK_BEGIN;
        for (unsigned fam = WORKER; fam < families.nops(); fam += JOBS) {
            logd("Preparing family {}", fam + 1);
            auto &&family  = families[fam];
            symvector xi(x.begin(), x.begin() + family.nops());
            auto uf = feynman_uf(family, loopmomenta, productrules, xi);
            auto g = uf.first + uf.second;
            auto br = bracket(g, xi);
            for (auto &&stemcoef : br) {
                auto &&it = coef2uniqid.find(stemcoef.second);
                if (it == coef2uniqid.end()) {
                    logd("Unique coefficient C{} = {}", coef2uniqid.size(), stemcoef.second);
                    coef2uniqid[stemcoef.second] = coef2uniqid.size();
                }
            }
            vector<pair<vector<int>, int>> gbr;
            for (auto &&stemcoef : br) {
                gbr.push_back(make_pair(stemcoef.first, coef2uniqid[stemcoef.second]));
            }
            gbrackets[fam] = gbr;
        }
        if (FORK) {
            logd("Exporting results");
            ofstream f(tmpdir::filename(WORKER));
            f << coef2uniqid.size() << "\n";
            for (auto &&coefid : coef2uniqid) {
                f << " " << coefid.second << " " << coefid.first << "\n";
            }
            f << gbrackets.size() << "\n";
            for (auto &&famgbr : gbrackets) {
                unsigned nx = families[famgbr.first].nops();
                f << " " << famgbr.first << " " << nx << " " << famgbr.second.size() << "\n";
                for (auto &&stemcoef : famgbr.second) {
                    f << "  " << stemcoef.second;
                    for (unsigned i = 0; i < nx; i++) {
                        f << " " << stemcoef.first[i];
                    }
                    f << "\n";
                }
            }
            logd("Export done");
        }
    FORK_END;
    if (FORK) {
        logd("Loading & combining worker outputs");
        for (int worker = 0; worker < JOBS; worker++) {
            ifstream f(tmpdir::filename(worker));
            int nc; f >> nc;
            map<int, int> cidmap;
            for (int i = 0; i < nc; i++) {
                int id; f >> id;
                string line; f >> ws; getline(f, line);
                ex coef = reader(line);
                auto it = coef2uniqid.find(coef);
                if (it == coef2uniqid.end()) {
                    logd("Unique coefficient C{} = W{}.C{} = {}",
                            coef2uniqid.size(), worker, id, coef);
                    coef2uniqid[coef] = cidmap[id] = coef2uniqid.size();
                } else {
                    logd("Known coefficient C{} = W{}.C{} = ", it->second, worker, id, coef);
                    cidmap[id] = it->second;
                }
            }
            int nfam; f >> nfam;
            for (int i = 0; i < nfam; i++) {
                int fam; f >> fam;
                int nx; f >> nx;
                int nterms; f >> nterms;
                vector<pair<vector<int>, int>> br;
                for (int t = 0; t < nterms; t++) {
                    int coefid; f >> coefid;
                    vector<int> stem(nx);
                    for (int x = 0; x < nx; x++) {
                        int p; f >> p;
                        stem[x] = p;
                    }
                    br.push_back(make_pair(stem, cidmap[coefid]));
                }
                gbrackets[fam] = br;
            }
        }
        logd("Done combining output");
        tmpdir::remove();
    }
    // Enumerate the sectors we're interested in.
    map<pair<unsigned, uint64_t>, uint64_t> sector2idx;
    for (unsigned fam = 0; fam < families.nops(); fam++) {
        uint64_t nsectors = 1ul << families[fam].nops();
        // TODO: these two loops waste time.
        for (unsigned ssize : familysizeset) {
            for (uint64_t sec = 1; sec < nsectors; sec++) {
                if ((bitcount(sec) == ssize)/* && !zeros[fam][sec]*/) {
                    sector2idx[make_pair(fam, sec)] = sector2idx.size();
                }
            }
        }
    }
    // Compute canonical permutations and hashes for top-level
    // sectors (those with largest number of propagators).
    uint64_t totalsectors = sector2idx.size();
    logd("Total interesting sectors: {}", sector2idx.size());
    hash_t *canonicalhashes = (hash_t*)shared_alloc(sizeof(hash_t) * totalsectors);
    uint8_t *canonicalperms = (uint8_t*)shared_alloc(maxfamilysize * sizeof(uint8_t) * totalsectors);
    uint8_t *hash_done = (uint8_t*)shared_alloc(sizeof(uint8_t) * totalsectors);
    uint8_t *fully_mapped = (uint8_t*)shared_alloc(sizeof(uint8_t) * totalsectors);
    atomic<int> *hashes_done = (atomic<int> *)shared_alloc(sizeof(atomic<int> *));
    FORK_BEGIN;
        logd("Precomputing canonical polynomials of each family");
        for (unsigned fam = WORKER; fam < families.nops(); fam += JOBS) {
            unsigned nx = families[fam].nops();
            uint64_t sector = (1ul << nx) - 1;
            uint64_t idx = sector2idx[make_pair(fam, sector)];
            assert(!hash_done[idx]);
            auto perm = canonical_variable_permutation(gbrackets[fam], nx);
            memcpy(canonicalperms + idx*maxfamilysize, &perm[0], nx);
            hash_t h = canonical_hash(gbrackets[fam], nx, perm);
            canonicalhashes[idx] = h;
            // "Release"-ordered fence. This prevents any preceding
            // reads or writes from being ordered past
            // any subsequent writes.
            // We need this so that any thread observing
            // hash_done to be 1 would also observe
            // canonicalhashes to be of correct value.
            // On x86 this is only an instruction for the
            // compiler not to reorder memory access across
            // this line; the hardware does not do this sort
            // of memory reordering.
            atomic_thread_fence(memory_order_release);
            hash_done[idx] = true;
            hashes_done->fetch_add(1, memory_order_relaxed);
            logd("Family {}, sector {} is {}", fam+1, sector, h);
        }
    FORK_END;
    // Wait for all top-level hashes to finish; then, for
    // each sector size walk through each sub-sector of
    // each family in order, compute and remember hashes of
    // their canonical polys, noting any duplicates along
    // the way.
    map<int, string> fam2mommap;
    if (FORK) {
        tmpdir::create("feynson");
    }
    FORK_BEGIN;
        // Pair each worker with its own set of family sizes.
        auto famsize = familysizeset.rbegin();
        // Skip WORKER sizes.
        for (int i = 0 ; (i < WORKER) && (famsize != familysizeset.rend()); i++, famsize++);
        for (; famsize != familysizeset.rend(); ) {
            logd("Computing symmetries for families with {} propagators", *famsize);
            map<hash_t, pair<int, uint64_t>> hash2sector;
            for (unsigned fam = 0; fam < families.nops(); fam += 1) {
                unsigned nx = families[fam].nops();
                if (nx != *famsize) continue;
                // Looking at family fam.
                uint64_t sector = (1ul << nx) - 1;
                uint64_t idx = sector2idx[make_pair(fam, sector)];
                assert(hash_done[idx]);
                if (1) {
                    // Check if the sector hash was already seen
                    // and put into canonicalhashes. This is an
                    // optional speedup.
                    auto secit = hash2sector.find(canonicalhashes[idx]);
                    if (secit != hash2sector.end()) {
                        int fam2 = secit->second.first;
                        uint64_t sec2 = secit->second.second;
                        uint64_t idx2 = sector2idx[secit->second];
                        fully_mapped[fam] = true;
                        logi("Family {} (top sector {}) is symmetric to family {}, sector {}",
                                fam+1, sector, fam2+1, sec2);
                        uint8_t *perm = canonicalperms + idx*maxfamilysize;
                        uint8_t *perm2 = canonicalperms + idx2*maxfamilysize;
                        logd("Fam {} is {}, sec {} perm: {}", fam+1, families[fam],
                                sector, vector<uint8_t>(perm, perm+nx));
                        logd("Fam {} is {}, sec {} perm: {}", fam2+1, families[fam2],
                                sec2, vector<uint8_t>(perm2, perm2+nx));
                        exvector src, dst;
                        for (unsigned i = 0; i < nx; i++) {
                            int p2i = bitposition(sec2, perm2[i]);
                            logd("  {} == {}", families[fam][perm[i]], families[fam2][p2i]);
                            src.push_back(families[fam][perm[i]]);
                            dst.push_back(families[fam2][p2i]);
                        }
                        auto mommap = find_momenta_map(src, dst, loopmomenta);
                        logd("Mom map: {}", mommap);
                        fam2mommap[fam] = to_string(mommap);
                        goto found;
                    }
                }
                for (auto &&famsecidx : sector2idx) {
                    unsigned fam2 = famsecidx.first.first;
                    int sec2 = famsecidx.first.second;
                    int idx2 = famsecidx.second;
                    if (fam2 >= fam) continue;
                    if (bitcount(sec2) != nx) continue;
                    if (fully_mapped[fam2]) continue;
                    if (!hash_done[idx2]) {
                        unsigned nx2 = families[fam2].nops();
                        auto br2 = subsector_bracket(gbrackets[fam2], nx2, sec2);
                        auto perm = canonical_variable_permutation(br2, nx);
                        memcpy(canonicalperms + idx2*maxfamilysize, &perm[0], nx);
                        hash_t h = canonical_hash(br2, nx, perm);
                        canonicalhashes[idx2] = h;
                        atomic_thread_fence(memory_order_release);
                        hash_done[idx2] = true;
                        hashes_done->fetch_add(1, memory_order_relaxed);
                        logd("Sector {}:{} is {}", fam2+1, sec2, h);
                    }
                    // "Acquire"-ordered memory fence. This prevents
                    // any following reads or writes to be ordered
                    // before any preceeding reads.
                    // We need this here so that the following reads from
                    // canonicalhashes would not be reordered before the
                    // read of hash_done above -- unlikely as it is.
                    atomic_thread_fence(memory_order_acquire);
                    if (canonicalhashes[idx] == canonicalhashes[idx2]) {
                        fully_mapped[fam] = true;
                        logi("Family {} (top sector {}) is symmetric to family {}, sector {}",
                                fam+1, sector, fam2+1, sec2);
                        uint8_t *perm = canonicalperms + idx*maxfamilysize;
                        uint8_t *perm2 = canonicalperms + idx2*maxfamilysize;
                        logd("Fam {} is {}, sec {} perm: {}", fam+1, families[fam],
                                sector, vector<uint8_t>(perm, perm+nx));
                        logd("Fam {} is {}, sec {} perm: {}", fam2+1, families[fam2],
                                sec2, vector<uint8_t>(perm2, perm2+nx));
                        exvector src, dst;
                        for (unsigned i = 0; i < nx; i++) {
                            int p2i = bitposition(sec2, perm2[i]);
                            logd("  {} == {}", families[fam][perm[i]], families[fam2][p2i]);
                            src.push_back(families[fam][perm[i]]);
                            dst.push_back(families[fam2][p2i]);
                        }
                        auto mommap = find_momenta_map(src, dst, loopmomenta);
                        logd("Mom map: {}", mommap);
                        fam2mommap[fam] = to_string(mommap);
                        goto found;
                    }
                }
                logi("Family {} (top sector {}) is unique", fam+1, sector);
                hash2sector[canonicalhashes[idx]] = make_pair(fam, sector);
            found:;
            }
            // Skip JOBS sizes.
            for (int i = 0 ; (i < JOBS) && (famsize != familysizeset.rend()); i++, famsize++);
        }
        if (FORK) {
            logd("Exporting results");
            ofstream f(tmpdir::filename(WORKER));
            f << fam2mommap.size() << "\n";
            for (auto &&fammap: fam2mommap) {
                f << fammap.first << " " << fammap.second << "\n";
            }
        }
    FORK_END;
    if (FORK) {
        logd("Importing worker results");
        for (int worker = 0; worker < JOBS; worker++) {
            ifstream f(tmpdir::filename(worker));
            int nentries; f >> nentries;
            for (int i = 0; i < nentries; i++) {
                int fam; f >> fam;
                string line; f >> ws; getline(f, line);
                fam2mommap[fam] = line;
            }
        }
        tmpdir::remove();
    }
    int ndone = 0;
    for (unsigned i = 0; i < totalsectors; i++) {
        ndone += !!hash_done[i];
    }
    logd("Canonized {} sectors out of {}, {}% of hashes wasted",
        ndone, totalsectors, (hashes_done->load() - ndone)*100/(hashes_done->load()));
    shared_free(canonicalperms, maxfamilysize * sizeof(uint8_t) * totalsectors);
    shared_free(canonicalhashes, sizeof(hash_t) * totalsectors);
    shared_free(hash_done, sizeof(uint8_t) * totalsectors);
    shared_free(fully_mapped, sizeof(uint8_t) * totalsectors);
    cout << "{";
    for (unsigned fam = 0; fam < families.nops(); fam++) {
        if (fam == 0) { cout << "\n "; }
        else { cout << ",\n "; }
        auto it = fam2mommap.find(fam);
        if (it != fam2mommap.end()) {
            cout << it->second;
        } else {
            cout << "{}";
        }
    }
    cout << "\n}\n";
}

#define IFCMD(name, condition) \
    if ((argc >= 1) && !strcasecmp(argv[0], name)) \
        if (!(condition)) { \
            cerr << "feynson: malformed '" << argv[0] \
                 << "' invocation (use -h to see usage)" << endl; \
            return 1; \
        } else
int
main(int argc, char *argv[])
{
    // Disable buffering from cerr to make logs from multiple
    // worker processes not overlap.
    std::ios_base::sync_with_stdio(false);
    cerr.unsetf(std::ios::unitbuf);
    bool SHORT = false;
    for (int opt; (opt = getopt(argc, argv, "hqsCVj:")) != -1;) {
        switch (opt) {
        case 'h': usage(); return 0;
        case 'V': cout << VERSION; return 0;
        case 'q': VERBOSE = false; break;
        case 'C': COLORS = true; break;
        case 's': SHORT = true; break;
        case 'j': JOBS = atoi(optarg); if (JOBS < 1) JOBS = 1; break;
        default: return 1;
        }
    }
    argc -= optind;
    argv += optind;
    IFCMD("ufx", argc == 2) {
        main_ufx(argv[1]);
    }
    else IFCMD("zero-sectors", argc == 2) {
        main_zerosectors(argv[1], SHORT);
    }
    else IFCMD("symmetrize", argc == 2) {
        main_symmetrize(argv[1]);
    }
    else if (argc == 0) {
        cerr << "feynson: no command provided (use -h to see usage)" << endl;
        return 1;
    }
    else {
        cerr << "feynson: unrecognized command '"
             << argv[0]
             << "' (use -h to see usage)" << endl;
        return 1;
    }
}
