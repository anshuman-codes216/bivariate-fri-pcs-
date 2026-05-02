#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <stdexcept>
#include <string>
using namespace std;

// ─────────────────────────────────────────────────────────────
//  Prime field arithmetic  (mod P)
//  P must be prime and P-1 must be divisible by every M you use.
//  97 works for M in {1,2,4,8,16,32} since 97-1 = 96 = 2^5 * 3.
//  Swap P for a larger prime (e.g. 769, 3329, a 64-bit prime) if
//  you need larger M or higher degree.
// ─────────────────────────────────────────────────────────────
const int P = 97;

int wrap(int a)          { return ((a % P) + P) % P; }
int add (int a, int b)   { return wrap(a + b); }
int sub (int a, int b)   { return wrap(a - b); }
int mul (int a, int b)   { return wrap(a * b); }

int power(int base, int exp) {
    int result = 1;
    base = wrap(base);
    for (; exp > 0; exp >>= 1) {
        if (exp & 1) result = mul(result, base);
        base = mul(base, base);
    }
    return result;
}
int inv(int a) { return power(a, P - 2); }   // Fermat's little theorem

// ─────────────────────────────────────────────────────────────
//  Polynomial type: poly[i] is the coefficient of X^i
// ─────────────────────────────────────────────────────────────
using Poly = vector<int>;

int evalPoly(const Poly& poly, int x) {
    int result = 0, xpow = 1;
    for (int coeff : poly) {
        result = add(result, mul(coeff, xpow));
        xpow   = mul(xpow, x);
    }
    return result;
}

// ─────────────────────────────────────────────────────────────
//  Roots of unity
//  Finds the smallest generator g of (Z/PZ)* then derives an
//  M-th primitive root.  Requires M | (P-1).
// ─────────────────────────────────────────────────────────────
int primitiveRoot() {
    // Brute-force search for a generator of (Z/PZ)*
    // Works fine for the small primes used here.
    for (int g = 2; g < P; g++) {
        bool ok = true;
        // Check that g^((P-1)/q) != 1 for every prime factor q of P-1
        int pm1 = P - 1, tmp = pm1;
        for (int q = 2; q * q <= tmp; q++) {
            if (tmp % q == 0) {
                if (power(g, pm1 / q) == 1) { ok = false; break; }
                while (tmp % q == 0) tmp /= q;
            }
        }
        if (tmp > 1 && power(g, pm1 / tmp) == 1) ok = false;
        if (ok) return g;
    }
    throw runtime_error("No primitive root found — is P prime?");
}

vector<int> buildRoots(int M) {
    if ((P - 1) % M != 0)
        throw runtime_error("M=" + to_string(M) +
                            " does not divide P-1=" + to_string(P-1) +
                            ". Choose M that divides " + to_string(P-1) + ".");
    int g = primitiveRoot();
    int w = power(g, (P - 1) / M);   // primitive M-th root of unity
    vector<int> roots(M);
    roots[0] = 1;
    for (int i = 1; i < M; i++)
        roots[i] = mul(roots[i-1], w);
    return roots;
}

// ─────────────────────────────────────────────────────────────
//  Lagrange basis evaluation:  L_i(y)  over the set `roots`
// ─────────────────────────────────────────────────────────────
int lagrange(int i, int y, const vector<int>& roots) {
    int M = roots.size(), num = 1, den = 1;
    for (int j = 0; j < M; j++) {
        if (j == i) continue;
        num = mul(num, sub(y, roots[j]));
        den = mul(den, sub(roots[i], roots[j]));
    }
    return mul(num, inv(den));
}

// ─────────────────────────────────────────────────────────────
//  FRI folding:  f_alpha(X) = f_even(X) + alpha * f_odd(X)
//  Halves the degree each round until a constant remains.
// ─────────────────────────────────────────────────────────────
Poly fold(const Poly& poly, int alpha) {
    Poly result;
    result.reserve((poly.size() + 1) / 2);
    for (int i = 0; i < (int)poly.size(); i += 2) {
        int even = poly[i];
        int odd  = (i + 1 < (int)poly.size()) ? poly[i+1] : 0;
        result.push_back(add(even, mul(alpha, odd)));
    }
    return result;
}

// Returns the final constant after repeated folding;
// records each random alpha in `alphasUsed`.
int friCommit(const Poly& poly, vector<int>& alphasUsed) {
    Poly cur = poly;
    alphasUsed.clear();
    while (cur.size() > 1) {
        int alpha = (rand() % (P - 1)) + 1;   // alpha in [1, P-1]
        alphasUsed.push_back(alpha);
        cur = fold(cur, alpha);
    }
    return cur[0];
}

// Reproduce the same folding sequence and return the final constant.
int friVerify(const Poly& poly, const vector<int>& alphas) {
    Poly cur = poly;
    for (int alpha : alphas)
        cur = fold(cur, alpha);
    return cur[0];
}

// ─────────────────────────────────────────────────────────────
//  Commit phase: fingerprint each row polynomial
// ─────────────────────────────────────────────────────────────
struct Commit {
    vector<int>         C;       // one fingerprint per row
    vector<vector<int>> alphas;  // folding randomness per row
};

Commit commitAll(const vector<Poly>& rows) {
    Commit cm;
    cout << "\n=== COMMIT ===\n";
    for (int i = 0; i < (int)rows.size(); i++) {
        vector<int> a;
        int c = friCommit(rows[i], a);
        cm.C.push_back(c);
        cm.alphas.push_back(a);
        cout << "  Row " << i << " fingerprint = " << c << "\n";
    }
    return cm;
}

// ─────────────────────────────────────────────────────────────
//  Degree check + fingerprint verification
// ─────────────────────────────────────────────────────────────
bool verifyPoly(const vector<Poly>& rows, const Commit& cm, int T) {
    cout << "\n=== VERIFY POLY ===\n";
    for (int i = 0; i < (int)rows.size(); i++) {
        // Degree check: highest non-zero coefficient must be <= T
        int deg = (int)rows[i].size() - 1;
        while (deg > 0 && rows[i][deg] == 0) deg--;
        if (deg > T) {
            cout << "  FAIL: Row " << i << " has degree " << deg
                 << " > T=" << T << "!\n";
            return false;
        }
        if (friVerify(rows[i], cm.alphas[i]) != cm.C[i]) {
            cout << "  FAIL: Fingerprint mismatch for row " << i << "!\n";
            return false;
        }
    }
    cout << "  PASS: All rows are valid polynomials of degree <= " << T << ".\n";
    return true;
}

// ─────────────────────────────────────────────────────────────
//  Opening:  evaluate the bivariate polynomial
//            F(X,Y) = sum_i  F_i(X) * L_i(Y)
//  at (alpha, beta) and produce a proof.
// ─────────────────────────────────────────────────────────────
struct OpenProof {
    vector<int> fi;   // F_i(alpha) for each row i
    vector<int> Ri;   // L_i(beta)  for each row i
    int z;            // claimed F(alpha, beta)
};

OpenProof openAt(const vector<Poly>& rows, int alpha, int beta,
                 const vector<int>& roots) {
    cout << "\n=== OPEN at alpha=" << alpha << ", beta=" << beta << " ===\n";
    OpenProof proof;
    int z = 0;
    int M = (int)rows.size();
    for (int i = 0; i < M; i++) {
        int fi = evalPoly(rows[i], alpha);
        int Ri = lagrange(i, beta, roots);
        proof.fi.push_back(fi);
        proof.Ri.push_back(Ri);
        z = add(z, mul(fi, Ri));
        cout << "  F_" << i << "(alpha)=" << fi
             << "  L_" << i << "(beta)=" << Ri << "\n";
    }
    proof.z = z;
    cout << "  Claimed F(alpha,beta) = " << z << "\n";
    return proof;
}

// ─────────────────────────────────────────────────────────────
//  Verification: recompute z from the proof and compare
// ─────────────────────────────────────────────────────────────
bool verifyProof(const OpenProof& proof) {
    cout << "\n=== VERIFY ===\n";
    int z_check = 0;
    for (int i = 0; i < (int)proof.fi.size(); i++)
        z_check = add(z_check, mul(proof.fi[i], proof.Ri[i]));

    if (z_check != proof.z) {
        cout << "  FAIL: z mismatch! Got " << z_check
             << " but proof says " << proof.z << "\n";
        return false;
    }
    cout << "  PASS: F(alpha,beta) = " << z_check << "\n";
    return true;
}

// ─────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────
void printDivisors(int n) {
    cout << "  (divisors of " << n << ":";
    for (int d = 1; d <= n; d++)
        if (n % d == 0) cout << " " << d;
    cout << ")\n";
}

// ─────────────────────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────────────────────
int main() {
    srand((unsigned)time(0));

    cout << "=== Generalised Bivariate FRI PCS  (mod " << P << ") ===\n";
    cout << "    (P-1 = " << P-1 << "; M must divide P-1)\n";
    printDivisors(P - 1);

    // ── Parameters ──────────────────────────────────────────
    int M, T;
    cout << "\nNumber of row polynomials M (must divide " << P-1 << "): ";
    cin >> M;
    if (M <= 0 || (P - 1) % M != 0) {
        cout << "Error: M=" << M << " does not divide P-1=" << P-1 << ".\n";
        return 1;
    }

    cout << "Degree bound T (each row poly has degree <= T): ";
    cin >> T;
    if (T < 0) { cout << "Error: T must be >= 0.\n"; return 1; }

    // ── Roots of unity ──────────────────────────────────────
    vector<int> roots;
    try { roots = buildRoots(M); }
    catch (const exception& e) { cout << "Error: " << e.what() << "\n"; return 1; }

    cout << "\nRoots of unity (M=" << M << "): ";
    for (int r : roots) cout << r << " ";
    cout << "\n";

    // ── Polynomial input ────────────────────────────────────
    vector<Poly> rows(M);
    cout << "\nEnter " << M << " polynomials, each with " << T+1
         << " coefficient(s)  [a0 a1 ... a" << T << "]"
         << "  for  a0 + a1*X + ... + a" << T << "*X^" << T << ":\n";
    for (int i = 0; i < M; i++) {
        cout << "F_" << i << ": ";
        rows[i].resize(T + 1);
        for (int j = 0; j <= T; j++) {
            cin >> rows[i][j];
            rows[i][j] = wrap(rows[i][j]);
        }
    }

    // ── Commit ───────────────────────────────────────────────
    Commit cm = commitAll(rows);

    // ── Verify polynomials ──────────────────────────────────
    if (!verifyPoly(rows, cm, T)) {
        cout << "VerifyPoly failed. Stopping.\n";
        return 1;
    }

    // ── Evaluation point ────────────────────────────────────
    int alpha, beta;
    cout << "\nEnter alpha (0-" << P-1 << "): ";
    cin >> alpha;
    alpha = wrap(alpha);

    cout << "Enter beta  (0-" << P-1 << ", ideally not a root of unity): ";
    cin >> beta;
    beta = wrap(beta);

    // ── Open & verify ────────────────────────────────────────
    OpenProof proof = openAt(rows, alpha, beta, roots);

    if (verifyProof(proof))
        cout << "\n>>> SUCCESS: Proof verified!  F(" << alpha << ","
             << beta << ") = " << proof.z << "\n";
    else
        cout << "\n>>> FAILED: Proof is invalid.\n";

    return 0;
}
