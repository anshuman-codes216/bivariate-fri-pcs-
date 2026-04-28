#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace std;

const int P = 97;
int wrap(int a) { return ((a % P) + P) % P; }
int add(int a, int b) { return wrap(a + b); }
int sub(int a, int b) { return wrap(a - b); }
int mul(int a, int b) { return wrap(a * b); }

int power(int a, int exp) {
    int result = 1;
    a = wrap(a);
    while (exp > 0) {
        if (exp % 2 == 1) result = mul(result, a);
        a = mul(a, a);
        exp /= 2;
    }
    return result;
}
int inv(int a) { return power(a, P - 2); }
using Poly = vector<int>;
int evalPoly(const Poly& poly, int x) {
    int result = 0, xpow = 1;
    for (int coeff : poly) {
        result = add(result, mul(coeff, xpow));
        xpow   = mul(xpow, x);
    }
    return result;
}

vector<int> buildRoots(int M) {
    int w = power(5, (P - 1) / M);
    vector<int> roots(M);
    roots[0] = 1;
    for (int i = 1; i < M; i++)
        roots[i] = mul(roots[i-1], w);
    return roots;
}

int lagrange(int i, int y, const vector<int>& roots) {
    int M = roots.size(), num = 1, den = 1;
    for (int j = 0; j < M; j++) {
        if (j == i) continue;
        num = mul(num, sub(y, roots[j]));
        den = mul(den, sub(roots[i], roots[j]));
    }
    return mul(num, inv(den));
}

Poly fold(const Poly& poly, int alpha) {
    Poly result;
    for (int i = 0; i < (int)poly.size(); i += 2) {
        int even = poly[i];
        int odd  = (i + 1 < (int)poly.size()) ? poly[i+1] : 0;
        result.push_back(add(even, mul(alpha, odd)));
    }
    return result;
}

int friCommit(const Poly& poly, vector<int>& alphasUsed) {
    Poly cur = poly;
    alphasUsed.clear();
    while (cur.size() > 1) {
        int alpha = (rand() % 96) + 1;
        alphasUsed.push_back(alpha);
        cur = fold(cur, alpha);
    }
    return cur[0];  
}

int friVerify(const Poly& poly, const vector<int>& alphas) {
    Poly cur = poly;
    for (int alpha : alphas)
        cur = fold(cur, alpha);
    return cur[0];  
}

struct Commit {
    vector<int>         C;
    vector<vector<int>> alphas;
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

bool verifyPoly(const vector<Poly>& rows, const Commit& cm, int T) {
    cout << "\n=== VERIFY POLY ===\n";
    for (int i = 0; i < (int)rows.size(); i++) {
        if ((int)rows[i].size() - 1 > T) {
            cout << "  FAIL: Row " << i << " degree too high!\n";
            return false;
        }
        if (friVerify(rows[i], cm.alphas[i]) != cm.C[i]) {
            cout << "  FAIL: Fingerprint mismatch for row " << i << "!\n";
            return false;
        }
    }
    cout << "  PASS: All rows are valid.\n";
    return true;
}

struct OpenProof {
    vector<int> fi;
    vector<int> Ri;
    int z;
};

OpenProof openAt(const vector<Poly>& rows, int alpha, int beta,
                 const vector<int>& roots) {
    cout << "\n=== OPEN at alpha=" << alpha << ", beta=" << beta << " ===\n";
    OpenProof proof;
    int z = 0, M = rows.size();
    for (int i = 0; i < M; i++) {
        int fi = evalPoly(rows[i], alpha);
        int Ri = lagrange(i, beta, roots);
        proof.fi.push_back(fi);
        proof.Ri.push_back(Ri);
        z = add(z, mul(fi, Ri));
        cout << "  F_" << i << "(alpha)=" << fi
             << "  R_" << i << "(beta)=" << Ri << "\n";
    }
    proof.z = z;
    cout << "  Claimed F(alpha,beta) = " << z << "\n";
    return proof;
}

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

int main() {
    srand(time(0));

    cout << "=== Simple Bivariate FRI PCS (mod " << P << ") ===\n";

    int M = 4;
    int T = 3;

    vector<int> roots = buildRoots(M);
    cout << "\nRoots of unity (M=" << M << "): ";
    for (int r : roots) cout << r << " ";
    cout << "\n";

    vector<Poly> rows(M);
    cout << "\nEnter " << M << " polynomials, each with " << T+1
         << " coefficients (a0 a1 a2 a3) for a0+a1*X+a2*X^2+a3*X^3:\n";
    for (int i = 0; i < M; i++) {
        cout << "F_" << i << ": ";
        rows[i].resize(T + 1);
        for (int j = 0; j <= T; j++) {
            cin >> rows[i][j];
            rows[i][j] = wrap(rows[i][j]);
        }
    }

    Commit cm = commitAll(rows);

    if (!verifyPoly(rows, cm, T)) {
        cout << "VerifyPoly failed. Stopping.\n";
        return 1;
    }

    int alpha, beta;
    cout << "\nEnter alpha (0-96): "; cin >> alpha; alpha = wrap(alpha);
    cout << "Enter beta  (0-96, not a root): "; cin >> beta; beta = wrap(beta);

    OpenProof proof = openAt(rows, alpha, beta, roots);

    if (verifyProof(proof))
        cout << "\n>>> SUCCESS: Proof verified! F(" << alpha << "," << beta
             << ") = " << proof.z << "\n";
    else
        cout << "\n>>> FAILED: Proof is invalid.\n";

    return 0;
}