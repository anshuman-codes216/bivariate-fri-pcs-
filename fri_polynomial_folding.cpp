#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;

int main() {
    int degree;
    cout << "Enter the degree of your polynomial: ";
    cin >> degree;

    vector<int> poly(degree + 1);
    cout << "Enter your " << degree + 1 << " coefficients: \n";
    for (int i = 0; i <= degree; i++) {
        cin >> poly[i];
        poly[i] = ((poly[i] % 97) + 97) % 97;
    }

    srand(time(0));
    int round = 1;

    while (poly.size() > 2) {
        int alpha = (rand() % 96) + 1;

        cout << "\n--- FRI ROUND " << round << " ---\n";
        cout << "Random Challenge (Alpha): " << alpha << "\n";

        vector<int> next_poly;
        for (int i = 0; i < poly.size(); i += 2) {
            int even_coeff = poly[i];
            int odd_coeff = 0;
            if (i + 1 < poly.size()) {
                odd_coeff = poly[i + 1];
            }

            int new_coeff = even_coeff + (alpha * odd_coeff);
            new_coeff = ((new_coeff % 97) + 97) % 97;
            next_poly.push_back(new_coeff);
        }

        poly = next_poly;
        cout << "New Polynomial: ";
        for (int i = 0; i < poly.size(); ++i) {
            cout << poly[i] << "x^" << i;
            if (i < poly.size() - 1) cout << " + ";
        }
        cout << "\n";

        round++;
    }

    cout << "\nFINISHED! Your polynomial is now low-degree.\n";
    return 0;
}