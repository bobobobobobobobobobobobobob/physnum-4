#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ConfigFile.tpp"


using namespace std;

const double PI=3.1415926535897932384626433832795028841971e0;
constexpr double epsilon0= 8.8541878e-12;
// Résolution d'un système d'équations linéaires par élimination de
// Gauss-Jordan:
template<class T>
vector<T>
solve(const vector<T>& diag,
      const vector<T>& lower,
      const vector<T>& upper,
      const vector<T>& rhs)
{
    vector<T> solution(diag.size());
    vector<T> new_diag(diag);
    vector<T> new_rhs(rhs);

    for (int i = 1; i < diag.size(); ++i) {
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution[diag.size() - 1] =
      new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i(diag.size() - 2); i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    return solution;
}

//TODO build the epsilon function
double epsilon(const double& r, const double& r1, const double& epsilon_a, const double& epsilon_b) {
    if (r < r1){return epsilon_a ;}
    else{return epsilon_b ;}
}

//TODO build the rho_epsilon function (rho_lib / epsilon_0)
double rho_epsilon(const double& r,const double& r1,const double& rho0,const bool& uniform_rho_case)
{
    if (uniform_rho_case) {
        return rho0 / epsilon0;
    } else {
        if ((r < r1)){return  rho0 * sin(PI * r / r1);}
        else{return 0.0;}
    }
}


int
main(int argc, char* argv[])
{
    

    // USAGE: Exercise4 [configuration-file] [<settings-to-overwrite> ...]

    // Read the default input
    string inputPath = "configuration.in";
    // Optionally override configuration file.
    if (argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    // Override settings
    for (int i = 2; i < argc; i++)
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Read geometrical inputs
    const double R  = configFile.get<double>("R");
    const double r1 = configFile.get<double>("r1");

    const double rho0 = configFile.get<double>("rho0");

    // For the analytical comparison
    const bool uniform_rho_case = configFile.get<bool>("uniform_rho_case");
        
    // Dielectric relative permittivity
    const double epsilon_a = configFile.get<double>("epsilon_a");
    const double epsilon_b = configFile.get<double>("epsilon_b");
    
    // Boundary conditions
    const double VR = configFile.get<double>("VR");
    
    // Discretization
    const int N1 = configFile.get<int>("N1");
    const int N2 = configFile.get<int>("N2");
    
    // Fichiers de sortie:
    string fichier = configFile.get<string>("output");
    string fichier_phi = fichier+"_phi.out";
    string fichier_E   = fichier+"_E.out";
    string fichier_D   = fichier+"_D.out";

    // Create our finite elements
    const int pointCount = N1 + N2 + 1; // Number of finite elements
    const double h1 = (r1 - 0) / N1;
    const double h2 = (R  - r1) / N2;

    // Position of elements
    vector<double> r(pointCount);

    //TODO build the nodes vector r
    for (int i = 0; i <= N1; ++i)
        r[i] = i * h1; // De 0 à r1

    for (int i = 0; i < N2; ++i)
        r[N1 + 1 + i] = r1 + (i + 1) * h2; // De r1 à R
    
    // Arrays initialization
    vector<double> h(pointCount-1); 	// Distance between grid points
    vector<double> midPoint(pointCount-1);  // Midpoint of each grid element
    
    vector<double> rho_libre(pointCount-1,0.0); // nécessaire pour la c calcul de pho_libre 
        
    // TODO build the h vector and midpoint vector
    for (int i = 0; i < pointCount - 1; ++i) 
    {
    h[i] = r[i + 1] - r[i];
    midPoint[i] = (r[i] + r[i + 1]) / 2.0;
    }

    // Construct the matrix and right-hand side
    vector<double> diagonal(pointCount, 0.0);  // Diagonal
    vector<double> lower(pointCount - 1, 0.0); // Lower diagonal
    vector<double> upper(pointCount - 1, 0.0); // Upper diagonal
    vector<double> rhs(pointCount, 0.0);       // Right-hand-side
    
    // Loop over the intervals: add the contributions to matrix and rhs   
    for (int k = 0; k < pointCount-1; ++k) {
        
                   
        // TODO build the vectors diagonal, lower, upper, rhs
        double eps_milieu = epsilon(midPoint[k], r1, epsilon_a, epsilon_b);
        // Contribution à la matrice (schéma standard 1D éléments finis)
        double coeff = eps_milieu*midPoint[k] / (h[k]);

        diagonal[k]     += coeff;
        diagonal[k + 1] += coeff;
        lower[k]        -= coeff;
        upper[k]        -= coeff;
        double rho_milieu = rho_epsilon(midPoint[k], r1, rho0, uniform_rho_case); // rho_lib/epsilon0

        rho_libre[k]=rho_milieu*epsilon0;//pour la question c)

        double b1 = (h[k] / (2.0)) * rho_milieu * midPoint[k];

        rhs[k]     += b1;
        rhs[k + 1] += b1;
    }    

    // TODO boundary condition at r=R (modify the lines below)
    lower[lower.size() - 1]       = 0.0;
    diagonal[diagonal.size() - 1] = 1.0;
    rhs[rhs.size() - 1] = VR;


    // Solve the system of equations
    vector<double> phi = solve(diagonal, lower, upper, rhs);

    // Calculate electric field E and displacement vector D
    vector<double> E(pointCount - 1, 0.0);
    vector<double> D(pointCount - 1, 0.0);
    for (int i = 0; i < E.size(); ++i) {
        // TODO calculate E and D
        E[i] = -(phi[i + 1] - phi[i]) / h[i];
        D[i] = epsilon0*E[i]*epsilon(midPoint[i], r1, epsilon_a, epsilon_b); 
    }


    vector<double> Div_D(pointCount-1,0.0);
    for (int i = 0; i < Div_D.size(); ++i) {
    Div_D[i] = (1.0 / midPoint[i]) * ((r[i+1] * D[i]) - (r[i] * D[i])) / h[i];
    }
   
    // Export data
    {
        // Electric potential phi
        ofstream ofs(fichier_phi);
        ofs.precision(15);

        if (r.size() != phi.size())
            throw std::runtime_error("error when writing potential: r and "
                                     "phi does not have size");

        for (int i = 0; i < phi.size(); ++i) {
            ofs << r[i] << " " << phi[i] << endl;
        }
    }

    {
        // Electric field E
        ofstream ofs(fichier_E);
        ofs.precision(15);

        if (r.size() != (E.size() + 1))
            throw std::runtime_error("error when writing electric field: size of "
                                     "E should be 1 less than r");

        for (int i = 0; i < E.size(); ++i) {
            ofs << midPoint[i] << " " << E[i] << endl;
        }
    }
    {
        // Displacement field D
        ofstream ofs(fichier_D);
        ofs.precision(15);
        if (E.size() != D.size() || D.size() != rho_libre.size() || D.size() != Div_D.size())
        throw std::runtime_error("error when writing displacement field: sizes of vectors do not match");

        for (int i = 0; i < D.size(); ++i) {
            ofs << midPoint[i] << " " << D[i] << " " << rho_libre[i] << " " << Div_D[i] << endl;
        }
        
    }

    return 0;
}

