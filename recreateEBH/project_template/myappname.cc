#include "itensor/all.h"

using namespace itensor;

int main()
{
    int L = 3;
    auto sites = Boson(L, {"MaxOcc=", 60, "ConserveQNs", true});
    auto sweeps = Sweeps(8);
    sweeps.maxdim() = 40, 80, 400, 400, 800, 800, 1000, 1000;
    sweeps.cutoff() = 1E-16;

    Real gamma = 0.0;
    Real J = 1.0;
    Real g = 0.7;
    Real C6 = -0.25*pow(0.000005,6);
    Real R = 0.000005;
    Real d = 1.0;
    
    Real U = C6/(pow(R,6.0));
    Real V = C6/(pow(d,6.0)+pow(R,6.0));
    Real W = C6/(pow(2.0,6.0)*pow(d,6.0)+pow(R,6.0));

    auto ampo = AutoMPO(sites);
    for (int j=1; j<L; j++)
    {
        ampo += -J, "A", j, "Adag", j + 1;
        ampo += -J, "Adag", j, "A", j + 1;
        ampo += V, "N", j, "N", j + 1;
    }
    ampo += -J, "A", 1, "Adag", L;
    ampo += -J, "Adag", 1, "A", L;
    ampo += V, "N", 1, "N", L;

    for (int j=1; j<=L; j++)
    {
        ampo += -(j-(L+1)/2)*gamma,"N",j;
        ampo += g/2,"N",j,"N",j;
        ampo += -g/2,"N",j;
        ampo += U,"N",j,"N",j;
    }

    for (int j=1; j<=L-2; j+=2)
    {
        ampo += W,"N",j,"N",j+2;
    }
    

    auto H = toMPO(ampo);

    auto state = InitState(sites);
    for (int i=1; i<=L; i++)
    {
        state.set(i, "20");
    }
    auto psi0 = randomMPS(state);

    auto [energy, psi] = dmrg(H, psi0, sweeps, "Quiet");
    
    // Measuring <a_j>
    println("\naa = ");
    for (int j=1; j<L; j++)
    {
        psi.position(j);

        auto ket = psi(j)*psi(j+1);
        auto bra = dag(prime(ket, "Site"));

        auto aa = op(sites, "Adag", j)*op(sites,"A",j+1);
        
        // take inner product
        auto aaij = elt(bra*aa*ket);
        printfln("%d %.12f",j, aaij);
    }

    return 0;
}


