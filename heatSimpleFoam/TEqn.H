{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    kappaEff.correctBoundaryConditions();

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();
    TEqn.solve();

}
