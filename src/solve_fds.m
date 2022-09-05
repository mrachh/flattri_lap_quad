function X = solve_fds(F,b,Area)
    rhs = b.*sqrt(Area);
    X = srskelf_sv_nn(F,rhs);
    X = X./sqrt(Area);

end