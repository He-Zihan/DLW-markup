
************************MARKUP_DLW_CD***************************
clear mata
set mem 400m

gen y=ln(Y)
gen m=ln(M)
gen l=ln(L)
gen k=ln(K)
gen wagebill=ln(W)
gen lva=ln(Valueadd)

*-----------------------------------------------------------*
* OLS REGRESSION FOR STARTING VALUES
xi: reg y l m k i.industry i.year
foreach var in l m k {
    qui gen OLS_`var' = _b[`var']
}
qui gen OLSConst = _b[_cons]

//GMM initial value setting
foreach var in l m k {
    qui gen initial_`var' = OLS_`var'
}
qui gen initialConst = OLSConst

*------FIRST STAGE USING EXP AS INPUT-----------------------*
xtset id year
xi: reg y l m k i.industry i.year
predict phi
predict epsilon, res
label var phi "phi_it" 
label var epsilon "measurement error first stage"
gen phi_lag=L.phi

*-------------------------------------------------------------*
gen l_lag = L.l
gen k_lag = L.k
gen m_lag = L.m

*--------------COMPUTE CORRECTED SHARES---------------------*
gen y_c=y-epsilon
gen va_c=exp(y_c)
gen alpha_m=M/va_c

*-----------------------------------------------------------------*
drop _I*
sort id year
gen const=1
drop if missing(y, l_lag, m_lag, k, phi, phi_lag)

*--------------BEGIN MATA PROGRAM -----------------------*
*void  matrix()
 *betas are matrices of coefficients，crit is the criterion function. g and H are the gradient and Hessian
 *st_data  gets the dep var "."
 *st_data(., ("mpg", "weight")) returns the values of variables mpg and weight, all observations.
 *invsym(.)is the inverse of the real symmetry matrix.
mata:

void GMM_DLW_CD(todo, betas, crit, g, H) {
    PHI = st_data(., "phi")
    PHI_LAG = st_data(., "phi_lag")
    Z = st_data(., ("const", "l_lag", "m_lag", "k_lag"))
    X = st_data(., ("const", "l", "m", "k"))
    X_lag = st_data(., ("const", "l_lag", "m_lag", "k_lag"))
    Y = st_data(., "y")
    C = st_data(., "const")

    OMEGA = PHI - X * betas'
    OMEGA_lag = PHI_LAG - X_lag * betas'
    OMEGA_lag_pol = (C, OMEGA_lag)
    g_b = invsym(OMEGA_lag_pol' * OMEGA_lag_pol) * OMEGA_lag_pol' * OMEGA
    XI = OMEGA - OMEGA_lag_pol * g_b
    crit = (Z' * XI)' * (Z' * XI)
}

void DLW_CD() {
    initialvalue = st_data(1, ("initialConst", "initial_l", "initial_m", "initial_k"))
    S = optimize_init()
    optimize_init_evaluator(S, &GMM_DLW_CD())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_technique(S, "nm")
    optimize_init_nmsimplexdeltas(S, 0.1)
    optimize_init_which(S, "min")
    optimize_init_params(S, initialvalue)
    p = optimize(S)
    p
    st_matrix("beta_dlw_cd", p)
}

end

*------------------END MATA PROGRAM------------------*
cap program drop dlw_cd
program dlw_cd, rclass
preserve
sort id year
mata DLW_CD()
end

*--------------COMPUTE MARKUPS -----------------------*
dlw_cd
gen betam_cd = beta_dlw_cd[1,3]
gen Markup_DLWCD = betam_cd / alpha_m
sum Markup_DLWCD
