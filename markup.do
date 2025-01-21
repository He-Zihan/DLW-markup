

************************MARKUP DLW***************************
clear mata
set mem 400m

gen y=ln(Y)
gen m=ln(M)
gen l=ln(L)
gen k=ln(K)
gen wagebill=ln(W)
gen lva=ln(Valueadd)
* higher order terms on inputs
local M=3
local N=3
forvalues i=1/`M' {
gen l`i'=l^(`i')
gen m`i'=m^(`i')
gen k`i'=k^(`i')
*interaction terms
forvalues j=1/`N' {
gen l`i'm`j'=l^(`i')*m^(`j')
gen l`i'k`j'=l^(`i')*k^(`j')
gen k`i'm`j'=k^(`i')*m^(`j')
}
}
gen lk=l*k
gen lm=l*m
gen km=k*m
gen lmk=l*m*k
gen l2k2m2=l2*k2*m2
gen l3k3m3=l3*k3*m3
*-----------------------------------------------------------*
* OLS REGRESSION FOR STARTING VALUES
xi:reg y l m k l2 m2 k2 lm lk km lmk i.nace2 i.year
for any l m k l2 m2 k2 lm lk km lmk : qui gen OLSX=_b[X]
qui gen OLSConst=_b[_c]
//GMM initial value setting
for any l m k l2 m2 k2 lm lk km lmk : qui gen initialX=OLSX
qui gen initialConst=OLSConst
*------FIRST STAGE USING EXP AS INPUT-----------------------*
xtset id year
xi: reg y l* m* k* i.nace2 i.year
predict phi
predict epsilon, res
label var phi "phi_it 
label var epsilon "measurement error first stage
gen phi_lag=L.phi
*-------------------------------------------------------------*
gen l_lag=L.l
gen k_lag=L.k
gen m_lag=L.m
gen l_lag2=l_lag^2
gen k_lag2=k_lag^2
gen m_lag2=m_lag^2
gen l_lagm_lag=l_lag*m_lag
gen l_lagk_lag=l_lag*k_lag
gen k_lagm_lag=k_lag*m_lag
gen l_lagm_lagk_lag=l_lag*m_lag*k_lag
gen l_lagk=l_lag*k
gen km_lag=k*m_lag
gen l_lagm_lagk=l_lag*m_lag*k
*--------------COMPUTE CORRECTED SHARES---------------------*
gen y_c=y-epsilon
gen va_c=exp(y_c)
gen alpha_m=M/va_c
*-----------------------------------------------------------------*
drop _I*
sort id year
gen const=1
drop if y==.
drop if l_lag==.
drop if m_lag==.
drop if k==.
drop if phi==.
drop if phi_lag==.
/* FOR SKETCH CODE BELOW WE USE AR1 PROCESS ON PRODUCTIVITY, IN PAPER WE USE POLYNOMIAL EXPANSION AND HAVE HIGHER ORDER TERMS IN OMEGA_LAG_POL */

*--------------BEGIN MATA PROGRAM -----------------------*
*void  matrix()
 *betas are matrices of coefficients，crit is the criterion function. g and H are the gradient and Hessian
 *st_data  gets the dep var "."
 *st_data(., ("mpg", "weight")) returns the values of variables mpg and weight, all observations.
 *invsym(.)is the inverse of the real symmetry matrix.
mata:

void GMM_DLW_TL(todo,betas,crit,g,H)
{
	PHI=st_data(.,("phi"))
	PHI_LAG=st_data(.,("phi_lag"))
	Z=st_data(.,("const","l_lag","m_lag","k","l_lag2","m_lag2","k2","l_lagm_lag","l_lagk","km_lag","l_lagm_lagk"))
	X=st_data(.,("const","l","m","k","l2","m2","k2","lm","lk","km","lmk"))
	X_lag=st_data(.,("const","l_lag","m_lag","k_lag","l_lag2","m_lag2","k_lag2","l_lagm_lag","l_lagk_lag","k_lagm_lag","l_lagm_lagk_lag"))
	Y=st_data(.,("y"))
	C=st_data(.,("const"))
	
	OMEGA=PHI-X*betas'
	OMEGA_lag=PHI_LAG-X_lag*betas'
	OMEGA_lag_pol=(C,OMEGA_lag)
	g_b = invsym(OMEGA_lag_pol'OMEGA_lag_pol)*OMEGA_lag_pol'OMEGA
	XI=OMEGA-OMEGA_lag_pol*g_b
	crit=(Z'XI)'(Z'XI)
}

void DLW_TRANSLOG()
{
	initialvalue=st_data(1,("initialConst","initiall","initialm","initialk","initiall2","initialm2","initialk2","initiallm","initiallk","initialkm","initiallmk"))
	S=optimize_init()
	optimize_init_evaluator(S, &GMM_DLW_TL())
	optimize_init_evaluatortype(S,"d0")
	optimize_init_technique(S, "nm")
	optimize_init_nmsimplexdeltas(S, 0.1)
	optimize_init_which(S,"min")
	optimize_init_params(S,initialvalue)
	p=optimize(S)
	p
	st_matrix("beta_dlwtranslog",p)
}

end
*------------------END MATA PROGRAM------------------*
cap program drop dlw_translog
program dlw_translog, rclass
preserve
sort id year
mata DLW_TRANSLOG()
end
*--------------COMPUTE MARKUPS -----------------------*
// OLS, DLW estimates for variable input (here materials) are used to compute markup distribution
*---------------------ACF estimates---------------------*
dlw_translog
gen betam_tl1=beta_dlwtranslog[1,3]
gen betam_tl2=beta_dlwtranslog[1,6]
gen betalm_tl=beta_dlwtranslog[1,8]
gen betakm_tl=beta_dlwtranslog[1,10]
gen betalmk_tl=beta_dlwtranslog[1,11]
gen betam_tl=betam_tl1+2*betam_tl2*m+betalm_tl*l+betakm_tl*k+betalmk_tl*l*k

//"const"1,"l"2,"m"3,"k"4,"l2"5,"m2"6,"k2"7,"lm"8,"lk"9,"km"10,"lmk"11
*Calculate the markup
gen Markup_DLWTL=betam_tl/alpha_m
sum Markup_DLWTL

