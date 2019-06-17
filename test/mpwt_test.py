import mpwt
mpwt.remove_pgbds('fatty_acid_beta_oxydation_icyc,tca_cycle_ecolicyc')
mpwt.cleaning_input('test')
mpwt.multiprocess_pwt('test', patho_inference=True, verbose=True)

print('a')